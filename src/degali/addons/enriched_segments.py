"""Guarded local continuation under explicit, conditional mixing scalings.

Not a validated LH2 turbulence closure, long-range solver, or default model.
Every RHS rebuilds the current phase mesh and solves all modal rates.
"""

import copy
import math
import numpy as np
from .buoyancy_profile import BuoyancyConstrainedEnthalpySection
from .transverse_mixing import ConservativeTransverseMixing
from .enriched_transport import PrescribedRadialMixing, EnrichedModalTransport, advective_moments
from .coupled_shape_initialization import FixedMeshCoupledTransport


POLICIES = ("equilibrium_velocity_width", "constant_geometric_diffusivity")


def encode_enriched(projection, parameters):
    r, y, a, theta, u, x, z = projection.mixing.state
    par = np.asarray(parameters, float)
    if par.shape != (projection.count,) or not np.all(np.isfinite(par)):
        raise ValueError("supply all finite scalar modes")
    return np.r_[math.log(r), math.log(r*y), math.log(a), theta, math.log(u), x, z, par]


def decode_enriched(encoded, count):
    value = np.asarray(encoded, float)
    if value.shape != (7+count,) or not np.all(np.isfinite(value)):
        raise ValueError("invalid finite gauge-free state")
    with np.errstate(over="raise", invalid="raise"):
        try:
            rho, c, area, u = np.exp(value[[0, 1, 2, 4]])
        except FloatingPointError as exc:
            raise ValueError("encoded physical state overflows") from exc
    state = np.array([rho, c/rho, area, value[3], u, value[5], value[6]])
    BuoyancyConstrainedEnthalpySection._physical(state)
    return state, value[7:].copy()


def shifted_field_view(template, encoded):
    """Current fields/geometry only; the old baseline remains an origin record."""
    state, par = decode_enriched(encoded, template.count)
    original = template.section
    section = BuoyancyConstrainedEnthalpySection(original.jetplume, original.thermodynamics,
        thermal_width_ratio=original.thermal_width_ratio, quadrature_points=original.quadrature_points)
    section._quadrature_cache = original._quadrature_cache
    mixing = ConservativeTransverseMixing(section, state)
    view = copy.copy(template)
    view.section, view.mixing = section, mixing
    view.hc = section.phase_inverse.enthalpy_and_slope(state[0], state[0]*state[1])[0]
    view._grids, view._edges = {}, {}
    # No fit/retraction/old edge method may be used on this view. In particular
    # template.target is the ORIGINAL six moments, not a downstream constraint.
    return view, par


def scaled_mixing(reference, initial_state, current_state, policy):
    if not isinstance(reference, PrescribedRadialMixing) or policy not in POLICIES:
        raise ValueError("supply a resolved mixing input and select an explicit update policy")
    BuoyancyConstrainedEnthalpySection._physical(initial_state)
    BuoyancyConstrainedEnthalpySection._physical(current_state)
    ratio_a = initial_state[2]/current_state[2]
    factor = ((current_state[4]/initial_state[4])*math.sqrt(ratio_a)
              if policy == "equilibrium_velocity_width" else ratio_a)
    if not math.isfinite(factor) or factor <= 0.:
        raise ValueError("mixing update is not finite/positive")
    result = copy.copy(reference)
    result.knots = reference.knots.copy()
    result.coefficients = reference.coefficients.copy()
    result.coefficients[:, :, 0] *= factor
    # Column1 is an origin flux scale, NOT a physical downstream mass flux.
    return result, float(factor)


def forward_shape_horizon(projection, parameters, rates, *, samples=129):
    parameters, rates = np.asarray(parameters, float), np.asarray(rates, float)
    if (parameters.shape != (projection.count,) or rates.shape != (7+projection.count,)
            or not np.all(np.isfinite(parameters)) or not np.all(np.isfinite(rates))):
        raise ValueError("finite shape and rates of matching size required")
    t = np.linspace(0., 1., samples)
    a, b = np.meshgrid(t, t, indexing="ij")
    psi = projection.prepare(a.ravel(), b.ravel())["psi"]
    size = projection.basis.size
    horizons = []
    for start in (0, size):
        value = psi@parameters[start:start+size]
        speed = psi@rates[7+start:7+start+size]
        if max(abs(value)) > .1000000001:
            raise ValueError("initial shape is already outside trust region")
        mask = abs(speed) > 1e-12
        horizons.append(float(np.min((.1-np.sign(speed[mask])*value[mask])/abs(speed[mask]), initial=np.inf)))
    return min(horizons)


def assert_transport_domain(out):
    scalars = [*out["edge_defects"].values(), out["minimum_chi_momentum"], out["maximum_outward_mass"],
        out["curvature_half_width"], out["linear_scaled_error"], out["weak_heat_scaled_error"], out["mass_boundary_error"]]
    if not np.all(np.isfinite(scalars)) or not np.all(np.isfinite(out["source"])):
        raise ValueError("nonfinite transport diagnostic")
    if max(out["edge_defects"].values()) > .05:
        raise ValueError(f"pointwise edge gate left: {out['edge_defects']}")
    if out["minimum_chi_momentum"] < 0.:
        raise ValueError(f"negative inferred momentum mixing: {out['minimum_chi_momentum']}")
    if out["maximum_outward_mass"] >= 0. or out["curvature_half_width"] >= .1:
        raise ValueError("inflow/curvature gate left")
    if max(out["linear_scaled_error"], out["weak_heat_scaled_error"]) > 1e-8:
        raise ValueError("local transport numerical gate left")
    if abs(out["mass_boundary_error"])/max(abs(out["source"][0]), 1.) > 1e-8:
        raise ValueError("local mass boundary balance gate left")


class EnrichedShortSegment:
    def __init__(self, projection, parameters, *, scalar_mixing, mixing_update,
                 thermal_species_ratio, mechanical_work):
        # Reuse the established validation of all constitutive assumptions.
        EnrichedModalTransport(projection, parameters, scalar_mixing=scalar_mixing,
            thermal_species_ratio=thermal_species_ratio, mechanical_work=mechanical_work)
        if mixing_update not in POLICIES:
            raise ValueError("an explicit supported mixing update is required")
        self.template, self.initial = projection, encode_enriched(projection, parameters)
        self.reference, self.policy = scalar_mixing, mixing_update
        self.ratio, self.work = thermal_species_ratio, mechanical_work
        self.count = len(self.initial)

    def context(self, encoded):
        p, par = shifted_field_view(self.template, encoded)
        supplied, factor = scaled_mixing(self.reference, self.template.mixing.state, p.mixing.state, self.policy)
        return p, par, supplied, factor

    def evaluate(self, encoded, *, order=8, angular_order=48, exact=False, enforce=True):
        p, par, supplied, factor = self.context(encoded)
        if exact:
            out = EnrichedModalTransport(p, par, scalar_mixing=supplied,
                thermal_species_ratio=self.ratio, mechanical_work=self.work).assemble(order=order, angular_order=angular_order)
        else:
            out = FixedMeshCoupledTransport(p, par, scalar_mixing=supplied,
                thermal_species_ratio=self.ratio, mechanical_work=self.work,
                order=order, angular_order=angular_order).evaluate(par)
        t = np.linspace(0., 1., 129)
        a, b = np.meshgrid(t, t, indexing="ij")
        fields = p.fields(par, p.prepare(a.ravel(), b.ravel()))
        if np.any(fields["rho"] <= 0.) or np.any(fields["y"] <= 0.) or np.any(fields["y"] >= 1.):
            raise ValueError("physical density/species domain left")
        out["mixing_amplitude_ratio"] = factor
        out["maximum_log_shape"] = max(float(max(abs(fields[k]))) for k in ("log_c", "log_h"))
        if enforce:
            assert_transport_domain(out)
        return out

    def flux_values(self, encoded, *, order=8, angular_order=8):
        p, par = shifted_field_view(self.template, encoded)
        return advective_moments(p, par, order=order, angular_order=angular_order)[:5]

    def integrate(self, length, steps, *, callback=None):
        if (not math.isfinite(length) or length <= 0. or isinstance(steps, bool)
                or int(steps) != steps or steps < 1):
            raise ValueError("positive finite length and integer steps required")
        count = self.count
        value = np.r_[self.initial, np.zeros(5)]
        ds = length/int(steps)
        history, calls, cached = [], 0, None
        def rhs(v):
            nonlocal calls
            out = self.evaluate(v[:count])
            calls += 1
            return np.r_[out["rates"], out["source"]], out
        for j in range(int(steps)):
            stage = "start"
            try:
                first, out1 = rhs(value) if cached is None else cached
                stage = "midpoint"
                middle, out2 = rhs(value+.5*ds*first)
                candidate = value+ds*middle
                stage = "endpoint"
                endpoint = self.evaluate(candidate[:count])
                calls += 1
                cached = np.r_[endpoint["rates"], endpoint["source"]], endpoint
                value = candidate
                row = dict(step=j+1, reached_m=(j+1)*ds,
                    edge_defects={k: max(v["edge_defects"][k] for v in (out1, out2, endpoint)) for k in endpoint["edge_defects"]},
                    minimum_chi_momentum=min(v["minimum_chi_momentum"] for v in (out1, out2, endpoint)),
                    mixing_amplitude_ratio=endpoint["mixing_amplitude_ratio"],
                    maximum_log_shape=max(v["maximum_log_shape"] for v in (out1, out2, endpoint)))
                history.append(row)
                if callback:
                    callback(row)
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                return dict(completed=False, reached_m=j*ds, requested_length_m=length, divisions=int(steps),
                    failure_stage=stage, failure=str(exc), parameters=value[:count],
                    cumulative_sources=value[count:], history=history, rhs_calls=calls)
        return dict(completed=True, reached_m=length, requested_length_m=length, divisions=int(steps),
            parameters=value[:count], cumulative_sources=value[count:], history=history, rhs_calls=calls)
