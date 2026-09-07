"""Explicit ambient-reservoir flux boundary for a reduced thermal moment model.

This is a weak two-thermal-moment approximation, not pointwise satisfaction
of the Gaussian gradient law at the edge. The reduced buoyancy-work ledger
is mandatory and is NOT a full pressure/gravity/TKE energy model.
"""

import math
import numpy as np

from .energy_crosswind import IndependentEnergyCrosswind
from .buoyancy_profile import BuoyancyConstrainedEnthalpySection
from .transverse_mixing import ConservativeTransverseMixing
from .shear_thermal import ReducedShearThermal, affine_root


def encode_section(state, beta):
    BuoyancyConstrainedEnthalpySection._physical(state)
    if not math.isfinite(beta) or beta <= 0.:
        raise ValueError("thermal width ratio must be positive and finite")
    r, y, a, theta, u, x, z = state
    return np.array([math.log(r), math.log(r*y), math.log(a), theta, math.log(u), x, z, math.log(beta)])


def decode_section(parameters):
    par = np.asarray(parameters, float)
    if par.shape != (8,) or not np.all(np.isfinite(par)):
        raise ValueError("supply eight finite section parameters")
    r, c, a, u, beta = np.exp(par[[0, 1, 2, 4, 7]])
    state = np.array([r, c/r, a, par[3], u, par[5], par[6]])
    BuoyancyConstrainedEnthalpySection._physical(state)
    if not math.isfinite(beta) or beta <= 0.:
        raise ValueError("decoded thermal width is outside the finite domain")
    return state, float(beta)


class _PhaseForceView:
    """Forward the established source law, replacing only force quadrature."""
    def __init__(self, mixing, force):
        self.mixing, self.force = mixing, force

    def __getattr__(self, name):
        return getattr(self.mixing.section, name)

    def buoyancy_force(self, state):
        if not np.array_equal(state, self.mixing.state):
            raise ValueError("fixed-section force cannot evaluate a different state")
        return self.force


class ReservoirThermalMoments:
    def __init__(self, mixing, *, thermal_species_ratio, mechanical_work):
        if mechanical_work != "reduced_buoyancy_work":
            raise ValueError("explicitly select the reduced buoyancy-work energy convention")
        if not math.isfinite(thermal_species_ratio) or thermal_species_ratio <= 0.:
            raise ValueError("specify a positive finite thermal/species diffusivity ratio")
        if mixing.section.energy_transport != "total" or mixing.section.ground_interaction != "free":
            raise NotImplementedError("reservoir moments require a free total-energy section")
        self.mixing, self.ratio = mixing, thermal_species_ratio
        self.shear = ReducedShearThermal(mixing, axial_force_density=self.axial_force_density)

    def axial_force_density(self, q, data):
        return 9.81*(self.mixing.section.rhoa-data["rho"])*math.sin(self.mixing.state[3])

    def source_ledger(self, order=16):
        m = self.mixing
        def values(q):
            d = m.local(q)
            force = 9.81*m.area*(m.section.rhoa-d["rho"])
            return np.column_stack([force, force*d["u"]*math.sin(m.state[3])])
        force, work = m.partition.integrate(values, order, square=True)
        base = IndependentEnergyCrosswind.source_terms(_PhaseForceView(m, float(force)), m.state)
        source = base.copy()
        source[4] += work
        return dict(original=base, sources=source, buoyancy_force=float(force), buoyancy_work=float(work))

    def edge_fields(self, q, family, *, order=16):
        d = self.shear.radial_fields(q, family, thermal_species_ratio=self.ratio, order=order)
        ea = .5*self.mixing.wind**2
        total = ea*d["mass_flux"]
        heat = total-d["kinetic_flux"]
        return dict(**d, reservoir_total=total, reservoir_enthalpy=heat,
                    enthalpy_gradient_defect=d["enthalpy_flux"]-heat,
                    species_gradient_defect=d["species_flux"],
                    axial_momentum_gradient_defect=d["momentum_flux"]
                        -self.mixing.wind*math.cos(self.mixing.state[3])*d["mass_flux"])

    def budgets(self, family, *, order=16):
        m, p = self.mixing, self.mixing.partition
        b = self.shear.equilibrium_budgets(family, thermal_species_ratio=self.ratio, order=order)
        spread = m.sy*m.sy+m.sn*m.sn
        def edge(q):
            d = self.edge_fields(q, family, order=order)
            return np.stack([d["reservoir_enthalpy"], spread*q[:, None]*d["reservoir_enthalpy"],
                             d["reservoir_total"], d["mass_flux"], d["species_flux"]], axis=1)
        eh, eh2, ee, em, ec = 16*p.q0*p.edge_integrate(edge, order)
        return dict(**b, reservoir_heat_outward=eh, reservoir_weighted_heat_outward=eh2,
            reservoir_energy_outward=ee, mass_outward=em, species_outward=ec,
            weak_zeroth_residual=b["enthalpy_axial"]+eh-b["production"],
            weak_second_residual=b["weighted_enthalpy_axial"]+eh2-b["weighted_transverse"]-b["weighted_production"])

    def evaluate(self, *, order=16, probes=257):
        m, p = self.mixing, self.mixing.partition
        ledger = self.source_ledger(order)
        family = m.tangent_family(ledger["sources"], order=order)
        b = self.budgets(family, order=order)
        gamma = affine_root(b["weak_second_residual"])
        at, rates = np.array([1., gamma]), family.at(gamma)
        q = np.unique(np.r_[np.linspace(0., p.qmax, probes), .5*(p.knots[:-1]+p.knots[1:])])
        d = self.shear.radial_fields(q, family, thermal_species_ratio=self.ratio, order=order)
        chi_c, chi_p = d["chi_species"]@at, d["chi_momentum"]@at
        edge = self.edge_fields(np.linspace(p.q0, p.qmax, max(129, (probes+1)//2)), family, order=order)
        fm = edge["mass_flux"]@at
        curvature = abs(rates[3])*math.sqrt(2*p.q0)*m.sn
        scalar = {key: float(value@at) for key, value in b.items()}
        scale0 = max(1., *(abs(scalar[k]) for k in ("enthalpy_axial", "reservoir_heat_outward", "production")))
        scale2 = max(1., *(abs(scalar[k]) for k in ("weighted_enthalpy_axial", "reservoir_weighted_heat_outward", "weighted_transverse", "weighted_production")))
        error = max(abs(scalar["weak_zeroth_residual"])/scale0, abs(scalar["weak_second_residual"])/scale2)
        incoming = bool(np.all(fm < 0.))
        positive = bool(np.all(chi_c >= 0.) and np.all(chi_p >= 0.))
        valid = bool(family.maximum_scaled_residual <= 1e-8 and error <= 1e-5
                     and incoming and positive and curvature < .1)
        iq = p.integrate(lambda qq: qq*m.area*m.local(qq)["h"]*m.local(qq)["u"], order, square=True)
        spread_rate = 2*m.sy*m.sy*(m.log_sy_partials@rates)+2*m.sn*m.sn*(m.log_sn_partials@rates)
        moment_rate = scalar["weighted_enthalpy_axial"]+spread_rate*iq
        # Defects are diagnostic: weak natural BC does NOT make pointwise
        # constitutive gradients at the edge correct.
        hscale = np.maximum(abs((edge["local"]["h"]/edge["local"]["rho"])*fm), 1.)
        cscale = np.maximum(abs(edge["local"]["y"]*fm), 1e-12)
        pscale = np.maximum(abs(edge["local"]["u"]*fm), 1.)
        defects = dict(heat=float(np.max(abs(edge["enthalpy_gradient_defect"]@at)/hscale)),
            species=float(np.max(abs(edge["species_gradient_defect"]@at)/cscale)),
            momentum=float(np.max(abs(edge["axial_momentum_gradient_defect"]@at)/pscale)))
        return dict(ledger=ledger, family=family, budgets=b, scalar_budgets=scalar, gamma=gamma,
            rates=rates, moment_rate=float(moment_rate), weak_budget_scaled_error=float(error),
            minimum_chi_species=float(min(chi_c)), minimum_chi_momentum=float(min(chi_p)),
            maximum_outward_mass=float(max(fm)), incoming=incoming, positive_diffusion=positive,
            curvature_half_width=float(curvature), edge_gradient_defects=defects, valid=valid)


class ReservoirShortSegment:
    """Opt-in research driver; not the default plume/receptor implementation."""
    def __init__(self, jetplume, thermodynamics, *, thermal_species_ratio, mechanical_work):
        self.jp, self.th = jetplume, thermodynamics
        self.ratio, self.work = thermal_species_ratio, mechanical_work

    def evaluate(self, parameters, *, order=8, probes=257):
        state, beta = decode_section(parameters)
        section = BuoyancyConstrainedEnthalpySection(self.jp, self.th, thermal_width_ratio=beta)
        mixing = ConservativeTransverseMixing(section, state)
        result = ReservoirThermalMoments(mixing, thermal_species_ratio=self.ratio,
                                         mechanical_work=self.work).evaluate(order=order, probes=probes)
        if not result["valid"]:
            raise ValueError("short segment left the valid weak closure / inflow / diffusion domain")
        result["mixing"] = mixing
        return result

    @staticmethod
    def flux_values(mixing, *, order=16):
        spread = mixing.sy**2+mixing.sn**2
        theta = mixing.state[3]
        def values(q):
            d = mixing.local(q)
            rho, c, h, u = [d[k] for k in ("rho", "c", "h", "u")]
            return mixing.area*np.column_stack([rho*u, c*u, rho*u*u*math.cos(theta),
                rho*u*u*math.sin(theta), h*u+.5*rho*u**3, spread*q*h*u])
        return mixing.partition.integrate(values, order, square=True)

    def integrate(self, initial, length, steps):
        if (not math.isfinite(length) or length <= 0. or isinstance(steps, bool)
                or int(steps) != steps or steps < 1):
            raise ValueError("positive finite segment length and positive integer steps required")
        state = np.r_[np.asarray(initial, float), np.zeros(6)]
        ds = length/int(steps)
        worst_error, worst_defect, smallest_chi = 0., 0., math.inf
        calls = 0
        def rhs(value):
            nonlocal calls, worst_error, worst_defect, smallest_chi
            out = self.evaluate(value[:8])
            calls += 1
            worst_error = max(worst_error, out["weak_budget_scaled_error"])
            worst_defect = max(worst_defect, out["edge_gradient_defects"]["heat"])
            smallest_chi = min(smallest_chi, out["minimum_chi_species"], out["minimum_chi_momentum"])
            return np.r_[out["rates"], out["ledger"]["sources"], out["moment_rate"]]
        for _ in range(int(steps)):
            first = rhs(state)
            state += ds*rhs(state+.5*ds*first)
        return dict(parameters=state[:8], cumulative_sources=state[8:], rhs_calls=calls,
                    maximum_weak_residual=worst_error, maximum_edge_heat_defect=worst_defect,
                    minimum_sampled_chi=smallest_chi)
