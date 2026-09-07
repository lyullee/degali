"""Liquid hydrogen releases: one call, with the uncertainty attached.

Everything needed to assess an LH₂ release is in this package, but it is
spread across a flashing source term, a jet model, a mixing table and a
buoyant plume, each with its own configuration.  This module puts them behind
one function and, as importantly, reports how far each answer has been
checked.

.. code-block:: python

    from degali.lh2 import assess

    result = assess(rate=0.285, orifice=0.0254, storage_pressure=6.0,
                    wind=2.5, height=0.5)
    result.distance_to_lfl        # metres to 4 mole per cent
    result.regime                 # 'aloft', 'low' or 'grounded'
    print(result.report())

The separately validated quiescent near-nozzle research path is exposed as
``lh2_source_from_measured_throat`` plus ``run_lh2_near_field_research``. It
can be joined to JETPLU for horizontal wind-aligned releases with
``run_lh2_crosswind_research``. That interface has been independently tested
on seven PRESLHY releases. Its local-shear mechanism is supported, but the
    complete research path is not promoted. The later independent-energy model
    closes all seven interfaces and improves variance, FAC2 and width, but
    worsens mean bias and centre-height error. It therefore does not silently
    replace the atmospheric assessment above.

Every number carries a confidence, because they are not equally well
established.  The confidences are not written here: they live in
:mod:`degali.evidence` and are read from it, so that this docstring cannot
drift away from what the code prints.

What this will not do
---------------------
It will not give a defensible number beyond the range it was checked in, and
that range is **different for a jet and for a pool**.  The jet evidence is the
PRESLHY E3.5 campaign: 0.08 to 0.29 kg/s through 6 to 25 mm orifices, wind 0.6
to 4.2 m/s, and an array reaching 6 m.  The pool evidence is the NASA Langley
spills: 9 to 10 kg/s from one 9.1 m pond, wind 1.6 to 6.3 m/s, measured at a
tower row 33.8 m away.  Treating those as one range -- which this module did
until the two campaigns were separated -- reports a jet asked about 30 m as
inside the checked range when nothing has checked it.

Outside the relevant range it still runs, and :attr:`Assessment.warnings` says
so.  The warning is raised against the furthest distance a *reported* answer
depends on, not against the integration limit: asking about 30 m and
integrating to 100 m is one question, and warning about the 100 m is an answer
to a question nobody asked.

It will also not give a defensible *concentration* for a release the wind
steers rather than its own momentum.  On the PRESLHY 1 barg trials, where the
exit velocity is 4 to 13 m/s against a 1.5 to 4 m/s wind, nominally identical
releases gave arc maxima of 83 % and 4 %: the plume went wherever the wind
was pointing, and a steady jet model has nothing to say about that. The
warning fires below a velocity ratio of ten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .addons.axisymmetric_jet import (
        AxisymmetricJetResult,
        AxisymmetricJetSource,
        ConservedGaussianJet,
    )
    from .addons.notional import MeasuredThroatExpansion
    from .core.jetplume import JetPlume, JetResult

from .evidence import (
    LIFTOFF_HEIGHT,
    MOMENTUM_RATIO,
    NEAR_FIELD_CORRECTED,
    NEUTRAL_BUOYANCY,
    check_range,
    RANGE as VALIDATED,
)
#: Hydrogen's flammable limits, mole fraction.
LFL = 0.04
UFL = 0.75
STOICHIOMETRIC = 0.295


@dataclass(frozen=True)
class LH2ExpandedSource:
    """Measured throat state and its conserved atmospheric source plane."""

    source: "AxisymmetricJetSource"
    expansion: "MeasuredThroatExpansion"


@dataclass
class LH2NearFieldResearchResult:
    """Dry-air conserved near-field result with numerical audit attached.

    This result belongs to the quiescent, momentum-dominated Raman research
    path. It is intentionally distinct from :class:`Assessment`, which is the
    atmospheric trajectory and consequence-assessment path.
    """

    source: "AxisymmetricJetSource"
    solution: "AxisymmetricJetResult"
    model: "ConservedGaussianJet" = field(repr=False)
    maximum_boundary_residual: float
    maximum_species_drift: float
    maximum_energy_drift: float
    warnings: list[str] = field(default_factory=list)
    notes: dict[str, str] = field(default_factory=dict)

    @property
    def conservative(self) -> bool:
        """Whether the run meets the frozen numerical conservation limits."""
        return (
            self.maximum_boundary_residual < 1.0e-8
            and self.maximum_species_drift < 2.0e-4
            and self.maximum_energy_drift < 2.0e-4
        )

    def report(self) -> str:
        lines = [
            "LH2 conserved axisymmetric near field (research)",
            f"  stations                    : {len(self.solution.S)}",
            f"  streamline extent           : {self.solution.S[-1]:.4f} m",
            f"  centre temperature range    : "
            f"{self.solution.temperature.min():.2f}--"
            f"{self.solution.temperature.max():.2f} K",
            f"  boundary residual           : "
            f"{self.maximum_boundary_residual:.3e}",
            f"  species-flux drift          : {self.maximum_species_drift:.3e}",
            f"  energy-flux drift           : {self.maximum_energy_drift:.3e}",
            f"  conservation screen         : "
            f"{'pass' if self.conservative else 'fail'}",
        ]
        if self.warnings:
            lines.append("  applicability warnings:")
            lines += [f"    - {warning}" for warning in self.warnings]
        return "\n".join(lines)


@dataclass
class LH2CrosswindHandoff:
    """Audited boundary between the conserved near field and JETPLU."""

    model: "JetPlume" = field(repr=False)
    state: np.ndarray
    streamline_distance: float
    target_fluxes: dict[str, float]
    projected_fluxes: dict[str, float]
    relative_residuals: dict[str, float]
    target_halfwidth: float
    projected_halfwidth: float
    halfwidth_residual: float
    target_temperature: float
    projected_temperature: float
    temperature_residual: float
    thermodynamic_closure: str
    virtual_source_temperature: float
    energy_quadrature_residual: float
    energy_quadrature_points: int
    crosswind_entrainment: str
    accepted: bool
    failure_reasons: list[str] = field(default_factory=list)

    def report(self) -> str:
        """Human-readable conservation and compatibility audit."""
        virtual_temperature = (
            f"{self.virtual_source_temperature:.3f} K"
            if math.isfinite(self.virtual_source_temperature) else "n/a"
        )
        lines = [
            "LH2 near-field to crosswind handoff",
            f"  streamline station          : {self.streamline_distance:.5f} m",
            f"  total-mass residual         : "
            f"{self.relative_residuals['total_mass']:.3e}",
            f"  hydrogen residual           : "
            f"{self.relative_residuals['hydrogen']:.3e}",
            f"  horizontal-momentum residual: "
            f"{self.relative_residuals['momentum_x']:.3e}",
            f"  vertical-momentum residual  : "
            f"{self.relative_residuals['momentum_z']:.3e}",
            f"  total-energy residual       : "
            f"{self.relative_residuals['energy']:.3e}",
            f"  H2 half-width residual      : {self.halfwidth_residual:.3e}",
            f"  centre-temperature residual : {self.temperature_residual:.3f} K",
            f"  thermodynamic closure       : {self.thermodynamic_closure}",
            f"  virtual source temperature  : "
            f"{virtual_temperature}",
            f"  energy-quadrature residual  : "
            f"{self.energy_quadrature_residual:.3e}",
            f"  energy-quadrature points    : "
            f"{self.energy_quadrature_points}",
            f"  crosswind entrainment       : {self.crosswind_entrainment}",
            f"  handoff screen              : "
            f"{'pass' if self.accepted else 'fail'}",
        ]
        if self.failure_reasons:
            lines.append("  failure reasons:")
            lines += [f"    - {reason}" for reason in self.failure_reasons]
        return "\n".join(lines)

    def run(
        self, *, distmx: float, tol: float = 1.0e-4,
        smax: float = 1.0e5,
    ) -> "JetResult":
        """Run JETPLU only after the frozen handoff screen has passed."""
        if not self.accepted:
            raise RuntimeError(
                "crosswind integration refused because the conservative "
                "handoff screen failed"
            )
        return self.model.run(self.state.copy(), distmx=distmx, tol=tol, smax=smax)


@dataclass
class LH2IndependentEnergyInterface:
    """Five-flux audit and entry point for the seven-state research model."""

    model: object = field(repr=False)
    state: np.ndarray
    streamline_distance: float
    target_fluxes: dict[str, float]
    projected_fluxes: dict[str, float]
    relative_residuals: dict[str, float]
    target_halfwidth: float
    projected_halfwidth: float
    halfwidth_residual: float
    target_temperature: float
    projected_temperature: float
    temperature_residual: float
    target_buoyancy_force: float
    projected_buoyancy_force: float
    buoyancy_force_ratio: float
    energy_quadrature_residual: float
    energy_quadrature_points: int
    energy_transport: str
    houf_width_mapping: str
    ground_interaction: str
    velocity_spreading_ratio: float
    accepted: bool
    failure_reasons: list[str] = field(default_factory=list)

    def report(self) -> str:
        """Human-readable five-flux and profile-compatibility audit."""
        lines = [
            "LH2 independent-energy crosswind interface",
            f"  streamline station          : {self.streamline_distance:.5f} m",
            f"  total-mass residual         : "
            f"{self.relative_residuals['total_mass']:.3e}",
            f"  hydrogen residual           : "
            f"{self.relative_residuals['hydrogen']:.3e}",
            f"  horizontal-momentum residual: "
            f"{self.relative_residuals['momentum_x']:.3e}",
            f"  vertical-momentum residual  : "
            f"{self.relative_residuals['momentum_z']:.3e}",
            f"  {self.energy_transport}-energy residual: "
            f"{self.relative_residuals['energy']:.3e}",
            f"  H2 half-width residual      : {self.halfwidth_residual:.3e}",
            f"  centre-temperature residual : {self.temperature_residual:.3f} K",
            f"  buoyancy force target/model : "
            f"{self.target_buoyancy_force:.6g}/"
            f"{self.projected_buoyancy_force:.6g} N/m",
            f"  buoyancy force ratio        : {self.buoyancy_force_ratio:.6g}",
            f"  energy-quadrature residual  : "
            f"{self.energy_quadrature_residual:.3e}",
            f"  energy-quadrature points    : {self.energy_quadrature_points}",
            f"  Houf width mapping         : {self.houf_width_mapping}",
            f"  ground interaction         : {self.ground_interaction}",
            f"  velocity/scalar spreading : "
            f"{self.velocity_spreading_ratio:.6g}",
            f"  interface screen            : "
            f"{'pass' if self.accepted else 'fail'}",
        ]
        if self.failure_reasons:
            lines.append("  failure reasons:")
            lines += [f"    - {reason}" for reason in self.failure_reasons]
        return "\n".join(lines)

    def run(
        self,
        *,
        maximum_distance: float,
        maximum_step: float = 0.05,
        relative_tolerance: float = 1.0e-5,
    ):
        """Integrate only an interface that passed every frozen gate."""
        if not self.accepted:
            raise RuntimeError(
                "independent-energy integration refused because the interface "
                "screen failed"
            )
        return self.model.solve(
            self.state.copy(),
            maximum_distance=maximum_distance,
            maximum_step=maximum_step,
            relative_tolerance=relative_tolerance,
        )


@dataclass
class LH2CoupledResearchResult:
    """End-to-end conserved near-field and crosswind research setup."""

    near_field: LH2NearFieldResearchResult
    handoff: LH2CrosswindHandoff
    local_wind: float
    source_velocity_ratio: float
    applicability_failures: list[str] = field(default_factory=list)

    @property
    def accepted(self) -> bool:
        return (
            self.near_field.conservative
            and self.handoff.accepted
            and not self.applicability_failures
        )

    def report(self) -> str:
        return (
            self.near_field.report()
            + f"\n  handoff local wind          : {self.local_wind:.3f} m/s\n"
            + f"  source velocity/wind ratio : "
            f"{self.source_velocity_ratio:.3f}\n"
            + self.handoff.report()
            + (
                "\n  coupled applicability failures:\n"
                + "\n".join(
                    f"    - {reason}" for reason in self.applicability_failures
                )
                if self.applicability_failures else ""
            )
        )

    def run(
        self, *, distmx: float, tol: float = 1.0e-4,
        smax: float = 1.0e5,
    ) -> "JetResult":
        """Run the accepted coupled crosswind trajectory."""
        if not self.accepted:
            raise RuntimeError(
                "crosswind integration refused because the coupled research "
                "path failed conservation or applicability screening"
            )
        return self.handoff.run(distmx=distmx, tol=tol, smax=smax)


def audit_lh2_independent_energy_interface(
    near_field: LH2NearFieldResearchResult,
    crosswind_model: "JetPlume",
    *,
    streamline_distance: float | None = None,
    energy_quadrature: int = 32,
    energy_transport: str = "total",
    houf_width_mapping: str = "velocity",
    ground_interaction: str = "free",
    fit_velocity_spreading: bool = False,
    downstream_thermodynamic_profile: str = "density",
) -> LH2IndependentEnergyInterface:
    """Project five near-field fluxes into the proposed seven-state section.

    Unlike :func:`handoff_lh2_near_field_to_crosswind`, density and hydrogen
    mass fraction are independent centre states, so total energy is an exact
    fifth constraint rather than a post-hoc diagnostic.  The returned object
    audits boundary feasibility before downstream integration is permitted.
    """
    from .addons.energy_crosswind import IndependentEnergyCrosswind
    from .core.jetplume import JetIntegralFluxes

    model_class = IndependentEnergyCrosswind
    if downstream_thermodynamic_profile == "enthalpy":
        from .addons.enthalpy_profile import GaussianEnthalpyCrosswind
        if fit_velocity_spreading:
            raise ValueError("enthalpy-profile screening keeps the fixed velocity spreading ratio")
        model_class = GaussianEnthalpyCrosswind
    elif downstream_thermodynamic_profile != "density":
        raise ValueError("downstream thermodynamic profile must be 'density' or 'enthalpy'")
    if not near_field.conservative:
        raise ValueError("near-field solution failed its conservation screen")
    solution = near_field.solution
    if streamline_distance is None:
        streamline_distance = float(solution.S[-1])
    streamline_distance = float(streamline_distance)
    velocity, width, density, fraction, theta, x, height = (
        solution.state_at_s(streamline_distance)
    )
    target_mass = float(np.interp(
        streamline_distance, solution.S, solution.mass_flux
    ))
    target_hydrogen = float(np.interp(
        streamline_distance, solution.S, solution.species_flux
    ))
    target_momentum = float(np.interp(
        streamline_distance, solution.S, solution.momentum_flux
    ))
    if energy_transport == "total":
        target_energy = float(np.interp(
            streamline_distance, solution.S, solution.energy_flux
        ))
    elif energy_transport == "enthalpy":
        target_energy = near_field.model.enthalpy_flux(np.array([
            velocity, width, density, fraction, theta,
        ]))
    else:
        raise ValueError("energy transport must be 'total' or 'enthalpy'")
    target_temperature = float(np.interp(
        streamline_distance, solution.S, solution.temperature
    ))
    target = JetIntegralFluxes(
        total_mass=target_mass,
        contaminant_mass=target_hydrogen,
        momentum_x=target_momentum * math.cos(theta),
        momentum_z=target_momentum * math.sin(theta),
        energy=target_energy,
    )
    target_fluxes = {
        "total_mass": target.total_mass,
        "hydrogen": target.contaminant_mass,
        "momentum_x": target.momentum_x,
        "momentum_z": target.momentum_z,
        "energy": target.energy,
    }

    model = model_class(
        crosswind_model, near_field.model,
        quadrature_points=energy_quadrature,
        energy_transport=energy_transport,
        houf_width_mapping=houf_width_mapping,
        ground_interaction=ground_interaction,
    )
    sysz = max(
        near_field.model.spreading_ratio**2 * width**2 / 2.0, 1.0e-14
    )
    wind_trial = np.array([
        density, fraction, sysz, theta, max(velocity, 1.0e-8), x, height,
    ])
    local_wind = model._wind(wind_trial)
    initial_state = wind_trial.copy()
    initial_state[4] = max(
        velocity - local_wind * math.cos(theta), 1.0e-8
    )
    target_halfwidth = (
        near_field.model.spreading_ratio * width * math.sqrt(math.log(2.0))
    )

    def project_at_spreading(spreading_ratio: float):
        model.velocity_shape_exponent = float(spreading_ratio) ** 2
        candidate = model.project(target, initial_state)
        candidate_halfwidth = math.sqrt(
            2.0 * math.log(2.0) * candidate.state[2]
        )
        signed_width = candidate_halfwidth / target_halfwidth - 1.0
        return candidate, candidate_halfwidth, signed_width

    default_spreading = near_field.model.spreading_ratio
    projection, projected_halfwidth, signed_width = project_at_spreading(
        default_spreading
    )
    if fit_velocity_spreading and abs(signed_width) > 0.05:
        evaluated = [(default_spreading, projection, projected_halfwidth, signed_width)]
        # Avoid evaluating the exact 1.0/1.5 limits.  Near those limits the
        # five-flux projection can become almost singular and make an otherwise
        # small source audit disproportionately expensive.  A directional
        # probe finds the improving side, followed by at most two bracketed
        # secant estimates.
        probe = 1.30
        probe_result = project_at_spreading(probe)
        evaluated.append((probe, *probe_result))
        if abs(probe_result[2]) < abs(signed_width):
            directional = 1.45
        else:
            directional = 1.08
        evaluated.append((directional, *project_at_spreading(directional)))

        finite = sorted(
            (item for item in evaluated if math.isfinite(item[3])),
            key=lambda item: item[0],
        )
        brackets = [
            (left, right)
            for left, right in zip(finite[:-1], finite[1:])
            if left[3] * right[3] < 0.0
        ]
        if brackets:
            low_item, high_item = min(
                brackets,
                key=lambda pair: abs(pair[0][3]) + abs(pair[1][3]),
            )
            for _ in range(2):
                low, low_residual = low_item[0], low_item[3]
                high, high_residual = high_item[0], high_item[3]
                root = low - low_residual * (high - low) / (
                    high_residual - low_residual
                )
                root_result = project_at_spreading(root)
                root_item = (root, *root_result)
                evaluated.append(root_item)
                if abs(root_result[2]) <= 0.005:
                    break
                if low_residual * root_result[2] <= 0.0:
                    high_item = root_item
                else:
                    low_item = root_item
        spreading, projection, projected_halfwidth, signed_width = min(
            (
                item for item in evaluated
                if item[1].success and math.isfinite(item[3])
            ),
            key=lambda item: abs(item[3]),
            default=evaluated[0],
        )
        model.velocity_shape_exponent = float(spreading) ** 2
    velocity_spreading_ratio = math.sqrt(model.velocity_shape_exponent)
    # The accepted finer order is part of the interface audit and is retained
    # by the downstream energy equation rather than silently falling back to
    # the initial coarse projection order.
    model.quadrature_points = projection.quadrature_points
    state = projection.state
    projected_fluxes = {
        "total_mass": projection.fluxes.total_mass,
        "hydrogen": projection.fluxes.contaminant_mass,
        "momentum_x": projection.fluxes.momentum_x,
        "momentum_z": projection.fluxes.momentum_z,
        "energy": projection.fluxes.energy,
    }
    halfwidth_residual = abs(signed_width)
    projected_temperature = model.centre_temperature(state)
    temperature_residual = abs(projected_temperature - target_temperature)
    target_buoyancy_force = near_field.model.buoyancy_force(np.array([
        velocity, width, density, fraction, theta,
    ]))
    projected_buoyancy_force = model.buoyancy_force(state)
    buoyancy_force_ratio = (
        projected_buoyancy_force / target_buoyancy_force
        if abs(target_buoyancy_force) > 1.0e-15 else math.nan
    )

    reasons = []
    if not projection.success:
        reasons.append(projection.message)
    for name, residual in projection.relative_residuals.items():
        if residual > 1.0e-8:
            reasons.append(f"{name} balance exceeds 1e-8")
    if projection.quadrature_residual > 1.0e-5:
        reasons.append("energy quadrature changes by more than 1e-5")
    if halfwidth_residual > 0.05:
        reasons.append("hydrogen half-width mismatch exceeds 5%")
    if temperature_residual > 2.0:
        reasons.append("centre-temperature mismatch exceeds 2 K")

    return LH2IndependentEnergyInterface(
        model=model,
        state=state,
        streamline_distance=streamline_distance,
        target_fluxes=target_fluxes,
        projected_fluxes=projected_fluxes,
        relative_residuals=projection.relative_residuals,
        target_halfwidth=target_halfwidth,
        projected_halfwidth=projected_halfwidth,
        halfwidth_residual=halfwidth_residual,
        target_temperature=target_temperature,
        projected_temperature=projected_temperature,
        temperature_residual=temperature_residual,
        target_buoyancy_force=target_buoyancy_force,
        projected_buoyancy_force=projected_buoyancy_force,
        buoyancy_force_ratio=buoyancy_force_ratio,
        energy_quadrature_residual=projection.quadrature_residual,
        energy_quadrature_points=projection.quadrature_points,
        energy_transport=energy_transport,
        houf_width_mapping=houf_width_mapping,
        ground_interaction=ground_interaction,
        velocity_spreading_ratio=velocity_spreading_ratio,
        accepted=not reasons,
        failure_reasons=reasons,
    )


def handoff_lh2_near_field_to_crosswind(
    near_field: LH2NearFieldResearchResult,
    crosswind_model: "JetPlume",
    *,
    streamline_distance: float | None = None,
    thermodynamic_closure: str = "phase_profile",
    energy_quadrature: int = 32,
    crosswind_entrainment: str = "source_momentum",
) -> LH2CrosswindHandoff:
    """Project a conserved near-field station into JETPLU.

    Total mass, hydrogen mass and both momentum components are solved as four
    exact integral constraints. The default ``phase_profile`` closure carries
    the accepted near-field N2/O2/H2O phase and component-enthalpy profile
    into the crosswind table. ``phase_manifold`` additionally carries the
    actually visited upstream centre states when a large mixed source needs a
    wider concentration domain. ``centre_state`` and ``flux_energy`` reproduce
    the first and second candidates retained for falsification. The fixed
    criteria are recorded in
    ``docs/prereg-conservative-nearfield-crosswind-handoff.md`` and its
    flux-energy addendum.

    On success the supplied JETPLU model receives the mixed-source dilution
    table and the continuous Houf--Schefer entrainment constants. On failure
    its original thermodynamic and coefficient state is restored, and
    :meth:`LH2CrosswindHandoff.run` refuses downstream integration.
    """
    from scipy.optimize import brentq, least_squares

    from .core.constants import ATM_TO_PA
    from .core.jetplume import (
        J_CC, J_SYSZ, J_THETA, J_UC, J_X, J_Z,
    )

    if not near_field.conservative:
        raise ValueError("near-field solution failed its conservation screen")
    closures = {
        "centre_state", "flux_energy", "phase_profile", "phase_manifold",
    }
    if thermodynamic_closure not in closures:
        raise ValueError(
            "thermodynamic closure must be 'centre_state', 'flux_energy' "
            "'phase_profile' or 'phase_manifold'"
        )
    if crosswind_entrainment not in {"source_momentum", "local_shear"}:
        raise ValueError(
            "crosswind entrainment must be 'source_momentum' or "
            "'local_shear'"
        )
    solution = near_field.solution
    if streamline_distance is None:
        streamline_distance = float(solution.S[-1])
    streamline_distance = float(streamline_distance)
    axis_state = solution.state_at_s(streamline_distance)
    velocity, width, density, fraction, theta, x, height = axis_state
    target_temperature = float(np.interp(
        streamline_distance, solution.S, solution.temperature
    ))
    target_mass = float(np.interp(
        streamline_distance, solution.S, solution.mass_flux
    ))
    target_hydrogen = float(np.interp(
        streamline_distance, solution.S, solution.species_flux
    ))
    momentum = float(np.interp(
        streamline_distance, solution.S, solution.momentum_flux
    ))
    target_energy = float(np.interp(
        streamline_distance, solution.S, solution.energy_flux
    ))
    target_values = np.array([
        target_mass,
        target_hydrogen,
        momentum * math.cos(theta),
        momentum * math.sin(theta),
    ])
    target_fluxes = {
        "total_mass": target_values[0],
        "hydrogen": target_values[1],
        "momentum_x": target_values[2],
        "momentum_z": target_values[3],
        "energy": target_energy,
    }

    ambient = crosswind_model.th.ambient
    model = near_field.model
    if not math.isclose(ambient.tamb, model.ambient_temperature, abs_tol=0.5):
        raise ValueError("near-field and crosswind ambient temperatures differ")
    if not math.isclose(
        ambient.pamb * ATM_TO_PA, model.ambient_pressure, rel_tol=0.01
    ):
        raise ValueError("near-field and crosswind ambient pressures differ")
    if not math.isclose(
        ambient.humid, model.ambient_absolute_humidity,
        # JETPLU retains DEGADIS's rounded 28.97/18.016 kg/kmol while the
        # near field uses CoolProp's current molar masses. The same RH input
        # therefore differs by 1.16e-4 at 277 K; larger differences still
        # indicate genuinely inconsistent ambient conditions.
        rel_tol=2.0e-4, abs_tol=1.0e-10,
    ):
        raise ValueError("near-field and crosswind ambient humidities differ")

    # Preserve the caller's model if any acceptance condition fails.
    original = {
        "table": crosswind_model.th.table,
        "gas_temp": crosswind_model.th.gas.temp,
        "rhoa": crosswind_model.rhoa,
        "sc": crosswind_model.k.sc,
        "momentum_beta": crosswind_model.k.momentum_entrainment_beta,
        "density_scaled": crosswind_model.k.density_scaled_entrainment,
        "houf": crosswind_model.k.houf_buoyant_entrainment,
        "rk": crosswind_model.rk,
        "source_entrainment": crosswind_model._source_momentum_entrainment,
        "houf_alpha": crosswind_model._houf_alpha_buoy,
        "houf_froude": crosswind_model.houf_source_froude,
    }

    spreading = float(model.spreading_ratio)
    crosswind_model.k.sc = spreading**2
    crosswind_model.rk = crosswind_model.k.profile_constants()
    if crosswind_entrainment == "source_momentum":
        crosswind_model.k.momentum_entrainment_beta = (
            model.momentum_entrainment_beta
        )
        crosswind_model.k.density_scaled_entrainment = False
        crosswind_model._source_momentum_entrainment = (
            model._momentum_entrainment
        )
    else:
        # The handoff marks the end of the source-dominated axisymmetric
        # region.  Resume JETPLU's local perimeter/shear entrainment with the
        # Ricou--Spalding density scaling used by the independently validated
        # atmospheric model, while retaining the same Houf buoyancy term.
        crosswind_model.k.momentum_entrainment_beta = 0.0
        crosswind_model.k.density_scaled_entrainment = True
        crosswind_model._source_momentum_entrainment = None
    crosswind_model.k.houf_buoyant_entrainment = True
    crosswind_model._houf_alpha_buoy = model._buoyancy_coefficient
    crosswind_model.houf_source_froude = model._source_froude

    dry_air_fraction = (1.0 - fraction) / (
        1.0 + model.ambient_absolute_humidity
    )
    state_template = np.zeros(6)
    state_template[J_X] = max(float(x), 1.0e-30)
    state_template[J_Z] = float(height)
    sysz_guess = max(spreading**2 * float(width)**2 / 2.0, 1.0e-12)
    sya = (
        crosswind_model.deltay
        * state_template[J_X] ** crosswind_model.betay
    )
    sza = (
        crosswind_model.deltaz
        * state_template[J_X] ** crosswind_model.betaz
        * math.exp(
            crosswind_model.gammaz * math.log(state_template[J_X]) ** 2
        )
    )
    _sy, sz = crosswind_model._split(sysz_guess, sya, sza)
    local_wind = crosswind_model._wind(float(height), sz, math.cos(theta))
    uc_guess = max(float(velocity) - local_wind * math.cos(theta), 1.0e-8)
    scales = np.array([
        max(abs(target_values[0]), 1.0e-12),
        max(abs(target_values[1]), 1.0e-12),
        max(abs(target_values[2]), 1.0),
        max(abs(target_values[3]), 1.0),
    ])
    upper_velocity = max(10.0 * (float(velocity) + local_wind + 1.0), 10.0)
    next_initial = np.array([
        math.log(max(float(density * fraction), 1.0e-12)),
        math.log(sysz_guess), float(theta), math.log(uc_guess),
    ])

    def transferred_phase_profile_table():
        """JETPLU table made from the near-field equilibrium radial state."""
        from .core.thermo import AdiabaticTable

        eta = np.linspace(0.0, model.radial_limit, 241)
        (
            _velocity_profile, density_profile, fraction_profile,
            molecular_weight, temperature_profile, rho_h_profile,
        ) = model._profiles_at_eta(axis_state, eta)
        concentration = density_profile * fraction_profile
        enthalpy = rho_h_profile / density_profile - model._ambient_enthalpy
        yc = molecular_weight / model.fuel_molecular_weight * fraction_profile
        if thermodynamic_closure == "phase_manifold":
            centre_concentration = float(concentration[0])
            centreline_concentration = (
                solution.density * solution.mass_fraction
            )
            upstream = np.flatnonzero(
                (solution.S <= streamline_distance + 1.0e-12)
                & (
                    centreline_concentration
                    > centre_concentration * (1.0 + 1.0e-10)
                )
            )
            if len(upstream):
                extra_concentration = []
                extra_density = []
                extra_yc = []
                extra_enthalpy = []
                extra_temperature = []
                for index in upstream:
                    upstream_state = np.array([
                        solution.velocity[index], solution.width[index],
                        solution.density[index], solution.mass_fraction[index],
                        solution.theta[index], solution.x[index],
                        solution.y[index],
                    ])
                    (
                        _upstream_velocity, upstream_density,
                        upstream_fraction, upstream_mw, upstream_temperature,
                        upstream_rho_h,
                    ) = model._profiles_at_eta(
                        upstream_state, np.array([0.0])
                    )
                    extra_concentration.append(
                        upstream_density[0] * upstream_fraction[0]
                    )
                    extra_density.append(upstream_density[0])
                    extra_yc.append(
                        upstream_mw[0] / model.fuel_molecular_weight
                        * upstream_fraction[0]
                    )
                    extra_enthalpy.append(
                        upstream_rho_h[0] / upstream_density[0]
                        - model._ambient_enthalpy
                    )
                    extra_temperature.append(upstream_temperature[0])
                concentration = np.concatenate((
                    concentration, np.asarray(extra_concentration),
                ))
                density_profile = np.concatenate((
                    density_profile, np.asarray(extra_density),
                ))
                yc = np.concatenate((yc, np.asarray(extra_yc)))
                enthalpy = np.concatenate((
                    enthalpy, np.asarray(extra_enthalpy),
                ))
                temperature_profile = np.concatenate((
                    temperature_profile, np.asarray(extra_temperature),
                ))
        order = np.argsort(concentration)
        concentration = concentration[order]
        density_profile = density_profile[order]
        yc = yc[order]
        enthalpy = enthalpy[order]
        temperature_profile = temperature_profile[order]
        threshold = max(float(concentration[-1]) * 1.0e-12, 1.0e-18)
        keep = np.concatenate((
            np.array([True]), np.diff(concentration) > threshold,
        ))
        concentration = concentration[keep]
        density_profile = density_profile[keep]
        yc = yc[keep]
        enthalpy = enthalpy[keep]
        temperature_profile = temperature_profile[keep]
        if concentration[0] > threshold:
            concentration = np.insert(concentration, 0, 0.0)
            density_profile = np.insert(
                density_profile, 0, model.ambient_density
            )
            yc = np.insert(yc, 0, 0.0)
            enthalpy = np.insert(enthalpy, 0, 0.0)
            temperature_profile = np.insert(
                temperature_profile, 0, model.ambient_temperature
            )
        else:
            concentration[0] = 0.0
            density_profile[0] = model.ambient_density
            yc[0] = 0.0
            enthalpy[0] = 0.0
            temperature_profile[0] = model.ambient_temperature
        return AdiabaticTable(
            yc=np.asarray(yc), cc=np.asarray(concentration),
            rho=np.asarray(density_profile), h=np.asarray(enthalpy),
            t=np.asarray(temperature_profile),
            humid=model.ambient_absolute_humidity,
            humsrc=0.0,
            gasmw=model.fuel_molecular_weight * 1000.0,
        )

    def project_at_temperature(endpoint_temperature: float):
        """Build one virtual mixing line and solve its four native fluxes."""
        nonlocal next_initial
        if thermodynamic_closure in {"phase_profile", "phase_manifold"}:
            crosswind_model.th.table = transferred_phase_profile_table()
        else:
            crosswind_model.th.gas.temp = float(endpoint_temperature)
            handoff_enthalpy = crosswind_model.th.enthalpy(
                float(fraction), float(dry_air_fraction),
                float(endpoint_temperature),
            )
            crosswind_model.th.table = crosswind_model.th.build_adiabatic_table(
                float(fraction), float(dry_air_fraction), handoff_enthalpy,
                exact_grid=True,
            )
        crosswind_model.rhoa = crosswind_model.th.table.rhoa
        cc_limit_local = float(crosswind_model.th.table.cc[-1])

        def decode(parameters: np.ndarray) -> np.ndarray:
            candidate = state_template.copy()
            candidate[J_CC] = math.exp(parameters[0])
            candidate[J_SYSZ] = math.exp(parameters[1])
            candidate[J_THETA] = parameters[2]
            candidate[J_UC] = math.exp(parameters[3])
            return candidate

        def residual(parameters: np.ndarray) -> np.ndarray:
            represented_flux = crosswind_model.integral_fluxes(
                decode(parameters), include_energy=False
            )
            represented = np.array([
                represented_flux.total_mass,
                represented_flux.contaminant_mass,
                represented_flux.momentum_x,
                represented_flux.momentum_z,
            ])
            return (represented - target_values) / scales

        lower = np.array([
            math.log(1.0e-14), math.log(1.0e-14),
            -math.pi / 2.0 + 1.0e-10, math.log(1.0e-10),
        ])
        upper = np.array([
            math.log(max(cc_limit_local, 1.0e-14)),
            math.log(max(100.0 * sysz_guess, 10.0)),
            math.pi / 2.0 - 1.0e-10, math.log(upper_velocity),
        ])
        initial = np.minimum(np.maximum(next_initial, lower + 1.0e-12), upper - 1.0e-12)
        fit_local = least_squares(
            residual, initial, bounds=(lower, upper),
            xtol=1.0e-13, ftol=1.0e-13, gtol=1.0e-13,
            max_nfev=1000,
        )
        next_initial = fit_local.x
        state_local = decode(fit_local.x)
        flux_local = crosswind_model.integral_fluxes(
            state_local, include_energy=True,
            energy_quadrature=energy_quadrature,
        )
        return fit_local, state_local, flux_local, cc_limit_local

    virtual_source_temperature = (
        math.nan
        if thermodynamic_closure in {"phase_profile", "phase_manifold"}
        else target_temperature
    )
    fit, state, flux, cc_limit = project_at_temperature(
        virtual_source_temperature
    )
    energy_closure_failure = None
    if thermodynamic_closure == "flux_energy":
        lower_temperature = max(14.1, target_temperature - 40.0)
        upper_temperature = min(ambient.tamb, target_temperature + 40.0)
        lower_trial = project_at_temperature(lower_temperature)
        lower_energy = lower_trial[2].energy - target_energy
        upper_trial = project_at_temperature(upper_temperature)
        upper_energy = upper_trial[2].energy - target_energy
        if lower_energy == 0.0:
            virtual_source_temperature = lower_temperature
        elif upper_energy == 0.0:
            virtual_source_temperature = upper_temperature
        elif lower_energy * upper_energy > 0.0:
            energy_closure_failure = (
                "energy target is not bracketed within the frozen "
                "virtual-temperature interval"
            )
        else:
            def energy_residual(endpoint_temperature: float) -> float:
                return (
                    project_at_temperature(endpoint_temperature)[2].energy
                    - target_energy
                )

            virtual_source_temperature = float(brentq(
                energy_residual, lower_temperature, upper_temperature,
                xtol=1.0e-8, rtol=1.0e-12,
            ))
        fit, state, flux, cc_limit = project_at_temperature(
            virtual_source_temperature
        )
    energy_quadrature_residual = 0.0
    energy_quadrature_points = 0
    if thermodynamic_closure in {"phase_profile", "phase_manifold"}:
        energy_quadrature_points = energy_quadrature
        # Sharp cryogenic phase profiles can require more than the original
        # 32/64-point pair.  Increase quadrature resolution until the already
        # frozen 1e-5 convergence screen is actually met; never accept a
        # coarse integral by relaxing that screen.
        for _iteration in range(4):
            next_points = 2 * energy_quadrature_points
            fine_flux = crosswind_model.integral_fluxes(
                state, include_energy=True,
                energy_quadrature=next_points,
            )
            energy_quadrature_residual = abs(
                fine_flux.energy - flux.energy
            ) / max(abs(fine_flux.energy), 1.0)
            flux = fine_flux
            energy_quadrature_points = next_points
            if energy_quadrature_residual <= 1.0e-5:
                break
    projected_fluxes = {
        "total_mass": flux.total_mass,
        "hydrogen": flux.contaminant_mass,
        "momentum_x": flux.momentum_x,
        "momentum_z": flux.momentum_z,
        "energy": flux.energy,
    }
    relative_residuals = {
        name: abs(projected_fluxes[name] - target_fluxes[name]) / scale
        for name, scale in zip(
            ("total_mass", "hydrogen", "momentum_x", "momentum_z"),
            scales,
        )
    }
    relative_residuals["energy"] = (
        abs(flux.energy - target_energy) / max(abs(target_energy), 1.0)
    )
    target_halfwidth = spreading * float(width) * math.sqrt(math.log(2.0))
    projected_halfwidth = math.sqrt(
        2.0 * math.log(2.0) * state[J_SYSZ]
    )
    halfwidth_residual = abs(
        projected_halfwidth / target_halfwidth - 1.0
    )
    projected_temperature = float(
        crosswind_model.th.table.from_concentration(state[J_CC]).temp
    )
    temperature_residual = abs(projected_temperature - target_temperature)

    reasons = []
    for name in ("total_mass", "hydrogen", "momentum_x", "momentum_z"):
        if relative_residuals[name] > 1.0e-8:
            reasons.append(f"{name} balance exceeds 1e-8")
    energy_limit = 1.0e-8 if thermodynamic_closure == "flux_energy" else 0.02
    if relative_residuals["energy"] > energy_limit:
        reasons.append(f"energy mismatch exceeds {energy_limit:g}")
    if energy_closure_failure is not None:
        reasons.append(energy_closure_failure)
    if energy_quadrature_residual > 1.0e-5:
        reasons.append("energy quadrature changes by more than 1e-5")
    if halfwidth_residual > 0.05:
        reasons.append("hydrogen half-width mismatch exceeds 5%")
    if temperature_residual > 2.0:
        reasons.append("centre-temperature mismatch exceeds 2 K")
    if not fit.success:
        reasons.append(f"four-flux projection did not converge: {fit.message}")
    if (
        not np.all(np.isfinite(state))
        or state[J_CC] <= 0.0 or state[J_CC] > cc_limit * (1.0 + 1.0e-12)
        or state[J_SYSZ] <= 0.0 or state[J_UC] <= 0.0
    ):
        reasons.append("projected state is non-physical or outside its table")
    accepted = not reasons

    result = LH2CrosswindHandoff(
        model=crosswind_model,
        state=state,
        streamline_distance=streamline_distance,
        target_fluxes=target_fluxes,
        projected_fluxes=projected_fluxes,
        relative_residuals=relative_residuals,
        target_halfwidth=target_halfwidth,
        projected_halfwidth=projected_halfwidth,
        halfwidth_residual=halfwidth_residual,
        target_temperature=target_temperature,
        projected_temperature=projected_temperature,
        temperature_residual=temperature_residual,
        thermodynamic_closure=thermodynamic_closure,
        virtual_source_temperature=virtual_source_temperature,
        energy_quadrature_residual=energy_quadrature_residual,
        energy_quadrature_points=energy_quadrature_points,
        crosswind_entrainment=crosswind_entrainment,
        accepted=accepted,
        failure_reasons=reasons,
    )
    if not accepted:
        crosswind_model.th.table = original["table"]
        crosswind_model.th.gas.temp = original["gas_temp"]
        crosswind_model.rhoa = original["rhoa"]
        crosswind_model.k.sc = original["sc"]
        crosswind_model.k.momentum_entrainment_beta = original["momentum_beta"]
        crosswind_model.k.density_scaled_entrainment = original[
            "density_scaled"
        ]
        crosswind_model.k.houf_buoyant_entrainment = original["houf"]
        crosswind_model.rk = original["rk"]
        crosswind_model._source_momentum_entrainment = original[
            "source_entrainment"
        ]
        crosswind_model._houf_alpha_buoy = original["houf_alpha"]
        crosswind_model.houf_source_froude = original["houf_froude"]
    return result


def run_lh2_crosswind_research(
    source: "AxisymmetricJetSource | LH2ExpandedSource",
    *,
    wind: float,
    height: float | None = None,
    ambient_temperature: float = 295.0,
    ambient_pressure: float = 101325.0,
    relative_humidity: float = 0.0,
    roughness: float = 0.001,
    stability: str = "D",
    averaging: float = 60.0,
    wind_reference_height: float = 10.0,
    maximum_nearfield_distance: float = 0.08,
    handoff_distance: float | None = None,
    handoff_search_start: float | None = None,
    handoff_search_step: float | None = None,
    minimum_mass_fraction: float = 7.0e-4,
    radial_points: int = 81,
    maximum_step: float = 0.00025,
    relative_tolerance: float = 5.0e-8,
    hydrogen_spin_isomer: str = "normal",
    equilibrium_dry_air_condensation: bool = True,
    thermodynamic_closure: str = "phase_profile",
    nearfield_establishment: str = "entrained_mass",
    nearfield_energy_transport: str = "total",
    crosswind_entrainment: str = "source_momentum",
    consistent_phase_ambient: bool = False,
) -> LH2CoupledResearchResult:
    """Build the conserved near-field-to-crosswind LH2 research path.

    The local logarithmic-profile wind used by JETPLU is first calculated at
    the release height and then supplied as the near-field co-flow. This
    avoids creating momentum or kinetic energy by changing ambient velocity
    at the handoff. Only a horizontal release aligned with the mean wind is
    currently accepted; a cross-axis wind needs a genuinely three-dimensional
    near-field model rather than an axisymmetric approximation.
    """
    from .core.jetplume import J_UC
    from .validation.nearfield import hydrogen_gas_jet

    if handoff_search_start is not None and handoff_distance is not None:
        raise ValueError(
            "set either a fixed handoff distance or a handoff search, not both"
        )
    if handoff_search_step is not None and handoff_search_step <= 0.0:
        raise ValueError("handoff search step must be positive")

    expanded = source if isinstance(source, LH2ExpandedSource) else None
    source_plane = expanded.source if expanded is not None else source
    if not math.isclose(source_plane.theta, 0.0, abs_tol=1.0e-10):
        raise ValueError(
            "the coupled research path currently requires a horizontal "
            "release aligned with the mean wind"
        )
    release_height = source_plane.y if height is None else float(height)
    if release_height <= 0.0:
        raise ValueError("crosswind release height must be positive")
    if not math.isclose(source_plane.y, release_height, abs_tol=1.0e-12):
        source_plane = replace(source_plane, y=release_height)
        source = (
            LH2ExpandedSource(source=source_plane, expansion=expanded.expansion)
            if expanded is not None else source_plane
        )
    crosswind, initial_state = hydrogen_gas_jet(
        rate=source_plane.fuel_mass_flow,
        diameter=source_plane.diameter,
        velocity=source_plane.velocity,
        wind=wind,
        height=release_height,
        source_temperature=source_plane.temperature,
        source_density=source_plane.density,
        source_mass_fraction=source_plane.mass_fraction,
        theta=source_plane.theta,
        relative_humidity=relative_humidity,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        roughness=roughness,
        stability=stability,
        averaging=averaging,
        wind_reference_height=wind_reference_height,
        sc=1.16**2,
        density_scaled_entrainment=False,
        momentum_entrainment_beta=0.28,
        houf_entrainment=True,
    )
    local_wind = source_plane.velocity - float(initial_state[J_UC])
    source_velocity_ratio = (
        source_plane.velocity / local_wind if local_wind > 0.0 else math.inf
    )
    near_field = run_lh2_near_field_research(
        source,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        relative_humidity=relative_humidity,
        ambient_coflow_velocity=local_wind,
        maximum_distance=maximum_nearfield_distance,
        minimum_mass_fraction=minimum_mass_fraction,
        radial_points=radial_points,
        maximum_step=maximum_step,
        relative_tolerance=relative_tolerance,
        hydrogen_spin_isomer=hydrogen_spin_isomer,
        equilibrium_dry_air_condensation=(
            equilibrium_dry_air_condensation
        ),
        establishment=nearfield_establishment,
        energy_transport=nearfield_energy_transport,
        consistent_phase_ambient=consistent_phase_ambient,
    )
    # Both branches received the same RH, but JETPLU's cached saturation
    # table and rounded legacy molecular weights otherwise differ from the
    # phase model by about 0.1%. The handoff replaces the original dilution
    # table, so carry the near-field value forward exactly.
    crosswind.th.ambient.humid = near_field.model.ambient_absolute_humidity
    if handoff_search_start is not None:
        search_start = float(handoff_search_start)
        search_step = (
            source_plane.diameter
            if handoff_search_step is None else float(handoff_search_step)
        )
        if search_step <= 0.0:
            raise ValueError("handoff search step must be positive")
        search_stop = float(near_field.solution.S[-1])
        if search_start > search_stop:
            raise ValueError("handoff search starts beyond the near-field run")
        stations = search_start + search_step * np.arange(
            int(math.floor((search_stop - search_start) / search_step)) + 1
        )
        handoff = None
        for station in stations:
            handoff = handoff_lh2_near_field_to_crosswind(
                near_field,
                crosswind,
                streamline_distance=float(station),
                thermodynamic_closure=thermodynamic_closure,
                crosswind_entrainment=crosswind_entrainment,
            )
            if handoff.accepted:
                break
        assert handoff is not None
    else:
        handoff = handoff_lh2_near_field_to_crosswind(
            near_field,
            crosswind,
            streamline_distance=handoff_distance,
            thermodynamic_closure=thermodynamic_closure,
            crosswind_entrainment=crosswind_entrainment,
        )
    return LH2CoupledResearchResult(
        near_field=near_field,
        handoff=handoff,
        local_wind=local_wind,
        source_velocity_ratio=source_velocity_ratio,
        applicability_failures=(
            [
                f"source velocity/wind ratio {source_velocity_ratio:.3f} is "
                f"below the momentum-dominated limit {MOMENTUM_RATIO:g}"
            ]
            if source_velocity_ratio < MOMENTUM_RATIO else []
        ),
    )


def lh2_source_from_measured_throat(
    *,
    throat_diameter: float,
    throat_pressure: float,
    throat_temperature: float,
    throat_density: float,
    throat_velocity: float,
    mass_flow: float | None = None,
    ambient_pressure: float = 101325.0,
    theta: float = math.pi / 2.0,
    x: float = 0.0,
    y: float = 0.0,
) -> LH2ExpandedSource:
    """Build the ambient-pressure source used by the Raman validation.

    All pressures are absolute Pa and SI units are used throughout. If
    ``mass_flow`` is omitted it is calculated from the reported throat
    density, velocity and diameter. The HyRAM+ measured-throat expansion then
    conserves mass, axial momentum including pressure thrust, and total
    specific energy.
    """
    from .addons.axisymmetric_jet import AxisymmetricJetSource
    from .addons.notional import expand_measured_throat_to_ambient

    throat_area = math.pi * throat_diameter**2 / 4.0
    if mass_flow is None:
        mass_flow = throat_density * throat_velocity * throat_area
    expansion = expand_measured_throat_to_ambient(
        fluid="Hydrogen",
        mass_flow=mass_flow,
        throat_diameter=throat_diameter,
        throat_pressure=throat_pressure,
        throat_temperature=throat_temperature,
        throat_density=throat_density,
        throat_velocity=throat_velocity,
        ambient_pressure=ambient_pressure,
    )
    source = AxisymmetricJetSource(
        diameter=expansion.diameter,
        velocity=expansion.velocity,
        density=expansion.density,
        temperature=expansion.temperature,
        theta=theta,
        x=x,
        y=y,
    )
    return LH2ExpandedSource(source=source, expansion=expansion)


def run_lh2_near_field_research(
    source: "AxisymmetricJetSource | LH2ExpandedSource",
    *,
    ambient_temperature: float = 295.0,
    ambient_pressure: float = 101325.0,
    relative_humidity: float = 0.0,
    ambient_coflow_velocity: float = 0.0,
    maximum_distance: float = 0.13,
    minimum_mass_fraction: float = 7.0e-4,
    radial_points: int = 81,
    maximum_step: float = 0.00025,
    relative_tolerance: float = 5.0e-8,
    radiative_absorptivity: float = 0.0,
    include_argon_phase: bool = False,
    equilibrium_dry_air_condensation: bool = True,
    hydrogen_spin_isomer: str = "normal",
    establishment: str = "scalar_peak",
    energy_transport: str = "total",
    consistent_phase_ambient: bool = False,
) -> LH2NearFieldResearchResult:
    """Run the recommended conserved cryogenic free-jet research model.

    The fixed physics are scalar-peak-constrained Gaussian establishment,
    equilibrium N2/O2 phase change, and component temperature-dependent
    enthalpy. The 4/4 Raman result applies to dry air at 295 K and 101325 Pa;
    humidity and co-flow remain explicit sensitivity inputs and generate an
    applicability warning. ``hydrogen_spin_isomer`` changes only the
    downstream component caloric table; ``para`` and ``ortho`` are explicit
    sensitivities because the validation release composition was not reported.
    ``consistent_phase_ambient`` is an opt-in correction that uses the same
    explicit-species ideal EOS for ambient gas as for the phase inversion.
    It is not included in the historical Raman validation claim.

    ``source`` must be an ambient-pressure gas-plane state. Use
    :func:`lh2_source_from_measured_throat` when throat measurements are the
    available boundary condition.
    """
    from CoolProp.CoolProp import PropsSI

    from .addons.axisymmetric_jet import ConservedGaussianJet, phase_ambient_from_rh

    if not 0.0 <= relative_humidity <= 100.0:
        raise ValueError("relative humidity must lie between 0 and 100 percent")
    spin_species = {
        "normal": "Hydrogen",
        "para": "ParaHydrogen",
        "ortho": "OrthoHydrogen",
    }
    try:
        spin_name = hydrogen_spin_isomer.lower()
        hydrogen_enthalpy_species = spin_species[spin_name]
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            "hydrogen spin isomer must be 'normal', 'para' or 'ortho'"
        ) from exc
    expanded = source if isinstance(source, LH2ExpandedSource) else None
    source_plane = expanded.source if expanded is not None else source

    # A saturated-liquid storage calculation can return a perfectly conserved
    # ambient-pressure *two-phase* H2 state. Density and temperature alone do
    # not advertise that quality, and treating it as the single gas phase of
    # this Gaussian model would be a large, silent source error.
    if math.isclose(source_plane.mass_fraction, 1.0, abs_tol=1.0e-12):
        saturation_temperature = float(PropsSI(
            "T", "P", ambient_pressure, "Q", 1, "Hydrogen"
        ))
        saturated_vapour_density = float(PropsSI(
            "D", "P", ambient_pressure, "Q", 1, "Hydrogen"
        ))
        if (
            source_plane.temperature < saturation_temperature - 0.1
            or (
                source_plane.temperature <= saturation_temperature + 0.5
                and source_plane.density > 1.1 * saturated_vapour_density
            )
        ):
            raise ValueError(
                "the near-field Gaussian model requires a single-phase gas "
                "source; the supplied ambient-pressure H2 state is liquid "
                "or two-phase"
            )

    dry_air_mw = float(PropsSI("M", "Air"))
    water_mw = float(PropsSI("M", "Water"))
    saturation_pressure = float(PropsSI(
        "P", "T", ambient_temperature, "Q", 0, "Water"
    ))
    vapour_pressure = relative_humidity / 100.0 * saturation_pressure
    if vapour_pressure >= ambient_pressure:
        raise ValueError("water vapour pressure must be below ambient pressure")
    absolute_humidity = (
        water_mw / dry_air_mw
        * vapour_pressure / (ambient_pressure - vapour_pressure)
    )
    humid_air_mw = (1.0 + absolute_humidity) / (
        1.0 / dry_air_mw + absolute_humidity / water_mw
    )
    dry_air_density = float(PropsSI(
        "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
    ))
    ambient_density = dry_air_density * humid_air_mw / dry_air_mw
    if consistent_phase_ambient:
        dry_air_mw, absolute_humidity, ambient_density = phase_ambient_from_rh(
            ambient_temperature, ambient_pressure, relative_humidity,
            include_argon=include_argon_phase,
        )

    model = ConservedGaussianJet(
        source_plane,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        ambient_density=ambient_density,
        fuel_molecular_weight=float(PropsSI("M", "Hydrogen")),
        ambient_molecular_weight=dry_air_mw,
        fuel_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure,
            "Hydrogen",
        )),
        ambient_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure, "Air"
        )),
        radial_points=radial_points,
        conservative_establishment=establishment,
        equilibrium_air_condensation=True,
        equilibrium_dry_air_condensation=(
            equilibrium_dry_air_condensation
        ),
        ambient_absolute_humidity=absolute_humidity,
        ambient_coflow_velocity=ambient_coflow_velocity,
        temperature_dependent_phase_enthalpy=True,
        hydrogen_enthalpy_species=hydrogen_enthalpy_species,
        radiative_absorptivity=radiative_absorptivity,
        equilibrium_argon_condensation=include_argon_phase,
        consistent_phase_ambient=consistent_phase_ambient,
        energy_transport=energy_transport,
    )
    solution = model.solve(
        maximum_distance=maximum_distance,
        minimum_mass_fraction=minimum_mass_fraction,
        maximum_step=maximum_step,
        relative_tolerance=relative_tolerance,
        method="LSODA",
    )
    boundary_residual = max(
        abs(value) for value in model.establishment_residuals.values()
    )
    species_drift = float(np.max(np.abs(
        solution.species_flux / solution.species_flux[0] - 1.0
    )))
    relative_energy = solution.energy_flux.copy()
    if energy_transport == "total":
        relative_energy = (
            relative_energy
            - 0.5 * ambient_coflow_velocity**2 * solution.mass_flux
        )
    energy_drift = float(np.max(np.abs(
        relative_energy / relative_energy[0] - 1.0
    )))
    if radiative_absorptivity:
        corrected_energy = relative_energy - solution.radiative_heat_added
        energy_drift = float(np.max(np.abs(
            corrected_energy / corrected_energy[0] - 1.0
        )))

    warnings = []
    if relative_humidity != 0.0:
        warnings.append(
            "humidity is a sensitivity only; the Raman experiment did not "
            "report RH and the accepted 4/4 result is the dry-air case"
        )
    if ambient_coflow_velocity != 0.0:
        warnings.append(
            "co-flow is represented conservatively but was not part of the "
            "accepted dry-air configuration"
        )
    if radiative_absorptivity != 0.0:
        warnings.append(
            "radiation is an effective-absorptivity sensitivity; no Raman "
            "optical depth was reported and the accepted model is adiabatic"
        )
    if include_argon_phase:
        warnings.append(
            "argon phase change is a conservative composition sensitivity "
            "but did not consistently reduce the four Raman errors"
        )
    if not equilibrium_dry_air_condensation:
        warnings.append(
            "dry-air condensation is suppressed as a metastable-gas "
            "research bound; it is not a validated default"
        )
    if energy_transport == "enthalpy":
        warnings.append(
            "Li equation-35 enthalpy transport is a research sensitivity; "
            "the validated Raman default transports total energy"
        )
    if spin_name != "normal":
        warnings.append(
            "hydrogen spin composition is a caloric sensitivity; the Raman "
            "release composition was not reported and the accepted baseline "
            "uses normal hydrogen"
        )
    if not math.isclose(ambient_temperature, 295.0, abs_tol=0.5):
        warnings.append("ambient temperature is outside the 295 K Raman condition")
    if not math.isclose(ambient_pressure, 101325.0, rel_tol=0.01):
        warnings.append("ambient pressure is outside the 1 atm Raman condition")
    if not math.isclose(source_plane.mass_fraction, 1.0, abs_tol=1.0e-12):
        warnings.append("source is not the pure-hydrogen boundary used in validation")
    if expanded is not None:
        throat = expanded.expansion
        checks = (
            ("throat diameter", throat.throat_diameter, 0.001, 0.00125, "m"),
            ("throat temperature", throat.throat_temperature, 37.4, 45.7, "K"),
            ("throat pressure", throat.throat_pressure, 0.972e5, 2.422e5, "Pa"),
            ("throat velocity", throat.throat_velocity, 498.2, 558.9, "m/s"),
        )
        for name, value, lower, upper, unit in checks:
            if not lower <= value <= upper:
                warnings.append(
                    f"{name} {value:g} {unit} is outside the Raman range "
                    f"{lower:g}--{upper:g} {unit}"
                )

    return LH2NearFieldResearchResult(
        source=source_plane,
        solution=solution,
        model=model,
        maximum_boundary_residual=boundary_residual,
        maximum_species_drift=species_drift,
        maximum_energy_drift=energy_drift,
        warnings=warnings,
        notes={
            "configuration": (
                "scalar_peak; "
                + (
                    "equilibrium N2/O2"
                    + ("/Ar" if include_argon_phase else "")
                    if equilibrium_dry_air_condensation
                    else "metastable gaseous N2/O2/Ar"
                )
                + "; equilibrium H2O"
                + f"; component h(T); {spin_name} H2"
                + f"; {energy_transport} energy transport"
            ),
            "validation": (
                "Hecht-Panda final 2019 journal fits, nine Table-1 dry-air "
                "Raman cases; provisional four of four printed metrics "
                "within 25% pending aggregate-fit membership/uncertainty"
            ),
        },
    )


@dataclass
class Assessment:
    """What an LH₂ release does, and how well that is known."""

    #: Downwind distance at which the centreline falls to the lower flammable
    #: limit, m.  ``nan`` if the cloud never reaches it within the range run.
    distance_to_lfl: float
    #: Distance to the stoichiometric mixture, where an ignition does most
    #: damage, m.
    distance_to_stoichiometric: float
    #: Height of the lowest flammable gas, m -- at ``at_distance`` if one was
    #: given, otherwise at the far end of the flammable envelope.
    lowest_flammable_height: float
    #: ``"grounded"`` below a metre, ``"low"`` to ten metres, ``"aloft"``
    #: above.  This is the best-established output.
    regime: str
    #: Mole fraction at which the cloud stops being denser than air.
    neutral_buoyancy: float
    #: The plume centre, ``(x, z, mole fraction)`` along the trajectory.
    trajectory: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    warnings: list[str] = field(default_factory=list)
    notes: dict[str, str] = field(default_factory=dict)

    def report(self) -> str:
        if self.notes.get("model path") == "JetPlume":
            concentration_evidence = (
                f"   [{NEAR_FIELD_CORRECTED.cite()},"
                f" grade {NEAR_FIELD_CORRECTED.grade}]"
            )
        else:
            concentration_evidence = (
                "   [no direct pool-concentration validation]"
            )
        lines = [
            "Liquid hydrogen release",
            f"  cloud regime                : {self.regime}"
            f"   [{LIFTOFF_HEIGHT.detail['regime']} on the NASA trials,"
            f" grade {LIFTOFF_HEIGHT.grade}]",
            f"  lowest flammable gas        : {self.lowest_flammable_height:.1f} m"
            f"   [{LIFTOFF_HEIGHT.cite()}, grade {LIFTOFF_HEIGHT.grade}]",
            f"  distance to 4 mol % (LFL)   : {self.distance_to_lfl:.1f} m"
            f"{concentration_evidence}",
            f"  distance to stoichiometric  : "
            f"{self.distance_to_stoichiometric:.1f} m",
            f"  buoyant below               : "
            f"{self.neutral_buoyancy * 100:.2f} mol %"
            f"   [{NEUTRAL_BUOYANCY.source}, grade {NEUTRAL_BUOYANCY.grade}]",
        ]
        if self.warnings:
            lines.append("  outside the validated range:")
            lines += [f"    - {w}" for w in self.warnings]
        return "\n".join(lines)


def _check(kind, rate, wind, distance, velocity_ratio=None, **rest) -> list[str]:
    """Warnings for a release outside what its own campaign checked.

    ``distance`` must be the furthest distance a *reported* answer depends on,
    not the integration limit.  Passing the limit -- which this did until the
    two were separated -- makes every call with the default ``max_distance``
    warn about a distance the caller never asked about, while the answers it
    did ask for sit comfortably inside the range.
    """
    out = []
    if velocity_ratio is not None and velocity_ratio < MOMENTUM_RATIO:
        out.append(
            f"exit velocity is only {velocity_ratio:.1f} times the wind; below "
            f"{MOMENTUM_RATIO:g} the wind steers the plume and the "
            f"concentrations are not defensible"
        )
    out += check_range(
        kind, rate=rate, wind=wind, distance=distance, **rest
    )
    return out


def _first_falling_crossing(points, level: float) -> float:
    """Interpolate the first downwind crossing of ``level``.

    ``points`` is an iterable of ``(x, concentration)`` pairs. A plume can
    briefly increase in concentration at a ground image or handover, so the
    crossing is deliberately searched in trajectory order rather than after
    sorting by concentration.
    """
    for (x0, c0), (x1, c1) in zip(points, points[1:]):
        if c0 >= level > c1:
            if c0 == c1:
                return float(x1)
            f = (c0 - level) / (c0 - c1)
            return float(x0 + f * (x1 - x0))
    return float("nan")


def _jet_lowest_flammable_height(trajectory, x: float) -> float:
    """Lowest height at which a jet is at or above the LFL at ``x``.

    The jet has a Gaussian vertical profile plus its ground image. Returning
    the centre minus a nominal radius, as the pool path does, would discard
    the profile the validated comparison actually uses. Search that same
    profile and refine the first LFL crossing by bisection.
    """
    state = trajectory.at(x)
    if state is None:
        return 0.0
    high = max(state.z + 4.0 * state.sz, 4.0 * state.sz, 1.0)
    heights = np.linspace(0.0, high, 401)
    concentrations = np.array([
        trajectory.concentration_at(x, 0.0, float(z)) / 100.0
        for z in heights
    ])
    indices = np.flatnonzero(concentrations >= LFL)
    if not len(indices):
        return 0.0
    first = int(indices[0])
    if first == 0:
        return 0.0
    lo, hi = float(heights[first - 1]), float(heights[first])
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if trajectory.concentration_at(x, 0.0, mid) / 100.0 >= LFL:
            hi = mid
        else:
            lo = mid
    return hi


def assess(
    *,
    rate: float,
    wind: float,
    height: float = 0.0,
    orifice: float | None = None,
    pool_diameter: float | None = None,
    storage_pressure: float = 1.013,
    ambient_temperature: float = 288.15,
    relative_humidity: float = 65.0,
    ambient_pressure: float = 101325.0,
    max_distance: float = 100.0,
    at_distance: float | None = None,
) -> Assessment:
    """Assess a liquid hydrogen release.

    Parameters
    ----------
    rate
        Contaminant mass rate, kg/s.
    wind
        Speed at the release, m/s.
    height
        Release elevation, m.  Zero for a pool.
    orifice
        Orifice diameter, m, for a pressurised release.
    pool_diameter
        Pool diameter, m, for a spill.  One of ``orifice`` or ``pool_diameter``
        is required.
    storage_pressure
        Absolute, bar.  Only used for a pressurised release.

    Notes
    -----
    A pressurised release goes through the flashing source term and the jet
    model; a pool goes straight to the buoyant plume.  Both then use the same
    hydrogen mixing line, so the concentration and buoyancy answers are
    consistent between them.
    """
    if (orifice is None) == (pool_diameter is None):
        raise ValueError("give exactly one of an orifice diameter or a pool diameter")
    if rate <= 0.0:
        raise ValueError("mass release rate must be positive")
    if wind <= 0.0:
        raise ValueError("wind speed must be positive")
    if max_distance <= 0.0:
        raise ValueError("maximum distance must be positive")
    if at_distance is not None and not (0.0 <= at_distance <= max_distance):
        raise ValueError("at_distance must lie between zero and max_distance")
    if orifice is not None:
        if orifice <= 0.0:
            raise ValueError("orifice diameter must be positive")
        if height <= 0.0:
            raise ValueError("a pressurised jet requires a positive release elevation")
    elif pool_diameter <= 0.0:
        raise ValueError("pool diameter must be positive")

    from CoolProp.CoolProp import PropsSI

    from .addons import LiftoffPlume
    from .core.atmosphere import absolute_humidity
    from .core.thermo import (
        AmbientConditions,
        CoolPropBackend,
        GasProperties,
        Thermo,
    )
    from .core.constants import PI
    from .validation.flashing import equivalent_source

    PI_4 = PI / 4.0

    backend = CoolPropBackend("Hydrogen")
    humid, _rh = absolute_humidity(
        ambient_temperature, ambient_pressure / 101325.0,
        backend.water_vapour_pressure, relhum=relative_humidity,
    )
    boil = PropsSI("T", "P", ambient_pressure, "Q", 0, "Hydrogen")
    gas = GasProperties(
        name="LH2", mw=2.016, temp=boil,
        rho=PropsSI("D", "T", boil, "Q", 1, "Hydrogen"),
        coolprop_name="Hydrogen", ulc=UFL, llc=LFL,
    )
    ambient = AmbientConditions(
        tamb=ambient_temperature, pamb=ambient_pressure / 101325.0,
        humid=humid, tsurf=ambient_temperature, ihtfl=1, iwtfl=1,
    )
    thermo = Thermo(
        gas=gas, ambient=ambient, backend=backend, legacy_numerics=False
    )
    thermo.reference_enthalpies()

    # source composition: flashed for a pressurised release, pure vapour for a
    # pool, because a pool has already flashed by the time it boils off
    source_wc = 1.0
    exit_velocity = None
    if orifice is not None and storage_pressure > 1.1:
        stored = PropsSI("T", "P", storage_pressure * 1e5, "Q", 0, "Hydrogen")
        flash = equivalent_source(
            "Hydrogen", storage_temperature=stored,
            ambient_temperature=ambient_temperature,
            ambient_pressure=ambient_pressure, molecular_weight=2.016,
        )
        source_wc = flash.mass_fraction
        area = PI_4 * orifice * orifice
        exit_velocity = rate / max(flash.orifice_density * area, 1e-30)
    thermo.ambient.humsrc = 0.0
    table = thermo.build_adiabatic_table(1.0, 0.0, thermo.hmrte)

    # where the cloud stops being dense: the only crossing on the mixing line
    a = table.as_array()
    yc, rho = a[:, 0], a[:, 2]
    neutral = 1.0
    for i in range(1, len(a)):
        if (rho[i - 1] - table.rhoa) * (rho[i] - table.rhoa) < 0.0:
            f = (rho[i - 1] - table.rhoa) / (rho[i - 1] - rho[i])
            neutral = float(yc[i - 1] + f * (yc[i] - yc[i - 1]))

    traj, d_lfl, d_stoich, lowest = [], float("nan"), float("nan"), 0.0
    model_path = "LiftoffPlume"

    if orifice is not None:
        # A pressurised release must use the same jet path that was validated
        # against PRESLHY and Spadeadam. The previous implementation computed
        # a flash and then sent the result straight through LiftoffPlume, so
        # its documented "jet" answer was actually a pool-plume calculation.
        from .validation.nearfield import STEP, Trajectory, hydrogen_jet

        jp, y0 = hydrogen_jet(
            rate=rate, diameter=orifice, wind=wind, height=height,
            ambient_temperature=ambient_temperature,
            relative_humidity=relative_humidity,
            ambient_pressure=ambient_pressure,
            storage_pressure_barg=max(
                storage_pressure - ambient_pressure / 1.0e5, 0.0
            ),
            # The thermodynamically consistent flashing source already makes
            # the high-wind Spadeadam plume remain close to the ground and
            # lets the low-wind plume lift, as observed.  The older
            # `ground_effect` option permanently held both down; it was only
            # beneficial while the source table was silently deleting the
            # air entrained by flashing.  Keep the physical free-to-detach
            # trajectory on the public assessment path.
            corrections=True, ground_effect=False,
        )
        run = jp.run(
            y0, distmx=STEP,
            smax=max(40.0, 2.5 * max(max_distance, at_distance or 0.0)),
        )
        jet = Trajectory(jp.th.table, run.rows)
        model_path = "JetPlume"
        if jet.ok:
            rows = run.rows[np.argsort(run.rows[:, 0])]
            for row in rows:
                x, z = float(row[0]), float(row[1])
                if x < 0.0 or x > max_distance:
                    continue
                centre = jet.concentration_at(x, 0.0, z) / 100.0
                traj.append([x, z, centre])

            centreline = [(r[0], r[2]) for r in traj]
            d_lfl = _first_falling_crossing(centreline, LFL)
            d_stoich = _first_falling_crossing(centreline, STOICHIOMETRIC)
            if at_distance is not None:
                lowest = _jet_lowest_flammable_height(jet, at_distance)
            else:
                flammable = [r[0] for r in traj if r[2] >= LFL]
                if flammable:
                    lowest = _jet_lowest_flammable_height(jet, flammable[-1])
    else:
        plume = LiftoffPlume(
            rho_ambient=table.rhoa, wind=wind, segment_length=pool_diameter,
            density_of=lambda c: table.from_mass_fraction(
                max(min(c, 1.0), 0.0), wa=0.0
            ).rho,
        )
        lfl_wc = table.from_mole_fraction(LFL).wc
        result = plume.run(
            rate=rate, concentration=source_wc, velocity=wind, height=height,
            max_distance=max_distance, min_concentration=lfl_wc * 0.1,
            max_height=200.0,
        )

        for s in result.states:
            y = table.from_mass_fraction(s.concentration, wa=0.0).yc
            traj.append([s.x, s.z, y])
            if math.isnan(d_lfl) and y <= LFL:
                d_lfl = s.x
            if math.isnan(d_stoich) and y <= STOICHIOMETRIC:
                d_stoich = s.x
            if y >= LFL:
                lowest = max(s.z - s.radius, 0.0)

        if at_distance is not None and traj:
            arr = np.array(traj)
            if at_distance <= arr[-1, 0]:
                i = int(np.searchsorted(arr[:, 0], at_distance))
                i = min(i, len(result.states) - 1)
                s = result.states[i]
                lowest = max(s.z - s.radius, 0.0)
            else:
                lowest = 0.0  # the flammable cloud does not reach that far

    regime = "grounded" if lowest < 1.0 else ("low" if lowest < 10.0 else "aloft")

    # The furthest distance any *reported* number leans on. `max_distance` is
    # how far the integration was allowed to run, which is a computational
    # choice and not a claim about anything; warning on it fired on every call
    # that took the default while the reported answers were all in range.
    relied_on = max(
        [d for d in (at_distance, d_lfl, d_stoich) if d is not None and not math.isnan(d)]
        or [0.0]
    )
    kind = "jet" if orifice is not None else "pool"
    return Assessment(
        distance_to_lfl=d_lfl,
        distance_to_stoichiometric=d_stoich,
        lowest_flammable_height=lowest,
        regime=regime,
        neutral_buoyancy=neutral,
        trajectory=np.array(traj) if traj else np.zeros((0, 3)),
        warnings=_check(
            kind, rate, wind, relied_on,
            velocity_ratio=(exit_velocity / wind) if exit_velocity else None,
            orifice=orifice, height=height if kind == "jet" else None,
            diameter=pool_diameter,
        ),
        notes={
            "source": (
                "flashing jet" if orifice is not None and source_wc < 1.0
                else "vapour jet" if orifice is not None
                else "pool or vapour"
            ),
            "source mass fraction": f"{source_wc:.3f}",
            "model path": model_path,
        },
    )
