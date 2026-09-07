"""Physics DEGADIS does not have.

Everything in :mod:`degali.core` is a port of DEGADIS 2.1 and is validated
against it. This package is different: it implements published theory the
original lacks, so nothing here can be checked against the Fortran and it is
kept separate for that reason.

:mod:`~degali.addons.liftoff`
    Buoyant lift-off of a ground-based plume, from the URAHFREP integral
    model, as a standalone plume.

:mod:`~degali.addons.unified`
    Dense slumping and buoyant rise in one set of equations, with no handover
    threshold between them: entrainment blends between the ground-layer law
    and the free-perimeter law by how much of the cross-section is still on
    the ground. This removes the one free parameter the staged approach needed.

:mod:`~degali.addons.buoyant`
    The same vertical momentum balance as a *closure* for the downwind model,
    so an existing DEGADIS run can be given the ability to lift its cloud off
    the ground by swapping one object. Dense behaviour is unchanged; the
    difference appears only where DEGADIS would have frozen the cloud in
    place.
"""

from .unified import UnifiedClosure
from .axisymmetric_jet import (
    AxisymmetricJetResult,
    AxisymmetricJetSource,
    ConservedGaussianJet,
    InitialEntrainmentResult,
    SourceEnthalpyBoundary,
    entrain_and_heat_initial_plug,
)
from .hydrogen_eos import HydrogenGasDeparture, hydrogen_gas_departure
from .energy_crosswind import (
    IndependentEnergyCrosswind,
    IndependentEnergyJetResult,
    IndependentEnergyProjection,
)
from .notional import (
    EnergyConservingNotionalNozzle,
    IsentropicThroat,
    MeasuredThroatExpansion,
    NotionalNozzle,
    choking_pressure,
    energy_conserving_notional_nozzle,
    expand_measured_throat_to_ambient,
    isentropic_throat,
    notional_nozzle,
)
from .lh2_droplets import (
    critical_droplet_diameter_for_phase_delay,
    critical_diffusivity_for_phase_delay,
    FlashingHydrogenDropletSource,
    gasflow_phase_relaxation_coefficient,
    HomogeneousEquilibriumHydrogenSource,
    flashing_hydrogen_droplet_source,
    homogeneous_equilibrium_hydrogen_source,
    minimum_heat_limited_hydrogen_evaporation_time,
)
from .cryogenic_air import (
    AirPhaseEquilibrium,
    CondensedAirSource,
    CondensedAirStep,
    HydrogenEvaporationEndpoint,
    MultiphaseHydrogenSourcePlane,
    LiZone3,
    air_saturation_pressure,
    equilibrium_air_phase_split,
    li2026_zone3,
    multiphase_hydrogen_evaporation_endpoint,
    multiphase_hydrogen_source_plane,
    minimum_heat_limited_sublimation_time,
    particle_relaxation_time,
    particle_terminal_velocity,
    ranz_marshall_transfer_number,
    transported_condensed_air_source,
)
from .buoyant import LAMBDA, BuoyantClosure, make_buoyant
from .liftoff import (
    ALPHA,
    BETA,
    GAMMA,
    RI_CLEAR,
    RI_ONSET,
    RI_SUBSTANTIAL,
    LiftoffPlume,
    LiftoffResult,
    LiftoffState,
    liftoff_regime,
    richardson_liftoff,
)

__all__ = [
    "UnifiedClosure",
    "AxisymmetricJetSource", "AxisymmetricJetResult", "ConservedGaussianJet",
    "InitialEntrainmentResult", "entrain_and_heat_initial_plug",
    "IndependentEnergyCrosswind", "IndependentEnergyProjection",
    "IndependentEnergyJetResult",
    "SourceEnthalpyBoundary", "HydrogenGasDeparture", "hydrogen_gas_departure",
    "NotionalNozzle", "notional_nozzle", "choking_pressure",
    "IsentropicThroat", "isentropic_throat",
    "EnergyConservingNotionalNozzle", "energy_conserving_notional_nozzle",
    "MeasuredThroatExpansion", "expand_measured_throat_to_ambient",
    "FlashingHydrogenDropletSource", "flashing_hydrogen_droplet_source",
    "gasflow_phase_relaxation_coefficient",
    "critical_droplet_diameter_for_phase_delay",
    "critical_diffusivity_for_phase_delay",
    "HomogeneousEquilibriumHydrogenSource",
    "homogeneous_equilibrium_hydrogen_source",
    "minimum_heat_limited_hydrogen_evaporation_time",
    "AirPhaseEquilibrium", "equilibrium_air_phase_split",
    "CondensedAirSource", "CondensedAirStep", "transported_condensed_air_source",
    "HydrogenEvaporationEndpoint", "multiphase_hydrogen_evaporation_endpoint",
    "MultiphaseHydrogenSourcePlane", "multiphase_hydrogen_source_plane",
    "LiZone3", "li2026_zone3", "air_saturation_pressure",
    "particle_relaxation_time", "particle_terminal_velocity",
    "ranz_marshall_transfer_number",
    "minimum_heat_limited_sublimation_time",
    "BuoyantClosure", "make_buoyant", "LAMBDA",
    "LiftoffPlume", "LiftoffResult", "LiftoffState",
    "richardson_liftoff", "liftoff_regime",
    "ALPHA", "BETA", "GAMMA", "RI_ONSET", "RI_SUBSTANTIAL", "RI_CLEAR",
]
