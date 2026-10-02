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
from .yawed_crosswind import YawedCrosswind, YawedTrajectory
from .wall_jet_transition import (
    FiniteWallJetCrosswind,
    FiniteWallJetResult,
    WallJetDiagnostics,
    WallJetTransitionConfig,
)
from .transient_receptor import (
    FixedReceptor,
    ReceptorTrace,
    SteadyPlumeTable,
    WindowStatistics,
    WindHistory,
    first_order_sensor_response,
    meteorological_from_to_math_radians,
    replay_fixed_receptors,
)
from .site_geometry import (
    AxisAlignedCuboid,
    ObstacleEncounter,
    ObstacleScreenResult,
    SegmentIntersection,
    TransverseWall,
    WindFrame,
    screen_trajectory,
)
from .obstacle_wake_validation import (
    DenseGasFenceBenchmark,
    DenseGasFenceObservation,
    NeutralObstacleWakeBenchmark,
    NeutralObstacleWakePoint,
    pair_smedis_fence_trials,
    read_aij_case_h,
    read_smedis_fence_pair,
)
from .transport_evidence import (
    PARTICLE_SLIP_BOUNDARY,
    TURBULENCE_BOUNDARY,
    TransportEvidenceBoundary,
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
from .cryogenic_blowdown import (
    CryogenicBlowdownConfig,
    CryogenicBlowdownResult,
    CryogenicBlowdownState,
    TransientPipeWall,
    TransientTankWall,
    blowdown_state_to_flashing_droplet_source,
    blowdown_state_to_homogeneous_evaporation_source,
    blowdown_state_to_lh2_source,
    hem_blowdown_state_to_flashing_droplet_source,
    hem_blowdown_state_to_homogeneous_evaporation_source,
    run_cryogenic_blowdown,
)
from .pool_evaporation import (
    PoolEvaporationResult,
    PoolEvaporationStep,
    SolidSubstrate,
    semi_infinite_heat_flux,
    substrate_conduction_evaporation,
)
from .dynamic_pool import (
    ConstantHeatFluxSurface,
    DynamicPoolNumerics,
    DynamicPoolResult,
    DynamicPoolStep,
    simulate_axisymmetric_spreading_pool,
)
from .droplet_rainout import (
    DropletClass,
    DropletClassOutcome,
    DropletPopulationResult,
    DropletTransportInput,
    RainoutPoolCouplingResult,
    concurrent_rainout_pool,
    dynamic_rainout_pool,
    droplet_transport_input_from_flash,
    post_release_rainout_pool,
    transport_droplet_population,
)
from .finite_release import (
    FiniteReleasePuffHandoff,
    SteadyWindApplicability,
    assess_steady_wind_applicability,
    finite_release_puff_handoff,
)
from .finite_puff import (
    GaussianPuffConfig,
    GaussianPuffResult,
    GaussianPuffState,
    PuffReceptorTrace,
    integrate_finite_gaussian_puff,
)
from .lh2_droplets import (
    critical_droplet_diameter_for_phase_delay,
    critical_diffusivity_for_phase_delay,
    FlashingHydrogenDropletSource,
    gasflow_phase_relaxation_coefficient,
    HomogeneousEquilibriumHydrogenSource,
    flashing_hydrogen_droplet_source,
    homogeneous_equilibrium_hydrogen_source,
    homogeneous_equilibrium_hydrogen_source_from_postflash,
    minimum_heat_limited_hydrogen_evaporation_time,
)
from .cryogenic_air import (
    AirPhaseEquilibrium,
    N2O2CondensedPhaseScope,
    CondensedAirSource,
    CondensedAirStep,
    HydrogenEvaporationEndpoint,
    MultiphaseHydrogenSourcePlane,
    TwoVelocityRelaxation,
    TwoVelocityPhaseStep,
    LiZone3,
    air_saturation_pressure,
    cunningham_corrected_particle_relaxation_time,
    davies_cunningham_slip_correction,
    equilibrium_air_phase_split,
    n2o2_condensed_phase_scope,
    li2026_zone3,
    multiphase_hydrogen_evaporation_endpoint,
    multiphase_hydrogen_source_plane,
    minimum_heat_limited_sublimation_time,
    particle_relaxation_time,
    particle_terminal_velocity,
    particle_knudsen_number,
    schiller_naumann_relaxation_time,
    relax_two_velocity_drag,
    relax_two_velocity_schiller_naumann_drag,
    advance_two_velocity_phase_step,
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
    "YawedCrosswind", "YawedTrajectory",
    "FiniteWallJetCrosswind", "FiniteWallJetResult",
    "WallJetDiagnostics", "WallJetTransitionConfig",
    "FixedReceptor", "ReceptorTrace", "SteadyPlumeTable",
    "WindowStatistics", "WindHistory", "first_order_sensor_response",
    "meteorological_from_to_math_radians", "replay_fixed_receptors",
    "WindFrame", "AxisAlignedCuboid", "TransverseWall",
    "SegmentIntersection", "ObstacleEncounter", "ObstacleScreenResult",
    "screen_trajectory",
    "DenseGasFenceBenchmark", "DenseGasFenceObservation",
    "NeutralObstacleWakeBenchmark", "NeutralObstacleWakePoint",
    "pair_smedis_fence_trials", "read_aij_case_h", "read_smedis_fence_pair",
    "TransportEvidenceBoundary", "TURBULENCE_BOUNDARY", "PARTICLE_SLIP_BOUNDARY",
    "SourceEnthalpyBoundary", "HydrogenGasDeparture", "hydrogen_gas_departure",
    "NotionalNozzle", "notional_nozzle", "choking_pressure",
    "IsentropicThroat", "isentropic_throat",
    "EnergyConservingNotionalNozzle", "energy_conserving_notional_nozzle",
    "MeasuredThroatExpansion", "expand_measured_throat_to_ambient",
    "CryogenicBlowdownConfig", "CryogenicBlowdownState",
    "CryogenicBlowdownResult", "TransientTankWall", "TransientPipeWall",
    "blowdown_state_to_lh2_source",
    "blowdown_state_to_flashing_droplet_source",
    "blowdown_state_to_homogeneous_evaporation_source", "run_cryogenic_blowdown",
    "hem_blowdown_state_to_flashing_droplet_source",
    "hem_blowdown_state_to_homogeneous_evaporation_source",
    "SolidSubstrate", "PoolEvaporationStep", "PoolEvaporationResult",
    "semi_infinite_heat_flux", "substrate_conduction_evaporation",
    "ConstantHeatFluxSurface", "DynamicPoolNumerics", "DynamicPoolStep",
    "DynamicPoolResult",
    "simulate_axisymmetric_spreading_pool",
    "DropletClass", "DropletClassOutcome", "DropletPopulationResult",
    "DropletTransportInput", "RainoutPoolCouplingResult",
    "transport_droplet_population", "concurrent_rainout_pool",
    "dynamic_rainout_pool",
    "droplet_transport_input_from_flash",
    "post_release_rainout_pool",
    "FiniteReleasePuffHandoff", "SteadyWindApplicability",
    "finite_release_puff_handoff", "assess_steady_wind_applicability",
    "GaussianPuffConfig", "GaussianPuffResult", "GaussianPuffState",
    "PuffReceptorTrace", "integrate_finite_gaussian_puff",
    "FlashingHydrogenDropletSource", "flashing_hydrogen_droplet_source",
    "gasflow_phase_relaxation_coefficient",
    "critical_droplet_diameter_for_phase_delay",
    "critical_diffusivity_for_phase_delay",
    "HomogeneousEquilibriumHydrogenSource",
    "homogeneous_equilibrium_hydrogen_source",
    "homogeneous_equilibrium_hydrogen_source_from_postflash",
    "minimum_heat_limited_hydrogen_evaporation_time",
    "AirPhaseEquilibrium", "N2O2CondensedPhaseScope",
    "equilibrium_air_phase_split", "n2o2_condensed_phase_scope",
    "CondensedAirSource", "CondensedAirStep", "transported_condensed_air_source",
    "HydrogenEvaporationEndpoint", "multiphase_hydrogen_evaporation_endpoint",
    "MultiphaseHydrogenSourcePlane", "multiphase_hydrogen_source_plane",
    "TwoVelocityRelaxation", "relax_two_velocity_drag",
    "relax_two_velocity_schiller_naumann_drag",
    "TwoVelocityPhaseStep", "advance_two_velocity_phase_step",
    "LiZone3", "li2026_zone3", "air_saturation_pressure",
    "particle_knudsen_number", "davies_cunningham_slip_correction",
    "cunningham_corrected_particle_relaxation_time",
    "particle_relaxation_time", "schiller_naumann_relaxation_time",
    "particle_terminal_velocity",
    "ranz_marshall_transfer_number",
    "minimum_heat_limited_sublimation_time",
    "BuoyantClosure", "make_buoyant", "LAMBDA",
    "LiftoffPlume", "LiftoffResult", "LiftoffState",
    "richardson_liftoff", "liftoff_regime",
    "ALPHA", "BETA", "GAMMA", "RI_ONSET", "RI_SUBSTANTIAL", "RI_CLEAR",
]
