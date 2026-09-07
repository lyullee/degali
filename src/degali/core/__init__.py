"""Physics and numerics of the DEGADIS model."""

from .atmosphere import (
    StabilityDefaults,
    absolute_humidity,
    ambient_density,
    fit_alpha,
    friction_velocity,
    psi,
    richardson_star,
    richardson_thermal,
    stability_defaults,
    wind_log,
    wind_power,
)
from .blanket import Blanket, BlanketParameters, BlanketResult, BlanketState
from .crfg import SourceVectors, build_source_vectors
from .dose import DoseHistory, DoseRun, Receptor, receptor_times
from .downwind import Downwind, DownwindConstants, DownwindState, series
from .driver import DriverParameters, SourceRun
from .entrainment import entrainment_velocity, phi, phi_hat, surface_exchange
from .jetplume import JetCoefficients, JetPlume, JetResult, JetState, ellipse
from .observer import (
    Observer,
    ObserverKinematics,
    ObserverState,
    crossing_downwind,
    crossing_upwind,
    release_times,
)
from .steady import Profile, SteadyStateRun
from .szf import LayerState, sigma_z_over_source
from .thermo import (
    AdiabaticTable,
    AmbientConditions,
    CoolPropBackend,
    GasProperties,
    LegacyBackend,
    MixtureState,
    Thermo,
    make_backend,
)
from .timesort import Snapshot, TimeSort, sort_times
from .transient import ObserverResult, TransientResult, TransientRun

__all__ = [
    # atmosphere
    "StabilityDefaults", "absolute_humidity", "ambient_density", "fit_alpha",
    "friction_velocity", "psi", "richardson_star", "richardson_thermal",
    "stability_defaults", "wind_log", "wind_power",
    # thermodynamics
    "AdiabaticTable", "AmbientConditions", "CoolPropBackend", "GasProperties",
    "LegacyBackend", "MixtureState", "Thermo", "make_backend",
    # closures
    "entrainment_velocity", "phi", "phi_hat", "surface_exchange",
    # source
    "Blanket", "BlanketParameters", "BlanketResult", "BlanketState",
    "DriverParameters", "SourceRun", "SourceVectors", "build_source_vectors",
    "LayerState", "sigma_z_over_source",
    # downwind
    "Downwind", "DownwindConstants", "DownwindState", "series",
    "Profile", "SteadyStateRun",
    # transient
    "Observer", "ObserverKinematics", "ObserverState", "crossing_downwind",
    "crossing_upwind", "release_times",
    "ObserverResult", "TransientResult", "TransientRun",
    "Snapshot", "TimeSort", "sort_times",
    "DoseHistory", "DoseRun", "Receptor", "receptor_times",
    # jet/plume
    "JetCoefficients", "JetPlume", "JetResult", "JetState", "ellipse",
]
