"""Hecht--Panda Raman validation for cryogenic gaseous-hydrogen jets.

Only quantities printed in the public papers are scored.  The active benchmark
uses the final 2019 journal manuscript; the 2017 conference fits are retained
in the data file as provenance.  The comparison does
not digitize plotted curves and therefore cannot manufacture precision that
the source does not provide.  See ``docs/prereg-hecht-panda-raman.md`` for the
protocol frozen before the first model run.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq

from ..addons.axisymmetric_jet import (
    AxisymmetricJetSource,
    ConservedGaussianJet,
    entrain_and_heat_initial_plug,
)
from ..addons.notional import expand_measured_throat_to_ambient
from .nearfield import hydrogen_gas_jet


DATA_PATH = (
    Path(__file__).resolve().parents[3]
    / "reference" / "lh2" / "hecht-panda-raman-2017.json"
)
SOURCE_HEIGHT = 1.0
SAMPLE_MM = np.arange(40.0, 100.0 + 0.5, 1.0)


def _reported_sample_mm(condition: dict) -> np.ndarray:
    """One-mm samples over the stitched range implied by Table 1."""
    end = 40.0 + 10.0 * int(condition["heights"])
    return np.arange(40.0, end + 0.5, 1.0)


@dataclass(frozen=True)
class RamanVariant:
    """One pre-registered JETPLU closure combination."""

    name: str
    alfa1: float
    sc: float = 1.42
    density_scaled_entrainment: bool = False
    momentum_entrainment_beta: float = 0.0


VARIANTS = {
    "legacy": RamanVariant("legacy", alfa1=0.057),
    "current_lh2": RamanVariant(
        "current_lh2", alfa1=0.0875, density_scaled_entrainment=True
    ),
    "alpha_only": RamanVariant("alpha_only", alfa1=0.0875),
    "hyram_momentum": RamanVariant(
        "hyram_momentum", alfa1=0.057, sc=1.16**2,
        momentum_entrainment_beta=0.28,
    ),
}


@dataclass(frozen=True)
class RamanMetrics:
    """Aggregate slopes and profile diagnostics for one model variant."""

    variant: str
    cases: int
    points: int
    centerline_mass_slope: float
    mass_half_width_slope_mm: float
    mass_radial_coefficient: float
    centerline_temperature_slope: float
    temperature_half_width_slope_mm: float
    temperature_radial_coefficient: float
    centerline_mass_trend_rmse: float
    centerline_temperature_trend_rmse: float
    relative_errors: dict[str, float]

    @property
    def slope_passes(self) -> int:
        """Number of the four printed slopes within the frozen 25% band."""
        return sum(abs(value) <= 0.25 for value in self.relative_errors.values())


@dataclass(frozen=True)
class RamanComparison:
    """The observations and all pre-registered model variants."""

    observed: dict[str, float]
    variants: dict[str, RamanMetrics]

    def report(self) -> str:
        headings = (
            "variant", "mass decay", "mass width", "T decay", "T width",
            "A_Y", "A_T", "passes",
        )
        lines = [" | ".join(headings), " | ".join("---" for _ in headings)]
        for name, metric in self.variants.items():
            lines.append(" | ".join([
                name,
                f"{metric.centerline_mass_slope:.5f}",
                f"{metric.mass_half_width_slope_mm:.5f}",
                f"{metric.centerline_temperature_slope:.5f}",
                f"{metric.temperature_half_width_slope_mm:.5f}",
                f"{metric.mass_radial_coefficient:.1f}",
                f"{metric.temperature_radial_coefficient:.1f}",
                f"{metric.slope_passes}/4",
            ]))
        return "\n".join(lines)


def load_data(path: str | Path = DATA_PATH) -> dict:
    """Load Table 1 and the final-journal Figures 6 and 8 transcription."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _forced_origin_slope(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.dot(x, y) / np.dot(x, x))


def _slope_with_intercept(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Slope and intercept; the paper prints the former but not the latter."""
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def _profile_half_width(table, cc: float, sigma: float, response: str) -> float:
    """Radius where mass fraction or temperature excursion falls to half."""
    centre = table.from_concentration(cc)
    if response == "mass":
        denominator = centre.wc

        def normalized(fraction):
            return table.from_concentration(cc * fraction).wc / denominator
    elif response == "temperature":
        ambient_temperature = float(table.t[0])
        denominator = centre.temp - ambient_temperature

        def normalized(fraction):
            state = table.from_concentration(cc * fraction)
            return (state.temp - ambient_temperature) / denominator
    else:
        raise ValueError(f"unknown profile response {response!r}")

    concentration_ratio = brentq(
        lambda fraction: normalized(fraction) - 0.5,
        0.0, 1.0, xtol=1.0e-12,
    )
    return sigma * math.sqrt(-2.0 * math.log(concentration_ratio))


def _case_arrays(
    condition: dict, variant: RamanVariant, data: dict,
    *, sample_mm: np.ndarray = SAMPLE_MM,
) -> dict:
    """Run and sample one Table-1 source over the published camera range."""
    from CoolProp.CoolProp import PropsSI

    diameter = float(condition["diameter_mm"]) / 1000.0
    throat_area = math.pi * diameter**2 / 4.0
    rate = (
        float(condition["rho_throat_kg_m3"])
        * float(condition["v_throat_m_s"])
        * throat_area
    )
    source = expand_measured_throat_to_ambient(
        fluid="Hydrogen",
        mass_flow=rate,
        throat_diameter=diameter,
        throat_pressure=float(condition["P_throat_bar_abs"]) * 1.0e5,
        throat_temperature=float(condition["T_throat_K"]),
        throat_density=float(condition["rho_throat_kg_m3"]),
        throat_velocity=float(condition["v_throat_m_s"]),
        ambient_pressure=float(data["ambient_pressure_Pa"]),
    )
    plume, y0 = hydrogen_gas_jet(
        rate=rate,
        diameter=source.diameter,
        velocity=source.velocity,
        wind=0.0,
        height=SOURCE_HEIGHT,
        source_temperature=source.temperature,
        source_density=source.density,
        theta=math.pi / 2.0,
        relative_humidity=0.0,
        ambient_temperature=float(data["ambient_temperature_K"]),
        ambient_pressure=float(data["ambient_pressure_Pa"]),
        alfa1=variant.alfa1,
        sc=variant.sc,
        density_scaled_entrainment=variant.density_scaled_entrainment,
        momentum_entrainment_beta=variant.momentum_entrainment_beta,
    )
    result = plume.run(y0, distmx=0.001, tol=1.0e-5, smax=0.13)
    if len(result.rows) < 3:
        raise RuntimeError("cryogenic free-jet integration returned no profile")

    rows = result.rows
    z_model = rows[:, 1] - SOURCE_HEIGHT
    z = sample_mm / 1000.0
    if z_model[0] > z[0] or z_model[-1] < z[-1]:
        raise RuntimeError(
            f"model range {z_model[0]:.4g}--{z_model[-1]:.4g} m does not "
            "cover the Raman camera range"
        )
    cc = np.interp(z, z_model, rows[:, 2])
    sy = np.interp(z, z_model, rows[:, 3])
    sz = np.interp(z, z_model, rows[:, 4])
    sigma = np.sqrt(sy * sz)

    mass_fraction = np.array([
        plume.th.table.from_concentration(value).wc for value in cc
    ])
    temperature = np.array([
        plume.th.table.from_concentration(value).temp for value in cc
    ])
    mass_half_width = np.array([
        _profile_half_width(plume.th.table, c, width, "mass")
        for c, width in zip(cc, sigma)
    ])
    temperature_half_width = np.array([
        _profile_half_width(plume.th.table, c, width, "temperature")
        for c, width in zip(cc, sigma)
    ])

    stagnation_density = float(PropsSI(
        "D", "T", float(condition["T_nozzle_K"]),
        "P", float(condition["P_nozzle_bar_abs"]) * 1.0e5, "Hydrogen"
    ))
    ambient_density = plume.rhoa
    zeta = z / (diameter * math.sqrt(stagnation_density / ambient_density))
    xi = z / (0.5 * diameter)
    return {
        "z": z,
        "zeta": zeta,
        "xi": xi,
        "inverse_mass": 1.0 / mass_fraction,
        "mass_half_width_mm": 1000.0 * mass_half_width,
        "mass_radial_coefficient": math.log(2.0) / (mass_half_width / z) ** 2,
        "inverse_temperature": float(data["ambient_temperature_K"]) /
            (float(data["ambient_temperature_K"]) - temperature),
        "temperature_half_width_mm": 1000.0 * temperature_half_width,
        "temperature_radial_coefficient": math.log(2.0) /
            (temperature_half_width / z) ** 2,
    }


def _conserved_energy_case_arrays(
    condition: dict, data: dict,
    *, conservative_establishment: bool | str = False,
    equilibrium_air_condensation: bool = False,
    radial_points: int | None = None,
    thermodynamic_spreading_ratio: float | None = None,
    thermodynamic_profile: str = "density",
    ideal_gas_enthalpy_fluids: tuple[str, str] | None = None,
    initial_heating_temperature: float | None = None,
    ambient_absolute_humidity: float = 0.0,
    sample_mm: np.ndarray = SAMPLE_MM,
    integration_maximum_step: float | None = None,
    integration_relative_tolerance: float | None = None,
    ambient_coflow_velocity: float = 0.0,
    temperature_dependent_phase_enthalpy: bool = False,
    hydrogen_enthalpy_species: str = "Hydrogen",
    equilibrium_argon_condensation: bool = False,
    radiative_absorptivity: float = 0.0,
) -> dict:
    """Run the published five-balance axisymmetric model for one test."""
    from CoolProp.CoolProp import PropsSI

    diameter = float(condition["diameter_mm"]) / 1000.0
    throat_area = math.pi * diameter**2 / 4.0
    rate = (
        float(condition["rho_throat_kg_m3"])
        * float(condition["v_throat_m_s"])
        * throat_area
    )
    outlet = expand_measured_throat_to_ambient(
        fluid="Hydrogen",
        mass_flow=rate,
        throat_diameter=diameter,
        throat_pressure=float(condition["P_throat_bar_abs"]) * 1.0e5,
        throat_temperature=float(condition["T_throat_K"]),
        throat_density=float(condition["rho_throat_kg_m3"]),
        throat_velocity=float(condition["v_throat_m_s"]),
        ambient_pressure=float(data["ambient_pressure_Pa"]),
    )
    ambient_temperature = float(data["ambient_temperature_K"])
    ambient_pressure = float(data["ambient_pressure_Pa"])
    jet_source = AxisymmetricJetSource(
        diameter=outlet.diameter,
        velocity=outlet.velocity,
        density=outlet.density,
        temperature=outlet.temperature,
    )
    initial_heating = None
    if initial_heating_temperature is not None:
        initial_heating = entrain_and_heat_initial_plug(
            jet_source,
            minimum_temperature=initial_heating_temperature,
            ambient_temperature=ambient_temperature,
            ambient_pressure=ambient_pressure,
            ambient_density=float(PropsSI(
                "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
            )),
            fuel_molecular_weight=float(PropsSI("M", "Hydrogen")),
            ambient_molecular_weight=float(PropsSI("M", "Air")),
            fuel_heat_capacity=float(PropsSI(
                "C", "T", ambient_temperature, "P", ambient_pressure,
                "Hydrogen"
            )),
            ambient_heat_capacity=float(PropsSI(
                "C", "T", ambient_temperature, "P", ambient_pressure, "Air"
            )),
        )
        jet_source = initial_heating.source
    ambient_molecular_weight = float(PropsSI("M", "Air"))
    water_molecular_weight = float(PropsSI("M", "Water"))
    humid_molecular_weight = (
        (1.0 + ambient_absolute_humidity)
        / (
            1.0 / ambient_molecular_weight
            + ambient_absolute_humidity / water_molecular_weight
        )
    )
    dry_ambient_density = float(PropsSI(
        "D", "T", ambient_temperature, "P", ambient_pressure, "Air"
    ))
    ambient_density = (
        dry_ambient_density * humid_molecular_weight / ambient_molecular_weight
    )
    model = ConservedGaussianJet(
        jet_source,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        ambient_density=ambient_density,
        fuel_molecular_weight=float(PropsSI("M", "Hydrogen")),
        ambient_molecular_weight=ambient_molecular_weight,
        fuel_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure, "Hydrogen"
        )),
        ambient_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure, "Air"
        )),
        conservative_establishment=conservative_establishment,
        equilibrium_air_condensation=equilibrium_air_condensation,
        thermodynamic_spreading_ratio=thermodynamic_spreading_ratio,
        thermodynamic_profile=thermodynamic_profile,
        ideal_gas_enthalpy_fluids=ideal_gas_enthalpy_fluids,
        ambient_absolute_humidity=ambient_absolute_humidity,
        ambient_coflow_velocity=ambient_coflow_velocity,
        temperature_dependent_phase_enthalpy=(
            temperature_dependent_phase_enthalpy
        ),
        hydrogen_enthalpy_species=hydrogen_enthalpy_species,
        equilibrium_argon_condensation=equilibrium_argon_condensation,
        radiative_absorptivity=radiative_absorptivity,
        radial_points=(
            radial_points if radial_points is not None
            else 81 if equilibrium_air_condensation else 241
        ),
    )
    result = model.solve(
        maximum_distance=0.13,
        maximum_step=(
            integration_maximum_step
            if integration_maximum_step is not None
            else 0.002 if equilibrium_air_condensation else 0.001
        ),
        relative_tolerance=(
            integration_relative_tolerance
            if integration_relative_tolerance is not None
            else 5.0e-6 if equilibrium_air_condensation else 1.0e-6
        ),
        method="LSODA" if equilibrium_air_condensation else "RK45",
    )
    z = sample_mm / 1000.0
    if result.y[0] > z[0] or result.y[-1] < z[-1]:
        raise RuntimeError(
            f"conserved model range {result.y[0]:.4g}--{result.y[-1]:.4g} m "
            "does not cover the Raman camera range"
        )
    states = np.array([result.state_at_y(value) for value in z])
    mass_fraction = states[:, 3]
    temperature = np.array([
        model.centreline_temperature(state) for state in states
    ])
    mass_half_width = np.array([
        model.half_width(state, "mass") for state in states
    ])
    temperature_half_width = np.array([
        model.half_width(state, "temperature") for state in states
    ])

    stagnation_density = float(PropsSI(
        "D", "T", float(condition["T_nozzle_K"]),
        "P", float(condition["P_nozzle_bar_abs"]) * 1.0e5, "Hydrogen"
    ))
    zeta = z / (
        diameter * math.sqrt(stagnation_density / model.ambient_density)
    )
    xi = z / (0.5 * diameter)

    # Deliberately impossible blackbody upper bound: every station's full 5B
    # envelope absorbs as though it were at the colder centre temperature.
    sigma = 5.670374419e-8
    radiative_power_per_length = (
        2.0 * math.pi * model.radial_limit * result.width * sigma
        * np.maximum(ambient_temperature**4 - result.temperature**4, 0.0)
    )
    blackbody_heat = float(trapezoid(radiative_power_per_length, result.S))
    initial_state = np.array([
        result.velocity[0], result.width[0], result.density[0],
        result.mass_fraction[0], result.theta[0], result.x[0], result.y[0],
    ])
    velocity_i, density_i, _fraction_i, _mw_i, _temp_i, rho_h_i = (
        model._profiles(initial_state)
    )
    area_weight_i = 2.0 * math.pi * initial_state[1]**2 * model._eta
    initial_enthalpy_deficit = abs(float(trapezoid(
        (
            rho_h_i - model._ambient_enthalpy * density_i
        ) * velocity_i * area_weight_i,
        model._eta,
    )))
    h_first = model._mixture_enthalpy(
        result.temperature[0], result.mass_fraction[0]
    )
    h_last = model._mixture_enthalpy(
        result.temperature[-1], result.mass_fraction[-1]
    )
    centreline_warming_scale = (
        jet_source.fuel_mass_flow * abs(h_last - h_first)
    )
    return {
        "z": z,
        "zeta": zeta,
        "xi": xi,
        "inverse_mass": 1.0 / mass_fraction,
        "mass_half_width_mm": 1000.0 * mass_half_width,
        "mass_radial_coefficient": math.log(2.0) /
            (mass_half_width / z) ** 2,
        "inverse_temperature": ambient_temperature /
            (ambient_temperature - temperature),
        "temperature_half_width_mm": 1000.0 * temperature_half_width,
        "temperature_radial_coefficient": math.log(2.0) /
            (temperature_half_width / z) ** 2,
        "diagnostics": {
            "maximum_initial_heating_residual": (
                0.0 if initial_heating is None else max(
                    initial_heating.fuel_mass_residual,
                    initial_heating.momentum_residual,
                    initial_heating.energy_residual,
                    initial_heating.geometry_mass_residual,
                )
            ),
            "maximum_boundary_residual": max(
                abs(value) for value in model.establishment_residuals.values()
            ),
            "maximum_species_drift": float(np.max(np.abs(
                result.species_flux / result.species_flux[0] - 1.0
            ))),
            "maximum_energy_drift": float(np.max(np.abs(
                (
                    result.energy_flux
                    - 0.5 * ambient_coflow_velocity**2 * result.mass_flux
                    - result.radiative_heat_added
                ) / (
                    result.energy_flux[0]
                    - 0.5 * ambient_coflow_velocity**2 * result.mass_flux[0]
                    - result.radiative_heat_added[0]
                ) - 1.0
            ))),
            "minimum_velocity_m_s": float(np.min(result.velocity)),
            "minimum_width_m": float(np.min(result.width)),
            "minimum_density_kg_m3": float(np.min(result.density)),
            "minimum_mass_fraction": float(np.min(result.mass_fraction)),
            "maximum_mass_fraction": float(np.max(result.mass_fraction)),
            "minimum_temperature_K": float(np.min(result.temperature)),
            "maximum_temperature_K": float(np.max(result.temperature)),
            "blackbody_heat_W": blackbody_heat,
            "blackbody_to_initial_enthalpy_deficit": (
                blackbody_heat / max(initial_enthalpy_deficit, 1.0e-30)
            ),
            "blackbody_to_centreline_warming_scale": (
                blackbody_heat / max(centreline_warming_scale, 1.0e-30)
            ),
        },
    }


def _aggregate_metrics(
    name: str, arrays: list[dict], data: dict
) -> RamanMetrics:
    """Apply the frozen regressions to one model's nine sampled cases."""
    take = lambda key: np.concatenate([case[key] for case in arrays])
    zeta, xi = take("zeta"), take("xi")
    inverse_mass = take("inverse_mass")
    inverse_temperature = take("inverse_temperature")
    mass_width = take("mass_half_width_mm")
    temperature_width = take("temperature_half_width_mm")

    mass_slope, mass_intercept = _slope_with_intercept(zeta, inverse_mass)
    mass_width_slope = _forced_origin_slope(xi, mass_width)
    temperature_slope, temperature_intercept = _slope_with_intercept(
        zeta, inverse_temperature
    )
    temperature_width_slope = _forced_origin_slope(xi, temperature_width)
    observed = data["published_fits"]
    fitted_mass = mass_intercept + mass_slope * zeta
    fitted_temperature = temperature_intercept + temperature_slope * zeta
    errors = {
        "centerline_mass": mass_slope /
            observed["centerline_inverse_mass_fraction_slope"] - 1.0,
        "mass_half_width": mass_width_slope /
            observed["mass_fraction_half_width_slope_mm"] - 1.0,
        "centerline_temperature": temperature_slope /
            observed["centerline_inverse_temperature_slope"] - 1.0,
        "temperature_half_width": temperature_width_slope /
            observed["temperature_half_width_slope_mm"] - 1.0,
    }
    return RamanMetrics(
        variant=name,
        cases=len(arrays),
        points=len(zeta),
        centerline_mass_slope=mass_slope,
        mass_half_width_slope_mm=mass_width_slope,
        mass_radial_coefficient=float(np.median(take("mass_radial_coefficient"))),
        centerline_temperature_slope=temperature_slope,
        temperature_half_width_slope_mm=temperature_width_slope,
        temperature_radial_coefficient=float(np.median(
            take("temperature_radial_coefficient")
        )),
        centerline_mass_trend_rmse=float(np.sqrt(np.mean(
            (inverse_mass / fitted_mass - 1.0) ** 2
        ))),
        centerline_temperature_trend_rmse=float(np.sqrt(np.mean(
            (inverse_temperature / fitted_temperature - 1.0) ** 2
        ))),
        relative_errors=errors,
    )


def evaluate_variant(
    variant: str | RamanVariant, *, data: dict | None = None
) -> RamanMetrics:
    """Evaluate one frozen closure combination on all nine release states."""
    data = load_data() if data is None else data
    if isinstance(variant, str):
        variant = VARIANTS[variant]
    arrays = [
        _case_arrays(condition, variant, data)
        for condition in data["conditions"]
    ]

    return _aggregate_metrics(variant.name, arrays, data)


def evaluate_conserved_energy(*, data: dict | None = None) -> RamanMetrics:
    """Evaluate the independently implemented published energy model."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(condition, data)
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("conserved_energy", arrays, data)


def evaluate_conservative_establishment(
    *, data: dict | None = None
) -> RamanMetrics:
    """Evaluate the four-flux conservative development-zone boundary."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data, conservative_establishment=True
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("conservative_establishment", arrays, data)


def evaluate_scalar_constrained_establishment(
    *, data: dict | None = None
) -> RamanMetrics:
    """Evaluate the scalar-peak-constrained conservative boundary."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data, conservative_establishment="scalar_peak"
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("scalar_constrained_establishment", arrays, data)


def evaluate_distributed_air_condensation(
    *, data: dict | None = None, radial_points: int = 81
) -> RamanMetrics:
    """Evaluate the frozen equilibrium/no-slip condensed-air bound."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            radial_points=radial_points,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("distributed_air_condensation", arrays, data)


def saturated_ambient_absolute_humidity(data: dict | None = None) -> float:
    """Return the 100%-RH water/dry-air mass ratio at the test ambient."""
    from CoolProp.CoolProp import PropsSI

    data = load_data() if data is None else data
    temperature = float(data["ambient_temperature_K"])
    pressure = float(data["ambient_pressure_Pa"])
    vapour_pressure = float(PropsSI("P", "T", temperature, "Q", 0, "Water"))
    return float(
        PropsSI("M", "Water") / PropsSI("M", "Air")
        * vapour_pressure / (pressure - vapour_pressure)
    )


def evaluate_humid_air_frost(
    relative_humidity: float,
    *, data: dict | None = None, radial_points: int = 81,
) -> RamanMetrics:
    """Evaluate equilibrium frost feedback at an explicit ambient RH."""
    if not 0.0 <= relative_humidity <= 100.0:
        raise ValueError("relative humidity must lie between 0 and 100 percent")
    data = load_data() if data is None else data
    humidity = (
        relative_humidity / 100.0 * saturated_ambient_absolute_humidity(data)
    )
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            ambient_absolute_humidity=humidity,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics(
        f"humid_air_frost_{relative_humidity:g}pct", arrays, data
    )


def evaluate_temperature_dependent_phase_enthalpy(
    *, data: dict | None = None, radial_points: int = 81
) -> RamanMetrics:
    """Evaluate component h(T) inside the dry equilibrium phase model."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics(
        "temperature_dependent_phase_enthalpy", arrays, data
    )


def evaluate_argon_phase_completeness_audit(
    *, data: dict | None = None, radial_points: int = 81
) -> dict:
    """Evaluate argon once and return both coverage protocols and diagnostics."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            equilibrium_argon_condensation=True,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    reported_arrays = []
    for condition, case in zip(data["conditions"], arrays):
        mask = case["z"] <= _reported_sample_mm(condition)[-1] / 1000.0
        reported_arrays.append({
            key: (
                value[mask]
                if isinstance(value, np.ndarray)
                and value.ndim >= 1 and len(value) == len(mask)
                else value
            )
            for key, value in case.items()
        })
    diagnostics = [case["diagnostics"] for case in arrays]
    return {
        "original_549": _aggregate_metrics(
            "argon_phase_completeness_549", arrays, data
        ),
        "reported_369": _aggregate_metrics(
            "argon_phase_completeness", reported_arrays, data
        ),
        "diagnostics": {
            "maximum_boundary_residual": max(
                value["maximum_boundary_residual"] for value in diagnostics
            ),
            "maximum_species_drift": max(
                value["maximum_species_drift"] for value in diagnostics
            ),
            "maximum_energy_drift": max(
                value["maximum_energy_drift"] for value in diagnostics
            ),
            "minimum_temperature_K": min(
                value["minimum_temperature_K"] for value in diagnostics
            ),
            "maximum_temperature_K": max(
                value["maximum_temperature_K"] for value in diagnostics
            ),
            "minimum_velocity_m_s": min(
                value["minimum_velocity_m_s"] for value in diagnostics
            ),
            "minimum_width_m": min(
                value["minimum_width_m"] for value in diagnostics
            ),
            "minimum_density_kg_m3": min(
                value["minimum_density_kg_m3"] for value in diagnostics
            ),
            "minimum_mass_fraction": min(
                value["minimum_mass_fraction"] for value in diagnostics
            ),
            "maximum_mass_fraction": max(
                value["maximum_mass_fraction"] for value in diagnostics
            ),
        },
    }


def evaluate_argon_phase_completeness(
    *, data: dict | None = None, radial_points: int = 81
) -> RamanMetrics:
    """Evaluate argon under the primary unequal 369-point coverage."""
    return evaluate_argon_phase_completeness_audit(
        data=data, radial_points=radial_points
    )["reported_369"]


def evaluate_hydrogen_spin_isomer_enthalpy_audit(
    *, data: dict | None = None, radial_points: int = 81,
    hydrogen_enthalpy_species: str = "ParaHydrogen",
) -> dict:
    """Run a fixed H2 spin-isomer caloric sensitivity on both coverages."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            hydrogen_enthalpy_species=hydrogen_enthalpy_species,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    reported_arrays = []
    for condition, case in zip(data["conditions"], arrays):
        mask = case["z"] <= _reported_sample_mm(condition)[-1] / 1000.0
        reported_arrays.append({
            key: (
                value[mask]
                if isinstance(value, np.ndarray)
                and value.ndim >= 1 and len(value) == len(mask)
                else value
            )
            for key, value in case.items()
        })
    diagnostics = [case["diagnostics"] for case in arrays]
    label = hydrogen_enthalpy_species.lower()
    return {
        "original_549": _aggregate_metrics(
            f"{label}_enthalpy_549", arrays, data
        ),
        "reported_369": _aggregate_metrics(
            f"{label}_enthalpy", reported_arrays, data
        ),
        "diagnostics": {
            "maximum_boundary_residual": max(
                value["maximum_boundary_residual"] for value in diagnostics
            ),
            "maximum_species_drift": max(
                value["maximum_species_drift"] for value in diagnostics
            ),
            "maximum_energy_drift": max(
                value["maximum_energy_drift"] for value in diagnostics
            ),
            "minimum_temperature_K": min(
                value["minimum_temperature_K"] for value in diagnostics
            ),
            "maximum_temperature_K": max(
                value["maximum_temperature_K"] for value in diagnostics
            ),
            "minimum_velocity_m_s": min(
                value["minimum_velocity_m_s"] for value in diagnostics
            ),
            "minimum_width_m": min(
                value["minimum_width_m"] for value in diagnostics
            ),
            "minimum_density_kg_m3": min(
                value["minimum_density_kg_m3"] for value in diagnostics
            ),
            "minimum_mass_fraction": min(
                value["minimum_mass_fraction"] for value in diagnostics
            ),
            "maximum_mass_fraction": max(
                value["maximum_mass_fraction"] for value in diagnostics
            ),
        },
    }


def evaluate_four_flux_phase_two_scalar_audit(
    *, data: dict | None = None, radial_points: int = 81,
) -> dict:
    """Audit the fully conservative phase/two-scalar closure."""
    data = load_data() if data is None else data
    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            thermodynamic_spreading_ratio=thermal_ratio,
            thermodynamic_profile="temperature_mass_fraction",
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    reported_arrays = []
    for condition, case in zip(data["conditions"], arrays):
        mask = case["z"] <= _reported_sample_mm(condition)[-1] / 1000.0
        reported_arrays.append({
            key: (
                value[mask]
                if isinstance(value, np.ndarray)
                and value.ndim >= 1 and len(value) == len(mask)
                else value
            )
            for key, value in case.items()
        })
    diagnostics = [case["diagnostics"] for case in arrays]
    return {
        "original_549": _aggregate_metrics(
            "four_flux_phase_two_scalar_549", arrays, data
        ),
        "reported_369": _aggregate_metrics(
            "four_flux_phase_two_scalar", reported_arrays, data
        ),
        "diagnostics": {
            "maximum_boundary_residual": max(
                value["maximum_boundary_residual"] for value in diagnostics
            ),
            "maximum_species_drift": max(
                value["maximum_species_drift"] for value in diagnostics
            ),
            "maximum_energy_drift": max(
                value["maximum_energy_drift"] for value in diagnostics
            ),
            "minimum_temperature_K": min(
                value["minimum_temperature_K"] for value in diagnostics
            ),
            "maximum_temperature_K": max(
                value["maximum_temperature_K"] for value in diagnostics
            ),
            "minimum_velocity_m_s": min(
                value["minimum_velocity_m_s"] for value in diagnostics
            ),
            "minimum_width_m": min(
                value["minimum_width_m"] for value in diagnostics
            ),
            "minimum_density_kg_m3": min(
                value["minimum_density_kg_m3"] for value in diagnostics
            ),
            "minimum_mass_fraction": min(
                value["minimum_mass_fraction"] for value in diagnostics
            ),
            "maximum_mass_fraction": max(
                value["maximum_mass_fraction"] for value in diagnostics
            ),
        },
    }


def evaluate_blackbody_radiation_upper_bound(
    *, data: dict | None = None, radial_points: int = 81
) -> dict[str, float]:
    """Return the pre-registered perfect-absorber radiation ratios."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    deficit = np.array([
        case["diagnostics"]["blackbody_to_initial_enthalpy_deficit"]
        for case in arrays
    ])
    warming = np.array([
        case["diagnostics"]["blackbody_to_centreline_warming_scale"]
        for case in arrays
    ])
    heat = np.array([
        case["diagnostics"]["blackbody_heat_W"] for case in arrays
    ])
    return {
        "cases": float(len(arrays)),
        "maximum_blackbody_heat_W": float(np.max(heat)),
        "maximum_ratio_to_initial_enthalpy_deficit": float(np.max(deficit)),
        "maximum_ratio_to_centreline_warming_scale": float(np.max(warming)),
    }


def evaluate_radiative_blackbody_ode(
    *, data: dict | None = None, radial_points: int = 81
) -> dict:
    """Run the full-black external-energy ODE and both coverage protocols."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment="scalar_peak",
            equilibrium_air_condensation=True,
            temperature_dependent_phase_enthalpy=True,
            radiative_absorptivity=1.0,
            radial_points=radial_points,
            integration_maximum_step=0.00025,
            integration_relative_tolerance=5.0e-8,
        )
        for condition in data["conditions"]
    ]
    reported_arrays = []
    for condition, case in zip(data["conditions"], arrays):
        mask = case["z"] <= _reported_sample_mm(condition)[-1] / 1000.0
        reported_arrays.append({
            key: (
                value[mask]
                if isinstance(value, np.ndarray)
                and value.ndim >= 1 and len(value) == len(mask)
                else value
            )
            for key, value in case.items()
        })
    diagnostics = [case["diagnostics"] for case in arrays]
    return {
        "original_549": _aggregate_metrics(
            "radiative_blackbody_549", arrays, data
        ),
        "reported_369": _aggregate_metrics(
            "radiative_blackbody", reported_arrays, data
        ),
        "diagnostics": {
            "maximum_boundary_residual": max(
                value["maximum_boundary_residual"] for value in diagnostics
            ),
            "maximum_species_drift": max(
                value["maximum_species_drift"] for value in diagnostics
            ),
            "maximum_energy_drift_after_external_heat": max(
                value["maximum_energy_drift"] for value in diagnostics
            ),
            "minimum_temperature_K": min(
                value["minimum_temperature_K"] for value in diagnostics
            ),
            "maximum_temperature_K": max(
                value["maximum_temperature_K"] for value in diagnostics
            ),
        },
    }


def radiation_optical_depth_sensitivity(
    blackbody_ratio: float,
    optical_depths: tuple[float, ...] = (0.001, 0.01, 0.1, 0.5, 1.0),
) -> dict:
    """Scale a blackbody energy ratio by ``1-exp(-optical_depth)``."""
    if blackbody_ratio < 0.0:
        raise ValueError("blackbody ratio cannot be negative")
    values = {
        str(depth): blackbody_ratio * (1.0 - math.exp(-depth))
        for depth in optical_depths
    }
    threshold_absorptivity = min(0.01 / max(blackbody_ratio, 1.0e-30), 1.0)
    threshold_depth = (
        math.inf if threshold_absorptivity >= 1.0
        else -math.log1p(-threshold_absorptivity)
    )
    return {
        "ratios": values,
        "absorptivity_for_one_percent": threshold_absorptivity,
        "optical_depth_for_one_percent": threshold_depth,
    }


def evaluate_humid_air_frost_upper_bound(
    *, data: dict | None = None, radial_points: int = 81
) -> RamanMetrics:
    """Evaluate equilibrium frost feedback at the fixed 100%-RH upper bound."""
    metric = evaluate_humid_air_frost(
        100.0, data=data, radial_points=radial_points
    )
    return RamanMetrics(
        **{**metric.__dict__, "variant": "humid_air_frost_upper_bound"}
    )


def evaluate_differential_scalar_spreading(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate independently fixed turbulent heat/species diffusivities."""
    data = load_data() if data is None else data
    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            radial_points=radial_points,
            thermodynamic_spreading_ratio=thermal_ratio,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("differential_scalar_spreading", arrays, data)


def evaluate_measured_coflow(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate the reported 0.3 m/s honeycomb co-flow boundary."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            ambient_coflow_velocity=0.3,
            radial_points=radial_points,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("measured_0p3_m_s_coflow", arrays, data)


def evaluate_independent_mass_fraction_temperature(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate independent Gaussian mass-fraction/temperature profiles."""
    data = load_data() if data is None else data
    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            radial_points=radial_points,
            thermodynamic_spreading_ratio=thermal_ratio,
            thermodynamic_profile="temperature_mass_fraction",
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics(
        "independent_mass_fraction_temperature", arrays, data
    )


def evaluate_temperature_dependent_enthalpy(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate the frozen two-scalar model with ideal-gas h(T)."""
    data = load_data() if data is None else data
    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            radial_points=radial_points,
            thermodynamic_spreading_ratio=thermal_ratio,
            thermodynamic_profile="temperature_mass_fraction",
            ideal_gas_enthalpy_fluids=("Hydrogen", "Air"),
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("temperature_dependent_enthalpy", arrays, data)


def evaluate_initial_entrainment_heating(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate initial plug heating to the fixed dry-air dew point."""
    data = load_data() if data is None else data
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            radial_points=radial_points,
            initial_heating_temperature=82.15083337712093,
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("initial_entrainment_heating", arrays, data)


def evaluate_heated_two_scalar_jet(
    *, data: dict | None = None, radial_points: int = 121
) -> RamanMetrics:
    """Evaluate the fixed dew-point heating plus independent Y/T profiles."""
    data = load_data() if data is None else data
    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    arrays = [
        _conserved_energy_case_arrays(
            condition, data,
            conservative_establishment=True,
            radial_points=radial_points,
            initial_heating_temperature=82.15083337712093,
            thermodynamic_spreading_ratio=thermal_ratio,
            thermodynamic_profile="temperature_mass_fraction",
        )
        for condition in data["conditions"]
    ]
    return _aggregate_metrics("heated_two_scalar_jet", arrays, data)


def compare_reported_case_coverage() -> RamanComparison:
    """Apply the source-corrected unequal stitched-image coverage."""
    data = load_data()
    observed = {
        key: float(value)
        for key, value in data["published_fits"].items()
        if isinstance(value, (int, float))
    }

    def conserved_arrays(**options) -> list[dict]:
        return [
            _conserved_energy_case_arrays(
                condition, data,
                sample_mm=_reported_sample_mm(condition),
                **options,
            )
            for condition in data["conditions"]
        ]

    thermal_ratio = 1.16 * math.sqrt(0.70 / 0.85)
    return RamanComparison(
        observed=observed,
        variants={
            "published_establishment": _aggregate_metrics(
                "published_establishment", conserved_arrays(), data
            ),
            "conservative_establishment": _aggregate_metrics(
                "conservative_establishment",
                conserved_arrays(conservative_establishment=True), data,
            ),
            "heated_two_scalar_jet": _aggregate_metrics(
                "heated_two_scalar_jet",
                conserved_arrays(
                    conservative_establishment=True,
                    initial_heating_temperature=82.15083337712093,
                    thermodynamic_spreading_ratio=thermal_ratio,
                    thermodynamic_profile="temperature_mass_fraction",
                ),
                data,
            ),
        },
    )


def compare() -> RamanComparison:
    """Run all pre-registered variants and return an auditable comparison."""
    data = load_data()
    observed = {
        key: float(value)
        for key, value in data["published_fits"].items()
        if isinstance(value, (int, float))
    }
    return RamanComparison(
        observed=observed,
        variants={
            **{
                name: evaluate_variant(spec, data=data)
                for name, spec in VARIANTS.items()
            },
            "conserved_energy": evaluate_conserved_energy(data=data),
            "conservative_establishment":
                evaluate_conservative_establishment(data=data),
            "scalar_constrained_establishment":
                evaluate_scalar_constrained_establishment(data=data),
            "distributed_air_condensation":
                evaluate_distributed_air_condensation(data=data),
            "differential_scalar_spreading":
                evaluate_differential_scalar_spreading(data=data),
            "measured_0p3_m_s_coflow": evaluate_measured_coflow(data=data),
            "independent_mass_fraction_temperature":
                evaluate_independent_mass_fraction_temperature(data=data),
            "initial_entrainment_heating":
                evaluate_initial_entrainment_heating(data=data),
            "heated_two_scalar_jet":
                evaluate_heated_two_scalar_jet(data=data),
        },
    )


__all__ = [
    "DATA_PATH", "RamanVariant", "RamanMetrics", "RamanComparison", "VARIANTS",
    "load_data", "evaluate_variant", "evaluate_conserved_energy",
    "evaluate_conservative_establishment",
    "evaluate_scalar_constrained_establishment",
    "evaluate_distributed_air_condensation",
    "saturated_ambient_absolute_humidity",
    "evaluate_humid_air_frost",
    "evaluate_temperature_dependent_phase_enthalpy",
    "evaluate_argon_phase_completeness",
    "evaluate_argon_phase_completeness_audit",
    "evaluate_hydrogen_spin_isomer_enthalpy_audit",
    "evaluate_four_flux_phase_two_scalar_audit",
    "evaluate_blackbody_radiation_upper_bound",
    "evaluate_radiative_blackbody_ode",
    "radiation_optical_depth_sensitivity",
    "evaluate_humid_air_frost_upper_bound",
    "evaluate_differential_scalar_spreading",
    "evaluate_measured_coflow",
    "evaluate_independent_mass_fraction_temperature",
    "evaluate_temperature_dependent_enthalpy",
    "evaluate_initial_entrainment_heating",
    "evaluate_heated_two_scalar_jet", "compare_reported_case_coverage",
    "compare",
]
