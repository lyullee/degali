from dataclasses import replace

import numpy as np
import pytest

from degali.addons.atmospheric_surface_layer_mixing import (
    ATMOSPHERIC_SURFACE_LAYER_MIXING_SCHEMA,
    NeutralSurfaceLayerMixing,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule


def test_neutral_surface_layer_profile_and_solid_mask_are_exact():
    closure = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=0.205,
        roughness_length_m=0.005,
        source_id="TI021 workbook u-star and z0",
    )
    z = np.asarray((0.0, 0.5, 1.5))
    solid = np.zeros((3, 2, 2), dtype=bool)
    solid[1, 0, 1] = True

    viscosity = closure.eddy_viscosity_m2_s(z, solid)

    assert viscosity[:, 1, 0] == pytest.approx((0.0, 0.041, 0.123))
    assert viscosity[1, 0, 1] == 0.0
    assert closure.record()["schema"] == ATMOSPHERIC_SURFACE_LAYER_MIXING_SCHEMA
    assert "no concentration-derived coefficient selection" in closure.record()[
        "limitations"
    ]


@pytest.mark.parametrize(
    "kwargs",
    (
        {"friction_velocity_m_s": 0.0},
        {"roughness_length_m": -0.1},
        {"von_karman_constant": float("nan")},
        {"source_id": ""},
    ),
)
def test_neutral_surface_layer_rejects_invalid_evidence_inputs(kwargs):
    values = {
        "friction_velocity_m_s": 0.205,
        "roughness_length_m": 0.005,
        "source_id": "declared evidence",
    }
    values.update(kwargs)
    with pytest.raises(ValueError):
        NeutralSurfaceLayerMixing(**values)


def _coupled_config() -> ResolvedObstacleFlow3DConfig:
    mixture = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=2.0,
        z_min_m=0.0, z_max_m=2.0, species_mole_fraction=0.32,
        temperature_k=mixture.ambient_temperature_k,
    )
    surface_layer = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=0.205,
        roughness_length_m=0.005,
        source_id="TI021 workbook u-star and z0",
    )
    return ResolvedObstacleFlow3DConfig(
        domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
        wind_speed_m_s=0.2, kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=1, minimum_flow_steps=1, convergence_check_interval=1,
        lattice_inlet_speed=0.02, duration_s=0.001,
        scalar_time_step_s=0.001, scalar_diffusivity_m2_s=1.0e-5,
        thermal_diffusivity_m2_s=3.0e-6,
        source_schedule=SourceRateSchedule(
            (0.0, 0.001), (0.0, 0.0), source_id="initial-cloud-only",
        ),
        source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
        obstacles=(), ambient_temperature_k=mixture.ambient_temperature_k,
        ambient_pressure_pa=mixture.ambient_pressure_pa,
        ambient_air_density_kg_m3=mixture.ambient_density_kg_m3,
        gravity_m_s2=mixture.gravity_m_s2,
        thermal_buoyancy_coupling=True,
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud, binary_mixture_closure=mixture,
        thermodynamic_corrector_maximum_iterations=32,
        require_thermodynamic_corrector_convergence=True,
        variable_density_projection_linear_solver="cached_direct",
        subgrid_model="none", smagorinsky_constant=0.0,
        background_surface_layer_mixing=surface_layer,
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.74,
        scalar_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
        thermal_turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
    )


def test_surface_layer_mixing_enters_momentum_species_and_energy_diagnostics():
    result = solve_resolved_obstacle_flow_3d(_coupled_config())
    diagnostics = result.diagnostics

    assert diagnostics["background_surface_layer_mixing"]["source_id"] == (
        "TI021 workbook u-star and z0"
    )
    assert diagnostics["background_eddy_viscosity_minimum_m2_s"] == pytest.approx(
        0.041
    )
    assert diagnostics["background_eddy_viscosity_maximum_m2_s"] == pytest.approx(
        0.615
    )
    assert diagnostics["applied_scalar_diffusivity_minimum_m2_s"] > 1.0e-5
    assert diagnostics["applied_thermal_diffusivity_minimum_m2_s"] > 3.0e-6
    assert diagnostics["maximum_relative_mass_residual"] < 1.0e-10
    assert diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1.0e-8
    assert any(
        "no stability correction" in item
        for item in diagnostics["warnings"]
    )


def test_surface_layer_mixing_requires_transport_numbers_and_excludes_second_sgs():
    config = _coupled_config()
    with pytest.raises(ValueError, match="Schmidt and Prandtl"):
        replace(config, effective_schmidt_number=None)
    with pytest.raises(ValueError, match="cannot be combined"):
        replace(config, subgrid_model="wale", wale_constant=0.5)


def test_surface_layer_profile_rejects_non_boolean_or_mismatched_masks():
    closure = NeutralSurfaceLayerMixing(0.205, 0.005, "evidence")
    with pytest.raises(ValueError, match="boolean"):
        closure.eddy_viscosity_m2_s(
            np.asarray((0.5, 1.5)), np.zeros((2, 2, 2), dtype=float),
        )
    with pytest.raises(ValueError, match="matching z_m"):
        closure.eddy_viscosity_m2_s(
            np.asarray((0.5, 1.5)), np.zeros((3, 2, 2), dtype=bool),
        )
