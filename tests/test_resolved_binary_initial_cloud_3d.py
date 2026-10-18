from dataclasses import replace

import numpy as np
import pytest

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DRequest,
)
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
)
from degali.addons.resolved_obstacle_flow_3d import (
    RESOLVED_OBSTACLE_FLOW_3D_SCHEMA,
    ResolvedObstacleFlow3DConfig,
    read_resolved_obstacle_flow_3d_restart,
    rebind_restart_for_corrector_tolerance_audit,
    solve_resolved_obstacle_flow_3d,
    write_resolved_obstacle_flow_3d_restart,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule


def test_resolved_solver_advances_full_binary_dense_initial_cloud():
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0),
        flat_to_flat_m=4.0,
        z_min_m=0.0,
        z_max_m=3.0,
        species_mole_fraction=0.32,
        temperature_k=closure.ambient_temperature_k,
        label="small TI021 analogue",
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(8.0, 8.0, 8.0),
        cells=(8, 8, 8),
        wind_speed_m_s=0.2,
        kinematic_viscosity_m2_s=0.1,
        flow_steps=5,
        minimum_flow_steps=5,
        convergence_check_interval=1,
        lattice_inlet_speed=0.02,
        duration_s=0.01,
        scalar_time_step_s=0.01,
        scalar_diffusivity_m2_s=0.01,
        thermal_diffusivity_m2_s=0.01,
        source_schedule=SourceRateSchedule(
            (0.0, 0.01), (0.0, 0.0), source_id="initial-cloud-only",
        ),
        source_position_m=(1.0, 1.0, 1.0),
        source_sigma_m=0.5,
        obstacles=(),
        ambient_temperature_k=closure.ambient_temperature_k,
        ambient_pressure_pa=closure.ambient_pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        gravity_m_s2=closure.gravity_m_s2,
        thermal_buoyancy_coupling=True,
        active_momentum_solver="finite_volume_projection",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud,
        binary_mixture_closure=closure,
        thermodynamic_corrector_maximum_iterations=32,
        thermodynamic_corrector_minimum_iterations=2,
        thermodynamic_corrector_density_tolerance=1e-6,
        thermodynamic_corrector_velocity_tolerance=1e-6,
        require_thermodynamic_corrector_convergence=True,
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.7,
        stratification_damping_alpha=5.0,
        stratification_shear_regularization_s_inv=0.041,
        store_scalar_fields=True,
        subgrid_model="none",
        smagorinsky_constant=0.0,
    )

    baseline = solve_resolved_obstacle_flow_3d(replace(
        config,
        stratification_damping_alpha=None,
        stratification_shear_regularization_s_inv=None,
    ))
    handoff_frames = []
    result = solve_resolved_obstacle_flow_3d(
        config,
        handoff_request=BinaryMixtureHandoffPlane3DRequest(
            x_face_m=6.0,
            sample_interval_s=0.01,
        ),
        handoff_observer=handoff_frames.append,
    )

    assert result.schema == RESOLVED_OBSTACLE_FLOW_3D_SCHEMA
    assert result.species_name == "R12"
    assert result.species_mole_fraction is not None
    assert result.maximum_species_mole_fraction is not None
    assert result.hydrogen_mole_fraction is None
    assert np.nanmax(result.maximum_concentration_vol_percent) == pytest.approx(
        32.0, rel=0.01,
    )
    assert result.domain_mass_kg[0] > 0.0
    assert result.mass_residual_kg[0] == 0.0
    assert result.diagnostics["initial_cloud_exact_geometry_volume_m3"] == (
        pytest.approx(cloud.exact_volume_m3)
    )
    assert abs(result.diagnostics["initial_cloud_relative_geometry_volume_error"]) < 1e-12
    assert result.diagnostics["stratification_damping_alpha"] == 5.0
    assert result.diagnostics["maximum_gradient_richardson_number"] > 0.0
    assert result.diagnostics["maximum_relative_mass_residual"] < 1e-10
    assert result.diagnostics["maximum_relative_energy_residual"] < 1e-10
    assert np.max(np.abs(
        result.concentration_kg_m3 - baseline.concentration_kg_m3
    )) > 1e-12
    assert result.diagnostics["maximum_relative_mixture_density_departure"] > 1.0
    assert result.diagnostics["active_scalar_out_of_range_cell_steps"] == 0
    assert result.diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
    assert result.diagnostics["maximum_relative_mass_residual"] < 1e-10
    assert result.diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1e-8
    assert [frame.time_s for frame in handoff_frames] == [0.0, 0.01]
    assert all(
        frame.species_mass_flux_kg_m2_s.shape == (8, 8)
        for frame in handoff_frames
    )
    assert result.diagnostics["binary_handoff_recording_enabled"]
    assert result.diagnostics["binary_handoff_x_face_m"] == 6.0
    assert result.diagnostics["binary_handoff_x_face_index"] == 6
    assert result.diagnostics["binary_handoff_frame_count"] == 2


def test_binary_fv_restart_matches_one_continuous_run(tmp_path):
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=4.0,
        z_min_m=0.0, z_max_m=3.0,
        species_mole_fraction=0.32,
        temperature_k=closure.ambient_temperature_k,
    )
    common = dict(
        domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
        wind_speed_m_s=0.2, kinematic_viscosity_m2_s=0.1,
        flow_steps=5, minimum_flow_steps=5,
        convergence_check_interval=1, lattice_inlet_speed=0.02,
        scalar_time_step_s=0.01, scalar_diffusivity_m2_s=0.01,
        thermal_diffusivity_m2_s=0.01,
        source_schedule=SourceRateSchedule(
            (0.0, 0.02), (0.0, 0.0), source_id="initial-cloud-only",
        ),
        source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
        obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
        ambient_pressure_pa=closure.ambient_pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        gravity_m_s2=closure.gravity_m_s2,
        thermal_buoyancy_coupling=True,
        active_momentum_solver="finite_volume_projection",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud,
        binary_mixture_closure=closure,
        thermodynamic_corrector_maximum_iterations=32,
        thermodynamic_corrector_minimum_iterations=2,
        thermodynamic_corrector_density_tolerance=1e-6,
        thermodynamic_corrector_velocity_tolerance=1e-6,
        require_thermodynamic_corrector_convergence=True,
        store_scalar_fields=False, subgrid_model="none",
        smagorinsky_constant=0.0,
    )
    continuous_config = ResolvedObstacleFlow3DConfig(
        duration_s=0.02, **common,
    )
    continuous = solve_resolved_obstacle_flow_3d(continuous_config)
    first = solve_resolved_obstacle_flow_3d(
        continuous_config, stop_time_s=0.01,
    )
    assert first.restart_state is not None
    checkpoint = write_resolved_obstacle_flow_3d_restart(
        tmp_path / "restart.npz", first.restart_state,
    )
    restored = read_resolved_obstacle_flow_3d_restart(checkpoint)
    assert restored.configuration_sha256 == (
        first.restart_state.configuration_sha256
    )
    assert restored.time_s == pytest.approx(0.01)
    resumed = solve_resolved_obstacle_flow_3d(
        continuous_config, restart_state=restored,
        stop_time_s=0.02,
    )

    assert resumed.time_s == pytest.approx((0.01, 0.02))
    assert resumed.diagnostics["restart_applied"] is True
    assert resumed.diagnostics["restart_time_s"] == pytest.approx(0.01)
    assert resumed.concentration_kg_m3 == pytest.approx(
        continuous.concentration_kg_m3, rel=1e-11, abs=1e-13,
    )
    assert resumed.sensible_enthalpy_density_j_m3 == pytest.approx(
        continuous.sensible_enthalpy_density_j_m3, rel=1e-11, abs=1e-10,
    )
    assert resumed.velocity_u_m_s == pytest.approx(
        continuous.velocity_u_m_s, rel=1e-11, abs=1e-12,
    )
    assert resumed.velocity_v_m_s == pytest.approx(
        continuous.velocity_v_m_s, rel=1e-11, abs=1e-12,
    )
    assert resumed.velocity_w_m_s == pytest.approx(
        continuous.velocity_w_m_s, rel=1e-11, abs=1e-12,
    )
    assert resumed.maximum_species_mole_fraction == pytest.approx(
        continuous.maximum_species_mole_fraction, rel=1e-11, abs=1e-13,
    )
    with pytest.raises(ValueError, match="configuration SHA-256"):
        solve_resolved_obstacle_flow_3d(
            ResolvedObstacleFlow3DConfig(
                **{
                    **continuous_config.__dict__,
                    "kinematic_viscosity_m2_s": 0.11,
                }
            ),
            restart_state=restored,
            stop_time_s=0.02,
        )

    strict_config = ResolvedObstacleFlow3DConfig(
        **{
            **continuous_config.__dict__,
            "thermodynamic_corrector_density_tolerance": 1.0e-7,
            "thermodynamic_corrector_velocity_tolerance": 1.0e-7,
        }
    )
    rebound = rebind_restart_for_corrector_tolerance_audit(
        restored, source_config=continuous_config, target_config=strict_config,
    )
    strict_resumed = solve_resolved_obstacle_flow_3d(
        strict_config, restart_state=rebound, stop_time_s=0.02,
    )
    assert strict_resumed.diagnostics["restart_applied"] is True
    assert rebound.time_s == restored.time_s
    assert rebound.concentration_kg_m3 == pytest.approx(
        restored.concentration_kg_m3,
    )
    with pytest.raises(ValueError, match="unexpected changes"):
        rebind_restart_for_corrector_tolerance_audit(
            restored,
            source_config=continuous_config,
            target_config=ResolvedObstacleFlow3DConfig(
                **{
                    **strict_config.__dict__,
                    "kinematic_viscosity_m2_s": 0.11,
                }
            ),
        )


def test_quiescent_binary_cloud_removes_its_initial_carrier_momentum():
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=4.0,
        z_min_m=0.0, z_max_m=3.0,
        species_mole_fraction=0.32,
        temperature_k=closure.ambient_temperature_k,
    )

    def run(mode):
        frames = []
        config = ResolvedObstacleFlow3DConfig(
            domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
            wind_speed_m_s=0.2, kinematic_viscosity_m2_s=0.1,
            flow_steps=5, minimum_flow_steps=5,
            convergence_check_interval=1, lattice_inlet_speed=0.02,
            duration_s=0.01, scalar_time_step_s=0.01,
            scalar_diffusivity_m2_s=0.01,
            thermal_diffusivity_m2_s=0.01,
            source_schedule=SourceRateSchedule(
                (0.0, 0.01), (0.0, 0.0), source_id="initial-cloud-only",
            ),
            source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
            obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
            ambient_pressure_pa=closure.ambient_pressure_pa,
            ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
            gravity_m_s2=closure.gravity_m_s2,
            thermal_buoyancy_coupling=True,
            source_temperature_k=closure.ambient_temperature_k,
            active_momentum_solver="finite_volume_projection",
            low_mach_mass_flux_correction=True,
            initial_binary_cloud=cloud,
            binary_mixture_closure=closure,
            initial_binary_cloud_momentum=mode,
            thermodynamic_corrector_maximum_iterations=32,
            thermodynamic_corrector_minimum_iterations=2,
            thermodynamic_corrector_density_tolerance=1e-6,
            thermodynamic_corrector_velocity_tolerance=1e-6,
            require_thermodynamic_corrector_convergence=True,
            store_scalar_fields=False, subgrid_model="none",
            smagorinsky_constant=0.0,
        )
        result = solve_resolved_obstacle_flow_3d(
            config,
            handoff_request=BinaryMixtureHandoffPlane3DRequest(
                x_face_m=4.0, sample_interval_s=0.01,
            ),
            handoff_observer=frames.append,
        )
        return result, frames[0]

    carried, carried_frame = run("ambient_carried")
    quiescent, quiescent_frame = run("quiescent")

    assert carried.diagnostics["initial_binary_cloud_momentum"] == (
        "ambient_carried"
    )
    assert quiescent.diagnostics["initial_binary_cloud_momentum"] == "quiescent"
    assert quiescent.diagnostics["initial_cloud_maximum_speed_m_s"] < (
        carried.diagnostics["initial_cloud_maximum_speed_m_s"]
    )
    assert np.max(np.abs(quiescent_frame.species_mass_flux_kg_m2_s)) < (
        np.max(np.abs(carried_frame.species_mass_flux_kg_m2_s))
    )


def test_binary_initial_cloud_configuration_accepts_a_continuing_source():
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(3.0, 4.0), flat_to_flat_m=2.0,
        z_min_m=0.0, z_max_m=2.0,
        species_mole_fraction=0.32,
        temperature_k=closure.ambient_temperature_k,
    )
    config = ResolvedObstacleFlow3DConfig(
            domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
            wind_speed_m_s=0.2, kinematic_viscosity_m2_s=0.1,
            flow_steps=5, minimum_flow_steps=5,
            convergence_check_interval=1, lattice_inlet_speed=0.02,
            duration_s=0.01, scalar_time_step_s=0.01,
            scalar_diffusivity_m2_s=0.01,
            source_schedule=SourceRateSchedule(
                (0.0, 0.005, 0.01), (0.01, 0.0, 0.0),
            ),
            source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
            obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
            ambient_pressure_pa=closure.ambient_pressure_pa,
            ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
            gravity_m_s2=closure.gravity_m_s2,
            thermal_buoyancy_coupling=True,
            source_temperature_k=closure.ambient_temperature_k,
            active_momentum_solver="finite_volume_projection",
            initial_binary_cloud=cloud,
            binary_mixture_closure=closure,
        )
    assert config.source_schedule.released_mass_kg == pytest.approx(5.0e-5)


def test_binary_mixture_source_release_conserves_species_and_enthalpy():
    closure = BinaryIdealGasMixtureClosure(
        species_name="ethylene",
        species_molar_mass_kg_mol=28.05316e-3,
        species_heat_capacity_j_kg_k=1_500.0,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=101_325.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=2.0,
        z_min_m=0.0, z_max_m=2.0,
        species_mole_fraction=0.0,
        temperature_k=closure.ambient_temperature_k,
        label="ambient carrier initialization",
    )
    total_release_kg = 1.0e-6
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
        wind_speed_m_s=0.2, kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=5, minimum_flow_steps=5, convergence_check_interval=1,
        lattice_inlet_speed=0.03, duration_s=0.01,
        scalar_time_step_s=0.01, scalar_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=1.0e-20,
        source_schedule=SourceRateSchedule(
            (0.0, 0.005, 0.01), (total_release_kg / 0.005, 0.0, 0.0),
            source_id="ethylene_point_source",
        ),
        source_total_mass_to_species_mass_ratio=20.0,
        source_position_m=(2.0, 4.0, 1.0), source_sigma_m=0.5,
        obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
        ambient_pressure_pa=closure.ambient_pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        gravity_m_s2=closure.gravity_m_s2,
        thermal_buoyancy_coupling=True, source_temperature_k=303.55,
        inlet_temperature_profile_z_m=(0.0, 8.0),
        inlet_temperature_profile_k=(300.0, 291.15),
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud, binary_mixture_closure=closure,
        transient_average_start_s=0.005,
        thermodynamic_corrector_maximum_iterations=40,
        thermodynamic_corrector_density_tolerance=1.0e-7,
        thermodynamic_corrector_velocity_tolerance=1.0e-7,
        require_thermodynamic_corrector_convergence=True,
        subgrid_model="none", store_scalar_fields=True,
    )

    result = solve_resolved_obstacle_flow_3d(config)

    assert result.species_name == "ethylene"
    assert result.diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
    assert result.diagnostics["maximum_relative_mass_residual"] < 1.0e-8
    assert result.diagnostics["maximum_relative_energy_residual"] < 1.0e-8
    assert result.injected_mass_kg[-1] == pytest.approx(total_release_kg)
    assert result.diagnostics["source_total_mass_to_species_mass_ratio"] == 20.0
    assert np.sum(result.concentration_kg_m3) > 0.0
    assert result.maximum_species_mole_fraction is not None
    assert np.max(result.maximum_species_mole_fraction) > 0.0
    assert result.temperature_fields_k is not None
    assert np.max(result.temperature_fields_k[0]) > np.min(
        result.temperature_fields_k[0]
    )
    assert result.diagnostics["inlet_temperature_boundary_model"] == (
        "conservative_x_min_sensible_enthalpy_inflow; linearly interpolated in height with endpoint hold"
    )
    assert result.time_average_sample_count == 1
    assert result.time_averaged_velocity_u_m_s == pytest.approx(
        result.velocity_u_m_s
    )


def test_projected_inlet_initializes_fv_at_physical_molecular_viscosity():
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(4.0, 4.0), flat_to_flat_m=2.0,
        z_min_m=0.0, z_max_m=2.0,
        species_mole_fraction=0.32,
        temperature_k=closure.ambient_temperature_k,
    )
    config = ResolvedObstacleFlow3DConfig(
        domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
        wind_speed_m_s=3.9, kinematic_viscosity_m2_s=1.5e-5,
        flow_steps=5, minimum_flow_steps=5,
        convergence_check_interval=1, lattice_inlet_speed=0.03,
        duration_s=0.001, scalar_time_step_s=0.001,
        scalar_diffusivity_m2_s=1.0e-5,
        thermal_diffusivity_m2_s=3.0e-6,
        source_schedule=SourceRateSchedule(
            (0.0, 0.001), (0.0, 0.0), source_id="initial-cloud-only",
        ),
        source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
        obstacles=(), ambient_temperature_k=closure.ambient_temperature_k,
        ambient_pressure_pa=closure.ambient_pressure_pa,
        ambient_air_density_kg_m3=closure.ambient_density_kg_m3,
        gravity_m_s2=closure.gravity_m_s2,
        thermal_buoyancy_coupling=True,
        active_momentum_solver="finite_volume_projection",
        carrier_flow_initialization="projected_inlet",
        low_mach_mass_flux_correction=True,
        initial_binary_cloud=cloud,
        binary_mixture_closure=closure,
        thermodynamic_corrector_maximum_iterations=32,
        require_thermodynamic_corrector_convergence=True,
        subgrid_model="none",
    )

    result = solve_resolved_obstacle_flow_3d(config)

    assert result.diagnostics["carrier_flow_initialization"] == "projected_inlet"
    assert result.diagnostics["flow_steps_completed"] == 0
    assert result.diagnostics["base_lattice_relaxation_time"] is None
    assert result.diagnostics["maximum_relative_mass_residual"] < 1e-10
    assert result.diagnostics[
        "maximum_variable_density_continuity_residual_kg_m3_s"
    ] < 1e-8


def test_projected_inlet_rejects_d3q19_active_momentum():
    with pytest.raises(
        ValueError, match="requires the finite_volume_projection",
    ):
        ResolvedObstacleFlow3DConfig(
            domain_m=(8.0, 8.0, 8.0), cells=(8, 8, 8),
            wind_speed_m_s=0.2, kinematic_viscosity_m2_s=0.1,
            flow_steps=5, minimum_flow_steps=5,
            convergence_check_interval=1, lattice_inlet_speed=0.02,
            duration_s=0.01, scalar_time_step_s=0.01,
            scalar_diffusivity_m2_s=0.01,
            source_schedule=SourceRateSchedule(
                (0.0, 0.01), (0.0, 0.0),
            ),
            source_position_m=(1.0, 1.0, 1.0), source_sigma_m=0.5,
            obstacles=(), carrier_flow_initialization="projected_inlet",
        )
