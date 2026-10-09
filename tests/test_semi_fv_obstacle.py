import numpy as np
import pytest
from dataclasses import replace

from degali.addons.semi_fv_obstacle import (
    DistributedScalarSource,
    RectangularObstacle2D,
    SemiFVConfig,
    SemiFVReceptor,
    SourceRateSchedule,
    run_semi_fv_refinement_study,
    solve_semi_fv_obstacle,
)


def test_obstacle_dimensions_and_distance_are_explicit():
    obstacle = RectangularObstacle2D(distance_m=20.0, width_m=4.0, height_m=5.0)
    assert obstacle.x_min_m == pytest.approx(18.0)
    assert obstacle.x_max_m == pytest.approx(22.0)
    assert obstacle.z_max_m == pytest.approx(5.0)


def test_semi_fv_numeric_boundaries_reject_boolean_coercion():
    with pytest.raises(ValueError, match="not boolean"):
        RectangularObstacle2D(True, 1.0, 1.0)
    with pytest.raises(ValueError, match="not boolean"):
        SourceRateSchedule((0.0, 1.0), (True, 0.0))
    with pytest.raises(ValueError, match="not boolean"):
        SemiFVConfig(source_rate_kg_s=True).validate()
    with pytest.raises(ValueError, match="not boolean"):
        SemiFVConfig(source_height_m=True).validate()
    with pytest.raises(ValueError, match="positive integers"):
        run_semi_fv_refinement_study(
            SemiFVConfig(),
            receptors=(SemiFVReceptor("sensor", 1.0, 0.5),),
            refinement_factors=(1, True),
        )


def test_positive_gravitational_settling_moves_the_scalar_downward():
    def centre_height(settling_m_s: float) -> float:
        config = SemiFVConfig(
            length_m=4.0, height_m=4.0, nx=40, nz=80,
            duration_s=0.5, time_step_s=0.001, wind_speed_m_s=0.2,
            diffusivity_m2_s=1.0e-5,
            gravitational_settling_m_s=settling_m_s,
            source_rate_kg_s=0.2, source_height_m=3.0, source_sigma_m=0.15,
        )
        result = solve_semi_fv_obstacle(config)
        field = np.nan_to_num(result.concentration_kg_m3)
        dx = config.length_m / config.nx
        dz = config.height_m / config.nz
        mass = float(field.sum() * dx * dz)
        return float((field * result.z_m[:, None]).sum() * dx * dz / mass)

    lower = centre_height(0.5)
    nominal = centre_height(0.0)
    upper = centre_height(-0.5)
    assert lower < nominal < upper


def test_semi_fv_without_obstacle_conserves_mass():
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=20.0, height_m=10.0, nx=80, nz=30,
        duration_s=2.0, time_step_s=0.02, source_rate_kg_s=0.2,
    ))
    assert result.diagnostics.applicability == "accepted"
    assert result.diagnostics.maximum_mass_residual_kg < 1.0e-10
    assert result.diagnostics.final_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert result.diagnostics.source_mass_ledger_residual_kg == pytest.approx(0.0, abs=1.0e-15)
    assert dict(result.diagnostics.source_mass_schedule_residual_kg) == {
        "primary": pytest.approx(0.0, abs=1.0e-15),
    }
    assert result.diagnostics.maximum_source_mass_schedule_residual_kg == pytest.approx(0.0, abs=1.0e-15)
    assert result.diagnostics.x_advection_courant <= 1.0
    assert 2.0 * result.diagnostics.vertical_diffusion_number <= 1.0
    assert np.nanmax(result.concentration_kg_m3) > 0.0


def test_semi_fv_diagnostics_reject_nonphysical_typed_values():
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=0.5, time_step_s=0.01, source_rate_kg_s=0.2,
    ))

    with pytest.raises(ValueError, match="mass_domain_kg"):
        replace(result.diagnostics, mass_domain_kg=-1.0)
    with pytest.raises(ValueError, match="diverted_mass_fraction"):
        replace(result.diagnostics, diverted_mass_fraction=1.1)
    with pytest.raises(ValueError, match="source_mass_injected_kg labels"):
        replace(
            result.diagnostics,
            source_mass_injected_kg=(("primary", 0.1), ("primary", 0.1)),
        )


def test_semi_fv_reports_obstacle_contact_and_conserves_routed_inventory():
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=30.0, height_m=12.0, nx=120, nz=48,
        duration_s=3.0, time_step_s=0.02, source_rate_kg_s=0.2,
        obstacle=RectangularObstacle2D(distance_m=10.0, width_m=3.0, height_m=4.0),
    ))
    assert result.diagnostics.obstacle_contact
    assert result.diagnostics.applicability == "conditional"
    assert 0.0 <= result.diagnostics.diverted_mass_fraction <= 1.0
    assert result.diagnostics.maximum_mass_residual_kg < 1.0e-9
    assert result.diagnostics.final_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert result.diagnostics.source_mass_ledger_residual_kg == pytest.approx(0.0, abs=1.0e-15)
    assert np.isnan(result.concentration_kg_m3).any()
    assert result.maximum_concentration_kg_m3 > 0.0


def test_semi_fv_routes_conservatively_through_the_union_of_multiple_obstacles():
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=30.0, height_m=12.0, nx=120, nz=48,
        duration_s=3.0, time_step_s=0.02, source_rate_kg_s=0.2,
        obstacle=(
            RectangularObstacle2D(8.0, 2.0, 3.0, label="wall-a"),
            RectangularObstacle2D(14.0, 2.0, 4.0, label="wall-b"),
        ),
    ))

    assert result.diagnostics.obstacle_contact
    assert result.diagnostics.obstacle_count == 2
    assert result.diagnostics.maximum_mass_residual_kg < 1.0e-9
    assert result.diagnostics.diverted_mass_fraction > 0.0

    with pytest.raises(ValueError, match="labels must be unique"):
        SemiFVConfig(
            obstacle=(
                RectangularObstacle2D(8.0, 2.0, 3.0, label="same"),
                RectangularObstacle2D(14.0, 2.0, 4.0, label="same"),
            ),
        ).validate()


def test_obstacle_cannot_touch_boundary():
    with pytest.raises(ValueError):
        solve_semi_fv_obstacle(SemiFVConfig(
            length_m=10.0, height_m=5.0,
            obstacle=RectangularObstacle2D(distance_m=0.5, width_m=1.0, height_m=2.0),
        ))


@pytest.mark.parametrize(
    ("distance_m", "width_m"),
    ((0.2, 1.0), (9.8, 1.0)),
)
def test_obstacle_cannot_extend_past_local_inlet_or_outlet(
    distance_m, width_m,
):
    with pytest.raises(ValueError, match="strictly inside the local domain"):
        SemiFVConfig(
            length_m=10.0, height_m=5.0,
            obstacle=RectangularObstacle2D(
                distance_m=distance_m, width_m=width_m, height_m=2.0,
            ),
        ).validate()


def test_obstacle_must_fit_inside_the_declared_vertical_domain():
    with pytest.raises(ValueError, match="inside the local domain"):
        solve_semi_fv_obstacle(SemiFVConfig(
            length_m=10.0, height_m=5.0,
            obstacle=RectangularObstacle2D(
                distance_m=5.0, width_m=1.0, height_m=2.0, base_height_m=4.0,
            ),
        ))


def test_semi_fv_rejects_an_explicitly_unstable_time_step():
    with pytest.raises(ValueError, match="stability limit"):
        solve_semi_fv_obstacle(SemiFVConfig(
            length_m=10.0, height_m=5.0, nx=40, nz=40,
            time_step_s=0.05, duration_s=1.0, diffusivity_m2_s=1.0,
        ))


def test_semi_fv_integrates_a_declared_atmospheric_source_schedule_exactly():
    schedule = SourceRateSchedule(
        time_s=(0.0, 0.1, 0.2), rate_kg_s=(2.0, 0.0, 0.0),
        source_id="direct-vapour-record",
    )
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=0.2, time_step_s=0.01, diffusivity_m2_s=0.05,
        source_rate_kg_s=0.0, source_schedule=schedule,
    ))

    assert schedule.released_mass_kg == pytest.approx(0.2)
    assert result.diagnostics.mass_injected_kg == pytest.approx(0.2)
    assert dict(result.diagnostics.source_mass_injected_kg) == {
        "primary": pytest.approx(0.2),
    }
    assert result.diagnostics.source_mode == "scheduled"
    assert result.diagnostics.final_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)


def test_linear_source_schedule_integrates_node_values_and_conserves_transport_mass():
    schedule = SourceRateSchedule(
        time_s=(0.0, 1.0, 2.0), rate_kg_s=(1.0, 3.0, 0.0),
        source_id="linear-history", rate_operator="linear",
    )
    assert schedule.released_mass_kg == pytest.approx(3.5)
    assert schedule.mass_between(0.25, 1.25) == pytest.approx(2.34375)
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=2.0, time_step_s=0.01, diffusivity_m2_s=0.05,
        source_rate_kg_s=0.0, source_schedule=schedule,
    ))

    assert result.diagnostics.mass_injected_kg == pytest.approx(3.5)
    assert result.diagnostics.final_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)


def test_source_schedule_rejects_an_unknown_rate_operator():
    with pytest.raises(ValueError, match="rate_operator"):
        SourceRateSchedule((0.0, 1.0), (1.0, 0.0), rate_operator="cubic")


def test_semi_fv_conserves_a_declared_internal_distributed_scalar_source():
    source = DistributedScalarSource(
        "in-flight-bin-01", x_m=1.0, height_m=0.8, sigma_m=0.2,
        schedule=SourceRateSchedule(
            (0.0, 0.1, 0.2), (0.0, 2.0, 0.0), source_id="droplet-ledger-bin-01",
        ),
    )
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=0.3, time_step_s=0.01, diffusivity_m2_s=0.05,
        source_rate_kg_s=0.0, distributed_sources=(source,),
    ))

    assert result.diagnostics.source_mode == "distributed"
    assert result.diagnostics.distributed_source_count == 1
    assert result.diagnostics.mass_injected_kg == pytest.approx(0.2)
    assert dict(result.diagnostics.source_mass_injected_kg) == {
        "primary": pytest.approx(0.0),
        "distributed:in-flight-bin-01": pytest.approx(0.2),
    }
    assert dict(result.diagnostics.source_mass_schedule_residual_kg) == {
        "primary": pytest.approx(0.0, abs=1.0e-15),
        "distributed:in-flight-bin-01": pytest.approx(0.0, abs=1.0e-15),
    }
    assert result.diagnostics.maximum_mass_residual_kg < 1.0e-10
    assert result.diagnostics.final_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)

    with pytest.raises(ValueError, match="inside a declared solid"):
        solve_semi_fv_obstacle(SemiFVConfig(
            length_m=4.0, height_m=2.0, nx=20, nz=10,
            duration_s=0.3, time_step_s=0.01, diffusivity_m2_s=0.05,
            source_rate_kg_s=0.0,
            obstacle=RectangularObstacle2D(1.0, 0.4, 1.0, label="wall"),
            distributed_sources=(source,),
        ))


def test_semi_fv_source_schedule_must_have_an_unambiguous_duration_and_rate():
    schedule = SourceRateSchedule((0.0, 1.0), (1.0, 0.0))
    with pytest.raises(ValueError, match="cannot exceed"):
        solve_semi_fv_obstacle(SemiFVConfig(
            duration_s=0.5, source_rate_kg_s=0.0, source_schedule=schedule,
        ))
    with pytest.raises(ValueError, match="must be zero"):
        SemiFVConfig(source_schedule=schedule).validate()


def test_semi_fv_rejects_a_nonzero_primary_schedule_endpoint():
    schedule = SourceRateSchedule((0.0, 1.0), (1.0, 0.2))
    assert not schedule.has_zero_endpoint
    with pytest.raises(ValueError, match="source schedule final endpoint rate must be zero"):
        SemiFVConfig(source_rate_kg_s=0.0, source_schedule=schedule).validate()


def test_semi_fv_rejects_a_nonzero_distributed_schedule_endpoint():
    source = DistributedScalarSource(
        "internal", x_m=1.0, height_m=0.5, sigma_m=0.2,
        schedule=SourceRateSchedule((0.0, 1.0), (0.1, 0.01)),
    )
    with pytest.raises(ValueError, match="distributed source 'internal' schedule final endpoint rate must be zero"):
        SemiFVConfig(source_rate_kg_s=0.0, distributed_sources=(source,)).validate()


def test_semi_fv_stops_a_short_source_schedule_but_continues_its_inventory():
    schedule = SourceRateSchedule(
        (0.0, 0.1), (2.0, 0.0), source_id="finite-release",
    )
    result = solve_semi_fv_obstacle(SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=0.3, time_step_s=0.01, diffusivity_m2_s=0.05,
        source_rate_kg_s=0.0, source_schedule=schedule,
    ))

    assert result.diagnostics.mass_injected_kg == pytest.approx(0.2)
    assert result.receptor_traces == ()


def test_semi_fv_refinement_study_reports_receptor_peak_and_dose_changes():
    config = SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=10, nz=10,
        duration_s=0.1, time_step_s=0.005, diffusivity_m2_s=0.05,
        source_rate_kg_s=0.2, source_height_m=0.5, source_sigma_m=0.2,
    )
    study = run_semi_fv_refinement_study(
        config, receptors=(SemiFVReceptor("sensor", 0.6, 0.5),),
        refinement_factors=(1, 2), relative_tolerance=0.05,
    )

    assert [factor for factor, _ in study.cases] == [1, 2]
    assert study.estimated_cell_steps == 34_000
    assert len(study.receptor_changes) == 1
    change = study.receptor_changes[0]
    assert change.coarse_peak_kg_m3 > 0.0
    assert change.fine_peak_kg_m3 > 0.0
    assert change.coarse_dose_kg_s_m3 > 0.0
    assert 0.0 < change.peak_relative_change <= 1.0
    assert 0.0 < change.dose_relative_change <= 1.0
    assert not study.converged
    assert "did not meet" in study.warnings[0]


def test_semi_fv_refinement_study_refuses_an_unbounded_work_request():
    with pytest.raises(ValueError, match="max_cell_steps"):
        run_semi_fv_refinement_study(
            SemiFVConfig(length_m=4.0, height_m=2.0, nx=10, nz=10, duration_s=0.1),
            receptors=(SemiFVReceptor("sensor", 0.6, 0.5),),
            max_cell_steps=1,
        )


def test_semi_fv_records_true_receptor_history_and_rejects_solid_receptor():
    config = SemiFVConfig(
        length_m=10.0, height_m=5.0, nx=80, nz=40,
        duration_s=1.0, time_step_s=0.005, source_rate_kg_s=0.2,
    )
    result = solve_semi_fv_obstacle(
        config, receptors=(SemiFVReceptor("sensor", 1.0, 0.5),)
    )
    trace = result.receptor_traces[0]
    assert trace.time_s[0] == pytest.approx(0.0)
    assert trace.time_s[-1] == pytest.approx(config.duration_s)
    assert len(trace.time_s) == len(trace.concentration_kg_m3)
    assert np.max(trace.concentration_kg_m3) > 0.0

    with pytest.raises(ValueError, match="inside declared solid"):
        solve_semi_fv_obstacle(
            SemiFVConfig(
                obstacle=RectangularObstacle2D(10.0, 2.0, 5.0),
            ),
            receptors=(SemiFVReceptor("in-wall", 10.0, 0.5),),
        )
