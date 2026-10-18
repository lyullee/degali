import numpy as np
import pytest

from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.transient_dense_gas_3d import (
    DenseGas3DReceptor,
    TransientDenseGas3DConfig,
    WindHistory,
    solve_transient_dense_gas_3d,
)
from degali.addons.source_state_ledger import SourceState, SourceStateLedger


def _config(*, source_density=10.0, wind_end=0.4, store_fields=True):
    schedule = SourceRateSchedule(
        (0.0, 0.2, 0.4), (1.0, 1.0, 0.0), source_id="event-3d",
    )
    wind = WindHistory(
        (0.0, 0.2, wind_end), (1.0, 1.0, 1.0), (270.0, 270.0, 180.0),
    )
    return TransientDenseGas3DConfig(
        domain_m=(12.0, 12.0, 8.0), cells=(12, 12, 8),
        duration_s=0.4, time_step_s=0.02,
        source_schedule=schedule, wind_history=wind,
        source_position_m=(2.0, 6.0, 4.0), source_sigma_m=0.8,
        source_density_kg_m3=source_density, source_h2_mass_fraction=0.9,
        ambient_density_kg_m3=1.2, diffusivity_m2_s=0.05,
        store_fields=store_fields,
    )


def test_source_and_wind_histories_are_coupled_on_the_common_clock():
    result = solve_transient_dense_gas_3d(_config(), receptors=(
        DenseGas3DReceptor("sensor", (3.0, 6.0, 4.0)),
    ))

    assert result.wind_vector_m_s[0] == pytest.approx((1.0, 0.0, 0.0))
    assert result.wind_vector_m_s[-1] == pytest.approx((0.0, 1.0, 0.0))
    assert result.receptor_traces_kg_m3["sensor"].shape == result.time_s.shape
    assert result.fields_kg_m3 is not None
    assert result.mean_velocity_m_s is not None
    assert result.mean_velocity_m_s[-1, 1] > result.mean_velocity_m_s[0, 1]


def test_three_dimensional_transport_closes_hydrogen_mass():
    result = solve_transient_dense_gas_3d(_config(), receptors=())

    assert result.injected_mass_kg[-1] == pytest.approx(0.4)
    assert result.total_mass_kg[-1] + result.outflow_mass_kg[-1] == pytest.approx(0.4)
    assert result.maximum_mass_residual_kg < 1.0e-10
    assert result.diagnostics["solver_scope"].startswith("3-D explicit")


def test_dense_source_generates_downward_buoyant_transport():
    result = solve_transient_dense_gas_3d(_config(source_density=10.0))
    assert float(np.max(result.maximum_buoyant_speed_m_s)) > 0.0
    initial = result.fields_kg_m3[1]
    final = result.fields_kg_m3[-1]
    z = result.z_m[:, None, None]
    initial_centroid = float(np.sum(initial * z) / np.sum(initial))
    final_centroid = float(np.sum(final * z) / np.sum(final))
    assert final_centroid < initial_centroid
    assert result.mean_velocity_m_s[-1, 2] < 0.0


def test_prognostic_velocity_can_relax_to_changing_wind_and_store_3d_velocity():
    config = _config(store_fields=False)
    config = TransientDenseGas3DConfig(
        **{
            **config.__dict__,
            "wind_relaxation_time_s": 0.2,
            "store_velocity_fields": True,
        }
    )
    result = solve_transient_dense_gas_3d(config)
    assert result.velocity_fields_m_s is not None
    assert result.velocity_fields_m_s.shape == (
        result.time_s.size, 3, config.cells[2], config.cells[1], config.cells[0],
    )
    # At the final time the wind target is northward, but finite relaxation
    # retains a non-zero memory of the earlier eastward forcing.
    assert result.mean_velocity_m_s[-1, 0] > 0.0
    assert result.mean_velocity_m_s[-1, 1] > 0.0


def test_wind_history_must_cover_entire_run():
    with pytest.raises(ValueError, match="cover the complete solver duration"):
        solve_transient_dense_gas_3d(_config(wind_end=0.3))


def test_source_state_ledger_handoff_runs_the_3d_solver():
    state = lambda time, rate: SourceState(
        time_s=time, h2_rate_kg_s=rate, h2_mass_fraction=0.9,
        temperature_k=90.0, density_kg_m3=10.0, area_m2=0.01,
    )
    ledger = SourceStateLedger(
        substance="hydrogen", stage="post_flash_atmospheric",
        states=(state(0.0, 1.0), state(0.2, 1.0), state(0.4, 0.0)),
        duration_s=0.4, metadata={"source_id": "ledger-event"},
    )
    config = _config()
    config = TransientDenseGas3DConfig(
        **{**config.__dict__, "source_schedule": ledger.to_source_rate_schedule()}
    )
    result = solve_transient_dense_gas_3d(config)
    assert result.diagnostics["source_id"] == "ledger-event"
    assert result.injected_mass_kg[-1] == pytest.approx(0.4)


def test_closed_axes_retain_scalar_inventory():
    base = _config(store_fields=False)
    config = TransientDenseGas3DConfig(
        **{**base.__dict__, "boundary_modes": ("closed", "closed", "closed")}
    )
    result = solve_transient_dense_gas_3d(config)

    assert result.outflow_mass_kg[-1] == pytest.approx(0.0)
    assert result.total_mass_kg[-1] == pytest.approx(result.injected_mass_kg[-1])
    assert result.diagnostics["boundary_modes"] == ["closed", "closed", "closed"]


def test_internal_solid_mask_blocks_flux_without_mass_loss():
    base = _config(store_fields=False)
    solid = np.zeros((base.cells[2], base.cells[1], base.cells[0]), dtype=bool)
    solid[:, 4:8, 5] = True
    config = TransientDenseGas3DConfig(
        **{**base.__dict__, "solid_mask_zyx": solid}
    )
    result = solve_transient_dense_gas_3d(
        config,
        receptors=(DenseGas3DReceptor("downstream", (8.0, 6.0, 4.0)),),
    )
    assert result.diagnostics["solid_cell_count"] == int(np.count_nonzero(solid))
    assert result.maximum_mass_residual_kg < 1.0e-10


def test_wall_receptor_samples_nearest_resolved_cell_centre():
    base = _config(store_fields=False)
    config = TransientDenseGas3DConfig(
        **{
            **base.__dict__,
            "source_position_m": (12.0, 12.0, 8.0),
            "boundary_modes": ("closed", "closed", "closed"),
        }
    )
    result = solve_transient_dense_gas_3d(
        config,
        receptors=(DenseGas3DReceptor("wall", (12.0, 12.0, 8.0)),),
    )

    assert float(np.max(result.receptor_traces_kg_m3["wall"])) > 0.0


def test_receptor_outside_physical_domain_is_rejected():
    with pytest.raises(ValueError, match="must lie inside the domain"):
        solve_transient_dense_gas_3d(
            _config(),
            receptors=(DenseGas3DReceptor("outside", (12.1, 6.0, 4.0)),),
        )


def test_diffusion_spreads_a_pulse_and_conserves_mass_next_to_obstacles():
    from degali.addons.transient_dense_gas_3d import _diffuse_neumann
    field = np.zeros((5, 5, 5))
    field[2, 2, 2] = 1.0
    solid = np.zeros_like(field, dtype=bool)
    solid[2, 2, 3] = True
    after = _diffuse_neumann(field, .1, .1, (1., 1., 1.), solid)
    assert after[2, 2, 1] > 0.0
    assert after[2, 1, 2] > 0.0
    assert after[1, 2, 2] > 0.0
    assert after[2, 2, 3] == 0.0
    assert after.sum() == pytest.approx(field.sum(), abs=1e-14)
    assert np.min(after) >= 0.0


@pytest.mark.parametrize("prognostic", [True, False])
def test_impermeable_ground_retains_falling_gas_and_preserves_open_top(prognostic):
    from dataclasses import replace
    base = replace(_config(), source_position_m=(6., 6., .2),
                   prognostic_velocity=prognostic)
    opened = solve_transient_dense_gas_3d(base)
    closed = solve_transient_dense_gas_3d(replace(base, impermeable_ground=True))
    assert opened.diagnostics['outflow_species_by_axis_kg'][2] > 0.0
    assert closed.diagnostics['outflow_species_by_axis_kg'][2] == pytest.approx(0., abs=1e-14)
    assert closed.total_mass_kg[-1] > opened.total_mass_kg[-1]
    assert closed.maximum_mass_residual_kg < 1e-12
    light = solve_transient_dense_gas_3d(replace(base, source_density_kg_m3=.1,
        source_position_m=(6., 6., 7.8), impermeable_ground=True))
    assert light.diagnostics['outflow_species_by_axis_kg'][2] > 0.0


def test_initial_species_inventory_includes_t0_trace_and_mass_ledger():
    from dataclasses import replace
    base = _config()
    initial = np.zeros((8, 12, 12))
    initial[3, 5, 2] = .2
    base = replace(base, boundary_modes=('closed', 'closed', 'closed'),
                   source_schedule=SourceRateSchedule((0., .4), (0., 0.)))
    result = solve_transient_dense_gas_3d(base, initial_species_kg_m3=initial,
        receptors=(DenseGas3DReceptor('initial', (2.5, 5.5, 3.5)),))
    assert result.total_mass_kg[0] == pytest.approx(.2)
    assert result.receptor_traces_kg_m3['initial'][0] == pytest.approx(.2)
    assert result.injected_mass_kg[-1] == 0.0
    assert result.total_mass_kg[-1] == pytest.approx(.2)
    assert result.maximum_mass_residual_kg < 1e-12


def test_frozen_stencils_match_physical_boundary_clamped_trilinear_samples():
    from degali.addons.transient_dense_gas_3d import _receptor_stencils, _sample_trilinear
    x, y, z = np.arange(4)+.5, np.arange(6)+.5, np.arange(3)+.5
    field = np.random.default_rng(931).random((3, 6, 4))
    receptors = tuple(DenseGas3DReceptor(str(i), p) for i,p in enumerate(
        ((0., 0., 0.), (4., 6., 3.), (1.2, 2.3, 1.1), (.4, 2.5, .4))))
    iz, iy, ix, weights = _receptor_stencils(x, y, z, receptors)
    sampled = np.sum(field[iz,iy,ix]*weights, axis=0)
    expected = [_sample_trilinear(field,x,y,z,r.position_m) for r in receptors]
    assert sampled == pytest.approx(expected, rel=1e-14, abs=1e-14)


def test_sparse_field_snapshot_callback_is_output_only_and_copies_state():
    config = _config(store_fields=False)
    captured = []

    def capture(time_s, field):
        captured.append((time_s, field.copy()))
        field[:] = 999.0

    result = solve_transient_dense_gas_3d(
        config,
        field_snapshot_interval_s=0.1,
        field_snapshot_callback=capture,
    )

    assert [item[0] for item in captured] == pytest.approx(
        [0.0, 0.1, 0.2, 0.3, 0.4]
    )
    assert result.fields_kg_m3 is None
    assert result.maximum_mass_residual_kg < 1.0e-10
    assert result.diagnostics["field_snapshot_count"] == 5


def test_full_cross_section_solid_wall_has_no_downstream_transport():
    cells = (30, 8, 8)
    solid = np.zeros((8, 8, 30), dtype=bool)
    solid[:, :, 15] = True
    config = TransientDenseGas3DConfig(
        domain_m=(30.0, 8.0, 8.0),
        cells=cells,
        duration_s=8.0,
        time_step_s=0.1,
        source_schedule=SourceRateSchedule((0.0, 7.999, 8.0), (0.1, 0.1, 0.0)),
        wind_history=WindHistory((0.0, 8.0), (1.0, 1.0), (270.0, 270.0)),
        source_position_m=(5.5, 4.0, 4.0),
        source_sigma_m=0.7,
        source_density_kg_m3=1.2,
        source_h2_mass_fraction=0.5,
        ambient_density_kg_m3=1.2,
        diffusivity_m2_s=0.1,
        prognostic_velocity=False,
        boundary_modes=("closed", "closed", "closed"),
        solid_mask_zyx=solid,
        store_fields=True,
    )
    result = solve_transient_dense_gas_3d(config)
    final = result.fields_kg_m3[-1]
    assert np.max(final[:, :, 15]) == 0.0
    assert np.max(final[:, :, 16:]) < 1.0e-40
    assert result.maximum_mass_residual_kg < 1.0e-12
