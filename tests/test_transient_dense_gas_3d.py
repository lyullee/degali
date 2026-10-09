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
