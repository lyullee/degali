import numpy as np
import pytest
from types import SimpleNamespace

from degali.addons.pool_evaporation import (
    SolidSubstrate,
    semi_infinite_heat_flux,
    substrate_conduction_evaporation,
)
import degali.lh2 as lh2


def _solid():
    return SolidSubstrate(
        conductivity_w_m_k=1.5, density_kg_m3=2200.0,
        heat_capacity_j_kg_k=900.0, initial_temperature_k=293.15,
        depth_m=1.0, cells=80,
    )


def test_conduction_pool_cools_and_evaporation_rate_falls():
    result = substrate_conduction_evaporation(
        _solid(), area_m2=0.25, duration_s=600.0, time_step_s=5.0,
    )
    assert result.total_evaporated_mass_kg > 0.0
    assert result.steps[-1].surface_adjacent_temperature_k < 293.15
    assert result.steps[-1].vapour_rate_kg_s < result.steps[0].vapour_rate_kg_s
    assert result.depth_is_effectively_semi_infinite


def test_finite_volume_flux_tracks_analytic_semi_infinite_early_time():
    substrate = SolidSubstrate(
        conductivity_w_m_k=1.5, density_kg_m3=2200.0,
        heat_capacity_j_kg_k=900.0, initial_temperature_k=293.15,
        depth_m=1.0, cells=800,
    )
    result = substrate_conduction_evaporation(
        substrate, area_m2=1.0, duration_s=100.0, time_step_s=1.0,
    )
    numerical = result.steps[-1].heat_flux_w_m2
    analytic = semi_infinite_heat_flux(substrate, elapsed_s=100.0)
    assert numerical == pytest.approx(analytic, rel=0.12)


def test_pool_inventory_caps_the_conduction_source():
    result = substrate_conduction_evaporation(
        _solid(), area_m2=1.0, duration_s=100.0, time_step_s=2.0,
        initial_liquid_mass_kg=1.0e-4,
    )
    assert result.total_evaporated_mass_kg == pytest.approx(1.0e-4)
    assert result.steps[-1].remaining_liquid_kg == pytest.approx(0.0)
    assert result.steps[-1].vapour_rate_kg_s == pytest.approx(0.0)


def test_declared_inflow_forms_pool_and_closes_inventory():
    result = substrate_conduction_evaporation(
        _solid(), area_m2=0.5, duration_s=20.0, time_step_s=1.0,
        initial_liquid_mass_kg=0.0, liquid_inflow_rate_kg_s=0.01,
        inflow_duration_s=5.0,
    )
    assert result.total_liquid_inflow_kg == pytest.approx(0.05)
    assert result.total_evaporated_mass_kg + result.steps[-1].remaining_liquid_kg \
        == pytest.approx(0.05)
    assert result.steps[0].liquid_inflow_rate_kg_s == pytest.approx(0.01)
    assert result.steps[-1].liquid_inflow_rate_kg_s == pytest.approx(0.0)


def test_invalid_pool_inputs_are_rejected():
    with pytest.raises(ValueError):
        SolidSubstrate(1, 1, 1, 293, 1, cells=2)
    with pytest.raises(ValueError):
        substrate_conduction_evaporation(_solid(), area_m2=0, duration_s=1, time_step_s=1)


def test_time_resolved_pool_path_keeps_transport_as_an_explicit_label(monkeypatch):
    evaporation = substrate_conduction_evaporation(
        _solid(), area_m2=.25, duration_s=3.0, time_step_s=1.0,
    )
    calls = []

    def fake_assess(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(trajectory=np.array([[0.0, 0.0, 1.0], [10.0, 0.0, 0.1]]))

    monkeypatch.setattr(lh2, "assess", fake_assess)
    history = lh2.assess_pool_history(
        evaporation, wind=2.0, pool_diameter=.5,
        receptor_distance=4.0, propagation_speed=1.0,
    )
    assert len(history.snapshots) == 3
    assert history.snapshots[0].source_interval_start_s == pytest.approx(0.0)
    assert history.snapshots[0].receptor_time_s == pytest.approx(5.0)
    assert history.snapshots[0].receptor_interval_start_s == pytest.approx(4.0)
    assert history.snapshots[0].receptor_centreline_mole_fraction == pytest.approx(.64)
    assert calls[0]["at_distance"] == pytest.approx(4.0)
    with pytest.raises(ValueError, match="propagation_speed"):
        lh2.assess_pool_history(evaporation, wind=2.0, pool_diameter=.5, receptor_distance=2.0)


def test_time_resolved_pool_optional_response_kernel_is_causal(monkeypatch):
    evaporation = substrate_conduction_evaporation(
        _solid(), area_m2=.25, duration_s=3.0, time_step_s=1.0,
    )

    def fake_assess(**kwargs):
        return SimpleNamespace(trajectory=np.array([[0.0, 0.0, 1.0], [10.0, 0.0, 0.1]]))

    monkeypatch.setattr(lh2, "assess", fake_assess)
    history = lh2.assess_pool_history(
        evaporation, wind=2.0, pool_diameter=.5,
        receptor_distance=4.0, propagation_speed=1.0, response_time_s=2.0,
    )
    values = [item.receptor_transient_mole_fraction for item in history.snapshots]
    assert history.transient_kernel_used
    assert history.response_time_s == pytest.approx(2.0)
    assert values[0] < values[1] < values[2] < 0.64
    with pytest.raises(ValueError, match="requires receptor_distance"):
        lh2.assess_pool_history(
            evaporation, wind=2.0, pool_diameter=.5, response_time_s=2.0,
        )
