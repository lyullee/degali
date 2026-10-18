import numpy as np
import pytest

from degali.addons.transient_enclosure import (
    EnclosureForcing,
    TransientEnclosureConfig,
    solve_transient_enclosure,
)


def _forcing(*, rate=0.01, ventilation=0.1):
    return EnclosureForcing(
        time_s=(0.0, 20.0),
        h2_mass_rate_kg_s=(rate, rate),
        ventilation_flow_m3_s=(ventilation, ventilation),
        source_id="test-enclosure",
    )


def _config(**updates):
    values = {
        "enclosure_volume_m3": 60.8,
        "effective_mixed_volume_fraction": 0.4,
        "effective_opening_area_m2": 0.6,
        "transport_delay_s": 0.0,
        "opening_height_m": 2.285,
    }
    values.update(updates)
    return TransientEnclosureConfig(**values)


def test_enclosure_volume_fraction_is_bounded_under_large_source():
    result = solve_transient_enclosure(
        _forcing(rate=10.0, ventilation=0.0),
        np.linspace(0.0, 20.0, 201),
        _config(effective_opening_area_m2=0.0),
    )

    assert np.all(result.hydrogen_volume_fraction >= 0.0)
    assert np.all(result.hydrogen_volume_fraction <= 1.0)
    assert result.hydrogen_volume_fraction[-1] > 0.99


def test_enclosure_update_closes_hydrogen_mass_balance():
    result = solve_transient_enclosure(
        _forcing(), np.linspace(0.0, 20.0, 201), _config(),
    )

    assert result.maximum_mass_residual_kg < 1.0e-12
    assert result.inventory_kg[-1] + result.cumulative_removed_kg[-1] == pytest.approx(
        result.cumulative_injected_kg[-1], abs=1.0e-12
    )
    assert result.diagnostics["source_displacement_included"] is True


def test_transport_delay_postpones_enclosure_response():
    time = np.arange(0.0, 20.1, 1.0)
    result = solve_transient_enclosure(
        _forcing(), time, _config(transport_delay_s=5.0),
    )

    assert np.all(result.hydrogen_volume_fraction[time <= 5.0] == 0.0)
    assert result.hydrogen_volume_fraction[time == 6.0][0] > 0.0


def test_larger_opening_exchange_reduces_concentration():
    time = np.linspace(0.0, 20.0, 201)
    closed = solve_transient_enclosure(
        _forcing(), time, _config(effective_opening_area_m2=0.0),
    )
    open_result = solve_transient_enclosure(
        _forcing(), time, _config(effective_opening_area_m2=0.8),
    )

    assert open_result.concentration_vol_percent[-1] < closed.concentration_vol_percent[-1]
    assert open_result.buoyant_opening_exchange_m3_s[-1] > 0.0


def test_forcing_must_cover_evaluation_clock():
    with pytest.raises(ValueError, match="cover the complete evaluation time"):
        solve_transient_enclosure(
            _forcing(), (0.0, 21.0), _config(),
        )
