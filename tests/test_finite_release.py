from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.finite_release import (
    assess_steady_wind_applicability,
    finite_release_puff_handoff,
)
from degali.addons.transient_receptor import WindHistory


class _Model:
    @staticmethod
    def section_widths(_state):
        return 2.0, 1.0


def _result(hydrogen_flux=0.1):
    arc = np.array([0.0, 10.0, 20.0])
    states = np.array([
        [1.0, 0.1, 2.0, 0.0, 2.0, x, 1.0] for x in arc
    ])
    fluxes = np.array([
        [2.0, hydrogen_flux, 4.0, 0.0, 100.0] for _ in arc
    ])
    return SimpleNamespace(arc_length=arc, states=states, fluxes=fluxes)


def test_finite_release_builds_full_section_conservative_inventory():
    handoff = finite_release_puff_handoff(
        _Model(), _result(), source_duration_s=5.0,
        source_hydrogen_mass_flow_kg_s=0.1,
    )
    assert handoff.status == "transition_ready"
    assert handoff.transition_arc_length_m == pytest.approx(10.0)
    assert handoff.centre_position_m == pytest.approx((10.0, 0.0, 1.0))
    assert handoff.total_mass_kg == pytest.approx(10.0)
    assert handoff.hydrogen_mass_kg == pytest.approx(0.5)
    assert handoff.momentum_kg_m_s == pytest.approx((20.0, 0.0, 0.0))
    assert handoff.longitudinal_half_length_m == pytest.approx(5.0)
    assert handoff.hydrogen_inventory_residual_kg == pytest.approx(0.0)
    assert handoff.requires_transient_puff_continuation


def test_finite_release_refuses_carrier_mass_as_hydrogen_clock():
    with pytest.raises(ValueError, match="carrier-mixture"):
        finite_release_puff_handoff(
            _Model(), _result(hydrogen_flux=0.2), source_duration_s=5.0,
            source_hydrogen_mass_flow_kg_s=0.1,
        )


def test_finite_release_reports_when_transition_is_beyond_domain():
    handoff = finite_release_puff_handoff(
        _Model(), _result(), source_duration_s=15.0,
        source_hydrogen_mass_flow_kg_s=0.1,
    )
    assert handoff.status == "transition_beyond_trajectory"
    assert not handoff.requires_transient_puff_continuation


def test_steady_wind_gate_handles_north_wrap_and_rejects_large_turn():
    wrap = WindHistory(
        time_s=[0.0, 5.0, 10.0],
        speed_m_s=[2.0, 2.1, 2.0],
        direction_from_deg=[359.0, 0.0, 1.0],
    )
    accepted = assess_steady_wind_applicability(wrap)
    assert accepted.applicable
    assert accepted.direction_span_deg == pytest.approx(2.0)

    turn = WindHistory(
        time_s=[0.0, 5.0, 10.0],
        speed_m_s=[2.0, 2.0, 2.0],
        direction_from_deg=[270.0, 290.0, 315.0],
    )
    rejected = assess_steady_wind_applicability(turn)
    assert not rejected.applicable
    assert "wind_direction_span_exceeds_steady_limit" in rejected.reasons
