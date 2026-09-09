"""Frozen target-selection guards for the six-trial flux extension."""

import pytest

from tools.audit_thermal_moment_flux_space_extension import first_downstream_profile_target


def test_first_well_constrained_profile_strictly_downstream_is_selected():
    trial = {"vertical_fits": [
        {"x": .79, "well_constrained": True},
        {"x": 1.78, "well_constrained": False},
        {"x": 4., "well_constrained": True},
    ]}
    assert first_downstream_profile_target(trial, .6) == pytest.approx(.79)
    assert first_downstream_profile_target(trial, .8) == pytest.approx(4.)


def test_missing_downstream_profile_fails_closed():
    with pytest.raises(ValueError, match="no well-constrained"):
        first_downstream_profile_target(
            {"vertical_fits": [{"x": .79, "well_constrained": False}]}, .5,
        )
