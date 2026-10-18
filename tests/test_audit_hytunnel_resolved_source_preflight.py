import runpy

import pytest


TOOL = runpy.run_path("tools/audit_hytunnel_resolved_source_preflight.py")


def test_preflight_config_is_explicitly_coarse_and_not_promotable():
    handoff = TOOL["_handoff"](0.06)
    config = TOOL["_case_config"](handoff)
    assert config.domain_m == (2.52, 0.90, 0.90)
    assert config.cells == (56, 20, 20)
    assert config.active_momentum_solver == "finite_volume_projection"
    assert config.nonstreamwise_wall_boundary == "all_pressure_open"
    assert config.flow_boundary_scheme == "characteristic_multidimensional"
    assert config.source_schedule.rate_operator == "piecewise_constant"
    assert config.source_schedule.released_mass_kg == pytest.approx(
        0.005 * config.duration_s
    )
    assert handoff.equivalent_diameter_m / 0.045 > 4.0


def test_preflight_classification_rejects_nominal_and_only_nominates_margin():
    candidates = [
        {"mean_handoff_density_departure_limit": 0.10, "passed": False},
        {"mean_handoff_density_departure_limit": 0.08, "passed": False},
        {"mean_handoff_density_departure_limit": 0.06, "passed": True},
    ]
    result = TOOL["classify_preflight"](candidates)
    assert result == {
        "status": "nominal_resolved_source_profile_preflight_failed",
        "nominal_0p10_passed": False,
        "least_conservative_passing_candidate": 0.06,
        "promotion_allowed": False,
    }


def test_preflight_classification_requires_candidates():
    with pytest.raises(ValueError, match="at least one"):
        TOOL["classify_preflight"]([])
