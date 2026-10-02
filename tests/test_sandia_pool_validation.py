import pytest

from degali.validation.sandia_pool import (
    PoolContourLowerBound,
    read_pool_contour_lower_bounds,
    screen_contour_lower_bound,
    screen_pool_contours,
    screen_report,
)


@pytest.fixture
def contour():
    return PoolContourLowerBound(
        identifier="example", figure="3-27", panel="b", contour_mole_fraction=.04,
        minimum_reach_m=3.0, spatial_resolution_m=.15,
        all_inflow_vaporization_supported=True,
        source_note="local-only reduction",
    )


def test_contour_is_a_one_sided_falsification_not_an_accuracy_target(contour):
    rejected = screen_contour_lower_bound(contour, 2.7)
    survives = screen_contour_lower_bound(contour, 4.0)
    assert rejected.outcome == "falsified_below_resolved_reach"
    assert survives.outcome == "not_falsified"
    assert not rejected.field_accuracy_claim_allowed
    assert not survives.default_promotion_allowed


def test_resolution_is_declared_not_a_hidden_fitted_tolerance(contour):
    assert screen_contour_lower_bound(contour, 2.85).outcome == "not_falsified"
    assert screen_contour_lower_bound(contour, 2.849).outcome == "falsified_below_resolved_reach"
    with pytest.raises(ValueError, match="finite"):
        screen_contour_lower_bound(contour, float("inf"))


def test_local_reader_and_screen_require_complete_unshared_identifiers(tmp_path, contour):
    source = tmp_path / "local.json"
    source.write_text(
        "[{'identifier': 'example'}]".replace("'", '"'), encoding="utf-8"
    )
    with pytest.raises(TypeError):
        read_pool_contour_lower_bounds(source)
    with pytest.raises(ValueError, match="match"):
        screen_pool_contours([contour], {"other": 3.0})
    report = screen_report(screen_pool_contours([contour], {"example": 2.0}))
    assert not report["third_party_observations_distributed"]
    assert not report["accuracy_score_calculated"]
    assert not report["all_lower_bounds_not_falsified"]
