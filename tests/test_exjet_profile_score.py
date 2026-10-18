from __future__ import annotations

import pytest
import numpy as np
from types import SimpleNamespace

from degali.validation.exjet import (
    sample_exjet_spatial_profile_predictions,
    score_exjet_profile_points,
)


def _observations():
    return [
        {
            "observation_id": "center",
            "station_x_m": 1.4,
            "radial_coordinate_m": 0.0,
            "coordinate_uncertainty_m": 0.03,
            "h2_molar_fraction": 0.05,
            "quantity_kind": "centreline_hydrogen_molar_fraction",
            "concentration_uncertainty_mole_fraction": 0.01,
        },
        {
            "observation_id": "radial",
            "station_x_m": 2.5,
            "radial_coordinate_m": 0.7,
            "coordinate_uncertainty_m": 0.03,
            "h2_molar_fraction": 0.03,
            "quantity_kind": "radial_hydrogen_molar_fraction",
            "concentration_uncertainty_mole_fraction": 0.015,
        },
    ]


def test_exjet_point_score_reports_vol_percent_errors_and_sparse_4vol_gate():
    result = score_exjet_profile_points(
        _observations(),
        [
            {
                "observation_id": "center", "h2_molar_fraction_predicted": 0.07,
                "h2_molar_fraction_coordinate_lower": 0.035,
                "h2_molar_fraction_coordinate_upper": 0.08,
            },
            {
                "observation_id": "radial", "h2_molar_fraction_predicted": 0.02,
                "h2_molar_fraction_coordinate_lower": 0.015,
                "h2_molar_fraction_coordinate_upper": 0.025,
            },
        ],
    )
    assert result["n"] == 2
    center = result["by_quantity"]["centreline_hydrogen_molar_fraction"]
    assert center["mean_bias_vol_percent_point"] == pytest.approx(2.0)
    assert center["mae_vol_percent_point"] == pytest.approx(2.0)
    radial = result["by_quantity"]["radial_hydrogen_molar_fraction"]
    assert radial["mean_bias_vol_percent_point"] == pytest.approx(-1.0)
    assert radial["pointwise_4_vol_percent_classification"]["indeterminate_observation_count"] == 1
    assert radial["pointwise_4_vol_percent_classification"]["confusion_matrix_n"] == 0
    assert center["pointwise_4_vol_percent_classification"]["indeterminate_prediction_count"] == 1
    assert result["dynamic_validation_supported"] is False
    assert result["coordinate_uncertainty_propagated"] is True
    assert result["quantitative_promotion_allowed"] is False


def test_exjet_score_rejects_missing_or_unexpected_prediction_ids():
    with pytest.raises(ValueError, match="IDs must match"):
        score_exjet_profile_points(
            _observations(),
            [{"observation_id": "center", "h2_molar_fraction_predicted": 0.05}],
        )


def test_exjet_score_rejects_nonphysical_mole_fraction():
    with pytest.raises(ValueError, match=r"in \[0, 1\]"):
        score_exjet_profile_points(
            _observations(),
            [
                {"observation_id": "center", "h2_molar_fraction_predicted": 1.1},
                {"observation_id": "radial", "h2_molar_fraction_predicted": 0.02},
            ],
        )


def test_resolved_field_sampler_uses_explicit_exjet_vertical_plane_transform():
    x = np.asarray((0.0, 1.0, 2.0))
    y = np.asarray((-1.0, 0.0, 1.0))
    z = np.asarray((0.0, 1.0, 2.0))
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    field = (xx + 2.0 * yy + zz) / 10.0
    result = SimpleNamespace(
        x_m=x,
        y_m=y,
        z_m=z,
        hydrogen_mole_fraction=field,
        solid_mask=np.zeros_like(field, dtype=bool),
        time_s=np.asarray((0.0, 4.0)),
    )
    sampled = sample_exjet_spatial_profile_predictions(
        result,
        [{
            "observation_id": "p", "station_x_m": 1.0,
            "radial_coordinate_m": 0.5, "coordinate_uncertainty_m": 0.03,
        }],
        nozzle_position_m=(0.0, 0.0, 0.5),
    )
    assert sampled["predictions"][0]["h2_molar_fraction_predicted"] == pytest.approx(0.2)
    assert sampled["predictions"][0]["h2_molar_fraction_coordinate_lower"] == pytest.approx(0.197)
    assert sampled["predictions"][0]["h2_molar_fraction_coordinate_upper"] == pytest.approx(0.203)
    assert sampled["predictions"][0]["sample_coordinate_solver_m"] == [1.0, 0.0, 1.0]
    assert sampled["sampled_time_s"] == pytest.approx(4.0)
    assert sampled["sampled_field_is_final_snapshot"] is True
    assert sampled["coordinate_uncertainty_propagated"] is True
    assert sampled["complete_pairing"] is True
    assert sampled["quantitative_promotion_allowed"] is False


def test_resolved_field_sampler_withholds_solid_stencil_instead_of_fallback():
    x = y = z = np.asarray((0.0, 1.0, 2.0))
    field = np.full((3, 3, 3), 0.2)
    solid = np.zeros_like(field, dtype=bool)
    solid[1, 1, 1] = True
    result = SimpleNamespace(
        x_m=x, y_m=y, z_m=z, hydrogen_mole_fraction=field,
        solid_mask=solid, time_s=np.asarray((0.0, 1.0)),
    )
    sampled = sample_exjet_spatial_profile_predictions(
        result,
        [{
            "observation_id": "p", "station_x_m": 1.0,
            "radial_coordinate_m": 0.0, "coordinate_uncertainty_m": 0.03,
        }],
        nozzle_position_m=(0.0, 1.0, 1.0),
    )
    assert sampled["predictions"] == []
    assert "intersects_solid" in sampled["failed_points"][0]["reason"]
    assert sampled["complete_pairing"] is False


def test_coordinate_uncertainty_is_not_claimed_when_an_endpoint_leaves_domain():
    x = y = z = np.asarray((0.0, 1.0, 2.0))
    field = np.full((3, 3, 3), 0.2)
    result = SimpleNamespace(
        x_m=x, y_m=y, z_m=z, hydrogen_mole_fraction=field,
        solid_mask=np.zeros_like(field, dtype=bool),
        time_s=np.asarray((0.0, 1.0)),
    )
    sampled = sample_exjet_spatial_profile_predictions(
        result,
        [{
            "observation_id": "near-top",
            "station_x_m": 1.0,
            "radial_coordinate_m": 1.49,
            "coordinate_uncertainty_m": 0.03,
        }],
        nozzle_position_m=(0.0, 0.0, 0.5),
    )
    assert sampled["complete_pairing"] is True
    assert sampled["complete_coordinate_uncertainty_propagation"] is False
    assert sampled["predictions"][0]["h2_molar_fraction_coordinate_lower"] is None


def test_sample_then_score_preserves_coordinate_bounds_end_to_end():
    x = np.asarray((0.0, 2.0, 4.0, 6.0))
    y = np.asarray((-1.0, 0.0, 1.0))
    z = np.asarray((0.0, 1.0, 2.0))
    field = np.full((3, 3, 4), 0.1)
    result = SimpleNamespace(
        x_m=x, y_m=y, z_m=z, hydrogen_mole_fraction=field,
        solid_mask=np.zeros_like(field, dtype=bool),
        time_s=np.asarray((0.0, 3.0)),
    )
    observations = _observations()
    sampled = sample_exjet_spatial_profile_predictions(
        result, observations, nozzle_position_m=(0.0, 0.0, 0.5),
    )
    assert sampled["complete_coordinate_uncertainty_propagation"] is True
    score = score_exjet_profile_points(observations, sampled["predictions"])
    assert score["coordinate_uncertainty_propagated"] is True
    assert score["by_quantity"]["radial_hydrogen_molar_fraction"][
        "combined_interval_comparison_n"
    ] == 1
