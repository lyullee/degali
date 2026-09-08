from tools.audit_downstream_residual_map import nearest, thermal_classification


def test_nearest_enforces_distance_tolerance():
    rows = [{"x": 1.0}, {"x": 2.0}]
    row, offset = nearest(rows, 1.005, "x", tolerance=0.01)
    assert row["x"] == 1.0
    assert offset == -0.004999999999999893
    assert nearest(rows, 1.05, "x", tolerance=0.01) is None


def test_classification_rejects_small_mean_ke_and_preserves_sign_change():
    observations = {
        "minimum": {
            "temperature_gap_K": 20.0,
            "conditional_specific_enthalpy_gap_J_kg": 200.0,
            "positive_enthalpy_gap_over_local_mean_ke": 20.0,
        },
        "p05": {
            "temperature_gap_K": 10.0,
            "conditional_specific_enthalpy_gap_J_kg": 100.0,
            "positive_enthalpy_gap_over_local_mean_ke": 15.0,
        },
        "median": {
            "temperature_gap_K": -5.0,
            "conditional_specific_enthalpy_gap_J_kg": -50.0,
            "positive_enthalpy_gap_over_local_mean_ke": 0.0,
        },
    }
    result = thermal_classification(observations, 0.2)
    assert result["mean_ke_primary_explanation_rejected"] is True
    assert result["minimum_positive_enthalpy_to_mean_ke_ratio"] == 15.0
    assert result["accounted_mechanical_temperature_scale_below_1K"] is True
    assert result["minimum_median_residual_sign_change"] is True
