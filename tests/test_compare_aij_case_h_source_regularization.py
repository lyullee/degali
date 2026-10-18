from __future__ import annotations

import csv
import json

import pytest

from tools.compare_aij_case_h_source_regularization import compare


def _write_run(directory, *, sigma, time_step=0.001, prediction=0.2):
    directory.mkdir()
    record = {
        "data": {"result_sha256": "same-observation-hash"},
        "configuration": {
            "source_sigma_m": sigma,
            "cell_size_m": 0.025,
            "scalar_time_step_s": time_step,
            "scalar_mean_window_s": [8.0, 12.0],
        },
        "runtime": {"wall_clock_s": 10.0},
        "numerics": {"maximum_relative_species_mass_residual": 1.0e-13},
    }
    (directory / "aij-case-h-tracer-score.json").write_text(
        json.dumps(record), encoding="utf-8",
    )
    with (directory / "aij-case-h-tracer-point-scores.csv").open(
        "w", encoding="utf-8", newline="",
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "point_id", "measured_ethylene_ppmv",
                "predicted_mean_ethylene_ppmv",
            ],
        )
        writer.writeheader()
        writer.writerows([
            {"point_id": 1, "measured_ethylene_ppmv": 0.1,
             "predicted_mean_ethylene_ppmv": prediction},
            {"point_id": 2, "measured_ethylene_ppmv": 0.3,
             "predicted_mean_ethylene_ppmv": prediction * 2.0},
        ])


def test_source_regularization_comparison_accepts_only_sigma_change(tmp_path):
    baseline = tmp_path / "baseline"
    narrow = tmp_path / "narrow"
    _write_run(baseline, sigma=0.025, prediction=0.2)
    _write_run(narrow, sigma=0.001, prediction=0.25)

    result = compare(baseline, narrow, tmp_path / "comparison")

    assert result["only_changed_configuration_field"] == "source_sigma_m"
    assert result["common_measurement_point_count"] == 2
    assert result["source_shape_resolution"]["narrow_underresolved"] is True
    assert result["promotion_allowed"] is False
    assert (tmp_path / "comparison" / "aij-case-h-source-regularization-sensitivity.png").exists()


def test_source_regularization_comparison_rejects_mixed_time_step(tmp_path):
    baseline = tmp_path / "baseline"
    narrow = tmp_path / "narrow"
    _write_run(baseline, sigma=0.025)
    _write_run(narrow, sigma=0.001, time_step=0.002)

    with pytest.raises(ValueError, match="only configuration change"):
        compare(baseline, narrow, tmp_path / "comparison")
