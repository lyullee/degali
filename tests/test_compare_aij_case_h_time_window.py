from __future__ import annotations

import csv
import json

import pytest

from tools.compare_aij_case_h_time_window import compare


def _write_run(directory, *, duration, window, released_mass, source_rate=1.0):
    directory.mkdir()
    record = {
        "data": {"result_sha256": "same-observation-hash"},
        "configuration": {
            "source_sigma_m": 0.025,
            "scalar_duration_s": duration,
            "scalar_mean_window_s": window,
            "integrated_released_species_mass_kg": released_mass,
            "source_mass_rate_kg_s": source_rate,
        },
        "runtime": {"wall_clock_s": duration * 2.0},
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
             "predicted_mean_ethylene_ppmv": 0.2},
            {"point_id": 2, "measured_ethylene_ppmv": 0.3,
             "predicted_mean_ethylene_ppmv": 0.4},
        ])


def test_time_window_comparison_accepts_only_declared_horizon_changes(tmp_path):
    short = tmp_path / "short"
    full = tmp_path / "full"
    _write_run(short, duration=12, window=[8, 12], released_mass=1.0)
    _write_run(full, duration=120, window=[0.001, 120], released_mass=10.0)

    result = compare(short, full, tmp_path / "comparison")

    assert result["common_measurement_point_count"] == 2
    assert result["promotion_allowed"] is False
    assert (tmp_path / "comparison" / "aij-case-h-time-window-sensitivity.png").exists()


def test_time_window_comparison_rejects_changed_source_rate(tmp_path):
    short = tmp_path / "short"
    full = tmp_path / "full"
    _write_run(short, duration=12, window=[8, 12], released_mass=1.0)
    _write_run(full, duration=120, window=[0.001, 120], released_mass=10.0,
               source_rate=1.1)

    with pytest.raises(ValueError, match="only duration"):
        compare(short, full, tmp_path / "comparison")
