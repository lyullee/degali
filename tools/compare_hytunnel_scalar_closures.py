"""Compare legacy, local-SGS and measured-mean-fan HyTunnel candidates.

The comparison fixes the 5% source handoff, Sc_t=0.7, grid, duration and
observation points.  It diagnoses the closure change; it is not a calibration
or permission to select a coefficient from the measured concentrations.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEGACY = ROOT / "outputs/hytunnel-source-mixing-sensitivity-2026-10-11-v1.json"
DEFAULT_LOCAL_SGS = ROOT / "outputs/hytunnel-resolved-sgs-diffusivity-handoff0p05-score-2026-10-11-v1.json"
DEFAULT_MEASURED_FAN = ROOT / "outputs/hytunnel-measured-fan-sgs-handoff0p05-score-2026-10-11-v1.json"
DEFAULT_OUTPUT = ROOT / "outputs/hytunnel-scalar-closure-comparison-2026-10-11-v1.json"
DEFAULT_PLOT = ROOT / "outputs/hytunnel-scalar-closure-comparison-2026-10-11-v1.png"


def compare_closures(
    legacy_path: str | Path = DEFAULT_LEGACY,
    local_sgs_path: str | Path = DEFAULT_LOCAL_SGS,
    measured_fan_path: str | Path = DEFAULT_MEASURED_FAN,
) -> dict[str, Any]:
    legacy_source = Path(legacy_path)
    candidate_source = Path(local_sgs_path)
    measured_fan_source = Path(measured_fan_path)
    legacy_record = json.loads(legacy_source.read_text(encoding="utf-8"))
    candidate_record = json.loads(candidate_source.read_text(encoding="utf-8"))
    measured_fan_record = json.loads(
        measured_fan_source.read_text(encoding="utf-8")
    )
    legacy = legacy_record["cases"]["handoff0p05_sc0p7"][
        "metrics_on_common_points"
    ]
    candidate_case = candidate_record["case_scores"]["published_box"]
    candidate = candidate_case["all_supported_points"]
    candidate_four = candidate_case["threshold_classification"]["four_vol_pct"]
    measured_fan_case = measured_fan_record["case_scores"]["published_box"]
    measured_fan = measured_fan_case["all_supported_points"]
    measured_fan_four = measured_fan_case["threshold_classification"]["four_vol_pct"]
    legacy_four = legacy["four_vol_confusion"]
    if int(legacy["pair_count"]) != int(candidate["pair_count"]):
        raise ValueError("closure comparison requires the same supported points")
    if int(legacy["pair_count"]) != int(measured_fan["pair_count"]):
        raise ValueError("fan-profile comparison requires the same supported points")

    cases = {
        "legacy_uniform_nu_base_over_sc": {
            "artifact": str(legacy_source),
            "diffusivity_model": "legacy_total_kinematic_viscosity_ratio",
            "mae_vol_pct_point": float(legacy["mae_vol_pct_point"]),
            "rmse_vol_pct_point": float(legacy["rmse_vol_pct_point"]),
            "mean_bias_vol_pct_point": float(legacy["mean_bias_vol_pct_point"]),
            "pearson_r": float(legacy["pearson_r"]),
            "fac2": float(legacy["fac2_floor_0p10"]),
            "four_vol_sensitivity": float(legacy_four["sensitivity"]),
            "four_vol_specificity": float(legacy_four["specificity"]),
            "false_positive": int(legacy_four["false_positive"]),
            "false_negative": int(legacy_four["false_negative"]),
        },
        "molecular_plus_local_sgs_over_sc": {
            "artifact": str(candidate_source),
            "diffusivity_model": "molecular_plus_subgrid_eddy_viscosity",
            "mae_vol_pct_point": float(candidate["mae_vol_pct_point"]),
            "rmse_vol_pct_point": float(candidate["rmse_vol_pct_point"]),
            "mean_bias_vol_pct_point": float(
                candidate["mean_bias_predicted_minus_observed_vol_pct_point"]
            ),
            "pearson_r": float(candidate["pearson_r"]),
            "fac2": float(
                candidate_case["hanna_statistics_floor_0p10_vol_pct"]["fac2"]
            ),
            "four_vol_sensitivity": float(candidate_four["sensitivity"]),
            "four_vol_specificity": float(candidate_four["specificity"]),
            "false_positive": int(candidate_four["false_positive"]),
            "false_negative": int(candidate_four["false_negative"]),
        },
        "measured_mean_fan_plus_local_sgs_over_sc": {
            "artifact": str(measured_fan_source),
            "diffusivity_model": "molecular_plus_subgrid_eddy_viscosity",
            "fan_inlet_profile": "measured_first_station",
            "mae_vol_pct_point": float(measured_fan["mae_vol_pct_point"]),
            "rmse_vol_pct_point": float(measured_fan["rmse_vol_pct_point"]),
            "mean_bias_vol_pct_point": float(
                measured_fan["mean_bias_predicted_minus_observed_vol_pct_point"]
            ),
            "pearson_r": float(measured_fan["pearson_r"]),
            "fac2": float(
                measured_fan_case["hanna_statistics_floor_0p10_vol_pct"]["fac2"]
            ),
            "four_vol_sensitivity": float(measured_fan_four["sensitivity"]),
            "four_vol_specificity": float(measured_fan_four["specificity"]),
            "false_positive": int(measured_fan_four["false_positive"]),
            "false_negative": int(measured_fan_four["false_negative"]),
        },
    }
    old = cases["legacy_uniform_nu_base_over_sc"]
    new = cases["molecular_plus_local_sgs_over_sc"]
    measured = cases["measured_mean_fan_plus_local_sgs_over_sc"]
    comparison_keys = (
        "mae_vol_pct_point", "rmse_vol_pct_point",
        "mean_bias_vol_pct_point", "pearson_r", "fac2",
        "four_vol_specificity",
    )
    return {
        "schema": "degali.hytunnel-scalar-closure-comparison.v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "local_sgs_and_measured_mean_fan_candidates_rejected_by_frozen_"
            "accuracy_gates"
        ),
        "promotion_allowed": False,
        "calibration_or_best_variant_selection_allowed": False,
        "fixed_basis": {
            "handoff_density_departure_limit": 0.05,
            "effective_schmidt_number": 0.7,
            "common_supported_point_count": int(legacy["pair_count"]),
            "duration_s": float(candidate_record["model_duration_s"]),
        },
        "cases": cases,
        "candidate_minus_legacy": {
            key: float(new[key] - old[key])
            for key in comparison_keys
        },
        "measured_fan_minus_uniform_local_sgs": {
            key: float(measured[key] - new[key])
            for key in comparison_keys
        },
        "measured_fan_minus_legacy": {
            key: float(measured[key] - old[key])
            for key in comparison_keys
        },
        "interpretation": (
            "Separating molecular and SGS diffusion removes the legacy numerical "
            "mixing compensation, but local WALE alone under-mixes this coarse "
            "carrier field. Transferring the independently measured first mean fan "
            "plane does not materially improve the frozen concentration metrics. "
            "The next gate is measured or independently resolved fan turbulence and "
            "scalar-flux evidence, not fitting Sc_t to this concentration sheet."
        ),
    }


def plot_comparison(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    old = record["cases"]["legacy_uniform_nu_base_over_sc"]
    new = record["cases"]["molecular_plus_local_sgs_over_sc"]
    measured = record["cases"]["measured_mean_fan_plus_local_sgs_over_sc"]
    labels = ("Legacy uniform", "Uniform + local SGS", "Measured mean fan + local SGS")
    colours = ("#6b7280", "#d97706", "#2563eb")
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    error_names = ("MAE", "RMSE", "Bias")
    keys = ("mae_vol_pct_point", "rmse_vol_pct_point", "mean_bias_vol_pct_point")
    x = range(len(keys))
    width = 0.25
    axes[0].bar([value - width for value in x], [old[key] for key in keys], width, label=labels[0], color=colours[0])
    axes[0].bar(list(x), [new[key] for key in keys], width, label=labels[1], color=colours[1])
    axes[0].bar([value + width for value in x], [measured[key] for key in keys], width, label=labels[2], color=colours[2])
    axes[0].set_xticks(list(x), error_names)
    axes[0].set_ylabel("vol%-point")
    axes[0].set_title("Concentration error (45 common points)")
    axes[0].legend(frameon=False)

    score_names = ("r", "FAC2", "4% specificity")
    score_keys = ("pearson_r", "fac2", "four_vol_specificity")
    axes[1].bar([value - width for value in x], [old[key] for key in score_keys], width, color=colours[0])
    axes[1].bar(list(x), [new[key] for key in score_keys], width, color=colours[1])
    axes[1].bar([value + width for value in x], [measured[key] for key in score_keys], width, color=colours[2])
    axes[1].axhline(0.90, color="black", linestyle="--", linewidth=1.0, label="4% target")
    axes[1].set_xticks(list(x), score_names)
    axes[1].set_ylim(0.0, 1.0)
    axes[1].set_title("Shape and threshold performance")
    axes[1].legend(frameon=False)
    figure.suptitle("HyTunnel scalar closure audit — no coefficient fitting")
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--local-sgs", type=Path, default=DEFAULT_LOCAL_SGS)
    parser.add_argument("--measured-fan", type=Path, default=DEFAULT_MEASURED_FAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args()
    record = compare_closures(args.legacy, args.local_sgs, args.measured_fan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    plot_comparison(record, args.plot)
    print(json.dumps({
        "status": record["status"],
        "candidate_minus_legacy": record["candidate_minus_legacy"],
        "measured_fan_minus_uniform_local_sgs": (
            record["measured_fan_minus_uniform_local_sgs"]
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
