"""Render the saved FFI strategy diagnostics; does not run or tune models."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--original-paired-csv", type=Path, required=True)
    args = parser.parse_args()
    with args.original_paired_csv.open(encoding="utf-8-sig", newline="") as stream:
        positions = {row["sensor"]: row for row in csv.DictReader(stream)}
    with (args.analysis_dir / "test4_replay_paired_predictions.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    analysis = json.loads((args.analysis_dir / "analysis.json").read_text(encoding="utf-8"))
    models = [
        ("observed_vol_pct", "Observed", "#17202a", "o", ""),
        ("DEGALI_static_high_10m", "DEGALI static", "#c05032", "s", "--"),
        ("DEGALI_replay_high_speed_direction_refined", "DEGALI measured-wind replay", "#2766a0", "o", "-"),
        ("SLABx_native", "SLABx native", "#8d9195", "^", ":"),
        ("SLABx_LH2_retained_research", "SLABx-LH2 research", "#2a7b55", "D", "-"),
    ]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    for axis, radius in zip(axes, (30.0, 50.0, 100.0)):
        grouped = {}
        for row in rows:
            if float(row["radius_m"]) != radius:
                continue
            position = positions[row["sensor"]]
            bearing = round(90.0 + math.degrees(math.atan2(float(position["y_crosswind_m"]), float(position["x_downwind_m"]))), 6)
            grouped.setdefault(bearing, []).append(row)
        bearings = sorted(grouped)
        for key, label, colour, marker, style in models:
            values = [np.mean([float(row[key]) for row in grouped[bearing]]) for bearing in bearings]
            axis.plot(bearings, values, label=label, color=colour, marker=marker,
                      linestyle=style, linewidth=1.6, markersize=5)
        axis.set_title(f"{radius:g} m arc")
        axis.set_xlabel("Sensor bearing (degrees)")
        axis.set_xticks(bearings)
        axis.set_ylim(bottom=0)
        axis.grid(axis="y", color="#e1e4e8", linewidth=0.6)
    axes[0].set_ylabel("Mean across three sensor heights (vol% H2)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.07))
    fig.suptitle("FFI Test 4: measured wind repairs the static core, but misses the lateral field", fontsize=13, y=0.97)
    text = "25-300 s wind bins; three measured heights per bearing. Lines join sampled positions only.\nReplay rotates steady fields instantaneously; source histories, travel-time memory and sensor lag remain unresolved."
    fig.text(0.5, 0.01, text, ha="center", va="bottom", fontsize=8, color="#555b61")
    fig.tight_layout(rect=(0.0, 0.23, 1.0, 0.92))
    fig.savefig(args.analysis_dir / "test4_strategy_diagnostic.png", dpi=180)
    fig.savefig(args.analysis_dir / "test4_strategy_diagnostic.pdf")
    plt.close(fig)
    replay_key = "DEGALI_replay_high_speed_direction_refined"
    summary = {
        "static_mae": analysis["scores"]["DEGALI_static_high_10m"]["mae_vol_pct_point"],
        "replay_mae": analysis["scores"][replay_key]["mae_vol_pct_point"],
        "research_mae": analysis["scores"]["SLABx_LH2_retained_research"]["mae_vol_pct_point"],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
