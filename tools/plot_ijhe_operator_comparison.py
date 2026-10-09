"""Render the peak-versus-common-window operator comparison for IJHE Fig. 3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def plot(baseline_manifest: Path, aligned_manifest: Path, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    baseline = json.loads(baseline_manifest.read_text(encoding="utf-8"))
    aligned = json.loads(aligned_manifest.read_text(encoding="utf-8"))
    peak = baseline["concentration_residual_map"]["frozen_arc_operator"]["statistics"]
    mean = aligned.get("aggregate", aligned.get("primary_aggregate_statistics"))
    if not isinstance(mean, dict):
        raise ValueError("aligned manifest does not contain aggregate statistics")
    labels = ["MG\n(observed/predicted)", "VG", "FAC2"]
    peak_values = [float(peak[key]) for key in ("mg", "vg", "fac2")]
    mean_values = [float(mean[key]) for key in ("mg", "vg", "fac2")]

    import numpy as np

    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.4, 4.9), constrained_layout=True)
    bars_peak = ax.bar(x - width / 2, peak_values, width, label="Peak / arc maximum", color="#0b4f71")
    bars_mean = ax.bar(x + width / 2, mean_values, width, label="Synchronised 20 s mean", color="#d28b20")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Metric value")
    ax.set_title("Observation operator changes the validation conclusion", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    for bars in (bars_peak, bars_mean):
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{bar.get_height():.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )
    ax.text(
        0.01,
        -0.23,
        "Peak score: 62 arcs from 9 trials; common-window score: 46 arcs from 7 trials.\n"
        "These are different observables and are not pooled.",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="svg", metadata={"Title": "IJHE observation operator comparison"})
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_manifest", type=Path)
    parser.add_argument("aligned_manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    plot(args.baseline_manifest, args.aligned_manifest, args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
