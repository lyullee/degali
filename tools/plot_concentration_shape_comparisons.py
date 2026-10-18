"""Build publication-ready concentration-shape and vol%-error figures.

The figures deliberately separate experimental concentration fields from
model-versus-observation comparisons.  FFI Test 4 maps use linear
interpolation only inside the sensor convex hull and are labelled as such;
they are not presented as CFD fields.  HyTunnel and USN figures preserve the
public workbook/MAT values, including values below zero in the audit source.
Display colours begin at zero because negative instrument offsets have no
physical flammability meaning; no source data or scored metrics are changed.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "concentration-shape-comparisons-2026-10-10"
FFI_REPLAY = (
    ROOT / "outputs" / "slabx-competitive-strategy-2026-10-08"
    / "test4_replay_paired_predictions.csv"
)
FFI_COORDINATES = (
    ROOT / "outputs" / "applied-energy-evidence-2026-10-03" / "test4-external"
    / "test4_slabx_degali_paired_sensors.csv"
)
FFI_ARCS = (
    ROOT / "outputs" / "applied-energy-evidence-2026-10-09c" / "paired_arcs.csv"
)
USN_AUDIT = ROOT / "outputs" / "usn-carpark-dynamic-audit-2026-10-10.json"
USN_EXP20 = (
    ROOT / "reference" / "datasets" / "usn-carpark-14405903" / "raw"
    / "Exp20.mat"
)
HYTUNNEL_WORKBOOK = (
    ROOT / "reference" / "datasets" / "hytunnel-ps-4668554" / "raw"
    / "HTE2440PS003MIXED_300321.xlsx"
)


def _plotting():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.lines import Line2D

    return plt, colors, Line2D


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _save(fig, output_dir: Path, stem: str) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.svg"]
    fig.savefig(paths[0], dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(paths[1], bbox_inches="tight", facecolor="white")
    return paths


def _style(plt) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "figure.titlesize": 17,
        "figure.titleweight": "bold",
        "svg.fonttype": "none",
    })


def _ffi_rows() -> list[dict[str, float | str]]:
    coordinate_by_sensor = {
        row["sensor"]: (float(row["x_downwind_m"]), float(row["y_crosswind_m"]))
        for row in _read_csv(FFI_COORDINATES)
    }
    rows: list[dict[str, float | str]] = []
    for row in _read_csv(FFI_REPLAY):
        x_m, y_m = coordinate_by_sensor[row["sensor"]]
        rows.append({
            "sensor": row["sensor"],
            "x_m": x_m,
            "y_m": y_m,
            "height_m": float(row["height_m"]),
            "observed": float(row["observed_vol_pct"]),
            "degali_static": float(row["DEGALI_static_high_10m"]),
            "degali_replay": float(row["DEGALI_replay_high_speed_direction_refined"]),
            "slabx_native": float(row["SLABx_native"]),
            "slabx_research": float(row["SLABx_LH2_retained_research"]),
        })
    return rows


def _mast_maxima(
    rows: Iterable[dict[str, float | str]], field: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grouped: dict[tuple[float, float], list[float]] = defaultdict(list)
    for row in rows:
        grouped[(float(row["x_m"]), float(row["y_m"]))].append(float(row[field]))
    ordered = sorted(grouped)
    return (
        np.asarray([key[0] for key in ordered]),
        np.asarray([key[1] for key in ordered]),
        np.asarray([max(grouped[key]) for key in ordered]),
    )


def plot_ffi_footprints(output_dir: Path) -> list[Path]:
    plt, colors, Line2D = _plotting()
    _style(plt)
    rows = _ffi_rows()
    panels = (
        ("observed", "Measured"),
        ("degali_replay", "DEGALI\nmeasured-wind replay"),
        ("slabx_research", "SLABx-LH2\nretained research path"),
    )
    bounds = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.5]
    cmap = plt.get_cmap("viridis")
    norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.7), sharex=True, sharey=True)
    for ax, (field, title) in zip(axes, panels):
        x_m, y_m, concentration = _mast_maxima(rows, field)
        triangulation = __import__("matplotlib.tri", fromlist=["Triangulation"]).Triangulation(x_m, y_m)
        image = ax.tricontourf(
            triangulation, concentration, levels=bounds, cmap=cmap, norm=norm,
            extend="max",
        )
        ax.scatter(
            x_m, y_m, c=concentration, cmap=cmap, norm=norm, s=56,
            edgecolor="white", linewidth=0.9, zorder=4,
        )
        if np.nanmin(concentration) < 4.0 < np.nanmax(concentration):
            ax.tricontour(
                triangulation, concentration, levels=[4.0], colors="#ef4444",
                linewidths=3.0, zorder=3,
            )
        ax.scatter([0.0], [0.0], marker="*", s=180, color="#111827", zorder=5)
        ax.annotate(
            "release", xy=(0.0, 0.0), xytext=(5.0, -4.0),
            arrowprops={"arrowstyle": "->", "color": "#111827"}, fontsize=9,
        )
        ax.set_title(title, pad=10, weight="bold")
        ax.set_xlim(-2.0, 104.0)
        ax.set_ylim(-26.0, 26.0)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Downwind distance x (m)")
        ax.grid(color="#d1d5db", linewidth=0.55, alpha=0.55)
    axes[0].set_ylabel("Crosswind distance y (m)")
    cbar = fig.colorbar(image, ax=axes, pad=0.018, fraction=0.032)
    cbar.set_label("H$_2$ concentration (vol%)")
    cbar.set_ticks(bounds)
    legend = [
        Line2D([0], [0], color="#ef4444", lw=3, label="4 vol% contour"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#64748b",
               markeredgecolor="white", markersize=7, label="sensor mast"),
    ]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.91),
               frameon=False, ncol=2)
    fig.suptitle("FFI Test 4 — flammable-footprint comparison")
    fig.text(
        0.5, 0.925,
        "275 s mean; maximum of three heights per mast. Linear interpolation only inside the sensor convex hull.",
        ha="center", color="#475569", fontsize=10,
    )
    fig.text(
        0.5, 0.02,
        "Visual diagnostic only: source mappings differ and this event was used during development; do not read as an independent model ranking.",
        ha="center", color="#7c2d12", fontsize=9.5,
    )
    fig.subplots_adjust(top=0.82, bottom=0.13, wspace=0.12)
    return _save(fig, output_dir, "ffi_test4_footprint_comparison")


def _aggregate_xy(points) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grouped: dict[tuple[float, float], list[float]] = defaultdict(list)
    for point in points:
        grouped[(float(point.x_m), float(point.y_m))].append(
            float(point.concentration_vol_percent)
        )
    ordered = sorted(grouped)
    return (
        np.asarray([key[0] for key in ordered]),
        np.asarray([key[1] for key in ordered]),
        np.asarray([float(np.mean(grouped[key])) for key in ordered]),
    )


def plot_hytunnel_conditions(output_dir: Path) -> list[Path]:
    from degali.validation.hytunnel_ps import read_hytunnel_concentration_workbook

    plt, colors, Line2D = _plotting()
    _style(plt)
    dataset = read_hytunnel_concentration_workbook(HYTUNNEL_WORKBOOK)
    wanted = (
        ("quiescent", 0.0, "Still air"),
        ("co-flow", 1.5, "Co-flow 1.5 m/s"),
        ("co-flow", 3.5, "Co-flow 3.5 m/s"),
        ("counter-flow", 3.5, "Counter-flow 3.5 m/s"),
    )
    selected = []
    for direction, speed, label in wanted:
        configuration = next(
            item for item in dataset.configurations
            if item.gas_label == "H2"
            and item.nozzle_diameter_mm == 4.0
            and item.nominal_mass_flow_g_s == 1.0
            and item.flow_direction == direction
            and math.isclose(item.nominal_ventilation_speed_m_s, speed)
        )
        points = [
            point for point in dataset.points
            if point.configuration_id == configuration.configuration_id
            and point.statistic == "average"
        ]
        selected.append((configuration, points, label))

    bounds = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 48.0]
    cmap = plt.get_cmap("turbo")
    norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 7.7), sharex=True, sharey=True)
    for ax, (configuration, points, label) in zip(axes.ravel(), selected):
        x_m, y_m, concentration = _aggregate_xy(points)
        triangulation = __import__("matplotlib.tri", fromlist=["Triangulation"]).Triangulation(x_m, y_m)
        image = ax.tricontourf(
            triangulation, concentration, levels=bounds, cmap=cmap, norm=norm,
            extend="max",
        )
        ax.scatter(x_m, y_m, s=8, c="#0f172a", alpha=0.5, zorder=4)
        if np.nanmin(concentration) < 4.0 < np.nanmax(concentration):
            ax.tricontour(
                triangulation, concentration, levels=[4.0], colors="white",
                linewidths=4.0, zorder=5,
            )
            ax.tricontour(
                triangulation, concentration, levels=[4.0], colors="#dc2626",
                linewidths=2.1, zorder=6,
            )
        ax.scatter([0.0], [0.0], marker=">", s=115, color="black", clip_on=False, zorder=7)
        if configuration.flow_direction == "co-flow":
            ax.annotate("air", xy=(1.75, 0.15), xytext=(1.25, 0.15),
                        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "white"},
                        color="white", weight="bold")
        elif configuration.flow_direction == "counter-flow":
            ax.annotate("air", xy=(1.25, 0.15), xytext=(1.75, 0.15),
                        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "white"},
                        color="white", weight="bold")
        ax.set_title(label, weight="bold")
        ax.set_xlim(0.0, 2.02)
        ax.set_ylim(-0.19, 0.19)
        ax.set_aspect("auto")
        ax.grid(color="white", linewidth=0.45, alpha=0.25)
    axes[1, 0].set_xlabel("Axial coordinate x (m)")
    axes[1, 1].set_xlabel("Axial coordinate x (m)")
    axes[0, 0].set_ylabel("Transverse coordinate y (m)")
    axes[1, 0].set_ylabel("Transverse coordinate y (m)")
    cbar = fig.colorbar(image, ax=axes, pad=0.018, fraction=0.028)
    cbar.set_label("Measured mean H$_2$ concentration (vol%)")
    cbar.set_ticks(bounds)
    legend = [Line2D([0], [0], color="#dc2626", lw=2.5, label="4 vol% contour")]
    fig.legend(handles=legend, loc="upper right", bbox_to_anchor=(0.89, 0.936), frameon=False)
    fig.suptitle("HyTunnel — ventilation reshapes the measured hydrogen plume")
    fig.text(
        0.5, 0.934,
        "Matched H$_2$ source: 4 mm nozzle, nominal 1 g/s. Published average matrices; no model fit.",
        ha="center", color="#475569", fontsize=10,
    )
    fig.text(
        0.5, 0.018,
        "Colour display is bounded at 0 vol%; small negative instrument offsets remain unchanged in the source audit.",
        ha="center", color="#475569", fontsize=9,
    )
    fig.subplots_adjust(top=0.84, bottom=0.1, wspace=0.12, hspace=0.27)
    return _save(fig, output_dir, "hytunnel_measured_plume_conditions")


def _side_view_samples(run, time_after_source_s: float):
    grouped: dict[tuple[float, float], list[float]] = defaultdict(list)
    for sensor in run.room_cloud_sensors:
        if sensor.position_raw_mm is None:
            continue
        raw_x, raw_height, _raw_lateral = sensor.position_raw_mm
        value = float(np.interp(
            run.source_t0_s + time_after_source_s,
            sensor.time_s,
            sensor.concentration_vol_percent,
        ))
        grouped[(-raw_x / 1000.0, raw_height / 1000.0)].append(value)
    ordered = sorted(grouped)
    return (
        np.asarray([key[0] for key in ordered]),
        np.asarray([key[1] for key in ordered]),
        np.asarray([float(np.mean(grouped[key])) for key in ordered]),
    )


def plot_usn_cloud_evolution(output_dir: Path) -> list[Path]:
    from degali.validation.usn_carpark import read_usn_carpark_mat

    plt, colors, Line2D = _plotting()
    _style(plt)
    run = read_usn_carpark_mat(USN_EXP20)
    times = (10.0, 30.0, 60.0, 360.0)
    bounds = np.arange(0.0, 10.0, 1.0)
    cmap = plt.get_cmap("plasma")
    norm = colors.BoundaryNorm(bounds, cmap.N, clip=True)
    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.9), sharex=True, sharey=True)
    for ax, time_s in zip(axes, times):
        x_m, height_m, concentration = _side_view_samples(run, time_s)
        triangulation = __import__("matplotlib.tri", fromlist=["Triangulation"]).Triangulation(x_m, height_m)
        image = ax.tricontourf(
            triangulation, concentration, levels=bounds, cmap=cmap, norm=norm,
            extend="max", alpha=0.9,
        )
        ax.scatter(
            x_m, height_m, c=concentration, cmap=cmap, norm=norm, s=85,
            edgecolor="white", linewidth=0.8, zorder=5,
        )
        if np.nanmin(concentration) < 4.0 < np.nanmax(concentration):
            ax.tricontour(
                triangulation, concentration, levels=[4.0], colors="white",
                linewidths=3.8, zorder=3,
            )
            ax.tricontour(
                triangulation, concentration, levels=[4.0], colors="#22c55e",
                linewidths=2.0, zorder=4,
            )
            status = f"{int(np.sum(concentration >= 4.0))}/{len(concentration)} positions ≥ 4%"
        else:
            status = "No averaged position ≥ 4%"
        ax.text(
            0.5, 0.04, status, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8.5, color="white", weight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#111827", "ec": "none", "alpha": 0.8},
        )
        ax.set_title(f"t = {time_s:.0f} s", weight="bold")
        ax.set_xlim(2.6, 12.0)
        ax.set_ylim(1.05, 2.29)
        ax.set_xticks([4.0, 6.0, 8.0, 10.0])
        ax.set_xlabel("Longitudinal position (m)")
        ax.grid(color="white", linewidth=0.45, alpha=0.25)
    axes[0].set_ylabel("Sensor height (m)")
    cbar = fig.colorbar(image, ax=axes, pad=0.035, fraction=0.025)
    cbar.set_label("Measured H$_2$ concentration (vol%)")
    cbar.set_ticks(bounds)
    legend = [Line2D([0], [0], color="#22c55e", lw=2.5, label="4 vol% contour")]
    fig.legend(handles=legend, loc="upper right", bbox_to_anchor=(0.89, 0.92), frameon=False)
    fig.suptitle("USN Experiment 20 — measured hydrogen-cloud evolution")
    fig.text(
        0.5, 0.91,
        "Room-array side view; values at repeated longitudinal/height coordinates are averaged across lateral offsets.",
        ha="center", color="#475569", fontsize=10,
    )
    fig.text(
        0.5, 0.018,
        "GH$_2$ confined-ventilation component evidence; this is not an LH$_2$ outdoor DEGALI validation result.",
        ha="center", color="#7c2d12", fontsize=9.2,
    )
    fig.subplots_adjust(top=0.8, bottom=0.14, wspace=0.1)
    return _save(fig, output_dir, "usn_exp20_cloud_evolution")


def _metrics(observed: np.ndarray, predicted: np.ndarray, threshold: float = 4.0) -> dict[str, float | int]:
    residual = predicted - observed
    observed_positive = observed >= threshold
    predicted_positive = predicted >= threshold
    tp = int(np.sum(observed_positive & predicted_positive))
    fn = int(np.sum(observed_positive & ~predicted_positive))
    fp = int(np.sum(~observed_positive & predicted_positive))
    tn = int(np.sum(~observed_positive & ~predicted_positive))
    return {
        "n": int(observed.size),
        "mae_vol_pct_point": float(np.mean(np.abs(residual))),
        "rmse_vol_pct_point": float(np.sqrt(np.mean(residual ** 2))),
        "bias_vol_pct_point": float(np.mean(residual)),
        "threshold_vol_pct": threshold,
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "sensitivity": None if tp + fn == 0 else tp / (tp + fn),
        "specificity": None if tn + fp == 0 else tn / (tn + fp),
        "accuracy": (tp + tn) / observed.size,
    }


def _error_records() -> tuple[list[dict], list[dict], dict]:
    ffi = _ffi_rows()
    observed = np.asarray([float(row["observed"]) for row in ffi])
    definitions = (
        ("DEGALI static", "degali_static", "development diagnostic"),
        ("DEGALI wind replay", "degali_replay", "development diagnostic"),
        ("SLABx native", "slabx_native", "native bridge"),
        ("SLABx retained research", "slabx_research", "retained research path"),
    )
    test4 = []
    for label, field, status in definitions:
        predicted = np.asarray([float(row[field]) for row in ffi])
        threshold_4 = _metrics(observed, predicted, threshold=4.0)
        threshold_2 = _metrics(observed, predicted, threshold=2.0)
        record = {
            "scope": "FFI Test 4, 30 sensors",
            "case": "Test 4",
            "model": label,
            "evidence_status": status,
            **threshold_4,
            "threshold_2_tp": threshold_2["tp"],
            "threshold_2_fn": threshold_2["fn"],
            "threshold_2_fp": threshold_2["fp"],
            "threshold_2_tn": threshold_2["tn"],
            "threshold_2_sensitivity": threshold_2["sensitivity"],
            "threshold_2_specificity": threshold_2["specificity"],
            "threshold_2_accuracy": threshold_2["accuracy"],
        }
        test4.append(record)

    arc_source = [row for row in _read_csv(FFI_ARCS) if row["dataset"] == "FFI"]
    arcs = []
    for case, chosen in (
        ("Tests 4+6", arc_source),
        ("Test 4", [row for row in arc_source if row["test_id"] == "4"]),
        ("Test 6", [row for row in arc_source if row["test_id"] == "6"]),
    ):
        obs = np.asarray([float(row["observed_vol_pct"]) for row in chosen])
        pred = np.asarray([float(row["corrected_vol_pct"]) for row in chosen])
        threshold_4 = _metrics(obs, pred, threshold=4.0)
        threshold_2 = _metrics(obs, pred, threshold=2.0)
        arcs.append({
            "scope": "FFI arc maxima",
            "case": case,
            "model": "DEGALI corrected baseline",
            "evidence_status": "far-field screen; only two releases",
            **threshold_4,
            "threshold_2_tp": threshold_2["tp"],
            "threshold_2_fn": threshold_2["fn"],
            "threshold_2_fp": threshold_2["fp"],
            "threshold_2_tn": threshold_2["tn"],
            "threshold_2_sensitivity": threshold_2["sensitivity"],
            "threshold_2_specificity": threshold_2["specificity"],
            "threshold_2_accuracy": threshold_2["accuracy"],
        })

    e35_source = [row for row in _read_csv(FFI_ARCS) if row["dataset"] == "E35"]
    e35 = []
    for label, field, status in (
        ("DEGALI corrected baseline", "corrected_vol_pct", "frozen primary validation"),
        ("Historical reconstruction", "historical_reconstruction_vol_pct", "reconstruction diagnostic"),
    ):
        obs = np.asarray([float(row["observed_vol_pct"]) for row in e35_source])
        pred = np.asarray([float(row[field]) for row in e35_source])
        threshold_4 = _metrics(obs, pred, threshold=4.0)
        threshold_2 = _metrics(obs, pred, threshold=2.0)
        e35.append({
            "scope": "PRESLHY E3.5 arc maxima",
            "case": "9 releases / 62 arcs",
            "model": label,
            "evidence_status": status + "; observed range 1.60–92.87 vol%",
            **threshold_4,
            "threshold_2_tp": threshold_2["tp"],
            "threshold_2_fn": threshold_2["fn"],
            "threshold_2_fp": threshold_2["fp"],
            "threshold_2_tn": threshold_2["tn"],
            "threshold_2_sensitivity": threshold_2["sensitivity"],
            "threshold_2_specificity": threshold_2["specificity"],
            "threshold_2_accuracy": threshold_2["accuracy"],
        })

    usn = json.loads(USN_AUDIT.read_text(encoding="utf-8"))
    return test4, arcs, e35, usn


def _metric_text(value: float | None) -> str:
    return "n/a" if value is None else f"{100.0 * value:.0f}%"


def plot_error_dashboard(output_dir: Path) -> tuple[list[Path], list[dict]]:
    plt, _colors, _Line2D = _plotting()
    _style(plt)
    test4, arcs, _e35, _usn = _error_records()
    palette = ["#94a3b8", "#2563eb", "#f59e0b", "#16a34a"]
    labels = [row["model"] for row in test4]
    short_labels = ["DEGALI\nstatic", "DEGALI\nwind replay", "SLABx\nnative", "SLABx\nresearch"]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 9.3))

    ax = axes[0, 0]
    width = 0.36
    mae = [row["mae_vol_pct_point"] for row in test4]
    rmse = [row["rmse_vol_pct_point"] for row in test4]
    ax.bar(x - width / 2, mae, width, color=palette, alpha=0.95, label="MAE")
    ax.bar(x + width / 2, rmse, width, color=palette, alpha=0.38,
           edgecolor=palette, linewidth=1.5, label="RMSE")
    for index, (left, right) in enumerate(zip(mae, rmse)):
        ax.text(index - width / 2, left + 0.06, f"{left:.3f}", ha="center", fontsize=8.5)
        ax.text(index + width / 2, right + 0.06, f"{right:.3f}", ha="center", fontsize=8.5)
    ax.set_xticks(x, short_labels)
    ax.set_ylabel("Error (vol% points)")
    ax.set_title("A. FFI Test 4 concentration error", loc="left", weight="bold")
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", color="#d1d5db", linewidth=0.65)

    ax = axes[0, 1]
    bias = [row["bias_vol_pct_point"] for row in test4]
    bars = ax.barh(x, bias, color=["#dc2626" if value > 0 else "#2563eb" for value in bias])
    ax.axvline(0.0, color="#111827", linewidth=1.0)
    ax.set_yticks(x, short_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Mean signed error (vol% points)")
    ax.set_title("B. Bias direction", loc="left", weight="bold")
    for bar, value in zip(bars, bias):
        ax.text(value + 0.025, bar.get_y() + bar.get_height() / 2,
                f"{value:+.3f}", ha="left", va="center", fontsize=9)
    ax.set_xlim(-0.52, 0.69)
    ax.grid(axis="x", color="#d1d5db", linewidth=0.65)

    ax = axes[1, 0]
    matrix = np.asarray([[row["tp"], row["fn"], row["fp"], row["tn"]] for row in test4])
    image = ax.imshow(matrix, cmap="Blues", aspect="auto", vmin=0, vmax=max(1, int(matrix.max())))
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center",
                    color="white" if matrix[i, j] > matrix.max() * 0.52 else "#111827", weight="bold")
        ax.text(4.0, i,
                f"sens {_metric_text(test4[i]['sensitivity'])}\nspec {_metric_text(test4[i]['specificity'])}",
                ha="left", va="center", fontsize=8.5)
    ax.set_xticks(range(4), ["TP", "FN", "FP", "TN"])
    ax.set_yticks(x, short_labels)
    ax.set_xlim(-0.5, 5.35)
    ax.set_title("C. 4 vol% decision counts", loc="left", weight="bold")
    ax.text(0.0, -0.14, "Only 3 measured positives: 100% sensitivity is not a broad confidence claim.",
            transform=ax.transAxes, color="#7c2d12", fontsize=9)

    ax = axes[1, 1]
    arc_rows = [row for row in _read_csv(FFI_ARCS) if row["dataset"] == "FFI"]
    positions = np.arange(len(arc_rows))
    residuals = [float(row["corrected_vol_pct"]) - float(row["observed_vol_pct"]) for row in arc_rows]
    arc_labels = [f"T{row['test_id']}\n{float(row['x_m']):.0f} m" for row in arc_rows]
    arc_colours = ["#2563eb" if value >= 0 else "#ef4444" for value in residuals]
    ax.axhline(0.0, color="#111827", linewidth=1.0)
    bars = ax.bar(positions, residuals, color=arc_colours)
    for bar, value in zip(bars, residuals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                value + (0.38 if value >= 0 else -0.38), f"{value:+.2f}",
                ha="center", va="bottom" if value >= 0 else "top", fontsize=8.5)
    ax.set_xticks(positions, arc_labels)
    ax.set_ylabel("Prediction − observation (vol% points)")
    ax.set_title("D. FFI far-field arc residuals", loc="left", weight="bold")
    ax.grid(axis="y", color="#d1d5db", linewidth=0.65)
    ax.set_ylim(-16.0, 1.25)
    overall = arcs[0]
    ax.text(
        0.98, 0.93,
        f"6-arc MAE {overall['mae_vol_pct_point']:.2f}\nRMSE {overall['rmse_vol_pct_point']:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "#cbd5e1"},
    )

    fig.suptitle("Hydrogen concentration error — one-page vol% summary")
    fig.text(
        0.5, 0.945,
        "vol% = concentration; vol% points = absolute difference between predicted and observed concentration",
        ha="center", color="#475569", fontsize=10,
    )
    fig.text(
        0.5, 0.018,
        "Sensor specification at a 4 vol% indication: Xensor ±3.00 vol%p (±3% full scale); Dräger ±0.08 vol%p (±2% of reading). Specifications are not event calibration errors.",
        ha="center", color="#334155", fontsize=9.2,
    )
    fig.subplots_adjust(top=0.88, bottom=0.105, hspace=0.42, wspace=0.31)
    paths = _save(fig, output_dir, "concentration_error_dashboard")
    return paths, test4 + arcs


def write_summary(output_dir: Path, records: list[dict], usn: dict) -> tuple[Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "concentration_error_summary.csv"
    fieldnames = [
        "scope", "case", "model", "evidence_status", "n",
        "mae_vol_pct_point", "rmse_vol_pct_point", "bias_vol_pct_point",
        "threshold_vol_pct", "tp", "fn", "fp", "tn",
        "sensitivity", "specificity", "accuracy",
        "threshold_2_tp", "threshold_2_fn", "threshold_2_fp", "threshold_2_tn",
        "threshold_2_sensitivity", "threshold_2_specificity", "threshold_2_accuracy",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    test4 = [row for row in records if row["scope"] == "FFI Test 4, 30 sensors"]
    arc = next(row for row in records if row["scope"] == "FFI arc maxima" and row["case"] == "Tests 4+6")
    e35 = [row for row in records if row["scope"] == "PRESLHY E3.5 arc maxima"]
    usn_events = []
    for event in usn["events"]:
        baseline = event["well_mixed_baseline"]
        usn_events.append({
            "experiment": event["experiment"],
            "observed_peak_vol_percent": baseline["observed_peak_vol_percent"],
            "well_mixed_peak_vol_percent": baseline["model_peak_vol_percent"],
            "peak_error_vol_pct_point": (
                baseline["model_peak_vol_percent"] - baseline["observed_peak_vol_percent"]
            ),
            "mae_vs_spatial_mean_vol_pct_point": baseline["mae_against_spatial_mean_vol_percent"],
            "rmse_vs_spatial_mean_vol_pct_point": baseline["rmse_against_spatial_mean_vol_percent"],
        })
    manifest = {
        "schema": "degali.concentration-shape-comparison-figures.v1",
        "date": "2026-10-10",
        "concentration_unit": "vol%",
        "error_unit": "vol% points",
        "sources": [
            {"path": str(path.relative_to(ROOT)), "sha256": _sha256(path)}
            for path in (FFI_REPLAY, FFI_COORDINATES, FFI_ARCS, USN_AUDIT, USN_EXP20, HYTUNNEL_WORKBOOK)
        ],
        "ffi_test4": test4,
        "ffi_arc_maxima_corrected": arc,
        "preslhy_e35_arc_maxima": e35,
        "usn_well_mixed_diagnostic": usn_events,
        "sensor_specification_at_4_vol_percent": {
            "Xensor_XEN_5320": {
                "declared_specification": "±3% of 0–100 vol% full scale",
                "absolute_width_vol_pct_point": 3.0,
                "specification_only_interval_vol_percent": [1.0, 7.0],
            },
            "Drager_Xam_5000_XXS_H2_HC": {
                "declared_specification": "±2% of reading over 0–4 vol%",
                "absolute_width_vol_pct_point": 0.08,
                "specification_only_interval_vol_percent": [3.92, 4.08],
            },
            "warning": "Declared specifications are not channel-level event calibration errors.",
        },
        "interpretation": {
            "ffi_test4": "Development-used single-event sensor comparison; not independent validation.",
            "ffi_arcs": "Far-field screen with only two releases and six arc maxima.",
            "usn": "No-fit 0-D well-mixed comparator, not the DEGALI free-field solver.",
            "hytunnel": "Experimental GH2 spatial component evidence; no DEGALI prediction is scored.",
        },
    }
    json_path = output_dir / "figure_manifest.json"
    json_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    def percent(value) -> str:
        return "n/a" if value is None else f"{100.0 * float(value):.1f}%"

    lines = [
        "# Hydrogen concentration error summary (vol% basis)",
        "",
        "`vol%` is concentration. `vol%p` is the absolute difference between two vol% values.",
        "The tables do not pool development, frozen validation and component-only evidence.",
        "",
        "## FFI Test 4 — 30 sensor rows (development-used event)",
        "",
        "| Model/path | MAE (vol%p) | RMSE (vol%p) | Bias (vol%p) | 2% sensitivity | 2% specificity | 4% sensitivity | 4% specificity |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in test4:
        lines.append(
            f"| {row['model']} | {row['mae_vol_pct_point']:.3f} | "
            f"{row['rmse_vol_pct_point']:.3f} | {row['bias_vol_pct_point']:+.3f} | "
            f"{percent(row['threshold_2_sensitivity'])} | {percent(row['threshold_2_specificity'])} | "
            f"{percent(row['sensitivity'])} | {percent(row['specificity'])} |"
        )
    lines.extend([
        "",
        "The 4 vol% result has only three measured-positive rows. Perfect sensitivity therefore has very low event-level breadth.",
        "",
        "## Frozen/far-field arc summaries", "",
        "| Scope | Model/path | n | MAE (vol%p) | RMSE (vol%p) | Bias (vol%p) |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for row in e35 + [arc]:
        lines.append(
            f"| {row['scope']} | {row['model']} | {row['n']} | "
            f"{row['mae_vol_pct_point']:.3f} | {row['rmse_vol_pct_point']:.3f} | "
            f"{row['bias_vol_pct_point']:+.3f} |"
        )
    lines.extend([
        "",
        "E3.5 contains much higher concentrations (observed 1.60–92.87 vol%), so its absolute vol%p errors must not be ranked directly against low-concentration FFI Test 4 rows.",
        "",
        "## USN no-fit well-mixed diagnostic (not DEGALI)", "",
        "| Experiment | Observed peak (vol%) | 0-D peak (vol%) | Peak error (vol%p) | MAE vs spatial mean (vol%p) | RMSE (vol%p) |",
        "|---:|---:|---:|---:|---:|---:|",
    ])
    for row in usn_events:
        lines.append(
            f"| {row['experiment']} | {row['observed_peak_vol_percent']:.3f} | "
            f"{row['well_mixed_peak_vol_percent']:.3f} | {row['peak_error_vol_pct_point']:+.3f} | "
            f"{row['mae_vs_spatial_mean_vol_pct_point']:.3f} | {row['rmse_vs_spatial_mean_vol_pct_point']:.3f} |"
        )
    lines.extend([
        "",
        "## Sensor specification boundary at a nominal 4 vol% indication", "",
        "| Instrument | Declared basis | Specification-only absolute width | Interval around 4 vol% |",
        "|---|---|---:|---:|",
        "| Xensor XEN-5320 | ±3% of 0–100 vol% full scale | ±3.00 vol%p | 1.00–7.00 vol% |",
        "| Dräger XXS H₂ HC | ±2% of reading over 0–4 vol% | ±0.08 vol%p | 3.92–4.08 vol% |",
        "",
        "These are declared specification widths, not measured event calibration errors and not model-residual corrections.",
    ])
    report_path = output_dir / "concentration_error_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return csv_path, json_path, report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()

    created: list[Path] = []
    created.extend(plot_ffi_footprints(output_dir))
    created.extend(plot_hytunnel_conditions(output_dir))
    created.extend(plot_usn_cloud_evolution(output_dir))
    dashboard_paths, records = plot_error_dashboard(output_dir)
    created.extend(dashboard_paths)
    _test4, _arcs, e35, usn = _error_records()
    created.extend(write_summary(output_dir, records + e35, usn))
    for path in created:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
