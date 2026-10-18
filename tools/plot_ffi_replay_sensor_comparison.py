"""Plot conditional FFI O1 replay traces against measured sensor histories."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from degali.validation.open_channel_h2_obstacle import load_open_channel_h2_obstacle_benchmark

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "reference/hydrogen/open_channel_usn/obstacle_benchmark_registry_v1.json"
DATA_ROOT = ROOT / "reference/hydrogen/open_channel_usn"

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay", type=Path, default=ROOT / "outputs/ffi-o1-conditional-replay.json")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/ffi-o1-conditional-replay.png")
    args = parser.parse_args()
    replay = json.loads(args.replay.read_text(encoding="utf-8"))
    solver_label = (
        "DEGALI resolved buoyant replay"
        if "resolved" in str(replay.get("schema", ""))
        else "DEGALI fast conditional replay"
    )
    benchmark = load_open_channel_h2_obstacle_benchmark(REGISTRY, DATA_ROOT)
    run = next(item for item in benchmark.cases if item.article_case == "O1").runs[0]
    selected = [label for label in ("sensor 1", "sensor 2", "sensor 6", "sensor 11", "sensor 14", "sensor 24")
                if label in replay["predicted_sensor_trace_vol_percent"]]
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    time_pred = np.asarray(replay["predicted_sensor_time_s"], dtype=float)
    for axis, label in zip(axes.flat, selected):
        measured = np.asarray(run.sensor_percent[label], dtype=float)
        measured_time = np.asarray(run.sensor_time_s, dtype=float)
        keep = (measured_time >= 60.0) & (measured_time <= 120.0)
        axis.plot(measured_time[keep] - 60.0, measured[keep], color="#111827", lw=1.0, label="measured")
        axis.plot(time_pred, replay["predicted_sensor_trace_vol_percent"][label], color="#d97706", lw=1.2, label=solver_label)
        axis.axhline(4.0, color="#dc2626", ls="--", lw=0.8)
        axis.set_title(label)
        axis.set_ylabel("H₂ (vol%)")
        axis.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("Elapsed time in 60–120 s window (s)")
    axes[-1, 1].set_xlabel("Elapsed time in 60–120 s window (s)")
    axes[0, 0].legend(loc="upper right", fontsize=8)
    fig.suptitle(f"FFI O1 conditional obstacle replay: measured vs {solver_label}", weight="bold")
    fig.text(0.5, 0.01, "Conditional diagnostic only: obstacle beam geometry, clock alignment, and calibration gates remain withheld.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(json.dumps({"output": str(args.output.resolve()), "sensor_count": len(selected)}, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
