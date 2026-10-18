"""Project live TI021 runtime from append-only progress samples."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


SCHEMA = "degali.thorney-ti021-runtime-projection.v1"


def project_runtime(
    progress_path: Path,
    *,
    duration_s: float,
    checkpoint_interval_s: float,
) -> dict[str, object]:
    path = progress_path.resolve()
    raw = path.read_bytes()
    samples = [
        json.loads(line)
        for line in raw.decode("utf-8").splitlines()
        if line.strip()
    ]
    if len(samples) < 3:
        raise ValueError("at least three progress samples are required")
    times = [float(item["model_time_s"]) for item in samples]
    elapsed = [float(item["elapsed_wall_clock_s_this_process"]) for item in samples]
    if any(not math.isfinite(value) for value in times + elapsed):
        raise ValueError("progress samples must be finite")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("model times must increase strictly")
    if any(right <= left for left, right in zip(elapsed, elapsed[1:])):
        raise ValueError("wall-clock samples must increase strictly")
    duration = float(duration_s)
    checkpoint = float(checkpoint_interval_s)
    if not math.isfinite(duration) or duration <= times[-1]:
        raise ValueError("duration_s must exceed the last model time")
    if not math.isfinite(checkpoint) or checkpoint <= times[-1]:
        raise ValueError("checkpoint interval must exceed the last model time")
    model_delta = times[-1] - times[-2]
    wall_delta = elapsed[-1] - elapsed[-2]
    seconds_per_model_second = wall_delta / model_delta
    latest_time = times[-1]
    latest_wall = elapsed[-1]
    return {
        "schema": SCHEMA,
        "status": "live_projection_from_latest_progress_interval",
        "progress_path": str(path),
        "progress_prefix_sha256": hashlib.sha256(raw).hexdigest(),
        "sample_count": len(samples),
        "latest_model_time_s": latest_time,
        "latest_elapsed_wall_clock_s": latest_wall,
        "latest_completed_fraction": latest_time / duration,
        "latest_interval": {
            "model_delta_s": model_delta,
            "wall_delta_s": wall_delta,
            "wall_seconds_per_model_second": seconds_per_model_second,
        },
        "projection": {
            "first_checkpoint_model_time_s": checkpoint,
            "first_checkpoint_elapsed_wall_clock_s": latest_wall + (
                checkpoint - latest_time
            ) * seconds_per_model_second,
            "full_duration_s": duration,
            "full_duration_projected_elapsed_wall_clock_s": latest_wall + (
                duration - latest_time
            ) * seconds_per_model_second,
            "full_duration_projected_wall_clock_h": (
                latest_wall + (duration - latest_time) * seconds_per_model_second
            ) / 3600.0,
            "method": "latest positive progress interval; not a completion result",
        },
        "observation_value_fields_accessed": [],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--progress-jsonl", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, required=True)
    parser.add_argument("--checkpoint-interval-s", type=float, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    record = project_runtime(
        args.progress_jsonl,
        duration_s=args.duration_s,
        checkpoint_interval_s=args.checkpoint_interval_s,
    )
    destination = args.output_json.resolve()
    destination.parent.mkdir(parents=True, exist_ok=False)
    destination.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
