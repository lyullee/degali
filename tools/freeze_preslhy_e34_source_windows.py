"""Freeze PRESLHY E3.4 release windows from corrected pool mass only.

This command intentionally opens only the first, source-history worksheet.  It
does not read hydrogen concentration, wind, trigger, temperature, model output
or an error metric.  Every qualifying post-fill decay returned by one declared
policy is retained, preventing a later analyst from choosing only favourable
sensor comparisons.

The frozen windows are evidence for the *timing* of candidate source events.
They do not promote liquid-mass loss to a measured atmospheric vapour flux.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from degali.validation.preslhy_e34 import (
    E34SourceWindowPolicy,
    read_e34_pool_history,
    source_only_evaporation_candidates,
)


SCHEMA = "degali.preslhy-e34-source-window-freeze.v1"
POLICY_ID = "preslhy-e34-corrected-pool-mass-all-post-fill-decays-v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def freeze_workbook(
    workbook_path: str | Path,
    *,
    policy: E34SourceWindowPolicy = E34SourceWindowPolicy(),
) -> dict[str, Any]:
    """Freeze all mass-only source windows in one public workbook."""

    workbook = Path(workbook_path).resolve()
    if not workbook.is_file():
        raise FileNotFoundError(workbook)
    history = read_e34_pool_history(workbook)
    candidates = source_only_evaporation_candidates(history, policy=policy)
    return {
        "trial_id": workbook.stem.removesuffix("-Final"),
        "workbook_path": str(workbook),
        "workbook_sha256": _sha256(workbook),
        "substrate": history.substrate,
        "source_sample_count": int(history.time_s.size),
        "source_time_range_s": [
            float(history.time_s[0]), float(history.time_s[-1]),
        ],
        "window_count": len(candidates),
        "windows": [asdict(candidate) for candidate in candidates],
    }


def build_source_window_freeze(
    data_directory: str | Path,
    output_directory: str | Path,
    *,
    policy: E34SourceWindowPolicy = E34SourceWindowPolicy(),
) -> dict[str, Any]:
    """Write one immutable mass-only selection record and review plot."""

    data = Path(data_directory).resolve()
    output = Path(output_directory).resolve()
    workbooks = sorted(data.glob("*Wind-Final.xlsx"))
    if not workbooks:
        raise ValueError("no E3.4 wind workbooks found")
    if output.exists():
        raise FileExistsError(output)

    trials = [freeze_workbook(path, policy=policy) for path in workbooks]
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY Experiment series 3.4",
        "source_doi": "10.35097/1319",
        "policy_id": POLICY_ID,
        "policy": asdict(policy),
        "selection_frozen_before_model_scoring": True,
        "all_qualifying_windows_retained": True,
        "selection_input_channels": ["Sync. Time [s]", "m(LH2) [g]"],
        "selection_inputs_explicitly_unused": [
            "hydrogen concentration",
            "wind",
            "trigger",
            "temperature",
            "model prediction",
            "validation error metric",
        ],
        "source_timing_evidence_status": "conditional",
        "atmospheric_vapour_flux_qualified": False,
        "qualification_holdbacks": [
            "liquid_mass_loss_is_not_a_direct_measurement_of_atmospheric_vapour_flux",
            "porous_substrate_mass_artefacts_require_review_where_applicable",
        ],
        "raw_measurements_copied": False,
        "trial_count": len(trials),
        "total_window_count": sum(trial["window_count"] for trial in trials),
        "review_plot": "source-window-freeze-review.png",
        "trials": trials,
    }
    output.mkdir(parents=True)
    plot_source_window_freeze(record, output / record["review_plot"])
    (output / "source-window-freeze.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def plot_source_window_freeze(
    record: dict[str, Any], output_path: str | Path,
) -> None:
    """Render a source-only review plot; no sensor or wind data are opened."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    trials = record["trials"]
    figure, axes = plt.subplots(
        len(trials), 1, figsize=(11.5, 2.5 * len(trials)), squeeze=False,
        constrained_layout=True,
    )
    for row, trial in enumerate(trials):
        axis = axes[row, 0]
        history = read_e34_pool_history(trial["workbook_path"])
        axis.plot(
            history.time_s, history.liquid_mass_kg,
            color="#1d4ed8", linewidth=0.8, label="corrected pool mass",
        )
        for index, window in enumerate(trial["windows"]):
            label = "all frozen post-fill windows" if index == 0 else None
            axis.axvspan(
                window["start_s"], window["end_s"],
                color="#f59e0b", alpha=0.22, label=label,
            )
            axis.axvline(window["start_s"], color="#b45309", linewidth=0.7)
        axis.set_title(
            f"{trial['trial_id']} — {trial['window_count']} source-only windows"
        )
        axis.set_ylabel("m(LH2) (kg)")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, loc="upper right")
    axes[-1, 0].set_xlabel("synchronized source clock (s)")
    figure.suptitle(
        "PRESLHY E3.4 source-window freeze — concentration and wind unused",
        fontsize=13,
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_directory", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record = build_source_window_freeze(
        args.data_directory, args.output_directory,
    )
    print(json.dumps({
        "schema": record["schema"],
        "policy_id": record["policy_id"],
        "trial_count": record["trial_count"],
        "total_window_count": record["total_window_count"],
        "source_timing_evidence_status": record[
            "source_timing_evidence_status"
        ],
        "atmospheric_vapour_flux_qualified": record[
            "atmospheric_vapour_flux_qualified"
        ],
        "output_directory": str(args.output_directory.resolve()),
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
