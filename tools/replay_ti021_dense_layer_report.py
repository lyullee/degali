"""Replay a TI021 dense-layer report with its recorded physical parameters."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from pathlib import Path

from score_ti021_dense_layer_hybrid import run


ALIASES = {
    "h": "horizontal_cell_size_m",
    "dt": "max_dt_s",
    "source_dt": "source_time_step_s",
    "entrainment": "dense_entrainment_model",
}


def replay_parameters(
    record: dict,
    *,
    component_diagnostics: bool = False,
    overrides: dict | None = None,
) -> dict:
    """Recover every recorded ``run`` parameter, failing on an unknown gap."""
    parameters = {}
    for name, parameter in inspect.signature(run).parameters.items():
        if name == "output":
            continue
        key = ALIASES.get(name, name)
        if name == "vertical_observation_component_diagnostics":
            parameters[name] = bool(component_diagnostics)
        elif key in record:
            parameters[name] = record[key]
        elif parameter.default is not inspect.Parameter.empty:
            parameters[name] = parameter.default
        else:
            raise ValueError(f"reference report does not record required parameter {name!r}")
    for name, value in (overrides or {}).items():
        if name not in parameters:
            raise ValueError(f"override is not a run() parameter: {name!r}")
        parameters[name] = value
    return parameters


def _parse_overrides(items: list[str]) -> dict:
    result = {}
    for item in items:
        name, separator, raw = item.partition("=")
        if not separator or not name:
            raise ValueError("--set requires NAME=JSON_VALUE")
        if name in result:
            raise ValueError(f"duplicate override {name!r}")
        try:
            result[name] = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"override {name!r} is not valid JSON") from exc
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--component-diagnostics", action="store_true")
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="NAME=JSON_VALUE",
        help="override a named run() parameter with an explicit JSON value",
    )
    args = parser.parse_args()
    reference = args.reference.resolve()
    output = args.output.resolve()
    if output.exists() or output.with_suffix(".npz").exists():
        raise FileExistsError("replay output already exists")
    record = json.loads(reference.read_text(encoding="utf-8"))
    overrides = _parse_overrides(args.set)
    parameters = replay_parameters(
        record,
        component_diagnostics=args.component_diagnostics,
        overrides=overrides,
    )
    run(output, **parameters)
    replay = json.loads(output.read_text(encoding="utf-8"))
    replay["replay_reference"] = {
        "path": str(reference),
        "sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        "parameter_source": "recorded run() fields plus explicit alias map",
        "component_diagnostics_override": bool(args.component_diagnostics),
        "parameter_overrides": overrides,
    }
    replay["reproduction_command"] = (
        f"python tools/replay_ti021_dense_layer_report.py --reference {reference} "
        f"--output {output}"
        + (" --component-diagnostics" if args.component_diagnostics else "")
        + "".join(
            f" --set {name}={json.dumps(value, separators=(',', ':'))}"
            for name, value in overrides.items()
        )
    )
    output.write_text(json.dumps(replay, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "dose_mape_percent": replay["dose_mape_percent"],
        "runtime_s": replay["prediction_runtime_s"],
    }, indent=2))


if __name__ == "__main__":
    main()
