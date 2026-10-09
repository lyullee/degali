"""Run the deterministic FFI source-state envelope without copying raw data.

Bounds use ``lower,nominal,upper``. A single value is an exact bound. The
result is a JSON residual map at the selected FFI sensor rows.
The report additionally contains corner-wise aggregate arc-max and
sensor-height observation-operator tables; raw source tables remain external.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

from degali.validation.ffi_source_state import (
    ffi_reference_provenance,
    FfiSourceState,
    ffi_source_state_envelope_report,
    run_ffi_source_state_envelope,
)
from degali.validation.spadeadam import load
from degali.addons.field_contracts import BoundedValue, CircularBoundedValue


def _bound(text: str | None, *, unit: str, name: str, circular: bool = False):
    if text is None:
        return None
    try:
        values = tuple(float(item.strip()) for item in text.split(","))
    except ValueError as error:
        raise ValueError(f"{name} must be one value or lower,nominal,upper") from error
    if len(values) == 1:
        lower = nominal = upper = values[0]
    elif len(values) == 3:
        lower, nominal, upper = values
    else:
        raise ValueError(f"{name} must be one value or lower,nominal,upper")
    source = f"CLI source-state envelope: {name}"
    if circular:
        return CircularBoundedValue(nominal, lower, upper, unit=unit, source=source)
    return BoundedValue(nominal, lower, upper, unit=unit, source=source)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", type=int, default=6)
    parser.add_argument("--reference-root", type=Path, default=None)
    parser.add_argument("--radius", type=float, default=None,
                        help="keep only readings within 0.5 m of this radius")
    parser.add_argument("--max-cases", type=int, default=64)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--no-corrections", action="store_true")
    for flag, dest, unit in (
        ("--rate-bounds", "rate_kg_s", "kg/s"),
        ("--orifice-bounds", "orifice_m", "m"),
        ("--release-height-bounds", "release_height_m", "m"),
        ("--storage-pressure-bounds", "storage_pressure_barg", "barg"),
        ("--wind-bounds", "wind_m_s", "m/s"),
        ("--ambient-temperature-bounds", "ambient_temperature_k", "K"),
        ("--humidity-bounds", "relative_humidity_pct", "%"),
        ("--ambient-pressure-bounds", "ambient_pressure_pa", "Pa"),
        ("--wind-reference-height-bounds", "wind_reference_height_m", "m"),
    ):
        parser.add_argument(flag, dest=dest, default=None,
                            help=f"{dest}: one value or lower,nominal,upper ({unit})")
    parser.add_argument("--wind-direction-bounds", default=None,
                        help="wind direction from in degrees: one value or lower,nominal,upper")
    args = parser.parse_args()

    trial = next(item for item in load(args.reference_root) if item.test == args.test)
    state = FfiSourceState.from_trial(trial)
    overrides = {}
    for dest, unit in (
        ("rate_kg_s", "kg/s"), ("orifice_m", "m"),
        ("release_height_m", "m"), ("storage_pressure_barg", "barg"),
        ("wind_m_s", "m/s"), ("ambient_temperature_k", "K"),
        ("relative_humidity_pct", "%"), ("ambient_pressure_pa", "Pa"),
        ("wind_reference_height_m", "m"),
    ):
        value = _bound(getattr(args, dest), unit=unit, name=dest)
        if value is not None:
            overrides[dest] = value
    direction = _bound(
        args.wind_direction_bounds, unit="deg", name="wind_direction_from_deg",
        circular=True,
    )
    if direction is not None:
        overrides["wind_direction_from_deg"] = direction
    state = replace(state, **overrides)

    readings = None
    if args.radius is not None:
        readings = tuple(
            reading for reading in trial.readings
            if abs(reading.radius - args.radius) < 0.5
        )
    result = run_ffi_source_state_envelope(
        trial, state, corrections=not args.no_corrections,
        readings=readings, max_cases=args.max_cases,
    )
    report = ffi_source_state_envelope_report(result)
    report["reference_provenance"] = ffi_reference_provenance(
        args.reference_root or "reference/spadeadam"
    )
    rendered = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        if args.output.exists():
            raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(rendered, end="")


if __name__ == "__main__":
    main()
