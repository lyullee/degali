"""Reproducible direct-vs-tabular field-LH2 flash timing check.

Run from the repository root with::

    python tools/benchmark_lh2_flash_table.py --count 1000 --nodes 161

It is a performance diagnostic, not a validation dataset or a replacement for
the regression comparison in ``tests/test_field_lh2.py``.
"""

from __future__ import annotations

import argparse
import json
import math
from time import perf_counter

from degali.addons.field_contracts import BoundedValue, ReleaseSource
from degali.addons.field_lh2 import (
    build_lh2_saturation_table_for_release,
    lh2_flash_source_from_release,
)


def _reference_release() -> ReleaseSource:
    return ReleaseSource(
        fluid="lh2",
        upstream_pressure=BoundedValue(0.4e6, unit="Pa"),
        upstream_temperature=BoundedValue(26.084, unit="K"),
        mass_flow_kg_s=BoundedValue(0.265, unit="kg/s"),
        opening_area_m2=BoundedValue(
            math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"
        ),
        discharge_coefficient=BoundedValue(0.8),
        liquid_fraction=BoundedValue(0.922),
        flash_model="homogeneous_equilibrium",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--nodes", type=int, default=161)
    args = parser.parse_args()
    if args.count < 1 or args.nodes < 4:
        parser.error("count must be positive and nodes must be at least four")

    release = _reference_release()
    reference = lh2_flash_source_from_release(release, ambient_temperature_k=293.15)
    started = perf_counter()
    table = build_lh2_saturation_table_for_release(release, nodes=args.nodes)
    table_build_s = perf_counter() - started

    started = perf_counter()
    for _ in range(args.count):
        lh2_flash_source_from_release(release, ambient_temperature_k=293.15)
    direct_s = perf_counter() - started

    started = perf_counter()
    tabular = None
    for _ in range(args.count):
        tabular = lh2_flash_source_from_release(
            release, ambient_temperature_k=293.15, property_table=table
        )
    tabular_s = perf_counter() - started
    assert tabular is not None
    quality_error = abs(
        tabular.flash.postflash_quality - reference.flash.postflash_quality
    ) / max(abs(reference.flash.postflash_quality), 1.0e-30)
    saving_per_case = (direct_s - tabular_s) / args.count
    payload = {
        "count": args.count,
        "nodes": args.nodes,
        "table_build_s": table_build_s,
        "direct_s": direct_s,
        "tabular_s": tabular_s,
        "per_case_speedup": direct_s / max(tabular_s, 1.0e-30),
        "quality_relative_error": quality_error,
        "break_even_cases": (
            table_build_s / saving_per_case if saving_per_case > 0.0 else None
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
