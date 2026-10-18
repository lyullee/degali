"""Reproduce the topology-only cache timing for the TI021 obstacle grid.

This microbenchmark times connected-component labels for the configured fluid
and blocked-face geometry. It does not time a projection solve or the complete
transient and must not be reported as an end-to-end solver speedup.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from degali.addons.resolved_obstacle_flow_3d import (
    _circular_arc_blocked_faces,
    _solid_mask,
)
from degali.addons.variable_density_projection_3d import (
    _projection_component_labels,
    _projection_geometry,
)
try:
    from tools.audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )
except ModuleNotFoundError:  # Direct ``python tools/<script>.py`` execution.
    from audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROPERTY_AUDIT = (
    ROOT / "outputs" / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def benchmark(
    *, property_audit: Path, cell_size_m: float, time_step_s: float,
    corrector_tolerance: float, warm_repeats: int,
) -> dict[str, object]:
    properties = json.loads(property_audit.read_text(encoding="utf-8"))
    config, _, _ = build_configuration(
        properties,
        cell_size_m=cell_size_m,
        time_step_s=time_step_s,
        active_steps=1,
        corrector_tolerance=corrector_tolerance,
        circular_arc_representation="zero_thickness_surface",
    )
    x, y, z, solid = _solid_mask(config)
    blocked_faces = _circular_arc_blocked_faces(config, x, y, z)
    spacing = (float(cell_size_m),) * 3
    open_faces = ("x_max",)
    cache_key = (
        solid.shape,
        solid.tobytes(),
        spacing,
        open_faces,
        tuple(mask.tobytes() for mask in blocked_faces),
    )

    _projection_component_labels.cache_clear()
    started = time.perf_counter()
    component_count, labels = _projection_component_labels(*cache_key)
    cold_s = time.perf_counter() - started
    if labels.flags.writeable:
        raise RuntimeError("cached topology labels must be immutable")

    warm_times_s: list[float] = []
    for _ in range(warm_repeats):
        started = time.perf_counter()
        cached_count, cached_labels = _projection_component_labels(*cache_key)
        warm_times_s.append(time.perf_counter() - started)
        if cached_count != component_count or cached_labels is not labels:
            raise RuntimeError("cache hit did not reuse the committed topology")

    # The production projection path serializes these masks to bytes before
    # entering the LRU caches. Time that key construction plus both cached
    # geometry lookups separately; the prebuilt-key hit above is not the
    # complete steady-state cache-access cost.
    full_lookup_times_s: list[float] = []
    for _ in range(warm_repeats):
        started = time.perf_counter()
        current_solid_bytes = solid.tobytes()
        current_blocked_bytes = tuple(mask.tobytes() for mask in blocked_faces)
        _projection_geometry(
            solid.shape, current_solid_bytes, spacing, open_faces,
            current_blocked_bytes,
        )
        current_count, current_labels = _projection_component_labels(
            solid.shape, current_solid_bytes, spacing, open_faces,
            current_blocked_bytes,
        )
        full_lookup_times_s.append(time.perf_counter() - started)
        if current_count != component_count or current_labels is not labels:
            raise RuntimeError("full cache path did not reuse the committed topology")

    return {
        "schema": "degali.ti021-projection-component-cache-benchmark.v1",
        "scope": (
            "topology and steady-state cache-access microbenchmark; not the "
            "linear projection solve or complete transient runtime"
        ),
        "geometry": {
            "cell_size_m": float(cell_size_m),
            "time_step_s": float(time_step_s),
            "grid_cells_xyz": list(config.cells),
            "fluid_cell_count": int((~solid).sum()),
            "solid_cell_count": int(solid.sum()),
            "blocked_face_counts_xyz": [
                int(mask.sum()) for mask in blocked_faces
            ],
            "open_boundary_faces": list(open_faces),
            "connected_fluid_components": int(component_count),
        },
        "timing": {
            "cold_topology_build_s": float(cold_s),
            "warm_cache_hit_mean_s": (
                float(sum(warm_times_s) / len(warm_times_s))
                if warm_times_s else None
            ),
            "warm_repeats": int(warm_repeats),
            "warm_times_s": warm_times_s,
            "full_cache_path_including_mask_byte_keys_mean_s": (
                float(sum(full_lookup_times_s) / len(full_lookup_times_s))
                if full_lookup_times_s else None
            ),
            "full_cache_path_times_s": full_lookup_times_s,
        },
        "provenance": {
            "property_audit_path": str(property_audit.resolve()),
            "property_audit_sha256": _sha256(property_audit),
            "benchmark_script_sha256": _sha256(Path(__file__).resolve()),
            "projection_source_sha256": _sha256(
                ROOT / "src" / "degali" / "addons"
                / "variable_density_projection_3d.py"
            ),
            "obstacle_flow_source_sha256": _sha256(
                ROOT / "src" / "degali" / "addons"
                / "resolved_obstacle_flow_3d.py"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--property-audit", type=Path, default=DEFAULT_PROPERTY_AUDIT)
    parser.add_argument("--cell-size-m", type=float, default=1.25)
    parser.add_argument("--time-step-s", type=float, default=0.02)
    parser.add_argument("--corrector-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--warm-repeats", type=int, default=10)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "outputs"
        / "ti021-projection-component-cache-benchmark-2026-10-12.json",
    )
    args = parser.parse_args()
    if args.warm_repeats < 1:
        parser.error("--warm-repeats must be at least one")
    record = benchmark(
        property_audit=args.property_audit.resolve(),
        cell_size_m=args.cell_size_m,
        time_step_s=args.time_step_s,
        corrector_tolerance=args.corrector_tolerance,
        warm_repeats=args.warm_repeats,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
