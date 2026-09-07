"""Merge disjoint PRESLHY field-result shards without rerunning solved cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def _merge_metrics(first: dict, second: dict) -> dict:
    n1, n2 = int(first["arcs"]), int(second["arcs"])
    total = n1 + n2
    if min(n1, n2) <= 0:
        raise ValueError("both field shards must contain concentration arcs")
    geometry = {
        name: math.nan
        for name in first["geometry"]
    }
    return {
        "arcs": total,
        # MG and VG definitions are exponential means of log ratio and its
        # square, so their disjoint pooled forms follow exactly from ln(MG)
        # and ln(VG). FAC2 and geometry are ordinary count-weighted means.
        "MG": math.exp(
            (n1 * math.log(first["MG"]) + n2 * math.log(second["MG"]))
            / total
        ),
        "VG": math.exp(
            (n1 * math.log(first["VG"]) + n2 * math.log(second["VG"]))
            / total
        ),
        "FAC2": (
            n1 * first["FAC2"] + n2 * second["FAC2"]
        ) / total,
        "geometry": geometry,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    first = json.loads(args.first.read_text(encoding="utf-8"))
    second = json.loads(args.second.read_text(encoding="utf-8"))
    for key in ("mode", "droplet_equilibrium_bound"):
        if first.get(key) != second.get(key):
            raise ValueError(f"field shards differ in {key}")
    for key, default in (("hydrogen_spin_isomer", "normal"),
                         ("source_energy_ledger", "legacy"),
                         ("phase_ambient_closure", "legacy_air_eos"),
                         ("downstream_thermodynamic_profile", "density"),
                         ("energy_quadrature", "unspecified")):
        if first.get(key, default) != second.get(key, default):
            raise ValueError(f"field shards differ in {key}")

    failed_first = {int(value) for value in first.get("failures", {})}
    failed_second = {int(value) for value in second.get("failures", {})}
    solved_first = set(map(int, first["selected_trials"])) - failed_first
    solved_second = set(map(int, second["selected_trials"])) - failed_second
    if solved_first & solved_second:
        raise ValueError("field shards contain overlapping solved trials")
    unresolved = {
        key: value for key, value in first.get("failures", {}).items()
        if int(key) not in solved_second
    }
    unresolved.update(second.get("failures", {}))

    vertical_total = int(first["vertical_rows"]) + int(second["vertical_rows"])
    baseline = _merge_metrics(first["baseline"], second["baseline"])
    candidate = _merge_metrics(first["candidate"], second["candidate"])
    for metrics_name in ("baseline", "candidate"):
        target = baseline if metrics_name == "baseline" else candidate
        a = first[metrics_name]["geometry"]
        b = second[metrics_name]["geometry"]
        for name in target["geometry"]:
            target["geometry"][name] = (
                first["vertical_rows"] * a[name]
                + second["vertical_rows"] * b[name]
            ) / vertical_total

    old_g, new_g = baseline["geometry"], candidate["geometry"]
    promoted = (
        not unresolved
        and abs(math.log(candidate["MG"])) < abs(math.log(baseline["MG"]))
        and candidate["VG"] < baseline["VG"]
        and candidate["FAC2"] >= baseline["FAC2"]
        and abs(new_g["width_ratio"] - 1.0)
        < abs(old_g["width_ratio"] - 1.0)
        and new_g["centre_mae"] <= old_g["centre_mae"]
    )
    interfaces = dict(first["interfaces"])
    interfaces.update(second["interfaces"])
    selected = sorted(solved_first | solved_second | set(map(int, unresolved)))
    payload = {
        "phase": "field",
        "mode": first["mode"],
        "droplet_equilibrium_bound": first["droplet_equilibrium_bound"],
        "hydrogen_spin_isomer": first.get("hydrogen_spin_isomer", "normal"),
        "source_energy_ledger": first.get("source_energy_ledger", "legacy"),
        "phase_ambient_closure": first.get("phase_ambient_closure", "legacy_air_eos"),
        "downstream_thermodynamic_profile": first.get("downstream_thermodynamic_profile", "density"),
        "energy_quadrature": first.get("energy_quadrature", "unspecified"),
        "aggregation": {
            "method": (
                "exact disjoint pooling: count-weighted ln(MG), ln(VG), "
                "FAC2 and geometry means"
            ),
            "source_files": [args.first.name, args.second.name],
            "source_sha256": {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (args.first, args.second)
            },
            "successful_trial_shards": [
                sorted(solved_first), sorted(solved_second)
            ],
        },
        "selected_trials": selected,
        "interfaces_accepted": (
            not unresolved
            and all(row["accepted"] for row in interfaces.values())
        ),
        "failures": unresolved,
        "baseline": baseline,
        "candidate": candidate,
        "vertical_rows": vertical_total,
        "promoted": promoted,
        "interfaces": interfaces,
    }
    for kind in ("baseline", "candidate"):
        for detail in ("pairs", "vertical_profiles"):
            if detail in first[kind] and detail in second[kind]:
                payload[kind][detail] = first[kind][detail] + second[kind][detail]
    args.output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"saved: {args.output}")
    print(json.dumps({
        "baseline": {k: v for k, v in baseline.items() if k not in {"pairs", "vertical_profiles"}},
        "candidate": {k: v for k, v in candidate.items() if k not in {"pairs", "vertical_profiles"}},
        "promoted": promoted,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
