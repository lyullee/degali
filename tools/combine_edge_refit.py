"""Combine disjoint, completed refit groups without rewriting their evidence."""

import argparse
import hashlib
import json
from pathlib import Path
from audit_transverse_mixing import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite combined evidence")
    rows, sources, hashes = [], [], {}
    for path in args.inputs:
        record = json.loads(path.read_text(encoding="utf-8"))
        if not record["completed"] or record["phase"] != "face_ray_phase_split_conservative_edge_refit":
            raise ValueError("only completed fixed-transport conservative refit groups are accepted")
        if record["selected_trials"] != [r["trial"] for r in record["rows"]]:
            raise ValueError("group does not contain exactly its selected trials")
        for relative, digest in record["sha256"].items():
            if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
                raise ValueError(f"source/input no longer matches: {relative}")
            if relative in hashes and hashes[relative] != digest:
                raise ValueError("groups used different source versions")
            hashes[relative] = digest
        rows.extend(record["rows"])
        sources.append(str(path.resolve().relative_to(ROOT)))
    if sorted(r["trial"] for r in rows) != [10, 11, 12, 22, 23, 24, 25]:
        raise ValueError("combined audit requires each of the seven frozen trials exactly once")
    rows.sort(key=lambda r: r["trial"])
    for path in [*[ROOT/s for s in sources], Path(__file__).resolve(), ROOT/"tests/test_edge_refit_retraction.py"]:
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    passed = [r["trial"] for r in rows if r["passed"]]
    payload = dict(phase="combined_face_ray_phase_split_conservative_edge_refit", completed=True,
        input_records=sources, rows=rows, passed_trials=passed,
        failed_trials=[r["trial"] for r in rows if not r["passed"]], all_passed=len(passed)==7,
        all_reference_quadrature_passed=all(r.get("quadrature_reference_passed", False) for r in rows),
        all_moments_passed=all(r.get("moment_passed", False) for r in rows),
        transport_reconstructed=False, shape_ode_solved=False, field_scored=False, promoted=False,
        sha256=hashes,
        limitations=["Passed means fixed-transport cross-section gates only, never downstream or observation validation.",
            "Unchanged earlier failures and disjoint current group records remain authoritative evidence."])
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    print(f"combined: passed={passed}, failed={payload['failed_trials']}, hashes={len(hashes)}")


if __name__ == "__main__":
    main()
