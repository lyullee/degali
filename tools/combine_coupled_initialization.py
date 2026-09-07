"""Combine the immutable two pilots and all five extension outcomes."""

import argparse
import hashlib
import json
from pathlib import Path
from audit_transverse_mixing import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("groups", type=Path, nargs="+")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite combined initialization evidence")
    rows, sha = [], {}
    for path in args.groups:
        group = json.loads(path.read_text(encoding="utf-8"))
        if not group["completed"] or group["phase"] not in {
            "simultaneous_shape_rate_initialization_pilots", "simultaneous_shape_rate_initialization_extension"}:
            raise ValueError("completed compatible initialization groups required")
        if group["selected_trials"] != [r["trial"] for r in group["rows"]]:
            raise ValueError("group selection is incomplete")
        for name, digest in group["sha256"].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest or (name in sha and sha[name] != digest):
                raise ValueError(f"inconsistent frozen dependency: {name}")
            sha[name] = digest
        rows.extend(group["rows"])
        sha[str(path.resolve().relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    if sorted(r["trial"] for r in rows) != [10, 11, 12, 22, 23, 24, 25]:
        raise ValueError("all seven trials required exactly once, including failures")
    sha[str(Path(__file__).resolve().relative_to(ROOT))] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows.sort(key=lambda r: r["trial"])
    output = dict(phase="combined_simultaneous_shape_rate_initialization", completed=True, rows=rows,
        selected_trials=[r["trial"] for r in rows], sha256=sha,
        passed_trials=[r["trial"] for r in rows if r.get("passed", False)],
        all_passed=all(r.get("passed", False) for r in rows),
        source_records=[str(p.resolve().relative_to(ROOT)) for p in args.groups],
        downstream_integrated=False, field_scored=False, promoted=False)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(dict(passed=output["passed_trials"], hashes=len(sha)), flush=True)


if __name__ == "__main__":
    main()
