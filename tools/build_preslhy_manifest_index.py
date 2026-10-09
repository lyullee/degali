"""Build a reviewer-facing index of conditional PRESLHY event manifests.

The index carries only derived metadata and hashes.  It deliberately omits the
third-party workbooks and absolute workstation paths, so it can be included in
an IJHE hand-off without redistributing the public archive or leaking local
filesystem layout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.preslhy-conditional-manifest-index.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.name


def _record(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    field = payload["field_evidence_manifest"]
    readiness = field["evidence_readiness"]
    return {
        "event_id": field["event_id"],
        "manifest_id": field["manifest_id"],
        "dataset_id": field["dataset_id"],
        "scope": field["scope"],
        "temporal_operator_id": field["temporal_operator_id"],
        "sensor_calibration_status": field["sensor_calibration_status"],
        "observed_row_count": field["observed_row_count"],
        "manifest_status": readiness["status"],
        "promotion_allowed": field["promotion_allowed"],
        "missing_requirements": list(readiness.get("missing_requirements", [])),
        "manifest_path": _relative(path),
        "manifest_sha256": _sha256(path),
    }


def build(output: Path, manifests: list[Path]) -> dict[str, object]:
    if not manifests:
        raise ValueError("at least one manifest is required")
    records = [_record(path.resolve()) for path in manifests]
    records.sort(key=lambda item: str(item["event_id"]))
    payload = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "public_dataset_doi": "10.35097/1481",
            "custodian_record_doi": "10.5445/IR/1000136281",
            "raw_workbooks_redistributed": False,
        },
        "promotion_rule": "All listed bundles are conditional traceability inputs; none is a promoted field-validation score.",
        "events": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    args = parser.parse_args(argv)
    payload = build(args.output, [path.resolve() for path in args.manifest])
    print(json.dumps({"output": str(args.output), "events": len(payload["events"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
