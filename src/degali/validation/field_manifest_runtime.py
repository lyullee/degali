"""Runtime gate for attaching a verified FieldEvidenceManifest to a solve."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..addons.field_evidence_manifest import read_field_evidence_manifest_json


def verify_runtime_manifest(
    path: str | Path, *, require_quantitative: bool = False,
) -> dict[str, Any]:
    """Read a saved execution manifest and expose its fail-closed readiness.

    The solver may produce a conditional physics witness when
    ``require_quantitative`` is false.  Production scoring can set it true; a
    conditional clock/calibration package then fails before the expensive
    solve rather than being mistaken for validated field evidence.
    """
    manifest_path = Path(path).resolve()
    manifest = read_field_evidence_manifest_json(manifest_path)
    quantitative = manifest.quantitative_readiness_record()
    if require_quantitative and quantitative.get("status") != "accepted":
        missing = quantitative.get("missing_requirements", [])
        raise ValueError(
            "quantitative FieldEvidenceManifest gate is not accepted: "
            + ", ".join(str(item) for item in missing)
        )
    return {
        "path": str(manifest_path),
        "manifest_id": manifest.manifest_id,
        "event_id": manifest.event_id,
        "scope": manifest.scope,
        "sensor_set_id": manifest.sensor_set_id,
        "operator_id": manifest.operator_id,
        "sensor_calibration_status": manifest.sensor_calibration_status,
        "quantitative_evidence_readiness": quantitative,
        "promotion_allowed": False,
    }


__all__ = ["verify_runtime_manifest"]
