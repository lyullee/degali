"""Explicit cross-file evidence manifest for LH2 field validation.

The readiness audit deliberately stops at channel discovery.  This module is
the next, manual boundary: an accountable reviewer selects one file for each
channel and declares the event/clock/operator identities that join them.  The
manifest fingerprints every selected file and can feed the existing
``FieldValidationEvidence`` record after its hashes are rechecked.  It never
marks a model or a design basis as promoted.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Literal, Mapping

from .field_contracts import FieldValidationEvidence
from .field_evidence_audit import (
    FieldEvidenceAudit,
    _CHANNELS,
    _json_object_pairs,
    _sha256_file,
)


FIELD_EVIDENCE_MANIFEST_SCHEMA = "degali.field-evidence-manifest.v1"
FIELD_EVIDENCE_MANIFEST_EXECUTION_SCHEMA = (
    "degali.field-evidence-manifest-execution.v1"
)
FIELD_EVIDENCE_MANIFEST_GATE_CODES = frozenset({
    "manifest_not_promoted",
    "manifest_files_verified",
    "manifest_metadata_complete",
    "manifest_metadata_conditional",
    "manifest_files_withheld",
})

_UNSPECIFIED_METADATA_VALUES = frozenset({
    "", "unspecified", "unknown", "not_recorded", "not_available",
    "missing", "n/a", "none", "legacy-unspecified",
})

SENSOR_CALIBRATION_STATUSES = frozenset({
    "certified",
    "specification_only",
    "missing",
    "legacy-unspecified",
})


def _string(value: object, name: str, *, reject_unspecified: bool = False) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    result = value.strip()
    if reject_unspecified and result.lower() == "unspecified":
        raise ValueError(f"{name} must not be unspecified")
    return result


def _sha(value: object, name: str) -> str:
    result = _string(value, name).lower()
    if len(result) != 64 or any(char not in "0123456789abcdef" for char in result):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return result


def _keys(
    value: Mapping[str, object], name: str, *, required: set[str],
    optional: set[str] | None = None,
) -> None:
    optional = set() if optional is None else optional
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required - optional)
    if missing or unknown:
        parts: list[str] = []
        if missing:
            parts.append("missing=" + ", ".join(missing))
        if unknown:
            parts.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(parts))


@dataclass(frozen=True)
class FieldEvidenceManifestArtifact:
    """One selected, hash-pinned evidence file in a cross-file manifest."""

    channel: Literal[
        "source_boundary", "weather", "obstacle_geometry",
        "receptor_observations", "common_clock",
    ]
    path: str
    sha256: str

    def __post_init__(self) -> None:
        if self.channel not in (*_CHANNELS, "sensor_registry"):
            raise ValueError("manifest artifact channel is unsupported")
        path = Path(self.path)
        if not path.is_absolute() or not self.path.strip():
            raise ValueError("manifest artifact path must be absolute")
        normalized_sha = _sha(self.sha256, "manifest artifact sha256")
        if self.sha256 != normalized_sha:
            raise ValueError("manifest artifact sha256 must be lowercase")

    def as_record(self) -> dict[str, str]:
        return {
            "channel": self.channel,
            "path": self.path,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class FieldEvidenceManifest:
    """Manual event join for the five field-evidence channels.

    This is an accountable qualification input, not validation itself.  The
    optional sensor-registry artifact pins calibration/response metadata that
    the receptor file alone cannot establish.  Missing registry or operator
    identity is retained as a conditional package status.
    A registry artifact with an explicit operator is still conditional unless
    ``sensor_calibration_status`` is ``certified``; nominal specifications are
    deliberately not treated as event calibration.  The
    ``promotion_allowed`` value is intentionally hard-coded to ``False`` in
    :meth:`as_record` and is never accepted as a constructor field.
    """

    manifest_id: str
    event_id: str
    audit_root: str
    channel_artifacts: tuple[FieldEvidenceManifestArtifact, ...]
    dataset_id: str
    observed_dataset_path: str
    observed_dataset_sha256: str
    observed_row_count: int
    source_boundary_id: str
    weather_id: str
    obstacle_geometry_id: str
    receptor_geometry_id: str
    temporal_operator_id: str
    common_clock_id: str
    scope: Literal["lh2_free_field", "lh2_obstacle_transport"] = (
        "lh2_obstacle_transport"
    )
    # A receptor file may contain many sensor IDs.  ``sensor_set_id`` points
    # to the registry/calibration export that defines those per-row IDs.
    # The legacy sentinel keeps older manifests readable, but makes their
    # readiness explicitly conditional.
    sensor_set_id: str = "legacy-unspecified"
    # Accountable person/system that assembled and reconciled the package.
    operator_id: str = "legacy-unspecified"
    # Optional for legacy manifests.  New packages should pin the registry or
    # calibration export that defines ``sensor_set_id``.
    sensor_registry_artifact: FieldEvidenceManifestArtifact | None = None
    # Distinguish an event-level calibration certificate from a nominal
    # instrument specification.  A registry file alone cannot establish the
    # former, so non-certified statuses remain conditional.
    sensor_calibration_status: Literal[
        "certified", "specification_only", "missing", "legacy-unspecified"
    ] = "legacy-unspecified"

    def __post_init__(self) -> None:
        for name, value in {
            "manifest_id": self.manifest_id,
            "event_id": self.event_id,
            "audit_root": self.audit_root,
            "dataset_id": self.dataset_id,
            "observed_dataset_path": self.observed_dataset_path,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "obstacle_geometry_id": self.obstacle_geometry_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
            "common_clock_id": self.common_clock_id,
        }.items():
            _string(value, f"manifest {name}", reject_unspecified=True)
        _string(self.sensor_set_id, "manifest sensor_set_id")
        _string(self.operator_id, "manifest operator_id")
        if self.sensor_calibration_status not in SENSOR_CALIBRATION_STATUSES:
            raise ValueError(
                "manifest sensor_calibration_status is unsupported"
            )
        if self.sensor_registry_artifact is not None:
            if not isinstance(
                self.sensor_registry_artifact, FieldEvidenceManifestArtifact
            ):
                raise TypeError(
                    "manifest sensor_registry_artifact must be a typed artifact"
                )
            if self.sensor_registry_artifact.channel != "sensor_registry":
                raise ValueError(
                    "manifest sensor_registry_artifact must use the sensor_registry channel"
                )
        root = Path(self.audit_root)
        if not root.is_absolute():
            raise ValueError("manifest audit_root must be absolute")
        observed_path = Path(self.observed_dataset_path)
        if not observed_path.is_absolute():
            raise ValueError("manifest observed_dataset_path must be absolute")
        if isinstance(self.observed_row_count, bool) or not isinstance(
            self.observed_row_count, int
        ) or self.observed_row_count <= 0:
            raise ValueError("manifest observed_row_count must be a positive integer")
        if self.scope not in {"lh2_free_field", "lh2_obstacle_transport"}:
            raise ValueError("manifest scope is unsupported")
        if (
            self.scope == "lh2_obstacle_transport"
            and self.obstacle_geometry_id.strip().lower() in {"none", "n/a", "not_applicable"}
        ):
            raise ValueError(
                "obstacle-transport manifest requires an obstacle geometry ID"
            )
        if not isinstance(self.channel_artifacts, tuple):
            raise TypeError("manifest channel_artifacts must be a tuple")
        if len(self.channel_artifacts) != len(_CHANNELS):
            raise ValueError("manifest must select exactly one artifact per channel")
        if not all(isinstance(item, FieldEvidenceManifestArtifact) for item in self.channel_artifacts):
            raise TypeError("manifest channel_artifacts must contain typed artifacts")
        channels = tuple(item.channel for item in self.channel_artifacts)
        if set(channels) != set(_CHANNELS) or len(set(channels)) != len(channels):
            raise ValueError("manifest channel_artifacts must contain each channel once")
        receptor_artifact = next(
            item for item in self.channel_artifacts
            if item.channel == "receptor_observations"
        )
        if Path(self.observed_dataset_path).resolve() != Path(receptor_artifact.path).resolve():
            raise ValueError(
                "manifest observed_dataset_path must equal the receptor-observation artifact"
            )
        if self.observed_dataset_sha256 != receptor_artifact.sha256:
            raise ValueError(
                "manifest observed_dataset_sha256 must match the receptor artifact"
            )
        normalized_sha = _sha(
            self.observed_dataset_sha256, "manifest observed_dataset_sha256"
        )
        if self.observed_dataset_sha256 != normalized_sha:
            raise ValueError("manifest observed_dataset_sha256 must be lowercase")

    @classmethod
    def from_audit(
        cls,
        audit: FieldEvidenceAudit,
        *,
        manifest_id: str,
        event_id: str,
        selected_paths: Mapping[str, str],
        dataset_id: str,
        observed_row_count: int,
        source_boundary_id: str,
        weather_id: str,
        obstacle_geometry_id: str,
        receptor_geometry_id: str,
        temporal_operator_id: str,
        common_clock_id: str,
        scope: Literal["lh2_free_field", "lh2_obstacle_transport"] = (
            "lh2_obstacle_transport"
        ),
        sensor_set_id: str = "legacy-unspecified",
        operator_id: str = "legacy-unspecified",
        sensor_registry_path: str | None = None,
        sensor_calibration_status: Literal[
            "certified", "specification_only", "missing", "legacy-unspecified"
        ] = "legacy-unspecified",
    ) -> "FieldEvidenceManifest":
        """Select and hash one candidate path per audit channel.

        The audit may be ``candidate_complete`` even when channels are split
        across files.  This constructor therefore requires an explicit path
        selection and never treats discovery order as event identity.
        """
        if not isinstance(audit, FieldEvidenceAudit):
            raise TypeError("audit must be a FieldEvidenceAudit")
        if audit.status != "candidate_complete" or not audit.scan_complete:
            raise ValueError(
                "manifest requires a candidate_complete audit with detected channels"
            )
        if not isinstance(selected_paths, Mapping):
            raise TypeError("selected_paths must be a mapping")
        selected = {
            _string(channel, "selected_paths channel"): _string(path, "selected_paths path")
            for channel, path in selected_paths.items()
        }
        if set(selected) != set(_CHANNELS):
            raise ValueError("selected_paths must contain exactly the five audit channels")
        path_by_channel = dict(audit.channel_paths)
        candidate_by_path = {candidate.path: candidate for candidate in audit.candidates}
        artifacts: list[FieldEvidenceManifestArtifact] = []
        for channel in _CHANNELS:
            path = str(Path(selected[channel]).resolve())
            if path not in path_by_channel[channel]:
                raise ValueError(
                    f"selected path is not an audited {channel} candidate: {path}"
                )
            candidate = candidate_by_path.get(path)
            if candidate is None or candidate.sha256 is None:
                raise ValueError(f"selected candidate lacks SHA-256 provenance: {path}")
            artifacts.append(
                FieldEvidenceManifestArtifact(
                    channel=channel, path=path, sha256=candidate.sha256,
                )
            )
        receptor = next(item for item in artifacts if item.channel == "receptor_observations")
        sensor_registry_artifact = None
        if sensor_registry_path is not None:
            registry_path = Path(
                _string(sensor_registry_path, "sensor_registry_path")
            ).resolve()
            if not registry_path.is_file():
                raise FileNotFoundError(registry_path)
            sensor_registry_artifact = FieldEvidenceManifestArtifact(
                channel="sensor_registry",
                path=str(registry_path),
                sha256=_sha256_file(registry_path),
            )
        manifest = cls(
            manifest_id=manifest_id,
            event_id=event_id,
            audit_root=str(Path(audit.root).resolve()),
            channel_artifacts=tuple(artifacts),
            dataset_id=dataset_id,
            observed_dataset_path=receptor.path,
            observed_dataset_sha256=receptor.sha256,
            observed_row_count=observed_row_count,
            source_boundary_id=source_boundary_id,
            weather_id=weather_id,
            obstacle_geometry_id=obstacle_geometry_id,
            receptor_geometry_id=receptor_geometry_id,
            temporal_operator_id=temporal_operator_id,
            common_clock_id=common_clock_id,
            scope=scope,
            sensor_set_id=sensor_set_id,
            operator_id=operator_id,
            sensor_registry_artifact=sensor_registry_artifact,
            sensor_calibration_status=sensor_calibration_status,
        )
        manifest.verify_files()
        return manifest

    def verify_files(self) -> None:
        """Recheck root containment and every selected file digest."""
        root = Path(self.audit_root).resolve()
        if not root.is_dir():
            raise NotADirectoryError(root)
        for artifact in self.channel_artifacts:
            path = Path(artifact.path).resolve()
            try:
                path.relative_to(root)
            except ValueError as error:
                raise ValueError(
                    f"manifest artifact path escapes audit_root: {path}"
                ) from error
            if not path.is_file():
                raise FileNotFoundError(path)
            current = _sha256_file(path)
            if current != artifact.sha256:
                raise ValueError(
                    f"manifest artifact changed since its provenance was recorded: {path}"
                )
        if self.sensor_registry_artifact is not None:
            artifact = self.sensor_registry_artifact
            path = Path(artifact.path).resolve()
            try:
                path.relative_to(root)
            except ValueError as error:
                raise ValueError(
                    f"manifest sensor registry path escapes audit_root: {path}"
                ) from error
            if not path.is_file():
                raise FileNotFoundError(path)
            current = _sha256_file(path)
            if current != artifact.sha256:
                raise ValueError(
                    "manifest sensor registry changed since its provenance was recorded"
                )
        observed = Path(self.observed_dataset_path).resolve()
        if _sha256_file(observed) != self.observed_dataset_sha256:
            raise ValueError(
                "manifest observed dataset changed since its provenance was recorded"
            )

    def as_validation_evidence(self) -> FieldValidationEvidence:
        """Return the existing score-input evidence after digest verification."""
        self.verify_files()
        return FieldValidationEvidence(
            dataset_id=self.dataset_id,
            path=str(Path(self.observed_dataset_path).resolve()),
            sha256=self.observed_dataset_sha256,
            row_count=self.observed_row_count,
            source_boundary_id=self.source_boundary_id,
            weather_id=self.weather_id,
            obstacle_geometry_id=self.obstacle_geometry_id,
            receptor_geometry_id=self.receptor_geometry_id,
            temporal_operator_id=self.temporal_operator_id,
            common_clock_id=self.common_clock_id,
            scope=self.scope,
        )

    def readiness_record(self) -> dict[str, object]:
        """Return the fail-safe readiness of the evidence package itself.

        A complete five-channel package with explicit sensor and accountable
        operator identities is ``accepted`` as an *evidence package* only.
        It is still never promoted to a validated model or design basis.
        Legacy/placeholder identities remain usable for traceability but are
        reported as ``conditional`` until the missing metadata is supplied.
        File disappearance or digest drift is handled by :meth:`verify_files`
        and therefore remains a hard verification failure (withheld).
        """
        missing: list[str] = []
        if self.sensor_set_id.strip().lower() in _UNSPECIFIED_METADATA_VALUES:
            missing.append("sensor_set_id")
        if self.operator_id.strip().lower() in _UNSPECIFIED_METADATA_VALUES:
            missing.append("operator_id")
        if (
            self.sensor_set_id.strip().lower() not in _UNSPECIFIED_METADATA_VALUES
            and self.sensor_registry_artifact is None
        ):
            missing.append("sensor_registry_artifact")
        if (
            self.sensor_set_id.strip().lower() not in _UNSPECIFIED_METADATA_VALUES
            and self.sensor_calibration_status != "certified"
        ):
            missing.append("sensor_calibration_certificate")
        status = "conditional" if missing else "accepted"
        return {
            "status": status,
            "missing_requirements": missing,
            "promotion_allowed": False,
            "gate_codes": [
                "manifest_metadata_conditional" if missing
                else "manifest_metadata_complete",
                "manifest_not_promoted",
            ],
        }

    def as_record(self) -> dict[str, object]:
        return {
            "schema": FIELD_EVIDENCE_MANIFEST_SCHEMA,
            "manifest_id": self.manifest_id,
            "event_id": self.event_id,
            "audit_root": self.audit_root,
            "channel_artifacts": [item.as_record() for item in self.channel_artifacts],
            "dataset_id": self.dataset_id,
            "observed_dataset_path": self.observed_dataset_path,
            "observed_dataset_sha256": self.observed_dataset_sha256,
            "observed_row_count": self.observed_row_count,
            "source_boundary_id": self.source_boundary_id,
            "weather_id": self.weather_id,
            "obstacle_geometry_id": self.obstacle_geometry_id,
            "receptor_geometry_id": self.receptor_geometry_id,
            "temporal_operator_id": self.temporal_operator_id,
            "common_clock_id": self.common_clock_id,
            "scope": self.scope,
            "sensor_set_id": self.sensor_set_id,
            "operator_id": self.operator_id,
            "sensor_calibration_status": self.sensor_calibration_status,
            "sensor_registry_artifact": (
                None if self.sensor_registry_artifact is None
                else self.sensor_registry_artifact.as_record()
            ),
            "evidence_readiness": self.readiness_record(),
            "promotion_allowed": False,
            "gate_codes": ["manifest_files_verified", "manifest_not_promoted"],
        }


def field_evidence_manifest_record(manifest: FieldEvidenceManifest) -> dict[str, object]:
    if not isinstance(manifest, FieldEvidenceManifest):
        raise TypeError("manifest must be a FieldEvidenceManifest")
    return manifest.as_record()


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def field_evidence_manifest_from_mapping(
    value: object, *, base_directory: str | Path | None = None,
) -> FieldEvidenceManifest:
    data = _mapping(value, "field evidence manifest")
    _keys(
        data,
        "field evidence manifest",
        required={
            "schema", "manifest_id", "event_id", "audit_root", "channel_artifacts",
            "dataset_id", "observed_dataset_path", "observed_dataset_sha256",
            "observed_row_count", "source_boundary_id", "weather_id",
            "obstacle_geometry_id", "receptor_geometry_id", "temporal_operator_id",
            "common_clock_id", "scope", "promotion_allowed", "gate_codes",
        },
        optional={
            "sensor_set_id", "operator_id", "sensor_registry_artifact",
            "sensor_calibration_status", "evidence_readiness",
        },
    )
    if data["schema"] != FIELD_EVIDENCE_MANIFEST_SCHEMA:
        raise ValueError("schema must be " + repr(FIELD_EVIDENCE_MANIFEST_SCHEMA))
    if data["promotion_allowed"] is not False:
        raise ValueError("promotion_allowed must be false")
    gate_codes = data["gate_codes"]
    if (
        not isinstance(gate_codes, list)
        or len(gate_codes) != 2
        or not all(isinstance(code, str) for code in gate_codes)
        or set(gate_codes) != {"manifest_files_verified", "manifest_not_promoted"}
    ):
        raise ValueError("manifest gate_codes must contain the fixed fail-safe codes")
    root = Path(_string(data["audit_root"], "manifest.audit_root"))
    if not root.is_absolute() and base_directory is not None:
        root = Path(base_directory) / root
    artifacts_data = data["channel_artifacts"]
    if not isinstance(artifacts_data, list):
        raise ValueError("manifest.channel_artifacts must be a JSON array")
    artifacts: list[FieldEvidenceManifestArtifact] = []
    for index, item in enumerate(artifacts_data):
        record = _mapping(item, f"manifest.channel_artifacts[{index}]")
        _keys(record, f"manifest.channel_artifacts[{index}]", required={"channel", "path", "sha256"})
        path = Path(_string(record["path"], f"channel_artifacts[{index}].path"))
        if not path.is_absolute() and base_directory is not None:
            path = Path(base_directory) / path
        artifacts.append(
            FieldEvidenceManifestArtifact(
                channel=_string(record["channel"], f"channel_artifacts[{index}].channel"),  # type: ignore[arg-type]
                path=str(path.resolve()),
                sha256=_sha(record["sha256"], f"channel_artifacts[{index}].sha256"),
            )
        )
    sensor_registry_data = data.get("sensor_registry_artifact")
    sensor_registry_artifact = None
    if sensor_registry_data is not None:
        record = _mapping(
            sensor_registry_data, "manifest.sensor_registry_artifact",
        )
        _keys(
            record,
            "manifest.sensor_registry_artifact",
            required={"channel", "path", "sha256"},
        )
        path = Path(_string(record["path"], "sensor_registry_artifact.path"))
        if not path.is_absolute() and base_directory is not None:
            path = Path(base_directory) / path
        sensor_registry_artifact = FieldEvidenceManifestArtifact(
            channel=_string(
                record["channel"], "sensor_registry_artifact.channel",
            ),  # type: ignore[arg-type]
            path=str(path.resolve()),
            sha256=_sha(
                record["sha256"], "sensor_registry_artifact.sha256",
            ),
        )
    observed_path = Path(_string(data["observed_dataset_path"], "manifest.observed_dataset_path"))
    if not observed_path.is_absolute() and base_directory is not None:
        observed_path = Path(base_directory) / observed_path
    row_count = data["observed_row_count"]
    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count <= 0:
        raise ValueError("manifest.observed_row_count must be a positive integer")
    manifest = FieldEvidenceManifest(
        manifest_id=_string(data["manifest_id"], "manifest.manifest_id", reject_unspecified=True),
        event_id=_string(data["event_id"], "manifest.event_id", reject_unspecified=True),
        audit_root=str(root.resolve()),
        channel_artifacts=tuple(artifacts),
        dataset_id=_string(data["dataset_id"], "manifest.dataset_id", reject_unspecified=True),
        observed_dataset_path=str(observed_path.resolve()),
        observed_dataset_sha256=_sha(data["observed_dataset_sha256"], "manifest.observed_dataset_sha256"),
        observed_row_count=row_count,
        source_boundary_id=_string(data["source_boundary_id"], "manifest.source_boundary_id", reject_unspecified=True),
        weather_id=_string(data["weather_id"], "manifest.weather_id", reject_unspecified=True),
        obstacle_geometry_id=_string(data["obstacle_geometry_id"], "manifest.obstacle_geometry_id", reject_unspecified=True),
        receptor_geometry_id=_string(data["receptor_geometry_id"], "manifest.receptor_geometry_id", reject_unspecified=True),
        temporal_operator_id=_string(data["temporal_operator_id"], "manifest.temporal_operator_id", reject_unspecified=True),
        common_clock_id=_string(data["common_clock_id"], "manifest.common_clock_id", reject_unspecified=True),
        scope=_string(data["scope"], "manifest.scope"),  # type: ignore[arg-type]
        sensor_set_id=_string(
            data.get("sensor_set_id", "legacy-unspecified"),
            "manifest.sensor_set_id",
        ),
        operator_id=_string(
            data.get("operator_id", "legacy-unspecified"),
            "manifest.operator_id",
        ),
        sensor_registry_artifact=sensor_registry_artifact,
        sensor_calibration_status=_string(
            data.get("sensor_calibration_status", "legacy-unspecified"),
            "manifest.sensor_calibration_status",
        ),  # type: ignore[arg-type]
    )
    if "evidence_readiness" in data:
        supplied_readiness = _mapping(
            data["evidence_readiness"], "manifest.evidence_readiness",
        )
        if dict(supplied_readiness) != manifest.readiness_record():
            raise ValueError(
                "manifest.evidence_readiness must match the declared metadata"
            )
    return manifest


def write_field_evidence_manifest_json(
    manifest: FieldEvidenceManifest, path: str | Path,
) -> FieldEvidenceManifest:
    """Write a hash-pinned manifest without overwriting an existing file."""
    if not isinstance(manifest, FieldEvidenceManifest):
        raise TypeError("manifest must be a FieldEvidenceManifest")
    manifest.verify_files()
    output = Path(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing evidence manifest: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": FIELD_EVIDENCE_MANIFEST_EXECUTION_SCHEMA,
        "input": {"audit_root": manifest.audit_root},
        "field_evidence_manifest": manifest.as_record(),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    try:
        with output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(rendered)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing evidence manifest: {output}"
        ) from error
    return read_field_evidence_manifest_json(output)


def read_field_evidence_manifest_json(path: str | Path) -> FieldEvidenceManifest:
    """Read and verify a saved evidence manifest and all selected files."""
    source = Path(path)
    payload = json.loads(
        source.read_text(encoding="utf-8-sig"), object_pairs_hook=_json_object_pairs,
    )
    wrapper = _mapping(payload, "field evidence manifest execution")
    _keys(wrapper, "field evidence manifest execution", required={
        "schema", "input", "field_evidence_manifest",
    })
    if wrapper["schema"] != FIELD_EVIDENCE_MANIFEST_EXECUTION_SCHEMA:
        raise ValueError(
            "schema must be " + repr(FIELD_EVIDENCE_MANIFEST_EXECUTION_SCHEMA)
        )
    input_data = _mapping(wrapper["input"], "manifest input")
    _keys(input_data, "manifest input", required={"audit_root"})
    manifest = field_evidence_manifest_from_mapping(
        wrapper["field_evidence_manifest"], base_directory=source.parent,
    )
    input_root = Path(_string(input_data["audit_root"], "input.audit_root"))
    if not input_root.is_absolute():
        input_root = source.parent / input_root
    if input_root.resolve() != Path(manifest.audit_root).resolve():
        raise ValueError("input.audit_root must match manifest.audit_root")
    manifest.verify_files()
    return manifest


__all__ = [
    "FIELD_EVIDENCE_MANIFEST_SCHEMA",
    "FIELD_EVIDENCE_MANIFEST_EXECUTION_SCHEMA",
    "FIELD_EVIDENCE_MANIFEST_GATE_CODES",
    "SENSOR_CALIBRATION_STATUSES",
    "FieldEvidenceManifestArtifact",
    "FieldEvidenceManifest",
    "field_evidence_manifest_record",
    "field_evidence_manifest_from_mapping",
    "write_field_evidence_manifest_json",
    "read_field_evidence_manifest_json",
]
