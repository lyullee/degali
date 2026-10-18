# Field evidence manifest

`field-audit` only discovers whether the five required channels are visible.
It must not infer that files describe the same event. When an accountable
reviewer has reconciled the files, `FieldEvidenceManifest.from_audit()` records
that decision explicitly:

- one selected file and SHA-256 for `source_boundary`, `weather`,
  `obstacle_geometry`, `receptor_observations`, and `common_clock`;
- an event ID, common-clock ID, source/weather/obstacle/receptor IDs and the
  temporal operator ID;
- a `sensor_set_id` pointing to the sensor registry/calibration export and an
  accountable `operator_id` for the person or system that reconciled the files;
- the observed dataset ID, positive row count, path and digest.

The manifest is still a qualification input, not validation. Its serialized
record always contains `promotion_allowed: false` and the fixed gate codes
`manifest_files_verified` and `manifest_not_promoted`. It can be converted to
the existing `FieldValidationEvidence` only after every selected file is
inside the audited root and its current bytes match the recorded digest.

The package also emits `evidence_readiness`. It is `accepted` only when
`sensor_set_id`, `operator_id`, a pinned registry artifact and an explicit
`sensor_calibration_status="certified"` are present. A nominal instrument
specification (`specification_only`), an unavailable certificate (`missing`) or
a legacy manifest remains traceable but is marked `conditional`; it must not be
used as evidence of a qualified model. Missing channels are stopped earlier by
`field-audit` as `partial` or `withheld`, and digest/path drift is a hard
verification failure.

For quantitative transient scoring the package now emits a separate
`quantitative_evidence_readiness` record.  The selected common-clock file must
also declare `reference_time_utc`, `synchronization_method`, `offset_s`,
`uncertainty_s`, and `covered_channels` containing at least
`source_boundary`, `weather`, and `receptor_observations`.  A legacy clock file
that contains only `common_clock_id` remains traceable but is automatically
`conditional` for scoring.  This prevents a shared label from being mistaken
for a verified time alignment.

The metadata may be a one-row CSV or a JSON object (optionally nested under
`clock_alignment`), for example:

```text
common_clock_id,reference_time_utc,synchronization_method,offset_s,uncertainty_s,covered_channels
clock-a,2026-09-21T00:00:00Z,PTP,0.012,0.004,source_boundary;weather;receptor_observations
```

## Python boundary

```python
from degali.addons import (
    FieldEvidenceManifest,
    audit_field_evidence,
    write_field_evidence_manifest_json,
)

audit = audit_field_evidence("/data/lh2-event")
manifest = FieldEvidenceManifest.from_audit(
    audit,
    manifest_id="event-2026-09-21-review",
    event_id="lh2-trial-2026-09-21",
    selected_paths={
        "source_boundary": "/data/lh2-event/source.csv",
        "weather": "/data/lh2-event/weather.csv",
        "obstacle_geometry": "/data/lh2-event/obstacle.json",
        "receptor_observations": "/data/lh2-event/receptors.csv",
        "common_clock": "/data/lh2-event/clock.csv",
    },
    dataset_id="receptors-2026-09-21",
    observed_row_count=320,
    source_boundary_id="source-boundary-a",
    weather_id="weather-a",
    obstacle_geometry_id="obstacle-a",
    receptor_geometry_id="receptors-a",
    temporal_operator_id="mean-275s",
    common_clock_id="clock-a",
    sensor_set_id="sensor-registry-2026-09-21-r1",
    operator_id="operator-jkim",
    sensor_registry_path="/data/lh2-event/sensor_registry.csv",
    sensor_calibration_status="certified",
)
write_field_evidence_manifest_json(manifest, "field-evidence-manifest.json")
```

`from_audit()` requires the audit status to be exactly `candidate_complete`; it
refuses a partial/withheld or incomplete scan, an unselected channel, a path
that was not detected in the audit, or a candidate without a SHA-256. The reader
rechecks all five files, so changing a source, weather, geometry, receptor or
clock file requires a fresh audit and manifest. No manifest is generated for
the current connected folders because their required channels remain missing.

For the public USN/FFI wide-table exports, the audit accepts the unlabeled
leading row-index column and records the raw source/sensor histories with
SHA-256, row counts and the separate sensor registry. The reproducible
`tools/audit_ffi_field_evidence_manifest.py` intake writes
`outputs/ffi-o1-field-evidence-intake-2026-10-14-v2.json` and the corresponding
O2 record. These are deliberately `withheld` intake records—not manifests:
weather, machine-registered beam geometry, a source/sensor synchronization
transform and event calibration certificates are still absent.

The same explicit boundary is available from the CLI. The command consumes a
verified `field-audit` execution artifact and requires one `--selected-path`
for each channel; it never infers the event join or overwrites an existing
manifest:

```console
degali field-evidence-manifest-create field-audit-execution.json \
  --selected-path source_boundary=/data/lh2-event/source.csv \
  --selected-path weather=/data/lh2-event/weather.csv \
  --selected-path obstacle_geometry=/data/lh2-event/obstacle.json \
  --selected-path receptor_observations=/data/lh2-event/receptors.csv \
  --selected-path common_clock=/data/lh2-event/clock.csv \
  --manifest-id event-2026-09-21-review \
  --event-id lh2-trial-2026-09-21 \
  --dataset-id receptors-2026-09-21 \
  --observed-row-count 320 \
  --source-boundary-id source-boundary-a \
  --weather-id weather-a \
  --obstacle-geometry-id obstacle-a \
  --receptor-geometry-id receptors-a \
  --temporal-operator-id mean-275s \
  --common-clock-id clock-a \
  --sensor-set-id sensor-registry-2026-09-21-r1 \
  --operator-id operator-jkim \
  --sensor-registry-path /data/lh2-event/sensor_registry.csv \
  --sensor-calibration-status certified \
  --output field-evidence-manifest.json
```

For automation, verify a saved artifact with:

```console
degali field-evidence-manifest-verify field-evidence-manifest.json
```

The command reports `selected_files_verified=true` and keeps
`promotion_allowed=false` only after all five digests match.
If a selected file is missing, changed, or unreadable, it exits with code 2
and emits `evidence_readiness.status="withheld"` with the failure reason; it
never reports a stale package as verified.

The resolved O1 and fast→resolved→fast hybrid runners can attach the same
execution manifest with `--field-evidence-manifest`. They copy its event,
sensor-set, operator and quantitative-readiness record into the solver output.
For a production score, add `--require-quantitative-manifest`; a conditional
clock or calibration package then fails before the expensive solve. Without
that flag, the run remains an explicitly conditional physics witness.

To bind the manifest into a strict matched-sensor case, add this optional
object beside `validation_evidence` in `degali.field-validation-input.v1`:

```json
"validation_manifest": {
  "path": "../field-evidence-manifest.json",
  "sha256": "<sha256 of the manifest execution JSON>"
}
```

The case reader verifies the manifest and requires its observed path, digest,
row count, scope, and source/weather/geometry/operator IDs to equal the case's
`validation_evidence` record.
