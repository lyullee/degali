# Strict named field batch input

`degali field-batch` runs several strict field cases into one new output
directory. It is a deployment wrapper around
`export_field_screening_batch()`; it does not infer source, weather, geometry,
or units from filenames.

The input schema is `degali.field-batch-input.v1`:

```json
{
  "schema": "degali.field-batch-input.v1",
  "options": {
    "include_refinement": true,
    "refinement_factors": [1, 2],
    "relative_tolerance": 0.05,
    "max_uncertainty_cases": 64,
    "table_nodes": 161
  },
  "cases": [
    {
      "label": "post-flash-event-A",
      "case_path": "field-case.json",
      "atmospheric_source_case_path": "source-case.json"
    }
  ]
}
```

`label` is also the emitted report filename stem. It must match
`[A-Za-z0-9][A-Za-z0-9_.-]{0,127}`; path separators, whitespace-only names and
other path-like values are rejected before the referenced case/source files are
loaded. The parsed Python `FieldBatchInput` keeps its cases as an immutable
tuple, including when constructed directly rather than through JSON.

The JSON reader rejects exact duplicate keys and keys that collide after
whitespace/case normalization before loading any referenced case or source
artifact. A last-value-wins JSON edit therefore cannot silently change a case
path, source schedule, or batch option.

`case_path` must point to a strict `degali.field-screening-input.v1` case. The
batch parser records its resolved path and SHA-256. An optional
`atmospheric_source_case_path` points to a strict
`degali.field-atmospheric-schedule-input.v1` file; its path, SHA-256 and full
source evidence are retained in `batch-execution.json`. The schedule is routed
through the complete operational uncertainty envelope, so its source-rate
corners are combined with the declared weather/sensor corners.

Historian and phase-routing cases are rejected by this wrapper. They must use
their dedicated operational envelopes so a batch cannot silently mix different
source physics. Source-schedule cases always require refinement, an in-plane
receptor and resolved uncertainty. If conditional screening is explicitly
enabled, the `conditional_review_authorizations` mapping must contain one
authorization for every named case; that case-level authorization is retained
in `manifest.json` while every deterministic envelope corner still has to pass
the non-negotiable gates independently.

Run:

```console
degali field-batch batch.json --output-directory field-batch-output
```

The output directory must be new or empty. Every requested case is computed and
every case/manifest JSON payload is rendered before that directory is created;
a late report serialization or finite-value failure therefore leaves no
partial batch directory. It contains one case report,
`manifest.json`, and `batch-execution.json`. Both top-level execution records
include a `degali.field-batch-summary.v1` `summary` with case count,
decision-status counts, typed gate-code counts, allowed and withheld counts, and an
`all_screening_allowed` fail-safe flag derived from the actual operational
decisions. Add `--require-screening` to exit with status 2
if any named case remains withheld; the reports are still kept for audit.
The Python `FieldBatchExport` wrapper re-reads the emitted manifest and reports
before accepting the export, checking schema, case identity, decision summary,
review/refinement metadata and report schema. Existing files with a valid path
but mismatched content are rejected.
`field-verify` independently recomputes that summary from the manifest case
decisions, including `gate_code_counts`, so changing both the manifest and
execution summary cannot hide a fabricated batch disposition. Legacy v3
manifests without typed gate codes remain integrity-checkable.
The direct Python `write_field_screening_report()` and
`write_field_refinement_report()` helpers also use exclusive creation, so an
existing report path is rejected rather than overwritten (including a
concurrent-writer race).
Operational-envelope case reports and the emitted `manifest.json` use the same
exclusive creation boundary after JSON serialization; the CLI's
`batch-execution.json` is also created exclusively. Thus a competing writer
cannot replace a typed batch artifact after the directory preflight.
The manifest also stores a SHA-256 for each case report; the CLI's
`batch-execution.json` stores those report hashes plus the manifest hash for
post-export integrity checks.
Its compact per-case execution summary also repeats the manifest's typed
`gate_codes`; the verifier rejects a status/allowance/code mismatch instead
of leaving automation to reopen the full manifest to explain a withheld case.
Operational envelope case reports include a
`deterministic_sensor_envelope` with per-detector extrema over completed
corners; withheld detector traces remain counted with their reasons rather than
being replaced by zero.
Every emitted `operational_decision` also carries machine-readable `gate_codes`
in addition to the explanatory `reasons` and `required_actions`. Codes identify
the stable gate category (for example `refinement_missing`,
`sensor_trace_missing`, `obstacle_unrepresented`, or
`uncertainty_unresolved`) without requiring automation to parse prose.
For uncertainty/history/sensor/phase aggregate reports, the top-level
`gate_codes` is the stable union of all corner decisions; a conditional
aggregate holdback includes `conditional_review_required`, while a legacy
withheld corner with no typed code is represented by `corner_withheld`.

`field-verify` also checks field-screening and batch execution records. For a
nominal or direct atmospheric-source screening record it reopens the strict
case, replays the recorded numerical options and recomputes the report. Complex
phase/history/sensor-array records retain `report_recomputed=false` when the
execution artifact does not contain all replay options. For a batch record it
reopens the strict batch input and checks the manifest plus every report SHA,
then parses each report and requires the typed screening or operational-envelope
report schema (including agreement with the manifest execution kind); a hashed
arbitrary JSON object therefore cannot stand in for a case result. Batch-level
promotion remains governed by the manifest decision summary. It also compares
the execution summary and each compact status/`screening_allowed` field against
the manifest, so a caller cannot rewrite the decision-facing summary while
leaving the report files untouched.
Automation that requires a fresh numerical replay can add
`degali field-verify batch-execution.json --require-recomputed`; batch records
currently return exit 2 under that flag because their manifest/report integrity
is verified without rerunning every named field case.
