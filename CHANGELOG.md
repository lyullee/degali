# Changelog

## Unreleased

## 0.3.0 - 2026-10-09

- Reserved and embedded the exact Zenodo release DOI
  `10.5281/zenodo.23256538` across package, citation and submission metadata.
- Completed the previously reserved v0.2.0 Zenodo archive at
  `10.5281/zenodo.23105451` before preparing this release, preserving the
  one-tag/one-DOI version chain.

- Rechecked the official KITopen/RADAR PRESLHY E3.5 record and selectively
  retrieved Trial 10. The 575,705-byte workbook matches the local controlled
  copy byte-for-byte; its five sheets and sampling intervals are now captured
  in a hash-pinned retrieval audit. The evidence remains conditional because
  the workbook has no separate channel calibration certificate, absolute UTC
  clock, surveyed obstacle sheet or event-level geometry manifest.
- Replayed the conditional common-clock extraction for accepted PRESLHY Trials
  20 and 21, added a three-event manifest index and gate audit, and expanded
  the IJHE hand-off bundle to 32 hash-pinned files. All three events remain
  conditional and non-promotable because event-level calibration and surveyed
  geometry are still absent.
- Retrieved the public RADAR Trial 20 and Trial 21 workbooks through their
  individual-file endpoints and verified byte-for-byte matches to the local
  replay inputs; the combined retrieval index now records Trial 10 as well.
  Added it to the IJHE hand-off;
  the refreshed bundle now contains 33 hash-pinned files. The calibration,
  clock and geometry promotion gates remain unchanged.
- Made the E3.5 selection rule explicit in the IJHE manuscript and supplement:
  all 24 public workbooks are acknowledged, nine deterministic
  horizontal/momentum-dominated trials form the peak score, and the seven raw
  clock-aligned workbooks are identified separately. The rebuilt Word/PDF
  artifacts were rendered and visually checked page-by-page.
- Added the external-evidence search, selective-retrieval audit and
  multi-event manifest index to the author-neutral IJHE hand-off bundle; the
  latest bundle contains 33 files, zero missing entries and remains
  final-upload-gated.
- Clarified the IJHE supplementary sensor boundary: a nominal 4 vol% Xensor
  indication maps to an approximately 1–7 vol% specification-only interval,
  whereas the nominal Dräger interval is 3.92–4.08 vol%. The text now states
  explicitly that this can make a binary 4 vol% Xensor decision indeterminate
  and is not a correction to the 11/18 model-operator result.
- Corrected the English and Korean user guides during release preparation so
  they no longer linked to an unpublished DOI; v0.2.0 is now archived at
  `10.5281/zenodo.23105451`, and this release is pinned to its own DOI.
- Corrected 69 relative links in the archived development log and verified all
  342 repository Markdown files with zero broken local links.
- Added explicit IJHE graphical-abstract minimum-dimension auditing and a
  regression test for the publication snapshot.
- Added an IJHE portal preflight matrix separating completed upload items from
  author, permission and live-portal confirmation gates.
- Clarified the manuscript's 11/18 decision-change count as binary 4 vol%
  classifications at a fixed 1.5 m receptor over the declared line/wind grid.
- Made the DOCX audit fail on local repository path markers in paragraphs or
  tables, not only in the builder regression test.
- Added a submission audit failure gate for accidental TODO/FIXME/TBD markers
  in journal-facing source files.
- Added `tools/audit_ijhe_submission.py` and
  `docs/ijhe-submission-readiness.md` for a machine-readable IJHE manuscript
  readiness check. The audit verifies claim-boundary language, numerical
  convergence evidence and public-snapshot safety while retaining an explicit
  warning until matched time-resolved field evidence is available.
- Local `outputs/` and `artifacts/` products are now ignored by the public
  snapshot. They can contain third-party workbook/PDF-derived material and
  remain reproducible local derivatives rather than distributable source.
- Updated `docs/publication-guide.md` to match the `0.3.0` package, release
  metadata and version-DOI workflow instead of the obsolete `0.1.0`
  instructions.
- Added separate IJHE upload drafts for highlights and a cover letter. The
  readiness audit checks highlight length and keeps author/declaration
  placeholders as an explicit manual warning.
- Reworked the manuscript abstract and data/code statements to foreground the
  evidence-bounded contribution, retain the FFI Test 6 limitation and identify
  the versioned public software snapshot.
- Added `IJHE_SUPPLEMENTARY_INFORMATION.md` with frozen configuration, lane
  separation, headline metrics, Test 6 limitation, 3-D numerical verification
  and evidence-gate scope exclusions.
- Added `tools/audit_open_channel_h2_dataset.py` and
  `docs/open-channel-h2-boundary.md` to hash and summarize the 22-file public
  open-channel H₂ source/sensor timing dataset without promoting it to
  atmospheric LH₂ validation. Raw negative flow samples remain visible rather
  than being clipped.
- Added an IJHE figure/table generation register and a standalone evidence-gated
  graphical abstract. The submission audit now parses the SVG and verifies the
  register's reproducibility fields.
- Added a reproducible 3-D transient concentration-error benchmark. It compares
  coarse/medium/fine grids against a declared high-resolution numerical
  reference, reports receptor MAE/RMSE/peak/dose/arrival errors, and separates
  numerical discretization error from experimental validation error.
- Added an explicit 3-D source/wind-coupled transient dense-gas transport
  path. `transient_dense_gas_3d` integrates atmospheric source schedules and
  meteorological wind histories on a common clock, applies density-driven
  vertical buoyancy through a prognostic velocity state, returns receptor time
  traces, and audits cumulative mass residuals. The scope is an explicit
  reduced-order finite-volume solver, not a variable-density LES/RANS
  replacement.
- Added an opt-in source-state ledger and conservative post-impact wall-jet
  operators for SLABx-LH2 parity work. The ledger can hand off gas-only H2
  schedules to the existing semi-FV obstacle transport, while radial,
  bifurcated, lateral-mixing and vertical-profile operators preserve scalar
  inventory and remain withheld from horizontal FFI Test4 by default.
- Field-evidence intake now carries explicit `sensor_set_id` and accountable
  `operator_id` metadata. A supplied sensor registry/calibration file can be
  pinned with its own SHA-256 as `sensor_registry_artifact`; missing registry
  or operator metadata is reported as `conditional`, while missing channels or
  digest drift remain fail-safe `partial`/`withheld`. Added the operator CLI
  flags, intake template, strict round-trip checks and focused regression
  coverage; the existing five audited channels and `promotion_allowed=false`
  boundary are unchanged.
- Re-ran the connected external-evidence audit against the current folders as
  v14 artifacts. SLABx now scans 16,139 structured files (four additional
  files, with the same 15 source-boundary and 33 receptor candidates); no
  weather, obstacle-geometry or common-clock channel appeared, so the audit
  remains `partial` and `promotion_allowed=false`. Support-centre and related
  material remain `withheld`. All three v14 artifacts pass `field-audit-verify`.
- Integrity-only `field-verify` now cross-checks each serialized operational
  allowance against its attached physical applicability, including nested
  envelope cases. A tampered `blocked` report can no longer retain an allowed
  decision without detection; the focused field target is now **490 passed**.
  The complete repository regression is **1,659 passed, 145 skipped, 0 failures**.
- Field sensor-branch superposition now requires every branch to expose the
  identical receptor time grid, not merely matching start/end times. A branch
  with a different internal grid is rejected before scalar traces are summed,
  preventing implicit interpolation from mixing temporal operators or hiding a
  source shutoff/peak.
- Historian CSV paired calibration columns now require each row to satisfy
  `lower <= nominal <= upper` and `lower <= upper`; reversed or nominal-excluding
  bounds are rejected before source flashing or uncertainty propagation rather
  than being sorted or clipped.
- Semi-FV obstacle inputs now reject rectangles that touch or extend past the
  local inlet/outlet instead of allowing the mesh mask to clip declared geometry;
  typed diagnostics also reject non-finite/negative mass and Courant values,
  invalid diverted fractions and duplicate source-ledger labels.
- Fixed-sensor model-comparison thresholds are now constrained to the physical
  mole-fraction range `[0, 1]` at case, function and typed-result boundaries;
  an out-of-range threshold can no longer classify every sensor as below-limit.
- Operational envelope case constructors now cross-check an allowed decision
  against the attached physical applicability status: `screening_allowed`
  requires `accepted`, while `conditional_allowed` requires `conditional`.
  A fabricated allowed decision can no longer be paired with a blocked or
  mismatched physical result before aggregate serialization.
- Field-evidence audit triage now recognizes explicit `source_bound`/
  `source_boundary`, reference-wind speed aliases (`wind_ref_m_s`, `u_ref_m_s`,
  `wind_10m_m_s`) and measured mean-concentration aliases (`mean_c_pct`,
  `observed_mean_c_pct`, `mean_vol_pct`). These remain candidate inventory
  signals only; old collection plans remain readable and no event join or
  validation promotion is inferred.
- Field-evidence audit artifacts now include a deterministic
  `collection_requirements` plan for all five channels: missing-vs-detected
  status, candidate paths, accepted field-alias groups, manual qualification
  boundary and the next collection action. The plan is derived from the
  channel inventory, is checked on strict read-back, and remains advisory;
  `promotion_allowed` stays false and older v1 artifacts without the optional
  field remain readable.
- Field-evidence common-clock detection now accepts explicit
  `synchronization_id`, `sync_id`, `time_sync_id`, and `clock_sync_id` aliases,
  plus synchronization metadata when paired with timestamps. Blank values
  remain non-evidence; this improves triage only and does not relax promotion
  gates. New collection plans advertise the expanded explicit aliases, while
  already-published v1 plans using the historical alias group remain readable.
- Added the explicit `PressureDrivenMeasuredHistory` source path.  Absolute
  historian pressure/temperature states are recomputed through the declared
  orifice boundary and flashed interval by interval; pressure trends are never
  treated as leak rates.  Coherent P/T/liquid-fraction source corners are
  available through `direct_vapour_schedule_envelope_from_pressure_driven_history()`
  and `run_field_pressure_driven_history_envelope()`, with coherent opening-area/
  discharge-coefficient corners, approved quality gates and fail-safe rejection of unresolved
  ambient-pressure uncertainty.  Added field-transport regression coverage and
  public exports; the joint pressure-driven history path now propagates source
  location/direction before wind-plane projection, plus ambient,
  weather, stability and detector corners and has an operational refinement/
  aggregate fail-safe wrapper.  Strict `pressure_driven_history` case JSON and
  CSV import/`field-verify` replay now preserve the pressure-only boundary
  without fabricating a flow channel.  The current field/semi-FV target is
  **489 passed**.  The latest complete repository regression is **1,658 passed,
  145 skipped, 0 failures**.
- The approved measured-history joint envelope now carries declared source
  location/direction corners into the field selection and fixes them before
  wind-plane projection, receptor applicability and obstacle transport. The
  operational measured-history wrapper and deterministic sensor summary retain
  every such corner; the field/semi-FV target remains **489 passed**.
- The measured-history source-only field envelope now applies the same declared
  source location/direction corner propagation before wind-plane projection.
  Historian P/T/flow/liquid-fraction selections remain separate from this
  geometry cross-product, with a max-case guard and serialized geometry
  selections; weather, sensor and obstacle uncertainty remain outside this
  source-only helper.
- Added `docs/field-goal-audit-2026-10-07.md`, a requirement-by-requirement
  evidence table that keeps implemented field capabilities separate from the
  still-missing external source/weather/obstacle/common-clock validation set.
- Field-evidence audit near-misses now recognize common `source_rate_kg_s`,
  `mass_rate_reported_kg_s`, `rate_kg_s`, `wind_m_s`, and mean-wind aliases.
  They remain explicit diagnostics when source/weather identity is missing and
  are never promoted to validation channels; refreshed v13 audits keep SLABx
  partial and support-centre evidence withheld. The v13 scan records
  16,135 structured files, 15 source-boundary candidate paths and 33
  receptor-candidate paths, with the same weather/obstacle/common-clock
  promotion holdback. Audit artifacts now include
  deterministic `diagnostic_counts` by candidate-file near-miss category and
  reject stale or tampered counts on read-back. Timestamped PT/TT/LT
  instrument histories are also retained as non-source diagnostics rather than
  inferred leak boundaries.
- Added the hash-pinned `FieldEvidenceManifest` boundary. It requires an
  explicit candidate path for every evidence channel plus event/clock/source/
  weather/geometry/operator IDs, rechecks all selected files, can feed the
  existing `FieldValidationEvidence`, and always serializes
  `promotion_allowed=false`. `degali field-evidence-manifest-verify` exposes
  the same fail-safe file recheck for automation. The new
  `field-evidence-manifest-create` command turns a verified audit plus five
  explicit `channel=path` selections into that manifest without inferring an
  event join or overwriting an existing artifact. Strict validation cases can
  now optionally pin that manifest and reject any evidence-ID/path drift.
- The manifest-create CLI regression brings the current field/semi-FV target
  to **449 passed**. The last complete repository regression immediately
  before this CLI-only addition was **1,619 passed, 145 skipped, 0 failed**.
- Added the executable `tools/audit_ffi_test6_transient_receptors.py` replay
  promised by the Test 6 observation-operator note. It strictly reads the
  caller-supplied wind CSV, fingerprints it, reuses one steady free-field
  LH2 trajectory per declared wind speed, rotates the fixed reference sensors,
  applies only the declared first-order response, and emits compact
  lower-bound statistics with promotion and validation qualification disabled.
  `--verify` rechecks the history/reference hashes and compact report without
  rerunning CoolProp. The new replay and existing transient-receptor tests pass
  (**7 passed**).
- A fresh complete repository regression after the manifest CLI and Test 6
  replay additions passes **1,622 tests with 145 skips and 0 failures**.
- Test 6 replay verification now binds the observation-operator discriminator,
  no-storage/no-time-varying-source flags, exact response-time metadata,
  window sample count, sensor coverage and recomputed peak summaries; tampered
  compact rows fail closed without rerunning the plume.
- Test 6 replay inputs now retain the history time range, interval range,
  uniform-sampling diagnostic and wind-speed range. The values are rechecked
  from the hash-pinned CSV but remain diagnostics, not an automatic quality or
  validation promotion.
- The manifest-create CLI now has an explicit partial-audit regression: an
  incomplete five-channel audit returns exit 1 and leaves no output artifact.
  The current field/semi-FV target is **450 passed**.
- `FieldEvidenceManifest.from_audit()` now requires the audit discriminator
  `candidate_complete` explicitly, before inspecting any selected path. A
  direct partial-audit regression keeps incomplete readiness from being
  mistaken for a channel-selection error; the focused field/semi-FV target is
  now **451 passed**.
- `field-verify` now binds non-replayed measured-history reports back to the
  currently imported historian CSV provenance (SHA, event identity and every
  nested field report). Changing the CSV without changing the case JSON now
  fails closed instead of silently accepting a stale envelope; the focused
  field/semi-FV target is **452 passed**.
- Measured-history `field-screen` executions are now deterministically
  replayed by `field-verify` for both nominal/refined and joint-history
  envelope modes. The verifier compares the complete report and operational
  decision after reimporting the fingerprinted CSV; a changed numeric report
  now fails on fresh recomputation as well as on provenance drift. The focused
  field/semi-FV target remains **452 passed** (the existing measured-history
  envelope regression now also covers replay and tamper rejection).
- Phase-routing transport, standalone sensor-array, and joint source-sensor
  envelope executions now use the same replay path when their serialized case
  retains the complete phase/calibration/source options. The source-sensor
  path also reopens and rechecks its optional atmospheric source schedule;
  the ordinary bounded field uncertainty (`nominal_field`) path now replays
  its stored corner options as well. Genuinely option-incomplete envelopes
  remain integrity-only rather than receiving fabricated replay settings.
- Field-evidence audit key normalization now preserves Unicode letters and
  digits in workbook/JSON/CSV inventories. Multi-sheet XLSX sampling skips an
  ambiguous auxiliary sheet while retaining valid sheets, but keeps an
  all-ambiguous workbook unreadable. This keeps multilingual workbooks
  readable for diagnostics without expanding the explicit ASCII channel
  aliases or promoting a channel by inference. The refreshed v4 audits scanned
  15,714 SLABx files (33 receptor candidates, still `partial`); the support
  centre retains 20 files with no required channel and the related-material
  folder remains `withheld`.
- A current free-field preservation subset covering the axisymmetric jet,
  cryogenic-air/blowdown, LH2 scope/yaw/rainout, reference-parity,
  model-comparison and screening-gate paths completed with **247 passed, 139
  skipped, 0 failures**. Skips are optional CoolProp/reference-environment
  cases; this is a focused preservation check, not a replacement for the full
  repository baseline.
- Read-only field-evidence audit artifacts now bind `channel_paths`, coherent
  candidates, status and scan completeness back to the serialized candidate
  records. A tampered channel index can no longer pass `field-audit-verify`
  merely because the referenced files and their hashes are still present.
- Field-evidence readiness scanning now recognizes explicit observed-volume
  percentage aliases such as `observed_time_mean_vol_pct` as quantitative
  fixed-receptor observations, while still requiring the independent source,
  weather, obstacle and common-clock channels before any promotion.
- The latest repository regression excluding tests marked `slow` completed with
  **1,556 passed, 134 skipped, 63 deselected**; no failures were observed.
  The 63 deselected cases are the explicitly marked long-running integration
  set, so this result is a preservation check rather than a replacement for a
  complete-repository baseline.
- Replayed the external SLABx Test 4 fixed-sensor comparison into a fresh v3
  execution artifact. Current `field-verify` recomputation passes; the result
  remains `conditional`, `promotion_allowed=false`, and model-selection impact
  `withheld` because native source and temporal operators differ.
- Integrity-only `field-verify` now binds pressure-derived source reports back
  to the strict case boundary, including nested sensor/envelope reports: typed
  ambient pressure, uncertainty provenance, source ID, and derivation metadata
  are checked even when the numerical field run cannot be replayed from the
  saved execution options.
- After the pressure-driven LH2 source-rate adapter extension, the field/semi-
  FV target is **448 passed** (including audit-artifact binding, Unicode-header
  inventory and observed-
  volume-percent readiness regressions).
  The last complete repository baseline before
  this opt-in adapter was **1,601 passed with 145 skips and 0 failures**; the
  adapter itself is covered by the focused `test_field_lh2.py` module (21
  passed) plus a strict-case parser regression. The adapter now rejects
  non-hydrogen or gauge-pressure boundaries before throat evaluation, and the
  typed boundary requires an explicit ambient source and matching rate source
  ID.
- `field-verify` now applies the semi-FV conservation gate to serialized
  transport diagnostics even when a history/phase/sensor envelope cannot be
  replayed numerically. It checks finite residuals, per-source schedule
  residual uniqueness and maximum consistency, and requires a matching typed
  blocked applicability reason and operational `gate_codes` for any residual
  beyond tolerance; nested batch reports receive the same fail-safe check.
  The same integrity-only path now validates compact transport concentration
  and sensor extrema/time/density records, including withheld-detector warning
  and result consistency. Serialized operational decisions are also checked for
  status/allowance/reason/action invariants and valid unique gate codes.
  Source-level injected-mass labels and totals are cross-checked against the
  ledger residual and schedule-residual source set.
- The same integrity-only verifier now binds the source ledger back to the
  serialized transport-input provenance: declared direct-vapour and
  distributed-source schedule masses must match `primary` and
  `distributed:<label>` injected entries, with duplicate/missing labels
  rejected before an operational decision can be accepted.
- Serialized transport diagnostics now also recompute the closed-form final
  inventory identity, `mass_injected = mass_domain + mass_outflow`, against
  the same numerical tolerance. A forged maximum residual can no longer hide
  inconsistent or negative domain/outflow totals in an integrity-only
  artifact.
- The semi-FV solver now records `final_mass_residual_kg` in its typed
  diagnostics and applies the same final-inventory gate during direct field
  screening, not only during post-execution artifact verification.
- The fingerprinted atmospheric-source boundary now accepts an explicitly
  declared `declared_atmospheric_vapour` source kind in addition to
  `post_flash_atmospheric_vapour`. It remains tied to the scenario's explicit
  source location and nominal flash contract; pool/droplet schedules still
  require their location-bearing distributed or phase-routing handoff.
- Direct Python release, weather, sensor and surface boundaries now reject a
  non-empty unit label that disagrees with the quantity's expected SI unit.
  Blank labels remain accepted for legacy direct `BoundedValue` callers, while
  strict JSON continues to require the exact unit label before construction.
- `SourceRateSchedule` keeps its legacy `piecewise_constant` default and now
  supports an explicit `linear` time-node operator. Exact clipped-interval
  mass integration feeds the same semi-FV injection and source-ledger gates;
  strict JSON/CSV provenance records and checks the operator across all rate
  corners.
- `pressure_driven_lh2_mass_flow()` provides an opt-in pressure-driven LH2
  source-rate adapter. It crosses upstream pressure/temperature, opening area,
  discharge coefficient, and optional ambient-pressure endpoints through the
  existing homogeneous-equilibrium throat closure, returning an explicit
  `kg/s` bound with source provenance. The companion release adapter refuses
  to overwrite a non-zero measured rate, records the derivation metadata, and
  fails closed when any thermodynamic corner cannot be evaluated. Strict field
  JSON now exposes the same behavior only through an explicit
  `pressure_driven_mass_flow` object; omitted rate and derivation remain
  invalid rather than triggering implicit inference.
- Pressure-derived source rates now retain a typed derivation boundary in the
  field contract. Deterministic uncertainty corners vary the driving inputs
  and recompute the rate at the selected ambient-pressure corner instead of
  independently crossing an already-derived mass-flow interval with its own
  inputs; the structured derivation is retained in `ReleaseSource.as_dict()`.
- The v1 field execution verifier treats a missing schedule operator in an
  older serialized record as the unambiguous legacy `piecewise_constant`
  default, while preserving an explicitly supplied operator for exact
  provenance comparison.
- Direct `FieldSemiFVRequest` ambient uncertainty bounds now enforce the same
  physical units as the strict case parser (`K`, `Pa`, and `kg/m3`), preventing
  an incorrectly labelled deterministic corner from entering flash or sensor
  propagation.
- `SourceRateSchedule` now exposes one shared zero-endpoint contract. Finite
  semi-FV primary and distributed-source transport rejects a non-zero retained
  endpoint before injection, while low-level schedule inspection remains
  available; field/CSV schedule boundaries use the same helper and error basis.
- Semi-FV diagnostics now independently compare each source's declared
  schedule mass with the mass actually injected by the solver. A mismatch is
  retained per source, exposed through the report, and blocks field screening
  with the typed `source_schedule_mass_mismatch` gate code.
- Measured-history and pool-vapour schedule wrappers now require a typed,
  zero-terminated `SourceRateSchedule`, validate warning/alignment metadata,
  and cross-check measured-flash or pool-ledger mass before promotion to a
  field source boundary.
- The joint measured-history envelope now propagates ambient temperature and
  pressure corners through each history re-flash and ambient air-density
  corners through detector conversion, retaining and clearing those selections
  per resolved field case instead of leaving them as an unresolved operational
  holdback.
- Direct Python field screening/refinement report writers now use exclusive
  file creation: an existing evidence path is rejected rather than silently
  replaced, including a race for the same destination. Batch operational
  reports, manifests, and CLI `batch-execution.json` use the same exclusive
  creation boundary after payload serialization; the new-or-empty output
  directory preflight remains in place. Report/batch/execution regression
  coverage is **36 passed**, with the full field target at **397 passed**.
- The remaining field artifact writers now close the same race: model-sensor
  CSV, comparison/validation/source case JSON, and evidence-audit JSON are
  serialized or fingerprinted before exclusive creation. FFI source-state,
  field-screen, comparison, validation, verification, and field-source CLI
  outputs use the same exclusive text boundary, so a competing writer cannot
  replace a provenance-pinned artifact after its preflight check.
- Added a multi-source operational-envelope regression covering two independent
  lower/nominal/upper distributed-source schedules, all nine cross-product
  corners, and each source's injected-mass ledger residual.
- Batch export now renders every case report and the manifest payload before
  creating the target directory. A late JSON/finite-value failure therefore
  leaves no partial batch directory; the batch/report/execution target is now
  **37 passed** and the full field target is **398 passed**.
- Operational screening decisions now carry optional typed `gate_codes` beside
  human-readable reasons/actions. The evaluator emits stable codes for blocked
  physical scope, missing/failed refinement, missing sensor traces, omitted
  obstacles, unresolved uncertainty, incomplete transport, and conditional
  review requirements; legacy direct constructors remain compatible with the
  default empty tuple.
- Cartesian uncertainty/history/sensor/phase aggregate decisions now preserve
  the union of every corner's `gate_codes`, add an explicit conditional-review
  code when aggregate opt-in is absent, and use `corner_withheld` for legacy
  manually constructed withheld corners. Empty sensor/source aggregates carry
  the same machine-readable withheld marker rather than an untyped status.
  The full repository regression after this change is **1,563 passed, 145
  skipped** with no failures.
- Batch `batch-execution.json` case summaries now repeat the manifest
  `gate_codes`, and the execution verifier cross-checks the typed list along
  with status and `screening_allowed`; a compact artifact cannot hide why a
  case was withheld.
- Batch summaries now include deterministic `gate_code_counts`, so a consumer
  can aggregate fail-safe causes without reparsing every case decision.
- Batch execution verification now recomputes the manifest summary from its
  case decisions (including gate-code counts), while accepting older v3
  manifests that predate typed gate codes.
- Fixed-sensor model comparison now emits stable `gate_codes` for sensor-set
  or geometry blocks, basis/operator/temporal mismatches, conditional or
  unresolved executions, and matched bases. Decision-impact records add
  explicit `classification_aligned`, `classification_disagreement`, and
  `model_selection_withheld` codes while retaining the explanatory reasons.
- Matched field-validation scores now carry those comparison codes into the
  validation artifact and add explicit obstacle-evidence and lower-bound
  disposition codes, including `validation_qualified` only for a score with
  no remaining holdback.
- Read-only field-evidence audits now emit typed readiness codes for missing
  or partial channels, incomplete scans, coherent versus cross-file candidate
  layouts, and the invariant `promotion_not_allowed` holdback; older audit
  artifacts without the optional field remain readable.
- Phase/pool operational transport corners now retain specific
  `phase_routing_incomplete` and `pool_transport_withheld` codes alongside
  the general `transport_incomplete` gate, so aggregate fail-safe reports can
  distinguish ledger failure from the pool-to-field handoff boundary.
- The operational gate now maps semi-FV diagnostic holdbacks to typed causes:
  `transport_mass_residual_exceeded`, `source_ledger_mismatch`, and
  `obstacle_transport_conditional`, while retaining the detailed residual and
  scope text in the report.
- Added `degali field-verify` and `verify_field_execution_artifact()` for
  source-schedule, model-comparison, matched-validation, field-screening and
  batch execution JSON. The verifier checks execution input/SHA, reopens
  referenced strict cases/CSVs, checks batch manifest/report hashes, and
  recomputes deterministic reports where the execution options are complete;
  nominal/direct atmospheric-source screening now records those options for
  replay. Complex phase/history/sensor-array screening and batch records
  explicitly report `report_recomputed=false` when only artifact integrity can
  be rechecked. Verification always keeps
  `promotion_allowed=false`; it is not a validation or approval decision.
- Batch execution verification now parses every hash-matched report and
  requires the typed screening or operational-uncertainty-envelope schema to
  agree with the manifest execution kind. Updating a digest cannot make an
  arbitrary JSON object look like a valid case report.
- The same verifier now cross-checks the execution summary and compact case
  status/`screening_allowed` fields against the signed manifest, and screening
  execution verification rejects unsupported report or operational-decision
  schemas even when an older artifact is not replayable.
- The post-hardening execution/input/scope regression set is **84 passed**;
  the complete field/evidence-audit/semi-FV/obstacle-wake target set remains
  **397 passed** with no failures.
- A complete repository regression run now passes **1,561 tests with 145
  skips and 0 failures**, covering both the preserved free-field/reference
  paths and the field deployment additions.
- A fresh 2026-10-07 read-only audit of the connected folders found 15,714
  structured files in `SLABx_LH2`; 33 files expose candidate quantitative
  fixed-receptor rows, but source/weather/obstacle/common-clock channels remain
  missing, so the audit is `partial` and no validation evidence is promoted.
  The support-centre folder still has 20 structured files with no required
  channel, and the related-material folder has none. The current v4 audit
  artifacts retain these separate dispositions; v2/v3 remain available as
  historical pre-multi-sheet-normalization scans.
- `degali field-verify --require-recomputed` now exits 2 when an execution
  record could only be integrity-checked (for example a complex legacy
  screening or batch artifact), allowing automation to fail closed instead of
  treating a non-recomputed report as a fresh deterministic replay.
- Strict field JSON readers now reject duplicate or whitespace/case-normalized
  object keys before schema mapping. The guard covers screening cases,
  comparison/validation cases, source schedules, and batch inputs/artifacts,
  preventing last-value-wins JSON parsing from changing a fingerprint, source
  boundary, applicability flag, or operational decision.
- Field-evidence XLSX sampling now rejects empty or duplicate normalized header
  fields before collecting channel candidates. Ambiguous auxiliary sheets are
  retained as sheet-level diagnostics while valid sheets remain inspectable;
  an all-ambiguous workbook stays unreadable, so no dictionary overwrite can
  create a readiness channel.
- Field source-schedule, model-sensor, and matched-validation CSV imports now
  use exact row-width checks rather than `DictReader`'s implicit missing/extra
  field behavior. Empty-data files, truncated rows, and extra-value rows fail
  with row numbers; headers are also rejected when their stripped names collide.
  This keeps source rates, fixed-sensor predictions, and observed concentrations
  from acquiring a different logical mapping through whitespace or malformed
  rows.
- Field-evidence CSV inventory now rejects empty/duplicate normalized headers
  and rows whose field count differs from the header. Such files remain
  unreadable diagnostic candidates, preventing silent `DictReader` column
  overwrite or extra-value loss from producing a readiness channel.
- Field-evidence CSV scanning now validates row widths in a streaming pass and
  retains the full row count while sampling at most the first 200 data rows for
  bounded value/co-occurrence checks; large historian exports no longer need to
  be materialized in memory for readiness triage.
- Field-evidence JSON scanning now rejects duplicate or normalization-colliding
  object keys before parsing, preventing a last-value-wins overwrite from
  changing a readiness channel silently.
- Every field-evidence candidate now carries the lowercase SHA-256 of the
  scanned file bytes, allowing stale audit output to be detected without
  promoting the inventory into validation evidence.
- Added strict `read_field_evidence_audit_json()` /
  `write_field_evidence_audit_json()` artifact helpers. They use the CLI
  execution schema, refuse overwrite, and recheck all candidate digests before
  accepting or returning a saved audit; unpinned or root-escaping candidates
  are refused.
- Audit execution records now preserve and cross-check the `max_files` scan
  limit, keeping partial-scan status semantics tied to the bound that produced
  them.
- Added `degali field-audit-verify`, which rechecks a saved audit artifact's
  candidate files and supports `--require-complete` exit-2 gating for
  automation.
- Measured-history CSV import now rejects empty/duplicate headers and rows with
  extra or truncated values before any P/T/flow/phase source-quality gate,
  preventing `DictReader` overwrite or silent field loss at the historian
  flash boundary.
- `field-audit` no longer promotes header-only CSVs or empty JSON row arrays
  into source, weather, obstacle, receptor, or clock readiness channels. They
  remain diagnostic candidates with an explicit no-data reason, preventing an
  empty export from producing `candidate_complete` triage.
- The bounded XLSX audit sample now requires a non-empty row after the first
  non-empty header row before promoting any readiness channel; header-only
  workbooks remain diagnostic and retain an unknown row count.
- Blank CSV rows and empty JSON records are also withheld from channel
  detection, so a nominal row count cannot by itself create readiness evidence.
- Common placeholder values (`N/A`, `unknown`, `unspecified`, and
  `not-recorded`) are excluded from identifier and measurement checks in the
  field audit.
- Field-audit channels now require the identifying and physical fields to
  co-occur in one bounded CSV/JSON/XLSX record; values from separate rows are
  never combined into a synthetic source or receptor observation.
- JSON event metadata is carried only into its own bounded observation records
  for this co-occurrence check; unrelated files and row groups remain
  separate candidates.
- Field-evidence audit records now expose immutable
  `coherent_candidate_paths`, identifying same-file candidates that contain all
  five channels while keeping `promotion_allowed=false` and manual identity
  reconciliation mandatory.
- Typed field model sensor sets can now be exported through
  `write_field_model_sensor_set_csv()`. The fixed-sensor artifact carries the
  declared averaging-time column when present, returns SHA-256/path/row-count
  provenance, and refuses to overwrite an existing file; comparison basis IDs
  remain explicit rather than being inferred from the CSV.
- Strict field model-comparison cases may now pin each CSV with paired
  `csv_sha256` and `csv_row_count` fields. The reader verifies both against the
  resolved artifact before comparison, rejecting post-authoring mutation or
  incomplete fingerprint declarations.
- Comparison reports now mark whether CSV provenance was integrity-pinned by
  the case (`integrity_pinned`), keeping an observed fingerprint distinct from
  an expected-hash verification.
- Fixed-sensor comparison CSV imports now reject empty/duplicate headers and
  rows with more values than the header, preventing `DictReader` column
  overwrite or silent extra-field loss at the model-comparison boundary.
- CSV role columns for sensor ID, coordinates, prediction and averaging time
  must now be distinct and single-line, preventing one physical field from
  being interpreted simultaneously as geometry and concentration metadata.
- Added `field_model_comparison_case_record()` and
  `write_field_model_comparison_case_json()`. Typed comparison cases now have
  a fail-safe artifact bridge that preserves custom CSV column mappings,
  relative paths, basis IDs and fresh CSV fingerprints, then re-reads the
  emitted JSON to verify the round trip.
- Added `write_field_source_schedule_json()`, which materializes a typed
  atmospheric source schedule as a strict relative-path case, revalidates its
  source CSV and SHA-256, and re-reads the emitted JSON before returning.
- Added `field_validation_case_record()` and
  `write_field_validation_case_json()`. Matched validation cases now preserve
  model/observed CSV column mappings, evidence identity and relative paths,
  revalidate both CSVs against their typed records, and round-trip through the
  strict parser before returning.
- Matched validation CSV imports now reject empty/duplicate headers, extra row
  values and collisions among sensor, geometry, time, concentration, clock and
  obstacle column roles before scoring.
- AIJ Case-H and paired SMEDIS obstacle-wake observations now enforce a
  fail-safe physical-row boundary: coordinates must be finite, missing
  sentinels remain explicit, and negative/infinite concentrations, RMS values,
  or turbulent kinetic energy are rejected before a neutral/dense-gas
  benchmark can be inspected. This improves evidence quality triage without
  promoting either benchmark to an LH2 wake closure.
- Direct fixed-sensor model-comparison records now require immutable tuple
  collections, typed comparison bases/predictions/rows, boolean-safe CSV row
  counts, and blocked comparisons with no fabricated rows or disagreement
  counts. Decision-impact reasons and alert identifiers follow the same
  immutable boundary.
- Direct matched field-validation datasets now require tuple observations,
  evidence row-count/path agreement, and a complete observed-model comparison
  side. Validation scores reject fabricated model identities and conditional
  statuses without reasons before serialisation.
- Generic distributed vapour sources and `FieldSemiFVRequest` now reject
  mutable list-valued position, warning, obstacle, sensor and source-corner
  collections at construction, preserving deterministic source/geometry
  coverage after a request has been created.
- Operational screening decisions now require immutable reason/action tuples,
  and a directly constructed `conditional_allowed` decision must retain a
  non-empty explanatory reason and review action before it can be serialized.
- Direct operational uncertainty, measured-history, sensor-array,
  source-sensor, and phase-routing envelope records now cross-check their
  aggregate decision against every corner. Fabricated allowed statuses,
  missing aggregate refinement/uncertainty requirements, and corner decisions
  that waive those requirements are rejected before report serialization;
  individual allowed corners also require completed screening and refinement
  results.
- Phase-routing operational aggregates now require tuple sensor/case
  collections and at least one retained corner; an empty sensor-array
  aggregate can only be withheld and cannot claim resolved evidence.
- `FieldBatchCase` now cross-checks physical applicability, operational status,
  refinement completion/availability and review metadata, so a manifest cannot
  pair an allowed decision with blocked/conditional physical scope or claim
  refinement that was not requested or completed.
- `FieldBatchExport` now verifies that manifest/report artifacts exist inside
  the declared output directory and that report filenames match their case
  labels before the export object can be returned or serialized.
- `FieldBatchExport` now parses those artifacts and cross-checks manifest schema,
  case order/identity, applicability, refinement, review, operational decision
  summary and report schema before accepting the typed export.
- Batch manifests now fingerprint every emitted case report, and
  `batch-execution.json` records both those report hashes and the manifest hash
  so later consumers can detect post-export artifact mutation.
- `ReleaseSource.metadata` and `StabilityScalarMixingClosure.diffusivity_m2_s`
  now snapshot caller mappings into immutable, type-checked mappings, so
  provenance and declared stability diffusivities cannot change after a case is
  audited.
- `FieldPhaseRoutingConfig` now snapshots gas/phase model-option mappings as
  immutable mappings, preserving phase-ledger controls after construction.
- Direct `ReleaseSource` construction now enforces `m` units for location
  uncertainty and dimensionless `1` units for direction uncertainty, matching
  the strict case-file parser before any corner is transported.
- Strict batch input now requires an immutable tuple of cases even for direct
  construction, and rejects unsafe case labels before loading referenced
  artifacts or deriving report paths.
- `FieldLH2FlashResult` and `FieldLH2SourcePreparation` now validate their
  nested flash/source types, phase flows, residual signs, effective geometry,
  warning strings and JSON-safe closure diagnostics at construction. A blocked
  preparation cannot expose a flash plane to downstream transport.
- Completed DEGALI field traces can now be exported directly into the matched
  fixed-sensor model-comparison contract with an explicit peak/final/full-trace
  operator on the true or indicated channel. Withheld sensors and blocked
  executions fail closed, and execution applicability/operator provenance keeps
  conditional runs from being reported as accepted model comparisons.
- Completed source superpositions can now use the same model-comparison export
  boundary, retaining their conditional source-interaction provenance. Sensor
  superposition also rejects steady/transient branch mixing before conservative
  traces are combined.
- Sensor superposition now inspects branch source/sensor/ambient/obstacle and
  distributed-source bounds plus stability, phase-routing and historian
  resolution flags; unresolved labels remain in warnings and set
  `uncertainty_complete=False` until the corresponding deterministic envelope
  is run.
- `FieldSemiFVEnvelope` now rejects duplicate corner selections before
  operational envelopes consume it, so deterministic uncertainty coverage
  cannot pass through a duplicate/missing-corner set comparison.
- Measured-history flash and field-screening envelope records now validate
  bound-selection keys, typed schedules/results, warning strings and unique
  source/joint corner selections before transport/report construction.
- Phase-routing results, phase/transport corner envelopes, droplet-handoff
  refinement outputs, historian imports and batch-export manifests now reject
  malformed types, non-finite values, duplicate selections and incomplete
  case coverage at construction, before ledger or report serialization.
- Fingerprinted atmospheric source-schedule envelopes now apply the same
  zero-endpoint contract to direct lower/upper rate-bound constructors as the
  CSV importer, preventing a bound history from continuing past its declared
  clock; bound-column names are validated before transport.
- Field screening results and deterministic envelopes now require tuple-based
  sensor/projection/corner collections at construction, so mutable list inputs
  cannot bypass deterministic case/report contracts.
- Operational uncertainty-envelope reports now include a deterministic sensor
  extrema summary (peak/final/time-average true and indicated signals) for
  every declared detector. Bounds are computed only from completed corner
  traces; missing or off-plane traces remain explicitly counted as withheld
  with their reasons, and the summary is never labelled a probability interval.
- The same deterministic sensor-extrema contract is now present in the joint
  measured-history, source-plus-sensor, and phase/pool operational reports;
  phase cases with no executable screening trace remain withheld instead of
  disappearing from the aggregate.
- Standalone detector-calibration (`--sensor-array-envelope`) operational
  reports now emit the same per-detector extrema and withheld accounting, so
  every CLI uncertainty branch has one auditable sensor summary shape.
- `ReleaseSource`, `WeatherState`, `FieldSemiFVRequest` and distributed-source
  geometry boundaries now reject boolean-valued numeric components at
  construction, preventing Python `True`/`False` coercion from becoming a
  physical coordinate, duration, ambient state, tolerance or source width.
- LH2 blowdown and measured-history adapters now fail closed on negative phase
  mass flows and on direct-vapour/liquid mass-ledger mismatch; the previous
  non-negative remainder fallback no longer hides an over-released schedule.
  Sensor and semi-FV boundaries also reject non-finite/boolean averaging or
  source-height inputs and non-positive detector gain.
- Operational screening no longer treats a nominal direct post-flash schedule as
  resolving unrelated release geometry uncertainty. Source location/direction
  bounds remain gated (and finite release timing remains explicit) while the
  schedule replaces only the upstream thermodynamic/rate boundary. A direct
  schedule with unpropagated source geometry is therefore withheld instead of
  being conditionally allowed from a nominal plane.
- Direct `FieldModelComparisonCase` construction now rejects boolean threshold
  and position-tolerance values instead of coercing them to `1.0`; parser and
  comparison-function numeric contracts are now consistent.
- Common `BoundedValue` and `CircularBoundedValue` Python boundaries now reject
  boolean values instead of coercing `True` to the physical number `1.0`.
- Direct field source, weather, surface, sensor, scenario and request records
  now reject untyped nested boundaries at construction instead of failing later
  through attribute access during screening.
- The low-level atmospheric schedule, obstacle, distributed-source and
  semi-FV numerical/refinement boundaries likewise reject boolean coercion
  before solving.
- `SemiFVConfig.gravitational_settling_m_s` now has an explicit positive-
  downward sign convention; the finite-volume face flux is corrected to move
  positive settling toward decreasing height, with a regression for the
  vertical centre-of-mass direction.
- Matched field-validation scoring and `FieldValidationCase` now apply the same
  boolean threshold/tolerance rejection, so validation gates cannot be changed
  by Python truth-value coercion.
- Direct comparison/validation case records now reject untyped model or dataset
  objects before report generation rather than relying on later attribute access.
- The local FFI Test 6 observation-operator audit now emits separate aggregate
  `arc_max_table` and `sensor_height_table` results. The first is the finite
  arc maximum at each selected radius; the second exposes exact-coordinate
  projection grouped by reported sensor height. Raw sensor rows remain
  external, lower-bound semantics are explicit, and the legacy aggregate
  30 m metrics are retained for compatibility.
- FFI source-state envelope reports now carry the same arc-max and
  sensor-height observation operators for every deterministic source/weather
  corner. Available and withheld sensor counts, lower-bound satisfaction and
  the exact corner selection remain attached; withheld predictions are never
  zero-filled, mixed groups remain `partial`, and operator diagnostics do not
  qualify validation.
- Field semi-FV requests now support explicit bounded ambient temperature,
  pressure and air-density boundaries. Ambient corners are propagated through
  flash, transport and sensor conversion; reports retain the boundary and
  operational screening withholds a nominal result until those bounds are
  resolved. The phase-routing/pool handoff repeats its phase ledger at each
  ambient temperature/pressure corner and carries air-density corners into
  field conversion; measured-history operational runs remain fail-safe until
  their time-aligned history and ambient corners are both resolved.
- Reused LH2 saturation tables now include both endpoint saturation
  temperatures implied by a declared ambient-pressure bound, preventing an
  otherwise valid pressure corner from being rejected by a nominal-only table
  domain.
- A fingerprinted `FieldValidationEvidence` record now keeps applicability
  `conditional` until the separate matched observed-vs-model validation score
  is produced; dataset identity alone can no longer look like model
  qualification.
- Matched field validation now accepts an explicit optional observation-kind
  column (`exact` or `lower_bound`). Censored lower-bound rows are never
  folded into MAE/RMSE/bias: reports expose satisfaction fraction and deficit,
  satisfied censored datasets remain `conditional`, and violated or unmatched
  bounds are `withheld`.
- The field LH2 flash handoff now exposes mass-partition, momentum and energy
  residuals in screening reports. A non-finite or out-of-tolerance residual
  blocks field transport instead of allowing the semi-FV scalar solver to
  consume an unclosed source; the post-flash liquid/vapour split remains
  explicit and unrouted liquid is not silently evaporated. Measured-history
  and cryogenic-blowdown schedule adapters apply the same closure gate before
  emitting a time-varying atmospheric source. The report records the exact
  mass/phase/momentum/energy tolerances used by that gate, and a blocked
  source retains JSON-safe residual diagnostics without exposing an unclosed
  flash plane to transport.
- Added explicit bounded release-duration uncertainty to the field source
  contract and strict case parser. Timing lower/upper corners now propagate
  through the field envelope, while fixed measured-history schedules reject an
  unresolved duration bound instead of reusing a nominal time axis.
- Field requests now reject anonymous direct post-flash schedules with the
  default `source_id="declared"`; measured-history and pool adapters retain
  their explicit provenance IDs.
- Added bounded Cartesian source-location coordinates to the field release
  contract. Each lower/upper geometry corner is reprojected and transported;
  mixed numeric/bounded coordinate declarations are rejected.
- Added explicit bounded geometry for global axis-aligned and oriented field
  obstacles. Every geometry corner is reprojected and transported with its
  declared evidence ID; unresolved nominal obstacle bounds remain withheld and
  never imply a probability distribution or wake closure.
- Added bounded source-direction components with explicit non-zero-corner
  validation. Direction corners are retained for jet/phase handoff paths but
  do not create an unvalidated momentum effect in the direct scalar closure.
- Added `FieldValidationEvidence` for matched LH2 validation datasets. The
  record fingerprints the external artefact and names its common clock,
  source/weather/obstacle/receptor geometry and temporal operator; a bare
  `lh2_validation_available=True` flag is now fail-closed in strict input and
  does not qualify a free-field record for obstacle transport.
- Validation datasets now reject source, weather, receptor-geometry or
  temporal-operator identifiers that disagree between the evidence record and
  the comparison basis. Physical obstacle IDs remain a separate namespace
  from model mask/wake representation IDs and are checked by the scoring gate.
- Semi-FV transport reports a source-level injected-mass ledger for the primary
  release and each evidenced distributed source, preserving the conservation
  audit when multiple source schedules are combined.
- Field screening now checks the source-ledger residual independently of the
  cell inventory residual; an inconsistent per-source bookkeeping total is
  retained diagnostically but fails closed before operational promotion.
- The read-only evidence audit now requires an explicit weather identity before
  promoting wind data to the weather channel; unjoined wind/temperature fields
  remain diagnostic notes and cannot make a candidate look complete.
- Model-comparison CSV provenance now records whether an averaging-time column
  was actually present. A declared averaging operator that cannot be checked
  against the imported CSV leaves numerical rows available but withholds the
  model-selection impact.
- Distributed atmospheric handoffs now reject unknown `source_kind` values;
  only declared post-flash, pool, droplet, and generic atmospheric categories
  can enter the local scalar transport.
- Model-sensor sets now reject CSV provenance whose row count or type does not
  match the supplied prediction set, preventing metadata-only fingerprints from
  being mistaken for a complete comparison input.
- Matched validation imports now require the evidence artefact path to resolve
  to the observed CSV path as well as matching its digest and row count; a
  copied file with the same bytes cannot silently acquire another dataset's
  provenance label.
- The conservative semi-FV obstacle solver now rejects obstacles whose top
  extends above the declared vertical domain instead of silently truncating
  the geometry in the solid-cell mask.
- Global obstacles crossing the local z=0 datum are now returned as an
  explicit blocked projection (and fully below-datum geometry as a diagnostic
  omission) instead of escaping the field workflow as an uncaught geometry
  exception or being silently clipped.
- Operational screening decision records now validate their status, boolean
  gate flags, non-empty reasons/actions, and required holdback explanation at
  construction time; malformed status objects cannot be serialized as an
  apparently safe screening result.
- Field applicability records now apply the same fail-safe type and
  non-empty-text checks to their uncertainty flag, reasons, and warnings
  before they enter reports or operational gates. A blocked record cannot claim
  complete uncertainty, and a non-blocked record cannot carry blocking reasons;
  accepted physical results may still retain explicit conditional warnings.
- Field screening and sensor-result records now validate their direct API
  structure at construction time: malformed typed subrecords, duplicate sensor
  labels, invalid obstacle projections, and trace-bearing results without a
  transport result are rejected while the operational gate retains its own
  missing-sensor defense.
- Direct post-flash field schedules now require a zero final endpoint at the
  Python request boundary as well as the strict CSV boundary, so a finite
  release cannot silently acquire an ambiguous post-duration source rate.
- Plain field uncertainty envelopes, sensor-calibration envelopes and
  refinement results now reject malformed direct records before report
  serialization, including empty corner sets, duplicate/invalid corner keys,
  non-finite selections and untyped nested results.
- The public model-comparison function now rejects non-comparison-set inputs
  and boolean threshold/tolerance arguments before they can be coerced into
  numeric comparison operators.
- Model-selection impact records now validate status, alert-sensor identity,
  monitored-distance finiteness, runtime ratio and consistency with the matched
  comparison before they can be serialized as a decision-facing result.
- Operational uncertainty aggregates now preserve the logical AND of every
  corner's uncertainty-resolution flag; a withheld corner cannot be reported as
  uncertainty-complete merely because the aggregate is already withheld.
- The phase-routing/pool operational aggregate applies the same corner-wise
  uncertainty flag rule, keeping refinement holdbacks distinct from unresolved
  physical uncertainty in phase-transport reports.
- The measured-history operational aggregate now preserves that same
  corner-wise uncertainty flag in withheld, conditional-review and allowed
  paths instead of hardcoding the aggregate as uncertainty-complete.
- Phase-routing operational corner records now validate finite/unique phase and
  sensor selections, typed nested screening/refinement/decision records, and
  reject an allowed decision when the field screening record is absent.
- General and measured-history operational corner records now validate unique,
  finite, report-safe selections and typed envelope membership; an allowed
  corner cannot omit its required refinement result.
- Phase-routing operational corners now apply the same required-refinement
  guard even when a screening result object is present.
- Phase-routing operational envelopes now validate every declared sensor
  selection (finite numeric values, unique keys and unique corners) before
  comparing the physical and calibration corner sets.
- Sensor-calibration uncertainty envelopes now reject malformed corner pairs,
  duplicate corners, and empty case sets when the shared field transport
  completed; completed reports cannot masquerade as calibration-free results.
- Added `run_field_operational_sensor_array_envelope()`, which reuses the
  shared scalar transport, attaches exact detector-calibration corners,
  applies the numerical-refinement/operational gates to each corner, and
  emits a versioned aggregate report without waiving source or weather
  uncertainty.
- `degali field-screen --sensor-array-envelope` now emits the same auditable
  detector-calibration operational envelope and rejects mixing it with source,
  historian or phase uncertainty modes.
- Added the joint source–sensor operational envelope and strict
  `--joint-source-sensor-envelope` CLI path. Physical corners are solved with
  nominal supplemental calibration, then crossed with exact calibration
  traces; every Cartesian corner still receives the ordinary refinement,
  receptor and resolved-uncertainty gates.
- Added the `degali.validation.ffi_source_state` source-state envelope for
  public FFI/Spadeadam horizontal releases. It crosses explicitly evidenced
  source, wind and ambient bounds, evaluates the existing free-field jet at
  reported sensor coordinates, and emits a lower-bound-aware residual map
  without fitting coefficients or claiming operational validation.
- Exposed the same FFI source-state audit through the strict
  `degali ffi-source-state` CLI, including bounded-corner parsing,
  versioned execution JSON and `--require-complete` fail-safe exit handling.
- FFI source-state execution records now fingerprint the `conditions.csv` and
  `sensors.csv` inputs and distinguish completed numerical corners from the
  always-false operational/validation qualification flags.
- FFI residual reports now add censor-aware one-sided lower-bound satisfaction
  and worst-deficit metrics, without replacing the descriptive symmetric
  residual fields or treating them as an accuracy score.
- Field uncertainty envelopes now preflight distributed-source schedule
  duration against every primary release-duration corner and reject a
  time-misaligned combination with an explicit remediation message, instead
  of failing midway through corner construction.
- Model-comparison basis, prediction-row, and comparison-result records now
  validate numeric/boolean types and internal concentration, ratio, geometry,
  sensor-coverage, disagreement-count, and mean-difference invariants before
  a report can serialize them; malformed direct API objects fail closed.
- Generic distributed atmospheric sources now accept paired same-clock lower
  and upper rate histories. Their deterministic rate corners are propagated
  with source geometry/width corners, retained in reports, and treated as
  unresolved by nominal operational screening until the envelope is run.
- Operational screening decision records now reject any directly constructed
  allowed status whose required refinement evidence is unavailable or whose
  required uncertainty remains unresolved; the gate cannot be bypassed by
  serializing a hand-built status object.
- Matched validation score records now validate dataset/model/comparison
  identity, finite and non-negative error metrics, row-derived MAE/RMSE/bias/
  maximum values, and qualified-status consistency before serialization.
- Sensor traces and multi-source superposition records now validate finite,
  strictly ordered common time axes, equal array lengths, non-negative scalar
  concentration, true mole-fraction consistency, and warning/type contracts
  before a combined field report can be serialized.
- The 2026-10-06 read-only evidence re-audit now records the stricter
  weather-identity rule in the handover/readiness notes: the connected
  `SLABx_LH2` directory remains `withheld` (15,604 structured files), as do
  the 20 support-centre files and the zero-file related-data folder; no
  validation evidence or site acceptance is promoted from these scans.
- Added the strict `degali.field-validation-input.v1` / `field-validate` path
  for one-per-sensor observed summaries. It verifies the observed CSV digest,
  row count, common clock, obstacle ID and averaging operator before reporting
  matched-sensor MAE/RMSE/bias; it never interpolates an unmeasured hazard
  distance.
- Added the read-only `field-audit` evidence readiness scanner. It inventories
  structured external CSV/JSON files for source, weather, obstacle, receptor,
  and common-clock channels while keeping promotion fail-safe; it never creates
  a validation evidence record from filenames or pool-radius outputs.
- A bounded field-audit scan now returns `scan_complete=false` with an explicit
  withheld/partial reason instead of promoting a truncated directory result.
- Sensor-summary CSVs with only receptor geometry or threshold crossing times
  are retained as non-quantitative audit notes, never as concentration
  observations.
- Added the strict `field-source` atmospheric schedule adapter. It fingerprints
  an already-post-flash/pool/droplet-vapour rate CSV, binds source and common
  clock IDs, requires a zero final endpoint, and refuses phase or leak-rate
  inference from upstream plant tags.
- The same adapter now preserves paired lower/upper rate columns as deterministic
  same-clock source corners instead of treating a nominal schedule as exact.
- `run_field_semi_fv_envelope(..., atmospheric_source_schedule=...)` now combines
  those source-rate corners with weather/surface/sensor corners, fixes replaced
  upstream source parameters at nominal values, and retains schedule provenance
  in the envelope report.
- `run_field_operational_uncertainty_envelope(...,
  atmospheric_source_schedule=...)` now carries those same source-rate corners
  through per-corner numerical refinement and fail-safe operational decisions.
- Atmospheric schedule envelopes now retain declared source location/direction
  uncertainty corners instead of collapsing all source geometry to nominal; a
  non-exact release-duration bound is rejected unless time-aligned schedule
  corners are supplied.
- Added an operational-envelope regression that verifies those source-geometry
  corners remain separate after atmospheric schedule replacement.
- Generic `FieldDistributedVapourSource` inputs now accept explicit bounded
  global position and vertical-width values. The semi-FV envelope expands these
  source-shape corners, operational screening detects unresolved nominal bounds,
  and reports retain the original evidence-backed intervals.
- Strict `field-screening-input.v1` cases can now carry those explicit
  distributed atmospheric sources with unit-labelled schedules; arbitrary
  imported schedules and phase inference remain rejected.
- `field-screen --atmospheric-source-case` composes a separately fingerprinted
  atmospheric schedule only at the complete uncertainty-envelope boundary;
  its input path/SHA and schedule provenance are retained, while historian and
  phase-routing source paths remain mutually exclusive.
- The direct atmospheric-source envelope now accepts only
  `post_flash_atmospheric_vapour`; pool/droplet schedules are rejected unless a
  location-bearing distributed or phase-routing handoff is supplied.
- Atmospheric source CSVs with zero integrated release mass are rejected at the
  import boundary instead of producing a misleading zero-concentration screen.
- `field-audit` now samples bounded XLSX worksheet headers/early rows while
  leaving workbook row counts and cross-sheet identity unresolved; promotion
  remains fail-safe and still never creates validation evidence.
- Named field batches can now route selected cases through
  `atmospheric_source_schedules`, exporting the complete operational
  uncertainty-envelope schema with source provenance and all deterministic
  source/weather/sensor corners. The batch path rejects attempts to waive its
  refinement, in-plane-sensor, or resolved-uncertainty gates. A named-case
  `conditional_review` authorization now applies consistently to this envelope
  path as well; every deterministic corner must still pass the non-negotiable
  gates before the aggregate can become `conditional_allowed`.
- Batch `manifest.json` and `batch-execution.json` now carry the same typed
  `degali.field-batch-summary.v1` fail-safe summary: case count,
  decision-status counts, screening-allowed and
  withheld counts, plus an `all_screening_allowed` flag derived from the
  operational decisions rather than physical applicability labels.
- A directly constructed batch case cannot claim `conditional_allowed` unless
  its accountable conditional-review authorization is present and marked as
  applied; the invariant is enforced before summary or manifest serialization.
- The read-only field-evidence audit now validates candidate/channel/reason
  tuple contracts, duplicate channel paths and contradictory `partial` states
  before a readiness record can be serialized; malformed audit objects fail
  closed instead of looking like a promotion candidate.
- Added the strict `degali.field-batch-input.v1` / `field-batch` CLI boundary.
  It fingerprints every referenced case and source JSON, emits a separate
  `batch-execution.json`, rejects historian/phase-routing mixing, and supports
  fail-safe exit 2 when any named case remains withheld.
- Added the opt-in `LH2_COUPLED_TRANSIENT` research path. It conservatively
  routes direct flash vapour through the finite jet/plume/puff model and
  residual liquid through droplet rainout and the concurrent pool model.
  In-flight and pool evaporation remain explicit time-resolved atmospheric
  source terms until defensible launch closures are supplied; they are not
  folded into the nozzle jet.
- Added global lateral source offsets to the yawed finite-release path without
  changing its default origin or the frozen DEGADIS 2.1 model.

## 0.2.0 - 2026-10-03

- Added a finite-duration release path with automatic steady-to-puff handoff,
  transient receptor histories and explicit mass-conservation diagnostics.
- Added an opt-in axisymmetric dynamic pool with substrate heat transfer,
  cryogenic evaporation, droplet rainout coupling and validation helpers for
  the JUEL and Sandia liquid-hydrogen pool experiments.
- Added finite wall-jet transition, yawed crosswind, site-geometry and obstacle
  wake screening paths while preserving the existing DEGADIS-compatible
  defaults.
- Added cryogenic blowdown, mixed-phase and thermal-property extensions,
  including a cached CoolProp lookup-table path for repeated evaluations.
- Expanded the public validation portfolio and documented the evidence
  boundaries for PRESLHY, ELVHYS, FFI Test 6 and open-channel hydrogen data.

- Added an identical-sensor Trial 10/23 audit of the coefficient-free direct-
  flux thermal-moment path. Both step levels reach 4 m and pass conservation,
  terminal and 40-sensor refinement gates, but thermal centre-amplitude error
  worsens while variance error improves. The preregistered observation gate
  therefore fails and no research or default path is promoted.
- Added an opt-in `radial_core` constitutive-validity domain after diagnosing
  that the full-circle cumulative diffusivity identity was being evaluated on
  square-clipped outer contours. It changes no full-square conserved integral,
  source or side flux, and does not clip a diffusivity. The existing
  `full_square` behavior remains default.
- Corrected the independent-width ground-image receptor to superpose direct
  and reflected species and enthalpy fields separately, avoiding a nonphysical
  image cross term when the thermal and species widths differ.
- Completed the preregistered direct-flux thermal-moment extension for PRESLHY
  Trials 11, 12, 22, 23, 24 and 25. All six linear-phase paths reach their
  frozen downstream stations at both step levels, retain positive sampled
  diffusivity, close independently reconstructed order-16 balances below
  `5.05e-11`, and pass the `0.005` terminal convergence gate. No downstream
  sensor score or default-model promotion is implied.
- Replaced the legacy `1e-4`-tolerance width root only inside the independent
  Gaussian-enthalpy research path with the cancellation-safe analytic positive
  root of the same JETPLU width constraints. The core DEGADIS-compatible
  splitter remains unchanged. This removes the numerical width plateau that
  stopped Trials 12 and 24.
- Added an opt-in, nodally exact C1 Hermite phase-property lookup. An independent
  8,000-point oracle screen reduced temperature RMS lookup error from
  `1.15e-4` to `4.19e-5 K` with no nonnegative density slopes. The package
  default remains the established bilinear lookup.
- Added a pre-registered classical-RK4 path that transports total mass,
  hydrogen mass, vector momentum, total energy and the physical thermal second
  moment directly, reconstructing a bounded physical section at every stage.
  Trial 10 reaches x=1.78 m at both 0.005/0.0025 m steps; independent order-16
  balance errors are below `4.35e-11`, and terminal state/flux differences are
  below `1.23e-4`. This removes the prior state-space accumulation error without
  fitting a new diffusivity or promoting the research closure to a default.
- Added a guarded adaptive downstream marcher for the opt-in conservative
  thermal second-moment state. Trial 10 reaches the first profile station and
  passes step-result convergence, but both independent long-range balance
  checks fail the preregistered limit. No receptor score or promotion is made;
  direct six-flux/moment marching is the next required numerical path.
- Added a pre-registered Trial 10 measured-wind displacement audit using the
  identical five-sensor thermal/H2 operator. Directional yaw reduces all four
  centre/load errors but increases both variance errors; magnitude-only and
  full-vector paths also fail the componentwise gate. The 0.02/0.01 m vector
  result is converged, so no wind/default path is promoted.
- Added an identical five-sensor model/measurement profile gate. Stored model
  trajectories are now compared with observed thermal and hydrogen centre
  amplitudes and finite-window variances without reintegration or coefficient
  fitting. The tested same-width enthalpy path improves variance slightly but
  worsens centre amplitude, so it is not promoted.
- Added a pre-registered, coefficient-free observation operator for collocated
  PRESLHY Trial 10/23 temperature-deficit and hydrogen profiles. The public
  code computes finite five-sensor zeroth/centroid/second moments and lag
  diagnostics, while the external workbooks and generated audit JSON remain
  local and excluded from distribution. No model coefficient or default was
  changed.
- Renamed the project, Python package and command from `degadisx` to `degali`.
  DEGALI means **Dense Gas Dispersion for Liquid Hydrogen**. This is a naming
  migration only; it does not change the frozen physical or numerical results.
- Prepared a compact public repository while retaining the full local research
  ledger. Third-party PDFs, raw workbooks, checkpoints and large generated
  fields are excluded without deleting them.
- Replaced the development-log front page with a concise status, installation,
  evidence and safety overview. The former front page is preserved under
  `docs/README-development-log-2026-09-06.md`.
- Added publication, contribution and security policies plus an automated
  repository size, restricted-file and credential check.
- Version 0.1.0 was the initial public alpha snapshot. Research-only paths in
  this release remain opt-in unless their documentation explicitly says
  otherwise.

## Research history

- Added joint PSD bounds showing that TKE and axial normal-stress work cannot
  both be assumed small for the frozen7 initial shear fields; all7 resolved.
  Added explicit normal momentum/work source allocation. Preserved fixed-center
  Q2 failure; new exact-geometry six-DOF mobile initialization passes both .02/2
  synthetic cases without changing source targets or fitting Q/r (max1.489e-11).
  Mobile source6 tests passed,105 hashes match. Added straight-frame stress/
  production transformation primitives; no full normal/pressure closure or
  observed/default adoption (2026-09-06,05:47 KST).

- Completed the explicit-Q source energy reallocation audit: six retained
  moments close to5.566e-12 while Q is paid from the unchanged total source
  energy. Original projection preserved;95 hashes match. Four new source
  tests passed separately AFTER full425. Synthetic-only, no actual Q selected
  or physical initial-boundary acceptance (2026-09-06,05:10 KST).

- Implemented opt-in finite Q=rho*k modal transport with thermal dissipation
  D, turbulent production P-D, and total H+meanKE+Q energy. All original scalar
  rows retained. Manufactured phase-refinement and actual energy gates pass
  (max1.992e-7), but its physical input gates FAIL and it is not adopted.
  Full non-slow425 passed. Added explicit-Q boundary source energy reallocation
  helper for separate testing. [Results](docs/finite-tke-transport-results.md)
  (2026-09-06,05:00 KST). Defaults/field scores unchanged.

- Completed exact-geometry conservative aligned witness re-solve and its new
  actual-energy differences (max7.66e-8); all retained equations and independent
  boundary gates pass. Extended necessary shear/TKE bounds to all7 initial
  cases (10.51–11.72% of mean axial KE flux), all numerically resolved.
  Next finite-TKE operator specified, not yet implemented. No physical closure,
  downstream observation or default promotion (2026-09-06,04:30 KST).

- Added opt-in exact-constraint transverse geometry with consistent tangents;
  actual full-energy directional differences pass all four registered checks,
  max7.18e-8. The old direction is not reused as a new conservative solution.
  Added sharp PSD shear/TKE lower-bound diagnostic (no fitted coefficient):
  minimum TKE flux10.56% of mean axial KE for the fixed prior aligned witness.
  Non-slow403 passed;3 later bound tests passed separately. No field/default
  adoption. [Results](docs/exact-geometry-and-tke-results.md) (2026-09-06).

- Added bounded-memory phase operators and a direction-preserving four-mode
  circulation feasibility pilot: failed trial10 edges3.543/4.000/1.963%,
  nonradial stress1.5e-15, all original weak/global rows retained. Not a closure.
  Preserved failed raw/paired/root-polished energy checks; new high-precision
  smooth enthalpy diagnosis separates cancellation and width/wind inconsistency.
  No field/default change. Non-slow385 passed;7 later unit tests separately pass.
  [Results](docs/aligned-flux-and-energy-conditioning-results.md) (2026-09-06).

- Completed seven guarded new-shape segment tests:6/7 pass short intervals,
  trial10's corner heat5.01189% failure independently reproduced. Added
  rate-free compatibility bounds, a rejected boundary-row prototype and
  conservative divergence-free flux modes retaining all original weak rows.
- Hard boundary-moment solenoidal pilots fail; separate small-amplitude
  inequality FEASIBILITY witness satisfies sampled gates but needs independent
  phase-volume verification and a physical flow/stress closure. No downstream
  adoption or observed score. Full non-slow371 passed;3 later LP tests pass.
  [Results](docs/conservative-flux-witness-results.md) (2026-09-06).

- Extended simultaneous initialization to all7 unchanged cases; all pass
  independent local gates. Added degree6/37-state operator coverage.
  New guarded research segment driver evolves fields without endpoint
  refitting and requires explicit mixing-amplitude scaling.12 guard/state
  tests pass; full non-slow355 passed,128 skipped,58 deselected578.82s.
  Downstream tests ongoing; trial10 heat boundary failure is preserved.
  No observed score/default promotion.
  [Results](docs/coupled-initialization-extension-results.md) (2026-09-06).

- Added research-only simultaneous shape/rate initialization: vectorized
  fixed-mesh iterations re-solve all modal rates for every candidate and
  difference sample, with hard six-moment constraints. Final candidates
  are retracted/verified by separate moving-phase quadrature and65 rays.
- Both preselected pilots25/11 pass local gates, including the newly added
  momentum boundary5% requirement. Trial25 scalar defects fall to3.24/2.30%;
  trial11 no longer requires negative inferred viscosity on checked samples.
  No fitted diffusivity, new flow redistribution, default or observed score.
  Added6 tests;342 passed,128 skipped,58 deselected. Other five pending.
  [Results](docs/coupled-shape-initialization-results.md) (2026-09-06).

- Removed redundant thermal-width/shape coordinates in a new research path;
  added23/37-state conditional modal transport with newly reconstructed
  mass/axial-momentum fluxes, stress work, angular scalar diffusion and weak
  reservoir boundaries. The supplied scalar mixing profile remains explicit.
- Preserved the original4/7 numerical screen and independently verified the
  fine mixing input for a combined7/7 numerical result. Constitutive adoption
  fails all seven: renewed edge mismatches and two negative-viscosity cases.
  Higher-order witnesses locate those signs near finite-section corners.
  Added12 tests; non-slow336 passed,128 skipped,58 deselected. No default,
  new downstream segment or field score. [Results](docs/enriched-shape-transport-results.md).

- Added face/phase-intersection angular splitting, rebuilt phase-ray
  quadrature, six-moment derivatives and hard-constraint scalar shape
  refitting. All seven frozen cases now pass independent section gates.
  Internal edge defects are0.95–3.63%, not observed prediction errors.
- Kept prior failed evidence, original shape degrees, log correction limit,
  transport, physical parameters and core/defaults unchanged. Added7 tests;
  non-slow:324 passed,128 skipped,58 deselected. No new shape-coefficient
  transport law, field score or default promotion.
  [Results](docs/edge-conservative-refit-results.md).

- Added opt-in square-symmetric C/H shape enrichment with fixed transport,
  six nonlinear moment constraints, exact phase inversion and independent
  edge/moment verification. Three of seven cases reach both edge gates;
  zero pass all conservation gates, so no candidate is adopted.
- Retained both basis degrees and failures; added a separate phase-cell
  polar integration diagnosis without changing original coefficients.
  Angular reference reproduction still needs refinement in five cases.
  Added16 tests; non-slow:317 passed,128 skipped,58 deselected.
  No new field score, downstream shape transport or default changes.
  [Results](docs/edge-profile-enrichment-results.md).

- Added opt-in ambient-reservoir total-energy flux boundaries and a reduced
  buoyancy-work ledger for weak thermal zeroth/second moments. Added a free
  thermal-width RHS and short research segment driver without changing the
  core, original boundary guards, source coefficients or default options.
- All seven initial weak balances and independent physical flux differences
  pass. Original 2/4/8-step short-segment checks pass 6/7; retained trial 10's
  fixed-step failures and separately verified it with adaptive DOP853.
  Combined short-segment checks pass 7/7, not pointwise/field validation.
  Edge species/heat gradient defects remain 12–25%. Added 14 tests; non-slow:
  301 passed, 128 skipped, 58 deselected. [Results](docs/reservoir-thermal-moment-results.md).

- Added reduced axial shear-stress/work reconstruction and an explicit
  thermal equilibrium compatibility screen. All seven numerical audits pass;
  unity thermal/species diffusion plus immediate shear heating fails the
  independent zeroth heat balance in all seven frozen boundary cases.
- Separated finite-edge advection/diffusion and ruled out merely reducing
  the heat-conversion fraction while keeping the current five source terms.
  Preserved failed physical candidates. Added 14 tests; non-slow regression:
  287 passed, 128 skipped, 58 deselected. No core/default or field-score change.
  See [the results](docs/shear-thermal-compatibility-results.md).

- Added conditional conservative transverse mass/species reconstruction,
  actual-source/coflow replay, a free thermal-width tangent family and an
  inferred anisotropic mixing tensor. Negative diffusivity is rejected;
  thermal/species ratio and width rate are never silently supplied.
- Added monotonic phase-cell-split quadrature and independent finite-change
  verification of all seven boundary tangent families. Retained failures of
  the dense unsplit gradient reference. Added 16 tests for diffusion,
  moving/curved sections, phase splitting and invalid/underdetermined inputs.
  Non-slow: 273 passed, 128 skipped, 58 deselected. No default promotion.

- Added research thermal-shape operators: mean advective enthalpy second
  moment and reduced derivative, planar curved/moving rectangular balance,
  and variable-density specific-enthalpy diffusion response. No diffusivity,
  mechanical heating distribution or downstream beta_H closure is invented.
- Preserved a failed seven-case fixed-order phase-gradient quadrature audit;
  added separate adaptive volume/edge integration with explicit error and
  nonconvergence reporting. Gaussian diffusion and curved manufactured
  balances plus negative input/closure tests are covered by 32 new tests.
  Adaptive verification passes 7/7 fixed sections, including the net response;
  no observed accuracy score or downstream transport closure is claimed.
  Non-slow suite: 257 passed, 128 skipped, 58 deselected. No default change.

- Added a boundary-only independent enthalpy width constrained by five fluxes
  and the near-field buoyancy moment. All seven fixed boundaries pass both
  conservation/force and existing temperature/width gates, with independent
  adaptive GK21 confirmation. No new field accuracy score or promotion.
- Added implicit phase/moment derivatives, fixed-flux warm starts and a
  force-cancellation rejection regression. Unclosed downstream transport and
  inappropriate five-constraint use fail explicitly. Pinned the recovered
  scientific environment in a project-local virtual environment.

- Added an opt-in Gaussian volumetric-enthalpy profile with strict local C,H
  phase inversion, a dedicated conservative downstream inverse and actual
  non-Gaussian density buoyancy integration. The seven-case / 38-arc /
  17-profile rerun is complete. Same-41-sensor minimum-temperature MAE drops
  about 5%, but height and median-temperature errors worsen; not promoted.
- Added exact polar integration over the original finite square, checking
  all five fluxes instead of energy alone for the new profile. Retained the
  original failed numerical checkpoint and rejected it for subsequent runs.
  Added high-order downstream audits, thermodynamic-profile/quadrature merge
  guards and hashed reuse of frozen temperature reductions.

- Added opt-in explicit-species ambient/phase EOS consistency, including RH
  conversion and an exact ambient table node. No-H2 ambient gas no longer
  acquires a spurious approximately -1.20 K temperature offset in this mode.
  The 7-interface / 38-arc / 17-profile rerun retains the existing gates and
  defaults. Centre-height MAE changes from 0.06359 to 0.06159 m; concentration
  changes are small and vertical width worsens slightly, so no promotion.
- Added frozen-trajectory replay with temperature, five-flux and sensor-arc
  identity checks before comparing raw thermocouple observations. Workbook
  hashes and fixed sampling windows are retained without modifying sources.
  Field merging now rejects incompatible ambient closures as well as source
  physics/isomers, and newly generated field shards include input/code hashes.

- Corrected the measured-LH2 HEM research source's pipe-to-Gaussian energy
  ledger: carry incoming pipe kinetic energy and explicit condensed-phase
  enthalpy with a checked component reference. Independent upstream checks
  remove a 5.99–6.22% energy mismatch on four measured nozzles; this is not an
  experimental accuracy claim. Historical HEM field results remain archived.
  The completed 38-arc/17-profile rerun reduces HEM centre MAE from 0.0700
  to 0.0636 m, but concentration and the joint promotion gate do not improve.
- Added stable pure-H2 HEOS volume/caloric departure screening, deliberately
  separate from a mixture EOS. The 20 K source requires nonideal treatment;
  multiplying the full H2/air mixture volume by pure-H2 Z is not adopted.
- Preserved the crosswind tensor quadrature while merging exactly duplicate
  radial exponents by symmetry. Added source-physics/isomer guards to field
  shard merging, pointwise field output and downstream trajectory storage.

- Added a spin-consistent LH2 source/near-field selector.  The pure-para
  PRESLHY end member preserves normal hydrogen as the default and passes all
  seven conservative interfaces. It reduces the fast-bound centre-height MAE
  from 0.0700 to 0.0545 m, but worsens VG, FAC2 and vertical width, so it is
  retained as an explicit composition uncertainty rather than promoted.
- Added a measured-pipe LH2 flash/droplet source that conserves mass,
  pressure-thrust momentum and total energy while retaining separate liquid
  and vapour rates.  It reproduces the published ESREL/TNO micrometre
  atomisation scale and fixes `C_ds=10/15/20` as a non-fitted sensitivity.
- Implemented PRESLHY D3.1 equation 45 for the physical GASFLOW-MPI phase
  relaxation coefficient and inverse critical-diffusivity/diameter screens. The
  measured-source droplets would require effective diffusivity below 0.071%
  of the 20.4 K NBS H2 reference value, or droplets 37.6--52.9 times larger,
  to survive the collective source transit, so no empirical phase-delay
  coefficient is introduced.
- Added a phase-correct homogeneous-equilibrium evaporation-source bound and
  connected it to the seven-trial PRESLHY conservative interface/field path.
  It passes all interfaces and improves VG and vertical width, but worsens MG,
  FAC2 and centre-height MAE, so the validated default is unchanged.
- Made the downstream five-flux state inversion retry strict failed solves
  from neighbouring physical width/velocity states without relaxing its
  tolerance.  This recovers the formerly fragile trial-23 path while leaving
  established solutions unchanged; the complete non-slow suite passes.

- Added a conservative axisymmetric-to-JETPLU LH2 handoff with transferred
  phase thermodynamics, adaptive energy quadrature, explicit crosswind
  entrainment selection, and reproducible PRESLHY coupled validation. The
  local density-scaled shear transition repairs the diagnosed vertical-width
  loss, while the fully four-flux path now fails closed where JETPLU needs an
  independent thermal-profile state.
- Pre-registered the coefficient-free seven-state crosswind extension that
  separates centre density from centre H2 fraction and transports total energy
  alongside mass, species and vector momentum.
- Implemented its first fail-closed gate: independent density/composition
  Gaussian profiles and an adaptive-quadrature five-flux boundary solve. All
  seven frozen PRESLHY 10D interfaces pass, including the two former
  single-scalar failures.
- Added the downstream five-balance ODE using conservative-flux RK4 with a
  positive state inversion at every stage. The 0.02/0.01 m PRESLHY runs are
  converged and close long-range balances below `1.04e-7`; VG, FAC2 and width
  improve, but MG and centre-height MAE fail the frozen promotion rule.
- Added a common PRESLHY geometry observation operator that fits modelled
  direct-plus-ground-image values at the actual sensor heights. It confirms
  that the independent-energy rise/width rejection is not a hidden-state
  comparison artefact.
- Added Li et al. (2026)'s enthalpy-only established-flow balance as an
  explicit alternative to HyRAM+ total energy. It passes all seven interfaces
  and conserved marches but changes centre MAE by less than 0.04 mm, ruling
  out resolved kinetic-energy thermalisation as the rise-error mechanism.
- Added a direct near-field versus projected crosswind buoyancy-moment audit.
  All seven interfaces preserve the force sign and differ by at most 6.761%,
  rejecting an interface density-shape discontinuity as the rise-error cause.
- Corrected the Houf local-Froude width conversion from scalar to velocity
  e-folding width (`B=sqrt(2 sigma_y sigma_z)/lambda`). The repeated 7-trial
  run improves MG/VG to 1.062/1.163 and internal width ratio to 0.976, but
  centre-height MAE remains above baseline, so the path stays research-only.
- Added a trial-10 vertical-momentum budget. It finds only 0.0657 N handoff
  momentum difference but 5.695 versus 2.734 N cumulative buoyancy by 6 m,
  locating the remaining rise error in the early thermal/density evolution.
- Added a coefficient-free metastable dry-air condensation bound. Trial 10
  remains warm and weakly buoyant and fails the unchanged 2 K interface gate,
  so delayed nucleation is rejected before any seven-trial field scoring and
  pre-10D mass entrainment becomes the next isolated mechanism.
- Added a four-flux `source_flux` Gaussian-establishment bound after tracing
  the source-zone definitions in Li et al. (2026). It removes a second
  application of Zone-V momentum entrainment inside Zone IV, passes 7/7 and
  improves concentration to MG/VG/FAC2 1.012/1.162/0.976, but does not beat
  baseline centre-height and width errors and remains research-only.
- Added coefficient-free `geometry` and `surface_layer` ground-contact bounds
  to the independent-energy crosswind model. Trial 10 remains 0.469 m too high
  at 6 m and fails the frozen all-section gate, so both remain off-default and
  no seven-trial score is calculated.
- Added an explicit Li-equation-35 enthalpy transport option to the conserved
  near field. Its `source_flux` plug-to-Gaussian boundary is incompatible on
  the representative source and now fails closed before downstream scoring.
- Added direct PRESLHY thermocouple validation from the original trial-10 and
  trial-23 workbooks, including local sensor-height temperature observation.
  The 42-point audit finds the independent-energy candidate systematically
  too warm, especially in trial 23, and locates the next gap in phase/thermal
  development rather than drag.
- Corrected the D3.6 Table A3 reader for Unicode minus signs and seven-character
  thermocouple serials.  Added an explicit sustained-mean/peak source selector;
  the peak bound improves temperature and geometry but fails the frozen VG
  gate, so the mean-rate baseline remains unchanged.
- Added a synchronized raw-workbook source/temperature audit. Trial 23 shows
  the expected colder-at-higher-flow sign at all eight centreline channels,
  but the common contrast is too small to explain the cold-core residual;
  trial 10 is dominated by continued apparatus cooling. This prevents a
  transient source correction from being mistaken for a phase-model repair.
- Separated LH2 caloric temperature from mechanical tanker pressure and
  pre-registered a boiling-liquid source bound.  The temperature-only form
  fails the source applicability gate in all seven PRESLHY trials because a
  pressure-thrust-free reconstruction loses the available pressure work.  It
  remains off-default while measured nozzle pressure, quality/density and
  effective area are reconstructed.  Added a range-selective RADAR TAR
  utility so individual workbooks can be obtained without the 11.3 GB archive.

### Finite-rate condensed-air applicability is now quantified

The cryogenic-air add-on now exposes Sandia's uncorrected Ranz--Marshall
transfer number and a heat-limited minimum N2/O2 particle sublimation time.
The latter assigns all convective heat to latent heat and omits Stefan
resistance, so it is an explicit fastest-transfer bound rather than a fitted
kinetic model.  Across the eligible LH2 sources, 1 um N2 particles disappear
within 7--40 mm on this bound while 100 um particles usually survive past the
equilibrium handoff.  No finite-rate default is adopted because the field
campaigns do not constrain particle size, number density or nucleation.

## Unreleased

First working version.

### The conserved LH2 near field now has an audited crosswind handoff

A new research path carries the accepted axisymmetric near field into JETPLU
by projecting total mass, hydrogen and horizontal/vertical momentum as exact
integral constraints. Total energy, hydrogen half-width and centre temperature
are independent fail-closed screens. JETPLU now exposes its exact integral
fluxes, including analytic Gaussian kinetic-energy factors.

Two pre-registered reduced thermodynamic mappings fail unchanged. The accepted
candidate transfers the full near-field N2/O2/H2O equilibrium and component
enthalpy radial profile. On the frozen 0.08 m representative boundary it
closes native balances below `6e-16`, energy to 0.0585%, H2 width to 1.57% and
centre temperature to 0.94 K on a strict-grid confirmation.
`run_lh2_crosswind_research` assembles the
consistent local-wind/near-field/handoff path for horizontal wind-aligned
releases; rejected handoffs cannot integrate downstream. Independent
PRESLHY/Spadeadam rescoring remains outstanding.

A pre-registered Spadeadam 4/6 pilot using the existing 1 micrometre conserved
source reconstruction stops before sensor scoring: test 6 passes a
phase-manifold table-domain extension, while test 4 fails fixed width and
temperature screens and lies below the momentum-dominated velocity ratio.
The coupled runner now makes that applicability condition fail-closed.

### The Raman benchmark now follows the final journal source

The active Hecht--Panda observations now use the final 2019 article's printed
mass-decay and mass-width fits, 0.2771 and 0.07069, while preserving the 2017
conference values as versioned provenance. Model predictions are unchanged.
The recommended dry conserved model remains inside all four 25% bands, but its
mass-centre error is now -24.86% and its temperature-width error +24.56%.

The validation record also no longer silently resolves a source contradiction:
Table 1 and radial panels contain nine conditions, whereas the aggregate-fit
legends list an additional 4 bar, 45 K, 1.25 mm series. Pass counts are marked
provisional pending the authors' fit membership, reduced data and uncertainty.

### Para-hydrogen calorics and ELVHYS source identifiability are audited

The conserved axisymmetric model now accepts an explicit normal-, para- or
ortho-hydrogen component enthalpy table without changing the species mass or
density relation. The public LH2 research entry point exposes the same
off-default spin-isomer sensitivity. A pre-registered nine-case para-hydrogen
run closes every conservation threshold, but worsens the corrected Raman
thermal-centreline error from -21.38% to -34.79%; normal hydrogen remains the
recommended dry baseline.

Checksum-verified ELVHYS Test-10/Test-11 raw records were independently
reduced. A quantitative model score is deliberately withheld: the public
archive contains no hydrogen mass-flow channel, the nominal repeats have
different pressure/temperature source histories, and official documents give
conflicting 200 and 250 mm nozzle elevations. The exact missing fields and a
ready-to-send HSE request are documented.

A second pre-registered candidate combines independent thermal/species radial
profiles with phase equilibrium through a fully conservative four-flux
boundary. It reduces the Raman temperature-width error from +24.56% to -2.27%
but worsens centreline temperature to -29.14% and is retained only as a
research closure. The former unsafe phase/two-scalar combination remains
blocked unless four-flux establishment is selected.

### The accepted Raman model now has a public research entry point

Measured throat conditions can now be expanded through the same conserved
HyRAM+ boundary used in validation and passed directly to the recommended
dry-air axisymmetric model. The returned object carries source state,
trajectory, boundary residual, species drift, external-energy-corrected drift
and applicability warnings. Existing atmospheric and DEGADIS defaults are
unchanged.

Atmospheric argon phase change was added as a pre-registered explicit
sensitivity. It remains 4/4 and grid-converged but worsens the thermal-centre
error, so it is not enabled in the recommended configuration. A grey
radiative-absorption energy source and separate external-heat ledger were also
implemented. Even the perfect-black `5B` upper bound slightly worsens both
mass metrics, so its absorptivity stays zero by default.

### Humid-air frost is now included in the conserved Raman model

The axisymmetric research model can now carry explicit ambient absolute
humidity and solve local H2/N2/O2/H2O phase equilibrium, including water
freezing, vaporisation-plus-fusion latent heat and condensed-water volume.
It remains off by default. A pre-registered saturated-air upper bound closes
all boundary invariants and keeps nine-case energy drift below `1e-4`, but
overpredicts centreline warming and is rejected. A 40%-RH sensitivity puts
all four Hecht--Panda slopes inside 25%; the experiment does not report RH,
so this is a diagnosed input uncertainty rather than a fitted improvement.

The Raman comparison now follows the unequal stitched-image coverage in the
source paper (369 rather than 549 samples). Fixed turbulent Prandtl/Schmidt
spreading, temperature-dependent ideal enthalpy, dry-air-dew-point initial
heating and their two-scalar combination were also tested prospectively and
rejected without altering production defaults.

The reported 0.3 m/s honeycomb co-flow can also be represented explicitly.
Entrained co-flow mass, vector momentum and kinetic energy are conserved and
the entrainment speed is relative to the moving ambient. Its effect on all
four Raman slopes is below 0.8%, so it is retained as a boundary option but
rejected as an explanation of the thermal-centreline error.

The dry equilibrium phase model now optionally replaces constant ambient
heat capacities with low-density Helmholtz ideal-gas enthalpy tables for H2,
N2, O2 and H2O. This candidate closes boundary/species/energy residuals at
`1.50e-14`, `1.53e-5` and `3.77e-5`, and passes all four Raman slopes under
both the 549-point audit and corrected 369-point protocol. It is the first
end-to-end conservative Raman candidate to do so and is the recommended
dry-air axisymmetric research configuration; atmospheric production defaults
remain unchanged pending humidity-conditioned data.

### A conserved-energy cryogenic free-jet model is Raman-validated

An independent implementation of the published axisymmetric Gaussian
mass/momentum/species/energy balances now covers quiescent cold-gas near
fields. It includes plug-to-Gaussian flow establishment, source-momentum and
buoyancy entrainment, kinetic energy, and separate density/species profiles.
On all nine Hecht--Panda Raman releases the official-establishment variant
passes all four pre-registered centreline and half-width slopes; its aggregate
metrics agree with a separate official HyRAM 6.1 oracle to within 0.75%.
An end-to-end audit then exposed 8.9--14.8% species loss in that published
boundary. Two exactly conservative alternatives are retained as rejected
research variants because both fail the centreline-temperature criterion.

### Source total energy and pressure thrust are explicit research options

The transported condensed-air source can now include kinetic energy in the
storage-to-evaporation balance, removing an energy-creation defect at its
first station.  A separate HyRAM+ 6.0 path finds the homogeneous-equilibrium
critical throat, infers `Cd` from measured flow, conserves geometric-orifice
pressure thrust with the Yuceil--Otugen convention and solves the atmospheric
state from total enthalpy.  Unit tests close throat energy and notional-nozzle
mass, momentum and energy.

Neither option changes the corrected default.  The energy-only candidate
worsens PRESLHY variance/centre height and independent Spadeadam bias.  The
pressure-thrust candidate improves one PRESLHY centre-height subset but
worsens campaign-level bias and width, fails Spadeadam, and has no physical
ambient-pressure state for three of nine eligible PRESLHY releases.  Those
incompatibilities are reported rather than clipped.

The source can also preserve the Li et al. equation-22 evaporation-zone
distance instead of placing an already air-loaded state at the orifice.  The
12--78 mm coordinate correction changes no handoff state and slightly lowers
the predicted plume, but does not improve the two validation campaigns enough
for adoption.

Li et al.'s stationary-condensate limit is now available as a second
off-default research bound with separate gas/particle endpoint kinetic energy.
Deposited solid N2/O2 removes its enthalpy and zero axial momentum, with mass,
momentum and energy closed explicitly.  It is rejected: PRESLHY variance and
centre/width geometry and independent Spadeadam bias are all worse than the
corrected default.  Together with the fully carried bound, this prevents an
unsupported intermediate slip coefficient from being fitted to field errors.

### Condensed-air transport is available as a rejected research bound

The LH2 evaporation endpoint now permits stable solid N2 and O2 below their
triple points, closes component partial pressures and enthalpy, and feeds an
off-default source march that conserves retained/dropped mass, axial momentum
and energy.  The march retains H2 non-ideality at its component partial
pressure, transports re-evaporating condensed air, applies
Schiller--Naumann settling, and preserves its physical distance when handing
off to JetPlume.  Unit and integration tests cover pressure, energy,
continuity and the explicit source option.

Fixed 1, 10 and 100 um cases were registered before comparison.  None is
adopted: PRESLHY common-arc VG becomes 2.058, 2.150 and 8.946, centre-height
MAE becomes 0.285, 0.295 and 0.482 m, and independent Spadeadam MG/VG becomes
1.353/1.598, 1.357/1.600 and 1.460/1.694.  Fine retained condensate reproduces
the too-buoyant warm-source result; coarse dropout makes the remaining source
more H2-rich and worse.  The public corrected default is unchanged.

### The flashing source is now conserved through the first ODE station

The equivalent-source flash calculation already included the air needed to
evaporate the remaining hydrogen, but the first jet lookup placed that state
on a pure saturated-hydrogen mixing curve. For Spadeadam test 6 this changed
the hydrogen mass fraction from about 0.400 to 0.985 and discarded most of the
entrained air. The corrected path rebuilds the thermodynamic table from the
mixed flash composition and enthalpy, and a regression test now requires mass
fraction, temperature and density to survive the first lookup.

Gas-branch CoolProp evaluation is scoped to corrected jet sources. Applying it
globally damaged two independently checked NASA pool regimes. The public LH2
jet preset no longer forces permanent ground contact: with the conserved
source, test 4 remains low while low-wind test 6 detaches, matching the
reported Spadeadam regime split.

The expanded source also used to manufacture total mixture momentum while
stationary air was entrained: about 3.28 times the incoming value at 5 barg
and 6.20 times at 1 barg.  The corrected source now enforces
`u_out = Y_H2 u_in`, matching the no-slip initial-entrainment balance in the
current Sandia HyRAM jet model.  Pressure thrust from an under-expanded zone
and condensed-phase slip remain separate, explicit uncertainties.

Seven 0.35/0.53 m PRESLHY arcs lie upstream of the enlarged established source
plane. On the 62 common downstream arcs, the complete corrected model has MG
1.047, VG 1.425 and FAC2 0.84. On 23 filtered vertical fits it raises the mean
spread ratio from 0.64 to 0.731 and leaves mean / median signed centre-height
errors of 0.036 / 0.002 m. The earlier table-only source reached a spread ratio
of 0.933 but did so with non-conserved momentum, so that number is superseded.
Houf–Schefer buoyancy entrainment and DEGADIS surface-layer top entrainment are
available as pre-registered, off-by-default experiments because neither gave a
campaign-wide improvement without a competing regression.

The one-call report now cites the corrected jet evidence (MG 1.047, n=62)
instead of the historical baseline. Pool LFL distance no longer borrows the
PRESLHY jet statistic; it states that no direct pool-concentration validation
is available.

### The LH2 statistics are now computed

Reduced tables of the PRESLHY E3.5 campaign and the REDIPHEM archive make the
field results reproducible without either archive, and the suite recomputes
them rather than quoting them. `validation/nearfield.py` holds the liquid
hydrogen jet configuration, which had never been committed: every LH2 result
in this project had been produced by assembling the thermodynamics, the
boundary layer, the flash and the coefficients by hand in a session.

The configuration is checked against a full parameter dump of a reference run
and matches it to four significant figures. Getting there needed four fixes,
one of which is worth stating on its own: **`distmx` is the integration step,
not the limit.** Passing the 40 m limit there integrated the whole plume in a
single stride and produced a trajectory that was monotone, smooth, plausible
and wrong by a factor of two. The averaging time (60 s) and roughness
(0.001 m) were nowhere in the records and had to be recovered by matching
`deltay` and `ustar` against the dump; the source flow is the window mean, not
the peak.

**The published results reproduce.** Provenance exactly — 24 workbooks to 9
trials to 69 arcs. Concentration to two per cent: MG 0.722, CI [0.574, 0.903],
VG 1.44, FAC2 0.82 against a published 0.738, [0.591, 0.918], 1.41, 0.83. The
vertical spread ratio at the band it was quoted over, 0.64, exactly.

The vertical spread reproduces exactly, 0.64 -> 0.97 against a published
0.64 -> 0.94, once three unrecorded choices are pinned: the basis is 23
momentum-filtered fits rather than the 42 used for the trajectory, the
statistic is the mean of per-fit ratios rather than the ratio of medians, and
there is no distance filter.

**One published result is withdrawn.** The trajectory figures `1.07 -> 0.19 m`
come from configurations that cannot be reconstructed, and were taken over
different subsets besides. From the two dumped configurations the pair is
**0.877 -> 0.658 m** on trial 10 at 6 m, and 1.00 -> 0.85 m over the 5-7 m
band: the corrections take about a quarter off the rise, not four fifths. The
fault itself is confirmed and is the largest open defect in the LH2 path --
the measurement puts the plume below the nozzle where the model puts it a
metre above.

Every discrepancy in this exercise traced to an aggregation reported without
its subset, so every statistic now carries its distance band, release heights,
filter and n. `docs/lh2-recomputed.md` has all of it, and superseded figures
are marked in place rather than deleted.

### Review pass: the reported numbers and the computed ones

A review compared what the package prints and what its documents claim against
what the suite actually produces. The physics was unaffected; the reporting was
not.

**Two defects in `degali.lh2`.** The scope warning checked `max_distance` —
the integration limit, which defaults to 100 m — instead of the furthest
distance any reported answer relies on, so every call taking the default warned
that 100 m was out of range while every number it returned was inside it. And a
single `VALIDATED` range conflated two campaigns: the wind bounds came from the
NASA pool spills and the distance bound from their 33.8 m tower row, but both
were applied to jets whose corrected evidence spans 0.79 to 6 m. A jet asked about 30 m
was reported as inside the checked range when nothing had checked it. The range
is now per release type, and a "range" that is really one point says so.

**`degali.evidence`, new.** Every validated number, with its sample size,
grade and source, in one module. `lh2.py` printed a superseded statistic —
`MG 0.74, VG 1.22, n=53` — from a hard-coded f-string, while every document
carried the current `MG 0.738, VG 1.41, n=69`;
the neutral-buoyancy concentration appeared as four different values across the
README, the handover, a test assertion and the program's own output. Both are
now read from one place, and the report carries each claim's grade so that an
n=4 result cannot be read as firmly as an n=69 one.

**Guards, so this does not recur.** A retired figure reappearing in prose fails
the suite; the printed evidence is checked against the module; the scope-warning
regression has a test naming it; and the absolute-path guard now walks the
tests as well as the package — it had only ever walked `src/`, and the tests had
accumulated eleven machine paths.

**Test data is located from the environment.** `DEGALI_E35_ROOT`,
`DEGALI_E35_REPORT` and `DEGALI_SMEDIS_ROOT` join `REDIPHEM_ROOT`. The
reduced PRESLHY tables are derived products of this work rather than the
dataset, so they are vendored under `reference/preslhy/` with their provenance;
two far-field results now reproduce on a clone with no third-party data at all.

**Citation metadata.** `CITATION.cff` had no `authors` block, which CFF 1.2.0
requires and Zenodo rejects. Added, with `version`, and a test that blocks the
placeholder once the version is no longer a development one. The REDIPHEM entry
cited a title that is not the title of anything: REDIPHEM is the CEC DG XII
ENVIRONMENT project, and the database report is Nielsen and Ott, *A collection
of data from dense gas experiments*, Risø-R-845(EN), 1995.

**Documentation corrections.** The NASA lift-off table in `docs/liftoff.md` was
stale by up to 4 m against what the code produces, and quoted a mean error of
−0.1 m where the four values give −1.0. The suite regenerates it. The
handover's neutral-buoyancy figure of 85–90 mol % is superseded; the mixing table
gives 99.97, and the reason the two are so far apart — one linear segment at the
top of the table, in a regime where air is treated as a non-condensing ideal gas
— is now recorded alongside it. Test counts have been removed from prose
entirely: they were quoted in four documents and had drifted in all four.

### The model

All six DEGADIS 2.1 programs are ported: the source blanket (`DEG1`), the
steady and transient downwind models (`DEG2S`, `DEG2`), the cloud snapshots
and receptor histories (`DEG3`, `DEG4`), and the jet/plume model with its
touchdown bridge (`JETPLU`, `DEGBRIDG`). All five EPA test cases run end to
end in Python.

### Validation against the original

The Fortran is built from source in this repository and run as part of the
test suite. Purpose-written probes extract its internal state at full double
precision, because printed output carries five significant figures and cannot
distinguish a correct port from a close one. Every module below the
integration drivers matches to round-off.

Eight portability patches were needed to build the 2012 sources with gfortran.
One was load-bearing: `COMMON /ERROR/` mixes `REAL*8` with an `INTEGER*4`, and
`ESTRT1`'s `EQUIVALENCE` overlay made gfortran pad it differently there than
elsewhere, so the air entrainment coefficient was silently read as zero.

### Three things that trip up a port

- `ADIABAT(ifl=1)` reads its `wa` argument rather than computing it; `SZLOCAL`
  passes an uninitialised local, so `SZF` runs its lookups as though the
  mixture contained no dry air.
- `SURFAC` receives a mix of heated and adiabatic layer properties, and `PSS`
  and `SSG` differ by one letter in which temperature they pass.
- `GAMINC` returns the *unregularised* incomplete gamma, and says so in the
  comment above the assignment. That is deliberate, not a defect — but
  mapping it onto `scipy.special.gammainc`, which is regularised, leaves the
  flammable-mass derivative low by Gamma(1/(1+alpha)).

### Field evaluation

Against the Risoe REDIPHEM database, with bootstrapped confidence intervals.

At the lowest instrumented height Burro gives MG = 0.81 with a 95 % interval
of [0.63, 1.02], which straddles Hanna's bias bound rather than satisfying it.
Over the whole cloud nothing passes: DEGADIS gets the ground-level centreline
roughly right and the vertical distribution badly wrong, so the conventional
single-height statistic is carried by the one elevation where the vertical and
lateral errors cancel. The lateral spread runs one to three times too wide.

The reader refuses to guess channel type numbers, which are per series. An
earlier version fell back on a plausible list for FLADIS, which ships no
definitions, and produced a clean-looking result computed on channels that
were not concentrations. That result is withdrawn.
- Field screening now fails closed on excessive semi-FV mass-inventory
  residuals. The transport diagnostics remain in the report for audit, but a
  non-conservative obstacle/source solve cannot be promoted to an operational
  screening result.
