# Field evidence readiness audit

`degali field-audit <directory>` performs a read-only inventory of structured
CSV/JSON/XLSX files. It looks for five conservative channels needed before a site
LH₂ obstacle-transport validation case can be assembled:

- an identified source boundary with a physical release rate;
- weather/meteorology with an explicit wind quantity;
- obstacle geometry or an obstacle-geometry identifier;
- fixed-receptor H₂ concentration observations with coordinates; and
- a common clock or synchronization identifier.

For the common-clock triage channel, explicit non-blank values in
`common_clock_id`, `clock_id`, `synchronized_clock_id`,
`synchronization_id`, `sync_id`, `time_sync_id`, or `clock_sync_id` are
accepted as synchronization identifiers. A timestamp paired with non-blank
`synchronization`, `synchronization_method`, `time_alignment`,
`common_time`, or `clock_offset_s` metadata is also retained as a clock
candidate. Blank values remain non-evidence, and these aliases only improve
inventory triage; they do not qualify an event join or relax promotion gates.

The command reports `withheld`, `partial`, or `candidate_complete`. The last
status means only that all channels appeared somewhere in the directory; it
does not prove that the files describe the same event. The report always sets
`promotion_allowed` to `false` and never creates a `FieldValidationEvidence`
record. Manual reconciliation of source identity, weather interval, obstacle
representation, sensor calibration, and clock is still required before
`field-validate` can be used.

That reconciliation can be captured as a separate, hash-pinned
`degali.field-evidence-manifest.v1` artifact with
`FieldEvidenceManifest.from_audit()`. The manifest requires one explicitly
selected candidate per channel and declares the event, IDs and temporal
operator; it still sets `promotion_allowed=false` and never bypasses the
matched-score gate. See [field-evidence-manifest.md](field-evidence-manifest.md).

Sensor registry/calibration data is intentionally not inferred as a sixth audit
channel. Provide it explicitly to manifest creation with `sensor_registry_path`
or `--sensor-registry-path`; the manifest then pins that file separately with
`sensor_set_id` and a SHA-256. Omitting the registry or accountable
`operator_id` leaves the package `conditional`; a registry with only nominal
specifications also remains `conditional` until
`sensor_calibration_status=certified` is declared. Missing audited channels
remain `partial`/`withheld`.

For automation, the report also emits deterministic `gate_codes`:
`channels_missing`, `channels_partial`, or `scan_incomplete` identify why
readiness is incomplete; `coherent_candidate_unqualified` and
`cross_file_candidate_unqualified` distinguish the two candidate-complete
shapes; and `promotion_not_allowed` is always present because inventory alone
never promotes validation evidence. The codes supplement, rather than replace,
the explanatory `reasons`.

The report also emits `coherent_candidate_paths`: files in which all five
channels were found together under the bounded record checks. An empty list
means the channels are necessarily split across candidates. Even a non-empty
list is only a same-file readiness signal; it does not authorize validation or
model promotion.

Each artifact also emits `diagnostic_counts`, a deterministic count of
candidate files by near-miss category (for example, release-rate fields without
source identity, wind fields without weather identity, operating pressure/
temperature/level or PT/TT/LT historian data without a source boundary,
missing concentration, missing clock, unreadable files, and bounded sampling).
A file is counted at
most once per category, so repeated rows cannot inflate the readiness summary.
When a saved artifact is read, these counts are recomputed from the immutable
candidate notes and a stale or tampered summary is rejected.

The artifact also emits `collection_requirements` for every required channel.
Each entry records whether a candidate was `missing` or merely `detected`, the
candidate paths, the accepted normalized field-alias groups, the remaining
qualification requirement, and a next collection action when the channel is
missing. `detected` never means qualified: the action changes to the manual
event/clock/geometry/calibration reconciliation still required before a
manifest can be selected. The requirements are derived from `channel_paths`
and are checked on read-back; a tampered or stale action plan is rejected.
Older v1 artifacts without this optional field remain readable.

The audit record itself is fail-safe: candidate fields, channel paths and
reasons are immutable, duplicate channel paths are rejected, and contradictory
states such as a `partial` audit with no detected channel are refused before
JSON serialization. On read-back, each channel path and coherent-candidate
path is reconstructed from the candidate `detected_channels`, so an index-only
artifact mutation is rejected even when the referenced files and SHA-256
digests are unchanged.

Header-only CSVs and empty JSON row arrays remain diagnostic candidates but do
not populate any readiness channel. Field names alone are not treated as a
measured source, weather record, or fixed-receptor observation; a non-empty
data record is required before the channel can contribute to
`candidate_complete`.
CSV audit inputs also require a non-empty, unique normalized header and exactly
the declared number of fields on every non-blank row. Duplicate/empty headers
and extra or truncated row values are retained as unreadable diagnostic
candidates rather than allowing `DictReader`-style silent column overwrite or
field loss to create a readiness channel.
XLSX sampling applies the same empty/normalized-duplicate header rule before
collecting bounded worksheet values. A malformed auxiliary sheet is skipped
with a sheet-specific diagnostic while valid sheets remain inspectable; a
workbook whose sampled sheets are all malformed remains an unreadable
diagnostic candidate. No channel is created from a skipped sheet, and the
XLSX sample is still not a full dataset validation.
For large CSVs, the audit validates every row's field count and retains the
complete `row_count`, but samples at most the first 200 data rows for value and
co-occurrence checks. This keeps readiness triage bounded without treating a
partial sample as a complete validation dataset.
JSON object keys are likewise required to be unique after normalization;
duplicate or normalization-colliding keys are retained as unreadable
diagnostic candidates instead of allowing the parser's last-value-wins rule to
alter a source, weather, geometry or observation field silently.
Each candidate also carries the lowercase SHA-256 of the bytes that were
scanned. The digest is provenance only: changing a file requires a fresh audit,
and it never promotes the candidate into `FieldValidationEvidence`.
The Python API `write_field_evidence_audit_json()` and
`read_field_evidence_audit_json()` use the same execution schema as the CLI,
refuse overwrite, and recheck every candidate file digest before accepting a
saved report. The saved JSON is created exclusively after candidate
fingerprints are rechecked, so a competing writer cannot replace the audit
artifact after preflight. Strict artifact writing also refuses an unpinned candidate with
no digest. A stale report therefore fails closed instead of silently being
reused after an input change. Candidate paths are required to resolve inside
the declared audit root; path-escaping artifacts are rejected before any file
is opened.
The execution record also preserves the `max_files` scan limit used to produce
the status, so a partial scan cannot be mistaken for a complete scan made with
a different bound.
Saved artifacts can be checked from automation with
`degali field-audit-verify audit.json`; the command rechecks candidate files and
returns exit code 2 under `--require-complete` when the stored status is not
`candidate_complete`.
The bounded XLSX sample applies the same rule after identifying its first
non-empty row as the header: a workbook with no later non-empty row remains a
diagnostic candidate, while its row count stays deliberately unknown.
Rows or JSON records whose values are all blank are treated the same way;
their presence alone does not count as an observation.
Common missing-value tokens such as `N/A`, `unknown`, `unspecified`, and
`not-recorded` are likewise not declared identifiers or measurements.
Channel fields must also co-occur in one bounded record: a source ID from one
row cannot be paired with a rate from another, and a receptor coordinate from
one row cannot be combined with a concentration from another.
For JSON event objects, top-level source/weather/geometry metadata is joined
only with the contained observation records from that same object; unrelated
files or separate row groups are not synthesized into one event.

If the structured-file limit is reached, the report sets `scan_complete` to
`false` and withholds `candidate_complete`; it does not silently treat a
partially scanned directory as complete.

Generated reports are new files and are never overwritten:

```console
degali field-audit C:\\data\\lh2-trial --output audit.json
```

Published pool-radius tables, CFD scenario inventories, and model-result files
are intentionally not treated as fixed-receptor H₂ observations merely because
their names mention a source, release, or concentration.

The separate neutral-gas Case-H and paired dense-gas fence readers apply the
same fail-safe principle at the row boundary: finite coordinates are required,
missing published sentinels remain missing, and negative or infinite
concentrations, RMS values, and turbulent kinetic energy are rejected. A clean
row therefore means only that the supplied observation is internally usable;
it does not change the benchmark's explicit non-LH₂ validation scope or allow
an obstacle-wake coefficient to be fitted automatically.

Sensor-summary files that contain only coordinates or threshold crossing times
are retained as candidate notes, but are explicitly labelled as non-quantitative
and do not populate the `receptor_observations` channel.
Explicit observed concentration columns, including volume-percent aliases such
as `observed_vol_pct`, `observed_mean_vol_pct` and
`observed_time_mean_vol_pct`, do populate that channel when sensor coordinates
and non-blank numeric values co-occur in the same record. This still does not
provide source, weather, obstacle or common-clock identity and therefore cannot
promote the audit.
Likewise, a release-rate or wind/temperature record without an explicit
source/weather identity is retained as a diagnostic note rather than treated
as a joinable boundary condition. Common source-rate aliases
(`source_rate_kg_s`, `mass_rate_reported_kg_s`, `rate_kg_s`) and wind aliases
(`wind_m_s`, `mean_wind_speed_m_s`, `mean_wind_from_deg`) are included in this
near-miss diagnostic, but do not relax the identity requirement. For legacy
collection artifacts, explicit `source_bound`/`source_boundary` fields,
reference-wind fields (`wind_ref_m_s`, `u_ref_m_s`, `wind_10m_m_s`) and
mean-concentration fields (`mean_c_pct`, `observed_mean_c_pct`, `mean_vol_pct`)
are also recognized as candidate aliases. They remain inventory signals only:
they do not create a cross-file event join or qualify a source, weather record,
or receptor observation for promotion. XLSX files are sampled from bounded
worksheet headers and
early rows; their `row_count` is deliberately left unknown and no workbook
identity, sheet relationship, or calibration is promoted automatically. If
optional `openpyxl` is unavailable, the workbook remains an unreadable note
rather than blocking the audit or creating evidence.
Timestamped operating pressure/temperature/level fields, including common
PT/TT/LT instrument-tag patterns, receive an explicit non-source diagnostic
when no atmospheric source-rate boundary is present. A pressure trend or valve
sequence therefore cannot become a leak rate by filename, tag name, or
co-occurrence alone.
