# Field evidence readiness review (2026-10-05)

## Purpose and rule

This review identifies what the connected local projects can support in the
DEGALI field-deployment path. It does not copy any operating or third-party
experimental raw data into this repository, and it does not convert an
ordinary operating transient into a release-dispersion validation case.

An observed pressure decrease or a vent-valve sequence is not, by itself, an
atmospheric release source. The source boundary, release geometry, phase
state, meteorology and receptor operator must remain separately evidenced.

## Reviewed local evidence

The read-only `field-audit` scan was rerun against the connected folders on
2026-10-07. The latest v14 `SLABx_LH2` audit scanned 16,139 structured files and was classified as
`partial`: 15 files expose an explicit `source_bound`/source-boundary field and
33 files expose candidate quantitative fixed-receptor rows through explicit
observed-volume-percentage columns, while wind/temperature cues without an
explicit weather identity remain diagnostic only. Reference-wind and mean-
concentration aliases improve triage but do not create an event join. The
audit's `diagnostic_counts` records these
categories by candidate file and
is recomputed from the immutable notes on read-back, so a stale summary cannot
silently change readiness. The v13 artifact also emits per-channel
`collection_requirements`: detected candidates still receive a manual
qualification action, while missing channels receive a concrete collection
action without enabling promotion.
The independent weather, `obstacle_geometry` and `common_clock` channels remain
missing, and the source-boundary candidates are model-result summaries rather
than a qualified site release boundary, so `promotion_allowed=false`. The fresh execution artifact
is `outputs/current-slabx-audit-2026-10-07-v14.json`.
With the bounded XLSX header/sample
reader enabled, `액화수소검사지원센터 데이터` exposed 20 structured CSV/XLSX
files (17 readable workbook candidates; 14 retain ambiguous auxiliary-sheet
diagnostics), but none of the five joinable validation channels; its status was
`withheld`. `액화수소센터 관련
자료` exposed no structured CSV/JSON/XLSX files and was also `withheld`. Fresh
support-centre and related-material execution artifacts are
`outputs/current-support-center-audit-2026-10-07-v14.json` and
`outputs/current-related-audit-2026-10-07-v14.json`.
These results do not inspect PDF/PPTX claims as quantitative observations and
do not create validation evidence.

### Static engineering documents checked separately

The connected document set was also checked as static engineering context,
without promoting diagram labels or scenario prose into event evidence. The
two-ton tank GA drawing provides vessel dimensions and support geometry, while
the tank P&ID and the A1-9106 instrument diagram show nozzle/vent topology,
instrument tags and the TK-1101/TK-1102 process arrangement. The scenario
description deck describes three operating scenarios and an automatic BOG vent
cycle from 6.5 bar to 6.0 bar. These documents are useful for preparing a
future coordinate/obstacle and source-boundary request, but they do not identify
one measured event, a non-negative atmospheric mass-rate history, simultaneous
weather, fixed-receptor H2 observations, calibration/operator metadata or a
common clock. They therefore remain context-only and cannot be selected into a
`FieldEvidenceManifest` without the missing event-linked exports.

| Local source | What is available | Appropriate use | Cannot support without more evidence |
|---|---|---|---|
| `SLABx_LH2/model-comparison` | FFI/PRESLHY-derived inputs, fixed-sensor comparison tooling, time-window and wind-history provenance, SLABx/DEGALI comparison ledgers | matched model-form studies; source/handoff and temporal-operator audits | validation of DEGALI's local obstacle/wake closure or transfer of a fitted FFI state to a plant |
| `SLABx_LH2` Test 4 material | Digitised 5 s anemometer history from the public report and declared sensor windows | transient-wind input/provenance workflow demonstration | raw high-frequency plant meteorology or a generic sensor response calibration |
| `액화수소검사지원센터 데이터` vent workbooks | Seven storage-vent sheets at mainly 1 s spacing; time, tank pressure channels and level channels are visible in the reviewed schema/quality inventory | historian sampling/continuity and operational-sequence plausibility checks | a direct vapour source history: no reviewed independent mass-flow, release-temperature, phase-fraction, release direction, local wind or downwind H2-concentration record |
| `액화수소검사지원센터 데이터` intake/supply workbooks | Several 1 min operational histories, pressure, temperature and level quality summaries | slow storage/transfer operating context; tag-quality review | short transient source reconstruction or sensor-arrival validation |
| `액화수소 저장탱크_밴트 시스템` parameter ledger | equipment dimensions, tank design pressure, documented nominal flow labels, operating sequences and source locators | scenario geometry and auditable engineering assumptions | measured discharge rate; one nominal normal-volume-flow label is not a mass flow until standard temperature, standard pressure, composition and meter/calibration basis are supplied |

The reviewed vent-data quality summary lists a 1 s baseline with occasional 2 s
gaps. It is useful evidence for sampling cadence, but a cadence is not sensor
response time or calibration uncertainty. The intake/supply histories are
mainly 60 s, so they cannot establish a seconds-scale release source boundary.

## DEGALI input readiness

The code-level promotion flag is now tied to `FieldValidationEvidence` rather
than a bare boolean. A strict case that sets
`lh2_validation_available: true` must fingerprint the external dataset and
declare a common clock plus source, weather, obstacle, receptor and temporal
operator identifiers. A free-field dataset remains conditional when an
obstacle is declared; this contract records evidence readiness but does not
create missing plant measurements.

| DEGALI field input | Current evidence disposition | Required promotion evidence |
|---|---|---|
| Upstream pressure | pressure tag time histories exist for some storage/vent events | pressure reference (absolute/gauge), tag calibration certificate, synchronised event window |
| Upstream temperature | some slower operating records include temperature-tag quality entries | exact upstream location, units/offset, response time and calibration uncertainty for the selected event |
| Atmospheric H2 mass flow | not established by the reviewed operational workbooks | calibrated flow meter or independently auditable mass balance; if normal volume is used, declared normal-state basis and composition |
| Liquid fraction / flash state | not established | measured/declared thermodynamic state and phase-fraction basis at the source boundary |
| Hole/vent geometry and direction | plant drawings can support a declared scenario after drawing review | release-point identifier, opening area/Cd evidence, orientation and elevation in the site coordinate system |
| Wind speed/direction/stability | not co-recorded in the reviewed plant event extracts | local met-mast record, sensor height, averaging period, clock alignment and stability/roughness evidence |
| H2 receptor concentrations | not present in the reviewed operating records | receptor locations/heights, H2 analyser calibration and t90, acquisition clock and averaging operator |
| Obstacle validation | no relevant LH2 obstacle concentration dataset identified | a controlled or field case with source, geometry, wind and H2 receptors all observed on a common clock |

## Safe integration path

1. Keep the existing `MeasuredHistoryQualityCriteria` gate. A selected
   historian interval may only become a `MeasuredReleaseHistory` after its
   pressure, temperature, mass-flow and optional liquid-fraction channels,
   units, response times, offsets and uncertainty bounds are explicitly
   supplied.
   `read_measured_history_csv()` is the controlled import boundary: it takes
   an explicitly mapped SI CSV plus a distinct phase/flash-state evidence ID
   and per-channel calibration/response evidence,
   rejects missing or nonphysical channels, retains the input SHA-256 and does
   not infer a release source from pressure-only Excel exports. The promoted
   field report retains the event ID, file fingerprint and each imported
   channel's unit, tag/calibration evidence and timing declaration.
   For long selected windows, the same route automatically builds one bounded
   LH2 saturation table over the stated temperature interval; the provenance
   records whether it was used. Small windows retain the direct reference
   path because table construction would not be a speed benefit.
   The same controlled map is available to `degali field-screen`; its exact
   JSON shape is documented in
   [field-measured-history-input.md](field-measured-history-input.md).
2. Use the present tank records first as *plausibility and data-readiness*
   evidence. They can establish candidate operating/vent windows and sampling
   gaps; they must not be labelled leak-rate calibration or dispersion
   validation.
3. Build one site scenario only after an engineering owner declares the vent
   location/orientation, opening/Cd boundary and normal-volume conversion
   basis. The resulting field report should remain conditional until the
   meteorology and sensor-receptor inputs above are also evidenced.
4. Keep FFI/SLABx comparisons as matched-basis model-form studies. Do not
   blend their public experimental sensors with the plant operating records or
   use them to assign an unvalidated obstacle-wake coefficient.

## Immediate collection request

For one named vent event, request a minimal, read-only export with:

- common-clock timestamps for upstream absolute pressure, upstream
  temperature, vent mass flow (or normal-volume flow with its exact standard
  basis), and any liquid/level measurement used to infer phase;
- tag units, calibration certificate/reference, accuracy interval, t90 (or
  response specification), and historian clock offset/quality statement;
- valve command/feedback, vent ID, opening diameter/Cd basis, elevation and
  global direction;
- local wind speed/direction at a declared height, plus stability/roughness
  evidence; and
- any H2 detector positions, calibration and averaging/response settings.

Without the last two groups, a release source can still be screened as a
declared engineering scenario, but it cannot validate field concentration or
obstacle behaviour.
