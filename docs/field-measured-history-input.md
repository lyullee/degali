# Quality-gated measured-history field input

`degali field-screen` can attach one selected SI historian export to a field
case through the optional `measured_history` object. This is a source-boundary
route, not Excel tag discovery: the CSV must already have a numeric,
strictly-increasing relative-seconds column and declared SI P/T/mass-flow
channels. A pressure drop, a valve signal or a nominal normal-volume flow is
not accepted as a substitute.

Add the following object to a `degali.field-screening-input.v1` case file.
`csv_path` is resolved relative to the case JSON file.

```json
{
  "measured_history": {
    "csv_path": "selected-vent-event.csv",
    "map": {
      "time_s_column": "time_s",
      "event_id": "vent-event-001",
      "event_evidence_id": "event-window-review-001",
      "phase_evidence_id": "source-phase-review-001",
      "pressure_pa": {
        "value_column": "pressure_pa",
        "unit": "Pa",
        "source_id": "PT-1105A",
        "calibration_evidence_id": "PT-1105A-calibration-001",
        "response_time_s": 0.2,
        "time_offset_s": 0.0,
        "absolute_half_width": 500.0
      },
      "temperature_k": {
        "value_column": "temperature_k",
        "unit": "K",
        "source_id": "TT-1106",
        "calibration_evidence_id": "TT-1106-calibration-001",
        "absolute_half_width": 0.25
      },
      "mass_flow_kg_s": {
        "value_column": "mass_flow_kg_s",
        "unit": "kg/s",
        "source_id": "FT-1101",
        "calibration_evidence_id": "FT-1101-calibration-001",
        "relative_half_width": 0.05
      },
      "liquid_fraction": {
        "value_column": "liquid_fraction",
        "unit": "1",
        "source_id": "LT-1101",
        "calibration_evidence_id": "LT-1101-phase-basis-001",
        "absolute_half_width": 0.0
      }
    },
    "quality_criteria": {
      "maximum_sample_interval_s": 1.0,
      "maximum_response_time_s": 1.0,
      "maximum_absolute_time_offset_s": 0.1,
      "maximum_relative_half_width": 0.1,
      "evidence_id": "historian-quality-procedure-001"
    }
  }
}
```

Each channel must use either paired `lower_column`/`upper_column` fields or
declared `absolute_half_width` and/or `relative_half_width`. The phase evidence
ID is required even if a liquid-fraction channel exists. The main case's
transient `source.duration_s` must equal the selected CSV duration. An
unresolved bounded `duration_uncertainty` cannot be attached to this schedule
because the historian fixes one time axis; use a separate time-aligned source
envelope instead of silently reusing the nominal schedule.

For paired calibration columns, every row must satisfy
`lower <= nominal <= upper` (and `lower <= upper`). Rows that reverse the
bound order or exclude the nominal reading are rejected before any flash or
uncertainty envelope is built; the importer never repairs them by sorting or
clipping.

The historian reader rejects empty or duplicate headers and any non-blank row
whose field count differs from the header. It does not allow `DictReader` to
overwrite a channel or silently discard an extra value before the pressure,
temperature, flow and phase quality gates run.

The resulting report retains the CSV SHA-256, event identifier, phase evidence
and channel calibration/timing declarations. A nominal history with non-exact
P/T/flow/phase bounds remains **withheld** for operational screening after its
normal calculation. Use the same CLI with `--uncertainty-envelope` to
propagate every global history selection with weather and every in-plane
detector calibration corner, refine every resulting case, and aggregate the
fail-safe decision:

```bash
degali field-screen measured-field-case.json --uncertainty-envelope \
  --allow-conditional --output measured-field-envelope.json
```

The programmatic counterpart is
`run_field_operational_joint_measured_history_envelope()`. This deterministic
envelope also propagates declared source location/direction and ambient
temperature/pressure bounds through
each history re-flash and ambient air-density bounds through detector
conversion. It is not a confidence interval or an obstacle-validation result.
Its operational report includes `deterministic_sensor_envelope` with
per-detector peak/final/time-average extrema from completed history/field
corners; missing history traces remain withheld with their reasons.

For a source-history-only sensitivity screen, the lower-level
`run_field_measured_history_envelope()` likewise crosses declared source
location/direction bounds before wind-plane projection, while leaving weather,
sensor and obstacle uncertainty to their dedicated field envelope helpers.

## Pressure-driven orifice history

When a historian has only an absolute upstream pressure and temperature trace,
do not convert the pressure trend into a leak rate.  The explicit alternative is
`PressureDrivenMeasuredHistory` together with a `ReleaseSource` whose
`opening_area_m2`, `discharge_coefficient`, and source identity are declared
and whose mass-flow boundary is a zero placeholder.  A single direct schedule
uses their nominal values; the envelope path crosses their lower/upper bounds
as source-boundary inputs.  The source helper
`direct_vapour_schedule_from_pressure_driven_history()` recomputes each interval
through the pressure-driven throat closure, flashes it, and retains the
post-flash liquid as an unrouted ledger term.  Its envelope counterpart
`direct_vapour_schedule_envelope_from_pressure_driven_history()` couples global
pressure/temperature/liquid-fraction and orifice-bound choices before
recomputing flow; it does not form independent per-sample rate crossings.

The same boundary can be selected from a strict case JSON with
`pressure_driven_history`.  Its `map` has the same `time_s_column`, event,
phase, source/calibration/timing and bound declarations as `measured_history`,
but contains only `pressure_pa`, `temperature_k`, and optional
`liquid_fraction`; it deliberately has no `mass_flow_kg_s` key.  The case
`scenario.source.mass_flow_kg_s` must be an exact zero placeholder. Opening
area and discharge coefficient may carry declared bounds; a nominal parse is
diagnostic and remains unresolved until the joint envelope crosses those
source-boundary values.
The imported CSV SHA-256 and mapping are carried into every field report.

For a fixed ambient boundary, use `run_field_pressure_driven_history_envelope()`.
It requires the same approved `MeasuredHistoryQualityCriteria` and carries the
derived source schedule into every semi-FV case.  When ambient
temperature/pressure uncertainty must be propagated with weather and detector
corners, use `run_field_joint_pressure_driven_history_envelope()`; each ambient
corner is passed into the throat closure before flash and the field request fixes
that corner so no uncertainty is silently dropped. Opening-area and
discharge-coefficient bounds are crossed in the same source envelope, while
declared source location/direction bounds are fixed before field wind-plane
projection. The raw historian file
still needs its own event, phase, clock, calibration and SHA-256 evidence; the
typed programmatic object is not a substitute for that evidence.

For refinement and the aggregate operational gate, wrap the same source with
`run_field_operational_joint_pressure_driven_history_envelope()`.  It requires
every joint corner to retain its own refinement, applicability status and
fail-safe decision; a single withheld corner withholds the aggregate result.

`degali field-verify` reimports the case CSV and compares its current SHA-256,
event identity and complete measured-history provenance against every nested
field report. For nominal/refined and joint-history executions it also reruns
the deterministic field calculation with the recorded options and compares
the complete report and operational decision. Editing the CSV or report
without regenerating the execution now fails closed; the case JSON fingerprint
alone is not treated as sufficient historian provenance. Phase-routing and
other envelope modes that do not serialize enough replay options remain
integrity-only.
