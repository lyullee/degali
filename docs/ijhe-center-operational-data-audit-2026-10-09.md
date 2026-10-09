# Inspection-centre operational-data audit — 2026-10-09

This is a read-only inventory of the local `액화수소검사지원센터 데이터`
workspace. It was checked as a possible source for the IJHE matched-transient
validation gate. No source workbook or raw historian export is copied into the
public DEGALI snapshot.

## Inventory

The workspace contains **17 XLSX workbooks and 3 CSV operating-sequence
tables**. The files fall into three evidence groups:

| Group | Observed content | Evidence role |
|---|---|---|
| Predicted/stepwise scenario workbooks | Ten-run sheets with `from_step`, `to_step`, `progress_pct`, `action`, PT, TT and DPT tags | Operating-state and source-sequence context |
| Recovered one-minute process histories | Date-labelled sheets with time, pressure (PT), temperature (TT), differential-pressure (DPT) and video-display columns | Possible source-state history after custodian review |
| Operating-sequence CSVs | Intake, storage/vent and supply actions with qualitative leak/fire/explosion classifications | Procedure and hazard-context metadata |

The channels are process instrumentation. They are not atmospheric receptor
measurements. The workbooks do not provide a spatial H₂ concentration field.

## Matched-event gate

| Required channel | Found in this workspace? | Decision |
|---|---:|---|
| Event ID and accountable operator | No stable cross-file event manifest | `withheld` |
| H₂ source boundary or approved mass-flow history | Pressure, temperature and DPT are present; no measured H₂ mass-flow channel or approved conversion was found | `withheld` |
| Wind speed, direction, height and stability | Not found | `withheld` |
| Fixed-receptor H₂ concentration time series | Not found | `withheld` |
| Sensor coordinates/heights and calibration certificate | Not found | `withheld` |
| Obstacle geometry and coordinate reference | Not found | `withheld` |
| Common clock across source, weather and receptors | Not established | `withheld` |

Consequently these files **cannot close the IJHE matched-transient warning**.
They must not be relabelled as a dispersion-validation event, used to produce
MG/VG/FAC2, or used to claim obstacle-wake validation. Pressure histories may
be retained as a source-state diagnostic only after a custodian supplies an
approved pressure-to-flow method and uncertainty.

## What would make the data usable

The custodian would need to provide a controlled derived package with:

1. a stable `event_id` and operator/custodian identity;
2. a calibrated source-rate history or an approved pressure-based conversion;
3. event-clock weather at a declared height, including wind direction and
   stability;
4. H₂ concentration time series at fixed, surveyed sensor coordinates;
5. sensor calibration certificates, response time and delay information;
6. surveyed obstacle geometry in a declared coordinate system; and
7. SHA-256 digests and clock-offset/drift records for every selected file.

Until those items are supplied, the appropriate use is **operational source
context and data-intake testing**, not quantitative model validation. This
classification is consistent with the controlled request in
[`docs/ijhe-field-evidence-request-2026-10-09.md`](ijhe-field-evidence-request-2026-10-09.md).
