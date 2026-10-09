# Matched LH2 field-evidence request

This request is the minimum export needed to close the IJHE matched-transient
evidence warning. It is intentionally written as a controlled-data request:
the custodian may provide derived CSV/JSON exports and SHA-256 digests without
redistributing proprietary raw historian files.

## Why the current inspection-centre workbooks are not sufficient

The transferred inspection-centre workspace contains process histories (PT,
TT and DPT channels, mostly at one-minute resolution) and operating-sequence
metadata. The repository field-evidence audit classified it as `withheld`:
it does not provide an atmospheric source boundary, event-linked wind,
fixed-receptor H₂ concentration, obstacle geometry or a common-clock ID.
Those process histories can be retained as a `source_history` supplement, but
they must not be relabelled as a matched dispersion event. The detailed
read-only inventory is recorded in
[`docs/ijhe-center-operational-data-audit-2026-10-09.md`](ijhe-center-operational-data-audit-2026-10-09.md).

The public PRESLHY E3.5 Trial 10 workbook has now been reconciled into a
conditional five-channel bundle for replay and manifest testing. It supplies
the source, weather, receptor and common-clock fields, but it does not supply a
channel-level calibration certificate or a measured three-dimensional site
survey. It therefore does not yet close the request below or justify an
operational-accuracy claim.

## Required event package

Provide one event or trial with a stable `event_id` and a single declared time
axis. Every file must include its units, timezone and source-system identity.

| Channel | Required fields | Acceptable evidence |
|---|---|---|
| Source boundary | `source_id`, release start/end, location and orientation, orifice/geometry, phase, mass-flow history in kg/s; if pressure-based, the approved conversion and uncertainty | Historian export, calibrated flowmeter export, or signed source-term calculation tied to the event |
| Weather | `weather_id`, timestamp, wind speed, wind direction, measurement height, stability/temperature and station location | Anemometer/met mast export on the event clock; include averaging interval |
| Obstacles | `obstacle_geometry_id`, coordinate reference, building/wall/barrier vertices or dimensions, heights and vertical datum | Survey/CAD/GIS export or a signed geometry table |
| Receptors and sensors | `sensor_id`, x/y/z, H₂ concentration time series, units, sample interval, alarm/threshold definition | Raw or controlled derived sensor export; include missing-value and saturation codes |
| Sensor metrology | calibration date and gas, calibration curve or certificate, accuracy/bias, resolution, response time, sampling/display delay and detection limit | Calibration certificate plus the sensor registry used during the event |
| Common clock | `common_clock_id`, timezone/UTC offset, synchronization method, start/end coverage, offset and drift estimate for each source | NTP/PTP/GPS record, historian clock record or an explicit signed offset table |
| Accountability | `operator_id`, data custodian, processing version, file names and SHA-256 digests | A manifest reviewed by the accountable operator |

## Submission format

```text
event-<id>/
  source_boundary.csv
  weather.csv
  obstacle_geometry.json
  receptor_observations.csv
  sensor_registry.csv
  common_clock.json
  FieldEvidenceManifest.json
```

The manifest must identify the five promotion channels (`source_boundary`,
`weather`, `obstacle_geometry`, `receptor_observations`, `common_clock`) and
the sensor registry. Missing channels remain `conditional` or `withheld`;
the model must not infer wind, geometry, calibration or a time alignment from
nearby files.

## Minimum path to an IJHE revision

One complete event with calibration and site-geometry evidence is enough to
replace the current conditional warning with a qualified case study. Three or
more events spanning low, nominal and high wind would support a stronger
operational claim. Until then, the manuscript should keep the current wording:
numerical transient results are verification only, and industrial use is
limited to evidence-bounded screening and routing of cases to SLABx/CFD.
