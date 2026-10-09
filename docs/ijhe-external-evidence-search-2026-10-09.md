# External matched-event search — 2026-10-09

This record documents the final local search for a time-resolved field event
that could close the IJHE operational-validation warning.

## Sources checked

- `reference/preslhy/raw/`: seven transferred E3.5 workbooks used by the
  common-window audit.
- `reference/hydrogen/open_channel_usn/extracted/`: 22 extracted open-channel
  files used by the USN timing audit.
- `reference/obstacle_wake/`: available obstacle/wake reference folders.
- FFI readiness outputs, including
  `outputs/ffi-test4-time-history-readiness-2026-10-08-v2.json`.
- The FFI Report 20/03101 site/instrumentation pages, which document the
  outdoor-pad containers, drum, instrument box and wind/concentration
  instrumentation; this is a geometry boundary, not a matched Trial 10 event.
- Existing field-audit manifests and legacy model-audit artifacts under the
  local `outputs/` directories.
- The transferred `액화수소검사지원센터 데이터` workspace, which contains 17
  Excel workbooks and three CSV operating-sequence files. The workbooks carry
  minute-scale process histories (PT/TT/DPT channels) and sequence metadata;
  they are operational/process histories, not atmospheric receptor data.
- The separate local `validation-data/smedis/` evidence directory, which
  contains 28 SMEDIS workbooks and reduced condition/receptor tables. These
  are generic dense-gas trials (ammonia, propane, R-12, SO2 and LNG/SF6/CO2
  records), not LH2 events. The EEC control/fence pair was also reviewed as a
  possible obstacle observation lane.
- The public [PRESLHY WP3 rainout dissemination record](https://hysafe.info/wp-content/uploads/sites/3/2021/05/3_4_Coldrick_Rainout_HSE_v2.pdf), which lists nominal
  ranges and accuracy specifications for Xensor, Dräger and Gill instruments.
  These specifications are useful metadata, but are not event-specific
  calibration certificates for the transferred Trial 10 channels.
- The National Laboratory of the Rockies H2SAFE data catalogue record
  ([DOI 10.7799/17118570](https://doi.org/10.7799/17118570)), which describes
  five indoor helium-release experiments with time-resolved concentrations,
  three-dimensional sensor/source coordinates and laboratory geometry. It is
  a useful sensor-operator and manifest-ingestion boundary case, not LH2
  validation evidence.
- The published [NREL HyWAM LH2 release description](https://www.aiche.org/sites/default/files/chs-member-files/manuscript_hydrogen_wide_area_monitoring_of_lh2_buttner_et_al2.pdf), which reports a measured
  three-dimensional sampling-frame layout and near-field LH2 monitoring. The
  publication provides geometry and instrument context, but not a complete,
  independently hashable Trial 10 time-series/calibration package suitable
  for promotion here.
- The official [KITopen/RADAR PRESLHY E3.5 research-data record](https://publikationen.bibliothek.kit.edu/1000136281),
  DOI `10.5445/IR/1000136281` (the archived RADAR package is versioned as
  `10.35097/1481`).  Its public metadata describes 25 elevated LH2 releases,
  H2 concentration and temperature channels, pipe/tank pressure, mass flow,
  near- and far-field weather channels, humidity, O2 depletion and video.
  The package is openly licensed under CC BY-SA 4.0, but its published
  archive is approximately 11.3 GB.  The associated public D3.6 report is
  also available as a selective report artifact and supplies nominal
  instrument accuracies, the shared logging architecture, receptor
  coordinates/heights, weather-station heights and the two obstruction-test
  offsets.  It does **not** supply channel-level calibration certificates,
  absolute-UTC offsets or a surveyed three-dimensional site/obstacle model.
  It is therefore a high-value custodian lead and a conditional
  geometry/uncertainty source, not silently promoted evidence.

## Result

The inspection-centre and generic SMEDIS folders do not satisfy the required
promotion channels.  A separate reconciliation of the public PRESLHY E3.5
Trial 10 workbook and D3.6 coordinate report now provides a **conditional
matched event package**.  The derived bundle is documented in
[`docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md`](ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md).
Its audit returns `candidate_complete`, and its `FieldEvidenceManifest`
returns `evidence_readiness=conditional` with `promotion_allowed=false`: the
registry contains nominal response metadata but no event-level calibration
certificate.

The conditional events were reconciled against the following five required
promotion channels:

1. event-linked source history;
2. weather speed, direction and stability on the same clock;
3. receptor coordinates and H₂ concentration time series;
4. sensor calibration/response metadata; and
5. an explicit common-clock identifier with offsets and covered channels.

The source, receptor and common-clock channels are represented in the derived
Trial 10, Trial 20 and Trial 21 bundles. The weather channel contains speed and
direction records, but stability and measurement-height metadata are not
provided in the workbook; it is therefore only a partial match. The D3.6
report provides a nominal specification table (including Xensor XEN-5320
`±3% of full scale`, Dräger Xam 5000 H2 `±2% of reading`, and an estimated
cryogenic Coriolis flow accuracy of `±3% of reading`), but these are device
specifications rather than event-level calibration certificates. The bundles
record the sensor registry and response interval without inventing a
calibration error. The report also provides source-relative receptor
coordinates/heights and documents obstruction trials at 0.16 m and 0.18 m
from the release; it does not provide a complete surveyed 3-D obstruction
mesh, so the promotion lane remains an open-pad/conditional geometry lane.

The additional replay was performed from locally controlled public-archive
copies using the same extraction command and D3.6 coordinate report. It yielded
two further verified manifests:

| Event | Receptor rows | Manifest status | Promotion |
|---|---:|---|---|
| Trial 10 | 5,053 | `conditional` | `false` |
| Trial 20 | 19,344 | `conditional` | `false` |
| Trial 21 | 2,232 | `conditional` | `false` |

The row counts are derived observation records, not independent validation
trials. The replay demonstrates that the evidence workflow is reproducible
across multiple accepted releases; it does not remove the calibration,
absolute-UTC and surveyed-geometry gaps.

Separately, the FFI report confirms that its outdoor pad contains known
structures. The six-arc FFI lane therefore remains a conditional transport
screen because DEGALI does not resolve that three-dimensional geometry; the
source and hash record are in
[`docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md`](ijhe-ffi-site-geometry-boundary-2026-10-09.md).

The newly checked inspection-centre workspace is an explicit partial-data case,
not a hidden matched event. Running the repository field-evidence audit over
that folder returns `status=withheld`; all five promotion channels are missing,
with `instrument_history_without_source_boundary=1` and `sampling_limited=17`.
The files do not expose an event-linked atmospheric mass-rate boundary, wind
speed/direction/stability, receptor coordinates with H₂ concentration, obstacle
geometry, sensor calibration/response metadata, or a common-clock ID. Their
timestamps are useful for a future process/source-history handoff only after an
operator supplies the missing release and atmospheric joins.

The remaining E3.5 workbooks still provide source and sensor histories for the
declared 20-second sensitivity lane, but their weather records have not been
individually reconciled into manifests. The FFI Test 4 readiness output retains
`common_clock_established=false`. The open-channel H₂ audit reports
`weather_geometry_common_clock_missing`. Obstacle folders do not add a matched
LH₂ concentration event. The SMEDIS files can support a separate generic-gas
transfer benchmark, but their non-LH2 thermodynamics, static reductions and
missing channel-level clock/calibration metadata prevent promotion into the
LH2 manuscript score. The machine-readable run is
`outputs/smedis-transfer-field-audit-2026-10-09.json`; it returns
`status=withheld`, `promotion_allowed=false` and all five promotion channels
as `missing`.

The additional public records do not change the promotion decision. The
PRESLHY D3.6 instrument table supports a declared nominal-specification
sensitivity and bounds the expected sensor contribution, but not a claim of
calibration-complete accuracy; H2SAFE is a helium indoor surrogate; and the HyWAM publication does not
provide the complete Trial 10 event package. None supplies the missing
custodian-confirmed channel calibration plus measured site-geometry pair.

The KITopen/RADAR record is stronger than a bibliographic lead because its
metadata explicitly names source, atmospheric and receptor-related channels,
and confirms a public data package.  It still cannot close the gate from the
landing record alone: the archive must be selectively retrieved and its
per-trial files reconciled against the five promotion channels, with the data
custodian's calibration and geometry provenance retained.  Until that audit is
performed, the record remains `conditional` and is not added to the headline
validation score.

As a selective retrieval check, the public RADAR file
`trial_10_13-09-2019alldata.xlsx` was downloaded on 2026-10-09 through the
record's individual-file endpoint.  It is 575,705 bytes with SHA-256
`17987bc4e2f0a8d7bec40900063e3a3c2ce7a9ad6f9ab4014a5715ab04924aa`, exactly
matching the locally controlled Trial 10 workbook.  The workbook contains
five visible sheets: `Flexlogger` (209 x 99), `Draeger` (209 x 121), `Xensor`
(693 x 33), `LocalWeather` (93 x 18) and `Flowmeter` (209 x 15).  The headers
show one-second Flexlogger/Draeger/Flowmeter records, approximately 0.3-second
Xensor records and five-minute local-weather records.  No separate calibration
certificate, surveyed obstacle sheet, absolute UTC clock or event-level
geometry manifest is present in that workbook.  This confirms that the public
archive can reproduce the current conditional lane and that the remaining gap
is provenance/custodian evidence, not file discoverability.

The same byte-level retrieval check was then completed for the public Trial 20
and Trial 21 endpoints. Trial 20 is 758,361 bytes with SHA-256
`8b6a4a68ed731a70f109cc2a4d2bb0cf969f611b51fbca27eaf0a71cd0e8e662`; Trial 21
is 1,840,887 bytes with SHA-256
`11db48177c17ed2181be5259de6c66b78167d78d8e42940b4b011c0cea84fe7e`. Both
public downloads match the locally controlled workbooks byte-for-byte. The
machine-readable record is
`outputs/preslhy-public-retrieval-index-2026-10-09.json` now records all three
byte-level matches and pins the public D3.6 coordinate-report artifact
(6,518,995 bytes; SHA-256
`865b2b9f966e05023fbb581a4ed68f95f2353d6d6db2ad137aae7a9cf17c4f67`) with
the pages used for nominal accuracy, receptor coordinates, weather locations
and obstruction offsets. The report itself is not redistributed. The index
retains the calibration, clock and geometry promotion boundary and does not
redistribute the raw workbooks.

## Decision

The IJHE manuscript should retain a **conditional matched-event warning**.
Trials 10, 20 and 21 are sufficient to demonstrate that the intended
source/weather/receptor/clock reconciliation workflow is repeatable across
multiple accepted releases, but they do not support an operational accuracy or
obstacle-wake claim. Numerical convergence, an open-pad geometry declaration,
or a sensor registry without calibration cannot be promoted to field
validation. A stronger claim still requires custodian-confirmed calibration and
an independently reviewed site-geometry package.
