# PRESLHY E3.5 Trial 10 conditional field-evidence package

The public PRESLHY E3.5 workbook for Trial 10 contains a Flowmeter source
history, a LocalWeather record and an approximately 3 Hz Xensor concentration
history.  The D3.6 report supplies the Xensor receptor coordinates and heights,
the 18 s sampling-line delay, the common logging architecture and the nominal
instrument accuracy table (Xensor XEN-5320: `±3% of full scale`; cryogenic
Coriolis flow: estimated `±3% of reading`).  These streams use the recorded
wall-clock fields in the workbook; the documented 18 s sampling-line delay is
kept as an explicit observation-operator input.  The D3.6 report also lists
obstruction trials at 0.16 m and 0.18 m from the release, but does not provide
a complete surveyed three-dimensional obstruction mesh.

The raw workbook and report are not copied into this repository.  The
reproducible extraction tool is
`tools/build_preslhy_trial_evidence_bundle.py`.  It writes a derived bundle
under `outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/` when the custodian
supplies the public files.  The bundle contains source, weather, receptor,
open-pad geometry, common-clock and sensor-registry records, plus a
hash-pinned `FieldEvidenceManifest`.

## Reproduction

```powershell
$env:PYTHONPATH = 'src'
python tools/build_preslhy_trial_evidence_bundle.py `
  --workbook '.../trial_10_13-09-2019alldata.xlsx' `
  --report '.../PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.22.pdf' `
  --output-dir outputs/preslhy-e35-trial10-evidence-2026-10-09-v2

python -m degali.cli field-audit `
  outputs/preslhy-e35-trial10-evidence-2026-10-09-v2 `
  --output outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-audit-execution.json `
  --require-complete

$bundle = 'outputs/preslhy-e35-trial10-evidence-2026-10-09-v2'
python -m degali.cli field-evidence-manifest-create `
  "$bundle/field-audit-execution.json" `
  --selected-path "source_boundary=$bundle/source_boundary.csv" `
  --selected-path "weather=$bundle/weather.csv" `
  --selected-path "obstacle_geometry=$bundle/obstacle_geometry.json" `
  --selected-path "receptor_observations=$bundle/receptor_observations.csv" `
  --selected-path "common_clock=$bundle/common_clock.json" `
  --manifest-id preslhy-e35-trial10-manifest-2026-10-09-v3 `
  --event-id preslhy-e35-trial-10-2019-09-13 `
  --dataset-id preslhy-e35-trial10-xensor-window `
  --observed-row-count 5053 `
  --source-boundary-id preslhy-e35-trial-10-2019-09-13-source `
  --weather-id preslhy-e35-trial-10-2019-09-13-weather `
  --obstacle-geometry-id preslhy-e35-trial-10-2019-09-13-open-pad-geometry `
  --receptor-geometry-id preslhy-e35-trial-10-2019-09-13-receptors `
  --temporal-operator-id xensor-release-window-plus-18s-0.3s `
  --common-clock-id preslhy-e35-trial-10-2019-09-13-shared-wall-clock `
  --sensor-set-id preslhy-e35-xensor-trial10-a1 `
  --operator-id degali-reproducible-bundle-builder-v1 `
  --sensor-registry-path "$bundle/sensor_registry.csv" `
  --sensor-calibration-status specification_only `
  --scope lh2_free_field `
  --output "$bundle/field-evidence-manifest.json"

python -m degali.cli field-evidence-manifest-verify `
  "$bundle/field-evidence-manifest.json"
```

The local run returned `status=candidate_complete`; the channels are
intentionally split across source, weather, geometry, receptor and clock
files, so the accountable manifest selection is explicit rather than inferred
from filenames. The subsequent manifest returned
`evidence_readiness=conditional`, `selected_files_verified=true` and
`promotion_allowed=false`: the registry records the sensor set and response
interval, and the explicit `sensor_calibration_status=specification_only`
field records that the public workbook does not provide an event-level
calibration certificate. The current manifest file is
`outputs/preslhy-e35-trial10-evidence-2026-10-09-v2/field-evidence-manifest.json`
with SHA-256
`0b59a6dc21515c84ff9043b1cdc638f12757ec8bfff1e7051e6a9892583adf84`.

| item | value |
|---|---|
| Event | `preslhy-e35-trial-10-2019-09-13` |
| Common clock | `preslhy-e35-trial-10-2019-09-13-shared-wall-clock` |
| Source rows | 208, with the recorded rate and a non-negative physical-rate field |
| Weather rows | 92, five-minute LocalWeather record |
| Xensor rows | 5,053 across 31 mapped receptors |
| Observation operator | source release window plus 18 s transport delay; 0.3 s Xensor sampling |
| Raw workbook SHA-256 | `17987bc4e2f8a0d7bebc40900063e3a3c2ce7a9ad6f9ab4014a5715ab04924aa` |
| Coordinate-report SHA-256 | `865b2b9f966e05023fbb581a4ed68f95f2353d6d6db2ad137aae7a9cf17c4f67` |

## Boundary

This is a **conditional evidence and replay package**, not a new headline
validation score.  The public material does not include a channel-level
calibration certificate in the extracted bundle; LocalWeather is five-minute
rather than high-frequency mast data; and the obstacle record is an explicit
open-pad declaration rather than a measured three-dimensional site survey.
Negative Xensor baselines are retained in a raw column and clipped only in the
physical concentration field. Consequently the manifest is conditional for
traceability and remains deliberately non-promotable. The manuscript must
continue to state that stronger operational or obstacle-resolved claims need
custodian-confirmed calibration and a matched site-geometry package.

The public source is the PRESLHY E3.5 dataset (DOI
[10.35097/1481](https://doi.org/10.35097/1481)); raw third-party files remain
outside the publication snapshot unless redistribution rights are confirmed.

The D3.6 source report records the nominal sensor specifications, not
event-level calibration certificates.  Therefore the specification is useful
for an uncertainty/sensitivity bound, while `sensor_calibration_status` remains
`specification_only` and the manifest remains non-promotable.

## Related accepted-trial replays

The same extraction command and D3.6 coordinate report were replayed against
the locally controlled public-archive workbooks for accepted Trials 20 and 21.
Their derived bundles and manifests are retained outside the journal-facing
source tree and indexed in
`outputs/preslhy-conditional-manifest-index-2026-10-09.json`:

| Event | Receptor rows | Manifest status | Promotion |
|---|---:|---|---|
| Trial 20 | 19,344 | `conditional` | `false` |
| Trial 21 | 2,232 | `conditional` | `false` |

These replays demonstrate workflow repeatability across multiple accepted
releases. They do not add headline validation points and retain the same
calibration, absolute-UTC and measured-geometry limitations described above.

The public RADAR workbooks for Trials 20 and 21 were also retrieved through
their individual-file endpoints and matched byte-for-byte to the locally
controlled copies. The combined retrieval audit is
`outputs/preslhy-public-retrieval-index-2026-10-09.json`. This strengthens
traceability of the replay inputs only; it does not convert the bundles into
calibration-complete field validation or permit promotion.
The same index pins the non-redistributed D3.6 report used for the replay
(6,518,995 bytes, SHA-256
`865b2b9f966e05023fbb581a4ed68f95f2353d6d6db2ad137aae7a9cf17c4f67`; pages
14–15, 43–46 and 50), preserving the coordinate and instrument-provenance
boundary without treating the report as event-level calibration.
