# Matched LH2 comparison and transfer-screen evidence (2026-10-03)

This is a derived, no-fit evidence package for the manuscript. The comparator
is DEGALI's local `corrections=False` historical reconstruction, **not** the
original DEGADIS Fortran code, HyRAM+, SLABx, or an independent implementation.
Both paths receive the same reported source, meteorology, sensor coordinates,
60 s model averaging setting and arc-maximum spatial observation operator.

## Reproduction

From the repository root:

```powershell
$env:PYTHONPATH = 'src'
.\.venv\Scripts\python.exe tools\audit_applied_energy_comparison.py `
  --output-dir tmp\applied-energy-evidence-reproduction `
  --repeats 3 --bootstrap 5000
```

The seed is fixed to `20261003`; `manifest.json` records input SHA-256 hashes,
method and timing environment. Scientific result values and percentile
intervals reproduce deterministically. Wall-clock timings are machine- and
load-dependent. The original 2026-10-03 snapshot is in
`outputs/applied-energy-evidence-2026-10-03/` as CSV and SVG. A regenerated
2026-10-09 snapshot (`outputs/applied-energy-evidence-2026-10-09/`) uses the
same fixed seed and adds trial-cluster VG intervals to `aggregate_metrics.csv`;
the headline scores are unchanged.

For the separate external-model lane, provide the accepted SLABx Test 4
`sensor_predictions_n640.csv` from the adjacent project:

```powershell
.\.venv\Scripts\python.exe tools\audit_test4_slabx_degalix.py `
  --slabx-sensor-csv '<slabx-repository>\model-comparison\outputs\test4_slabx_matched_sensor_gate_full\sensor_predictions_n640.csv' `
  --output-dir tmp\test4-external-reproduction
```

The audit also writes `test4_field_model_comparison.json`. It fingerprints the
paired CSV and explicitly labels the source-boundary and time-operator
mismatches, so its numerical differences cannot be mistaken for a direct
common-source validation score.

## Matching and limitations

| Dimension | Matched? | Treatment |
|---|---|---|
| Source and wind | Yes, between model paths | E3.5 window-mean flow; FFI published low-mast wind with stated 10 m reference-height assumption |
| Sensor geometry | Yes | Predict at each E3.5 sensor; FFI at reported arc heights, then take arc maximum |
| Model averaging | Yes, 60 s | Same setting for both model paths |
| Observation time operator | **No** | Field concentration is a reported release-window peak, not a synchronized 60 s average. Raw time-aligned channels are unavailable in this package. |

Thus “matched” below means a **paired model-path comparison** and a common
spatial observation operator. It does not mean the model and experimental time
operators are identical. The field scores are conditional on that mismatch.

## Results

| Data / path | Tests | Arcs | MG (observed/predicted) | VG | FAC2 | 95% cluster-bootstrap MG | 95% cluster-bootstrap VG | 95% cluster-bootstrap FAC2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| E3.5 historical reconstruction | 9 | 62 | 1.156 | 23.101 | 0.774 | 0.721–2.528 | 1.149–7,206 | 0.563–0.932 |
| E3.5 corrected | 9 | 62 | 1.047 | 1.425 | 0.839 | 0.771–1.416 | 1.132–1.882 | 0.677–0.983 |
| FFI historical reconstruction | 2 | 6 | 12.600 | 4283.065 | 0.000 | 3.921–40.495 | 7.307–2,510,435 | 0.000–0.000 |
| FFI corrected | 2 | 6 | 1.245 | 1.373 | 0.833 | 1.031–1.503 | 1.075–1.754 | 0.667–1.000 |

These are percentile bootstrap intervals with the *trial* as resampling unit
(5,000 draws), not confidence bounds on the true atmospheric source or model
form. The wide E3.5 historical VG interval reflects heterogeneous trial
residuals; the corrected interval is narrower but remains a sensitivity
diagnostic. The FFI interval is especially fragile with only two trials. Trial
rows also include within-trial arc-resampling intervals; arcs are correlated
and those are descriptive, not independent experimental uncertainty.

Improvement is not uniform. Corrected-path trial MG remains 2.366 on E3.5
Trial 20 and 0.473 on Trial 21, and FFI Test 6 at 30 m remains 21.0 observed
versus 6.03 vol % predicted (3.48-fold low). The full trial and arc tables
must accompany any aggregate claim.

## Runtime

Wall-clock timing includes source thermodynamics, plume integration and
receptor extraction, with each variant warmed once before measurement. Three
runs per trial were made sequentially in one process; the table is the median
of trial-level medians, with minimum–maximum across trials. It does not claim
a speedup over an external code or transient CFD.

| Data / path | Median per trial | Across-trial range |
|---|---:|---:|
| E3.5 historical reconstruction | 0.373 s | 0.362–0.387 s |
| E3.5 corrected | 0.351 s | 0.344–0.359 s |
| FFI historical reconstruction | 1.666 s | 1.505–1.827 s |
| FFI corrected | 1.567 s | 1.538–1.595 s |

These small within-session differences are not treated as a statistically
established speed advantage. A prior run in the same workspace gave medians
of 0.159 s (E3.5 corrected) and 0.699 s (FFI corrected), showing that
background load materially affects wall time; the current snapshot's exact
environment and per-trial records are in `manifest.json` and
`runtime_by_trial.csv`. A dedicated controlled benchmark is still needed
for a journal speed claim.

## Separate SLABx Test 4 sensor comparison

The adjacent `SLABx_LH2/model-comparison` project contains an accepted
640/64-refinement SLABx Test 4 calculation at 30 physical sensors. Its
published setup uses 0.828 kg/s, 25.4 mm, 6.7 m/s wind at 10 m, 3.3 °C,
75% relative humidity, 0.03 m roughness, neutral stability and a 275 s
sensor-mean window. The DEGALI corrected path was rerun with those same
bulk inputs and sensor coordinates; it was **not** combined with the
low-mast/peak-concentration six-arc comparison above.

| Model | MAE (vol %-point) | RMSE (vol %-point) | LFL TP / FN / FP |
|---|---:|---:|---:|
| SLABx (accepted 640/64 run) | 0.648 | 1.211 | 3 / 0 / 0 |
| DEGALI corrected, quasi-steady | 1.539 | 2.679 | 3 / 0 / 3 |

This is a **conditional** common-receptor comparison, not a strictly common
temporal response: SLABx reports a 275 s time mean, while DEGALI's 275 s
parameter sets a steady meteorological averaging closure. The codes also
use their own native two-phase/flash source mappings. The SLABx input file
and its SHA-256 fingerprint, exact paired rows and DEGALI script are recorded
in `test4-external/`. The strict general comparison execution record,
`field-model-comparison-execution-v3.json`, also fingerprints the manifest
that documents the missing CSV averaging-time column and native-model limits;
it intentionally returns `conditional` and model selection `withheld`. The
external run did not record comparable wall time,
so no cross-code speed ratio is claimed. The three extra DEGALI false
positives matter for an operational LFL detection task despite a small
30-sensor sample.

## LH2 transfer consequence case

FFI Test 6 is a physical horizontal LH2 bunkering/transfer-like release:
0.833 kg/s through a 25.4 mm orifice, 2.3 m/s published low-mast wind, 0.5 m
release height. At a fixed 1.5 m-high downwind centreline receptor, the
computed 4 vol % LFL exit is **12.6 m** with the historical reconstruction and
**34.3 m** with the corrected path (+21.7 m, 2.72×). At 20 m the respective
model values are approximately 0.46 and 16.4 vol %, so a hypothetical 20 m
screening line would change category. This is a model-choice demonstration,
**not** a safety setback or authorization to operate: the corrected model's
observed Test 6 30 m arc error, time-operator mismatch and source/wind
uncertainty make an operational clearance decision unsupported.

## Figure and table inventory

- `ffi_paired_arc_errors.svg`: paired FFI observed/predicted ratios at 30,
  50 and 100 m; each test labelled.
- `e35_trial_bias.svg`: trial MG and within-trial descriptive intervals.
- `transfer_lfl_screen.svg`: both Test 6 model profiles at 1.5 m, and 4% line.
- `paired_arcs.csv`: exact common-arc source for every plotted/model statistic.
- `trial_metrics.csv`: per-test metrics, descriptive within-test intervals and
  median runtime.
- `aggregate_metrics.csv`: pooled scores and trial-cluster intervals.
- `runtime_by_trial.csv`, `transfer_test6_lfl_profile.csv`, `manifest.json`:
  timing, case trace and provenance.
- `LH2_validation_evidence.xlsx`: reader-facing copy of summary and tables.
  The artifact renderer crashed in this Windows runtime; the exported XLSX
  was independently reopened and cell dimensions, key values and formula
  were inspected. Visual layout in Excel remains to be checked before journal
  circulation.

## Remaining paper gate

A genuinely independent same-input/geometry/**time-operator** comparison
requires either a runnable external model plus a defined model observation
operator and synchronized field time series, or a narrower explicitly
time-averaged subset. Neither is supplied by the current files. Do not call
the local historical branch an independent external benchmark.
