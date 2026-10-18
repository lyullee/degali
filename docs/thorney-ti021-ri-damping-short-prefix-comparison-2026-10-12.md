# TI021 Richardson-damping short-prefix comparison

Date: 2026-10-12  
Purpose: observation-blind numerical/physical sensitivity check of the opt-in
stable-stratification vertical-diffusivity candidate. This is not a validation
or promotion result.

## Frozen pair

Both runs use the same TI021 R12 property audit, workbook, 5 m grid, 0.02 s
time step, 2 s horizon, 60 coordinate-only receptors, and nonlinear tolerance
`1e-4`. The candidate changes only the explicitly opted-in damping parameters:
`alpha=5`, exponent `1`, and shear regularization `0.041 s^-1`. Observation
values are not used to configure either run.

| Run | Output | Wall time | Status | Max relative species-mass residual |
|---|---|---:|---|---:|
| Neutral baseline | `outputs/thorney-ti021-ri-damping-baseline-dx5m-2s-2026-10-12-v1/` | 14.90 s | `accepted_inner_transient_qualification` | `2.7643e-16` |
| RI-damped candidate | `outputs/thorney-ti021-ri-damping-alpha5-dx5m-2s-2026-10-12-v1/` | 15.32 s | `accepted_inner_transient_qualification` | `2.7643e-16` |

Both run records report all their computational gates passing. This only
establishes successful execution and conservation for this short prefix.

## Paired results

The 60 labels and all 101 time samples match exactly. Comparing the candidate
minus baseline sensor concentration traces gives:

- pooled trace RMSE: `0.002865 vol%-points`;
- pooled trace MAE: `0.0003974 vol%-points`;
- largest absolute sensor/time difference: `0.04951 vol%-points`;
- endpoint trace RMSE: `0.01076 vol%-points`;
- mean per-sensor peak change: `+0.00346 vol%-points` (15 higher, 3 lower,
  42 unchanged);
- maximum of the 60 sensor peaks: `2.6724` to `2.7219 vol%`;
- receptors reaching 4 vol%: `0/60` in both runs;
- maximum field concentration: `32.0 vol%` in both runs.

The candidate therefore does not show an improvement in this window. It
slightly raises the maximum receptor concentration and changes only a small
part of the sensor signal. Since this is an early prefix with no 4 vol%
receptor arrival, it cannot adjudicate the dense-cloud front, obstacle bypass,
recirculating wake, peak/dose timing, hazard volume, or full-event convergence.

## Decision

Keep Richardson damping **opt-in and non-default**. Do not cite this pair as
TI021 agreement, hydrogen validation, evidence of improved accuracy, or grounds
for promotion. The existing 0--20 s three-grid audit remains the controlling
grid-sensitivity evidence and continues to withhold full-event quantitative
use because hazard-volume and receptor convergence gates fail. A longer
paired run is justified only after a predeclared cost/benefit decision; it must
retain the frozen parameters, score observation data only after execution, and
reassess all obstacle-region and hazard-volume metrics rather than selecting
the candidate from peak concentration alone.

