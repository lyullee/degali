# FFI paired obstacle-effect diagnostic — 2026-10-14

This audit runs the same fast three-dimensional transport path for the public
FFI hydrogen channel control case A and obstacle case O1 over 60–120 s, then
compares the sensor-wise change `O1 − A`.  The purpose is to separate a common
absolute source/concentration error from the obstacle response.  It is not a
facility-validation score because the published O1 beam assembly, leading-edge
datum and source/sensor clock operator remain unresolved.

## Result

The machine-readable result is
`outputs/ffi-o1-control-obstacle-delta-2026-10-14.json`.

| Quantity | Result |
|---|---:|
| Receptor pairs retained | 24 |
| Cell size / grid | 0.05 m / 116×18×16 |
| Runtime (control + obstacle) | 11.96 s |
| Predicted Δ sensor MAE | 0.328 vol%-points |
| Predicted Δ sensor RMSE | 0.537 vol%-points |
| Predicted Δ bias | −0.062 vol%-points |
| Signed correlation of predicted vs measured Δ | −0.098 |
| Maximum mass residual | 1.21×10⁻¹⁷ kg |
| Promotion | false |

The numerical conservation ledger is clean, but the conditional beam surrogate
does not reproduce the measured sensor-wise obstacle effect: its signed pattern
has effectively no correlation with the paired observation, and its predicted
changes are generally much smaller than measured changes at several top
receptors.  This is stronger evidence than the absolute MAPE alone that the
remaining error is dominated by unresolved geometry/clock/near-roof transport,
not merely a uniform source-rate scale.

The next quantitative step therefore requires the machine-registered beam
cross-sections and datum, followed by a common-clock replay.  Parameter fitting
to the 24 deltas is deliberately not performed.

