# Observed-wind four-cell contrast (same keys) — 2026-09-06

Scope: trial 10 with fixed SOURCE/HANDLE parameters and same 19-temperature rows/5-concentration arcs/5-geometry comparisons.
Key files:

- `reference/preslhy/yawed_trial10_step002_2026-09-06/complete.json`
- `reference/preslhy/observed_wind_magnitude_step002_2026-09-06/complete.json`
- `reference/preslhy/observed_wind_vector_step002_2026-09-06/complete.json`
- `reference/preslhy/observed_wind_vector_step001_2026-09-06/complete.json`

Baseline (CONTROL 19T mean MAE): `27.0312 K`

| Case | temp mean MAE (K) | Δ vs control | temp minimum MAE (K) | temp p05 MAE (K) | MG | VG | width ratio | centre MAE |
|---|---:|---:|---:|---:|---:|---:|---:|
| yaw .02 | 19.0453 | -29.54% | 23.3639 | 20.2822 | 0.932706 | 1.016145 | 1.235037 | 0.177467 |
| magnitude .02 | 29.2721 | +8.29% | 19.1587 | 17.6668 | 0.774623 | 1.142112 | 1.093047 | 0.188774 |
| vector .02 | 23.7324 | -12.20% | 20.2293 | 18.0680 | 0.871162 | 1.050781 | 1.292207 | 0.200676 |
| vector .01 | 23.7332 | -12.20% | 20.2313 | 18.0700 | 0.871186 | 1.050776 | 1.292243 | 0.200680 |

Findings:

- `measured_vector .01` and `.02` are nearly identical for this case; no meaningful step-size signal.
- `yaw` gives best temperature MAE but worsens near-core trajectory/width structure.
- `magnitude` hurts temperature at this stage and does not recover concentration.
- `vector` improves concentration relative to control and does not improve temperature enough to justify geometry/trajectory side-effect alone.

No adoption has been made from this four-cell contrast alone.
