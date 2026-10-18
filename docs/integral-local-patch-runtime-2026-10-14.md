# Grid-free exterior and local obstacle patch: measured checkpoint

The new executor in `src/degali/addons/integral_obstacle_patch.py` advances
only a local obstacle patch. It accepts an interval-mean species face flux
from an exterior transport callback. The benchmark supplies that flux with a
Gaussian cross-section integral and exact advective delay. This is a passive
isothermal exterior model, not the DEGADIS dense-gas integral closure.

The patch carrier is the existing viscous D3Q19 obstacle solver. Its startup is
included in the measured runtime. A time-averaged carrier is retained during
the subsequent scalar transient, including its reverse-flow cells. Molecular
plus resolved SGS diffusivity is applied with a declared Schmidt number.
This does not solve time-varying buoyancy or pressure feedback from the source.

The inlet packet is imposed in the finite-volume face-flux operator. Inlet,
outlet and other boundary fluxes are independently integrated and reconciled
with the patch field inventory. No concentration overwrite, clipping, or
post-run mass renormalization closes the balance. Outlet face distributions
are continued to downstream receptors using a grid-free Green function.
Lateral/top exits are inventoried but their downstream concentration is not
yet continued. A reverse-flow outlet is rejected and requires a longer patch.

## Measured case

- 600 s horizon, 0.02 kg/s atmospheric H2 for the first 120 s, then shutdown.
- Uniform wind 2 m/s; ambient 293.15 K, 101325 Pa.
- Reference extent 100 x 24 x 20 m.
- Patch x=30..50 m, y=6..18 m, z=0..10 m.
- Block x=36..40 m, y=10..14 m, z=0..4 m.
- Exterior time-marched grid cells: zero.
- Both refinement cases have a 15 s physical carrier spin-up and a final
  3 s velocity averaging window.

| Quantity | 1 m local mesh | 0.5 m local mesh |
|---|---:|---:|
| Patch cells | 2,400 | 19,200 |
| Cold-start compute time | 3.333 s | 60.163 s |
| Including NPZ output | 3.463 s | 61.046 s |
| Maximum relative species ledger residual | 2.11e-14 | 6.14e-14 |
| Reverse-flow lee cells | 17 | 240 |

The 1 m case meets the 30 s runtime target on this workstation. The 0.5 m
case does not. These are actual 600 s runs, not extrapolations from a short
prefix. They are a new test case; the TI021 runtime is not a directly matched
speed comparison.

The 1 m sensor-peak changes relative to 0.5 m are 9.32% front, 15.62% roof,
19.75% side, 21.68% wake and 24.60% downstream. Thus the all-sensor 20% peak
refinement gate fails. This is numerical sensitivity, not experimental error.
Neither result has a recalculated independent experimental MAPE.

The obstacle-free local-grid diagnostic adds 5.12..16.32% peak error relative
to the analytic solution at 1 m. The proposed 1% obstacle-free interface gate
is not passed. This is evidence to improve reconstruction/sampling and to
dispatch an entirely obstacle-free case directly to the exterior integral
model, rather than claiming that conservative mass balance implies accurate
point concentration.

## Reproduction

From the repository with its `.venv` and `PYTHONPATH=src`:

```powershell
.\.venv\Scripts\python.exe tools/benchmark_integral_obstacle_patch.py --duration-s 600 --flow-steps 1000 --output outputs/integral-local-patch-2026-10-14/coarse-600s.json
.\.venv\Scripts\python.exe tools/benchmark_integral_obstacle_patch.py --duration-s 600 --spacing-m 0.5 --flow-steps 2000 --output outputs/integral-local-patch-2026-10-14/fine-600s.json
.\.venv\Scripts\python.exe tools/compare_integral_obstacle_patch.py --coarse outputs/integral-local-patch-2026-10-14/coarse-600s.json --fine outputs/integral-local-patch-2026-10-14/fine-600s.json --output outputs/integral-local-patch-2026-10-14/comparison.json
```

The JSON and NPZ artifacts include concentration fields, sensor traces,
face outlet rates, inventories, runtime and carrier residuals. The figure is
`outputs/integral-local-patch-2026-10-14/comparison.png`.
