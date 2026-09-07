# Pre-registration: conservative-handoff numerical convergence

Date frozen: 2026-09-04

## Baseline

The accepted measured-throat phase-profile handoff at `S = 0.080 m` used 41
radial points, 1 mm maximum axial step and `2e-6` relative ODE tolerance. It
gave 0.0585% energy mismatch, 1.571% H2 half-width mismatch, 1.045 K centre
temperature mismatch and `1.25e-6` doubled-order energy-quadrature change.

## Confirmatory run

Repeat the identical physical source and atmosphere with the public strict
settings:

- 81 radial points,
- 0.25 mm maximum step,
- `5e-8` relative tolerance,
- the same `S = 0.080 m` handoff,
- default 32/64-point energy quadrature.

## Frozen decision rule

The strict run must still pass every original interface threshold. In
addition, relative to the baseline its:

- energy mismatch may change by at most 0.2 percentage point,
- H2 half-width mismatch by at most 1 percentage point,
- centre-temperature mismatch by at most 0.2 K, and
- local wind by at most `1e-9` relative.

No parameter or handoff position is changed if this confirmation fails.

## Results

The strict run completed in 255.6 s and passed:

| diagnostic | coarse | strict | change | limit |
|---|---:|---:|---:|---:|
| energy mismatch | 0.05852% | 0.05853% | +0.00001 percentage point | 0.2 point |
| H2 half-width mismatch | 1.571% | 1.571% | below displayed precision | 1 point |
| centre-temperature mismatch | 1.045 K | 0.943 K | -0.102 K | 0.2 K |
| energy quadrature residual | 1.25e-6 | 1.916e-6 | — | original 1e-5 screen |
| maximum native balance residual | 5.04e-16 | 2.09e-16 | — | original 1e-8 screen |

The local wind is identical at `1.687380256997 m/s`. The phase-profile
handoff is therefore accepted as numerically converged at this station.
