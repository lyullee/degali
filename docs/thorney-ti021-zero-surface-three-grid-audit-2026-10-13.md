# TI021 zero-thickness obstacle: matched three-grid impact audit

Machine-readable report and figure:

- `outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v1/ti021-three-grid-impact-convergence.json`
- `outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v1/ti021-three-grid-impact-convergence.png`

Matched 3-view × 3-grid maximum-concentration contours:

- `outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v2/ti021-three-grid-maximum-concentration-contours.png`
- `outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v2/ti021-three-grid-maximum-concentration-contours.json`

The contour panels use shared physical planes (`z=2.5 m`, `y=100 m`, `x=60 m`), a common 0–32 vol% scale, and a black 4 vol% contour. Values are each cell's maximum over the full 0–20 s interval, not a simultaneous snapshot; the plotted 4 vol% contour is therefore an “ever exceeded” slice and must not be presented as an instantaneous cloud boundary or a validated hazard contour. Cell-centred values are interpolated to the same plane on all three grids.

## Matched runs

The audit reuses finalized, observation-blind TI021 initial-cloud runs at 5.0, 2.5, and 1.25 m. Each covers 0–20 s with `dt=0.02 s`, corrector tolerance `1e-4`, the same zero-thickness semi-circular surface, identical workbook/property SHA-256 inputs, quiescent initial-cloud momentum, and exact 2:1 grid refinement. All run numerical gates pass. This isolates spatial-grid sensitivity for the idealized zero-thickness model; the real scaffold/tarpaulin member geometry and tarp gauge remain unidentified.

## Result: important outputs are not grid-converged

- Common-fluid concentration-field L1/RMSE errors contract under refinement.
- Maximum-ever volume at or above 4 vol% increases from **14,875 m³** (5 m) to **15,000 m³** (2.5 m) to **15,431.6 m³** (1.25 m). The sequence is monotonic but nonconvergent, so no defensible GCI is available.
- Only **6/25 (24%)** eligible receptor peak triplets are monotonic-convergent (median fine-grid GCI **56.1%**, p90 **178.9%**). For integrated dose, only **9/25 (36%)** converge (median GCI **244.0%**, p90 **1,383.8%**). Arrival time is less sensitive: **12/25 (48%)** converge (median GCI **7.68%**, p90 **51.9%**), still short of the declared 70% gate.
- At the side-tip bypass region, the fine grid first resolves a localized maximum of **5.35 vol%** and a 4 vol% volume of **19.53 m³**; both coarser runs report zero there. This is a concrete warning that the bypass/hazard contour remains unresolved at operationally tempting coarse spacings.
- Status is `impact_window_grid_sensitivity_not_converged_full_event_withheld`. The field-average contraction is not enough to promote a 4 vol% contour or sensor peak/dose estimate.

This is a 20 s dense-gas R12 obstacle-surrogate numerical audit, not a hydrogen accuracy result and not the 600 s full-event grid-convergence audit. Next: refine the obstacle-tip and lee-wake regions or use a defensible adaptive/local mesh, then rerun the frozen impact audit without tuning against the scores. The reported physical fence geometry still needs the scaffold and tarp construction registered or its omission justified.
