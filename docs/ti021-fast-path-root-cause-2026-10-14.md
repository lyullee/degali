# TI021 fast-path failure: measured causes and unresolved hypotheses

## Main finding

The 30-second path is not a faithful reduced version of the dense-gas solver.
It transports a passive species on a frozen air carrier and continues only
the downstream-face packets with a Gaussian kernel. Source collapse, density
feedback, gravity-current spreading and the evolving cloud inventory/momentum
are not carried into that exterior kernel. In addition, the applied air and
species ground boundaries disagree.

Instrumentation did not alter either baseline run. An explicitly labelled
counterfactual changed only the ground impermeability treatment in memory.
The production solver has not been repaired by this analysis.

## 1. Air penetrates the ground in the failed fast configuration

`tools/score_ti021_integral_patch_30s.py` selects `all_pressure_open`, whose
air-flow operator includes `z_min`. `integral_obstacle_patch.py` applies the
species open faces x-min, x-max, y-min, y-max and z-max, leaving z-min closed.
Thus zero species ground flux does not imply zero air normal velocity.

Measured maximum ground normal speed:

- 1,000 carrier steps: 2.4022 m/s.
- 5,000 carrier steps: 5.7273 m/s.
- Impermeable-ground diagnostic: exactly zero.

The air field is projected to small divergence with its open ground face
still present. Removing that face only in the species equation changes the
local carrier-volume divergence used by transport. This can cause artificial
compression/dilution even while integrated species mass closes.

## 2. The bad carrier lifts the dense species out of the patch

For 5,000 steps, the species centroid rises from 6.73 m at t=0 to 22.14 m at
t=10.10 s. The fraction of remaining species in the lowest 5 m falls from
38.5% to 3.5%. Ultimately 2,617.09 kg, or 79.5439% of the initial R12, exits
the 40 m roof. Only 19.6592% exits x-max.

The continuation reads only x-max face packets. Roof/side packets appear in
the conservation ledger but never enter downstream receptor concentration.
The bookkeeping is conservative; the exterior field prediction is incomplete.

## 3. Ground-only counterfactual separates this error from missing physics

The diagnostic excludes z-min from the pressure-open LBM faces, uses its
existing normal-population reflection, and fixes ground staggered normal
velocity at zero before projection. Sides and roof remain pressure-open.
This is an impermeable **free-slip** ground diagnostic, not a qualified
no-slip production patch.

| Quantity | Failed open ground, 5,000 steps | Impermeable ground, 5,000 steps |
|---|---:|---:|
| Carrier statistical residual | 0.57997 | 0.00936 |
| Declared 0.01 residual gate | fail | pass |
| Initial mass exported at roof | 79.5439% | 0.1651% |
| Dose MAPE | 98.8205% | 97.8668% |
| Predicted/observed positive-sensor dose sum | 0.5256% | 2.9597% |

A separate no-slip diagnostic overflowed during carrier startup and produced
no score. It is retained as `diagnostic-sealed-ground-flow5000-failure.json`;
no accuracy improvement is inferred from that failed intervention.

The successful impermeability probe demonstrates the boundary sensitivity,
but also shows that removing pathological roof export does not restore
measured sensor exposure. Ground boundaries are not the sole cause.

## 4. Dense-cloud dynamics are a leading missing mechanism in the fast path

For the declared initial 32 mole% R12/air mixture at 291.15 K and 100500 Pa:

- Mixture density: 2.4241 kg/m3.
- Ambient air: 1.2025 kg/m3.
- Density ratio: 2.0158.
- Vertical acceleration magnitude scale, g*(rho-rho_air)/rho: 4.9418 m/s2.
- Gravity/inertia velocity scale over 13 m: 8.015 m/s, versus wind 3.9 m/s.
- Idealized height-collapse timescale sqrt(2H/a): 2.294 s.
- Fast model's applied density-dependent acceleration: zero.

These are order-of-magnitude scales, not reconstructed experimental collapse
times. They establish that density forces are not small perturbations. Missing
collapse and horizontal hydrostatic-pressure-driven spreading also miss the
cloud's height evolution, lateral footprint and near-ground residence. Neutral
isotropic mixing cannot replace those dynamics.

The literature supports using layer-averaged mixture mass, horizontal
momentum, pollutant mass and energy equations for such clouds. DISPLAY-2
explicitly evaluated Thorney trials 8 and 21 with obstacles:
https://www.sciencedirect.com/science/article/abs/pii/S0304389403000116
Related two-layer gravity-current work:
https://www.mdpi.com/2311-5521/6/1/2

## 5. Failure occurs before the obstacle, and near-ground sampling aliases

S001 is at release-local (x,y,z)=(-14,14,0.4) m, radius 19.8 m from the source,
well inside the 50 m fence. Its observations and predictions are:

| S001 quantity | Observation | Old 84.9% run | Fast impermeable-ground probe |
|---|---:|---:|---:|
| Arrival (s) | 12 | 1.75 | 2.0 |
| Departure (s) | 305 | 27.625 | 10.5 |
| Dose (vol% s) | 1,832 | 72.990 | 7.278 |

The diagnostic fast detection window uses an unfiltered 0.1 vol% threshold;
the old timing score uses its 0.6 s averaging protocol. That observation-
operator difference cannot by itself account for the hundreds of seconds of
missing residence, but timing metadata still needs independent verification.

S001 at 0.4 m and S002 at 2.4 m share a horizontal location. Their observed
doses differ by 5.020x (1832 versus 364.936 vol% s). Both models produce
identical traces for the pair, because their lowest cell centre is 2.5 m and
both receptors clamp to that cell. The source/gravity-current layer and the
sensor observation operator need vertically resolved or layer-aware treatment.

The old 84.9% run's authoritative record actually has **5 m** inner cells,
not the 2.5 m previously stated in chat. Its outer cells are also 5 m.
The old solver does contain variable-density momentum and buoyancy, so
one must not claim that its entire error has the same cause as the newly
passive path. Shared sampling alias and severely short pre-obstacle residence
are confirmed; stratification-dependent mixing and collapse dynamics in the
old run remain hypotheses to isolate.

## Excluded gross explanations

- The holder volume is integrated with exact polygon/cell overlap:
  2048.2036 m3; initial R12 inventory 3290.1158 kg.
- Initial conversion uses R12 0.120913 kg/mol, not hydrogen molecular weight,
  and gives the intended 32 vol% initial cell concentration.
- Same workbook SHA-256, 600 s horizon and 84 positive-dose scoring mask.
- Global species ledger residual is about 2e-15 relative in the fast runs.

No factor-of-100 unit error or vanished source inventory was identified.
This does not rule out source-clock/sensor calibration issues.

## Repair order supported by this diagnosis

1. Make air, species and energy ground-face ownership/impermeability identical;
   stop scoring a non-converged or leaking carrier as a usable prediction.
2. Restore the actual dense-gas integral/shallow-layer source model and its
   cloud height, density, momentum and residence before the obstacle.
3. Transfer conserved state and all significant exiting face packets; do not
   count exported mass without continuing its concentration.
4. Resolve near-ground vertical structure with an anisotropic/stretched FV
   patch or an explicit layer reconstruction. Cubic 5 m LBM cells cannot
   distinguish these observation heights.
5. Re-score the identical sensor set and separate R12 mechanism evidence
   from H2 and cryogenic density/buoyancy-transition validation.

Artifacts: `outputs/integral-local-patch-2026-10-14/failure-cause-summary.json`
and `.png`; detailed unmodified and counterfactual instrumentation reports
are in the same directory.
