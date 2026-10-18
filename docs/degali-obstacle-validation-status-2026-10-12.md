# DEGALI obstacle-resolved validation status

**Snapshot:** 2026-10-12  
**Scope:** transient, obstacle-resolved Thorney Island TI021 impact-window run

This document records the current evidence boundary. A numerical gate passing is
not treated as field validation unless the corresponding observation and
convergence gates also pass.

The AIJ Case H passive-ethylene diagnostic was extended with a same-grid,
same-window SGS scalar-diffusion sensitivity. On 582 common sensor means over
8--12 s, adding Smagorinsky gradient diffusion with `Sc_t=0.7` changes RMSE
from 2624.5 to 2280.7 ppmv and bias from +745.6 to +688.0 ppmv; relative
species-mass residual remains below `5.5e-14`, and runtime is 12.04 s. This
13.1% RMSE decrease is a useful mechanism sensitivity only: both predictions
remain strongly high, the diffusivity is uncalibrated, the D3Q19 comparison
does not apply source/nozzle momentum, and the 120 s observation window is not matched. It
does not validate ethylene concentration or hydrogen dispersion.

A generic ethylene--air binary EOS was then exercised with the FV momentum /
mass-flux projection and the exit momentum inferred from the reported flow
rate and 4 mm port (`0.46394 m/s`). A matched 0.5--1.0 s startup control with
and without that momentum closes the maximum variable-density continuity,
global mass, and momentum-balance residuals at `1.42e-13 kg m-3 s-1`,
`1.72e-16 kg s-1`, and `1.35e-16 kg m s-1`, respectively. The paired mean
sensor-field change is only `-0.031 ppmv` (paired RMSE difference `0.222
ppmv`) over 582 points. This short model startup is not comparable to the
experiment's 120 s mean, and the exit velocity itself is inferred rather than
measured; it proves that the coupled code path executes conservatively, not
that the source boundary is physically resolved or validated. The FV path's
coarse runtime is about 30 s per simulated second on this machine.

For a stronger sensor comparison, both solver tiers were then run on the same
25 mm mesh, 1 ms scalar time step, 8--12 s averaging window, source rate and
Smagorinsky/Schmidt scalar settings. Over 582 common sensors, the D3Q19 carrier
with SGS scalar diffusion scored RMSE/bias `2283.0/+688.9 ppmv` (14.91 s);
the generic ethylene--air EOS/FV momentum path with inferred nozzle momentum
scored `1007.8/+263.6 ppmv` (270.62 s), with Pearson `r=0.743` versus `0.556`.
Finite-volume continuity and global-mass residuals remained `1.55e-13 kg
m-3 s-1` and `2.05e-16 kg s-1`. This is a matched-grid component improvement
but changes EOS/pressure/momentum as a package. A second matched FV pair then
turned source momentum off while holding the other settings fixed: RMSE changed
only from `1008.095` to `1007.762 ppmv` (0.033%), and the paired sensor-field
RMSE difference was `0.373 ppmv`. Thus the inferred momentum term does not
explain the large tier gain in this setup; remaining candidates include
variable-density EOS/pressure coupling and interactions with source/grid
regularization. The FV mean concentration is still about `3.16x` the measured
mean, and its 8--12 s model window does not match the experiment's 120 s
acquisition. Promotion and hydrogen claims remain withheld.

## Evidence matrix

| Requirement | Current evidence | Status |
|---|---|---|
| Conservation-based 3-D obstacle solver | `src/degali/addons/resolved_obstacle_flow_3d.py`; 52 related regression tests pass | **pass** |
| Obstacle-aware stratification diagnostic | Fine-grid 20 s gradient-Richardson audit passes with solid-cell exclusion; diagnostic-only, no flux damping | **pass for computed window** |
| Mass, momentum, species and sensible-enthalpy residuals | 0.02 s run checkpoints at 1 and 2 s; residuals are below the declared numerical gates | **pass for computed window** |
| Stability and corrector convergence | Checkpoints report zero non-converged corrector steps and stability numbers below one | **pass for computed window** |
| Pressure-outlet momentum backflow state | `variable_density_momentum_3d.py` now uses ambient-at-rest velocity on high-x backflow and interior velocity only on outflow; paired boundary tests and the momentum/obstacle/binary-corrector regression group pass (74 tests) | **numerical boundary fix verified; no new facility score attributed** |
| Closed-rear zero-ambient FV startup | New explicit zero-wind mode separates physical ambient wind (0 m/s) from positive lattice velocity scale, requires projected inlet and variable-density FV momentum, and passes a one-step H₂-source/positive-buoyancy smoke with continuity, species-mass and energy gates | **component smoke only; no FFI channel-scale score** |
| Updated embedded-arc wall shear closure | New 5 m, 20 s zero-surface transient activates smooth-wall Reichardt stress (max 0.00718 Pa); all conservation gates pass. Same-grid pre/post comparison gives concentration-field RMSE 6.93e-6 and sensor-trace RMSE 6.12e-6 vol%-points | **execution verified only; calibration and wake validation withheld** |
| Front/top/side concentration fields | 2 s segment field and 3-view contour at `outputs/thorney-ti021-dt0p02-two-second-contours-2026-10-12/`; 5 s maximum-field contour at `outputs/thorney-ti021-dt0p02-five-second-contours-2026-10-12/` | **pass for computed window** |
| 4 vol% geometry | 2 s field gives 5,023.4 m³; full-event geometry is not yet accepted | **conditional** |
| Three-grid spatial convergence | Matched-dt 0--20 s audit holds the numerical zero-thickness surface fixed and field norms contract, but hazard/receptor gates fail. Source describes tarpaulin on scaffold; tarp gauge and scaffold sections are unreported, and scaffold is omitted from the model | **withheld** |
| Two-time-step sensitivity | Same 0--20 s, dx=1.25 m comparison at dt=0.02 and 0.015625 s; field/receptor RMSE and input SHA-256 values recorded | **pass for computed window** |
| Hazard-volume convergence | Zero-surface matched-dt maximum-envelope volumes are 14,875 / 15,000 / 15,431.64 m³; qualifying monotonic GCI is absent | **withheld** |
| Receptor peak/dose convergence | Spatial GCI ratios are 2:1; monotonic fractions are 0.24 (peak) / 0.36 (dose), median fine GCI 56.1% / 244.0% | **withheld** |
| Independent TI021 observation comparison | The 2.5 m resolved screening run covers the complete 0--600 s record; the 1.25 m finest impact-window run covers 0--20 s only | **conditional / finest-grid withheld** |
| Cross-file field evidence manifest | PRESLHY Trial 10 five-channel manifest re-verified; all selected file digests match | **conditional** |
| Full event promotion | A 600 s, 2.5 m obstacle-resolved screen exists, but there is no accepted 600 s fine-grid run or hazard/sensor convergence package | **withheld** |

## Current runs

### Pressure-open boundary correction (2026-10-12)

Inspection of the finite-volume momentum face flux found that the high-x
pressure outlet multiplied negative (inward) mass flux by the adjacent cell's
velocity. That implicitly imported interior momentum from the exterior during
backflow, whereas this outlet is coupled to quiescent ambient air. The boundary
now uses zero exterior velocity for inward mass flux and retains the adjacent
cell velocity for outward flux. Focused tests exercise both signs, and the
variable-density momentum, resolved-obstacle-flow and binary-mixture-corrector
regression group passes (74 tests total). The TI021 nested outer transient below
uses a separate flow implementation, so this code correction is not counted as
a change to that run's results; it must be included in any subsequent run using
`advance_variable_density_momentum_3d` before its facility metrics are scored.

### Zero-ambient-wind finite-volume startup (2026-10-12)

The resolved solver previously required positive `wind_speed_m_s` even when a
case had a physically quiescent ambient and a closed rear boundary. It now
allows `wind_speed_m_s=0` only with finite-volume momentum, projected-inlet
initialization, and an explicit positive `velocity_reference_scale_m_s` used
solely for numerical velocity/time and pressure-coefficient normalization. Any
declared physical inlet profile must remain zero; forced-flow inlet histories
are rejected. A one-step 3-D test with a source and an obstacle verifies zero
rear-face velocity, positive H₂ buoyancy, finite-volume continuity, species
mass, and energy residual gates. The resolved-solver plus variable-density-
momentum suite passes 68 tests. This is startup-path evidence only: it does not
yet qualify the FFI channel outlet, turbulence closure, beam geometry, clocks,
or concentration score.

The updated zero-surface wall-stress closure has a short, full-gate execution
witness at
`outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-2s-dt0p02-2026-10-12-v1/`.
It covers 2 s on the 5 m / 12,800-cell grid in 16.53 s wall-clock; the maximum
reported shear is 0.00284 Pa, maximum friction velocity 0.0486 m/s, relative
species-mass residual `2.76e-16`, enthalpy residual `0`, and momentum-balance
residual `2.43e-9 kg m s^-1`. This confirms execution and conservation with the
new closure, not wall-stress calibration or a separated-wake benchmark. The
20 s, 1.25 m field run predates the wall-shear code change and must not be
described as evidence for that closure.

The same closure has now been run through 20 s at 5 m in
`outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-20s-dt0p02-2026-10-12-v1/`.
Runtime was **108.28 s**; maximum shear was **0.00718 Pa**, maximum friction
velocity **0.0747 m/s**, and wall impulse **(-4.736, 0.342, -8.706) N s**.
Relative species-mass residual was `5.53e-16`, enthalpy residual `0`, maximum
continuity residual `5.71e-13 kg m^-3 s^-1`, and momentum-balance residual
`3.08e-9 kg m s^-1`; all run gates pass. The materialized prefix is
`outputs/thorney-ti021-zero-surface-wall-shear-updated-dx5m-prefix-20s-2026-10-12-v1/`.

The hash-pinned pre/post comparison is
`outputs/thorney-ti021-wall-shear-closure-20s-comparison-2026-10-12-v1.json`.
Over common cells, the peak-ever concentration RMSE is `6.93e-6 vol%-points`,
sensor trace RMSE `6.12e-6 vol%-points`, and 4-vol% maximum-envelope volume
does not change (both `14,875 m3`). The largest velocity-u difference is
`2.36e-4 m/s`. Thus this closure has a measurable but very small effect on this
coarse 0--20 s case; no claim of experimental accuracy follows, and roughness
or separated-wake validation remains absent.
Regional comparison locates the largest velocity response in the windward
inner-face shell (`u` RMSE `8.88e-5 m/s`, 54 cells) and lee-wake shell
(`5.16e-5 m/s`, 150 cells); side-tip bypass response is smaller
(`6.20e-6 m/s`, 24 cells). Regional concentration RMSE remains at most
`4.70e-5 vol%-points`. These sparse coarse-grid shells are only a sensitivity
diagnostic; they do not validate stagnation, separation or recirculation.

The independent AIJ Case H instantaneous neutral-wake audit is hash-pinned at
`outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_resolved_3d_grid_audit.json`.
It compares 176 common center-plane U/W points on three grids. The finest grid
has normalized U/W correlations of **0.9483/0.8943** and recirculation-sign
accuracy **0.9253**, but the quantitative gate remains withheld because the
finest-pair U error changes by 11.28%, W-RMSE changes by 12.08%, the lee-speed
change is 23.13%, and the flow residual misses its tolerance. This is useful
neutral-wake component evidence only; it is not H₂ concentration validation.

The newer obstacle-face-aligned, time-averaged 3-D AIJ audit is
`outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_time_averaged_aligned_grid_audit.json`.
It is `accepted_component_grid`: 176 common points, all three mean flows pass
the 1% statistical-stationarity tolerance, obstacle alignment error is
effectively zero, and finest-pair U-MAE / W-RMSE changes are 0.94% / 6.15%.
This is stronger evidence for the neutral mean-flow component than the
instantaneous audit, but it still withholds scalar concentration validation,
TKE validation, and cryogenic-H₂ prediction.

The first 3-D ethylene-tracer scoring pilot is recorded at
`outputs/aij-case-h-tracer-resolution-sensitivity-2026-10-12-v1/` and is
reproduced by `tools/score_aij_case_h_resolved_tracer.py` plus
`tools/compare_aij_case_h_tracer_resolution.py`. It maps the README's 4 mm
floor nozzle and 5.83×10⁻⁶ m³/s pure-ethylene source into an ethylene mass rate
at 293 K/1 atm, then scores model mass concentration as ideal-gas ethylene ppmv
without fitting to the sensor values. On the same 8–12 s interval and 582
common points, the 25 mm and 12.5 mm predictions have RMSE **2,624.5** and
**5,449.1 ppmv**, respectively; bias is **+745.6/+906.9 ppmv**, against a
measured common-point mean of **121.9 ppmv**. Both carrier-flow stationarity
gates pass and species-mass residuals are **5.41×10⁻¹⁴ / 2.14×10⁻¹³**. This is
an adverse scalar result, not an accepted concentration validation. It is not
a pure grid-convergence test because scalar dt and Gaussian source width scale
with dx. The D3Q19 passive-scalar runs also have zero vertical source-exit
momentum and no validated turbulence-scalar covariance; those source/wake
physics must be addressed before rescoring.

A controlled medium-grid spin-up diagnostic is recorded at
`outputs/aij-resolved-spinup-residual-2026-10-12.json` and reproduced by
`tools/diagnose_aij_case_h_spinup_residual.py`. Extending the same nominal grid
from 1,500 to 4,500 D3Q19 steps (0.836 s physical flow time) reduced the
instantaneous relative velocity-update residual from 0.002454 to 0.001716;
over the last 2,500 steps it changed only 1.30%. The projected continuity
divergence remained 9.66×10⁻¹⁰ s⁻¹, while the solver's 10⁻⁶ instantaneous
velocity-update threshold remained unmet. This residual is not the continuity
error. The separate time-averaged aligned-grid audit is the relevant mean-flow
comparison and now passes its declared component-grid gates. Do not spend
runtime on unbounded instantaneous spin-up alone as a presumed fix; concentration
and cryogenic-H₂ validation remain separate unmet requirements.

The complete-horizon resolved screening member is
`outputs/thorney-ti021-coarse-600s-dx2p5m-dt0p0625-checkpoint1s-tol1e3-2026-10-11/`.
It covers 0--600 s on a 2.5 m grid, has 60 receptor traces, passes its numerical
gates, and reports 7,120.52 s wall-clock (6,783.20 s solver-reported). Its
full-horizon score is retained as a screening-tier result; it is not treated as
finest-grid validation because the matched 1.25 m run currently covers only the
0--20 s impact window.

Its review-ready front/top/side/lee-wake package is
`outputs/thorney-ti021-600s-screening-package-2026-10-12-v1/`. The package
reports all concentration values in vol%, the 4-vol% maximum-ever volume,
60-receptor traces, and the fixed obstacle-region metrics. The fourth panel is
masked to the exterior lee-wake shell rather than being inferred from the plan
view alone.

The linked post-run TI021 score records dose MAE **265.93 vol% s**, dose MAPE
**86.54%**, Pearson **r=0.8474**, Spearman **ρ=0.9141**, arrival-time MAE
**22.55 s**, departure-time MAE **180.55 s**, and six-arc peak/dose MAPEs of
**52.09%/79.48%**. These are screening-tier comparison metrics; they are not
finest-grid promotion evidence.

The screening package also exposes pressure, velocity, density and temperature
state fields plus conservation residuals. In this particular 600-s TI021
screening run, temperature stayed at 291.15 K and the relative density span was
only **3.04×10⁻¹²**; therefore the low-temperature density/buoyancy transition is
implemented in the solver but explicitly marked **not exercised in this run**.

The separate active-thermal development witness is
`outputs/obstacle-physics-rebuild-2026-10-10/obstacle_thermal_buoyancy_comparison.json`.
With an 80 K source it reaches **21.14 vol%**, minimum **274.49 K**, and
**5.21%** maximum mixture-density departure; maximum mass and energy residuals
are **6.07×10⁻¹⁸ kg** and **6.63×10⁻¹⁶**. This demonstrates the sign-changing
density/buoyancy path, but remains a development witness rather than facility
validation because its boundaries and active-scalar closure are not externally
qualified.

The completed fine impact-window run is:

`outputs/thorney-ti021-three-grid-impact-dx1p25m-20s-dt0p015625-checkpoint1s-tol1e4-2026-10-11-v3`

The earlier mixed-time-step audit is retained as a historical artifact; the
spatially isolating, same-`dt` three-grid audit is:

`outputs/thorney-ti021-matched-dt-grid-convergence-2026-10-12-v3/ti021-three-grid-impact-convergence.json`

The newly run 5 m and 2.5 m members cover 0--20 s at `dt=0.02 s`, matching the
existing 1.25 m member. They completed in **95.65 s** and **862.01 s** wall time,
respectively, and all per-run numerical gates passed. Their record, concentration
field and receptor-trace hashes are pinned in the three-grid audit. The audit
isolates spatial sensitivity but still rejects hazard-volume and receptor
convergence; see the corrected convergence interpretation below.

The additional stable-time-step run has completed:

`outputs/thorney-ti021-three-grid-impact-dx1p25m-20s-dt0p02-checkpoint1s-tol1e4-2026-10-12-v1`

It uses `dx=1.25 m`, `dt=0.02 s`, 20 s duration and 1 s checkpoints. All 20
segments pass the conservation and stability diagnostics with no stderr output.
The authoritative runtime is 22,620.93 s wall-clock (22,517.76 s solver-reported).
The final field maximum is 32.1019 vol%, and the final 4-vol% volume is 654.30 m³
(not a converged or full-event hazard result).

A separate zero-thickness-surface sensitivity run at the same 1.25 m grid and
0.02 s step has now committed its first 0--1 s segment at
`outputs/thorney-ti021-zero-surface-dx1p25m-20s-dt0p02-2026-10-12-v1/`.
The segment and restart hashes are recorded in its checkpoint journal; the
materialized prefix at
`outputs/thorney-ti021-zero-surface-prefix-1s-dt0p02-2026-10-12-v1/` passes
the computational prefix gates (complete 51-frame time axis, 60 receptors,
mass/enthalpy conservation, continuity and corrector checks). The maximum is
32.0709 vol%; the final-field and any-time 4 vol% volumes are 3,220.70 m³ and
3,572.27 m³. The three-view field is
`ti021-zero-surface-prefix-1s-contours.png` in that prefix directory. These are
one-second, computed-window sensitivity outputs—not observation agreement,
converged hazard volumes, or full-event results. At the latest inspection the
20 s run remained active; only the 1 s checkpoint was committed.

The matched 0--1 s wall-representation comparison is recorded at
`outputs/thorney-ti021-geometry-representation-prefix-comparison-2026-10-12.json`
and reproduced by
`tools/compare_thorney_ti021_geometry_representation_prefixes.py`. The voxel
wall occupies 536 solid cells while the zero-thickness surface occupies none;
over 818,664 shared fluid cells, endpoint concentration RMSE is only 0.00206
vol%-points and pooled 60-sensor trace RMSE is 8.60×10⁻¹¹ vol%-points. No
receptor reached 4 vol% by 1 s. Endpoint 4-vol% volumes are 3,216.80/3,220.70
m³ (voxel/surface), and historical-envelope volumes are 3,562.50/3,572.27 m³.
This near-equality is an early-time non-arrival result; it does not establish
equivalent obstacle-wake transport. The later committed prefixes are needed to
compare windward impingement, roof/side bypass and lee-wake behavior.

The current checkpoint runner reconstructs the steady resolved momentum field at
the start of every checkpointed segment before restoring the scalar restart state.
This is conservative and restart-safe, but it adds repeated flow-solve overhead to
short checkpoint intervals. The measured 1-s rate should not be treated as the
intrinsic cost of a single uninterrupted transient. A production execution-time
audit must separately report checkpoint interval and, after a flow-state restart
extension is implemented and tested, the uninterrupted-flow rate.

An independently materialized, immutable 5-s prefix is available at
`outputs/thorney-ti021-checkpoint-prefix-5s-dt0p02-2026-10-12/`. Its record is
`ti021-inner-transient.json`; all prefix gates pass, including the 60-receptor
trace axis, species-mass and sensible-enthalpy conservation, continuity residual,
and zero non-converged corrector steps. This is a qualified computed-window
artifact only and does not promote the full event or hazard-volume convergence.
The independently regenerated three-view figure is
`outputs/thorney-ti021-checkpoint-prefix-5s-dt0p02-2026-10-12/ti021-prefix-5s-final-contours.png`.

The completed run also has immutable prefixes through 6 s. Its immutable record and field
are in `outputs/thorney-ti021-checkpoint-prefix-6s-dt0p02-2026-10-12/`, with
the three-view figure at
`outputs/thorney-ti021-checkpoint-prefix-6s-dt0p02-2026-10-12/ti021-prefix-6s-final-contours.png`.
The 6-s checkpoint passes all numerical gates; the historical maximum remains
32.1019 vol%, and the 4-vol% historical geometry is 11,212.89 m³ at that
checkpoint. These are computed-window values only and remain withheld from
full-event promotion.

The reproducible 5--6 s evolution summary is
`outputs/thorney-ti021-prefix-evolution-2026-10-12/ti021-5s-to-6s.json`.
Over that additional second, the final-field 4-vol% volume increased from
5,048.83 to 5,751.95 m³, while the historical maximum-envelope volume increased
from 9,849.61 to 11,212.89 m³. The 60-receptor label hash is unchanged, so the
sensor axis is identical across prefixes. This is a computed-window evolution
diagnostic, not a validation or promotion result.

The zero-thickness-surface 20-s sensitivity run has now committed through 18 s.
The immutable materialization is
`outputs/thorney-ti021-zero-surface-prefix-18s-dt0p02-2026-10-12-v2/` and passes
all prefix gates (901 frames, 60 receptors, mass/enthalpy conservation,
continuity and corrector convergence). Its time-maximum three-view field is
`ti021-zero-surface-prefix-18s-contours.png`. The 4-vol% endpoint volume is
626.95 m³, while the any-time envelope is 15,431.64 m³; the receptor peak is
23.27 vol% and endpoint maximum is 3.27 vol%. These are computed-window values,
not independent validation or a converged full-event hazard estimate. The
reproducible same-geometry evolution from 2 to 18 s is
`outputs/thorney-ti021-zero-surface-prefix-evolution-2s-to18s-2026-10-12-v2/ti021-2s-to18s.json`.
The source solver remains live after this committed-prefix extraction.
At the same stop time, a no-fit voxel-wall vs zero-thickness comparison is
recorded in `outputs/thorney-ti021-geometry-representation-prefix-18s-2026-10-12-v2.json`.
On 818,664 shared fluid cells its endpoint field RMSE is **0.1633 vol%-points**;
the endpoint 4-vol% volume is **1,525.39 m³** (voxel) versus **626.95 m³**
(zero-thickness), while any-time envelopes are **14,724.61** and **15,431.64 m³**.
This shows material endpoint-hazard sensitivity to wall representation, not which
representation is more accurate; the current model horizon is truncated.
The 19-s prefix is separately hash-pinned at
`outputs/thorney-ti021-zero-surface-prefix-19s-dt0p02-2026-10-12-v1/`; its
matched wall-representation comparison is
`outputs/thorney-ti021-geometry-representation-prefix-19s-2026-10-12-v1.json`.
At 19 s, common-fluid endpoint RMSE is **0.1607 vol%-points**. The endpoint
4-vol% volume is **783.20 m³** (voxel) vs **205.08 m³** (zero-thickness),
while any-time volumes remain **14,728.52** vs **15,431.64 m³**. This
computed-window gap persists but does not identify the more accurate geometry.
The zero-surface run has now finalized and its full 20-s checkpoint was
independently materialized at
`outputs/thorney-ti021-zero-surface-prefix-20s-dt0p02-2026-10-12-v1/`.
All materialization gates pass: 1,001 frames, 60 receptors, relative species
mass residual `4.15e-16`, enthalpy residual `0`, maximum continuity residual
`3.67e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The 20-s
field figure is `ti021-zero-surface-prefix-20s-contours.png`; its maximum is
`32.10 vol%`. The matched 20-s voxel-vs-surface comparison is
`outputs/thorney-ti021-geometry-representation-prefix-20s-2026-10-12-v1.json`:
common-fluid endpoint RMSE `0.1575 vol%-points`, voxel/surface 4-vol% endpoint
volumes `654.30/167.97 m3`, and historical envelopes `14,730.47/15,431.64 m3`.
This is boundary-representation sensitivity, not evidence that either geometry
is more accurate.

The blind score against TI021 is at
`outputs/thorney-ti021-zero-surface-prefix-20s-score-2026-10-12-v1/` and is
explicitly `withheld_as_horizon_limited`: only 23 of 53 positive sensors cross
the 0.1-vol% detection threshold; normalized dose bias is `-0.955`, dose MAPE
`0.978`, arrival MAE `26.93 s`, and six-arc peak MAPE `0.621`. Observed arc
peaks are at `86.5--133.5 s`, so this 20-s score is not final accuracy. All six
observed arc peaks are below 4 vol%; that is true-negative behavior only, not
a 4-vol% sensitivity test. The zero-surface run's recorded wall time is
`21,321.90 s` (5.92 h); the separate voxel-wall run's `22,620.93 s` timing
below must not be attributed to it.
An observation-blind score against TI021 was also run on this 18-s prefix at
`outputs/thorney-ti021-zero-surface-prefix-18s-score-2026-10-12-v1/`. It is
explicitly `withheld_as_horizon_limited`: 53 sensors have positive observations,
but only 21 cross the model's 0.1 vol% detection threshold and 31 observed
arrivals are missed. Observed six-arc peaks occur at 86.5--133.5 s, well beyond
the 18-s model horizon. Its apparent dose bias/MAPE therefore describes a
truncated, non-comparable window—not final model accuracy. This confirms that
finishing the remaining 2 s of the current run cannot close TI021 validation;
the required next benchmark run must span the observation event horizon.

For context, the existing full-horizon 600-s / 2.5-m run was scored separately
in `outputs/thorney-ti021-inner-transient-score-dx2p5m-2026-10-11-v1/`. It
covers TI021 time but remains withheld from quantitative validation: 42/53
positive sensors are detected, 10/52 arrival times are missed, normalized dose
bias is **-0.888**, dose MAPE is **0.865**, and arrival MAE is **22.56 s**.
Six-arc peaks are predicted at **0.52--0.72 vol%** versus observed
**1.16--1.47 vol%** (peak MAPE **0.519**). The grid aliases low-height
receptors; all observed six-arc peaks are below 4 vol%, so this does not test
4-vol% sensitivity. Full event coverage alone is insufficient for promotion.

The active run subsequently committed a 7-s prefix. The immutable record, field,
receptor traces, and three-view contour figure are in
`outputs/thorney-ti021-checkpoint-prefix-7s-dt0p02-2026-10-12-v3/`.
The 6--7 s evolution report is
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v2/ti021-6s-to-7s.json`.
The final-field 4-vol% volume increased from **5,751.95 to 6,335.94 m³** and
the historical maximum-envelope volume from **11,212.89 to 12,355.47 m³**.
The historical peak remained **32.1019 vol%**; the endpoint maximum decreased
from **10.4329 to 8.4941 vol%** as the cloud advected through the receptor grid.
All prefix numerical gates pass, but the artifact remains a computed-window
evolution only: full-event promotion and independent-field validation remain
withheld.

The active run has now committed an 8-s prefix at
`outputs/thorney-ti021-checkpoint-prefix-8s-dt0p02-2026-10-12-v1/`, including
the regenerated contour figure. The reproducible 7--8 s evolution report is
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v3/ti021-7s-to-8s.json`.
The final-field 4-vol% volume changed from **6,335.94 to 6,490.23 m³**, while
the historical maximum-envelope volume changed from **12,355.47 to
13,261.72 m³**. The historical peak stayed at **32.1019 vol%**; the endpoint
maximum decreased from **8.4941 to 6.6510 vol%**. All prefix numerical gates
pass, but this remains a computed-window evolution and is not an independent
field-validation or full-event promotion result.

The live runtime projection, based on the latest positive 1-s progress interval,
is recorded at
`outputs/runtimeaudit-2026-10-12-v2/ti021-runtime-live.json`. At 5 s it estimates
7.30 h wall time for the 20-s checkpointed run. This is explicitly a live
projection, not a completion result; it includes repeated flow reconstruction at
each checkpoint and must not be generalized to an uninterrupted production run.

After the 7-s checkpoint, the refreshed projection is
`outputs/runtimeaudit-live-2026-10-12-v5/runtime-v5.json`. The latest interval
was 1 model second per 1,219.78 wall seconds, giving a provisional **7.18 h**
for the 20-s run. This remains a live projection only and will be superseded by
the next positive checkpoint interval.

The next live projection, after the 8-s checkpoint, is
`outputs/runtimeaudit-live-2026-10-12-v6/runtime-v6.json`. Its latest interval
is 1 model second per **1,119.44 wall seconds**, giving a provisional **6.82 h**
for the 20-s run.

The same 0--5 s window was compared at `dt=0.02 s` versus `dt=0.015625 s` using
`tools/compare_thorney_ti021_temporal_prefix.py`. The reproducible report is
`outputs/thorney-ti021-temporal-prefix-comparison-2026-10-12/ti021-temporal-prefix-comparison.json`.
The field MAE/RMSE are 0.000399/0.007383 vol%-points, and the interpolated
60-receptor time-series MAE/RMSE are 0.002520/0.016003 vol%-points. The 4-vol%
volumes differ by 9.77 m³ in this computed window. These results support local
time-step sensitivity only; they do not establish full-event or independent-field
validation.

The completed 0--20 s, 1.25-m runs were also compared at `dt=0.02 s` and
`dt=0.015625 s` using the same temporal-prefix comparator. The hash-pinned
report is `outputs/thorney-ti021-temporal-prefix-comparison-20s-2026-10-12-v2/ti021-temporal-prefix-comparison-20s.json`.
Across 818,664 maximum-field cells the RMSE is **0.00618 vol%-points**; across
76,860 interpolated receptor samples the RMSE is **0.00876 vol%-points**. The
two 4-vol% maximum-envelope volumes are **14,724.61** and **14,730.47 m³**
(difference **5.86 m³**, about **0.040%** of the smaller volume). This is
evidence of timestep insensitivity for this 20-s computed window only; it does
not repair the non-monotonic spatial hazard-volume GCI or establish full-event
sensor/hazard convergence. The provenance record contains SHA-256 digests for
both input records, fields and receptor traces. The updated comparator and the
temporal/grid-audit focused regression tests pass (**10 tests**).

The temporal comparison tool and obstacle-resolved validation path pass the focused
regression suite: 55 tests passed in 64.08 s, including temporal-prefix comparison,
checkpoint equivalence, stratification diagnostics, grid-convergence audits, and
the 3-D resolved obstacle solver tests.

The fine-grid stratification diagnostic is reproducible at
`outputs/thorney-ti021-stratification-fine-2026-10-12-v3/ti021-stratification-diagnostic.json`
with a companion figure in the same directory. During this audit, resolved solid
cells were found to carry NaN density sentinels. The diagnostic now validates
finite/positive values only in fluid cells and sanitizes solid entries before
finite-volume face arithmetic; the solid-cell result remains zero and excluded
from summaries. The targeted stratification and audit tests pass (`10 passed`).
This diagnostic reports local buoyancy/shear structure only and does not modify
transport coefficients or promote a turbulence closure.

At the 20 s endpoint, the region means/p95 of gradient Richardson number were:
windward face **1.378 / 6.028**, roof shear band **0.353 / 0.914**, side-tip
bypass **0.072 / 0.570**, and lee-wake shell **0.0146 / 0.136**. These values
quantify the expected stable near-face stratification and weakly stratified wake
without treating the diagnostic as concentration validation.

The real-event evidence join is hash-pinned in
`outputs/IJHE_UPLOAD_AUTHOR_APPROVED_2026-10-10_DEGALI_DYNAMIC_V11/PRESLHY_Trial10_field_evidence_manifest.json`.
Its verification record is
`outputs/preslhy-trial10-manifest-verification-2026-10-12.json`: source boundary,
weather, obstacle geometry, receptor observations, common clock, and the sensor
registry all verify by SHA-256. The package remains **conditional** because the
event-level sensor calibration certificate is not present; `promotion_allowed` is
therefore `false`. This manifest is evidence provenance, not a claim that the
5-second TI021 computed prefix has been independently validated.

The 9-s checkpoint was then atomically committed from the live run at
`outputs/thorney-ti021-checkpoint-prefix-9s-dt0p02-2026-10-12-v3/`. Its record
passes all numerical gates: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The field
maximum over the prefix remains **32.1019 vol%**. A reproducible contour figure
is `outputs/thorney-ti021-checkpoint-prefix-9s-dt0p02-2026-10-12-v3/ti021-prefix-9s-final-contours.png`.

The 8-to-9-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v4/ti021-8s-to-9s.json`.
The final 4-vol% volume decreased from **6,490.23 to 6,150.39 m^3** while the
historical envelope increased from **13,261.72 to 13,781.25 m^3**; the endpoint
maximum changed from **6.6510 to 5.9744 vol%** and the endpoint mean from
**0.7638 to 0.6657 vol%**. Receptor-label hashing is unchanged. This is still
computed-window evolution only (`promotion_allowed=false`), not independent
field validation or full-event promotion.

The refreshed live runtime projection is
`outputs/runtimeaudit-live-2026-10-12-v7/runtime-v7.json`. Using the latest
positive 1-s interval (1,115.34 wall seconds per model second), the projected
20-s wall time is **6.81 h**. This is a live projection, not a completion
measurement, and includes checkpoint flow reconstruction.

The live runtime projection was refreshed after the 9-s progress record at
`outputs/runtimeaudit-live-2026-10-12-v8/runtime-v8.json`. The latest positive
1-s interval remains **1,115.34 wall seconds per model second**, giving a
projected 20-s wall time of **6.81 h**. This is still a live projection rather
than a completion measurement; the solver remains at the 9-s checkpoint while
the 9--10-s segment is being computed.

The 10-s checkpoint was then atomically committed at
`outputs/thorney-ti021-checkpoint-prefix-10s-dt0p02-2026-10-12-v1/`. All
numerical gates pass again: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**. Its final contour figure is
`outputs/thorney-ti021-checkpoint-prefix-10s-dt0p02-2026-10-12-v1/ti021-prefix-10s-final-contours.png`.

The 9-to-10-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v5/ti021-9s-to-10s.json`.
The final 4-vol% volume decreased from **6,150.39 to 5,388.67 m^3** while the
historical envelope increased from **13,781.25 to 14,126.95 m^3**; endpoint
maximum changed from **5.9744 to 5.4595 vol%** and endpoint mean from
**0.6657 to 0.5991 vol%**. Receptor-label hashing is unchanged. This remains
computed-window evolution only (`promotion_allowed=false`).

The refreshed live runtime projection is
`outputs/runtimeaudit-live-2026-10-12-v9/runtime-v9.json`. The latest positive
1-s interval is **1,128.49 wall seconds per model second**, projecting the
20-s run to **6.85 h**. This is a live projection, not a completion result.

The 11-s checkpoint was then materialized at
`outputs/thorney-ti021-checkpoint-prefix-11s-dt0p02-2026-10-12-v1/`. All
numerical gates pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**. Its final contour figure is
`outputs/thorney-ti021-checkpoint-prefix-11s-dt0p02-2026-10-12-v1/ti021-prefix-11s-final-contours.png`.

The 10-to-11-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v7/ti021-10s-to-11s.json`.
The final 4-vol% volume decreased from **5,388.67 to 4,984.38 m^3** while the
historical envelope increased from **14,126.95 to 14,339.84 m^3**; endpoint
maximum changed from **5.4595 to 5.0522 vol%** and endpoint mean from
**0.5991 to 0.5529 vol%**. Receptor-label hashing is unchanged. This remains
computed-window evolution only (`promotion_allowed=false`).

The refreshed live runtime projection is
`outputs/runtimeaudit-live-2026-10-12-v11/runtime-v11.json`. The latest positive
1-s interval is **1,129.66 wall seconds per model second**, projecting the
20-s run to **6.85 h**. This is a live projection, not a completion result.

The 12-s checkpoint was then materialized at
`outputs/thorney-ti021-checkpoint-prefix-12s-dt0p02-2026-10-12-v1/`. All
numerical gates still pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**. Its final contour figure is
`outputs/thorney-ti021-checkpoint-prefix-12s-dt0p02-2026-10-12-v1/ti021-prefix-12s-final-contours.png`.

The 11-to-12-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v8/ti021-11s-to-12s.json`.
The final 4-vol% volume decreased from **4,984.38 to 4,664.06 m^3** while the
historical envelope increased from **14,339.84 to 14,476.56 m^3**; endpoint
maximum changed from **5.0522 to 4.7153 vol%** and endpoint mean from
**0.5529 to 0.5200 vol%**. The final-field p99 changed from **1.1809 to
1.1876 vol%**, and the historical peak remains **23.2845 vol%** at the
receptors. Receptor-label hashing is unchanged. This remains computed-window
evolution only (`promotion_allowed=false`).

The refreshed live runtime projection is
`outputs/runtimeaudit-live-2026-10-12-v12/runtime-v12.json`. The latest positive
1-s interval is **1,173.17 wall seconds per model second**, projecting the
20-s run to **6.96 h**. This is a live projection, not a completion result.

The 13-s checkpoint was then materialized at
`outputs/thorney-ti021-checkpoint-prefix-13s-dt0p02-2026-10-12-v1/`. All
numerical gates pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**. Its final contour figure is
`outputs/thorney-ti021-checkpoint-prefix-13s-dt0p02-2026-10-12-v1/ti021-prefix-13s-final-contours.png`.

The 12-to-13-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v9/ti021-12s-to-13s.json`.
The final 4-vol% volume increased from **14,476.56 to 14,576.17 m^3** in this
comparison output, while the reported historical envelope was not used as a
promotion criterion. The field MAE is **0.00230 vol%**, RMSE **0.03080 vol%**,
and maximum absolute field difference **2.7380 vol%**; receptor comparison is
zero over the common 0--12-s frame because the 13-s trace is interpolated to
the shorter reference. This is computed-window evolution only
(`promotion_allowed=false`), not independent validation or full-event
promotion.

The 14-s checkpoint was subsequently materialized at
`outputs/thorney-ti021-checkpoint-prefix-14s-dt0p02-2026-10-12-v1/`. The same
numerical gates pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**.

The 13-to-14-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-2026-10-12-v10/ti021-13s-to-14s.json`.
The final 4-vol% volume increased from **14,576.17 to 14,638.67 m^3**. The
field MAE is **0.00207 vol%**, RMSE **0.02413 vol%**, and maximum absolute
field difference **2.2682 vol%**; receptor comparison is zero over the common
0--13-s frame because the 14-s trace is interpolated to the shorter reference.
This is still computed-window evolution only (`promotion_allowed=false`), not
independent validation or full-event promotion.

The 17-s checkpoint was independently materialized at
`outputs/thorney-ti021-checkpoint-prefix-17s-dt0p02-2026-10-12-v1/`. All
numerical gates pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%**. Its final 4-vol% volume is **2,273.44 m^3**
and its historical maximum-envelope volume is **14,718.75 m^3**.

The 14-to-17-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-14-17s-2026-10-12-v1/summary.json`.
The final 4-vol% volume decreased by **1,519.53 m^3** while the historical
maximum-envelope volume increased by **80.08 m^3**; the endpoint maximum changed
from **4.1685 to 3.5640 vol%**. These are computed-window diagnostics only and
do not establish independent field validation or full-event promotion.

The 18-s checkpoint was independently materialized at
`outputs/thorney-ti021-checkpoint-prefix-18s-dt0p02-2026-10-12-v1/`. All
numerical gates pass again: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%** and the final 4-vol% volume is **1,525.39 m^3**.

The 17-to-18-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-17-18s-2026-10-12-v1/summary.json`.
The final 4-vol% volume decreased by **748.05 m^3**, while the historical
maximum-envelope volume increased by **5.86 m^3**; endpoint maximum changed
from **3.5640 to 3.4057 vol%**. This remains computed-window evidence only
(`promotion_allowed=false`).

The live runtime projection was refreshed after the 18-s progress record at
`outputs/runtimeaudit-live-2026-10-12-v14/runtime-v14.json`. The latest
positive interval is **768.62 wall seconds per model second**, giving a
provisional full-run projection of **6.24 h**. This is a live projection only,
not a completion measurement.

The 19-s checkpoint was independently materialized at
`outputs/thorney-ti021-checkpoint-prefix-19s-dt0p02-2026-10-12-v1/`. All
numerical gates still pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum remains **32.1019 vol%** and the final 4-vol% volume is **783.20 m^3**.

The 18-to-19-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-18-19s-2026-10-12-v2/summary.json`.
The final 4-vol% volume decreased by **742.19 m^3**, while the historical
maximum-envelope volume increased by **3.91 m^3**; endpoint maximum changed
from **3.4057 to 3.2649 vol%**. This remains computed-window evidence only
(`promotion_allowed=false`).

The full 20-s run completed and was independently materialized at
`outputs/thorney-ti021-checkpoint-prefix-20s-dt0p02-2026-10-12-v1/`. All
numerical gates pass: maximum relative species-mass residual
`4.1465e-16`, sensible-enthalpy residual `0`, maximum continuity residual
`4.6221e-11 kg m^-3 s^-1`, and zero non-converged corrector steps. The prefix
maximum is **32.1019 vol%** and the final 4-vol% volume is **654.30 m^3**;
the historical maximum-envelope volume is **14,730.47 m^3**.

The 19-to-20-s computed-window evolution is recorded at
`outputs/thorney-ti021-prefix-evolution-19-20s-2026-10-12-v1/summary.json`.
The final 4-vol% volume decreased by **128.91 m^3**, while the historical
maximum-envelope volume increased by **1.95 m^3**; endpoint maximum changed
from **3.2649 to 3.1389 vol%**. This remains computed-window evidence only
(`promotion_allowed=false`).

The separate voxel-wall run record reports **22,620.93 s wall-clock** and
**22,517.76 s solver-reported runtime** for the 20-s, 1.25-m, dt=0.02-s
resolved obstacle case. The run produced the HDF5 handoff, 60-receptor
time-series, final field, checkpoint history, and progress log. The runtime
is recorded for reproducibility; it is not an industrial screening runtime.

## Convergence audit interpretation

The matched-time-step three-grid audit is
`outputs/thorney-ti021-matched-dt-grid-convergence-2026-10-12-v4/ti021-three-grid-impact-convergence.json`.
All three runs use `dt=0.02 s`, identical 20-s duration and 2:1 spatial
refinement, with the same workbook and property-audit digests. Final and maximum
field L1/RMSE differences contract on common fluid cells, but this does **not**
qualify fixed-geometry spatial convergence: the wall thickness is set to `dx`
(5 / 2.5 / 1.25 m), so physical geometry changes along the refinement path.
The 4-vol% hazard volumes
(14,625.00 / 14,781.25 / 14,730.47 m³) are non-monotonic. Receptor convergence
is worse when GCI uses the physically correct spatial refinement ratios:
peak/dose monotonic fractions are 0.217/0.208, with median fine GCI
52.7%/92.2%. Therefore
`full_event_grid_convergence_promoted=false`; the fine-grid field is a qualified
computed-window result, while hazard geometry and sensor metrics remain withheld.
The completed fixed zero-thickness numerical-surface 0--20 s audit is
`outputs/thorney-ti021-zero-surface-matched-dt-grid-convergence-20s-2026-10-12-v2/`;
v2 supersedes the metrics-identical v1 with an explicit record of the published
tarpaulin-on-scaffold construction and the model's scaffold omission.
Matched `dt=0.02 s` runs at 5/2.5/1.25 m contract field L1/RMSE/L∞ norms
(final-field RMSE 0.1106 then 0.0878 vol%-points on the 5-m comparison grid),
but the 4-vol% volume sequence 14,875/15,000/15,431.64 m³ has no qualifying
monotonic GCI. Only 6/25 eligible receptor peaks and 9/25 doses are monotonic
convergent (0.24/0.36); median fine-grid GCI is 56.1% for peaks and 244.0%
for dose. The idealized zero-thickness numerical boundary is held fixed across
grids, so this is a valid convergence diagnostic for that numerical idealization;
however, the predeclared hazard-volume and receptor gates fail. Separately, the
published Phase II description identifies the real barrier as tarpaulin sheets
hung over a scaffold frame, arranged on a 180-degree arc 50 m from the source.
The tarpaulin gauge and scaffold member sections are not specified in the source
available here, and the numerical model omits the scaffold. Thus the calculation
is not yet a physical-geometry qualification of the facility fence. See Davies
and Singh (1985), [DOI 10.1016/0304-3894(85)85043-3](https://doi.org/10.1016/0304-3894(85)85043-3).
The available component time audits similarly accept only the explicitly tested
component cases; no broad full-event time-convergence claim is made.

This matched-dt audit also corrects a diagnostic defect in the earlier
`thorney-ti021-three-grid-impact-convergence-2026-10-12-v1` report: that report
used temporal refinement ratios (2 and 4) inside the *spatial* receptor GCI
calculation. Its reported receptor GCI numbers are superseded. The audit now
accepts equal time steps to isolate spatial effects, requires time steps not to
get coarser with grid refinement, and always computes receptor GCI with the
2:1 spatial ratios. Both the matched-dt runs and corrected audit still fail the
predeclared hazard/receptor convergence promotion gates; no favorable claim is
being inferred from the method correction.

## Fixed-surface short-prefix diagnostic

Existing zero-thickness-surface runs at 5 m and 2.5 m cover 20 s. The 1.25 m
zero-surface run has since advanced beyond the earlier 3 s prefix: the
authoritative live process (PID 25588) has a committed checkpoint through 19 s
and progress through 19.8 s (99%); a 20 s checkpoint has not yet been committed.
Do not treat the separate voxel-wall 20 s package above as completion of this
zero-surface run. To make an honest like-for-like diagnostic without relabeling
it as a 20 s result, the committed 0--3 s prefixes were materialized on all
three grids and audited at matched `dt=0.02 s`. The updated materializer now
preserves the obstacle-boundary representation from the checkpoint journal;
earlier materialized records had omitted it, which made downstream audits
silently assume the finite-voxel default. The new audit is
`outputs/thorney-ti021-zero-surface-matched-3s-convergence-2026-10-12-v2/`.
The fixed numerical surface and field L1/RMSE contraction gates pass, but the
4-vol% volume GCI is 52.0%, only 33.3% of eligible receptor peaks are
monotonic-convergent, and the dose gate has 0% monotonic-convergent receptors.
The 20 s gate is explicitly false; this is diagnostic evidence only. The
1.25 m zero-surface run is continuing from committed checkpoints, not restarted
from zero; its terminal 20 s result remains pending.

The earlier v3 audit recorded voxelized semi-circular arc sampling: at 5 / 2.5 /
1.25 m the cell-band estimate was 140 / 150 / 167.5 m versus the analytical
157.08 m, with non-monotonic relative errors of 10.87% / 4.51% / 6.63%. That
one-cell voxel-wall path remains a geometry-regularisation sensitivity and is
not a measured fence thickness. A later matched-dt 0--20 s audit (v2) holds
the ideal zero-thickness immersed surface fixed across all three transport
grids. Its field norms contract, but hazard-volume GCI and receptor peak/dose
convergence still fail; the numerical-geometry comparison is therefore valid
for the stated idealization but not numerically converged by the predeclared
engineering gates.

The physical facility fence is better described than the older audit implied:
the Phase II paper says it consisted of tarpaulin sheets hung over a scaffold
frame, on a 180-degree arc 50 m from the source. It does not give tarp gauge or
scaffold member sections. The resolved model's zero-thickness impermeable
surface is a deliberate thin-sheet idealization and omits the frame; it is not
the fully registered physical assembly. The next steps are to improve the
hazard/receptor convergence results for the idealized surface and, separately,
obtain or justify the omitted frame/tarp details before facility-level
physical interpretation. The workbook porosity field is `1` for the described
impermeable fence, while report and workbook wind/height values also conflict;
these remain explicit rather than silently changing the geometry.

Exact segment/face geometry for a zero-thickness circular arc is implemented in
`src/degali/addons/site_geometry.py`, including staggered x/y blocked-face
masks. Those masks are now wired into the finite-volume pressure projection,
variable-density projection, and resolved/fast H2 and enthalpy transport paths;
the binary-mixture corrector also receives them. The focused regression batch
for geometry, transport, projections, thermal/binary correctors, and fast and
resolved obstacle flow passed before the D3Q19 surface-mode integration
(**98 tests**); the integrated focused batch passed (**113 tests**) before
the TI021 fixed-surface runner option was added. After that option, the
affected TI021-config/geometry/resolved-solver batch passed (**72 tests**).
The fast 3-D screening path now fails closed when a finite-thickness obstacle
occupies no cell centres, rather than silently omitting it; zero-thickness arc
walls remain represented by their embedded blocked-face treatment. The new
sub-grid crosswind-rod regression is included in the 113-test batch.
This verifies code paths and blocked-face/link invariants, not a field
prediction.

The resolved solver now has an explicit `circular_arc_representation` option.
Its default, `finite_thickness_voxel`, preserves prior behavior. The opt-in
`zero_thickness_surface` requires `CircularArcWall.thickness_m=0`, excludes
circular arcs from the solid-cell mask,
uses the Bouzidi-Firdaouss-Lallemand linear interpolated bounce-back with the
geometric wall-hit fraction on crossing D3Q19 links, and applies the exact
arc-crossed finite-volume face masks to projection and scalar/enthalpy fluxes.
The config rejects mismatched zero/finite thickness and representation modes.
If the upstream interpolation node is solid due to another voxel obstacle, it
falls back to halfway bounce-back for that link.

The TI021 configuration and checkpointed runner now expose the same explicit
mode, `--circular-arc-representation zero_thickness_surface`. In this mode the
report-priority fence stays at radius 50 m and height 5 m with zero declared
thickness on every grid, rather than silently setting thickness equal to `dx`.
A 0.25 s, 5 m-grid smoke run used 320 fractional D3Q19 wall links and closed
with relative species-mass residual `1.38×10⁻¹⁶` and continuity residual
`3.55×10⁻¹² kg/m³/s`; it is too short for plume/wake validation. The first
20 s fixed-surface impact-window run at 5 m completed all per-run numerical
gates in 146.35 s wall time. Its output is
`outputs/thorney-ti021-zero-surface-dx5m-20s-dt0p02-2026-10-12-v1/`.
The matched 2.5 m member has also completed all per-run numerical gates over
the same 0--20 s window in 1,155.27 s wall time; its species-mass residual is
`4.15×10⁻¹⁶` and continuity residual is `3.46×10⁻¹² kg/m³/s`. Its output is
`outputs/thorney-ti021-zero-surface-dx2p5m-20s-dt0p02-2026-10-12-v1/`.
The active 1.25 m, 20 s checkpointed member has now committed **2.0 s** in two
50-step segments. Segment wall times were **2,110.54 s** and **2,136.41 s**;
the second-segment rate extrapolates to about **10.7 h** for the remaining
18 simulated seconds. That remains a planning estimate, not a completed-run
runtime. Both committed segments pass the recorded mass, energy, continuity,
momentum and corrector gates; the live process remains active.
For context, the matched 2.5 m member completed the same 20 s window in
1,155.27 s (57.8 s per simulated second), so the first two 1.25 m segments
average about **36.7× slower per simulated second** than the coarser full-run
average. This is not a steady throughput ratio and flags a material runtime
risk that must be revisited with later segments and a completed fine-grid run
before claiming practical high-resolution use.
The checkpointed runner's progress-log default is now 0.1 simulated seconds
(five frames at `dt=0.02 s`) rather than 30 simulated seconds, so future
short-horizon runs expose intermediate progress at a useful cadence. This is a
runner-default change only; it does not alter the already-running process or
its recorded configuration.

The verified 0--2 s prefix is materialized at
`outputs/thorney-ti021-zero-surface-prefix-2s-dt0p02-2026-10-12-v1/` with all
prefix gates passing and 60 receptor traces. In this short computed window the
maximum historical concentration is **32.1019 vol%**, the endpoint field
maximum is **31.9927 vol%**, the endpoint 4 vol% volume is **3,605.47 m³**, and
the historical 4 vol% envelope is **5,021.48 m³**. The largest receptor value
reaches **13.1185 vol%** by 2 s. Relative to the 1 s prefix, the endpoint
hazard volume grows **384.77 m³**, the historical envelope grows **1,449.22
m³**, and the maximum receptor rises from effectively zero to 13.12 vol%.
These are computed-window dynamics, not a full-event hazard extent or
independent experimental validation. The three-view contour is
`outputs/thorney-ti021-zero-surface-prefix-2s-dt0p02-2026-10-12-v1/ti021-zero-surface-prefix-2s-contours.png`.
The companion
`outputs/thorney-ti021-zero-surface-prefix-2s-dt0p02-2026-10-12-v1/ti021-zero-surface-prefix-2s-endpoint-contours.png`
shows the endpoint field (the existing contour is the time-maximum field).
Both figures now overlay the configured idealized zero-thickness circular arc
and its front/side projections: the cell-centred solid mask is empty for this
representation, so a mask-only plot had misleadingly omitted the obstacle.
The plotting utility accepts `--field-view endpoint|time-maximum`; these are
two views of the same computed prefix, not independent evidence or validation.
The projection routine was also found to recompute connected components of the
same fixed fluid/obstacle graph on every time step. A bounded read-only cache
now reuses those labels by grid/solid/boundary geometry; 15 projection tests
cover fixed-geometry reuse and changed solid/blocked-face topologies, and the
resolved-obstacle suite passed in the preceding combined run. The
currently running PID predates this source change and therefore does not
exercise the optimization. The repeatable microbenchmark
`tools/benchmark_ti021_projection_component_cache.py` reconstructs the exact
TI021 1.25 m surface geometry (819,200 fluid cells; 352 blocked x-faces and
352 blocked y-faces for the idealized arc). The latest cold topology build was
6.77 s. The 0.46 μs mean cache-hit figure is the lookup with a prebuilt key;
the production-shaped path, including serializing masks to form cache keys and
looking up both projection geometry and labels, averaged **1.62 ms** over ten
repeats. This corrects the earlier comparison, which used only the prebuilt-key
lookup to describe the runtime path. Even the 1.62 ms figure is not a
projection solve or end-to-end benchmark. If the old topology traversal would
otherwise recur once per each of 50 steps, the cold traversal is about 15.8%
of the second segment's measured 2,136 s, while cache-key/lookup overhead is
negligible at that scale. The machine-readable record is
`outputs/ti021-projection-component-cache-benchmark-2026-10-12.json` and pins
SHA-256 for the property audit, benchmark script, projection source, and
obstacle-flow source. End-to-end benefit still requires a fresh matched run
after the live run reaches a safe terminal state.

Important geometry limits remain: masks block whole staggered faces and do not
provide fractional cut-cell areas. The D3Q19 arc rule follows the published
Bouzidi et al. curved-wall interpolation [DOI](https://doi.org/10.1063/1.1399290),
but this implementation and its coupling to DEGALI's forcing/projection path
have not yet been validated against a curved-wall flow benchmark. This mode
therefore remains a numerical treatment under development. A one-second
paired synthetic run after the interpolation change reports **0.479
vol%-points** endpoint field RMSE and maximum absolute difference **5.74
vol%-points** between voxel and ideal-surface geometries; this is a geometry
sensitivity result, not an accuracy score.
One-step smoke cases exercise both D3Q19 link bounce-back and finite-volume
projected-inlet surface masks; the D3Q19 case satisfies its declared global
species-mass residual gate. These are implementation checks, not flow-field
accuracy evidence.
The existing 5/2.5/1.25-m TI021 concentration and convergence artifacts were
not regenerated with this mode. Consequently the v4 fixed-geometry convergence
failure and all field-validation/promotion limitations above remain in force;
no updated concentration or convergence claim follows from these software
tests.

The reproducible paired implementation-sensitivity run is recorded at
`outputs/curved-arc-boundary-mode-comparison-2026-10-12-v1.json`, with a
common-scale plan/vertical/wake figure at the matching `.png` path. It compares
a 0.2-m voxel band to a true zero-thickness surface on a 0.1-m grid for a
synthetic 1-s, 0.012-kg release. Endpoint concentration RMSE is **0.479
vol%-points**, maximum absolute common-fluid difference is **5.74 vol%-points**,
and the historical 4-vol% volume differs by **0.023 m³**. Both runs' maximum
species-mass residual is **1.91×10⁻¹⁷ kg** (voxel) and **1.73×10⁻¹⁷ kg**
(surface). The record pins setup, runner and
figure SHA-256 values; it is explicitly labelled numerical sensitivity only,
not independent flow or dispersion validation.

The 3-D obstacle geometry also now includes `CrosswindCylinder` for finite
spanwise horizontal rods, with a circular cross-section mask in x-z. This
supports the shape family in the ExJet obstacle-array paper, but exact rod
count/centres and the direct release source remain unknown; no ExJet score has
been promoted using assumed geometry registration.

## Independent non-obstacle component evidence

The H2 tunnel resolved-transient replay at
`outputs/hytunnel-resolved-transient-score-2026-10-11-v2.json` provides an
independent transport-component check, not an obstacle-wake validation. On the
54-point published-box comparison, the qualified endpoint proxy gives
`MAE=2.1792 vol%`, `RMSE=2.4905 vol%`, `r=0.7334`, and a mean bias of
`+2.1697 vol%`; the 4-vol% classification has sensitivity 1.0 but specificity
0.5217. The score is deliberately not promoted because the workbook contains
steady-window averages, the common domain ends at 2.52 m, and no physical
obstacle is present. This evidence supports the scalar transport handoff only;
it cannot support the front/top/side/wake claim.

## FieldEvidenceManifest status

The manifest verifier has been run against the PRESLHY Trial 10 package at
`outputs/preslhy-trial10-manifest-verification-2026-10-12.json`. The five core
artifacts (source boundary, weather, obstacle geometry, receptor observations,
and common clock) and their SHA-256 digests verify, with 5,053 observation rows.
The result remains `conditional` and `promotion_allowed=false` because only
specification-level sensor calibration evidence is available; a calibration
certificate is still missing. The E3.4 obstacle-transport manifests have the
same limitation, plus report-figure-only, non-registered obstacle geometry and
unresolved atmospheric vapour flux. These are evidence-package gates, not
silently filled assumptions.

The external-data triage and import order are recorded in
`docs/public-obstacle-evidence-triage-2026-10-12.md`. The ExJet obstructed
free-jet study is the highest-priority new H2 obstacle benchmark; DataverseNO
adds public obstacle/source timing and pressure records, while H2SAFE is kept as
an explicitly labelled helium-surrogate sensor/facility test. The H2SAFE archive
itself is now SHA-256 verified locally (five events, STEP/STL geometry and sensor
CSV members), but its audit remains conditional because concentration units,
sensor calibration, release-clock alignment and some HVAC boundary mappings are
not declared.

## Newly imported ELVHYS WP4.2 component evidence

DataverseNO DOI `10.18710/JXJP0H` is a newly identified public CC0 dataset with
48 cryogenic-H2 transfer-connection-space tests. A minimal paired subset for
vent-restriction sensitivity (tests 29 and 31; nominally 2 barg and a 25.4 mm
nozzle) is SHA-256 pinned under the ignored local directory
`reference/hydrogen/elvhys_tcs_wp4_2/`. Both concentration histories contain 16
H2 channels; their sample intervals are 0.05 s (20 Hz, 18,764 records) and
0.002 s (500 Hz, 273,649 records). The README declares Xensor uncertainty of
±1 vol% over 0–100 vol% and five-point calibration, and the paired pressure
files include LH2 release-rig and nozzle pressure. The imported `FLMT` channels
are ventilation-fan flow, not hydrogen release mass flow. Test 29 uses four VAV
plates; test 31 uses three blank plates and one VAV plate. This is promising
real-H2 evidence for confined vent-restriction, pressure and sensor-history
components, but not an open-field obstacle-wake benchmark. Detailed registered
plate apertures and event-specific calibration certificates are not in the
subset. Visual review of the sensor-details drawing established a nominal
1 m³ enclosure, 3 mm steel walls, four rear vents, and the source/sensor
coordinate convention. It also exposed a geometry ambiguity: the README states
400 × 200 mm vents (0.08 m² by simple rectangular area), while the drawing has a
0.02 m² passive-vent callout without stating whether it is per vent or the
effective opening. The README separately gives pressure-peaking total open
areas of 0.0416/0.0208/0.0104 m² for tests 29/30/31. These are recorded in a
hash-pinned geometry registry, but the aperture shapes and exact opening
coordinates remain unresolved, so quantitative geometry stays withheld. The
four time vectors (concentration, pressure, ventilation flow and
ambient conditions) match exactly for both tests. One conditional
`FieldEvidenceManifest` is emitted per event, preserving the shared clock and
the sensor-calibration-certificate gate. No concentration score or promotion
is claimed. The repeatable audit is
`tools/audit_elvhys_tcs_component.py`, with machine-readable result
`outputs/elvhys-tcs-component-audit-2026-10-12.json`.

The experimental 1.25 m TI021 fine-grid transient has since reached its full
20 s window (20 committed segments; all 60 receptors present). The extracted
prefix audit reports mass/enthalpy conservation, continuity residual below
1e-7, converged correctors, and no use of observations. This is a computed
impact window, not the 600 s event: do not promote it as full-event validation.
The run took 22,620.93 s wall clock, which is why further long runs were not
used to close the remaining gates.

The ELVHYS registry join also exposed a fail-closed mismatch: CSV channels
include H2Bottom2/7/9 and H2Side7 although the README table labels those
positions temperature-only, while README concentration positions Bottom6 and
Bottom8 are absent from the CSV. Coordinate labels exist for the observed
channels, but sensor type identity is not reconciled. The audit records the
exact lists; no ELVHYS concentration score is allowed until the dataset owner
or an authoritative channel map resolves them.

## TI021 nested transient and interface-state correction (12 October)

The checkpointed 600 s outer R12 run has completed after resuming from its
hash-verified 370 s restart. It covered 88 × 128 × 8 cells at 5 m and 0.125 s,
with 4,800 conservative handoff queries. The conservation/numerical gates
passed (maximum continuity residual 1.24e-12 kg m-3 s-1, cumulative relative
species residual 9.47e-16, zero nonconverged correctors); the reported total
wall time, including the failed first attempt, was 2,430 s. This is only a
coarse R12 solver execution, not hydrogen validation.

The independent 145-sensor no-fit score failed quantitatively: normalized dose
bias -0.905, dose MAE 186.54 vol% s, 84 positive observations versus 53 model
detections at 0.1 vol%, and 30 missed arrivals among 82 observations with
arrival times. Nine nominal overlap receptors all lie exactly on the inner/outer
cut plane, so this is not a shared-volume convergence check. Three active
cut-plane sensors (S047-S049) also showed 32.9-35.2% peak and 20.1-25.1% dose
differences; their inner peaks occur about 11-12 s before the outer peaks.
Unconditional validation and overlap coupling remain withheld.

The audit isolated a thermodynamic-state error in the old handoff: the outer
inlet inferred primitive mass fraction and sensible enthalpy from *net* species
and energy fluxes, even though those fluxes include molecular diffusion and
conduction. Handoff schema v2 now transports the advective species and
sensible-enthalpy flux components separately while retaining the complete net
fluxes for conservation. The outer equation of state uses the advective
components; the species and energy updates still use the full net fluxes.
Version-1 HDF5 histories remain readable but are explicitly marked as using
the legacy state-from-net-flux fallback. A matched 60 s, 5 m TI021 prefix
the legacy state-from-net-flux fallback. A matched 60 s, 5 m checkpointed inner
and outer prefix was then run with the v2 handoff. Both runs passed their
numerical gates. At t=38 s near S047, the inner receptor was 0.595 vol%; the two
adjacent advective-state inlet cells were 0.567 and 0.601 vol%, averaging to
0.584 vol% at the cut-plane crosswind location (about 1.9% difference). Across
all nine cut-plane sensors, after the same 0.6 s sensor-response filter, the
maximum peak change was 4.56% and maximum dose change was 8.95%; all nine meet
the 10% coupling gate, although six sensors were zero in this 60 s prefix and
only three exercised the plume. This is direct evidence that the new
advective/net-flux split corrects the prior active-sensor boundary discrepancy
on the impact prefix. It does not yet re-score the 600 s full event, resolve the
145-sensor dose underprediction, or qualify H2. Relevant regression coverage:
19 focused handoff/outer tests and 63 resolved-flow/TI021/scoring tests passed;
the final legacy-reader compatibility test also passes in the focused suite.

## Updated v2 full-event no-fit result (13 October)

The v2 handoff was carried through the complete 600 s TI021 outer transient in
`outputs/thorney-ti021-outer-handoff-v2-600s-dx5m-dt0p125-2026-10-12-v1/`.
The 88 × 128 × 8 grid (5 m, 0.125 s) completed 4,800 steps in a solver-reported
2,752.06 s. Its numerical gates passed: maximum continuity residual
1.2373e-12 kg m-3 s-1, cumulative relative species-mass residual 9.85e-16,
zero nonconverged correctors, and maximum conservative-remap residual
4.31e-16. The final record is `accepted_completed_outer_transient`. These are
numerical/component gates, not proof of agreement with measurements.

The post-run, observation-blind 145-sensor no-fit score is preserved in
`outputs/thorney-ti021-nested-145-sensor-score-v2-600s-2026-10-13-v3/` and
remains withheld as unconditional quantitative validation. All 145 sensor
locations are covered (60 inner, 94 outer, 9 overlapping; inner preferred in
the nominal overlap), but a geometry check shows that all nine shared sensors
lie exactly on the y=100 m coupling plane. There are zero shared sensor
locations more than one 5 m cell inside both domains, so this sensor set does
not establish an interior overlap solution. The separate dose comparison also
fails: the maximum peak change is 4.56%, while maximum dose change is 12.96%
(S047/S048; S049 is 11.11%), above the preregistered 10% tolerance. Thus the
v2 advective/net-flux split fixes the large prefix peak discrepancy but does
not establish a grid-independent time-integrated interface signal. The score
now fails closed unless at least three distinct horizontal sensor locations
are within both domains beyond the handoff-plane buffer; a real spatial
overlap run is still needed.

### Short interior-overlap diagnostic (13 October)

The inner runner now accepts a declared domain extent and physical handoff
location, and the outer runner derives its origin and downwind domain length
from that recorded interface. A 60-s, observation-blind matched 5-m/0.125-s
diagnostic used a 300×200×40 m inner domain with handoff at site y=125 m. Both
solvers passed their conservation/corrector gates. Their 84 and 77 sensor
sets share 24 receptors at six distinct horizontal positions, all at least one
cell inside the overlap; no shared sensor lies on the coupling plane. The
maximum inner/outer change is 3.83% for 0.6-s peaks and 4.22% for integrated
dose over this 60-s window. Runtime was 104.8 s inner and 215.3 s outer.

This is a real interior coupling diagnostic, not the experimental 145-sensor
score: the chosen y=125 m boundary leaves eight sensors uncovered (S100–S103,
S124–S127), and 60 s is shorter than the full observed event. Therefore it
does not validate the concentration physics or promote any hazard map. The
reproducible JSON/plot are in
`outputs/thorney-ti021-interior-overlap-60s-dx5m-handoff125-2026-10-13-v1/audit/`;
the audit tool explicitly keeps full-event and promotion gates false.

Across the 84 sensors with positive measured dose, dose MAE is 185.53 vol% s,
normalized mean bias is -0.900, and only 3.57% fall within a factor of two.
Model detection at the 0.1 vol% criterion is 54/145 versus 84 positive
observations. For the 82 sensors with observed arrivals, 29 are missed and the
mean absolute arrival error among matched detections is 31.90 s; 24/68 observed
departures are missed, with 190.17 s mean absolute error among matches. The
maximum predicted 0.6 s sensor peak is 9.268 vol%, with 9 sensor peaks at or
above 4 vol%. These are model predictions only: the dose workbook does not
provide 145 independent measured peak concentrations, so the predicted LFL
sensor count is not an observed hazard-footprint validation. Large dose errors
are concentrated at the near-source sensors (e.g. S004: 2060.7 measured vs
69.48 predicted vol% s). These findings indicate severe dose underprediction
and short plume-duration mismatches; correlation alone (Pearson r=0.716,
Spearman rho=0.858) is not an adequate validation claim.

The resulting figures are `ti021-inner-time-maximum-concentration.png` and
`ti021-outer-time-maximum-concentration.png`. They are explicitly labelled as
R12 tracer fields, not hydrogen or flammability results. The inner view
contains the resolved obstacle; the outer view is the downstream
post-handoff domain and therefore intentionally has no obstacle overlay. The
plotter now ignores NaN solid cells when taking projection maxima, so the
4 vol% reference contour is not silently suppressed by obstacle cells.
Regression check: `tests/test_score_thorney_ti021_nested_145_sensors.py` (5
passed); the two full-event scoring and plotting commands are recorded by
their output artifacts above.

This full-event v2 score supersedes the earlier v1 numerical values above; the
earlier failure diagnosis and matched 60 s interface investigation remain
useful historical evidence. R12 is a tracer-gas component test. It does not
validate hydrogen thermodynamics, cryogenic buoyancy, independent obstructed
H2 experiments, or full-event finest-grid convergence.

### Failure localization: dose deficit is mainly persistent-exposure loss

The no-fit traces distinguish this from a simple non-arrival case. At S001,
the model predicts a 9.146 vol% filtered peak and 72.71 vol% s dose, while the
workbook reports 1832.0 vol% s and an arrival/departure window of 12--305 s.
The model's 0.1 vol% window is only 1.75--26.75 s. At S004 the corresponding
model values are 9.114 vol%, 69.48 vol% s and 1.25--32.125 s, versus 2060.7
vol% s and an observed 10--258 s window. Thus the computed short-lived peaks
do not sustain the measured long exposure; increasing peak amplitude alone
has not been established as the remedy. The dataset does not provide
independent measured peak time series for direct peak-to-peak validation.

The coarse full-event vertical grid is also visibly inadequate for these
sensors: its cell centres are z=2.5, 7.5, ... m, so every 0.4 m and 2.4 m
receptor pair at a common horizontal position uses the bottom-cell value. The
machine-readable score diagnostic reports 91/145 receptors at or below the
bottom-cell centre and all 42 of 42 such pairs as having bit-identical
predicted traces. In the 22 pairs with positive measured dose at both heights,
the measured-dose ratio has median
1.39 and maximum 5.65 (S001/S002 differ by a factor of 5.02: 1832.0 versus
364.936 vol% s). The 5 m-high report-priority fence is likewise represented
by only one vertical cell. This is a demonstrated resolution limitation, not
proof that interpolation is the sole source of the overall dose bias.

The selected frozen prior uses report meteorology (3.9 m/s and -6.1 degrees
relative to the array); the workbook alternative is 3.4 m/s and -15.8 degrees.
The report-priority fence is 5 m while the workbook gives 6 m. These
concentration-blind discrepancies are already registered as sensitivities,
not values to select by minimizing the sensor score. The initial condition
remains the report-scale instantaneous R12 cloud with quiescent cloud
momentum. A higher-resolution, observation-blind persistence/vertical-profile
test and the registered source/weather/geometry sensitivities are therefore
required before attributing the residual dose deficit to one physical closure.

## Promotion rule

Until the full event, hazard-volume convergence, receptor convergence, and
time-window-matched independent observations are available, the model may be
used for numerical development and conditional engineering screening only. It
must not be described as a fully validated field predictor.

### Follow-up checkpoint-prefix diagnostic (13 October)

The active 2.5 m TI021 inner-domain R12 run is a 60 s, 0.0625 s-step
observation-blind diagnostic using a zero-thickness obstacle surface. Its
20 s checkpoint is committed and its segment/restart hashes match the journal.
At that checkpoint the species relative residual is 4.15e-16, continuity
residual 9.54e-13 kg m-3 s-1, momentum residual 2.03e-9 kg m s-1, with zero
nonconverged corrector steps. This establishes numerical checkpoint integrity
only; the run is not complete and is not a hydrogen validation.

The same 2.5 m process has since advanced to a committed 50 s checkpoint
(current stdout has reached 54 s). The 50 s segment and restart hashes match
the journal; its species residual is 4.89e-16, continuity residual
2.18e-13 kg m-3 s-1, momentum residual 2.12e-9 kg m s-1, and corrector
nonconvergence count remains zero. Advective, species-diffusion, and thermal
diffusion numbers are 0.237, 0.272, and 0.258, respectively.

The existing
`tools/materialize_thorney_ti021_checkpoint_prefix.py` can extract a committed
prefix without waiting for the 60 s run to finish. The 5 m, zero-thickness,
20 s prefix and active 2.5 m, zero-thickness, 20 s prefix both pass their
numerical gates and use the same 1e-4 corrector tolerance. Their matched
physics configuration (apart from cell/time discretization) is confirmed by
the updated comparator. The earlier experimental v1 60 s reconstruction
failed its species-mass gate because segment-local outflow accumulators reset
at each boundary; it is retained only as a failed diagnostic artifact and
must not be used. The corrected v2 package passes its recomputed conservation
gates, but the established materializer is the reproducible workflow.

The checkpoint-prefix work also exposed a confound in the planned 5 m vs 2.5 m
comparison: the 5 m source uses `finite_thickness_voxel` and corrector
tolerance 1e-6, while the active 2.5 m diagnostic uses
`zero_thickness_surface` and tolerance 1e-4. These are physical/numerical
configuration changes, not grid/time-step changes. The comparison gate now
compares all configuration fields other than explicit grid/time discretization
keys, records mismatching fields, and labels the result as a physics-mismatch
diagnostic rather than grid convergence. A clean grid study still requires a
matched-physics 2.5 m run and a third 1.25 m representative window.

An existing 5 m zero-thickness 20 s prefix then enabled a valid matched-physics
two-level sensitivity using the active 2.5 m run's committed 20 s prefix.
There are 60 receptors, 12,800 common fully-fluid coarse cells, no measured
values accessed, and the comparator reports no non-grid/time configuration
mismatches. The spatial spacing halves from 5 to 2.5 m, while the timestep
changes from 0.02 to 0.0625 s; therefore this is a coupled space/time
sensitivity, not isolated spatial convergence. Final-field relative L1
difference is 44.8%; maximum-ever-field relative L1 difference is 29.9%. Mean
receptor peak and dose relative differences are 60.5% and 77.0%, respectively;
arrival differs by 0.86 s on average. The 4 vol% maximum-ever volume differs by
only 1.75% (14,875 vs 15,140.625 m3), showing that hazard volume can appear
stable while local concentration fields and sensor peaks are not. The record
is explicitly `diagnostic_two_grid_comparison_pending_1p25m_window`, with
`grid_convergence_promoted=false`; it is not a convergence or measurement
validation claim. Figure and machine-readable record are in
`outputs/thorney-ti021-inner-grid-diagnostic-20s-zero-surface-2026-10-13-v1/`.

### Completed 60 s, 2.5 m resolved R12 inner run (13 October)

The checkpointed run has now finalized at 60.0 s after 960 steps and 1,010.67
s wall time. It records 60 receptors × 961 time frames. All numerical gates
pass: species relative residual 2.76e-16, sensible-enthalpy residual 0,
continuity residual 6.77e-12 kg m-3 s-1, global mass residual
7.41e-11 kg s-1, momentum residual 2.45e-9 kg m s-1, and zero nonconverged
correctors. The largest advective/diffusion/thermal-diffusion numbers are
0.500/0.272/0.258. The final journal is `finalized`; the 6 segment hashes
and all 6 restart hashes were independently recomputed and matched. The solver
process exited successfully after writing its final record. Run record status is
`accepted_inner_transient_qualification`, which is a numerical qualification
for this 60 s R12 window only.

The time-maximum field reaches 32.046 vol% and the largest receptor peak is
18.114 vol%; model-only 4 vol% time-maximum volume is 15,140.625 m3. The
field image (plan/front/side projections; R12 and reference contour labelled)
is `outputs/thorney-ti021-inner-r12-fine-60s-dx2p5m-dt0p0625-zero-thickness-2026-10-13-v2/ti021-inner-60s-time-maximum-concentration-v2.png`.
The temperature remains exactly 291.15 K; density spans 1.203--1.235 kg/m3
at the final frame, so this remains an isothermal R12-tracer run, not a
cryogenic-H2 test. In the final lee-wake shell mean u is +2.085 m/s and only
2.10% of cells have u < 0. The weak reverse-flow fraction does not establish
a robust recirculation bubble.

A post-run no-fit score against TI021 is explicitly horizon-limited and
withheld: 60 sensors, 53 positive observed doses, and 43 model detections at
0.1 vol%. Predicted dose over only 0--60 s has MAE 269.63 vol% s, MAPE 0.8825
and normalized bias -0.900, but the workbook dose covers the longer published
event; those dose metrics are not a matched-window error. Nine of 52 available
arrivals are not detected by 60 s; the matched-arrival mean absolute error is
23.23 s. The six reported arc peak times lie at 86.5--133.5 s, outside this
run horizon, so peak-time and peak-amplitude comparisons are also provisional.
The score artifact is
`outputs/thorney-ti021-inner-r12-fine-60s-score-2026-10-13-v1/` and retains
the decision `withheld_as_horizon_limited`.

Runtime extrapolation to a 600 s equivalent at the measured cell-step rate is
2.80 h; this is an estimate, not a completed full-horizon run. No independent
obstructed-H2 data were accessed or validated by this exercise.

To isolate temporal discretization from the coupled grid/time comparison,
the same 2.5 m, zero-thickness, 1e-4-tolerance, quiescent-cloud physics was
checked over the common 0--20 s prefix at dt=0.02 and 0.0625 s. Workbook and
property-audit hashes are identical, both records pass their numerical gates,
and both runs declare no observation values used to configure them. On 102,400
fluid cells the time-maximum-field RMSE is 0.0293 vol% (maximum pointwise
difference 0.862 vol%); the 4 vol% volume differs by 140.625 m3 (0.94%).
Interpolated receptor histories have 0.0120 vol% MAE and 0.0540 vol% RMSE,
with a 1.227 vol% maximum pointwise difference; maximum receptor peaks are
17.643 versus 18.114 vol%. This suggests small aggregate 0--20 s temporal
sensitivity for this pair, with a notable localized peak discrepancy. Two
time steps do not establish temporal convergence or later-event accuracy.
The observation-blind comparison is recorded at
`outputs/thorney-ti021-temporal-2p5m-20s-dt0p02-vs0p0625-2026-10-13-v1/ti021-temporal-comparison.json`.

The predefined local sampling shells reinforce that caution: maximum-ever
4 vol% volume in the windward inner-face shell is 1,125 m3 (5 m) versus
2,156.25 m3 (2.5 m), and in the exterior lee-wake shell 250 m3 versus
812.5 m3. The roof shear band and side-tip shell remain below 4 vol% in both
runs, but their peaks also differ. These are fixed geometric samples of the
R12 tracer field; they do not by themselves prove stagnation, separation, or
recirculation topology. Near-equality of total hazard-threshold volume is an
aggregation effect that masks substantial local redistribution.

As a first flow-topology check, final 20 s checkpoint velocity fields were
sampled in the same fixed shells (positive x is the measured ambient-flow
direction). Ambient mean u is +4.068 m/s at 5 m and +4.059 m/s at 2.5 m. In
the exterior lee-wake shell the mean u is +1.827 and +2.250 m/s, and the
fraction of cells with u < 0 is 0% and 1.21%, respectively. The windward
inner-face shell has mean u +1.502/+1.640 m/s; the roof-shear band +2.975/
+3.240 m/s. This is consistent with deceleration and roof bypass, but does
not establish a robust recirculation bubble: the small, grid-dependent
reverse-flow fraction requires a resolved time-averaged streamline/vorticity
and wall-shear audit before assigning separation or wake-recirculation
physics. This velocity check uses model fields only, not observations.

The matched 20 s R12 runs also carry a nontrivial composition-density signal:
mixture density spans 1.203--1.299 kg/m3 on 5 m and 1.203--1.346 kg/m3 on
2.5 m, at a uniform 291.15 K. Thus the current surrogate case exercises
isothermal variable-density coupling, but it does not exercise cold-H2
thermal contraction or buoyancy reversal. Mean gauge pressure is positive in
the windward shell (+2.42/+2.35 Pa) and near zero in the roof/wake shells
(roof -0.07/-0.22 Pa; wake +0.42/+0.06 Pa). These are uncalibrated model-field
trends, not independent pressure validation.

### FFI source/sensor elapsed-origin audit (13 October)

The obstacle-present O1/O2 raw repeats were re-audited without reading H2
concentration values for synchronization. All four files have source-flow
and H2-sensor elapsed-time columns that each begin at exactly 0 s. Their
median recorded sample intervals are 0.01 s and 0.63--0.64 s, respectively.
That is consistent with the article's release-relative time convention, but
does not prove subsequent clock drift, acquisition latency, or calibration.

The detailed coverage audit found a more immediate blocker: the flow-meter
record ends 191.79--254.14 s before the corresponding sensor record ends, and
the last recorded flow is still non-zero in all four repeats (about 13--54% of
each file's own recorded maximum). Thus the source history does not establish
the full release shutoff and cannot drive or score the full sensor transient.
Positive-part trapezoidal integration of the recorded window gives 12.02 and
12.30 g for O1, and 56.68 and 56.44 g for O2; these are window integrals, not
event-total released masses. Common-clock and complete-event promotion remain
withheld. However, the fully covered 60--120 s source window is useful: mean
measured rates are 0.0883/0.0854 g/s for O1 and 0.4365/0.4381 g/s for O2,
within -2.9/-6.1% and -3.0/-2.7%, respectively, of the companion LES's fixed
0.091/0.45 g/s inputs. This supports no-fit transient-source replay through
120 s, but does not repair the later source-history gap. The detailed audit
and SHA-256 fingerprints are at
`outputs/ffi-shared-elapsed-origin-audit-2026-10-13-v4/`. Seven targeted tests
pass in the project `.venv`.

This result does not resolve the higher-priority blockers: the exact member
cross-sections and placement of the two-beam obstacle assembly, the closed-
rear/open-end facility boundary in a channel-scale run, the truncated source
history, or the missing sensor calibration certificate. No FFI model score is
produced or promoted.

### Matched zero-surface three-grid impact audit (13 October)

The already-completed observation-blind 5/2.5/1.25 m TI021 runs were
re-audited as a controlled grid family: all use the same 20 s window,
`dt=0.02 s`, zero-thickness ideal fence, property/workbook hashes, and
quiescent cloud. The new result verifies that the earlier mixed-dt/mixed-wall
audit was not the only available evidence. Common-fluid field L1/RMSE errors
contract, but the 4 vol% maximum-ever volume (14,875/15,000/15,431.6 m3) is
monotonic nonconvergent. Only 6/25 receptor peaks and 9/25 doses are
monotonic-convergent; side-tip bypass is absent at 5/2.5 m but the 1.25 m run
resolves a localized 5.35 vol% maximum and 19.53 m3 above 4 vol%. Thus field
average agreement masks unresolved threshold geometry and sensor metrics.
Status remains `impact_window_grid_sensitivity_not_converged_full_event_withheld`;
the fine-grid 20 s audit cannot substitute for 600 s convergence or H2
validation. Full numerical evidence and plot:
`outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v1/`.
Matched common-plane plan/streamwise/crosswind maximum-concentration contours
are in `outputs/thorney-ti021-three-grid-zero-surface-matched-dt-2026-10-13-v2/`;
these are per-cell time maxima over 0--20 s, not simultaneous cloud snapshots.
The audit and scorer regression tests pass (13 tests total); the plotting tool
now explicitly uses a noninteractive backend so the test is reproducible in a
headless Windows environment. The shared-plane contour sampler has two
additional unit tests.

### Conserved-mixture EOS inversion optimization (13 October)

The binary ideal-gas closure's inversion from species mass concentration and
sensible-enthalpy density previously used 64 vectorized bisection iterations.
It now solves the equivalent carrier-concentration quadratic using a
cancellation-resistant root and explicitly handles the zero-species and
pure-species limits. A deterministic set of 512 physical states across the
closure temperature/composition bounds, the limiting states, and the
3-D-corrector tests pass (16 targeted tests).

On the identical 2-step FFI control startup, wall time fell from 23.407 s to
11.408 s (2.05x); maximum H2 concentration changed by 4.44e-16 vol%-points
and the maximum mass residual remained 1.06e-22 kg. This is still a 0.01 s,
no-obstacle solver-path smoke and is not an experimental score. A separate
2 s, 5 m TI021 zero-surface transient completed all declared numerical gates
in 17.09 s wall time; compared with its pre-optimization 2 s artifact, common
field differences remain below 5e-11 in absolute units (the maximum is gauge
pressure), with species-mass residual 2.76e-16 and zero non-converged
corrector steps. The TI021 runtime did not materially improve, indicating
that this inversion was not its dominant cost. Artifacts are in
`outputs/ffi-control-transient-smoke-2026-10-13-v4/` and
`outputs/thorney-ti021-analytic-eos-smoke-dx5m-2s-2026-10-13-v1/`.
These execution checks do not upgrade obstacle, concentration, or full-event
validation status.

The subsequent FFI profile showed that Python per-cell sparse-matrix assembly
in the 3-D face-mass projection was the largest remaining startup cost. That
assembly is now vectorized over the six neighbor directions while retaining
the same solid-cell, embedded blocked-face, fixed-outlet and Jacobi-CG
operators. The 106-test obstacle/projection/EOS regression group passes; an
additional irregular-solid test closes x/y/z blocked faces and keeps the
projected divergence below `1e-8 s-1`. On the same FFI 0.01 s, 2-step smoke,
runtime decreased from 11.408 s (analytic-EOS-only v4) to 9.263 s (v5), with
maximum field difference `9.94e-14` and mass residual `1.06e-22 kg`. Relative
to the original 23.407 s v3 smoke, the combined optimizations are 2.53x
faster. The case remains no-obstacle, 0.01 s startup-only evidence; it does
not estimate runtime linearly for a 60--120 s facility validation.

### FFI obstacle-sensor uncertainty clarification (13 October)

The open-channel experiment paper reports XEN-5320 accuracy as ±1%, while the
manufacturer's preliminary XEN-5320 datasheet specifies 1% full-scale
inaccuracy. If that full scale is 0--100 vol%, its conditional interpretation
is ±1 vol%-point, or 50% of a 2 vol% reading and 25% of a 4 vol% reading. The
article's shorthand does not state the reference basis, and no event-specific
calibration certificate is included with the local raw package; therefore the
registry preserves both attributed claims and does not assign event-level
accuracy. This is distinct from the ±3% full-scale value cited for the
separate FFI Test 4 dataset. The FFI obstacle score remains withheld. Sources:
[Henriksen et al. (2025)](https://doi.org/10.1016/j.jlp.2025.105669),
[Xensor XEN-5320 datasheet](https://xensor.nl/images/pdf-files/XEN-5320/Thermal_Conductivity_Gas_Sensor_XEN-5320-HP-USB-A.pdf).

### TI021 nested-domain handoff correction (13 October)

The earlier nested 600-s score's nine shared sensors were all exactly on the
handoff plane, so their agreement did not verify shared-fluid coupling. The
corrected handoff at site `y=100 m` plus a 300-m inner crosswind domain yields
145/145 sensor coverage and 24 true interior overlap sensors at six locations.
The 60-s diagnostic has peak/dose changes of 4.86%/5.61%; the corrected 600-s
run retains the peak agreement (4.65%) but integrated-dose change reaches
19.62%, failing the 10% coupling limit. Inner and outer runs passed all their
own numerical gates, with runtimes 414.4 s and 1573.1 s. The full TI021 score
still materially underpredicts exposure: dose MAPE 85.0%, normalized bias
-89.85%, and arrival MAE 31.76 s (29 of 82 available arrivals missed). A
piecewise trace audit at S058 finds inner/outer doses 7.402/7.725 vol%·s
(0--60 s), 19.057/15.567 (60--120 s), and 2.914/0.343 (120--200 s); most of
the mismatch accumulates after 60 s. This points to a long-horizon solution
or coupling inconsistency to investigate, but is not yet a causal diagnosis.
Full-event evidence is now in
`outputs/thorney-ti021-nested-145-sensor-score-v3-600s-2026-10-13-v1/`; no
experimental validation or model promotion is allowed. The earlier full-run
configuration remains preserved as a baseline.

### Inner crosswind-span diagnostic (13 October)

The first 200 s, 640 m-wide pilot is invalidated: receptor solver-y remained
referenced to the old 200 m center and handoff `q` edges remained at +/-100 m.
The scorer's coordinate hardening caught 57 bad inner samplers and the origin
error; do not cite that pilot as model evidence. Its corrected 120 s successor
uses a 320 m crosswind origin and +/-320 m handoff edges and has zero sampler
mismatches. Against the valid 200 m-wide, 450 m-long baseline over the same
0--120 s prefix and 36 common interior sensors, the 12 active sensors show
maximum peak/dose coupling changes of 2.30%/7.20%, versus 4.87%/12.97% at
200 m. The 640 m run adds valid interior coverage (54 sensors, 15 locations),
but its all-sensor maxima remain 17.40% peak and 116.33% dose due to very weak
signals; the existing acceptance gate therefore still fails. This supports a
crosswind-domain-sensitivity diagnosis for active-plume coupling only. It does
not supersede the 600 s result (19.62% dose coupling change), establish full
event convergence, or validate against hydrogen measurements. See
`outputs/thorney-ti021-wide-domain-overlap-score-120s-centered-2026-10-13-v1/`.
