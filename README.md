# DEGALI — Dense Gas Dispersion for Liquid Hydrogen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22646258.svg)](https://doi.org/10.5281/zenodo.22646258)
[![PyPI](https://img.shields.io/pypi/v/degali.svg)](https://pypi.org/project/degali/)

DEGALI stands for **Dense Gas Dispersion for Liquid Hydrogen**. It is a modern
Python model built from a verified reimplementation of the US EPA DEGADIS 2.1
dense-gas dispersion model and extended for cryogenic hydrogen releases.

> **Alpha research software.** The DEGADIS 2.1 compatibility path is strongly
> regression-tested. The liquid-hydrogen extensions are suitable for research,
> scenario comparison and sensitivity studies, but are not independently
> certified for regulatory separation distances or safety-critical design.

한국어 요약: 현재 공개본은 원본 DEGADIS 재현 경로와 액화수소 연구 확장을
함께 제공합니다. PRESLHY 수평 야외제트 범위에서 정량 비교를 마쳤지만,
미완성 TKE/압력/열폭 폐쇄와 제한된 독립 검증 때문에 설비 인허가 판단의
단독 근거로 사용하면 안 됩니다. 자세한 판정은
[1차 결과](docs/stage1-results-2026-09-06.md)를 참고하십시오.

처음 설치하거나 실제 계산을 시작하려면 상세 사용 가이드를 참고하십시오:
**[한국어](USER_GUIDE_KO.md) · [English](USER_GUIDE.md)**. LH₂ 간편
평가, 결과 해석, CSV 저장, 민감도 분석, 기존 DEGADIS 입력 덱과 CLI 사용법을
단계별로 설명합니다.

## What is included

- Python ports of all six DEGADIS 2.1 programs. The independently obtained
  Fortran oracle is retained locally and is not redistributed.
- Steady, transient, receptor-dose and jet-to-ground workflows.
- Legacy and modern thermodynamic backends.
- Opt-in LH2 source flashing, buoyant trajectory, air condensation/freezing,
  component enthalpy, ground interaction and thermal-profile research paths.
- Experimental finite-TKE, independent thermal-width and Reynolds-stress
  operators. These are deliberately marked as research-only until their
  physical closure inputs and downstream field performance are validated.

## Current evidence

The port reproduces the five EPA reference cases and was checked locally
against a source-built Fortran implementation. The first frozen LH2 assessment uses
seven horizontal outdoor PRESLHY trials, 38 concentration sections, 17
vertical profiles and 41 temperature sensors.

| Quantity | Current provided path |
|---|---:|
| Concentration MG (observed/predicted; ideal 1) | 1.102 |
| Concentration VG (ideal 1) | 1.233 |
| Concentration FAC2 | 36/38 |
| Vertical width, predicted/measured | 1.091 |
| Plume-centre MAE | 0.054 m |
| Minimum-temperature MAE | 27.25 K |

An opt-in thermal-profile candidate lowers minimum-temperature MAE to 16.47 K
but worsens other acceptance metrics, so it is not the default. These samples
are correlated observations from one campaign, not 38 or 41 independent
validation experiments.

### Current IJHE manuscript baseline

The manuscript freeze used for the IJHE submission track is reported separately
from the historical first-pass table above. On nine momentum-dominated
PRESLHY E3.5 releases (62 arc maxima), the pre-specified peak operator gives
MG = 1.047, VG = 1.425 and FAC2 = 0.839. A synchronised 20 s mean operator
on seven trials gives MG = 0.544, VG = 4.601 and FAC2 = 0.761, while the
independent six-arc FFI/DNV screen gives MG = 1.245, VG = 1.373 and FAC2 =
0.833. These are separate evidence lanes, not one pooled accuracy score; the
low-wind FFI Test 6 residual remains an explicit application limit.

The accompanying 3-D transient solver is a numerical research extension. Its
grid-convergence and mass-conservation results are documented in the
[IJHE submission-readiness record](docs/ijhe-submission-readiness.md) and do
not constitute measured field validation.

A separate finite-block Tier-3 prototype now couples D3Q19 momentum, bounded
MUSCL hydrogen/enthalpy transport, ideal-gas mixture density and a
variable-density mass-flux pressure correction. Its manufactured conservation
checks pass, but the active three-level time and grid studies remain withheld:
roof/side low-concentration metrics are not grid independent. It is therefore
not enabled as a quantitative facility-design claim. See the
[obstacle-physics rebuild record](docs/obstacle-physics-rebuild-plan-2026-10-10.md).
The product target is a dynamic local-facility consequence model that
complements, rather than duplicates, a free-field slab model. An opt-in WALE
closure is implemented, but a fixed-grid AIJ comparison found no material
accuracy gain over either Smagorinsky or no SGS and nearly doubled runtime;
the default has therefore not been changed.
The carrier-flow observation operator can now report a post-spin-up time mean.
On AIJ Case H its 1000-step window passes a predeclared 1% first-half/second-
half stationarity test (`0.8015%`), while grid convergence and hydrogen
concentration validation remain explicitly withheld. An opt-in non-equilibrium
extrapolation outlet was also tested and rejected because it increased AIJ
errors and the instantaneous residual.
The subsequent AIJ grid audit identified and removed a geometry-comparison
error: the old `0.020 m` grid voxelized the nominal `0.10 m` block as
`0.12 m`. With all block faces aligned on `0.025/0.01667/0.0125 m` grids,
model-field changes decrease monotonically and the carrier-flow component grid
gate passes. Obstacle voxelization bounds/errors are now explicit diagnostics;
this does not promote the still-unconverged hydrogen roof/side concentrations.
The Gaussian release is now integrated over each control volume, so nested
fine-cell source weights exactly recover the corresponding coarse-cell weight.

An opt-in unified finite-volume active solver now removes the incompatible
stacking of thermodynamic pressure correction onto isothermal BGK pressure
dynamics. It advances momentum, hydrogen and sensible enthalpy with one shared
face mass flux and one variable-density projection. Four-level time refinement
is first order and the aligned three-grid numerical component gate passes;
global momentum residuals are below `9e-15 kg m/s`. Release momentum velocity
is an explicit input. This candidate is still research-only: independent
obstacle hydrogen experiments have not yet qualified it, so the legacy default
and the prohibition on facility-design claims remain unchanged.
The earlier 26--36% roof/side refinement failure belongs to the superseded
isothermal-BGK/thermal-pressure feedback route; the unified finite-volume path
passes the current numerical component gates. That is not an experimental
accuracy claim.

Thorney Island Trial 021 is now represented by a native semicircular-wall
solid rather than a rectangular surrogate. A three-grid geometry-only audit
reduces arc-length error from 5.68% to 1.22% to 0.27% while preserving wall
connectivity, source clearance and the report-priority 5 m height. The local
evidence audit retains all 145 dose sensors (84 positive) and six arc
summaries. A concentration replay is intentionally withheld because Trial 021
is an initially finite, relative-density-2.02 R12/air cloud. Schema v20 now
adds a separate conservative full binary-mixture initial-cloud path: the
twelve-sided holder is integrated exactly into finite volumes, species mass is
transported by total projected mass flux with mass-fraction diffusion, and EOS
density/buoyancy feed the momentum solve. Its 32 mol% initial state reconstructs
58.389 g/mol and relative density 2.0158, while the 14 m across/13 m fill gives
2048.204 m3 versus the reported 2050 m3. The molecular property gate is
frozen separately: the measured CFC12-air correlation gives
`9.78865e-6 m2/s` at 291.15 K with reported +/-2% uncertainty, and the exact
CoolProp 8.0.0 revision supplies an ambient R12 cp of 606.0749 J/kg/K. This
does not qualify or tune turbulent scalar mixing. A separate
`projected_inlet` initialization now lets the finite-volume momentum path keep
that physical molecular viscosity without forcing a metre-scale atmospheric
flow through an invalid near-0.5 BGK relaxation time. The 1 m, 291,200-cell
full-site single-step preflight represents the 5 m wall exactly and passes
species-mass, continuity, boundedness, Mach and nonlinear-corrector gates in
32.4 s. Its cached symmetric smoothed-aggregation AMG pressure preconditioner
cuts the matched 2 m projection iteration count 40.31-fold and total runtime
2.87-fold relative to Jacobi-CG without changing the declared mixture-state
or continuity gates. The neutral surface-layer closure is now frozen without
reading TI021 concentration, dose or arc observations. It applies
`nu_t=kappa*u_star*z` using the reported `u_star=0.205 m/s` and
`z0=0.005 m`; all six combinations of `Sc_t={0.2,0.7,1.3}` and
`Pr_t={0.74,1.0}` are mandatory sensitivities, so the primary 0.7/0.74 member
cannot be selected after scoring. The 1 m, 291,200-cell primary-closure
preflight passes all numerical and geometry gates in 34.0 s. This is still
only an execution/numerical gate. A subsequent sensor-coverage audit found
that the preflight box contains only 33 of 145 sensors; the network reaches
500 m downwind and +/-300 m crosswind. Direct 1.25 m simulation of the
monolithic 600 x 640 x 40 m sensor-complete domain is estimated from measured
DEGALI cell-step rates at no less than 1,730 h per prior member. The coordinate-
consistent nested plan instead uses a 200 x 200 x 40 m resolved inner domain
containing the complete fence and all six 75/100 m arc centreline observations.
Its physical site-local y=100 m plane is inner solver x=160 m and maps exactly
to x-min of a 440 x 640 x 40 m outer domain spanning y=100--540 m. The inner
and outer cover 60 and 94 sensors with a 9-sensor overlap; their union covers
all 145. The larger 600 m domain is retained only as a monolithic feasibility
reference. The first handoff operator now
conservatively remaps signed total-mass anomaly, species
mass, sensible enthalpy and lateral/vertical momentum fluxes from a 1.25 m
inner plane to a 5 m outer plane with relative integral residuals below
`1e-12`, without clipping or renormalisation. The time-interpolated outer binary boundary is now
implemented and independently passes species (`1e-10 kg`), enthalpy (`1e-6 J`)
and continuity (`1e-7 kg/m3/s`) gates on an 8,192-cell coupling audit. It
reconstructs boundary composition, temperature, density and velocity from
mass/species/enthalpy flux ratios and holds the prescribed x-min total mass
flux fixed through pressure projection. The coupled outer execution step now
also applies reconstructed transverse/vertical velocity to the finite-volume
momentum predictor and feeds the projected face state into the next step. A
48-step, observation-blind same-grid audit reconstructs both tangential
momentum fluxes to round-off and gives maximum common-receptor peak/dose
changes of `7.35e-9`/`2.15e-8`, below the 10% handoff gate. This deliberately
does not promote the 5 m mesh: a separate 2.5-to-5 m diagnostic changes peak
by 18.81% and dose by 45.12%. The actual 200 x 200 x 40 m inner solver now also
streams the physical x=160 m plane to append-only HDF5 in the shared
`q=-site_x` coordinate and conserves all five remapped flux integrals to
`1.94e-16`. The HSE report's stationary-cloud description (p. 150) is now an
explicit `quiescent` initial-momentum condition rather than silently assigning
the ambient wind to the dense cloud. With that report-based condition, the
primary 5 m inner run completes 600 s, all 60 inner sensors and 4,801 handoff
frames in 501.1 s. Species mass closes to `3.59e-15` relative, continuity stays
below `3.0e-12 kg/m3/s`, and no thermodynamic corrector step fails. Its first
observation-blind score is intentionally not promoted: the six 75/100 m arc
peaks are underpredicted by 47.8% MAPE, dose by 81.5% MAPE, and detected
arrivals are 25.2 s early. At 5 m, the fence has one cell and 0.4/2.4 m sensors
share the first vertical cell. The evidence status is therefore
`conditional_pending_outer_transient_and_grid_convergence`; 2.5/1.25 m checks
and the outer 145-sensor transient remain.

The checkpointed 2.5 m inner member is now complete through 600 s in 6,085.8 s
wall time, with all numerical gates passing: maximum relative species-mass
residual is `6.91e-16`, continuity residual is `9.83e-12 kg/m3/s`, and no
nonlinear-corrector step fails. Resolution does not close the validation gap.
Across the six published arc summaries it underpredicts peak concentration by
51.9% MAPE and dose by 79.5% MAPE, while peak time is 60.7 s early in absolute
mean. The 60-sensor dose rank structure remains strong (`r=0.847`,
Spearman `rho=0.914`) but dose MAPE is 86.5%. The full-event 5-to-2.5 m
comparison changes maximum-field relative L1 by 28.9%, receptor peak and dose
by 36.5%/28.1% MAPD, despite only a 1.67% change in maximum-ever 4 vol% volume.
It is therefore a non-converged diagnostic, not quantitative validation. A
model-only 0--10 s peak-window replay also passes all predeclared `1e-4`
versus `1e-6` corrector gates (field relative L1 `2.91e-4`, receptor Linf
`0.00499 vol%`), excluding production tolerance as the source of the large
physical bias. The observation-blind 1.25 m impact-window member is the next
required grid-falsification result; full-event finest-grid convergence and the
outer 145-sensor transient remain withheld.
The consolidated schema-v12 evidence record hash-pins both 600 s members,
both no-fit scores, the two-grid comparison and the strict replay in
`outputs/thorney-ti021-evidence-audit-2026-10-11-v16/`.
The first 1.25 m impact-window attempt is retained as a rejected numerical
member: `dt=0.03125 s` reached species/thermal diffusion numbers
`0.553501`/`0.523581`, above the explicit 0.5 limit. The revised observation-
blind plan uses `dt=0.0625/0.03125/0.015625 s` on the 5/2.5/1.25 m grids,
preserving an exact 2:1 accepted time-refinement path while keeping the failed
attempt visible. The stable-`dt` v2 execution was deliberately stopped at
`t=0.328125 s` (22 frames) before its first 10 s atomic checkpoint after
measured throughput showed an approximately seven-hour recovery exposure.
It remains unscored, its closed HDF5/progress hashes are preserved, and no
physics was changed. The replacement v3 execution uses the identical grid,
time step, tolerances and boundaries with a 1 s atomic checkpoint. The
reproducible revised plan is in
`outputs/thorney-ti021-three-grid-impact-window-plan-2026-10-11-v5/`; schema
v3 hash-binds both the rejected stability member and the unscored operational
interruption. It
projects the fine-grid advective/species/thermal diffusion numbers as
`0.2497/0.2724/0.2576`, passing both explicit stability gates. No
result from either superseded member is scored.
The replacement's live progress samples reached model time 0.25, 0.50 and
0.75 s in 581.47, 1126.54 and 1694.89 wall-clock seconds. This confirms active
advancement. The latest positive interval projects the first 1 s checkpoint at
`2263.24 s` and the full 20 s window at about `12.63 h`; the machine-readable
projection is in `outputs/thorney-ti021-fine-runtime-projection-2026-10-11-v2/`.
This is runtime evidence only, not a completed physics result.
The next physical-bias investigation is preregistered in
[the stratification audit plan](docs/thorney-ti021-stratification-bias-hypothesis.md):
the current neutral surface-layer eddy diffusivity lacks local dense-layer
stratification damping. It will first be diagnosed without altering fluxes;
no Richardson damping coefficient may be selected from TI021 observations.
The diagnostic-only finite-volume `Ri_g` operator and fixed-region summary are
now implemented and tested, including explicit shear regularisation and a
separate undamped neutral-diffusivity exposure. They are not coupled into the
running baseline and do not alter the production closure.
The required regularisation sensitivity was frozen before diagnostic output at
`0.0205/0.041/0.082 s^-1 = 0.5/1/2 * u*/H_wall`; all three members are mandatory
and non-promoting. The plan is in
`outputs/thorney-ti021-stratification-diagnostic-plan-2026-10-11-v1/`.
The envelope driver executes and retains all three plan members; it reports
regional and receptor-cell ranges without a best-member selector.

An energy-conserving high-pressure-source adapter now prevents choked or
supersonic notional-nozzle states from being injected directly into the
low-Mach tier. It entrains declared ambient co-flow until both the physical
Mach and active-scalar density-departure gates pass, then maps the conserved
H2, momentum and sensible-enthalpy increments to a cell-integrated anisotropic
source. An opt-in tunnel mode applies no-slip sidewalls and roof for genuinely
confined facilities. It is not the HyTunnel boundary: the official HyTunnel
CFD prescribed the inlet at `x=0` and used open non-reflecting conditions on
the other five faces. DEGALI now provides an opt-in five-face pressure-open
candidate: adjacent-plane velocity/non-equilibrium stress, ambient pressure,
open mass-flux projection and scalar outflow are coupled consistently. It has
passed the declared one-dimensional linear acoustic tests: reflected amplitude
falls from 91.72% for fixed pressure to 1.64% for the characteristic form at
both the streamwise and lateral-normal faces. A 30° single-mode test passes at
4.188%, but that one-dimensional-impedance form fails the broader sweep at
45°/60° long waves (8.268%/12.774%). A separate angle-aware characteristic
candidate uses the local normal/tangential acoustic velocity ratio rather than
a fitted incidence angle; all seven resolved angle/wavelength cases pass, with
a worst reflection of 3.907% at 60°. A localised 45° packet aimed at the
x-max/y-max intersection also passes a conservative post-direct-pulse gate at
2.594%. The first angle-aware candidate was wavelength-sensitive: its scaled
8-cell x-y/x-z/y-z/xyz matrix passed at 4.142%, 4.154%, 4.142% and 3.644%,
but the 10-cell matrix failed at 6.076%, 6.073%, 6.075% and 5.999%. A separate
second-order-extrapolated angle-aware candidate removes that failure without a
fitted coefficient. Its corresponding 8-cell ratios are 2.418%, 2.439%,
2.420% and 1.846%, and its 10-cell ratios are 3.192%, 3.179%, 3.191% and
3.152%; both declared 3-D matrices pass the 5% gate for that raw-extrapolation
candidate. It is not a unified boundary solution: the same candidate becomes
unstable in a finite-vortex audit. A separate TVD-limited, local
dilatation/curl-separated candidate remains stable and keeps vorticity error
to 0.943%, but its velocity/density contamination is 19.825%/243.230% and its
10-cell acoustic matrix is 5.502--5.564%, so it also fails. Vortical outflow
remains unqualified. A
same-grid 0.120343 s plume audit comparing the published-size box with a box
expanded downstream and laterally passes: the common-field normalised L1 is
1.092%, the worst primary-output change is 3.865%, and the 4 vol% volume changes
0.235%. An explicit cached-direct pressure projection reduces the expanded
63,000-cell run from 477.6 s to 220.4 s and now passes the separate 5 min
development target without changing those concentration metrics. The official
0.707246 s/3,056-step audit is now complete: worst primary-output change is
0.9901%, common maximum-field L1 is 1.2447%, and the enlarged domain clears
the lateral boundary gate. Its 1,402.2 s runtime fails the 5 min development
target. The solver remains opt-in because the large sparse factorization is
memory intensive. The official CFD also uses a different notional-nozzle
definition, so unified vortical boundary qualification and source-model
sensitivity remain required.

A deliberately simpler no-fit check has already been completed: the same
conservative source/handoff was connected to the existing axisymmetric
free-jet model and compared with 108 published HyTunnel average-concentration
points downstream of the handoff. Although Pearson `r = 0.937`, it
overpredicts by `+2.688 vol%p` on average, has `FAC2 = 0.519`, and predicts a
3.413 m 4 vol% centreline length instead of the measured 1.65 m. This rejects
the free-jet shortcut for this forced-flow case; it is not a score for the
resolved 3-D open-box model. A separate supported-overlap score now exists:
the 6% handoff and 2.52 m common box support 54 of 153 average points without
near-field or far-field extrapolation. A deterministic replay reproduces the
frozen maximum slice exactly and its last-20% receptor peak-to-peak drift is
only 0.00687%, qualifying the final 0.707 s slice as a steady-endpoint proxy.
Against the averages it gives MAE/RMSE 2.179/2.491 vol%p, `r = 0.733`, FAC2
0.556, 4 vol% sensitivity 1.0 and specificity 0.522. The accuracy gates still
fail, so this score is withheld.
Reproduce the free-jet route with
`python tools/validate_hytunnel_free_jet_handoff.py`. Reproduce the resolved
overlap score with `python tools/score_hytunnel_resolved_transient.py`.

A frozen, no-selection sensitivity audit separates the source handoff from the
effective scalar-mixing closure on 45 common points. At the 6% handoff,
`Sc_eff = 0.5/0.7/1.0` gives MAE 1.670/2.113/2.529 vol%p; tightening the
handoff to 5% at `Sc_eff = 0.7` gives 2.185 vol%p. Thus the current error is
more sensitive to the unvalidated mixing closure than to this admissible
handoff change. Even the lowest-error perturbation has FAC2 0.689 and 4 vol%
specificity 0.675, so observations select no coefficient or default. Reproduce
with `python tools/audit_hytunnel_source_mixing_sensitivity.py`.

That four-case result is now labelled a **legacy uniform-diffusivity
sensitivity**, not a qualified turbulent-mixing study. In schema v16 the
resolved solver can instead transport species with the local field
`D_molecular + nu_sgs / Sc_t` and enthalpy with
`alpha_molecular + nu_sgs / Pr_t`. Cellwise coefficients use conservative
face-averaged diffusive fluxes, and the finite-volume active-momentum route
updates its SGS viscosity from the evolving velocity field. The previous
uniform `nu_base/Sc` route remains selectable only to reproduce frozen v14
artifacts. In a 100-step HyTunnel preflight, the local-SGS closure forced the
source admissibility decision from the former 6% handoff (10.46% density
departure, failed) to the predeclared 5% candidate (9.25%, passed); neither
preflight is an external accuracy claim.

The full 3,056-step 5%-handoff local-SGS run is also complete. Its endpoint is
stationary and its conservation/applicability gates pass, but the common
45-point score is worse: MAE/RMSE `3.425/3.976 vol%p`, bias `+3.380 vol%p`,
`r = 0.396`, FAC2 `0.378`, 4 vol% sensitivity `1.0` and specificity `0.25`.
The matched legacy 5%-handoff surrogate had MAE `2.185 vol%p` and specificity
`0.45`; that is not a reason to restore it as physics. It shows that uniform
numerical diffusion was compensating for unresolved mixing. Reproduce the
fixed comparison with `python tools/compare_hytunnel_scalar_closures.py`.

The active finite-volume route can now impose a complete nonuniform `z-y`
mean-velocity inlet plane. For HyTunnel, the public 3.5 m/s co-flow workbook's
first average plane (`x = 0.25 m`) is shifted to the computational inlet with
its SHA-256, sheet, statistic and coordinates recorded. Its five measured
velocities span 2.902--4.241 m/s; interpolation on the published-size grid
produces a 2.907--4.231 m/s inlet with mean 3.525 m/s. The workbook provides
no velocity time series, RMS turbulence, integral scale, spectrum or
dissipation, so none is invented. The full 3,056-step measured-mean-fan
candidate is stationary and conservative, but it does not improve accuracy:
MAE/RMSE are `3.432/3.986 vol%p`, bias is `+3.383 vol%p`, `r=0.392`, FAC2 is
`0.356`, and 4 vol% sensitivity/specificity remain `1.0/0.25`. Relative to the
uniform-in-y local-SGS candidate, MAE worsens by only `0.0073 vol%p`; the mean
fan profile is therefore not the missing physics. Its runtime is 385.7 s and
misses the 5 min development target. Reproduce with
`--fan-inlet-profile measured_first_station`; the three-way comparison is
`outputs/hytunnel-scalar-closure-comparison-2026-10-11-v2.json/.png`. The
standalone evidence gate is reproduced by
`python tools/audit_hytunnel_fan_flow_evidence.py`.

The public PRESLHY E3.4 pool-spill archive now supplies a separate dynamic
component path. Four wind workbooks contain corrected pool mass plus
synchronized 10 Hz three-component sonic wind and local H2 at 0.35/0.45/0.55 m.
The workbook flags CH35 defective in all four trials and flags the Sand02 wind
channels defective, leaving three trials for conditional dynamic replay. A
mass-only preregistration step freezes every post-fill decay before any H2
value or model result is read: 13 windows are retained across the four files.
Run `python tools/freeze_preslhy_e34_source_windows.py ...`; the reviewed record
is `outputs/preslhy-e34-source-window-freeze-2026-10-11-v2/`. These windows fix
event timing but do not turn liquid-mass loss into a measured atmospheric
vapour flux. Report-figure-only obstacles, nominal rather than event sensor
calibration, and the unresolved sonic-to-site coordinate transform keep model
qualification withheld.

The first concentration-sealed Sand03 replay is now complete. Window 1
(`437.62--706.12 s`) uses the frozen mass-loss proxy, a 10 s block-mean Young
U/V/W inlet in its instrument-mean-wind frame, and CH45/CH55 as vertical
pool-centre receptors. The complete `268.5 s`, `26,850`-step run takes `241.4 s`
on the current machine and passes its numerical preflight and full-run gates:
maximum lattice Mach `0.1424`, maximum mixture-density departure `1.184%`, and
relative H2-mass/energy residuals `1.50e-13`/`1.93e-13`. These numerical results
do not imply accuracy. After both model and scoring policy were hash-frozen,
the sealed observations were opened and the primary accuracy gate failed:
median peak error `21.53 vol%p`, median relative dose error `124.4%`, and pooled
4 vol% sensitivity/specificity `0.594/0.787`. CH45/CH55 pointwise MAE are
`4.280/2.305 vol%p`.

The failed primary is retained rather than tuned away. The observation already
contains `15.86/3.52 vol%` at event time zero while the numerical field starts
from zero, because the mass-only scoring window begins after plume generation
during filling. Both model peaks also occur at `109.96 s`, in the lowest inlet-U
quartile: model concentration correlates with inverse applied inlet U at
`0.889/0.858`, versus only `0.085/0.057` in the observations. This diagnoses a
missing source prehistory and an over-strong whole-inlet response to a local
sonic lull; it does not authorize a time shift, source rescaling or closure
fit. Reproduction and immutable records are under
`outputs/preslhy-e34-sand03-window01-preflight-2026-10-11-v1/` via
`tools/preflight_preslhy_e34_dynamic_component.py`,
`tools/score_preslhy_e34_dynamic_component.py`, and
`tools/diagnose_preslhy_e34_dynamic_component_score.py`.

A post-primary fill-prehistory diagnostic then uses corrected mass and wind
only. A fixed 10% baseline-to-peak mass crossing starts the model at
`308.12 s`, `129.5 s` before the retained score window; the frozen post-fill
rate is extended backward without H2 fitting. The 398 s/39,800-step diagnostic
runs in `365.5 s` and passes its numerical gates. It improves median peak error
from `21.53` to `14.63 vol%p` and normalized arrival error from `150%` to `55%`,
but median dose error remains `116.3%` and pooled 4 vol% sensitivity/specificity
are only `0.590/0.792`. Prehistory is therefore a real but secondary error
source. The diagnostic cannot replace the failed primary or serve as
independent validation. Reproduction:
`tools/build_preslhy_e34_fill_prehistory_diagnostic.py`,
`tools/run_preslhy_e34_fill_prehistory_diagnostic.py`, and
`outputs/preslhy-e34-sand03-window01-fill-prehistory-2026-10-11-v1/`.

Schema v17 also accepts a synchronized three-component inlet history on the
finite-volume active momentum path. It is linearly interpolated on the solver
clock and passed through momentum, pressure correction and conservative scalar
transport. The measured vector is applied over the declared spatial inlet
shape; this is a temporal replay capability, not a synthetic spatial-turbulence
inlet. E3.4 is not sent through that boundary until its coordinate transform
and any x-min flow reversals have a declared treatment.

The first resolved-source preflight also found a separate numerical gate:
passing the 10% density-departure limit on the uniform handoff plane does not
guarantee that its moment-matched Gaussian source profile passes locally. On a
declared 45 mm coarse grid, the nominal mapping reaches 13.342% and fails; an
8% mean-plane candidate reaches 11.562% and also fails, while a 6% candidate
reaches 9.664%. This uses the characteristic pressure-open candidate. No
concentration data enter this check. The 6% case is only a
candidate for the subsequent grid/time/domain study, not a fitted or validated
replacement. Reproduce it with
`python tools/audit_hytunnel_resolved_source_preflight.py`.
Reproduce the boundary reflection audit with
`python tools/audit_pressure_open_boundary_reflection.py`.
Reproduce the 30° audit with
`python tools/audit_pressure_open_boundary_oblique_reflection.py`.
Reproduce the angle/wavelength audit with
`python tools/audit_pressure_open_boundary_sweep.py`.
Use `--boundary-scheme characteristic_multidimensional_extrapolated` to
reproduce the linear-acoustic multi-wavelength candidate; it is not qualified
for vortical outflow. The default preserves the failed 1-D-impedance baseline.
Reproduce the two-face corner audit with
`python tools/audit_pressure_open_boundary_corner_reflection.py`.
Reproduce the finite-vortex audit with
`python tools/audit_pressure_open_vortical_outflow.py`.
Reproduce the short-horizon plume-domain audit with
`python tools/audit_hytunnel_plume_domain_sensitivity.py`.

A separate bounded enclosure-average dynamic model is evaluated on the public
USN ventilated-container GH2 record. Exp19 identifies three closure
coefficients; frozen Exp20/21/23 parameter holdouts give pooled MAE 0.282
vol%p, RMSE 0.332 vol%p and Pearson r = 0.980 over 900 s spatial-mean
histories. Their spatial-mean 4 vol% first-arrival MAE is 3.70 s and duration
MAE is 27.1 s; these are not any-sensor alarm scores. This is retrospective
GH2 enclosure evidence, not prospective
blind, individual-sensor, outdoor LH2 or CFD validation. Reproduce it with
`python tools/validate_usn_enclosure_dynamic.py`.

The public H2SAFE record (DOI `10.7799/17118570`) is also available as an
independent full-scale indoor geometry and forced-ventilation evidence lane.
DEGALI hash-pins the official archive/document, verifies all extracted members,
joins five one-second helium traces to their 24/37-sensor registries, checks
source/sensor coordinates against the STEP/STL frame, and inventories named
HVAC assemblies. It explicitly maps the published vertical `Y` axis to
DEGALI's solver `z` axis. The typed HVAC audit preserves ten Lab-1 boundaries:
the reported inlet and outlet totals are 1.61767 and 1.61683 m3/s (0.052%
relative imbalance). Lab 2 reports eight supply boundaries totalling 1.42906
m3/s and two recirculation fans, but no external outlet-flow boundary, so that
closure remains withheld. The downloaded CSVs do not declare a concentration unit or
map relative seconds to the documented release clock; event calibration,
helium-to-hydrogen equivalence and STEP opening-face normals are also absent.
Consequently this lane is conditional component evidence only: direct vol%
scoring and quantitative H2/LH2 promotion are disabled. Reproduce the current
audit by running `tools/audit_h2safe_dataset.py` against the raw dataset root
and providing a new `--output` path.

Run `python tools/audit_ijhe_submission.py` before circulating a manuscript
snapshot. It checks the claim boundaries, numerical-verification artifact and
public-snapshot safety, while retaining an explicit warning until a matched
time-resolved field event is available.

The draft upload bundle is [the manuscript](MANUSCRIPT_LH2_VALIDATION.md),
[Supplementary Information](IJHE_SUPPLEMENTARY_INFORMATION.md),
[Highlights](IJHE_HIGHLIGHTS.md) and a
[cover-letter draft](IJHE_COVER_LETTER_DRAFT.md).
The [graphical abstract](IJHE_GRAPHICAL_ABSTRACT.svg) and
[figure/table register](docs/ijhe-figure-table-register.md) keep the upload
artifacts reproducible.

An author-neutral assembled draft is generated with
`tools/build_ijhe_upload_bundle.py` in
`outputs/IJHE_UPLOAD_DRAFT_LATEST/`. Its `BUNDLE_STATUS.md` and
`bundle-manifest.json` keep the remaining author and external-evidence gates
visible; the directory must not be treated as a final upload until those gates
are closed.

The separately audited [open-channel H₂ timing boundary](docs/open-channel-h2-boundary.md)
is retained as supporting evidence only; it is not promoted to atmospheric
LH₂ validation because its weather, geometry and common-clock evidence are
incomplete.

The research-only six-flux thermal-moment marcher now passes its preregistered
numerical extension across PRESLHY Trials 11, 12, 22, 23, 24 and 25: all coarse
and refined paths reach their frozen stations, with maximum independent balance
residual `5.05e-11` and maximum terminal flux difference `1.57e-4`. This is a
conservation, solvability and step-convergence result—not a new field-accuracy
score or safety-design validation. See the
[thermal-moment extension results](docs/thermal-moment-flux-extension-results.md).

A subsequent identical-sensor audit for Trials 10 and 23 reaches 4 m with
strong conservation and step convergence under an opt-in complete-radial-core
constitutive check. It slightly improves thermal spread but worsens centre
temperature-deficit amplitude, so the field-acceptance gate fails and neither
the candidate nor the new domain option is promoted. See the
[thermal-moment observation result](docs/thermal-moment-observation-results.md).

See [model status](docs/stage1-results-2026-09-06.md),
[claim grading](docs/claim-grading.md), and
[data/reproduction notes](docs/DATA_AND_REPRODUCTION.md). The current
availability and limitations of the LH2 measurements are summarized in the
[LH2 data-readiness register](docs/lh2-data-readiness.md).

## Field-deployment extension (preliminary)

An opt-in local 2-D wind-plane, conservative semi-finite-volume screening
path now makes the release boundary, flash state, weather, surface, sensors,
obstacle geometry and deterministic uncertainty bounds explicit. A selected
historian event can enter only through an SI CSV map with calibrated P/T/mass
flow channels, a separate phase/flash-state evidence identifier and a retained
file fingerprint; an ordinary pressure trend or nominal volume flow is not
converted into a leak source. Long historian windows can reuse a bounded LH2
property table instead of calling CoolProp for every interval.

This path is suitable only for conditional operational screening and research.
It is not a CFD solver, a building-wake closure, a design-basis calculation or
an approval decision; local LH2 obstacle validation remains outstanding. See
the [deployment roadmap](docs/field-deployment-roadmap-2026-10-05.md) and
[connected-project evidence review](docs/field-evidence-readiness-2026-10-05.md).

For a declared nominal site case, start from the synthetic
[field-case template](docs/field-screening-case.example.json), replace every
evidence ID and value under engineering control, then run:

```bash
degali field-screen field-case.json --output new-field-execution.json
```

The command runs its grid/time refinement by default and records its input
SHA-256, numerical report and fail-safe disposition. `--allow-conditional` is
an explicit accountable-review opt-in for screening only after relevant
source/weather/detector uncertainty has been propagated; it never enables a
design-basis or approval result. The flag requires the strict case's
[`conditional_review` record](docs/field-conditional-review.md), which is
retained in the execution JSON. The nominal template does not include a
historian or a post-flash schedule. A selected SI historian CSV may be added
only through the [quality-gated measured-history mapping](docs/field-measured-history-input.md);
post-flash schedules still retain their dedicated API/provenance gate and must
carry an explicit non-default `source_id` when attached through the Python
field request API.

For a pressure/temperature-only historian, the programmatic
`PressureDrivenMeasuredHistory` path is the only supported source derivation:
`direct_vapour_schedule_from_pressure_driven_history()` recomputes an explicit
orifice boundary interval by interval, and
`run_field_pressure_driven_history_envelope()` propagates its coherent source
corners into the semi-FV screen. Use
`run_field_joint_pressure_driven_history_envelope()` when ambient, weather and
detector corners must be propagated with that source. A pressure trend is never interpreted as a
leak rate; declared opening area/discharge coefficient (including any bounds
that must be crossed by the joint envelope) and source location/direction
corners are fixed before wind-plane projection. Independent event, phase,
clock and calibration evidence remain mandatory.
The strict case-file equivalent is `pressure_driven_history`; its CSV map
contains P/T (and optional liquid fraction) but no mass-flow column, so the
zero source placeholder cannot be mistaken for a measured rate. `field-screen`
and `field-verify` retain the pressure-driven history kind and its CSV digest
through the joint uncertainty/refinement envelope.

Ambient temperature, absolute pressure and air-density bounds may also be
declared at the top level as `ambient_temperature_uncertainty_k`,
`ambient_pressure_uncertainty_pa` and `ambient_air_density_uncertainty_kg_m3`.
They are deterministic corners, not probability intervals; a nominal result
with unresolved ambient bounds remains withheld by the operational gate. The
phase-routing/pool transport adapter carries the same ambient corners through
the phase ledger and clears them only on the corresponding field corner;
measured-history operational runs remain withheld until their time-aligned
history and ambient boundaries are both propagated.
When an ambient-pressure bound is present, the reusable LH2 saturation table
also spans both endpoint saturation temperatures before any corner is run.
Likewise, a fingerprinted validation dataset is evidence identity only: the
field gate remains conditional until a matched observed-vs-model validation
score is generated.
Direct ambient uncertainty objects enforce `K`, `Pa`, and `kg/m3` units before
their corners enter flash, transport, or detector conversion, matching the
strict case-file boundary.
The joint measured-history envelope carries the same ambient temperature and
pressure corners through each history re-flash and applies ambient air-density
corners to detector conversion; a completed joint case clears only the ambient
corner it actually evaluated.

When the case contains relevant bounded source, weather or primary-detector
inputs, use `--uncertainty-envelope`. It runs every deterministic corner and
its refinement, then allows screening only if every corner passes the same
gate. With a quality-gated `measured_history` object, the same flag instead
runs the joint P/T/flow/phase-history, weather and detector envelope:

```bash
degali field-screen field-case.json --uncertainty-envelope --allow-conditional \
  --output new-field-envelope.json
```

The envelope is still a sensitivity set rather than a confidence interval.
For a finite transient release, `scenario.source.duration_s` may itself be a
bounded object, or a numeric duration may be paired with
`duration_uncertainty`; the nominal duration and every lower/upper timing
corner are retained in the source contract. The transport duration must cover
the upper declared release duration. No timing probability or unrecorded
source schedule is inferred.
The source `location_m` may likewise use three SI bounded coordinate objects;
their Cartesian lower/upper corners are propagated through the wind-plane
projection and transport case. Mixed numeric/bounded coordinate arrays are
rejected, and this remains a deterministic geometry sensitivity rather than a
probability distribution.
The same three-component bounded form is available for `direction_m`; every
corner must remain a non-zero vector and is preserved for jet/phase-routing
handoff calculations. The direct scalar screen records these direction
corners even where its post-flash source closure does not use momentum.
An explicit base case may also declare `distributed_vapour_sources` for
already-atmospheric pool/droplet or other ledger branches. Each source carries
its global position, a unit-labelled vertical scalar width, a finite
zero-endpoint rate schedule and evidence ID; bounded position/width values are
expanded into separate source-shape corners and are not treated as a
probability distribution. The operational gate withholds a nominal result
until those corners are propagated.
For categorical meteorological stability uncertainty, declare evidence-backed
alternatives and a diffusivity for every class; the local solver does not
invent a turbulence correlation. The [stability-alternative protocol](docs/field-stability-alternatives.md)
adds those cases to both nominal and historian joint envelopes.

`scenario.weather.direction_deg` is a circular degree interval, rather than a
linear number range. A wind record crossing north can therefore declare
`lower: 350`, `nominal: 0`, `upper: 10`; the envelope evaluates 350° and 10°
and records the short 20° sector. It does not reinterpret that input as a
340° sweep through south or average a changing wind into one wind plane.

Obstacle footprints may be global-axis-aligned or explicitly rotated; both are
reprojected for each wind-direction corner and retained in the report. Linear
dimensions (and an oriented bearing when needed) may be declared as explicit
bounded geometry values with an `uncertainty_evidence_id`. The field envelope
then recomputes every obstacle-mask corner, while the nominal operational
screen remains withheld until those corners are resolved. See the [obstacle
input contract](docs/field-obstacle-input.md) for the exact JSON forms,
required site-coordinate reference and the limits of the conservative 2-D
mask.

`lh2_validation_available: true` is not a free-standing claim. Strict case
files must also carry `validation_evidence` with a SHA-256 fingerprint,
positive row count, common-clock ID, source/weather/obstacle/receptor geometry
IDs and temporal-operator ID. A free-field evidence record cannot qualify an
obstacle case; an unbacked validation flag is rejected before transport.

For an actual fixed-sensor validation replay, `degali field-validate` accepts
the separate `degali.field-validation-input.v1` contract. It checks every
observed row's clock/obstacle IDs and file fingerprint before calculating
matched-sensor MAE, RMSE and bias. See
[field-validation-input.md](docs/field-validation-input.md); no unmeasured
hazard distance is interpolated. If an observed CSV declares an
`observation_kind` column, `lower_bound` rows are treated as censored
one-sided constraints, excluded from symmetric accuracy metrics, and prevent
`qualified` status even when satisfied.

`degali field-audit <directory>` scans an external CSV/JSON/XLSX directory read-only
for the five evidence channels needed to assemble that validation case. It is
readiness triage only: even `candidate_complete` keeps
`promotion_allowed=false` and cannot create `FieldValidationEvidence`; a scan
limit sets `scan_complete=false` and also prevents promotion. The saved record
includes candidate-file `diagnostic_counts` for near-miss categories and
recomputes them on read-back. See
[field-evidence-audit.md](docs/field-evidence-audit.md).

After an accountable reviewer reconciles the five files to one event and common
clock, `FieldEvidenceManifest.from_audit()` can record the selected paths,
SHA-256 digests and explicit source/weather/geometry IDs. New packages should
also provide `sensor_set_id`, a SHA-pinned `sensor_registry_artifact` (sensor
calibration/response metadata), and an accountable `operator_id`; missing
metadata is reported as `evidence_readiness=conditional`. The manifest
rechecks its files and remains `promotion_allowed=false`; see
[field-evidence-manifest.md](docs/field-evidence-manifest.md).
The package layout and minimum CSV/JSON fields are listed in
[field-evidence-intake-template.md](docs/field-evidence-intake-template.md).
`field-validation-input.v1` may carry the same manifest as
`validation_manifest`; the parser then rejects any path, digest or explicit-ID
mismatch between the manifest and `validation_evidence`.
The equivalent CLI boundary is
`degali field-evidence-manifest-create <audit-execution.json>` with five
repeated `--selected-path channel=path` options and explicit event/identity
metadata, plus `--sensor-set-id`, `--operator-id` and optionally
`--sensor-registry-path`; it refuses incomplete audits and existing output
paths. Verification emits machine-readable `conditional` or `withheld`
readiness instead of silently treating missing evidence as valid.

For the public FFI/Spadeadam horizontal-release lane, use
`degali ffi-source-state`. It crosses explicitly declared source, wind and
ambient endpoints through the unchanged free-field jet and evaluates every
reported sensor coordinate. The output is a deterministic residual envelope,
not a confidence interval or site-acceptance record:

```bash
degali ffi-source-state --test 6 --radius 30 \
  --rate-bounds 0.82,0.833,0.85 \
  --wind-bounds 2.2,2.3,2.4 --max-cases 4 \
  --require-complete --output ffi-test6-source-state.json
```

See [the FFI source-state envelope note](docs/ffi-source-state-envelope.md).

Public component-level dynamic, spatial and CFD evidence can be audited with
`tools/audit_usn_carpark.py`, `tools/audit_hytunnel_ps.py` and
`tools/audit_multifidelity_cfd_evidence.py`. The first reproduces a published
4 vol% dynamic sensor operator, the second preserves HyTunnel concentration
and ventilation matrices, and the third pins the MultHyFuel obstructed-jet
CFD context. None silently promotes GH2/confined/plot-only evidence into the
outdoor LH2 headline score; see
[the multi-fidelity intake note](docs/ijhe-multifidelity-evidence-2026-10-10.md).

When a reviewer has an owner-exported Test 6 wind history, the separate
observation replay can be run with
`tools/audit_ffi_test6_transient_receptors.py`. It requires the exact
`time_s,wind_speed_ms,wind_direction_from_deg` CSV contract, a declared
steady-wind table and window, fingerprints the history, and reports compact
true/indicated sensor statistics. It does not add plume storage or infer a
time-varying source; its execution artifact always keeps promotion disabled.
Use the same tool with `--verify <replay.json>` to recheck the history and
reference-table hashes without rerunning the plume.

`degali field-source` validates a separate, already-atmospheric source-rate
CSV. Its strict evidence record fingerprints the file and binds it to a source
boundary and common clock; the final endpoint must be zero. It does not infer a
flash, liquid fraction, pool evaporation, or leak rate from pressure/level data.
See [`field-source-input.md`](docs/field-source-input.md).
The same parsed schedule can be passed to
`run_field_operational_uncertainty_envelope(...,
atmospheric_source_schedule=...)` so each deterministic source corner receives
its own refinement and fail-safe screening decision.
The resulting report also contains `deterministic_sensor_envelope`: per-detector
peak/final/time-average extrema over completed corners, with missing or withheld
traces counted explicitly rather than replaced by zero. These extrema are
deterministic sensitivity bounds, not confidence intervals.
The same field is emitted by the measured-history, source-plus-sensor and
phase/pool operational envelope reports, as well as the standalone
`--sensor-array-envelope` calibration path.
The CLI keeps that source contract in a separate file and composes it only at
the complete-envelope boundary:

```bash
degali field-screen field-case.json --uncertainty-envelope \
  --atmospheric-source-case source-case.json \
  --output field-source-execution.json
```

The source-case path and SHA-256 are retained in the execution input record.

For named batch exports, pass the parsed schedule through
`export_field_screening_batch(...,
atmospheric_source_schedules={"case-name": schedule})`. Those cases are
written as complete operational-uncertainty-envelope reports and keep the
refinement, in-plane receptor and resolved-uncertainty gates enabled; the
batch API rejects attempts to waive them. Conditional-review opt-in requires
one authorization for every named case and is retained in the manifest; the
aggregate can become `conditional_allowed` only after every deterministic
source/weather/sensor corner passes those gates.

For a file-based named batch, use the strict
[`field-batch-input.md`](docs/field-batch-input.md) contract:

```bash
degali field-batch batch.json --output-directory field-batch-output
```

The execution record fingerprints the batch file, every case JSON and every
optional source-schedule JSON. Its top-level
`degali.field-batch-summary.v1` `summary` is derived from the typed
operational decisions and reports case/status counts plus the
fail-safe `all_screening_allowed` flag. Historian and phase-routing cases are kept on
their dedicated envelope paths rather than mixed into this batch.
For a measured-history `field-screen` execution, `field-verify` also reimports
the current historian CSV and compares its SHA-256/event provenance with every
nested report, so changing the CSV without regenerating the execution fails
closed.
Complete bounded `nominal_field`, phase-routing, standalone sensor-array and
joint source-sensor executions are replayed from their serialized numerical
options; the joint source-sensor path also reopens its optional atmospheric
schedule. Only genuinely option-incomplete envelope or batch records remain
integrity/semantic checks without a fabricated replay.
`degali field-verify` also checks serialized semi-FV conservation diagnostics
inside non-replayable envelope and batch reports: per-source schedule residuals
must be finite and consistent with their recorded maximum, and any residual
outside the numerical gate must carry a matching blocked applicability reason
and operational `gate_codes` entry. Compact transport concentration and sensor
extrema/time/density records are also checked, including the invariant that a
withheld detector has warnings and no fabricated result. The serialized
operational decision itself must retain consistent status, allowance,
reason/action, and unique gate-code fields. Source-level injected-mass labels
and totals are cross-checked against the ledger and schedule residual records.
It cannot be combined with the measured-history or phase-routing source paths,
and it is rejected unless every source corner is refined and gated.
The CLI direct-envelope composition accepts only
`post_flash_atmospheric_vapour`; pool or droplet schedules require their
location-bearing distributed/phase-routing handoff and are not placed at the
nozzle by inference.

For an explicitly routed liquid/rainout branch, the [phase-routing transport
handoff](docs/field-phase-routing-transport.md) preserves direct flash and
adds only a conservative time-resolved pool-evaporation ledger as a separate
field source. It withholds missing pool source histories rather than filling
them with a fitted evaporation or wake coefficient. The dedicated
[`degali.field-phase-routing-screening-input.v1` example](docs/field-phase-routing-case.example.json)
can be executed with `degali field-screen ... --uncertainty-envelope`; the CLI
requires the complete phase/pool/weather/surface/sensor envelope for this path.
If evidence-backed bounds are available for phase/pool scalars, the strict
`phase_routing.uncertainty` object propagates their lower/upper corners (with
the upper post-release bound enforced in the transport duration); no phase
range or droplet-population distribution is inferred. The pool launch scalar
width can also carry an explicit bounded `m` interval, which is evaluated as
separate source-shape corners. Droplet population alternatives are accepted
only as complete evidence-backed class corners whose mass fractions sum to one.

To audit DEGALI against supplied SLABx/CFD predictions—or to expose a
steady-versus-transient mismatch—use the separate fixed-sensor comparison
contract. It preserves the two CSV input fingerprints, units, sensor geometry,
source/weather IDs, temporal operator and (when applicable) the obstacle-mask or
wake-closure representation rather than inferring equivalence from model names:

```bash
degali field-compare comparison.json --output new-comparison.json --require-comparable
```

See the [comparison template](docs/field-model-comparison-case.example.json)
and [comparison protocol](docs/field-model-comparison.md). A mismatch retains
numerical rows but withholds a model-selection conclusion; even an aligned
comparison reports only monitored receptors, not an interpolated risk distance.

## Installation

Install Python 3.10 or newer. For LH2 work, install DEGALI with the CoolProp
extra:

```bash
python -m pip install "degali[coolprop]"
```

Check the installed version and command-line interface:

```bash
python -c "import degali; print(degali.__version__)"
degali --help
```

CoolProp is required for LH2 calculations. The traditional DEGADIS
compatibility path can be installed with `python -m pip install degali`.

## Quick manual

### Liquid-hydrogen release without an input deck

The simplest public interface is `degali.lh2.assess`. The example below
represents a horizontal 0.1 kg/s release through a 10 mm orifice, 0.5 m above
ground, from saturated LH2 stored at 6 bar absolute:

```python
from degali.lh2 import assess

result = assess(
    rate=0.10,                  # kg/s
    wind=2.0,                   # m/s at release height
    height=0.50,                # m
    orifice=0.010,              # m; use pool_diameter instead for a pool
    storage_pressure=6.0,       # bar(a), not gauge pressure
    ambient_temperature=288.15, # K
    relative_humidity=65.0,     # percent
    ambient_pressure=101325.0,  # Pa
    max_distance=30.0,          # m
    at_distance=10.0,           # m
)

print(result.report())
for warning in result.warnings:
    print("WARNING:", warning)
```

The principal outputs are:

- `distance_to_lfl`: centreline distance to the hydrogen LFL of 4 mol%, m;
- `distance_to_stoichiometric`: centreline stoichiometric distance, m;
- `lowest_flammable_height`: lowest flammable-gas height, m;
- `regime`: `grounded`, `low`, or `aloft`;
- `trajectory`: NumPy columns `[distance, centre height, mole fraction]`;
- `warnings`: validation-range and applicability warnings. Do not discard them;
- `screening_scope`: `qualified`, `conditional`, or `out_of_scope`.

For an automated screening gate, add `strict_scope=True`; an out-of-range
request raises `ApplicabilityError` instead of returning an extrapolated value.
For explicit source/wind sensitivity, use `assess_envelope(rates=[...],
winds=[...], ...)`. The supplied values are evaluated without fitting a
correction or claiming a statistical confidence interval.

To diagnose a measured centreline point such as the FFI Test 6 residual, use
`assess_observation_envelope(...)`. It returns every declared rate/wind
hypothesis within the requested factor of the observation. Multiple matches
are reported as non-identifiable; DEGALI never chooses a fitted correction.
Sensor-height or arc-maximum values require a separate validated observation
operator and cannot be inferred from this centreline diagnostic.
For exact mast coordinates, `project_lh2_jet_to_sensors(points_m=[...])`
evaluates the Gaussian vertical/lateral profile directly and returns mole
fraction and temperature at each point.

Horizontal releases at a non-zero bearing can be explored with
`run_lh2_yawed_crosswind_research()`. This uses the conserved six-flux yaw
kernel, rejects reverse-axial flow, and is explicitly marked unvalidated; it
does not change the primary `assess()` path.

For a defensible comparison with another consequence model, use the comparison
gate rather than comparing screenshots or separately reduced tables. It
requires identical source, wind vector, receptor operator, averaging window,
phase closure and geometry; a fitted prediction is reported but cannot win a
ranking:

```python
from degali.validation import ComparisonCase, ModelPrediction, compare_models

case = ComparisonCase(
    source_mode="horizontal_jet", source_rate_kg_s=0.12,
    wind_speed_m_s=2.0, wind_direction_rad=0.0, release_height_m=1.0,
    receptor_operator="arc_max_1s_peak", averaging_time_s=1.0,
    phase_closure="measured_throat", geometry="open_horizontal",
)
report = compare_models(
    observations,
    [
        ModelPrediction("degali", degali_values, case),
        ModelPrediction("other_model", other_values, case),
    ],
    case=case,
)
print(report.ranking_allowed, report.results)
```

This protocol does not claim that DEGALI is generally more accurate than
HyRAM, PHAST or EFFECTS. It makes that claim testable when independently
generated predictions are available and otherwise reports that no ranking was
performed.

For a workflow-level safety check, `degali.validation.evaluate_screening()`
combines the result scope with an optional obstacle-geometry screen. It can
approve a qualified clear-path *screening* run, but it always returns
`design_basis_allowed=False` and `approval_allowed=False`.

Specify exactly one source geometry: `orifice` for a pressurised jet or
`pool_diameter` for a pool/evaporation source. `storage_pressure` is bar(a),
while `ambient_pressure` is Pa. All temperatures are K.

### Save the LH2 trajectory

```python
import numpy as np

np.savetxt(
    "lh2_trajectory.csv",
    result.trajectory,
    delimiter=",",
    header="distance_m,height_m,centreline_mole_fraction",
    comments="",
)
```

### Existing DEGADIS input decks

```python
from degali import run_steady, run_transient, run_jet_to_ground

profile, source = run_steady("CASE.INP", backend="legacy")
print(profile.distance_to(0.05))

transient = run_transient("TRANSIENT.INP")
profile, jet, source = run_jet_to_ground("JET.INO", "GROUND.IN")
```

The equivalent command-line workflows are:

```bash
degali steady CASE.INP
degali transient TRANSIENT.INP --snapshot 60 --snapshot 120
degali dose TRANSIENT.INP --at 50 --at 100 --at 200
degali jet JET.INO --bridge GROUND.IN
```

DEGADIS decks are positional files: a missing value changes the meaning of
every field that follows. Their pressure convention also differs from
`assess()`. Preserve the original deck, verify units, and change one field at
a time. Input decks and licensed third-party data are not bundled.

For complete installation, input, output, plotting, parameter-study, deck,
validation-range, and troubleshooting instructions, see the root-level
**[English user guide](USER_GUIDE.md)** or
**[한국어 사용자 가이드](USER_GUIDE_KO.md)**.

## Verification

Fast development checks:

```bash
python -m pytest -m "not slow" -q
```

The original Fortran oracle and raw REDIPHEM, SMEDIS and PRESLHY files are not
redistributed. Tests that need
them skip with an explicit message unless the documented environment variables
are configured. Their provenance, reduction method and aggregate validation
results are documented. See [publication scope](docs/publication-scope.md).

The consolidated [technical reference](docs/technical-reference.md) describes
the governing physics, LH2 extensions, validation process, quantitative
results and remaining model-form limitations.

For the IJHE manuscript track, use the
[submission-readiness record](docs/ijhe-submission-readiness.md) and run
`python tools/audit_ijhe_submission.py` before circulating a draft.

## Repository layout

```text
src/degali/  model and validation code
tests/       self-contained numerical, physical and regression tests
tools/       research and release verification utilities
docs/        derivations, provenance, audits, limitations and results
```

## Scope and safety

The current validated domain does not cover arbitrary equipment conditions,
obstacles, indoor releases, pool spreading, downward/strongly wind-steered
jets, or general transient source behaviour. Always compare safety-critical
results with independent experiments and an accepted consequence-analysis
workflow. See [security and safety reporting](SECURITY.md).

## Citation and license

Citation metadata, including the published concept DOI and the current
version-specific DOI, are provided in [`CITATION.cff`](CITATION.cff). The
README DOI badge uses the concept DOI so it continues to resolve to the latest
archived DEGALI release.

The exact archive for DEGALI `v0.3.0` is
[`10.5281/zenodo.23256538`](https://doi.org/10.5281/zenodo.23256538).

The complete GitHub, Zenodo DOI and PyPI release sequence is documented in the
[publication guide](docs/publication-guide.md).

The DEGALI Python implementation is available under the MIT License. Original
Fortran implementations and external experimental datasets are not bundled or
redistributed.

The chronological development log that previously occupied this front page is
preserved in [docs/README-development-log-2026-09-06.md](docs/README-development-log-2026-09-06.md)
and [CHANGELOG.md](CHANGELOG.md).
