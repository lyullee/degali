# DEGALI validation snapshot — 2026-10-14

This snapshot records the strongest current evidence without promoting
diagnostic or surrogate runs to facility validation.

| Evidence path | Result | Promotion boundary |
|---|---:|---|
| AIJ Case H resolved 3-D obstacle flow, time audit | `accepted`; 3000→6000 step U MAE change 2.39%, W RMSE change 2.92% | obstacle-flow/time convergence only; not hydrogen validation |
| AIJ Case H resolved 3-D grid audit | `withheld`; finest-pair U MAE/W RMSE and lee velocity changes exceed 10% | no grid-independent resolved claim |
| AIJ Case H passive-tracer score | Pearson `0.802`; RMSE `196.1 ppmv`; runtime `9.65 s` | ethylene tracer component evidence only |
| AIJ Case H lateral-domain FV tracer sensitivity (`outputs/aij-case-h-lateral-domain-tracer-comparison-2026-10-14-v1/aij-case-h-lateral-domain-tracer-comparison.json`) | paired 582-point 8–12 s FV/SGS/nozzle-momentum comparison; lateral width `0.30→0.50 m` changes RMSE `1007.8→763.5 ppmv` (`24.2%` reduction) but Pearson `0.743→0.293` | side-boundary physics sensitivity identified; no H₂/facility promotion |
| FFI obstacle/control baseline | 29 sensors, paired O1/O2 trials, 60–120 s windows, SHA-pinned | scoring withheld until geometry/clock/calibration gates clear |
| FFI paired O1 obstacle-effect diagnostic (`outputs/ffi-o1-control-obstacle-delta-2026-10-14.json`) | same fast path on control A and conditional O1: 24 paired receptors, 0.05 m grid, 11.96 s total; predicted Δ MAE/RMSE `0.328/0.537 vol%-points`, signed correlation `−0.098`; mass residual `1.21×10⁻17 kg` | conservation passes but the conditional beam surrogate misses the measured obstacle-effect pattern; geometry/clock/calibration remain withheld |
| FFI FieldEvidenceManifest intake (`outputs/ffi-o1-field-evidence-intake-2026-10-14-v2.json`, `outputs/ffi-o2-field-evidence-intake-2026-10-14.json`) | raw source/sensor CSVs and registry are fingerprinted; source/receptor channels are conditional; weather, obstacle geometry and common clock remain withheld | fail-closed intake only; no five-channel manifest or promotion until event weather, machine-registered beam geometry, clock transform and calibration certificate are supplied |
| FFI machine-registered geometry adapter | `load_registered_ffi_obstacle_geometry()` and both resolved/hybrid runners accept only `degali.ffi-obstacle-geometry.v1` with `registration_status=machine_registered_exact`; bounds and coordinate frame are checked before mesh construction | executable geometry path is ready, but no such machine-registered O1 export is present in the public package, so current FFI runs remain conditional |
| FFI manifest runtime gate | both O1 runners accept `--field-evidence-manifest`; `--require-quantitative-manifest` rejects conditional clock/calibration packages before solving, while default mode records the manifest as a conditional witness | prevents an unverified field package from being reported as quantitative validation |
| FFI resolved startup smoke | max `3.243 vol%`; mass residual `1.06e-22 kg`; runtime `8.17 s` | regularised source smoke, not facility validation |
| FFI O1 conditional replay | 26 sensor traces, max `9.23 vol%`, mass residual `1.21e-17 kg`, runtime `0.63 s`; time-integrated all-sensor MAE `4.15 vol%` (diagnostic), observed ≥4 vol% MAE `7.93 vol%` | rectangular beam proxy; geometry/clock/calibration withheld |
| FFI O1 conditional crossbeam replay | 24 sensor traces after rejecting receptor stencils touching/adjacent to the surrogate, max `9.47 vol%`, mass residual `1.21e-17 kg`, runtime `5.88 s` on a `116×18×16` (`0.05 m`) grid with only `x_max` open; diagnostic MAE/RMSE `2.93/4.89 vol%-points` | two-beam surrogate preserves the reported `36.8%` blockage ratio, but registered beam CAD, clock alignment and calibration remain withheld; no facility score |
| FFI O1 thin-top-beam shape hypothesis | same conservative 24-receptor stencil set, max `9.47 vol%`, mass residual `1.30e-17 kg`, runtime `5.77 s`; diagnostic MAE/RMSE `2.94/4.89 vol%-points` | `0.15 m` upper-beam height inferred from the reported blockage and near-full span; shape hypothesis does not improve the common-sensor diagnostic and remains withheld |
| FFI O1 crossbeam transport sensitivity | controlled exploratory runs with wind `1e-6–0.10 m/s`, diffusivity `0.002–0.01 m²/s`, source width `0.03–0.08 m`; runtime `4.9–6.9 s`; sensor-set-dependent MAE `3.52–4.08 vol%-points` | changing scalar transport parameters changes peak concentration (`9.47–84.39 vol%`) but does not remove the spatial underprediction; this points to missing quiescent/buoyant roof-layer and open-end boundary physics, not a basis for parameter fitting |
| FFI resolved-thermal preflight | rejected before integration: the blockage-preserving conditional beams leave only `0.05 m` lateral fetch, below the resolved solver's two-cell side-fetch contract | no resolved FFI score until beam section/placement is machine-registered or a separately declared lower-blockage sensitivity geometry is used |
| FFI O1 resolved buoyant prefix (`outputs/ffi-o1-resolved-buoyant-prefix-60s-2026-10-14.json`) | 60 s replay on a `58×9×8` (`0.1 m`) grid in `44.28 s`; composition-density buoyancy enabled; 27 sensor histories; mass/energy residuals `6.07e-18 kg` / `2.49e-14 J`; diagnostic MAE/RMSE `3.99/5.75 vol%-points`; maximum 4 vol% volume `0.021 m³` | lower-blockage resolved sensitivity geometry, numerical carrier speed `0.1 m/s`, and two receptors omitted for solid stencil; field/contour artifacts are physics witnesses only, no facility promotion |
| FFI O1 resolved closure/boundary sensitivity (`outputs/ffi-o1-resolved-closure-sensitivity-3s-2026-10-14.json`) | matched 3 s prefix: baseline RMSE `4.27823` and MAPE-floor `63.315%`; Smagorinsky `C_s=0.12`, `Sc_t=Pr_t=0.7` changes RMSE by only `−5.14×10⁻5 vol%-points` (runtime `2.34→3.40 s`); characteristic-vortical outlet changes RMSE by `−8.35×10⁻5`; `w=-60 m/s` source-jet and `σ=0.02 m` variants do not materially reduce error | SGS, outlet and source-shape changes are physics-path sensitivities only; conditional beam CAD, event clock and calibration still prevent promotion |
| FFI O1 resolved source-jet sensitivity | 60 s replay with flow-scaled downward source momentum (reference `w=-20 m/s`) and relaxed Mach gate `0.16` completed in `21.45 s`; diagnostic MAE/RMSE `3.96/5.68 vol%-points`, maximum `9.03 vol%`, 4 vol% volume `0.030 m³`; reducing Gaussian source width `0.08→0.05 m` gives MAE/RMSE `3.99/5.75` and maximum `12.04 vol%` | small numerical improvement only; inferred 4.6 mm-orifice velocity is roughly `0–66 m/s`, so the capped closure is not a calibrated nozzle model; strict Mach `0.15` gate fails over 60 s and the run is not promoted |
| Time-resolved source-jet interface | `ResolvedSourceMomentumHistory3D` now interpolates source momentum on the same event clock as the measured mass-flow schedule; 3 s dynamic prefix passes mass/energy gates with `w` scaled to the public O1 flow trace | API/solver capability validated; no quantitative facility claim until nozzle pressure/temperature/discharge data and a resolved or subgrid jet closure are available |
| Fast-versus-resolved O1 common-sensor comparison | 22 common sensors; resolved-minus-fast mean absolute difference `0.135 vol%-points`, RMSE `0.184`, bias `+0.022`, maximum absolute difference `0.475`; fast/resolved runtime `5.88/21.45 s` | useful interface consistency check only: the two geometries and source closures are not machine-identical, so this is not an accuracy ranking |
| FFI O1 resolved grid prefix | identical 3 s buoyant prefix at `0.1 m` (`58×9×8`, `3.06 s`) and `0.05 m` (`116×18×16`, `135.52 s`); common 27-sensor mean absolute grid change `0.00247 vol%-points`; maximum field `6.72→7.20 vol%` | short-prefix grid change is small for receptor means but runtime grows `44.3×`; full-window fine-grid convergence is not claimed |
| FFI O1 hybrid fast→resolved patch→fast coupling (`outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14.json`) | non-overlapping `1.8/1.6/2.4 m` blocks on a `58×9×8` (`0.1 m`) composite grid; 3,000 incoming and 3,000 outgoing bridge exchanges; active thermal/buoyancy runtime `56.26 s`; diagnostic MAE/RMSE `4.08/5.89 vol%-points`; maximum `7.82 vol%`; maximum partitioned mass residual `1.56×10⁻17 kg`; energy residual `0 J`; continuity residual `1.97×10⁻15 kg m⁻³ s⁻¹`; maximum pressure correction `0.147 Pa`; resolved patch Mach `0.125` | characteristic packet velocity `p_in/m_in` and enthalpy flux callbacks are active at the patch boundary without long-window Mach accumulation. Conditional beam surrogate and sensor score remain withheld |
| FFI hybrid patch bounds | `tools/run_ffi_o1_hybrid_patch_prefix.py` now accepts explicit `--patch-x-min-m/--patch-x-max-m` bounds; a `1.9–3.3 m` one-second isothermal smoke completed in `0.227 s` with `50/50` bridge exchanges and mass residual `2.71×10⁻20 kg` | obstacle-local grid extent is now configurable; the smoke is interface evidence only and does not validate the conditional beam geometry |
| FFI hybrid patch SGS interface | same `1.9–3.3 m` one-second thermal bridge with `C_s=0.12`, `Sc_t=Pr_t=0.7` and molecular-plus-SGS transport; runtime `0.481 s`, `50/50` exchanges, mass residual `2.71×10⁻20 kg`, continuity residual `1.67×10⁻15 kg m⁻³ s⁻¹` | closure selection is now explicit at the fast↔resolved interface; short prefix remains a coupling witness, not a facility accuracy result |
| FFI O1 hybrid 3 s thermal visualization (`outputs/ffi-o1-hybrid-patch-prefix-2026-10-14-field.png`) | assembled `58×9×8` composite maximum H₂ field with the resolved patch obstacle mask and 4 vol% contour; active thermal runtime `1.49 s`; bridge `150/150` exchanges; maximum mass residual `8.13×10⁻20 kg` | visual/ledger witness only; the 3 s prefix is too short for a facility accuracy claim |
| FFI O1 fast/resolved/thermal-hybrid comparison (`outputs/ffi-o1-fast-resolved-hybrid-thermal-comparison-2026-10-14.json`) | common-sensor thermal-hybrid-minus-fast mean absolute difference `0.110 vol%-points`; thermal-hybrid-minus-resolved `0.180`; runtimes `5.88/21.45/56.26 s` | consistency and cost diagnostic only; no ranking because geometry, source closure and momentum treatment differ |
| FFI O1 hybrid temporal convergence (`outputs/ffi-o1-hybrid-temporal-convergence-3s-2026-10-14.json`) | declared isothermal audit at `dt=0.02/0.01 s`; common 27-sensor mean absolute change `2.40×10⁻4 vol%-points` (RMSE `7.67×10⁻4`, max `0.00377`); field RMSE `0.00357 vol%`; runtime ratio `1.92×`; mass residuals `8.13×10⁻20/2.17×10⁻19 kg` | numerical time/interface convergence only; active thermal path is constrained to the D3Q19 physical step and is reported separately |
| FFI O1 hybrid flow-iteration sensitivity (`outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-flow20-2026-10-14.json`) | active-thermal 60-s rerun with resolved flow iterations `40→20`; runtime `56.26→53.66 s`; common 27-sensor mean change `3.83×10⁻7 vol%-points` RMS `9.88×10⁻7`; field RMS `2.25×10⁻8 vol%` | optional speed mode; same conditional geometry and not a validation score |
| FFI O1 proxy-width sensitivity | widths `0.4/0.6/0.8 m`; corrected all-sensor MAE `4.157/4.150/4.150 vol%`; ≥4 vol% MAE `7.96/7.93/7.93 vol%` | width alone does not remove underprediction; no promotion |
| Low-temperature obstacle witness | minimum `274.87 K`, density departure `4.98%`, continuity residual `7.97e-15` | physics-path witness, not measured validation |
| Resolved active-thermal obstacle comparison (`outputs/obstacle-physics-rebuild-2026-10-10/obstacle_thermal_buoyancy_comparison.json`) | D3Q19 resolved obstacle path with conservative H2/enthalpy transport and density buoyancy; active run `4.17 s`; maximum buoyancy acceleration `0.494 m/s²`; minimum `274.87 K`; density departure `4.98%`; lee recirculation cells `6`; mass residual `6.94×10⁻18 kg`; energy residual `1.82×10⁻11 J` | resolved physics witness; flow is not externally calibrated and remains withheld from facility prediction |
| Multifidelity diagnostic | Legacy/Tier-1/Tier-2 runtimes `2.58/0.32/2.15 s`; Tier-2 lee velocity negative | synthetic diagnostic, not experimental validation |
| Quantitative obstacle transport contract | Fast 3-D and resolved 3-D diagnostics explicitly report `nonlocal_scalar_routing=false`; resolved path records local conservative species transport and optional sensible-enthalpy transport | contract evidence; accuracy remains gated by external validation |
| Hybrid transactional coupling smoke | initial/final field mass `0.125/0.125 kg` | interface ledger only; no facility accuracy claim |
| Conservative field-packet coupling contract | shared-clock packet exchange preserves total mass, H2 species mass, sensible enthalpy and momentum to `1e-10`; donor overdraw is rejected | callback/ledger evidence; production outer↔patch wiring still required |
| Solver callback insertion points | fast and resolved 3-D paths accept momentum callback; resolved thermal path accepts conserved species/enthalpy callback before its corrector | hook-level evidence; face-resolved production packet driver still required |
| Face-flux packet extraction | Cartesian interface utility integrates outward total/species mass, enthalpy and momentum over a declared plane and injects incoming packets into a boundary strip; 18 tests pass | utility-level conservation evidence; concurrent outer↔patch execution still not promoted |
| Concurrent callback bridge smoke (`outputs/hybrid-callback-smoke-2026-10-14.json`) | fast outer + resolved patch share 2 steps; `2/2` packets exchanged; fast mass residual `0`, resolved relative mass residual `0`; runtime `0.024 s` | overlapping toy domains and nominal/normal momentum closure; conditional hook validation only |
| Interface packet ledger callbacks | bridge separates state application from packet mass/enthalpy ledger accounting; fast/resolved callbacks accept explicit ledger hooks; bridge/smoke regression passes | prevents hidden packet injection but does not by itself validate facility accuracy |
| Fast TI021 outer screening (`outputs/thorney-ti021-fast-outer-screening-600s-2026-10-14.json`) | 600 s simulated in `15.47 s` on an 88×128×8 (5 m) obstacle grid; HDF5 handoff callback called 2,400 times; maximum mass residual `2.46×10⁻11 kg`; peak R12 proxy `0.549 vol%` | runtime/ledger target met for screening; potential-flow closure, R12 proxy, and one-way handoff mean this is not an H2 facility validation or wake-accuracy claim |
| H2 fast obstacle screen (`outputs/h2-obstacle-fast-screening-600s-2026-10-14.json`) | declared H2 release, 2 m grid, 600 s in `1.73 s`; peak `5.125 vol%`; 4 vol% region `176 m³`; mass residual `1.83×10⁻15 kg` | H2 field/threshold/sensor-output path is operational; release and obstacle are synthetic, external facility manifest remains withheld |
| PRESLHY E3.4 Sand03Wind H2 replay (`outputs/preslhy-e34-sand03wind-h2-fast-obstacle-replay-2026-10-14.json`) | geometry-derived anisotropic pool source plus horizontal-wind/turbulent-mixing surrogate; 329.5 s replay in `3.63 s`; peak `100 vol%` (clipped); 4 vol% region `0.483 m³`; mass residual `8.39×10⁻13 kg`; overall MAE/RMSE `2.332/4.588 vol%-points`, MAPE floor-1 `62.4%` | dynamic score remains withheld: source is liquid-mass proxy, report obstacle is not registered, calibration is specification-only, and buoyant wake is unresolved |
| PRESLHY E3.4 Concrete03Wind H2 replay (`outputs/preslhy-e34-concrete03wind-h2-fast-obstacle-replay-2026-10-14.json`) | same frozen source/geometry policy; 175.0 s replay in `5.98 s`; adaptive accepted step `0.01 s`; peak `96.18 vol%`; 4 vol% region `0.410 m³`; mass residual `2.26×10⁻13 kg`; overall MAE/RMSE `2.532/3.824 vol%-points`, MAPE floor-1 `63.1%` | independent wind/material case remains conditional/withheld for the same source-flux, obstacle-registration and calibration gates |
| PRESLHY E3.4 Gravel03Wind H2 replay (`outputs/preslhy-e34-gravel03wind-h2-fast-obstacle-replay-2026-10-14.json`) | same policy; 220.5 s replay in `1.90 s`; accepted step `0.02 s`; peak `83.60 vol%`; 4 vol% region `0.372 m³`; mass residual `2.58×10⁻13 kg`; overall MAE/RMSE `1.188/2.159 vol%-points`, MAPE floor-1 `43.9%` | independent wind/material case remains conditional/withheld; no facility promotion |
| Quantitative manifest clock gate | common-clock file now requires origin, synchronization method, offset, uncertainty and covered channels; legacy ID-only clocks remain conditional | prevents transient scores from using an unverified shared label |
| TI021 local resolved patch, 2.5 m cells, 120×120×20 m, 2 s | 18,432 cells; 26 in-domain sensors; all physical conservation/continuity gates passed; wall time `21.21 s`; linear 600 s estimate `1.76 h` | accepted as a local-patch qualification only; one-way handoff export, no dynamic outer/patch coupling, and not a full-event score |

## Current promotion rule

The next quantitative hydrogen result must use a common-clock
`FieldEvidenceManifest`, registered obstacle geometry, source history, sensor
calibration metadata, and the resolved or coupled solver path. Until those
inputs are complete, outputs remain `conditional` or `withheld` by design.

## Latest FFI replay clarification

The O1 conditional replay now computes both predicted and reference 4 vol%
exceedance durations from the full 60–120 s time histories.  The reference
duration is a diagnostic description of the measured record; it is not a
validation score because the published obstacle assembly, source/sensor clock
alignment, and calibration certificate remain unresolved.  The corrected
time-integrated diagnostic errors are:

- all 26 paired receptors: MAE `4.1497 vol%`, RMSE `5.8364 vol%`, bias
  `−3.7950 vol%`, MAPE with a `0.1 vol%` floor `151.8%`;
- observed receptors at or above `4 vol%`: 13 receptors, MAE/RMSE
  `7.9331/8.2396 vol%`, bias `−7.9331 vol%`.

These values expose the current proxy's underprediction and are retained as a
failure-diagnostic input for the next geometry/clock reconstruction step, not
as evidence of facility accuracy.

The follow-up sensitivity sweep is stored in the corresponding
`outputs/ffi-o1-crossbeam-sensitivity-*-2026-10-14.json` records.  It varied
carrier speed, scalar diffusivity and source width without changing the
crossbeam geometry.  Peak values moved substantially, but the paired sensor
MAE stayed in the `3.52–4.08 vol%-points` range.  In particular, sensors above
the ceiling source remained high in the measurements while downstream roof
receptors stayed low in the fast potential-flow replay.  This is a diagnostic
of the tier's missing quiescent/buoyant roof-layer transport and open-end
boundary representation; it is not evidence that more parameter tuning will
close the gap.

## TI021 runtime reduction and remaining coupling work

The first reduced-cost resolved run uses a `120 × 120 × 20 m` local patch at
`2.5 m` spacing (`48 × 48 × 8 = 18,432` cells).  It completes a 2 s prefix in
`21.21 s`, compared with `54.28 s` for the `200 × 200 × 20 m` patch and
`106.24 s` for the `200 × 200 × 40 m` patch.  The reduction is real, but a
linear extrapolation still gives about `1.76 h` for 600 s and therefore does
not meet a 30 s facility-screening target.

The run is now labelled `local_patch_sensor_subset` with
`coupling_mode=one_way_handoff_frame_export` and
`dynamic_outer_patch_coupling_applied=false`.  This prevents the output from
being mistaken for a coupled full-event prediction.  The next implementation
step is a shared-clock outer/patch exchange of species mass, momentum and
sensible enthalpy at the patch boundary; only after that exchange passes its
conservation ledger can a 600 s nested score be attempted.

## Fast outer screening path (new)

The runtime bottleneck now has a usable conditional path.  The fast outer
solver reads the durable TI021 HDF5 handoff lazily, remaps the inner plane to
the outer boundary, and applies the positive species mass flux directly to the
finite-volume x-min face, rather than overwriting a concentration field after
the solve.  Interface injection is included in the same mass ledger as source
and outflow terms.  On the 5 m, 88×128×8 grid the 600 s run completed in
15.47 s with maximum residual `2.46×10⁻11 kg`, meeting the screening runtime
target by a wide margin.

This closes the requested fast-versus-resolved branch: use the fast path for
geometry/obstacle triage and threshold maps, and reserve the variable-density
resolved path for cases that pass triage.  The fast result remains conditional:
its flow is inviscid potential flow, it has no separation/recirculation or
buoyancy feedback, and the current TI021 source is an R12 handoff rather than
an H2 calibration dataset.  Those limitations are serialized in the output
diagnostics and prevent accidental promotion.

The H2 branch also writes a compact maximum-field artifact and a contour image
(`outputs/h2-obstacle-fast-screening-600s-2026-10-14.png`) showing the obstacle
footprint and the 4 vol% isoline.  Full time-resolved fields remain opt-in so a
600 s run does not create multi-gigabyte artifacts by default.

The PRESLHY replay now exercises the same H2 output path against a real public
workbook.  The source window is selected from corrected liquid mass only; the
wind speed is projected to the measured horizontal sonic magnitude; and the
H2 channels are loaded only after the solver configuration is frozen.  The
pool footprint is represented by a cell-integrated anisotropic source, with a
  declared 0.20 m near-field mixing-layer surrogate and a transparent
wind/length-scale turbulent-diffusivity rule.  This reduced the Sand03Wind
replay's overall sensor MAE from `102.18` to `2.33 vol%-points` compared with
the original point-like source, while mass conservation remained below
`9×10⁻13 kg`.  Repeating the same frozen policy without changing parameters
gave `2.53` vol%-points for Concrete03Wind and `1.19` vol%-points for
Gravel03Wind.  The Concrete case triggered the explicit scalar-stability
safeguard, halving the step to `0.01 s` rather than accepting an unstable
result.  These remain conditional replays, not a promotion: the source flux
is still a liquid-mass proxy, the report obstacle is not registered, and
buoyant/separated wake physics and calibration certificates are absent.
