# Public obstacle-evidence triage

The following sources were checked for independent obstacle or facility evidence.
They are kept separate by physical relevance; no surrogate or combustion result
is silently promoted to an unignited H2 concentration validation.

| Source | What is available | DEGALI use | Gate |
|---|---|---|---|
| [USN/FFI open-channel H₂ experiments and companion LES](https://doi.org/10.1016/j.ijhydene.2024.10.038), [public data DOI 10.23642/USN.26117989](https://dataverse.no/dataset.xhtml?persistentId=doi%3A10.23642%2FUSN.26117989), and [experimental article](https://doi.org/10.1016/j.jlp.2025.105669) | 29 H₂ sensors and separate source/sensor clocks; 5.8 × 0.9 × 0.8 m channel; ceiling source at y=0.5 m. The article reports XEN-5320 accuracy as ±1%; the manufacturer datasheet gives a preliminary 1% full-scale inaccuracy (equivalent to ±1 vol%-point only if applied to the 0–100 vol% range); event-specific calibration certificates are absent. The article Fig. 2(b) identifies a wooden cross-beam assembly; configuration 1 has two top cross beams and reported blockage 0.37 (36.8% in the companion LES paper); configurations 2/3 add beams (0.45/0.53). The companion CFD is OpenFOAM 10 `buoyantReactingFoam` LES with WALE; obstacle mesh 542,466 cells; O1/O2 sensor-mean scatter and FB/NMSE/MG/VG summaries use 29 positions and a 60–120 s averaging window. These CFD outputs are article/figure-level, not a machine-readable field archive. The obstacle prose labels both 0.8 m width and 0.75 m height on z, and does not specify individual beam sections or exact registration. Treating 0.8 × 0.75 m as a solid transverse rectangle would imply 83.3%, but that is only a rejected simplification—not evidence against the reported assembly blockage. | Highest-priority direct H₂ obstacle benchmark so far: repeat histories, source-flow traces, obstacle/no-obstacle controls, and a published independent LES comparator. Sensor coordinates and eight control/obstacle CSV identities are machine-registered and SHA-256 verified. | conditional: do not instantiate quantitative beam geometry until individual sections, x placement, and y-position datum are resolved. The confined ceiling-jet/channel case is not an outdoor LH₂ plume. No DEGALI-vs-FFI score yet. |
| [ExJet obstructed free-jet model-validation paper](https://hysafe.info/uploads/papers/2023/252.pdf), HySafe 2023, ID252 | 38.0±0.1 barg and 8±1 °C are reported; array described as 1.5-m-long, 0.1-m-diameter cylinders spaced approximately 0.35 m vertically; centreline and radial H2 molar-fraction profiles at 1.4, 2.5 and 4.4 m are plotted. The obstructed direct-orifice diameter is not reported in this source. Table 1 gives four model-specific pseudo-source conditions, not a measured source truth. | Best external obstacle-array concentration-profile benchmark; source-term spread is itself a sensitivity axis. Paper states concentration over-prediction grows downstream for some models and reports geometry-related velocity dips/asymmetry. | conditional: only plot-level observations and schematic geometry; no common-clock history or measured direct source flow. |
| [DataverseNO replication dataset](https://doi.org/10.18710/WSKBIJ), DOI `10.18710/WSKBIJ` | API-verified 52-file set: 51 experiment rows, 39 rows with complete pressure columns, 4.8/9.4 mm nozzles, 0.36–221.52 bar reservoir pressure, four pressure sensors; high-speed `.cine` and raw text files | Source/pressure/obstacle-timing component benchmark; not a direct unignited concentration-field validation because no H₂ concentration channel or machine-registered obstacle geometry is supplied | component-ready; concentration promotion withheld |
| [ELVHYS WP4.2 transfer-connection-space dataset](https://dataverse.no/dataset.xhtml?persistentId=doi%3A10.18710%2FJXJP0H), DOI `10.18710/JXJP0H` | Public CC0 dataset of 48 tests in a 1 m³ LH₂ transfer enclosure. Local subset pins tests 29 and 31: same nominal 2 barg/25.4 mm release, differing vent restrictions (four VAV plates vs three blanks and one VAV plate), 16 H₂ channels at 20 Hz and 500 Hz, release-rig/nozzle pressure, and ventilation flow. README declares Xensor sensors ±1 vol% over 0–100 vol% and five-point calibration. | Fast independent component benchmark for cryogenic-H₂ sensor histories, pressure coupling, and restricted-vent sensitivity. The 500 Hz test gives strong timing resolution. | conditional: no measured H₂ mass-flow channel in the imported FLMT files (they contain fan flow), detailed VAV geometry/report registration and event-specific calibration certificates are not included; confined vent restriction is not an open-field obstacle-wake validation. No DEGALI score yet. |
| [OEDI/NREL H2SAFE dataset](https://data.openei.org/submissions/8799), DOI `10.7799/17118570` | SHA-256-verified archive with CSV concentration histories, 3-D sensor positions, release coordinates, STEP/STL laboratory geometry and HVAC metadata; helium surrogate | Strong facility geometry/sensor/time-series package for transport and sensor-operator testing; not hydrogen density validation | conditional surrogate |
| [GASFLOW-MPI parking CFD dataset, Parts I–IV](https://doi.org/10.5281/zenodo.18740541), [II](https://doi.org/10.5281/zenodo.18745400), [III](https://doi.org/10.5281/zenodo.18754599), [IV](https://doi.org/10.5281/zenodo.18754802) | The records describe a 50 × 30 × 3 m garage, 6 × 6 m column grid, ceiling jet fans, 12 source locations × 5 release rates (1/30/50/100/150 g/s) × 3 ventilation rates (ACH 3/6/10 h⁻¹), and 60 s of 1 s concentration fields. Described solver: compressible Navier–Stokes + k–ε + Boussinesq buoyancy + H₂ species transport; 288,000 cells, with <3.6% deviation from a fine-mesh reference claimed. These are CFD simulations, not experiments. | Potential independent CFD-to-CFD comparison of transient obstacle bypass, ventilation and interpolated sensor histories. Useful only as a numerical comparator, not independent experimental validation. | **unavailable**: all four Zenodo parts mark files embargoed until 2027-05-01; metadata only, no mesh, boundaries, fields or score can presently be checked. |
| [Water-flume dense-fluid cubic-obstacle experiments](https://doi.org/10.3390/fluids10080194), Fluids 2025 | Center-plane concentration fields and downstream wake time series for a Thorney-Island-like cubic obstacle | Neutral/dense-gas wake physics and wake-region regression test; not a hydrogen sensor validation | conditional neutral-fluid component |
| [BOS obstacle-shape H2 concentration study](https://doi.org/10.1016/j.est.2024.110764), Journal of Energy Storage 83 (2024) 110764 | Hydrogen concentration distributions measured around flat-plate, sidewall and dome obstacles; reported Gaussian-like/W/W-like radial signatures and 3.81% BOS inversion error | Shape-signature and contour-regression target for front/side/roof transport; useful for testing whether obstacle geometry changes concentration topology as expected | withheld until raw BOS coordinates, geometry, clock and calibration are machine-registered |

## Immediate import priority

1. Resolve the USN/FFI obstacle projection and placement convention against the
   published figure and blockage ratio, then add a closed-rear, zero-ambient-
   flow source boundary before comparing frozen paired histories and companion
   LES without fitting observations. The resolved 3-D finite-volume path now
   has an explicit zero-ambient-wind mode with a separate numerical velocity
   scale and projected zero inlet; a one-step source/buoyancy smoke passes, but
   it is not yet a channel-scale FFI run. Its `all_no_slip_tunnel` wall option
   alone still does not define the closed-rear/open-end pressure boundary.
2. Obtain the ExJet obstructed concentration profile tables and machine-register
   the cylinder-array geometry; create a new manifest with source, geometry,
   release state, receptors and common clock.
3. Use the locally hash-pinned DataverseNO README/Excel package for
   source/pressure/obstacle-timing checks. The audit confirms 51 experiment
   rows but only 39 complete pressure rows; use raw `.txt`/`.cine` files only
   after selecting an event and registering its geometry. Direct concentration
   promotion remains withheld.
4. Use the locally verified H2SAFE archive as a sensor-operator and
   facility-geometry test with an explicit helium-surrogate label. Its current
   audit reports five events, verified archive members and geometry bounds, but
   withholds quantitative H2/vol% promotion because concentration units,
   calibration and release-clock/boundary mappings are incomplete. Rechecking
   the current OEDI catalog description confirms helium and time-resolved
   concentration histories, but it does not declare the trace unit or provide
   a calibration certificate. Do not infer vol% from the numeric range of the
   CSV values.
5. Keep GASFLOW-MPI in the future-CFD comparison queue because the record is
embargoed.
6. Use the hash-pinned ELVHYS WP4.2 tests 29/31 for a separate LH₂ confined-
   enclosure/vent-restriction component comparison. Do not label this an open-
   obstacle wake benchmark, and do not infer release mass flow from fan flow or
   fit it to the H₂ sensor histories. The local fail-closed audit is
   `outputs/elvhys-tcs-component-audit-2026-10-12.json` and raw files live in
   the ignored `reference/hydrogen/elvhys_tcs_wp4_2/` directory.

## Newly registered FFI matched cases

The primary CFD article's Table 2 labels O1 as `test-02-1a/b` (0.091 g/s)
and O2 as `test-05-1a/b` (0.45 g/s); the local public README maps these to
T00027/T00028 and T00029/T00030. Its unobstructed controls are A
(`test-02-0a/b`, 0.086 g/s; T00012/T00013) and E (`test-05-0a/b`,
0.48 g/s; T00003/T00004). All eight local files have 29 sensor channels.

The existing hash-pinned 4 vol% timing audit gives a descriptive signal that
the obstacle matters: the median of each run's sensor-peak median is 9.17 vs
10.58 vol% (A vs O1), and 8.96 vs 10.33 vol% (E vs O2). Median sensor arrival
is 54.1 vs 54.6 s and 52.1 vs 53.4 s, respectively. These are summaries of
two runs per condition, not sensor-paired uncertainty estimates or DEGALI
errors. The flow classes are close but not identical; do not interpret the
peak shifts as a controlled causal effect until per-sensor pairing and source
histories are checked. This makes FFI a concrete next transient obstacle
benchmark rather than a currently passed validation.

`degali.validation.open_channel_h2_obstacle.load_open_channel_h2_obstacle_benchmark`
now consumes the registered coordinates and the eight raw CSVs together. It
checks all sensor IDs against the 29 article receptors, verifies each CSV's
SHA-256, and returns source and sensor clocks independently. The local intake
verified 8/8 mapped files; the result correctly remains
`geometry_quantitative_allowed=false` and `degali_score_allowed=false`.
Synthetic tamper tests verify that a changed CSV hash fails closed:

### Hash-pinned 60--120 s paired sensor baseline (2026-10-13)

`tools/compare_ffi_obstacle_control_window.py` extracts the published
60--120 s averaging window from all eight control/obstacle CSVs, interpolating
the exact window edges and calculating time-weighted sensor means plus the
source-flow integral on the same elapsed-time interval. The resulting
`outputs/ffi-obstacle-control-baseline-2026-10-13-v1/ffi-obstacle-control-window-60-120s.json`
contains all 29 receptor values and SHA-256 digests for each raw file. All four
paired source records cover 60--120 s; they do not cover the later end of each
complete sensor record, so this is a window-specific baseline, not a full
transient replay.

In raw measured sensor means, O1 minus A is +0.292 and +0.382 vol%-points for
the two repeats; O2 minus E is +1.067 and +0.185 vol%-points. Mean source rates
also differ within those pairs (0.0813/0.0883, 0.0808/0.0854, 0.3947/0.4365,
and 0.4797/0.4381 g/s for control/obstacle). The artifact retains both
unadjusted sensor deltas and an explicitly non-validated linear source-rate
normalization sensitivity; neither is a causal obstacle-only effect. Sensor
clock drift, calibration certificates, exact beam geometry/registration, and
the expanded open-end boundary remain unresolved. It explicitly keeps
`promotion_allowed=false` and contains no DEGALI predictions. Eight targeted
tests pass, including rejection of a requested window beyond source-history
coverage.

Version 2 adds descriptive summaries by registered sensor height/location
class (top, mid-high, mid-low, bottom), without pooling those groups as
independent uncertainty samples. Raw mean top-sensor obstacle-minus-control
differences are +0.666/+0.843 vol%-points for O1 repeats; O2 gives
+1.315/+0.696. O2 repeat 1 also has +1.685 at bottom and +0.683 at mid-low,
while O2 repeat 2 has -0.743 at mid-high. This repeat-to-repeat spatial
variation is a useful wake/stratification target for the model, but not a
stable causal signature or grounds to tune geometry. The v2 artifact is
`outputs/ffi-obstacle-control-baseline-2026-10-13-v2/ffi-obstacle-control-window-60-120s.json`.

### FFI control A direct-H2 3-D startup smoke (2026-10-13)

The first v1 direct-nozzle smearing trial exposed a poor source mapping: a
0.1 m Gaussian on the nozzle plane reduced the 61.7 m/s inferred port speed to
a 0.0040 m/s maximum after 0.01 s. It is retained as historical evidence at
`outputs/ffi-control-transient-smoke-2026-10-13-v1/`, not as the current
recipe.

The runner now connects FFI A to the existing conservative ambient-pressure
H2 low-Mach handoff. Given the article's 0.086 g/s control source rate and
4.6 mm port diameter, the inferred nozzle speed is 61.7 m/s. At a 0.30 Mach
limit and 10% mixture-density-departure limit, the handoff model mixes in
0.01027 kg/s of ambient air and places the source plane 0.459 m below the
ceiling; its equivalent diameter is 0.154 m, mixture velocity 0.513 m/s, and
mass/momentum/energy closure residuals are zero. This is a modelled plane, not
a measured station. At 0.025 m, the equivalent diameter spans 6.16 cells.

The v2 solver smoke advances 10 steps over 0.01 s in 72.54 s wall-clock; v3
uses 2 steps over the same interval in 23.41 s. Both are 232×36×32, unobstructed
control-only runs. The v3 maximum H2 field is 3.243 vol%, mass residual
`1.06e-22 kg`, energy residual zero, and continuity residual
`6.46e-12 kg m-3 s-1`. These two smoke settings are not a time-convergence
pair: they only demonstrate the source-handoff and species/EOS/pressure path.
The full v3 fields and diagnostics are at
`outputs/ffi-control-transient-smoke-2026-10-13-v3/`; v2's 10-step record is
preserved at `outputs/ffi-control-transient-smoke-2026-10-13-v2/`.

The v3 diagnostic still flags the one-cell-thick axial Gaussian as
under-resolved, uses a constant-property ideal-gas mixture, and has no obstacle.
It stops 60 s before the paper's 60--120 s averaging window, so no sensor score
or facility validation is claimed. The 23.4 s for 0.01 s simulated also makes
a direct fine-grid full-window run expensive; next work should test a valid
faster screening tier or streamline the pressure projection before attempting
that horizon. FFI beam geometry remains an independent unresolved obstacle gate.

The current Tier-1 3-D potential-flow model is **not** that shortcut for FFI:
it requires positive imposed wind and its declared scope excludes buoyancy,
separation, recirculation and turbulent mixing. A zero-wind FFI channel must
not be forced through it merely to obtain a faster 60--120 s contour. A focused
regression test now records that the fast-tier configuration rejects the
quiescent FFI condition; this case remains assigned to the resolved tier until
a buoyancy-capable fast model is independently established.

```powershell
$env:PYTHONPATH = 'src'
.\.venv\Scripts\python.exe -m pytest -q tests\test_open_channel_h2_obstacle_registry.py
```

5. Treat the BOS study as a literature-aggregate shape constraint only. Its
reported 3.81% inversion error is an observation-processing error, not a DEGALI
prediction error; do not use it to claim field accuracy.

The current repository therefore has a concrete external benchmark path, but
does not claim that any of these sources already closes the obstacle-H2 field
validation gate.

## FFI O1 crossbeam conditional replay (2026-10-14)

The public FFI/USN package and its companion article describe the O1 obstacle
as two top cross beams and report a 36.8% blockage ratio.  The repository now
constructs a deliberately conditional two-beam surrogate from that published
ratio and routes it through the fast 3-D obstacle solver.  The registry
coordinates are transformed explicitly into the solver's streamwise/lateral/
vertical frame; receptors falling inside the surrogate are omitted rather than
silently sampled through solid cells.

The reproducible run is
`outputs/ffi-o1-crossbeam-conditional-replay-2026-10-14.json`.  At 0.05 m
spacing it uses a 116×18×16 grid, completes in 5.88 s with only the physical
downstream `x_max` face open, and preserves mass to
`1.21×10^-17 kg`.  Five receptors whose interpolation stencils touch the
conditional beam are omitted, leaving 24 paired traces.  Its descriptive
time-window errors are MAE/RMSE `2.93/4.89 vol%-points`.  These values are retained as
diagnostics showing the effect of a physically shaped obstacle surrogate, not
as a claim that the solver matches the facility.

A second conditional geometry hypothesis,
`outputs/ffi-o1-thin-top-crossbeam-conditional-replay-2026-10-14.json`, uses
two near-ceiling beams of inferred height `0.15 m` and near-full lateral span.
It preserves the same reported blockage ratio and gives MAE/RMSE
`2.94/4.89 vol%-points` on the same 24-receptor stencil set.  The negligible
difference is evidence against selecting a beam section by diagnostic error
alone; exact beam CAD remains necessary.

This branch remains `conditional_geometry_replay_not_validation_score` because
the public download does not contain machine-registered beam CAD and datum,
the source and sensor event clocks are not certified as synchronized, and the
sensor calibration certificates are not included in the imported package.
The reported 36.8% blockage is therefore preserved as a traceable constraint,
not treated as sufficient geometry evidence.  The next promotion step is to
obtain the beam assembly coordinates and clock/calibration manifest, then rerun
the same script without changing numerical parameters.

A small, explicitly exploratory sensitivity set is also retained under
`outputs/ffi-o1-crossbeam-sensitivity-*-2026-10-14.json`.  Varying carrier
speed, diffusivity and source width changed the peak field but left the paired
sensor MAE between `3.52` and `4.08 vol%-points`.  The persistent pattern—high
measured roof concentrations upstream and over the obstacle with a much lower
fast-tier downstream field—cannot be repaired by treating those scalar knobs
as calibration parameters.  It identifies the next physics task: a
quiescent/open-end formulation with buoyant roof-layer transport and an
obstacle wake closure.  Until that path is independently checked, the FFI
replay remains a fail-closed diagnostic rather than a promoted score.

A resolved-thermal preflight was also attempted on the same 36.8%-preserving
surrogate.  The resolved solver rejected it before time integration because
its conservative fetch contract requires at least two lateral cell widths of
fluid on each side, while the conditional beams use a 0.05 m edge margin to
preserve the published blockage.  This is a useful fail-closed result: the
surrogate cannot be silently reused as resolved CFD geometry.  A resolved FFI
run therefore requires either the actual beam section/placement or a separately
declared, lower-blockage resolved sensitivity geometry.

The follow-on sensitivity path is now executable in
`tools/run_ffi_o1_resolved_buoyant_prefix.py`.  With a separately declared
0.21 m edge-margin placement, a 0.1 m resolved grid, the measured O1 source
window and composition-density buoyancy, the 60 s replay completed in 44.28 s.
It produced 27 boundary-aware sensor histories, mass residual
`6.07×10^-18 kg`, energy residual `2.49×10^-14 J`, diagnostic MAE/RMSE
`3.99/5.75 vol%-points`, and a maximum 4 vol% volume of `0.021 m³`.  The
maximum-field NPZ and center-plane contour are retained beside the JSON.  This
does not close the FFI validation gate: the geometry is intentionally
lower-blockage sensitivity geometry, the positive carrier speed is a numerical
condition rather than a measured wind, and two receptors were excluded because
their interpolation stencils intersect the conditional solid.  It does,
however, establish the requested composition-buoyant resolved path and its
conservation ledger on the public source history.

It produced 27 boundary-aware sensor histories, mass residual
`6.07×10^-18 kg`, energy residual `2.49×10^-14 J`, and diagnostic
MAE/RMSE `3.99/5.75 vol%-points`.  This does not close the FFI validation gate:
the geometry is intentionally lower-blockage sensitivity geometry, the
positive carrier speed is a numerical condition rather than a measured wind,
and two receptors were excluded because their interpolation stencils intersect
the conditional solid.  It does, however, establish the requested
composition-buoyant resolved path and its conservation ledger on the public
source history.

The resolved prefix also has a first grid check.  The same 3 s physical prefix
at 0.1 m (`58×9×8`) and 0.05 m (`116×18×16`) gives a common 27-sensor mean
absolute change of `0.00247 vol%-points`; the maximum field changes from
`6.72` to `7.20 vol%`.  Runtime rises from `3.06` to `135.52 s` (44.3×), so
the result supports receptor-mean grid stability for this short prefix but not
a full-window fine-grid production run.  The artifact is
`outputs/ffi-o1-resolved-buoyant-prefix-3s-dx005-2026-10-14.json`.

The resolved path was then used to isolate the missing source physics.  The
FFI paper reports a 4.6 mm ceiling orifice; applying the O1 source-flow range
to an ideal-gas H₂ density implies a diagnostic jet-speed envelope of roughly
0–66 m/s.  The previous resolved replay injected the correct mass history but
zero source momentum.  A bounded sensitivity with a flow-scaled downward
source-momentum closure (reference `w=-20 m/s`; finite-volume projection, Mach
gate relaxed from 0.15 to 0.16) completed the 60–120 s window in 21.45 s and changed the
diagnostic MAE/RMSE from `3.99/5.75` to `3.96/5.68 vol%-points`.  A narrower
Gaussian source (`σ=0.05 m` instead of `0.08 m`) changed MAE/RMSE only to
`3.99/5.75` while raising the maximum field to `12.04 vol%`.  These are
physics sensitivities, not validation scores: the machine-registered nozzle
jet is unresolved, the strict Mach-0.15 gate fails for the 60 s momentum run,
and the lower-blockage beam geometry remains conditional.  Artifacts are
`outputs/ffi-o1-resolved-buoyant-prefix-60s-dynamic-momz20-gate016-2026-10-14.json`
and `outputs/ffi-o1-resolved-buoyant-prefix-60s-sigma005-2026-10-14.json`.

For an explicit interface check, `tools/compare_ffi_fast_resolved.py` pairs the
conditional fast crossbeam report with the flow-scaled resolved report on the
common sensor labels.  There are 22 common sensors; resolved-minus-fast mean
absolute difference is `0.135 vol%-points` (RMSE `0.184`, bias `+0.022`, maximum
absolute difference `0.475`), with runtimes `5.88/21.45 s`.  This is a useful
cross-path consistency diagnostic, not an accuracy ranking, because the two
geometries and source closures are not machine-identical.  The reproducible
summary is `outputs/ffi-o1-fast-resolved-comparison-2026-10-14.json`.

### FFI O1 non-overlapping hybrid patch execution (2026-10-14)

The requested “fast outside / resolved only near the obstacle” path is now
executed rather than left as a callback smoke.  `tools/run_ffi_o1_hybrid_patch_prefix.py`
advances three non-overlapping streamwise blocks: fast `0–1.8 m`, resolved
conditional crossbeam patch `1.8–3.4 m`, and fast `3.4–5.8 m`.  The two
`SynchronizedSolverBridge` instances exchange boundary species packets on the
common `0.02 s` clock.  The active thermal/buoyancy 60 s O1 window completed
in `56.26 s` with 3,000
incoming and 3,000 outgoing exchanges, composite maximum `7.82 vol%`, and
partitioned mass residual `1.56×10⁻17 kg`.  The diagnostic sensor MAE/RMSE was
`4.07/5.89 vol%-points`; it is not a validation score because the beam CAD and
clock/calibration gates remain unresolved.  The characteristic-momentum report
is `outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14.json`, with the
assembled field in `outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14-fields.npz`
and its contour image in
`outputs/ffi-o1-hybrid-patch-60s-thermal-characteristic-2026-10-14-field.png`.

The three-way diagnostic in `outputs/ffi-o1-fast-resolved-hybrid-thermal-comparison-2026-10-14.json`
gives thermal-hybrid-minus-fast mean absolute difference `0.110 vol%-points`
and thermal-hybrid-minus-resolved `0.180` on their common
receptors.  These values are useful for regression/interface monitoring only;
they do not establish that the hybrid model is more accurate than either
single-fidelity path.

The patch's standalone resolved carrier remained at maximum lattice Mach
`0.125`.  The bridge now maps incoming characteristic velocity as
`u_in=p_in/m_in` to the patch boundary cells and normal face, eliminating the
previous direct cell-mass impulse that accumulated above Mach `0.30` over 60 s.
The successful conservative species-mass/enthalpy coupling and the new
momentum mapping remain conditional because the beam geometry and carrier boundary closure are
not machine-registered facility data.  The same report records zero thermal
ledger residual, `1.97×10⁻15 kg m⁻³ s⁻¹` continuity residual, and maximum
pressure correction `0.147 Pa` for the resolved patch.

### Hybrid time-step audit

Because the active thermal resolved solver requires its D3Q19 physical step,
time refinement is audited separately in an explicitly isothermal mode.  The
same 3 s hybrid prefix at `dt=0.02` and `0.01 s` changes the 27-sensor mean by
`2.40×10⁻4 vol%-points` (RMSE `7.67×10⁻4`, maximum `0.00377`) and the maximum
field by RMSE `0.00357 vol%`; runtime increases `1.92×`.  Both mass ledgers
remain below `2.2×10⁻19 kg`.  This is numerical/interface convergence, not a
thermal validation score.  The reproducible record is
`outputs/ffi-o1-hybrid-temporal-convergence-3s-2026-10-14.json`.

## ExJet source-term correction (2026-10-12)

The source PDF is now retained at
`reference/obstacle_wake/exjet_obstructed_free_jet_id252.pdf` (SHA-256
`2522d10547225dc09ae14b5748125118bba0e2e5ab1a593c5c172d39c1c0b1f8`). A
source audit found and corrected a consequential descriptor error: 12 mm is
the nozzle diameter for the *unobstructed* Daubech et al. case discussed on
the preceding page, not a reported direct-orifice diameter for the obstructed
case. The obstructed case's direct orifice is therefore recorded as unknown.
The four diameters/velocities/temperatures/compositions in Table 1 are retained
as model-specific pseudo-source alternatives (FLACS, CFX, OpenFOAM 1812/FRED,
OpenFOAM 1912+/EXORIS); they are CFD inputs, not experimental measurements and
must not be represented as four source truths.

Figures 11 and 12 have now been manually digitized into
`reference/obstacle_wake/exjet_obstructed_h2_plot_digitisation_v1.csv`. This is
14 uncertainty-tagged points: 3 centreline measurements and 11 radial points
at 1.4, 2.5, and 4.4 m. Plot-reading uncertainty is recorded separately from
physical sensor uncertainty (±0.010 centreline, ±0.015 radial H2 mole fraction;
±0.03 m radial position). It is a raster graph transcription, not author data.
Figure 11's digitized centreline values are 0.463, 0.289, and 0.207 at the
three stations; Figure 12 provides the corresponding partial radial profiles.
The paper supplies no time-resolved observation clock, so these profiles may
support spatial/steady-profile comparison only and cannot validate DEGALI's
transient response. The geometry description is still insufficient to claim
exact cylinder-array registration. Quantitative promotion remains withheld
until precise geometry, direct source characterization or defensible source
interval, and observation-state definition are obtained.

The BOS descriptor and fail-closed audit are
`reference/obstacle_wake/bos_obstacle_concentration_benchmark.json` and
`outputs/bos-obstacle-benchmark-audit-2026-10-12.json`. The audit records the
reported shape signatures and 3.81% inversion error while keeping
`promotion_allowed=false` because the paper does not supply a machine-registered
concentration grid or common-clock facility manifest.

The fail-closed ExJet readiness record is
`outputs/exjet-benchmark-audit-2026-10-12.json`; it currently reports
`withheld_missing_machine_data` and lists the missing geometry, observation,
clock and source-term artifacts. This output is intentionally generated even
before the ExJet raw package arrives so a future import cannot be mistaken for a
validated score. The separate H2SAFE readiness record is
`outputs/h2safe-evidence-audit-2026-10-12-v1/h2safe-audit.json`.

The DataverseNO component audit is
`outputs/dataverse-no-wskbij-audit-2026-10-12.json`; its small metadata package
is retained under `reference/obstacle_wake/dataverse_no_wskbij/`. It is marked
`component_source_pressure_ready` while keeping `promotion_allowed=false`.

The ExJet paper also gives a concrete modelling lesson for DEGALI: pseudo-source
selection materially changes the concentration result, and resolved versus
porous obstacle representation changes local velocity drops. That comparison is
now a required source-term/geometry sensitivity case, not an optional figure.

Figure 1 also gives a nozzle height of 0.75 m and a schematic array envelope
beginning 1.1 m downstream, spanning 3.5 m streamwise and 1.5 m vertically.
These reported dimensions are captured in
`reference/obstacle_wake/exjet_obstructed_schematic_geometry.json`, explicitly
tagged `schematic_envelope_only`. Cylinder count, individual row coordinates,
spanwise registration, and sensor-plane alignment remain unspecified. The
readiness audit now refuses to treat this envelope-only schematic as an exact
machine-registered geometry package.

The ExJet gate now distinguishes steady spatial-profile comparison from
transient validation. A common release clock is not a prerequisite for a
steady spatial profile score, but time-stamped measurements and a common event
clock are both mandatory for transient validation. The current digitized
package is not spatial-ready either: exact cylinder registration and measured
or independently bounded source terms remain absent. This separation avoids
both rejecting valid steady comparisons for lack of a clock and accidentally
promoting a steady graph as a transient benchmark.

The 3-D geometry layer now includes `CrosswindCylinder`, a finite horizontal
rod with circular x-z section and y-axis span, matching the qualitative ExJet
description of cylinders perpendicular to the streamwise jet. The paper reports
1.5 m rod length, 0.1 m diameter and approximately 0.35 m vertical separation;
its schematic also gives an array envelope. It does not supply the registered
rod count/centres or the direct obstructed-orifice diameter. The new primitive
therefore enables a faithful shape representation once those placements are
bounded, but does not make the current ExJet package machine-registered or
remove its source-term uncertainty gate.

`degali.validation.score_exjet_profile_points` now provides the spatial scoring
operator for a future paired DEGALI field export: exact observation-ID matching,
signed bias/MAE/RMSE in vol%-points, comparison against plot-reading uncertainty,
and a 4 vol% pointwise classification that marks observations indeterminate
when their digitization interval crosses the threshold. It never interpolates
or turns sparse points into a hazard volume, and it always labels these scores
as non-transient diagnostics. When paired with field-sampler bounds, the scorer
reports concentration-reading and coordinate-induced prediction ranges
separately and tests whether the combined residual interval contains zero.
Predicted 4 vol% classification is marked indeterminate when the coordinate
sensitivity range straddles the threshold. If callers omit those bounds, the
output explicitly indicates point-estimate-only classification.

`degali.validation.sample_exjet_spatial_profile_predictions` pairs a resolved
3-D final H2 mole-fraction field with these stations using an explicit vertical
measurement-plane transform. It uses trilinear cell-centred sampling only;
out-of-domain points and stencils touching solid cells are recorded as failed
rather than silently moved to a neighbouring fluid cell. Both observation and
source PDF hashes are checked by the readiness audit. The sampler uses only
the final field snapshot; it propagates radial-coordinate uncertainty with
endpoint and intervening cell-centre samples, and reports incomplete
propagation if any such sample is outside the domain or touches a solid
stencil. This bounded sensitivity is not a probabilistic confidence interval.
