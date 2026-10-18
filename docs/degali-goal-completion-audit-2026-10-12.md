# DEGALI reconstruction completion audit

This is a requirement-by-requirement audit of the active reconstruction goal.
`pass` means the current repository and output artifacts provide direct
evidence. `conditional` means the mechanism exists but the evidence is not
strong enough for an operational claim. `withheld` means the requirement is
not yet met.

| Requirement | Evidence | Status |
|---|---|---|
| Fail closed on legacy/non-validated obstacle concentration fields | Fast Tier-1 uses local scalar fluxes (`nonlocal_scalar_routing=false`) but explicitly reports `quantitative_prediction_allowed=false`; the resolved 3-D prototype also reports false; field decision rejects legacy non-local semi-FV routing | pass for the use-prohibition gate; validated obstacle prediction remains conditional/withheld |
| Resolve front, roof, side and lee-wake geometry | 3-D resolved solver, solid mask, front/roof/side/wake sampling in TI021 convergence audit | conditional: computed-window only |
| Conserve mass, momentum, H2 species and sensible enthalpy | 20-s prefix: species residual `4.1465e-16`, enthalpy residual `0`, continuity residual `4.6221e-11`, zero failed correctors | pass for computed window |
| Include density and buoyancy transition | active binary-mixture closure, temperature recovery and buoyancy feedback are implemented; a separate 80 K development witness reaches 21.14 vol%, 274.49 K minimum and 5.21% density departure with mass/energy residuals below 1e-15 | conditional: witness is a development case with unqualified boundaries/closure; TI021 facility screening itself has negligible density span (`3.04×10⁻¹²`) |
| Include turbulent mixing | molecular-plus-subgrid viscosity/diffusivity paths and diagnostics exist | conditional: no independent obstacle turbulence calibration |
| Provide screening and resolved tiers | fast local obstacle transport and resolved finite-volume/LBM paths are separate modules | pass at software-architecture level |
| Validate against AIJ Case H / neutral wake data | Hash-pinned aligned, time-averaged 3-D mean-flow audit over 176 U/W points is `accepted_component_grid`; no-fit ethylene-tracer scores cover 582 common concentration points | conditional mean-flow component; a matched 8–12 s D3Q19→generic ethylene-air FV tier comparison (same 25 mm mesh, 1 ms step, source and SGS scalar settings) lowers RMSE `2283.0→1007.8 ppmv` (55.9%) and bias `688.9→263.6 ppmv`. A second paired FV run isolates the configured nozzle-momentum term: RMSE `1008.09→1007.76 ppmv` (0.033% change), so it does not explain the tier improvement in this setup. FV continuity/global-mass residuals are `1.55e-13 kg m-3 s-1` / `2.05e-16 kg s-1`; mean concentration remains 3.16× measured and runtime is 270.6 s (18.1× D3Q19). Model average 8–12 s does not reproduce the 120-s measurement record; the mean-flow/scalar closure and source regularization still need work. No validation or H₂ promotion |
| Validate against SMEDIS/Thorney Island | TI021 600-s, 2.5-m score detects 42/53 positive sensors but has normalized dose bias `-0.888`, dose MAPE `0.865`, six-arc peak MAPE `0.519`; completed zero-surface 20-s, 1.25-m score detects 23/53 and has dose bias `-0.955`, dose MAPE `0.978`, six-arc peak MAPE `0.621` | conditional: full-horizon coarse grid aliases low sensors and underpredicts peaks/dose; fine run remains horizon-limited; independent hazard convergence remains withheld |
| Validate against public H2 obstacle/facility measurements | PRESLHY E3.4 manifests and H2 tunnel replay exist; H2SAFE archive is SHA-256 verified with STEP/STL geometry and five helium events; BOS obstacle-shape study is registered as a literature aggregate | withheld for obstacle-H2 qualification: geometry/flux/calibration/raw-coordinate gaps; H2SAFE remains helium-surrogate only |
| Provide front/top/side/wake vol% fields | 600-s screening package renders plan/top, front, side and an explicitly masked lee-wake shell; region metrics are hash-pinned in JSON | pass for computed screening window; finest-grid promotion withheld |
| Provide 4 vol% hazard region | Matched 20-s, 1.25-m wall comparison: voxel-wall endpoint/envelope `654.30/14730.47 m3`; zero-thickness endpoint/envelope `167.97/15431.64 m3` | conditional: computed-window geometry sensitivity is material; neither hazard field is independently validated or converged |
| Provide sensor time series and error metrics | Corrected full-event nested 600-s score covers all 145 sensors; 84 have positive dose, 54 predicted detections, dose bias `-0.899`, dose MAPE `0.850`, arrival MAE `31.76 s`; 24 shared interior sensors at six locations | conditional: strong underprediction, 29 missed available arrivals, 19.62% interior dose-coupling discrepancy, unresolved near-ground heights |
| Demonstrate grid and time convergence | matched-dt 0--20 s v2 audit holds the numerical zero-thickness surface fixed and contracts field L1/RMSE/L∞ norms; 4-vol% volumes are 14,875/15,000/15,431.64 m³; peak/dose monotonic convergence is 0.24/0.36 | withheld: hazard GCI and receptor convergence fail for the idealized surface; the physical tarpaulin-on-scaffold frame is incompletely registered and omitted |
| Reproducible evidence package | checkpoint hashes, manifest SHA-256, common-clock and artifact records exist | conditional: sensor calibration certificate missing |
| Report execution time | Corrected TI021 600-s inner/outer runs: `414.4 s` / `1573.1 s` wall clock; earlier 20-s grid runtimes retained above | pass for these reproducible calculations |

## Corrected TI021 nested overlap and full-event rerun (2026-10-13)

The earlier 600-s nested score used nine shared sensors located exactly on its
coupling plane. Those sensors cannot establish interior-domain coupling; that
score remains a useful historical 145-sensor baseline, but its overlap gate is
not a valid in-domain test. A corrected 60-s run moves the handoff to site
downwind `y=100 m` and expands only the inner crosswind span to 300 m, so the
full sensor union and a shared fluid region are both represented. It covers
all 145 sensors, including 24 true interior overlap sensors at six horizontal
locations; all matched coordinates agree. Maximum 0.6-s peak change is 4.86%
and integrated-dose change is 5.61% using the declared 0.1 vol% floors. The
nine sensors on the interface are reported separately, not used to claim
interior coupling.

The corrected 60-s model-versus-experiment score is deliberately still
withheld: it captures only the first 60 s of an event whose sensor records
extend later, and the dose score remains poor (84 positive observations,
normalized bias `-0.918`, MAPE `0.915`; 48 sensors predicted above 0.1 vol%).
These are not full-event accuracy estimates. A checkpointed 600-s rerun with
the corrected `300×200×40 m`, 5-m, `dt=0.125 s` inner setup and the same
`y=100 m` handoff has now completed. All internal gates passed: relative
species-mass residual `4.15e-16`, relative sensible-enthalpy residual `0`,
maximum continuity residual `2.45e-12 kg m-3 s-1`, zero nonconverged corrector
steps, and `414.4 s` wall time. The matched outer calculation also completed
600 s with all solver, conservative-remap and handoff gates passing (species
residual `7.11e-16`, enthalpy residual `0`, continuity residual
`1.11e-12 kg m-3 s-1`, remap residual `5.63e-16`, wall time `1573.1 s`).
The corrected run covers all 145 sensors; the shared interior overlap now has
24 sensors over six locations. Peak coupling passes at 4.65%, but integrated-
dose coupling fails at 19.62%. Against TI021, 84 sensors have positive dose;
the model detects 54, has dose MAPE `0.850`, normalized bias `-0.899`, and
arrival MAE `31.76 s` with 29 missed arrivals. Thus the full-horizon score is
now available but quantitative validation remains withheld. For example,
sensor S058's inner/outer dose is 7.402/7.725 vol%·s in 0--60 s, 19.057/15.567
in 60--120 s, and 2.914/0.343 in 120--200 s, localizing much of the coupling
mismatch to the later plume passage; this split is diagnostic, not a proven
physical cause. The prior v2 full-event run remains preserved for comparison.

## Current conclusion

DEGALI is a reproducible, conservation-checked obstacle-resolved research and
conditional screening system. It is not yet a fully validated field predictor:
the remaining decisive evidence is an independently calibrated, machine-registered
obstacle experiment with a time-synchronized H2 sensor network, plus hazard and
receptor convergence that passes the predeclared gates.

Newly imported ELVHYS WP4.2 tests 29/31 add a real-LH2, high-rate sensor and
restricted-vent component package with exact cross-file time vectors and
conditional event manifests. They do not close the independent open-field
obstacle-wake gate: the selected tests are confined vent-restriction cases,
their hydrogen mass-flow history and registered plate geometry are absent, and
event-specific calibration certificates were not supplied.

The hash-pinned evidence index is
`outputs/degali-validation-evidence-index-2026-10-12.json`; it currently reports
all listed artifacts present while keeping `promotion_allowed=false`.

The new FFI hybrid-patch witness closes an implementation gap between the fast
screening and resolved obstacle paths.  It uses fast transport upstream and
downstream of a resolved obstacle patch, synchronized conserved mass/species,
momentum and (in the production-shaped run) sensible-enthalpy packets.  The
60-s active thermal run completed in 56.26 s for 3,000 bridge exchanges in each
direction, with a maximum Mach number of 0.125, species-mass residual
`1.56×10^-17 kg`, enthalpy residual `0 J`, continuity residual
`1.97×10^-15 kg m^-3 s^-1`, and maximum pressure-correction residual `0.147 Pa`.
Its 24 shared FFI sensors have MAE/RMSE `4.081/5.892 vol%-points`; this is an
interface/conservation witness, not an independent facility-validation score.
The characteristic momentum mapping and enthalpy callback are covered by
regression tests.  A separate explicitly isothermal 3-s time-step audit gives
`dt=0.02→0.01 s` sensor mean absolute change `2.40×10^-4 vol%-points`, field
RMSE `0.00357 vol%`, and runtime `0.584→1.124 s`.  Because the active thermal
solver enforces its D3Q19 physical step, this audit is numerical/interface
convergence evidence only and does not promote the FFI geometry or close the
independent H2 obstacle gate.

For execution-time sensitivity, the same 60-s active-thermal case was rerun
with resolved-patch flow iterations reduced from 40 to 20 per scalar step.  The
wall time changed from `56.26` to `53.66 s` (4.62% reduction), while the 27
common predicted sensor means changed by only `3.83×10^-7 vol%-points` RMS
(`4.38×10^-6` maximum) and the assembled-field RMS change was
`2.25×10^-8 vol%`.  Species-mass, enthalpy and continuity residuals remained
at the same numerical scale.  This is an execution/solver-iteration
sensitivity result, not evidence that the conditional FFI geometry is
validated.

The FFI evidence boundary is now machine-readable at intake.  The wide USN
CSV parser accepts the unlabeled row-index column, records complete row counts
and SHA-256 values, and flags the separate source/sensor histories without
silently treating them as co-located receptors.  The O1/O2 intake records
(`outputs/ffi-o1-field-evidence-intake-2026-10-14-v2.json` and
`outputs/ffi-o2-field-evidence-intake-2026-10-14.json`) mark source and
receptor channels conditional, while weather, machine-registered beam
geometry and common clock remain withheld.  This is an explicit
`FieldEvidenceManifest` intake boundary, not a promotion or a substitute for
the missing event-linked files.

The AIJ neutral-wake evidence now separately indexes the newer aligned-grid
time-averaged 3-D audit at
`outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_time_averaged_aligned_grid_audit.json`.
This component-grid gate passes (176 common points, all three mean-flow runs
statistically stationary, obstacle faces aligned); the earlier instantaneous
3-D grid audit remains withheld. This does not change concentration validation
or hydrogen-prediction promotion, both of which remain false/withheld.

An observation-blind passive-ethylene component comparison now uses the AIJ
Case H 4 mm floor source and reported 5.83×10⁻⁶ m³/s source rate. The matched
8–12 s score samples 582 common concentration points: mean measured
concentration is 121.9 ppmv, versus 867.5 ppmv at 25 mm and 1,028.8 ppmv at
12.5 mm. RMSE is 2,624.5 and 5,449.1 ppmv, respectively. Both carrier-flow
means pass the 1% stationarity check and species-mass residuals are below
2.2×10⁻¹³, so conservation does not cure this scalar-model error. This is not a
pure grid-convergence test: timestep and Gaussian source width change with
resolution. The current source uses zero vertical exit momentum and omits a
validated turbulent-scalar covariance closure; therefore this is an adverse,
non-promoting component diagnosis.

The machine-readable counterpart is
`outputs/degali-goal-completion-audit-2026-10-12.json`. It currently reports
**2 pass, 6 conditional and 3 withheld** requirements; the withheld gates are
embedded-wall shear validation, full-event grid/time convergence, and
independent public obstacle-H₂ validation.

An additional observation-blind AIJ boundary audit expands the lateral domain
from `0.30` to `0.50 m` at the same `25 mm` spacing and the same `0.75 s`
carrier-flow horizon, translating the measurement origin so the obstacle and
source remain registered.  The wider domain lowers center-plane wake-velocity
RMSE from `0.06643` to `0.05523` (about 17%) and raises W correlation from
`0.881` to `0.893`, while U RMSE changes from `0.15665` to `0.17009`; both
cases are explicitly marked non-converged at this diagnostic horizon.  This
identifies side-boundary sensitivity as a real model uncertainty, not a reason
to fit concentration data.  It is retained as component evidence only in
`outputs/obstacle-physics-rebuild-2026-10-10/aij_case_h_lateral_domain_sensitivity.json`.

The same boundary change was then carried through the 8--12 s ethylene-air FV
candidate (SGS and inferred nozzle momentum unchanged) and scored on the same
582 point IDs.  Paired RMSE falls from `1007.8` to `763.5 ppmv` (`24.2%`), and
bias falls from `263.6` to `176.7 ppmv`; however, Pearson correlation falls
from `0.743` to `0.293`.  Thus the wider domain reduces one gross amplitude
error but does not yet provide a validated concentration field.  The paired
comparison is preserved at
`outputs/aij-case-h-lateral-domain-tracer-comparison-2026-10-14-v1/` and
remains explicitly non-promoting.

The requirement audit now also records the latest completed 600-s nested
145-sensor score (`outputs/thorney-ti021-nested-145-sensor-score-v2-600s-2026-10-13-v3/`).
It covers 84 positive-dose observations but predicts threshold arrival at only
54 sensors; normalized dose bias is -0.900 and dose MAPE is 85.4%. The nine
independent inner/outer overlap sensors fail the dose-coupling gate (maximum
relative dose change 13.0%, against the 10% limit). A geometric re-audit found
all nine nominal overlap sensors lie exactly on the y=100 m coupling plane;
there are zero shared sensor locations more than one cell inside both domains.
Thus this sensor comparison cannot validate an interior overlap solution,
even apart from its dose mismatch. In addition, 91 receptors
are at or below the inner-grid bottom-cell centre (2.5 m), so their vertical
positions are not resolved by the 5 m inner cells. This is a more complete
full-horizon diagnostic than the earlier 60-sensor summary, but it remains
explicitly withheld from quantitative promotion. The scorer now fails closed
unless at least three distinct horizontal sensor locations lie inside both
domains beyond the coupling-plane buffer; the regenerated score and generator
regression test retain this geometry diagnostic.

A short observation-blind spatial-overlap experiment then moved the handoff to
y=125 m and extended the resolved inner domain to 300×200×40 m. The completed
60-s inner and outer runs share 24 receptor locations at six distinct
horizontal positions, all more than one 5 m cell beyond the handoff plane.
Across this diagnostic window the maximum 0.6-s peak change is 3.83% and the
maximum integrated-dose change is 4.22%, both below the 10% component limit.
The inner solver took 104.8 s and the outer solver 215.3 s. This is the first
actual shared-interior sensor comparison rather than a comparison solely at
the inlet plane. It covers only 137/145 unique receptors and only 60 s, so it
does not replace the full-event 145-sensor score, establish grid convergence,
or validate DEGALI against the experiment. Hash-pinned details and the plot are
in `outputs/thorney-ti021-interior-overlap-60s-dx5m-handoff125-2026-10-13-v1/audit/`.

The matched-dt spatial audit supersedes the previous receptor-GCI calculation;
the earlier audit mistakenly used time-step ratios in a spatial GCI formula.
The corrected machine-readable grid/time item records the SHA-256-backed 20-s
temporal sensitivity report at
`outputs/thorney-ti021-temporal-prefix-comparison-20s-2026-10-12-v2/`; its
positive time-step result does not override the still-failed spatial hazard
GCI or receptor convergence gates.
The added matched fixed-surface 0--3 s diagnostic is at
`outputs/thorney-ti021-zero-surface-matched-3s-convergence-2026-10-12-v2/`.
It shows contracting field norms, but a 52.0% hazard-volume GCI, 33.3% peak
monotonic convergence and 0% dose monotonic convergence; its 20 s gate was
explicitly false for that short-prefix audit. The 1.25 m zero-surface run has
since completed and been independently materialized over 20 s; the matched
fixed-numerical-geometry multi-grid audit is also complete (v2), but hazard and
receptor gates fail and the real scaffold frame remains unrepresented.
The matched-dt 3-grid data and geometry metrics are at
`outputs/thorney-ti021-matched-dt-grid-convergence-2026-10-12-v4/`. It measures
the one-cell surface band's sampled fence arc lengths as 140/150/167.5 m versus
157.08 m analytically (10.87%/4.51%/6.63% error); this is non-monotonic and a
plausible contributor, not a proven sole cause, of the failed spatial gates.
The separate 1/0.5/0.25 m geometry-only audit passes but cannot qualify these
transport grids. The v4 gate explicitly rejects fixed-geometry spatial
promotion because represented wall thickness follows `dx`; field, hazard, and
receptor comparisons are joint grid/regularisation sensitivity.

## Follow-up evidence check (2026-10-13)

The current local regression check for the obstacle-resolved conserved-flow
path completed with **91 passed** across
`test_variable_density_momentum_3d.py`,
`test_variable_density_projection_3d.py`,
`test_binary_mixture_corrector_3d.py`, and
`test_resolved_obstacle_flow_3d.py` (61.67 s). This establishes tested building
block behavior only; it does not upgrade experimental validation or any
full-event gate.

The complementary reduced-order/field-gate regression check completed with
**80 passed** across `test_fast_obstacle_transport.py`,
`test_fast_obstacle_transport_3d.py`, `test_field_decision.py`, and
`test_field_workflow.py` (55.79 s). The 3-D potential-flow result is asserted
non-quantitative, and a legacy non-local obstacle-routing result is asserted
withheld even when conditional use is requested. This confirms separation of
the fast diagnostic tier from quantitative obstacle prediction at the tested
API/gate level; it does not qualify the resolved tier as an industrial
predictor.

A direct-H₂ unignited cylinder-impingement paper was added to the evidence
triage in
`docs/direct-h2-cylinder-impingement-benchmark-audit-2026-10-13.md`. The
publisher's accessible article record supports a 1.5 mm nozzle, 0.15--0.30 MPa
release pressure, and 0.20--0.40 m nozzle-to-cylinder spacing, but the article
page returned 403 and no machine-readable raw concentration/time-series package
was located in the scoped search. Exact geometry, registered sensor operator,
source history, and synchronized observations remain missing; this is
conditional evidence only, not a DEGALI score. The direct-H₂ obstacle
validation gate therefore remains withheld, and the overall reconstruction
goal remains incomplete.

## AIJ source-regularization sensitivity (2026-10-13)

The AIJ Case H scorer now accepts an explicit `--source-sigma-m` independent of
grid spacing, and the matched-comparison utility rejects any pair where that is
not the sole changed configuration field. At fixed 25 mm grid, 1 ms step,
8--12 s averaging interval, source rate, inferred source momentum, and SGS
settings, the 1 mm Gaussian sigma (D/4 radial-second-moment equivalent for a
4 mm circular opening) was compared against the prior 25 mm sigma case. Both
sample 582 identical observations and conserve source species to relative
residuals below `2.3e-13`.

The narrower under-resolved source increased mean prediction from 385.5 to
469.3 ppmv, raised RMSE from 1,007.8 to 1,930.3 ppmv, and increased the maximum
from 9,499 to 22,348 ppmv. It decreased prediction at 381 points and increased
it at 201 points, so the response is strongly spatially redistributed rather
than a uniform scale change. This is source-regularization sensitivity, not a
resolved nozzle model, observation-fitted result, or validation improvement.
The reproducible pair score and figure are
`outputs/aij-case-h-source-regularization-sensitivity-2026-10-13-v2/`; the
runner changes and tests are `tools/score_aij_case_h_resolved_tracer.py`,
`tools/compare_aij_case_h_source_regularization.py`, and
`tests/test_compare_aij_case_h_source_regularization.py` (6 targeted tests
passed including the scorer tests). The concentration residual remains large;
this sensitivity isolates source-kernel dependence but does not identify the
remaining transport-closure error.

In response, the resolved 3-D solver now reports Gaussian source FWHM in cells
by axis and emits an explicit warning when any axis has fewer than two cells
across the FWHM. This threshold is a diagnostic resolution flag, not a
validation or acceptance criterion. The 1 mm sensitivity kernel is 0.0942 cells
FWHM at 25 mm spacing and is flagged under-resolved; the 25 mm baseline kernel
is 2.355 cells FWHM. The matched report v2 records this distinction. The
resolved-flow, AIJ scoring, and pair-comparison regression tests pass (64
tests), including rejection of comparisons with any second changed setting.

### TI021 inner crosswind-span controlled prefix (13 October)

A valid width-640 m rerun was made after discovering that the earlier 200 s
wide-domain pilot retained the width-200 m receptor origin and handoff edges.
That pilot's score is explicitly invalidated for physical interpretation; its
corrected-coordinate rescore flags 57 inner sampler mismatches and a misplaced
crosswind origin. The corrected 120 s run centers the inner source at solver
`y=320 m` and maps the handoff plane over `q=-320..+320 m`; all 145 receptors
are covered, with 54 interior overlap sensors at 15 locations and zero inner
or outer sampler-coordinate mismatches.

To isolate the crosswind-span effect, the valid 200 m and corrected 640 m
runs were compared over the same 0--120 s window and the 36 common interior
labels. In the 12 sensors whose peak reaches 0.1 vol% in either member of
each pair, maximum peak coupling change decreases from 4.87% to 2.30%, and
maximum dose change decreases from 12.97% to 7.20% (dose p95: 12.68% to
7.20%). This is evidence that the previous inner lateral truncation was
contributing to active-plume coupling error. It is a short-horizon diagnostic,
not a full-event improvement claim: the predeclared all-sensor metric still
fails badly for near-zero signals (on the 640 m run: peak maximum 17.40%, dose
maximum 116.33% using the existing 0.1 vol% / 0.1 vol% s floor). The full
comparison is reproducible with `tools/compare_thorney_ti021_crosswind_span_prefixes.py`
and hash-pinned at `outputs/thorney-ti021-crosswind-span-prefix-comparison-2026-10-13-v1/`.
The official
gate and low-signal sensors are not redefined or excluded. The 600 s coupled
run still fails the dose gate (19.62%) and remains the authoritative full-event
coupling result. Artifacts: `outputs/thorney-ti021-wide-domain-overlap-score-120s-centered-2026-10-13-v1/`,
`outputs/thorney-ti021-long-domain-overlap-score-200s-2026-10-13-v1/`, and
the hash-pinned evidence index. Do not use the invalidated 200 s width-640
pilot as physical evidence.
