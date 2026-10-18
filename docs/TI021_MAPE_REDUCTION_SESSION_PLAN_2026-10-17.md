# TI021 released-mixture dose MAPE reduction session

## Latest status: species-volume observer repair, 28.592% diagnostic

The latest complete 600 s candidate uses both wall-aware diffusion and a
species-volume-conservative vertical observation operator. Its dose MAPE is
**28.5921439512%**, independently reproduced bit-for-bit in **25.1612 s**.
It improves the wall-repaired 29.0965% total but worsens the 2.4 m and 6.4 m
height groups and two distance groups. **Strict promotion is rejected and the
10% goal is not achieved.** The 29.0965% control also contains the now-known
legacy vertical observation defect; it is a historical comparison, not a
qualified model. See the final section for the complete iteration evidence.

The goal remains active. The preceding review was progress because it
reproduced an observation-operator defect and changed the next action from
parameter tuning to reconstruction repair. The current iteration implements
that repair and measures its full-event accuracy, resolution and runtime.

## Prior status: wall-diffusion defect found and repaired

This section supersedes the historical promotion labels below. The 28.0407%
result is reproducible but has a newly identified physical boundary defect:
horizontal turbulent diffusion ignored finite-height wall faces, although
advection respected them. Global conservation alone did not detect it.
The repaired diagnostic gives **29.0965% MAPE in 22.2871 s** for the complete
600 s event. It is not an accuracy improvement or a qualified prediction.
The 10% target remains unmet. The historical 28.0407% model now fails the
known-wall-defect eligibility gate, and the repaired model fails the accuracy
and group-nondegradation gates. Neither is a newly accepted benchmark.

All MAPE figures here compare time-integrated released-mixture concentration
(vol% s), not instantaneous hydrogen concentration. TI021 is an R12/N2
surrogate experiment. Repeated selection on its observations is calibration,
not independent external validation, even when the solver does not read dose
observations while advancing its state.

The fixed-prediction optimizations below are diagnostic bounds only within
their specified pair/tower inventories or finite convex hull. They are not
proof that a revised physical model, vertical transport, or untested
parameter choice cannot attain 10%.

## Objective and acceptance gate

Reduce the TI021 0--600 s dose MAPE over all 84 positive observations to at
most 10%. The production metric is the time-integrated **released R12/N2
mixture concentration** in vol% s. A result is accepted only when it:

- recomputes all 600 s and all 145 receptors without sensor-specific factors;
- passes species-mass, mixture-mass, and momentum conservation checks;
- does not materially worsen a height or downwind group;
- is reproducible from a saved configuration and result file; and
- completes in at most 30 s on the current machine.

The loop is: diagnose residual structure, change one physical mechanism, run
focused tests, recompute the full event, and adopt or reject from the gates.

## Metric correction found on 2026-10-17

The scanned field report states that GAS and HGAS readings are concentration
of the **released gas mixture**. The former scorer instead compared its R12
component mole fraction with those readings. Consequently the historical
29.5441456% and the 29.223--29.481% variants are R12-component regression
diagnostics, not valid scores against the reported observable.

The source input was also audited. The report gives a 14 m-across, 13 m-high
initial tent. `RegularPolygonPrismInitialCloud` interprets the legacy
`source_initial_radius_m=14` input as a full polygon diameter and stores a 7 m
apothem, so the geometry itself was correct; only the parameter name was
misleading. The CLI now passes this input through and every output records the
diameter, apothem, and legacy-name semantics explicitly.

The scorer now supports two explicit modes:

- `released_mixture`: production TI021 comparison, pseudo-species molar mass
  0.0585 kg/mol and mass-weighted heat capacity;
- `r12_component`: backward regression only.

The run wrapper and residual auditor refuse to compare candidates whose
observable modes differ.

## Verified results

| Run | Dose MAPE | NMB | Solver runtime | Decision |
|---|---:|---:|---:|---|
| Correct-observable untuned run | 176.569% | +1.264 | 14.03 s | reject |
| Entrainment 0.30, horizontal diffusivity 2.5 m2/s | 41.993% | -0.040 | 16.21 s | diagnostic best without elevated branch |
| Same plus wall-impact elevated branch (0.5 transfer, 0.02 s-1 return, 6 m centre, 1.5 m sigma) | 40.869% | -0.050 | 25.40--27.66 s | overall improvement, but reject strict promotion because the 6.4 m group worsens 5.88 percentage points |
| Radially relaxing cavity vertical closure (exponent 12, radial power 1) | 39.591% | +0.015 | 15.52 s | all height and distance groups nondegraded |
| Same plus exterior entrainment relaxation 0.5/100 m and vertical-tail diffusivity 0.01 m2/s | **33.597%** | -0.030 | 16.40 s | current strict candidate; every height and distance group improves |
| Same with source crosswind scale 1.0 | 33.375% | -0.026 | 17.60 s | lower overall MAPE, but reject strict promotion because three groups regress slightly |
| Same plus conservative elevated-branch slip transport | 34.218% | -0.049 | 26.75 s | reject; 6.4 m group worsens 4.01 percentage points and total MAPE rises 0.62 points |
| Spatially local conserved wake reservoir, 0.005/0.005 s-1 | **31.893%** | -0.054 | 20.47 s | diagnostic best; reject strict promotion because the 6.4 m group worsens 0.27 percentage points |
| Local visible reservoir, 0.005/0.1 s-1 | 32.613% | -0.015 | 15.56 s | reject; 0.4 m, 6.4 m, and far groups regress |
| Species-only local delay, 0.005/0.1 s-1 | 32.762% | -0.017 | 18.07 s | reject; preserved carrier mass and momentum but did not widen the pulse enough |
| Wall-elevated branch with 18 m shear-layer development ramp | 33.341% | not promoted | 28.25 s | reject; total improves slightly but 2.4 m, 6.4 m, and middle-distance groups regress |
| Local reservoir 0.005/0.01 s-1 + spatially advected independent wind + unexpanded 14 m source footprint | **31.724%** | -0.046 | **14.43 s** | **promoted current strict candidate; every height and distance group improves, conservation passes** |
| Same plus longitudinal-only diffusivity 5 m2/s | 32.277% | not promoted | 16.88 s | reject; pulse broadening worsens dose allocation |
| Same plus moving local reservoir at 0.5 ambient-speed fraction | 31.674% | not promoted | 16.14 s | lower total, but reject because 2.4 m and 6.4 m groups regress |
| Same plus fully visible Gaussian reservoir at z=3 m, sigma=2 m | 48.111% | not promoted | 16.98 s | reject; stored carrier volume grossly overstates 4.4 m dose |
| Age-transported wall branch (3.5 m centre, 2 m initial sigma, settling 0.005 m/s, Kz 0.05 m2/s, 50 m development) | **31.554%** | -0.019 | **24.82 s** | diagnostic best, but reject strict promotion because the <150 m group worsens 0.124 percentage points |
| Transported conservative wake reservoir, 200 s residence / full local-wind speed | 32.798% | -0.099 | 16.18 s | reject; 105,607 kg remains hidden and dose is underpredicted |
| Transported conservative wake reservoir, 100 s residence / half local-wind speed | 31.936% | -0.067 | 18.64 s | reject; 2.4 m group worsens 2.508 percentage points |
| Coupled transported wake reservoir + age-resolved wall branch | 31.791% | not promoted | 29.98 s | reject; 2.4 m and near-distance groups regress |
| Global friction-velocity sweep, 0.16/0.18/0.23 m/s | 38.811/34.346/34.359% | not promoted | 23.59--24.43 s | reject; changing the whole plume speed cannot repair the path-specific residual |
| Crosswind diffusivity factor sweep, 0.5--3.0 | 32.057--37.848% | not promoted | 21--24 s | reject; one scalar width cannot fit centreline and outer-arc doses together |
| Conserved wall-elevated crosswind branch, +0.15 m/s, 150 m development | 30.981% | not promoted | 41.33 s | lower total and six of seven groups improve; reject because 6.4 m worsens 1.159 points and runtime exceeds 30 s |
| Same branch with conservative 1 s transport subcycle | **30.967%** | -0.018 | **28.95 s** | best diagnostic and under runtime gate; reject because 6.4 m worsens 1.185 points |
| Historical R12-component exact regression | 29.544% | -0.265 | 15.82 s | exact backward regression only |

The correct-observable sweep covered entrainment coefficients 0.18--0.50 and
horizontal diffusivities 1--4 m2/s. The 0.30/2.5 combination removed nearly
all global bias. A geometry-only vertical closure then concentrated the
collapsed source core near the ground and relaxed that concentration toward
the 50 m fence. Outside the fence, the wake-mixing coefficient relaxes over a
100 m length to half its near-wall value. Reducing the unresolved vertical
tail diffusivity from 0.28 to 0.01 m2/s removed the remaining 6.4 m regression.
None of these closures contains a sensor ID or observed dose.

Rejected checks include mean wind angles -12.2, -6.1, 0, and +12.2 degrees;
measured and direction-only wind histories; crosswind diffusion factors 2 and
4; time-relaxing global entrainment; bilinear receptors; and a four-layer
height-resolved solve. The four-layer run produced 73.716% MAPE in 40.37 s,
failing both accuracy and runtime gates.

All reported runs cover the full 600 s. Species relative residuals are below
5.9e-15 and the momentum gate passes. The focused regression suite has 50
passing tests. The default historical path remains bit-for-bit at 29.5441456%
in `r12_component` mode. The elevated branch now has a separately transported,
bounded species state whose flux is limited by donor inventory and receiver
mixture capacity. Its default is disabled. Its first full-event result was
rejected, so it is implementation evidence rather than a promoted model.

The new crosswind branch is a bounded subset of wall-impacted elevated
species, not a receptor multiplier. It moves perpendicular to the local wind,
carries its temporal first moment, and is limited by donor inventory and
receiver mixture capacity. A one-second conservative transport subcycle cut
its full-event runtime from 41.33 s to 28.95 s while changing MAPE by only
0.013 points. Its direction is supported by the residual: positive deflection
improves MAPE whereas the equal negative deflection worsens it to 31.747%.
It is nevertheless not promoted because a single elevated Gaussian still
forces the low/mid and 6.4 m paths to share one lateral trajectory. The next
mechanism must split the wall-elevated inventory into a lower separated-shear
branch and an upper over-wall branch with independent conservative paths.

## Residual evidence and next mechanism

The former 33.597% strict candidate was also compared with all workbook timing
summaries and eight independently digitized gas traces. Its mean peak-time
error is -32.9 s, its geometric-mean peak ratio is 1.335, and its median
equivalent-duration ratio is 0.558. It is therefore too early, too peaked,
and roughly half as wide as the measurements. A local wake reservoir improves
the integrated score to 31.893%, but a fully visible reservoir overfills the
4.4 m layer and a species-only delay leaves the timing-width error nearly
unchanged. These rejections isolate the missing mechanism as transported,
height-resolved wake age and meander rather than a scalar residence-time
constant.

The promoted 31.724% candidate combines three independent mechanisms: local
conservative cavity residence, spatial advection of the independently
digitized mast-direction history, and removal of the prior 1.5x crosswind
expansion of the registered 14 m source footprint. Relative to 33.597%, the
height-group improvements are 0.984, 1.381, 4.406, and 0.526 percentage
points, and the three distance-group improvements are 1.865, 1.150, and
2.670 points. Its latest-code rerun is bitwise stable in score and completes
in 14.43 s. The timing residual remains material: mean peak-time error
-31.6 s, peak ratio 1.418, and equivalent-duration ratio 0.553.

The next implemented diagnostic carries a conservative first temporal moment
with the wall-elevated species marker. Its local mean age lowers the vertical
centre by a declared settling velocity and grows the Gaussian variance by
`2*Kz*age`; capture, local storage, return, base advection, and relative
advection all carry the same age moment. The best full-event setting lowered
MAPE to 31.554% and improved all four height groups, but the <150 m distance
group increased from 30.341% to 30.465%, so it was not promoted. Its temporal
audit changed mean duration error only from -48.26 s to -46.20 s and median
equivalent-duration ratio only from 0.5525 to 0.5550. This falsifies vertical
age projection as the main timing remedy.

A second new state allows locally captured conserved material to advect out of
the cavity through a geometry-defined wake corridor before first-order return.
It passes conservation and moves stored mass outside the cavity in unit tests,
but both 100 s and 200 s residence candidates worsened the full-event score.
The 100 s run improved 0.4, 4.4, and 6.4 m groups while worsening 2.4 m by
2.508 points. Therefore a transported reservoir cannot return into the same
column-mean vertical state: its transported age and vertical moments must be
coupled before another residence sweep. A latest-code disabled-feature rerun
reproduced 31.724296939970355% exactly, proving no default-path regression.

Updated optimistic residual bounds for 31.724% are 27.326% after perfect
conservative left/right redistribution, 21.672% after perfect conservative
vertical redistribution within each tower, and 18.906% with one fitted
scalar per horizontal tower. Since all remain above 10%, the next admissible
step is not another single coefficient: it requires jointly transported
ground, over-wall, and wake-age states or a resolved local 3-D patch.

For the current strict 33.597% candidate, optimistic diagnostic bounds are:

- perfect conservative left/right redistribution: 28.739%;
- perfect conservative height redistribution within each tower: 22.349%;
- one fitted scalar for every exact downwind arc: 25.446%;
- one fitted scalar for every horizontal tower: 18.513%.

These are fitted lower-bound diagnostics for fixed predictions, not models.
The tested redistributions do not meet 10%; they do not exclude a physical
change that also changes those inventories. They motivate investigating
distance-dependent transport jointly with height allocation:

1. add independent but conservative ground, over-wall, and separated-wake
   transport states so vertical classes can have different horizontal paths;
2. drive their exchange by wall pressure, buoyancy, and wake-residence time,
   not by distance- or sensor-specific dose factors;
3. replace the rejected ambient-minus-layer slip law with an over-wall state
   whose momentum is initialized from blocked-face impulse and evolves with
   drag and settling before it exchanges with the ground/wake state;
4. validate timing and residence width against the digitized independent gas
   traces before the full dose run;
5. cross-check the same closure on an independent obstacle/wake case before
   any field-qualification claim.

## Reproducible artifacts

- `outputs/ti021-released-mixture-observable-600s-2026-10-17.json`
- `outputs/ti021-mixture-e0p3-d2p5-600s-2026-10-17.json`
- `outputs/ti021-mixture-e0p3-d2p5-wallbranch-f0p5-r0p02-c6-s1p5-paired-600s-2026-10-17.json`
- `outputs/ti021-mixture-vexp12-radial1-600s-2026-10-17.json`
- `outputs/ti021-mixture-strict-vexp12-extentr0p5-vtail0p01-600s-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-residual-audit-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-wallslip-f0p5-u1-r0p02-c6-s1p5-600s-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-temporal-residual-audit-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-localres0p005-600s-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-scalardelay-cap0p005-ret0p1-600s-2026-10-17.json`
- `outputs/ti021-mixture-scalardelay32p762-temporal-residual-audit-2026-10-17.json`
- `outputs/ti021-mixture-strict33p597-wallbranch-ramp18-f0p5-u1-r0p02-c6-s1p5-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-latest-code-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-temporal-residual-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-residual-structure-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-wakeage-f0p5-u0-r0p005-c3p5-s2-set0p005-kz0p05-ramp50-600s-2026-10-17.json`
- `outputs/ti021-mixture-wakeage31p554-temporal-residual-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-transported-wake-cap0p005-ret0p005-u1-l500-w50-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-transported-wake-cap0p005-ret0p01-u0p5-l500-w50-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-post-wake-age-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p724-post-crossbranch-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-wakeage-x0p15-ramp150-subcycle1s-600s-2026-10-17.json`
- `outputs/ti021-mixture-wakeage-x0p15-ramp150-600s-2026-10-17.json`
- `outputs/ti021-mixture-wakeage-crossm0p05-600s-2026-10-17.json`
- `outputs/ti021-frozen-r12-exact-regression-after-spatial-closures-2026-10-17.json`

The 10% target is not reached. No candidate is promoted as validated or
field-qualified.

## Local 3-D response coupling iteration (latest)

The local obstacle patch was corrected to score the registered released
R12/N2 mixture rather than the R12 component.  A paired obstacle/control run
was then completed with identical 5 m cells, 600 s horizon, source and wind.
Both runs satisfy the species ledger and the end-to-end 30 s runtime gate, but
their raw sensor fields are not suitable as a replacement model:

| local patch candidate | MAPE | runtime | decision |
|---|---:|---:|---|
| obstacle, integral-to-patch | 83.970% | 24.65 s | reject |
| no-obstacle control | 70.174% | 18.01 s | reject |

Even an optimistic single global amplitude rescaling leaves 66.42% and
59.98% MAPE, respectively.  The patch nevertheless provides conservative
obstacle-response moments.  At the downstream face the obstacle changes the
mean exit time from 14.06 to 24.41 s, the time standard deviation from 5.31
to 15.57 s, the height centroid from 4.15 to 11.21 m, and the fraction above
5 m from 0.242 to 0.892.  Splitting the obstacle packets shows a 5 m-and-below
wake centroid `(y,z)=(340.97,2.5) m` and an above-wall centroid
`(327.08,12.27) m` at essentially the same 24.4 s mean exit time.

The dense-layer solver now supports the following conservative, default-off
states derived from those moments:

- finite-duration elevated-branch advection scale;
- elevation marking of only the species flux that actually crosses a finite
  wall face;
- conservative local-reservoir crosswind convergence toward a wake axis;
- boundary packet space/time/height partitions in the local 3-D diagnostic.

The focused regression suite is **53 passed**. Full-event candidates were:

| coupling candidate | MAPE | runtime | decision |
|---|---:|---:|---|
| all impacted material, local-3D vertical/lag response | 33.013% | 28.76 s | reject |
| crossed-wall flux, vertical/lag response | 32.137% | 27.11 s | reject |
| crossed-wall flux, vertical response without lag | 32.146% | 24.78 s | reject |
| local-3D overflow residence mixture | 31.684% | 24.13 s | reject: 6.4 m and near group regress |
| advected/focused wake reservoir | 32.181% | 16.25 s | reject |
| cavity-only focused reservoir | 31.886% | 16.16 s | reject |
| height-split crosswind response | 32.919% | 32.85 s | reject: accuracy and runtime |

The best strict accepted score therefore remains **31.72429694%**.  The
physically derived overflow residence candidate lowers the global number by
0.040 percentage point but is not promoted because the 6.4 m and <150 m
groups regress.

An observable-safe model-span audit was also added.  It excludes legacy R12
component results and admits only the released-mixture observable.  Across
118 unique 600 s predictions, the best individual is 30.967%; even the
optimistic non-negative convex combination of complete model runs has a
29.205% MAPE lower bound within that finite convex hull. Blending those saved
predictions cannot reach 10%; this does not bound untested parameters or
physical models. It motivates testing height-dependent horizontal transport,
validated carrier response, or synchronized raw concentration history.
The frozen-carrier local patch is not promoted as that carrier
because its height-split crosswind response worsens the field score.

Latest reproducible artifacts:

- `outputs/ti021-integral-to-patch-released-mixture-flow1000-2026-10-17.json`
- `outputs/ti021-integral-to-patch-released-mixture-no-obstacle-flow1000-2026-10-17.json`
- `outputs/ti021-mixture-local3d-response-f0p858-z11p214-s5p629-a0p576-t24p41-r20-600s-2026-10-17.json`
- `outputs/ti021-mixture-wallflux-local3d-response-z11p214-s5p629-a0p576-t24p41-r20-600s-2026-10-17.json`
- `outputs/ti021-mixture-walloverflow-local3d-delay-f0p667-r0p0644-600s-2026-10-17.json`
- `outputs/ti021-mixture-local3d-cavity-focus-k0p0016-y1p24-600s-2026-10-17.json`
- `outputs/ti021-mixture-wallflux-local3d-vertical-crossm0p569-z12p267-s5p629-r20-600s-2026-10-17.json`
- `outputs/ti021-correct-observable-model-span-audit-2026-10-17.json`

## Height-dependent path and wind-propagation iteration

The fixed-height state previously carried a vertical concentration inventory,
but every layer inherited the same depth-mean horizontal path.  A conservative
relative-transport operator was added so fixed-height layers can now follow
independent streamwise and obstacle-wake crosswind velocities while exchanging
mass through the existing vertical diffusion operator.  All four conserved
components are moved together, finite-height walls admit only the portion of a
layer above the wall, and the closed relative-transport boundary preserves
global mixture, species, and momentum.

The work also exposed a pre-existing discrete defect: independent Rusanov
mixture/species splitting could leave an admissible column with a layer-local
species inventory greater than its mixture inventory.  A conservative
within-column water-fill now moves only the invalid excess to another height
layer instead of clipping it.  The focused regression suite is **55 passed**.

Paired 600 s results showed that the fully fixed-height observer is not yet a
replacement for the accepted continuous-column observer:

| candidate | MAPE | runtime | decision |
|---|---:|---:|---|
| two fixed layers, no relative path | 53.309% | 21.02 s | reject |
| two fixed layers, local-3D cross-height shear -0.058 s-1 | 54.744% | 23.67 s | reject |

The new path state is conservative and executable, but the fixed-layer volume
occupancy operator underpredicts the lower three sensor heights.  It remains
default-off and experimental rather than being tuned against sensor gains.

A separate residual diagnosis found that the spatial wind-history propagation
speed is consequential.  Reducing the externally advected mast-record speed
from 3.4 to 3.2 m/s lowered MAPE from 31.724% to 31.533%, but 6.4 m and the
150--300 m group regressed slightly.  Coupling that change with the existing
wall-impact elevated branch removes both regressions.  The final retained
candidate uses:

- 0.5 s maximum layer step;
- 3.2 m/s spatial wind-history propagation;
- unit wall-impact transfer intensity with 0.005 s-1 return;
- a 3.5 m centre, 2.0 m sigma, 0.005 m/s settling, and 0.05 m2/s vertical
  diffusion for the elevated observation branch;
- a 150 m development ramp and 0.025 m/s crosswind slip, updated every 1 s.

Its full-event result is **31.18804585114051% MAPE** in **14.63 s** on the
independent reproduction.  Every preregistered height and distance group
improves relative to 31.724%: deltas are -0.274, -0.725, -0.990, and -0.049
percentage point for 0.4/2.4/4.4/6.4 m, and -0.219, -0.719, and -1.464 points
for <150/150--300/>=300 m.  Species residual is below 1.3e-15 relative,
mixture residual is below 4e-8 kg, and the momentum gate passes.  A second
run reproduces all 145 predicted doses bit-for-bit (`max abs difference=0`).

This **31.188%** artifact replaces 31.724% as the current strict accepted
benchmark.  It does not satisfy the 10% session target, so the goal remains
active.  The mirrored wind-coordinate hypothesis (74.443%) and a -2 degree
history rotation (34.405%) are rejected.

Reproducible artifacts:

- `outputs/ti021-mixture-promoted31p188-post-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted31p188-post-regression-600s-2026-10-17.configuration.json`
- `outputs/ti021-mixture-multilayer2-control-vs31p724-600s-2026-10-17.json`
- `outputs/ti021-mixture-multilayer2-crossshear-m0p058-d70-vs31p724-600s-2026-10-17.json`
- `outputs/ti021-mixture-wind-advection3p2-vs31p724-600s-2026-10-17.json`
- `outputs/ti021-mixture-wind-history-mirrored-m12p2deg-vs31p724-600s-2026-10-17.json`
- `outputs/ti021-correct-observable-model-span-audit-after31p188-2026-10-17.json`

After adding the new runs, the observable-safe span contains 130 unique
released-mixture predictions.  Its optimistic whole-run convex MAPE remains
29.205%, so the newly promoted branch improves the strict engineering path but
does not remove the structural residual that prevents the 10% target.

## Instantaneous-release and residence-time recheck

The authoritative workbook was reread before introducing any additional
source history.  TI021 is explicitly a `puff`; its release-duration field is
1 s and its registered mass proxy is 4,962 kg.  The code intentionally retains
the workbook rate/duration fields as evidence rather than silently treating
them as a continuous atmospheric boundary.  A long finite source is therefore
not an admissible explanation for the measured pulse width.

Four conservative local-reservoir candidates were rerun from the promoted
31.188% state.  The best scalar residence setting (`capture=0.0025 s-1`,
`return=0.01 s-1`) reached 31.062%, but it regressed preregistered height and
distance groups and is not promoted.  More importantly, its temporal audit
worsened the actual structural residual:

| temporal diagnostic | promoted 31.188% | scalar residence 31.062% |
|---|---:|---:|
| mean arrival error | -22.43 s | -22.48 s |
| mean departure error | -72.91 s | -78.70 s |
| mean duration error | -51.50 s | -57.25 s |
| digitized mean peak-time error | -32.13 s | -32.38 s |
| median equivalent-duration ratio | 0.552 | 0.536 |

Negative errors mean the model is early.  The scalar reservoir therefore
slightly improves the integrated score while making the pulse still narrower;
it is the wrong mechanism for the observed time history.  Additional tuning of
one capture/return pair is stopped.

The next admissible iteration is a conservative two-path carrier model rather
than another dose multiplier:

1. split only physically wall-impacted inventory into a prompt over-wall path
   and a slower separated-wake path;
2. carry mixture mass, R12 mass, both horizontal momenta, and temporal moments
   on both paths, with exact exchange closure;
3. give the slow path an independently transported wake-age distribution and
   crosswind spreading, rather than returning it into the same column state;
4. reject a candidate first if the digitized trace peak time and equivalent
   duration do not improve, before spending a full promotion comparison;
5. promote only after the same 600 s dose MAPE, conservation, runtime, and all
   height/distance non-regression gates pass.

New diagnostic artifacts:

- `outputs/ti021-memory-cap0p0025-ret0p005-vs31p188-600s-2026-10-17.json`
- `outputs/ti021-memory-cap0p0025-ret0p01-vs31p188-600s-2026-10-17.json`
- `outputs/ti021-memory-cap0p005-ret0p005-vs31p188-600s-2026-10-17.json`
- `outputs/ti021-memory-cap0p01-ret0p005-vs31p188-600s-2026-10-17.json`
- `outputs/ti021-memory-cap0p0025-ret0p01-temporal-audit-2026-10-17.json`

## Wall-elevated path iteration and current benchmark

The separate slow wall-impact wake state was implemented as a default-off,
conservative four-component store.  It captures only blocked-face incident
flow, transports independently, and returns with a declared first-order rate.
It improved the mean departure error from -72.91 s to -55.27 s and the mean
duration error from -51.50 s to -33.95 s, but worsened the integrated dose
MAPE to 31.342%.  Gaussian observation of that state also worsened the score.
The timing improvement is real, but this path is not promoted because it puts
the delayed mass in the wrong horizontal/vertical sensor locations.

The accepted mechanism instead retains the existing wall-elevated inventory
and changes only its obstacle-induced horizontal path and normalized vertical
projection.  A bounded sweep found a clear crosswind-path optimum: MAPE falls
through 0.09, 0.15, and 0.30 m/s, then rises at 0.40--0.70 m/s.  The retained
setting is 0.30 m/s crosswind velocity, 1.5 m initial vertical centre, 0.9 m
initial sigma, 0.005 m/s settling, 0.05 m2/s vertical diffusivity, and a 150 m
development ramp.  The crosswind speed is 8.8% of the 3.4 m/s reference wind.

The independent 600 s reproduction gives **29.4918243491251% MAPE** in
**15.35 s**.  All 145 predicted doses reproduce bit-for-bit.  Relative to the
previous 30.5616508% benchmark, height-group deltas are -1.744, -1.582,
-0.294, and -0.004 percentage point for 0.4/2.4/4.4/6.4 m; distance-group
deltas are -0.038, -4.307, and -1.264 points for <150/150--300/>=300 m.
The maximum relative species residual is 7.4e-16, the mixture residual is
1.5e-8 kg, and the momentum gate passes.  This result is the new strict
benchmark, but it is not field-qualified and does not meet the 10% target.

Two falsification checks were added.  A conservative elevated-only crosswind
diffusion operator spreads the two obstacle shear layers without changing the
ground carrier, and the existing conserved tangential-wall bypass was exposed
to the run wrapper.  Elevated diffusivities 0.5--4 m2/s produced
29.532--29.873% MAPE; bypass fractions 0.05--0.40 produced 29.532--29.859%.
Both are worse than the zero-value 29.492% result and are rejected.  The new
diffusion path remains disabled by default.  The focused regression suite is
**54 passed**.

Residual structure remains decisive: at the new benchmark, the x<=-75 m
group still has 40.50% MAPE and -28.72% mean bias, while the two x>=75 m
positive sensors have 43.33% MAPE and -40.21% bias.  More one-sided drift
cannot repair both tails.  The next iteration must infer a geometry-based
left/right split of wall-impact inventory or a time-varying meander state,
then validate that state against timing and the independent obstacle evidence.

Current reproducible artifacts:

- `outputs/ti021-mixture-promoted30p562-post-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted29p492-post-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted29p492-post-regression-600s-2026-10-17.configuration.json`
- `outputs/ti021-dualpath-wc0p02-wr0p01-u0p3-temporal-audit-2026-10-17.json`
- `outputs/ti021-cross0p30-c1p5-s0p9-d0p5-vs30p562-600s-2026-10-17.json`
- `outputs/ti021-cross0p30-c1p5-s0p9-b0p05-vs30p562-600s-2026-10-17.json`

## Geometry-split shear layer and full-carrier wake iteration

A geometry-only lateral split was added to the wall-elevated path.  The local
crosswind coordinate is computed from the obstacle centre and instantaneous
carrier direction; material on either side is transported outward with a
smooth `tanh` transition.  The operator contains no receptor coordinates or
observations.  A symmetric unit test confirms that it preserves the lateral
centroid, increases the second moment, and retains exact species and momentum
closure.  An optional age-decay timescale represents downstream wake recovery.

The retained setting adds 0.095 m/s of outward motion with a 50 m transition
scale to the existing +0.30 m/s asymmetric branch.  Its independently
reproduced full-event result is **29.1906550801644% MAPE** in **14.34 s**.
All 145 doses reproduce bit-for-bit on the latest code.  Relative to 29.492%,
all four height groups and all three distance groups improve; the distance
deltas are -0.448, -0.078, and -0.021 percentage point.  Conservation and
momentum gates pass.  At 0.10 m/s the far group regresses, so stronger
permanent splitting is rejected.  Age-decaying 0.2--0.3 m/s candidates also
regress low-height or far groups and are not promoted.

The current temporal residual remains structural: mean arrival/departure/
duration errors are -22.48/-69.93/-48.49 s, the digitized peak is 32.25 s
early on average, and the median equivalent-duration ratio is 0.568.  The
model is still too early and too narrow even though its integrated spatial
score improved.

The updated optimistic bounds are 24.545% for perfect conservative mirrored
left/right redistribution, 19.279% for perfect vertical redistribution within
each tower, 23.070% for one scalar per downwind arc, and 17.305% for one scalar
per horizontal tower.  None reaches 10%; simultaneous horizontal, vertical,
and temporal transport remains necessary.

The meteorology coordinate audit confirmed that solver `x` is the physical
downwind coordinate and the wind history is advected along that axis.  The
stored alpha-0.1 history reduces the raw 5.64-degree heading standard
deviation to 0.56 degrees.  Reinstating 0.25/0.5/1.0 of the independently
digitized fluctuation worsened MAPE to 29.339/30.430/34.724%, so suppressed
meander is not the missing correction for this model.

Finally, the separate wall-impact wake reservoir was extended to carry all
four conserved components along independent streamwise, asymmetric
crosswind, and geometry-outward paths with an optional normalized vertical
observer.  Full-event candidates passed conservation and stayed below 30 s,
but produced 29.389--30.403% MAPE and regressed the 2.4/6.4 m or middle-
distance groups.  They are rejected and default-off.  The focused regression
suite is now **55 passed**.

New artifacts:

- `outputs/ti021-mixture-promoted29p191-post-fullwake-regression-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted29p191-temporal-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted29p191-residual-structure-audit-2026-10-17.json`
- `outputs/ti021-base0p30-out0p095-s50-vs29p492-600s-2026-10-17.json`
- `outputs/ti021-outward095-wind-alpha1p0-vs29p191-600s-2026-10-17.json`
- `outputs/ti021-fullwake-wc01-r01-u03-c02-o01-z1-s1-vs29p191-600s-2026-10-17.json`

## Species-weighted full-carrier wake promotion

The first full-carrier wake candidates captured about 58 t of total mixture to
move only 182 kg of released-mixture marker.  The branch was repeatedly
capturing highly diluted ambient carrier at the wall.  A new
`species_weighted` capture mode multiplies blocked-face capture rate by the
local released-mixture mass fraction while still transferring total mass,
species mass, and both momenta together.  The mode is geometry/field based;
it does not use a receptor value.

At `0.05 s-1` capture, `0.01 s-1` return, 0.3 ambient transport fraction,
0.2 m/s asymmetric and 0.1 m/s outward wake velocities, the branch captures
915.75 kg of mixture and 58.72 kg of released gas.  Its normalized vertical
observer is centred at 0.5 m with 0.8 m sigma.  Coupling this branch with an
elevated-path centre of 1.55 m, initial sigma 0.6 m and vertical diffusivity
0.045 m2/s produces the new strict result:

- full 600 s dose MAPE: **28.952452724147292%**;
- independent runtime: **21.88 s**;
- all 145 predicted doses reproduced bit-for-bit;
- height-group deltas versus 29.191%: -0.679, -0.013, -0.004, -0.049 points;
- distance-group deltas: -0.204, -0.423, -0.162 points;
- mixture residual `2.24e-8 kg`, relative species residual `1.28e-15`, and
  momentum gate passed.

The branch improves mean departure error from -69.93 to -64.02 s and mean
duration error from -48.49 to -42.61 s.  It does not yet repair the digitized
peak timing or equivalent-duration ratio, so it is accepted only as the new
engineering benchmark, not as field validation.  Return-rate sweeps from
0.0025 to 0.015 s-1 all worsened the promoted score or a preregistered group;
100 s remains the retained residence timescale.

The new optimistic bounds are 23.947% for mirrored-pair redistribution,
19.417% for vertical redistribution, 22.911% for arc scaling and 17.156% for
tower scaling.  The 10% goal therefore remains unmet and active.

Artifacts:

- `outputs/ti021-mixture-promoted28p952-speciesweighted-fullwake-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted28p952-speciesweighted-fullwake-600s-2026-10-17.configuration.json`
- `outputs/ti021-mixture-promoted28p952-temporal-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted28p952-residual-structure-audit-2026-10-17.json`

## Age-resolved wake and full-carrier prompt-path falsification

The slow wall-wake state now optionally transports a conservative species-age
moment.  Its normalized vertical observer can therefore settle and diffuse as
a function of local residence age rather than applying one fixed Gaussian to
all wake material.  Four full-event vertical diffusivities from 0.005 to 0.04
m2/s increased MAPE from 28.952% to 29.012--29.265% and regressed preregistered
height/distance groups.  The mechanism is implemented default-off but is not
promoted.  The remaining stored wake species had a 570 s mean age at 600 s,
showing that old retained inventory, not missing vertical diffusion alone,
controls this branch.

A second optional mode changes the prompt over-wall path from a species marker
inside the column carrier to an independently transported four-component
state.  It carries mixture mass, released-mixture mass and both horizontal
momenta, plus a species-age moment; capture, transport and return all enter the
same mass and momentum ledgers as the slow wake.  The legacy `species_marker`
mode remains the default.

At the previous coefficients the conserved prompt state captured 431,753 kg
of total mixture and 1,636 kg of released mixture over the event.  It improved
the 0.4 m and 4.4 m group MAPEs by 3.677 and 0.104 percentage points, but
worsened 2.4 m and 6.4 m by 1.750 and 12.753 points; total MAPE was 30.133%.
Crosswind, geometry-outward speed and capture-fraction sweeps did not remove
the vertical/lateral conflict.  The branch is therefore retained only as a
conservative structural experiment, not promoted.

An updated convex-span diagnostic over 272 unique compatible 600 s runs gives
27.620% MAPE even when non-negative whole-run weights are optimized globally.
This is not a validation model, but it proves that mixtures of the existing
response shapes cannot reach 10%; a new spatial/temporal structure is needed.
The next admissible branch must partition prompt wall-contact material between
the column-following shear layer and full-carrier overflow without duplicating
species, then transport the two portions with separate age and vertical
moments.

The post-change legacy reproduction remains exactly 28.952452724147292% for
all 145 sensor doses and runs in 17.87 s.  The focused regression suite is
**57 passed** and the conservation/momentum gates pass.

Artifacts:

- `outputs/ti021-mixture-promoted28p952-post-conserved-branch-regression-600s-2026-10-17.json`
- `outputs/ti021-conserved-prompt-path-vs28p952-600s-2026-10-17.json`
- `outputs/ti021-conserved-prompt-f0p25-vs28p952-600s-2026-10-17.json`
- `outputs/ti021-conserved-prompt-f0p5-vs28p952-600s-2026-10-17.json`
- `outputs/ti021-conserved-prompt-f0p75-vs28p952-600s-2026-10-17.json`
- `outputs/ti021-correct-observable-model-span-audit-2026-10-17-v2.json`

## Conservative prompt/marker hybrid and independent vertical observer

The prompt wall-contact inventory can now be partitioned between the legacy
column-following species marker and an independently conserved four-component
carrier state without duplicating released species.  A 5% conserved split
lowered the full-event MAPE from 28.952% to 28.837%, but was rejected because
some height and distance groups regressed.

The conserved carrier now has independent centre, width, settling velocity,
vertical diffusivity, and minimum-centre settings.  All five settings fall
back to the marker values when omitted.  A full 600 s fallback regression
reproduced all 145 sensor doses bit-for-bit.  The best independent single-mode
candidate used a 0.7 m centre, 0.35 m width, 0.002 m/s settling and zero
vertical diffusion.  It reached **28.469% MAPE**, but regressed the 2.4 m,
6.4 m, and far-distance groups and was therefore not promoted.

The vertical-mode sweep isolated a structural conflict.  A 2.5 m-centred
mode improved 2.4/4.4 m and all distance groups but regressed 0.4 and 6.4 m
by 0.067 and 0.058 percentage points.  A low reattachment mode repaired the
0.4 m group but worsened the far group, while a 6 m mode strongly worsened
6.4 m.  A constrained non-negative mixture of the baseline and all tested
vertical response shapes could not improve the objective while preserving
every preregistered group; the optimizer returned 100% baseline weight.
Therefore another single Gaussian or scalar mixture is falsified.  The next
admissible mechanism is a genuinely height-resolved wall branch whose lower
reattachment and upper overflow inventories have independent horizontal
transport paths.

New diagnostic artifacts:

- `outputs/ti021-hybrid-prompt-f0p05-vertical-split-regression-600s-2026-10-17.json`
- `outputs/ti021-hybrid-f0p05-conserved-kz0-600s-2026-10-17.json`
- `outputs/ti021-hybrid-f0p05-conserved-c0p7s0p35-600s-2026-10-17.json`
- `outputs/ti021-hybrid-f0p05-conserved-c2p5s0p5-600s-2026-10-17.json`
- `outputs/ti021-hybrid-f0p05-conserved-c4s0p8-600s-2026-10-17.json`
- `outputs/ti021-hybrid-f0p05-conserved-c6s0p7-600s-2026-10-17.json`

## Two-path conserved wall-contact transport

The full-carrier prompt state was extended into two independent conserved
wall-contact paths: an upper overflow state and a lower reattachment state.
Both carry mixture mass, released-mixture mass, two horizontal momentum
components, and a species-age moment.  The lower state has independent
ambient-advection blending, asymmetric crosswind speed, geometry-outward
speed and scale, first-order return, centre, width, settling velocity, and
vertical diffusivity.  The upper conserved state also now has horizontal
transport settings independent of the legacy species marker.  Defaults keep
the lower split disabled and fall back to the historical upper-path settings.

The first sweep broke the previous single-Gaussian tradeoff.  A 25% lower
split improved six of seven preregistered groups and missed the 6.4 m gate by
only 0.033 percentage point.  A 75% split improved six groups and missed the
2.4 m gate by 0.063 point.  Subsequent height/return/path sweeps produced
several near-gate candidates around 28.81--28.85% MAPE, but none improved all
seven groups simultaneously.  No candidate was promoted.

The search also exposed an adaptive-step phase artifact: moving an arbitrary
accumulated interval every time a nominal one-second transport became due
caused discontinuous responses to sub-percent parameter changes.  When the
two-path branch is active, upper and lower transport now use separate clocks,
move exactly one declared interval, and retain the residual time.  The
response is smooth after this correction.  The best deterministic two-path
diagnostic reached **28.754% MAPE**, improving six groups but worsening the
6.4 m group by 0.183 percentage point.  Sensor-level inspection shows the
remaining high-level regression is dominated by originally underpredicted
high sensors losing column-following marker inventory, while adding an upper
Gaussian first reaches already overpredicted sensors.  The next admissible
mechanism is therefore a field-based wall-impact partition using local face
height/incident momentum, not another global split coefficient.

The final default-off 600 s regression reproduces all 145 hybrid-control
doses bit-for-bit at 28.8372098664% in 18.77 s.  Relative species residual is
1.46e-15, mixture residual is 2.05e-8 kg, and the momentum gate passes.  The
strict promoted benchmark remains 28.9524527241%; the 10% target remains
active and unmet.

Representative artifacts:

- `outputs/ti021-two-path-final-zero-split-regression-600s-2026-10-17.json`
- `outputs/ti021-two-path-split0p25-u0p1-c0-o0p03-600s-2026-10-17.json`
- `outputs/ti021-two-path-split0p75-u0p1-c0-o0p03-600s-2026-10-17.json`
- `outputs/ti021-two-path-deterministic-c1p5-kz0p001-split0p748-600s-2026-10-17.json`
- `outputs/ti021-two-path-upper-independent-c0p15-o0p095-600s-2026-10-17.json`

Reproduction pattern (an artifact contains the complete accepted run
settings, and unknown report fields are ignored):

```powershell
.\.venv\Scripts\python.exe -c "import json,inspect; from pathlib import Path; from tools.score_ti021_dense_layer_hybrid import run; c=json.loads(Path(r'outputs/ti021-mixture-promoted28p403-two-path-exact-reproduction-600s-2026-10-17.json').read_text()); a=inspect.signature(run).parameters; p={k:v for k,v in c.items() if k in a}; p.update(h=c['horizontal_cell_size_m'],dt=c['max_dt_s'],source_dt=c['source_time_step_s'],entrainment=c['dense_entrainment_model']); run(Path(r'outputs/ti021-promoted28p403-reproduction.json'),**p)"
```

The lower/upper split was then changed from a spatially uniform fraction to a
field-based partition.  At each wall-impact cell the lower share is bounded
by the adjacent wall-height/local-layer-height ratio; the remaining captured
four-component state follows the upper path.  This is a geometry and solved-
field closure and contains no sensor coordinate or observation.  Its best
full-event candidate reached **28.669% MAPE** (with a small slow-wake vertical
diffusivity), improving six preregistered groups but still worsening 6.4 m by
0.236 percentage point.  The effective lower capture share was about 0.29 for
a declared 0.75 maximum, confirming that the new rule materially separated
deep overflow from shallow reattachment.  It was not promoted.

The updated compatible-run span now contains 376 unique 600 s released-
mixture predictions.  Its unconstrained non-negative whole-run convex lower
bound improved from 27.620% to **26.915%**, but remains far above 10% and is
not a validation model.  This proves the two-path family added a genuinely
new response direction while also showing that the remaining error cannot be
removed by recombining the currently available whole-field solutions.

Additional evidence:

- `outputs/ti021-two-path-fieldpartition-split0p75-600s-2026-10-17.json`
- `outputs/ti021-two-path-fieldpartition-wakekz0p001-600s-2026-10-17.json`
- `outputs/ti021-correct-observable-model-span-audit-2026-10-17-v3.json`

## Strict promotion to 28.403% MAPE

The group-constrained whole-run span identified the column-following marker's
vertical diffusivity as complementary to the new two-path transport.  The old
benchmark used `Kz=0.045 m2/s`; after the wall-contact inventory was split,
that inherited value was no longer appropriate.  A full-event sweep at
0.01/0.02/0.03 m2/s found that **0.03 m2/s** is the first value to improve
all four height and all three distance groups simultaneously.

The independently reproduced promoted result is:

- 600 s positive-observed dose MAPE: **28.403478592377457%**;
- runtime including output: **21.01396 s**;
- all 145 predicted doses bit-for-bit identical to the discovery run;
- height-group deltas versus 28.952%: **-0.524, -0.940, -0.015, -0.804**
  percentage point at 0.4/2.4/4.4/6.4 m;
- distance-group deltas: **-0.150, -2.086, -0.317** percentage point;
- relative species residual: **2.38e-15**;
- mixture residual: **-2.24e-8 kg**;
- momentum conservation gate: **passed**.

This result replaces 28.9524527241% as the strict engineering benchmark.  It
does not satisfy the 10% session target and is not field qualification.

Promotion artifacts:

- `outputs/ti021-two-path-fieldpartition-markerkz0p03-600s-2026-10-17.json`
- `outputs/ti021-mixture-promoted28p403-two-path-exact-reproduction-600s-2026-10-17.json`
- `outputs/ti021-group-constrained-model-span-audit-2026-10-17-v1.json`

The promotion decision is now machine-audited by
`tools/audit_ti021_internal_promotion.py`.  The audit deliberately separates
an internal engineering-benchmark promotion from external validation.  It
checks the identical 145-sensor observable, bitwise independent reproduction,
all seven preregistered group gates, species/mixture/momentum conservation and
the 30 s runtime limit in one record:

- `outputs/ti021-promoted28p403-post-trajectory-internal-promotion-audit-2026-10-17.json`

## Post-promotion vertical-trajectory and observer audit

The residual-structure audit on the 28.403% benchmark gives three important
optimistic bounds.  Conservatively optimizing only the left/right split of 19
mirrored pairs gives **23.447% MAPE**; conservatively optimizing only the
height allocation within each horizontal tower gives **19.151%**; and even one
fitted scalar per horizontal tower gives **16.853%**.  These are diagnostic
bounds with the specified predictions and inventories fixed, not admissible
models. Those particular redistributions do not reach 10%, indicating a
horizontal-tower mismatch as well as height allocation error. They do not
rule out a revised transport closure that changes the inventories.

A common decaying post-impact vertical velocity was implemented for both the
conserved upper reservoir and the column-following wall-impact marker:

`z(age) = z0 + w0*tau*(1-exp(-age/tau)) - ws*age`.

Both are default-off, conserve the existing inventory, and preserve the
28.4034785924% benchmark bit-for-bit when `w0=0`.  The focused dense-layer test
suite remains **58 passed**.  Seven full 600 s conserved-branch trials reached
28.351% at best but degraded the 2.4 m group; eight main-marker trials all
failed to beat the zero-rise control.  In particular, lifting the main marker
improved 4.4 m while worsening 0.4 and 6.4 m, so uniform post-impact ascent was
rejected rather than promoted.

The 5 m nearest-cell observer was also replaced by the already implemented
bilinear observer in an otherwise identical replay.  MAPE increased to
**29.191%**, so grid-position quantization is not the missing mechanism.

Finally, a bounded shear-layer development-length sweep changed the common
wall-impact observation ramp from 150 m to 50--200 m.  A 100 m ramp reduced
overall MAPE to **28.181%**, improved the 0.4, 2.4 and 6.4 m groups and the
150--300 m band, but degraded 4.4 m by 0.126 percentage point and the near band
by 0.044 point.  A follow-up vertical-diffusivity sweep could not remove both
regressions.  This candidate is retained as a diagnostic and is not promoted.

Artifacts:

- `outputs/ti021-promoted28p403-residual-structure-audit-2026-10-17.json`
- `outputs/ti021-vertical-trajectory-sweep-2026-10-17.json`
- `outputs/ti021-marker-vertical-trajectory-sweep-2026-10-17.json`
- `outputs/ti021-two-path-bilinear-observer-600s-2026-10-17.json`
- `outputs/ti021-elevated-ramp-sweep-2026-10-17.json`
- `outputs/ti021-ramp100-markerkz-sweep-2026-10-17.json`

The next admissible model change is therefore **spatially selective** rather
than another global vertical coefficient: preserve the successful 100 m
development response only on wall-impact paths that feed the underpredicted
high towers, while retaining the current ground/reattachment path elsewhere.
That selector must be derived from wall-face geometry, incident momentum and
branch age (not sensor coordinates), and must pass the same full-event gates.

## Strict promotion to 28.041% and temporal-structure re-audit

The geometry-outward wall-impact path was narrowed without changing its
conserved inventory.  Reducing the outward transition scale from 30 m to
10 m first produced a strictly admissible 28.099651297843835% candidate.
Every preregistered height and distance group improved, the independent
reproduction was bitwise identical, and the end-to-end runtime including
output was 20.00 s.  A second bounded refinement used a 1 m transition scale
and changed the elevated observation centre from 1.65 to 1.625 m.  The
independently reproduced result was the strict engineering benchmark before
the later wall-diffusion defect was discovered:

- 600 s released-mixture dose MAPE: **28.040686504627093%**;
- runtime including output: **21.02670789998956 s**;
- all 145 predicted doses bitwise identical to the discovery run;
- all four height and all three distance groups improve relative to 28.100%;
- maximum relative species residual: `4.03e-15`;
- mixture residual: `9.31e-9 kg`;
- momentum conservation gate: passed.

This is an internal benchmark promotion only.  It is not external field
validation and remains 18.041 percentage points above the session target.

The updated optimistic diagnostics remain far above 10%: perfect conservative
mirrored-pair redistribution gives 23.210%, perfect conservative vertical
redistribution within each tower gives 18.801%, one fitted scalar per exact
downwind arc gives 22.072%, and one fitted scalar per horizontal tower gives
16.685%. These fixed-prediction redistributions do not meet the target;
they do not establish a lower bound for untested physical transport models.

The new time-history audit identifies the next structural error.  Across 82
timed positive sensors, mean arrival error is -22.43 s.  Across 67 complete
arrival/departure pairs, mean departure and duration errors are -55.13 and
-33.72 s.  The eight independently digitized traces have a mean peak-time
error of -32.13 s, a geometric-mean peak ratio of 1.452, and a median
equivalent-duration ratio of 0.561.  The model is still too early, too peaked,
and too narrow.

Falsification checks after this promotion rejected the following global
remedies:

- lowering the effective layer wind from 0.50 to 0.35--0.45 m/s worsened MAPE
  to 28.104--28.193%;
- increasing isotropic horizontal diffusivity from 2.5 to 3.5--7.5 m2/s
  worsened MAPE to 29.418--35.692%;
- using measured wind-speed fluctuations with direction, or speed alone,
  worsened MAPE to 29.147% and 29.441%;
- reducing the mast-record propagation speed below 3.2 m/s is not admissible
  with the authoritative 0--850 s record because the complete 600 s spatially
  shifted interval is no longer covered.

The next admissible mechanism is a conservative multi-path longitudinal
velocity distribution attached only to wall-impacted material.  Upper
overflow, lower reattachment, and slow separated-wake inventories must carry
mixture, released species, both horizontal momenta, and age moments without
duplication.  The iteration is accepted only if it improves the digitized
peak time and equivalent duration first, then lowers the full 600 s dose MAPE
without degrading any preregistered group, while retaining the 30 s runtime
and conservation gates.  Sensor IDs and observed-dose multipliers remain
forbidden.

Promotion and audit artifacts:

- `outputs/ti021-mixture-promoted28p100-outscale10-reproduction-600s-2026-10-17.json`
- `outputs/ti021-promoted28p100-internal-promotion-audit-2026-10-17.json`
- `outputs/ti021-mixture-promoted28p041-outscale1-center1p625-reproduction-600s-2026-10-17.json`
- `outputs/ti021-promoted28p041-internal-promotion-audit-2026-10-17.json`
- `outputs/ti021-promoted28p041-residual-structure-audit-2026-10-17.json`
- `outputs/ti021-promoted28p041-temporal-residual-audit-2026-10-17.json`

## Wall-aware horizontal diffusion: repair and full-event verification

The horizontal diffusion stencil formerly exchanged mixture mass, released
species and horizontal momentum across all adjacent cells. Unlike advection,
it did not consult the wall-face heights. A dedicated solver regression
places gas on one side of a 100 m solid wall with gravity, drag and
entrainment disabled: the legacy stencil puts species on the other side;
the repaired stencil does not. Thus this is a reproduced transport defect,
not a conclusion inferred merely from the sensor score.

`horizontal_diffusion_respect_walls=True` uses shared nonnegative face
conductances and antisymmetric conserved updates. For a finite-height wall,
the open fraction is the cloud height above the wall divided by the local
maximum column height, including the configured blockage/porosity factor.
The same conductance is used for mixture, species and both signed horizontal
momenta. This is a reduced uniform-height geometric closure, not resolved
3-D overflow. The legacy default remains false solely for exact historical
replay; all new repaired-path trials must explicitly set it true.

| Complete 600 s replay | Dose MAPE | Output-inclusive runtime | Decision |
|---|---:|---:|---|
| Legacy stencil, latest source | 28.0406865% | 24.17 s (concurrent) | historical reproduction; known boundary gate fails |
| Wall-aware diffusion, independent isolated replay | 29.0965237% | 22.2871 s | physically necessary repair, not accuracy promotion |
| Wall-aware diffusion, cavity exponent 16 | 28.9310914% | 25.00 s | reject: ground-height group worsens 3.107 points versus repaired control |

The independent repaired replay exactly reproduces all 145 doses and the
saved concentration-history NPZ SHA-256 of the discovery run. Maximum
relative species residual is `2.0112e-15`, mixture residual is
`-1.86265e-8 kg`, and the momentum conservation gate passes. The reported
2335.77 kg blocked diffusive throughput sums absolute hypothetical prevented
face transfers; it is not a unique amount of released gas and can count
material more than once.

The repair improves time-history diagnostics but reveals a near-field
vertical mismatch. Mean arrival error changes from -22.43 to -20.29 s,
mean departure error from -55.13 to -34.00 s, and mean duration error from
-33.72 to -14.86 s. Digitized peak-time MAE decreases from 32.88 to 30.25 s;
the median equivalent-duration ratio rises from 0.561 to 0.620, still too
narrow. Arrival/departure estimates use a diagnostic 5% relative threshold,
not a confirmed experimental threshold definition.

Height-group MAPEs after repair are 29.954%, 22.360%, 32.085% and 33.822%
at 0.4, 2.4, 4.4 and 6.4 m. The 4.4 m group worsens by 5.319 percentage
points, mainly at near-source towers. This directs the next investigation
to the near-field vertical distribution and branch residence, not to
reintroducing unphysical diffusion through the fence for a lower score.

Three bounded full-event sweeps were evaluated:

- Entrainment 0.25, 0.275, 0.325, 0.35 and 0.40 all worsen the repaired
  0.30 control (29.667--36.003%).
- A diagnostic Richardson-dependent suppression of the vertical observation
  tail at coefficients 0.1, 1, 5 and 20 also worsens MAPE
  (29.114--29.317%). This speculative option was removed from source; its
  archived results are rejection evidence, not currently replayable options.
- Cavity exponents 10, 16, 20 and 24 trade high-level error against ground
  error. Exponent 16 is the lowest total at 28.931%, but fails the group gate.

### Strengthened acceptance audit

`tools/audit_ti021_internal_promotion.py` schema v2 now verifies the original
workbook sensor identities/observations, separate candidate/reproduction
files, all 145 complete 0--600 s traces, bounded concentrations, reintegrated
doses, matching observable, frozen baseline route groups, and the known
wall-diffusion boundary requirement in addition to the previous mass,
momentum, height/distance, runtime and bitwise reproduction tests.

The known-boundary check only covers this demonstrated defect; passing it
does not certify all physical boundary conditions. The earlier
`ti021-promoted28p041-full-trace-audit-2026-10-17.json` predates this extra
gate and is superseded by the legacy boundary re-audit below.

Focused verification: **65 tests passed** across the dense-layer solver and
promotion-audit suites; Python compilation and `git diff --check` passed.

Authoritative artifacts:

- `outputs/ti021-legacy28p041-boundary-reaudit-2026-10-17.json`
- `outputs/ti021-wall-diffusion-repair-promotion-audit-2026-10-17.json`
- `outputs/ti021-wall-diffusion-independent-reproduction-600s-2026-10-17.json`
- `outputs/ti021-wall-diffusion-temporal-audit-2026-10-17.json`
- `outputs/ti021-wall-diffusion-entrainment-sweep-2026-10-17.json`
- `outputs/ti021-wall-diffusion-stable-tail-sweep-2026-10-17.json`
- `outputs/ti021-wall-diffusion-core-profile-sweep-2026-10-17.json`

Reproduce the repaired control from the repository directory (use a new
output filename to preserve the independently audited artifact):

```powershell
.\.venv\Scripts\python.exe -c "import json;from pathlib import Path;from tools.sweep_ti021_vertical_trajectory import _base_parameters,_execute;b=json.loads(Path('outputs/ti021-wall-diffusion-finite-height-600s-2026-10-17.json').read_text());f='outputs/ti021-wall-diffusion-local-replay-600s.json';_execute((f,_base_parameters(b)));r=json.loads(Path(f).read_text());print(r['dose_mape_percent'],r['runtime_including_output_s'])"
.\.venv\Scripts\python.exe -m pytest tests/test_dense_layer_patch.py tests/test_audit_ti021_internal_promotion.py -q
```

### Next bounded iteration

1. Keep the wall-aware repair enabled and use the 29.0965% repaired control
   for mechanism ablations, while retaining 28.0407% only as a historical
   comparison. Do not discard the necessary boundary repair to chase MAPE.
2. Audit how the near-source collapsing core, wall-impact marker and upper
   branch allocate a single conserved species inventory in height. In
   particular, test the existing non-mass-normalized observation tail and
   concentration clipping before adding another tunable profile.
3. If that audit shows a real allocation defect, repair it first. Otherwise
   implement geometry/momentum/age-driven exchange between ground,
   overflow and slow-wake material; do not use sensor coordinates or
   observed-dose multipliers in that exchange.
4. Check isolated invariants, then recompute the complete 600 s event.
   Report global and height/distance/route MAPE, timing diagnostics,
   conserved inventories and output-inclusive runtime together. Reject
   regressions explicitly rather than declaring a failed run validated.
5. Accept a new numerical candidate only after an independent replay and
   the updated gates pass. The final target stays <=10%; reaching it on
   this repeatedly used experiment alone would still require a separate
   held-out case before claims of predictive field validation.

## Species-volume-conservative observer and bounded correction cycle

### Implemented repair and scope

At fixed temperature and pressure, reconstructed mole fraction must satisfy
`integral X(z) dz = species_mass_per_area / pure_species_density`.
Multiplying a mass fraction by a normalized height profile and then applying
a nonlinear mass-to-mole conversion does not satisfy this identity. Clipping
the final mole fraction at one can also discard reconstructed inventory.
The old solver's mass ledger cannot see either observation-only error.

`src/degali/addons/conservative_column_observation.py` now integrates
normalized polynomial, finite exponential/erfc tail and positive-half-space
Gaussian branch kernels into vertical bins. It combines their conserved
species volumes before capacity-constrained redistribution. Concentration
is bounded between zero and one with the same species-volume integral;
infeasible support capacity raises an error instead of silently deleting gas.
The discretization depends on column/branch geometry, not sensor elevation.
Near-coincident kernel edges are merged at roundoff scale, and limited
cell-linear sampling preserves each bin's mean and its concentration bounds.

The dense-layer solver exposes this through
`vertical_observation_mode='conservative_volume'` and
`vertical_observation_bins`. The former mode remains explicitly available as
`legacy_mass_fraction` for historical replay and is still the backward API
default. New conservative trials select the repaired mode explicitly through
the saved-configuration Python API. This does not silently change every old
runner or claim that all review findings have been fixed.

The main marker's upper share is subtracted from its lower share, not counted
twice. Conserved upper/reattached material stays in the lower profile before
its development ramp reaches the upper branch; the ramp no longer makes that
inventory invisible. Existing optional hidden-store visibility fractions
remain explicit: the new volume gate covers represented columns, not an
independent proof that every global reservoir is physically visible.

This remains an isothermal reduced-column observation model. It does not add
a vertical momentum equation or establish resolved ascent, roof overflow or
3-D wake validation. The source trajectory and horizontal transport state
were not changed by this repair. Direct tests verify identical solver states
with the old and new observer while independently integrating sampled
concentrations back to species inventory.

### Full-event results

| Path | Dose MAPE | End-to-end runtime | Decision |
|---|---:|---:|---|
| Wall-aware solver, legacy observer, latest-code replay | 29.096524% | 20.80 s | bitwise historical replay; observer defect retained |
| New observer, original transport parameters, 32 bins | 29.979793% | 28.19 s | correction control; no accuracy promotion |
| New observer, entrainment 0.25, vertical-tail Kz 0.02 m2/s | **28.592144%** | **25.16 s independent replay** | total improves, group gate fails |
| Same entrainment, Kz 0.04 m2/s | 29.395580% | 26.60 s | reject |
| Entrainment 0.275, Kz 0.02 m2/s | 29.619476% | 24.51 s | reject |

All runs compute the full 0--600 s event and all 145 receptors; MAPE uses the
same 84 positive observations. The best candidate changes no sensor-specific
coefficient or observed-dose multiplier. Relative to the wall-aware legacy
observer control, its height-group MAPE changes are -3.410, +2.765, -1.207 and
+1.355 percentage points at 0.4, 2.4, 4.4 and 6.4 m. The near, middle and far
distance-group changes are -1.172, +0.363 and +0.930 points. These are explicit
reasons for rejection, despite the 0.504-point total improvement.

The independent best-candidate replay has relative species mass residual
`1.8283e-15`, mixture residual `-1.86265e-9 kg`, a passing horizontal momentum
gate, and represented-column species-volume residual `6.2379e-16`.

### Numerical refinement and discarded prototype rankings

The initial piecewise-constant prototype had a 12.535% maximum sensor dose
change (normalized by the positive observed dose) from 64 to 128 bins, even
though global MAPE barely moved. Its coarse-profile and entrainment sweeps
are retained as diagnostic artifacts, not authoritative candidate rankings.
The edge merge and limited-linear reconstruction remove most of this
sampling sensitivity:

- original-parameter control, 32 to 64 bins: maximum normalized sensor dose
  change 0.355%; MAPE change 0.00731 percentage point;
- control, 64 to 128 bins: maximum normalized sensor dose change 0.149%;
  MAPE change 0.00127 point;
- best coupled candidate, 32 to 64 bins: maximum normalized sensor dose
  change 0.227%; MAPE 28.592144% versus 28.587501%.

These are observation-resolution checks, not horizontal-grid or full CFD
convergence. Finer observation bins are not claimed to satisfy the 30 s gate
on every run; the independently replayed 32-bin candidate does.

### Tests, audit and authoritative artifacts

**144 tests passed** across the conservative-column, dense-layer,
promotion-audit, binary-puff, finite-prism, binary-thermodynamic,
HDF5-to-layer and TI021 suites. Compilation and `git diff --check` passed.
The broader 2,338-test review result from the previous turn is not presented
as a rerun of all changes in this iteration.

Promotion audit schema v3 now requires explicit conservative observation mode
and finite represented-volume residual <=1e-9, in addition to the existing
wall, full-event trace, original workbook identity, reproduction,
mass/momentum, group and runtime gates. A zero solver mass residual on a
legacy observation path cannot satisfy the new requirement.

- `outputs/ti021-conservative-volume-iteration-summary-2026-10-17.json`
- `outputs/ti021-conservative-volume-linear32-control-600s-2026-10-17.json`
- `outputs/ti021-conservative-volume-linear32-e025-k002-600s-2026-10-17.json`
- `outputs/ti021-conservative-volume-e025-k002-reproduction-600s-2026-10-17.json`
- `outputs/ti021-conservative-volume-e025-k002-refinement64-600s-2026-10-17.json`
- `outputs/ti021-conservative-volume-e025-k002-promotion-audit-2026-10-17.json`
- `outputs/ti021-legacy-observer-post-volume-repair-regression-600s-2026-10-17.json`

Reproduce the best diagnostic without overwriting its audited artifacts:

```powershell
.\.venv\Scripts\python.exe -c "import json;from pathlib import Path;from tools.sweep_ti021_vertical_trajectory import _base_parameters,_execute;b=json.loads(Path('outputs/ti021-conservative-volume-linear32-e025-k002-600s-2026-10-17.json').read_text());_execute(('outputs/ti021-volume-candidate-local-replay.json',_base_parameters(b)))"
.\.venv\Scripts\python.exe -m pytest tests/test_conservative_column_observation.py tests/test_dense_layer_patch.py tests/test_audit_ti021_internal_promotion.py tests/test_binary_integral_puff.py tests/test_dense_prism_integral.py tests/test_binary_mixture_thermodynamics.py tests/test_hdf5_dense_layer_packets.py tests/test_thorney_ti021.py -q
```

### Next action, without changing the goal

Keep wall-aware diffusion and species-volume-conservative reconstruction
enabled. The new 28.592% result is a diagnostic starting point, not a
nondegrading promoted benchmark. Reinspect the 2.4/6.4 m residuals with branch
age/height contributions separated, then change the physically transported
ground/overflow/wake partition or entrainment history where that evidence
supports it. Do not restore the old concentration inflation to recover a
lower number. Repeat the same complete event and group/reproduction gates.
The other review findings (common statistics, exterior conversion/velocity,
legacy wall defaults, multilayer initialization and inherited core time
inversion) were not modified in this iteration and must not be reported fixed.

## Whole-model audit and corrected-wall benchmark promotion

The complete review found that finite-wall openness for a vertically
ground-loaded column had been raised to `p` instead of the analytically
integrated `p+1`.  The corrected law is covered by exact tests at `p=0,2,3`.
At the former `p=3` setting, this repair changed the best candidate from
28.592144% to 29.123884%, demonstrating that the old score had benefited from
excess overtopping.

The configuration provenance was then corrected to the final-report priority
inputs: wall height 5 m and source wind 3.9 m/s rather than the secondary
workbook alternatives 6 m and 3.4 m/s.  With the conservative observer and
wall-aware diffusion retained, this produced a current-code benchmark of
28.201528% in 24.70 s.

Because fixing the wall integral changes the physical meaning of the exponent,
the ground-loading exponent was re-evaluated without sensor-specific factors.
Using `p=2` (integrated openness power 3) reduced MAPE to 27.695521%.  Its
2.4/6.4 m regressions were 0.016/0.040 percentage point.  Matching the
wall-impact marker's initial sigma to the conserved upper branch at 0.50 m
removed those regressions.  The independently reproduced strict result is:

- dose MAPE: **27.6831474975817%**;
- output-inclusive replay runtime: **27.9596969 s**;
- all four height and all three distance groups nondegraded;
- complete 0--600 s trace and dose reintegration: pass;
- species/mixture/momentum and conservative-column gates: pass;
- candidate/replay time and concentration arrays: bitwise identical;
- implementation hashes: current and matching.

This is promoted only as the next internal engineering benchmark.  It remains
17.683 percentage points above the <=10% objective and does not establish
external validation.  The component audit of the immediately preceding
benchmark attributes 97.31% of predicted dose to the unresolved core, so the
next useful mechanism must change conservative core timing/vertical transport
rather than amplify the small upper/wake observation branches.

Artifacts:

- `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-600s-2026-10-17.json`
- `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-independent-600s-2026-10-17.json`
- `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-promotion-audit-2026-10-17.json`
- `docs/degali-whole-model-error-audit-2026-10-17.md`

## Wall-flux provenance and upper-carrier rejection

The current implementation can now combine true finite-wall crossing flux
with the hybrid conserved impact split.  A child species/age marker remains a
subset of the ordinary elevated marker through base FV transport, relative
transport, local cavity capture/return and wake capture.  It therefore
separates physical over-wall provenance from the impact-lift Gaussian without
creating mass or using receptor gains.  Eighty-seven focused tests pass.

The mechanism is active (419--427 kg marked in the full event), but all tested
600 s closures were rejected: the shared crosswind/Gaussian candidate reached
27.2962% but regressed four preregistered groups; a 5.5 m Gaussian gave
31.6872%; the core-kernel/shared-crosswind version gave 27.7536%; and the
core-kernel/5 m neutral-log carrier gave 28.2162% in 31.28 s.  Consequently the
promoted configuration keeps wall-flux marking disabled.

After these code changes the promoted path was independently replayed.  The
sensor arrays remain bitwise identical at 27.6831474975817%, the independent
output-inclusive runtime is 28.2565 s, all conservation and height/distance
gates pass, and the refreshed audit is:

- `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-currentcode2-promotion-audit-2026-10-17.json`

The current residual bounds are 22.782% for mirrored-pair redistribution,
19.245% for tower-wise conservative vertical redistribution, 21.460% for one
scale per downwind arc and 15.914% for one scale per tower.  Thirteen of 21
multi-height positive-dose towers have non-monotone arrival ordering.  The
next mechanism is therefore a conservative time-lagged vertical
recirculation/core state, not another global vertical kernel or overtop speed.

## Whole-model audit follow-up: rejected trajectory and repaired defects

An isolated age-dependent vertical trajectory for actual wall-crossing flux
was run for the complete 600 s event using baseline-preserving parameters.  It
scored 30.11864678% MAPE in 30.9594 s and worsened three of four height groups,
so the experimental mode was removed.  Earlier trajectory files produced by
`run_ti021_residence_candidate.py` are excluded because that CLI overwrites
multiple baseline fields with parser defaults; only the `_base_parameters`
control reproduced the promoted doses exactly.

Increasing and broadening the existing conserved impact split gave a useful
but non-promotable diagnostic.  The 20%, 3.5 m centre, 1.5 m sigma candidate
scored 27.05109672% in 28.0146 s, but regressed the 2.4 and 6.4 m groups by
0.1337 and 0.3390 percentage point.  The 40%, 4 m centre, 2 m sigma candidate
scored 27.11384369% and regressed 6.4 m by 2.8572 points.  This confirms that a
single broader Gaussian cannot supply the missing height-dependent timing.

The whole-model review also found that the cavity-backflow target was written
in global x coordinates and used `wind[0]`, while `cavity_align_with_wind`
rotated only the mask.  A transpose/90-degree witness was exact when backflow
was off but differed by 12.41--13.85% in mixture and 29.73--30.68% in momentum
when backflow was active.  The closure now uses a consistent local basis and
passes an exact rotation-covariance regression.  The corrected current-code
full-event score is 27.75955662% in 29.586 s, slightly worse than the
27.68314750% historical frame, so it is the default physical implementation
but not a promoted accuracy candidate.  Explicit `legacy_x` replay remains
available only to reproduce the old benchmark bit-for-bit.

The residence-candidate runner's implicit parser-default overrides were also
removed.  A repaired-runner 600 s control exactly reproduces baseline times,
concentrations and doses with zero parameter disagreements.  The legacy
transient supervisor now retains observer failure stage/type/message and the
CLI warns instead of allowing a caught numerical failure to look like a valid
zero-concentration result.

## Multilayer execution, entrainment and handoff audit

The current scorer's conservative-volume default was incorrectly rejected by
the multilayer input guard even though the multilayer branch has its own EOS
volume-occupancy observer.  The guard and regression test are repaired.  The
multilayer solver also documented top-surface ambient entrainment but placed
all entrained mass in the ground layer.  Corrected runs use the highest
materially occupied layer, with `legacy_ground` available only for historical
A/B replay.

The scorer then exposed a third integration defect: the integral model had
already collapsed the 13 m finite prism to 6.57593994 m at handoff, but the
multilayer initializer was given the original 0--13 m source bounds and
silently re-expanded it.  Initialization now derives its vertical geometry
from the actual handoff state and records it in the report.  The old
79.510594% (13-layer) and 78.196174% (three-layer) results are therefore
pre-repair history, not evidence about the current path.

The corrected seven-layer uniform-handoff 600 s run scores **56.95036658%**
in **24.7616 s** and passes species, mixture, horizontal momentum and vertical
momentum conservation.  Gaussian initialization scores 63.95--67.94%, while
simple global settling/lift variants worsen the uniform result to
59.42--68.27%.  The run also cumulatively projects 318.8958 kg of species
between height layers to restore local `species <= mixture` bounds.  Total
species is conserved, but that correction is too large to interpret as mere
round-off.  The current single-layer physical result remains bitwise unchanged
at **27.75955662%**.  Thus global multilayer replacement is rejected: the next
admissible route is a conservative local obstacle/wake vertical exchange
coupled back into the fast single-layer outer carrier, not more global
vertical bins or fitted sensor gains.

## Final whole-model audit reproduction

After all repairs, the current source was used to regenerate the 600 s
single-layer baseline, corrected cavity-frame candidate and an independent
candidate reproduction.  The candidate reproductions are bitwise identical
and score 27.75955662% dose MAPE.  The second run took 28.9963 s including
output, and species, mixture, momentum and conservative vertical-volume gates
pass.  Relative to the fresh 28.20152797% baseline, aggregate MAPE improves by
0.44197134 percentage point, but the 6.4 m group regresses by 0.31672472 point.
The strict audit therefore rejects promotion; the 10% objective remains
unmet.  Evidence is saved in
`outputs/ti021-cavity-frame-cavity-currentcode6-current-baseline-promotion-audit-2026-10-17.json`.

After the handoff-geometry repair and its two regression tests, the final
repository-wide suite completed with **2,430 passed and 145 skipped in
1,134.88 s**, with no failures.  This establishes implementation consistency
for the tested surface, not external physical validity.

## Persistent local wake and vertical-observer iteration

The local wake store was found to stop advecting when its return rate was
zero. Transport has been separated from return, and a height-based neutral
log-wind speed is available for a declared wake height. The repaired state
passes 78 focused tests and exact conservation checks.

The complete-event sweep did not close the target. The lowest dose result
that improves every preregistered height and distance group is
**26.77419422%**, obtained with a 0.05 s-1 below-wall/core capture, 0.01 s-1
return, and a visible local wake transported at 0.2 of the reduced layer wind.
It is rejected because the mean departure error worsens to -41.45 s and two
current-code replays took 30.84 and 30.80 s. Making the store persistent and
independently mobile repairs
mean departure timing to -1.66 s, but the conservative layer observer scores
27.56005890% and regresses the 2.4 m and near-distance groups. A 66.7 s return
time reaches 26.89297295% at 30.01 s but leaves small 2.4 m and far-distance
regressions. Neutral-log transport at 2.5 m does not remove the tradeoff.

This falsifies another single-state closure: one vertical distribution cannot
fit dose allocation and pulse duration together. The 10% goal remains active;
no sensor-specific factor, hidden mass or timing-only promotion is adopted.

## 2026-10-17 component-span and sequential-return checkpoint

- Exact component attribution of the 26.77419422% state closes to
  `2.274e-13 vol%-s`; the core contributes 96.31% of positive-sensor dose.
- Global 0--5x component amplitudes have an exact lower bound of 24.8962%, so
  branch strength tuning cannot reach 10%.
- A conservative local-reservoir -> reattached-carrier split was implemented.
  The 0% control reproduces the baseline exactly.  A 20% split with a 1.7 m
  centre gives a repeatable 26.68333974%, improves every registered height and
  distance group, and closes species to `2.011e-15` relative error.
- The candidate remains conditional because its isolated output-inclusive
  runtime is 32.3572 s and the <=10% goal is not reached.
- Full-return fast-carrier, wind-history speed 3.9 m/s and +/-3 degree wind
  direction candidates were rejected.  The next iteration must alter the
  conservative core's spatially varying vertical/crosswind transport, not add
  another global branch multiplier.

## 2026-10-18 full-horizon follow-up checkpoint

- A broad obstacle-wake meander candidate (`wake_length_m=150`,
  `wake_base_width_m=30`, `wake_meander_strength=0.30`) scored **28.7741%**
  in **27.59 s** and was rejected.
- A two-layer vertical-shear candidate conserved species and mixture mass but
  scored **60.7351%** in **24.35 s**; every registered height and distance
  group degraded.  Global multilayer replacement is therefore not admissible.
- Finite-prism continuous-release candidates scored **32.7954%** for 60 s and
  **38.0132%** for 10 s, so release-duration broadening does not explain the
  residual.
- Increasing the common diffusive-tail extension from 10 m to 20 m produced
  **26.15855%** in **29.21 s**, but degraded the 4.4 m and 150--300 m groups;
  the 15 m trial scored **26.16564%** and was also rejected.  The best
  individual remains 26.15868% and is not promoted.

These complete-event tests further narrow the missing mechanism to a coupled
near-field 3-D source/obstacle split and spatially varying core transport.
No candidate in this checkpoint uses sensor-specific gains, and none is
promoted unless the conservation, group, runtime and 10% gates all pass.

## 2026-10-18 conservative core-meander closure

A new optional carrier closure was implemented for obstacle-shed lateral
momentum across the dense core (`core_meander_strength`,
`core_meander_length_m`, `core_meander_width_m`).  It changes only the
carrier velocity; mixture, species and momentum still advance through the
same conservative FV fluxes.

- Strength 0.10, length 150 m, width 50 m: **26.15233%**, 28.87 s;
  species conservation passed, but the 6.4 m and far-distance groups
  degraded.
- Strength 0.05 with the same geometry: **26.14946%**, 33.65 s; the runtime
  gate failed and the 6.4 m group still degraded.
- Strength 0.10 with 180 degree phase: **26.18378%**, 28.91 s; rejected.

The closure remains available as an explicitly optional diagnostic, but none
of these runs is promoted.  The best accepted-path evidence remains above the
10% target; the next candidate must couple the near-field 3-D source/obstacle
split to the outer core without violating the group and runtime gates.

## 2026-10-18 named 3-D packet coupling checkpoint

The existing named donor packet was coupled through the conservative interface
at three downstream locations using the same 600 s replay and no sensor
calibration.  Relative to the 26.15868% reference, the results were:

| interface x | MAPE | runtime | decision |
|---:|---:|---:|---|
| 100 m | 44.3493% | 39.63 s | reject |
| 160 m | 37.7435% | 52.54 s | reject |
| 220 m | 35.1191% | 51.72 s | reject |

The donor packet conserves its named mixture/species/horizontal-momentum
ledger, but explicitly discards vertical impulse and sensible enthalpy and
does not couple lateral/top exports.  It therefore cannot pass the full
3-D handoff gate.  A wide-inner packet without the named donor schema was
also rejected at load time rather than silently reused.

An explicit schema-conversion audit on the wide-inner packet then produced
**50.4307%** when the full ambient normal flux was restored (the marked
mixture ledger became **5.86e7 kg**) and **56.2863%** when only positive
anomaly flux was retained.  Both conversions are rejected: the former
double-counts ambient carrier mass and the latter omits the height-resolved
energy and vertical-momentum coupling required by the receiver.
