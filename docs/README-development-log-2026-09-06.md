# degali

> **1차 완료 — 2026-09-06:** 실제 실행 경로와 성능 판정을 확정했습니다.
> 지면 반사 중복을 제거한 별도 `BASE_SINGLE_IMAGE` 실행 경로를 제공하며,
> 7개 시험의 새 적분·38/17/41 비교를 재현했습니다. 최신 열분포 후보는
> 온도 개선에도 농도 편향·높이 오차가 악화되어 기본판으로 승격하지 않습니다.
> 최종 134개 파일/환경 검증 통과. 아래의 진행 중 표시는 과거 기록입니다.
> [1차 결과](docs/stage1-results-2026-09-06.md) ·
> [실행 방법](docs/stage1-reproduction-2026-09-06.md).
> 전체 설비 검증/미완성 난류 연구가 끝난 것은 아니며 2차를 자동 시작하지 않습니다.

> **05:47 KST (2026-09-06):** source matching now includes prescribed Q,
> axial normal momentum and stress work. Allowing source centers/area/velocity
> and scalar widths to move recovers both registered synthetic cases with
> unchanged original6targets (max1.489e-11);105 hashes match. The failed old
> fixed-center large-Q case is preserved. Actual turbulence inputs, pressure/
> full stress transport and measured-distance validation remain open. Full
> non-slow regression is running. [Latest results](docs/finite-tke-transport-results.md).

> **05:00 KST (2026-09-06):** opt-in finite-TKE transport now separates shear
> production, turbulent storage and thermal dissipation. All retained equations
> and independent actual H+meanKE+Q differences pass a manufactured audit.
> Its synthetic closure inputs FAIL physical gates and are NOT adopted. Full
> non-slow425 passed; no observed accuracy/default claim. Initial source energy
> reallocation also passes its manufactured audit (05:10);4 later source tests
> passed separately. Next physical initial/boundary closure. [Results](docs/finite-tke-transport-results.md).
> Older blocks below are historical snapshots.

> **04:30 completion checkpoint (2026-09-06):** the re-solved exact-geometry
> local conservative witness passes boundaries, all retained equations and
> NEW-direction actual energy checks. All7 initial shear/TKE bound diagnoses
> resolve (minimum flux10.51–11.72% of mean axial KE). These are NOT a physical
> circulation/TKE closure or observed accuracy. Next: finite-TKE modal operator.
> [Current results](docs/exact-geometry-and-tke-results.md). All audits below ended.

> **Latest numerical/physics result (2026-09-06,04:10 KST):** coherent exact
> width/wind geometry passes independent actual full-energy differences
> (max7.18e-8). Re-solving the local conservative flow is in progress.
> A coefficient-free shear-realizability bound requires TKE flux at least
>10.56% of mean axial KE for the prior aligned trial10 witness; it is not
> an actual k value or a closure. Non-slow403 passed; no new field/default claim.
> [Current results](docs/exact-geometry-and-tke-results.md). Older blocks are history.

> **Latest overnight result (2026-09-06):** a direction-preserving conservative
> circulation witness passes trial10's local checks, but is NOT a turbulence
> closure. Independent energy differences expose cancellation and width/wind
> derivative consistency problems; all failures retained. No observed score or
> default change. Full non-slow385 passed. [Current results](docs/aligned-flux-and-energy-conditioning-results.md).
> Entries below are historical snapshots, including then-running audits.

> **Overnight conservative-flow research (2026-09-06):** local initialization
>7/7, guarded short segments6/7; trial10 heat boundary failure independently
> reproduced. Two hard boundary prototypes were rejected. A small solenoidal
> flux FEASIBILITY witness reaches H2/heat/P3.54/4.00/3.02% at the failed
> field, but is not a stable closure or observed validation. Independent
> phase-volume verification is running. Non-slow371 passed.
> [Current results](docs/conservative-flux-witness-results.md),
> [seven-case segments](docs/enriched-short-segment-results.md).

> **Seven-case initialization extension (2026-09-06):** all7 now pass
> independent local conservation, H2/heat/momentum boundary and sampled
> positive-viscosity gates.85 frozen hashes match; non-slow355 passed.
> Guarded downstream segments are under test: trial10 leaves the heat
> boundary gate at about0.44 mm, so local success is NOT downstream success.
> No new field score or default change. [Results](docs/coupled-initialization-extension-results.md).
> Older entries below are historical snapshots.

> **Simultaneous shape/transport initialization (2026-09-06):** both
> preregistered pilots25/11 pass independent local conservation, numerics,
> scalar/momentum boundary and sampled positive-viscosity gates. Trial25
> H2/heat defects fall16.88/16.54% to3.24/2.30%; trial11's negative inferred
> viscosity becomes positive. No new physical coefficient, downstream run,
> observed score or default promotion. Non-slow342 passed.
> [Results](docs/coupled-shape-initialization-results.md). Other five not rerun.

> **Coupled enriched-shape transport (2026-09-05):** gauge-free modal rates
> now reconstruct mass flux, stress work and scalar transport on the new shapes.
> Numerics pass7/7 with separately verified fine mixing input, but no case
> passes adoption gates: renewed boundary defects and negative inferred viscosity in
> two cases. No downstream run, observed accuracy or default change.
> Non-slow:336 passed. [Results](docs/enriched-shape-transport-results.md).

> **Conservative edge refit (2026-09-05):** all seven fixed-transport square
> sections now pass independent conservation, edge and physical checks.
> Face/ray phase splitting and constrained refitting reduce internal H2/heat
> edge defects from12–25% to0.95–3.63%. Earlier failures remain unchanged.
> This is not new downstream transport or measured accuracy; defaults are
> unchanged. Non-slow:324 passed. [Results](docs/edge-conservative-refit-results.md).

> **Square edge-profile enrichment (2026-09-05):** new two-dimensional H2/heat
> shapes reduce both fixed-transport edge defects below5% in3/7 cases, but
> **0/7 pass independent conservation checks**. Candidates are not adopted.
> Phase-cell diagnosis identifies under-resolved optimization integrals;
> further angular error control and conservative refitting are required.
> Non-slow:317 passed. Defaults/field scores remain unchanged.
> [Results](docs/edge-profile-enrichment-results.md).

> **Ambient-reservoir thermal moments (2026-09-05):** an explicit inflow-energy
> boundary now closes a research weak thermal-width transport law. Seven
> millimetre-scale segments have passed conservation/convergence checks,
> including separate adaptive verification of one failed fixed-step case.
> Pointwise species/heat edge defects remain 12–25%; no field accuracy or
> default promotion is claimed. [Results](docs/reservoir-thermal-moment-results.md).

> **Shear/thermal closure screen (2026-09-05):** reduced shear-work and heat
> moment coupling are implemented. Numerical checks pass 7/7, but unity scalar
> diffusion plus immediate shear heating fails the independent heat balance
> in all seven cases. The finite cold edge does not match the existing ambient
> energy-inflow condition. No new downstream run, accuracy claim or default.
> [Results and next boundary correction](docs/shear-thermal-compatibility-results.md).

> **Conservative transverse mixing (2026-09-05):** the independent thermal
> section now reconstructs transverse mass/species transport conditional on
> its unknown width rate. All seven boundaries pass independent flux-change
> verification; no width rate or thermal/species diffusivity ratio is fitted
> or selected. This is not yet a closed downstream thermal ODE.
> [Results and remaining closure](docs/transverse-conservative-mixing-results.md).

> **Thermal-shape transport operators (2026-09-05):** mean enthalpy second
> moments, moving/curved free-section balances and specific-enthalpy-gradient
> response are implemented. Finite-edge terms are material; phase-gradient
> quadrature needs separate error control. Adaptive repair passes all seven
> independent fixed-section checks. These are operators, not a closed
> downstream width law or a new accuracy claim. Defaults remain unchanged.
> [Results and remaining physics](docs/thermal-moment-operator-results.md).

> **Buoyancy-constrained thermal width (2026-09-05):** a separate enthalpy
> width now matches the five fluxes AND the near-field buoyancy moment in all
> seven boundary cases. Maximum centre-temperature mismatch is 0.581 K;
> independent adaptive integration also passes. This is boundary-only:
> downstream width transport and experimental error reduction are not yet
> demonstrated. [Results](docs/buoyancy-constrained-enthalpy-width-results.md).

> **Enthalpy-profile screen (2026-09-05):** an opt-in Gaussian volumetric
> enthalpy profile now has its own conservative inverse and actual-density
> buoyancy integration. All seven cases and 38 arcs / 17 profiles completed.
> On the same 41 thermocouples, minimum-temperature MAE falls about 5%, but
> median-temperature error and centre-height MAE worsen. Do not promote;
> defaults remain unchanged. [Results](docs/gaussian-enthalpy-profile-results.md).

> **Ambient consistency audit (2026-09-05):** an opt-in research correction
> removes a spurious cold/enthalpy offset from hydrogen-free ambient gas.
> The full 38-arc/17-profile comparison shows only small concentration/height
> changes and does not pass the joint promotion gate. No default was changed.
> See [ambient and downstream profile results](docs/phase-ambient-consistency-results.md).

> **Latest research correction (2026-09-05):** the measured-LH2 HEM branch now
> carries pipe kinetic energy and condensed-phase source enthalpy into the
> conservative near field. This removes an independently detected ~6% energy
> transfer defect; it is not a claim of ~6% better experimental accuracy.
> The default DEGADIS path is unchanged. See
> [EOS/source-ledger results](docs/preslhy-source-eos-ledger-results.md).

> **Liquid-hydrogen model update (2026-09-03):** the flashing-source mass,
> enthalpy and total momentum are now conserved into the jet ODE, and the
> public pressure-jet path is allowed to detach from the ground. Current
> validation numbers and the remaining condensed-air limitation are in
> [`docs/lh2-model-improvements-2026-09-03.md`](docs/lh2-model-improvements-2026-09-03.md).
> The current crosswind research extension now has a 7/7 conservative
> independent-energy boundary and downstream ODE result. It improves
> concentration variance, FAC2 and vertical width, but worsens mean bias and
> centre-height error, so it remains research-only. See
> [`docs/preslhy-independent-energy-interface-results.md`](docs/preslhy-independent-energy-interface-results.md).
> A common sensor-height ground-image fit confirms that the geometry rejection
> is not an observation-operator artefact. A second published energy equation,
> Li et al. (2026)'s enthalpy-only balance, changes centre MAE by less than
> 0.00004 m and therefore rules out mean-kinetic-energy thermalisation as the
> missing rise physics. Direct section-integrated buoyancy also agrees across
> all seven 10D interfaces within 6.761%, so the remaining diagnosis is the
> downstream accumulation of buoyancy rather than an interface-force jump.
> A subsequent HyRAM profile audit corrects the Houf velocity-width mapping;
> MG/VG improve to 1.062/1.163 and internal width ratio to 0.976, but centre
> MAE remains worse than baseline, so the corrected model is still research-only.
> Completely suppressing dry-air condensation was then tested as a
> coefficient-free non-equilibrium bound. It barely changes the 10D thermal
> state and fails the unchanged interface-temperature gate, so no field score
> was issued and pre-10D mass entrainment is now the isolated target.
> That source audit found a second application of established-flow momentum
> entrainment inside the 6.2D Gaussian-development zone. Removing it as a
> four-flux lower bound improves concentration to MG 1.012, but geometry still
> misses the all-metric gate, so it remains a research option.
> A following ground-contact bound lowers trial 10's 6 m centre by only 0.042 m
> and is not promoted. Applying Li equation 35 before 10D also fails the
> conservative Gaussian boundary before field scoring. The remaining rise
> residual therefore requires an additional distributed/turbulent-energy or
> phase-slip state; it cannot be repaired by a ground or Zone-IV multiplier.

A modern Python reimplementation of **DEGADIS 2.1**, the dense gas dispersion
model of Spicer and Havens (University of Arkansas), released by the US EPA in
1989 and still named in 49 CFR 193.2059 as an acceptable means of determining
LNG vapour-gas dispersion exclusion zones.

There is no other open Python implementation. A GitHub-wide search returns one
repository (`diegomat87/DegadisMIO`), a C#/WPF front end that shells out to the
original DOS executables; the dispersion physics is not ported. `degadis`,
`pydegadis` and `degali` are all unclaimed on PyPI.

## Status

The model is complete. All five EPA test cases run end to end in Python — a
steady LNG pool spill, the same spill as a transient, a buoyant jet that never
lands, and two pressurised jets that touch down and continue as ground-level
plumes. The suite checks the port against the original Fortran, built from
source and running in the same repository.

```python
from degali import run_steady, run_transient, run_jet_to_ground

profile, source = run_steady("B9.INP")
profile.distance_to(0.05)      # 497 m to the lower flammable limit
profile.mass_above_lfl         # 6732 kg

out = run_transient("B9T.INP")
out.snapshots[7].mass_above_llc          # the cloud at t = 79 s
out.dose([Receptor(x=400.0)])[0].peak    # what a receptor at 400 m sees

profile, jet, src = run_jet_to_ground("EX2.INO", "EX2.IN")
jet.distance                   # 562 m to touchdown
```

or from a shell:

```
degali steady    B9.INP
degali transient B9T.INP --snapshot 60 --snapshot 120
degali dose      B9T.INP --at 200 --at 400 --at 800
degali jet       EX2.INO --bridge EX2.IN
```

| Component | Fortran | Status |
|---|---|---|
| Reference oracle (build + run + probe) | all | done, 5/5 test cases reproduce |
| Constants | `DEGIN.PRM` | done |
| Boundary layer | `ATMDEF`, `PSIF`, `ALPH`, `RHOA` | done, parity |
| Numerics | `ZBRENT`, `RKGST`, `GAMMA`, `INCGAMMA`, `AFGEN` | done, parity |
| Thermodynamics | `TPROP`, `ENTHAL`, `SETDEN`, `ADIABAT`, `ADDHEAT` | done, parity |
| Entrainment & surface exchange | `RIPHIF`, `SURFAC` | done, parity |
| Input decks | `IO.FOR`, `DEGINP` | done, parity |
| Parameter files | `ESTRT1` (`.ER1`) | done, parity |
| Secondary source blanket | `DEG1`, `SRC1`, `SRC1O` | done, parity |
| Source vector thinning | `CRFG` | done, parity |
| Blanket-free source layer | `SZF` | done, parity |
| Blanket-free source driver | `NOBL` | provisional (no test case reaches it) |
| Handoff file | `TRANS`/`STRT2` (`.TR2`) | done, parity |
| Downwind dispersion | `PSS`, `SSG` | done, parity |
| Downwind driver | `DEG2S`, `PSSOUT`, `SSGOUT`, `SSOUT` | done; profile within 3 % |
| Observer kinematics | `UIT`, `XIT`, `T0OB`, `TUPF`, `TDNF` | done, parity |
| Observer layer growth | `OB`, `OBOUT` | done, parity |
| Transient supervisor | `SSSUP`, `DEG2` | done; seeding at parity |
| Cloud snapshots | `DEG3`, `GETTS`, `SORTS`, `SORTS1` | done; masses within 3.5 % |
| Receptor time history | `DEG4`, `GETTD`, `DOSOUT` | done |
| Jet/plume (Ooms) | `JETPLU`, `SETJET`, `ELLIPS` | done, parity |
| Jet decks | `JETPLUIN` (`.INO`), `.IND` | done |
| Jet-to-ground bridge | `DEGBRIDG` | done, parity |
| One-call API and CLI | — | done |
| Real-fluid properties | — | CoolProp backend, all stages |
| Transient observers | `DEG2`, `OB`, `UIT`, `TUPF`, `DEG3` | later |

## How correctness is established

Printed DEGADIS output carries five or six significant figures, which is not
enough to distinguish a correct port from one that is merely close. So the
original Fortran is built and run inside this repository, and a purpose-written
probe dumps its internal state at full double precision. Every Python module is
compared against that.

Two layers of checking:

1. **The oracle is validated first.** The patched Fortran must still reproduce
   EPA's 2012 golden listings. It does: `B9`, `EX1`, `EX2` and `EX3` match line
   for line apart from timestamps, and `B9T` differs on 8 of 1086 lines in the
   seventh significant figure, which the EPA readme calls out as expected
   between compilers.
2. **The Python is compared to the oracle**, not to the listings, at a relative
   tolerance of 1e-12.

```
pytest -m "not slow"         # normal development regression
pytest -m slow               # whole-model/campaign/research regression
```

The Fortran oracle is built on the first compatible Linux run. On Windows,
the pure-Python and stored-reference checks still run, while the complete
gfortran oracle regression remains a Linux release check.

The field-trial tests locate their data from the environment and skip without
it, because none of it is ours to redistribute:

```
export REDIPHEM_ROOT=...        # the Risoe database
export DEGALI_E35_ROOT=...    # PRESLHY E3.5 workbooks, DOI 10.35097/1481
export DEGALI_E35_REPORT=...  # PRESLHY D3.6, for the sensor positions
export DEGALI_SMEDIS_ROOT=... # the SMEDIS spreadsheets
```

The count is deliberately not quoted here. It was quoted in four places, drifted
in all four, and `test_no_retired_figure_survives_in_the_prose` now fails the
suite if a superseded figure reappears in the prose.

## Burro 9, end to end

The steady-state LNG spill test case runs from the input deck to the
concentration profile in Python alone:

| quantity | DEGADIS 2.1 | degali |
|---|---|---|
| profile rows recorded | 35 | 35 |
| mole fraction, interpolated to the reference distances | — | within 1.1 % |
| concentration | — | within 2.5 % |
| sigma_z | — | within 0.5 % |
| distance to 5 mol % (the LFL) | 474–514 m (bracket) | 497 m |
| mass above the LFL | 6740.2 kg | 6731.7 kg (0.13 %) |
| mass between UFL and LFL | 3521.4 kg | 3509.7 kg (0.33 %) |

The residual comes from one place. `ADDHEAT` resolves the heated temperature
with `ZBRENT` at a 1e-3 K tolerance and then forms the heat capacity as
`dh/(temp − amt)`; early in the dense phase that denominator is only a few
hundredths of a kelvin, so the quotient carries percent-level uncertainty,
which feeds the ground heat flux and hence `dh` itself. At the first output
point the distance, the bisection count and the other three integrated states
all agree to 1e-6 while the added-heat state differs by 1.6 %, which is what
identifies it. Consequence: output points land at slightly different
distances, so profiles are compared by interpolation rather than row by row.

## Burro 9 transient, end to end

The same spill declared unsteady runs from the input deck to cloud snapshots.
DEGADIS handles this by releasing 30 imaginary observers over the source at
staggered times and following each downwind; a snapshot is then assembled by
finding where every observer is at that instant.

| snapshot | points (DEGADIS / degali) | cloud extent | mass above the LFL |
|---|---|---|---|
| t = 13 s | 5 / 5 | 78.0 / 78.0 m | 1031 / 1021 kg |
| t = 46 s | 16 / 16 | 252 / 252 m | 4531 / 4580 kg |
| t = 79 s (peak) | 28 / 28 | 439 / 439 m | 6516 / 6612 kg |
| t = 123 s | 30 / 30 | 699 / 700 m | 1964 / 2056 kg |
| t = 211 s | 30 / 30 | 1247 / 1247 m | 0 / 0 kg |

All 19 snapshots have exactly the point count the Fortran found, which means
the same observers were in range at the same instants — a strict check on the
release times, both edge crossings, the downwind integration and the
arrival-time mapping together. Extent agrees to 0.3 %; the flammable mass to
3.5 % near its peak.

## Pressurised releases

`JETPLU` is an integral jet/plume model in the style of Ooms: four coupled
balances — contaminant mass, total mass, and momentum in x and z — solved as a
4×4 linear system at every step, with an elliptical cross-section whose aspect
ratio follows the ambient dispersion and a wind speed averaged over it.

All three EPA jet cases reproduce against the full-precision `.VEL` trajectory
dump:

| case | release | trajectory | touchdown |
|---|---|---|---|
| EX1 | MIC, buoyant, never lands | 5.4e-6 over 4.7 km | correctly reports no landing |
| EX2 | ammonia pipeline, 56 kg/s | 4e-4 downstream of the stiff first 10 m | 562.463 m vs 562.479 m |
| EX3 | as EX2, simplified density | 4e-4 over 3 km | 3005.489 m vs 3005.473 m |

## Against measurement

Reproducing DEGADIS says nothing about whether DEGADIS is right.
`degali.validation` reads the Risø **REDIPHEM** database of heavy-gas
release trials, turns a trial into a deck while recording every assumption it
had to make, and reduces the comparison to the measures the published
evaluations report — with bootstrapped confidence intervals, because at these
sample sizes a point estimate hides whether a criterion is met.

`MG` is observed over predicted, so below 1 means the model reads high.

**At the lowest instrumented height**, which is what the literature reports:

| series | n | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|---|
| Burro, 8 LNG spills | 61 | 0.811 | [0.63, 1.02] | 2.53 | 0.56 |
| Desert Tortoise, 4 NH₃ jets | 6 | 1.839 | [1.63, 2.04] | 1.48 | 0.83 |

**Over the whole cloud**, using every instrumented height:

| series | n | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|---|
| Burro | 173 | 25.6 | [11.2, 65.0] | 4e19 | 0.29 |
| Desert Tortoise | 17 | 1.97 | [1.64, 2.38] | 1.85 | 0.53 |

The gap between the two is the result. **DEGADIS gets the ground-level
centreline roughly right and the vertical distribution badly wrong**: at 50 to
140 m downwind the measured concentration at 3 m is 0.6 to 1.0 of the value at
1 m, the model gives 0.02 to 0.18, and at 8 m the model is four orders of
magnitude down. The lateral spread runs the other way, one to three times too
wide against Desert Tortoise. The conventional single-height statistic is
carried by the one elevation where those two errors cancel.

Burro's conventional interval straddles Hanna's 0.7 bound, so even that is not
a pass. Anything depending on the cloud's shape rather than its centreline
concentration — a flammable volume, an exposure at head height, a cloud over a
bund wall — inherits both errors.

See [docs/field-validation.md](docs/field-validation.md) for the method, what
had to be withdrawn, and why the way measurements are reduced matters more
here than most modelling choices do.

## Liquid hydrogen

```python
from degali.lh2 import assess

r = assess(rate=9.2, pool_diameter=9.1, wind=1.6, at_distance=30.0)
print(r.report())
```
```
Liquid hydrogen release
  cloud regime                : aloft   [4 of 4 on the NASA trials, grade C]
  lowest flammable gas        : 15.0 m   [RMS 2.9 m, n=4, grade C]
  distance to 4 mol % (LFL)   : 16.3 m   [no direct pool-concentration validation]
  distance to stoichiometric  : 4.2 m
  buoyant below               : 99.98 mol %   [mixing line, grade A]
```

Each number is quoted with how far it has been checked, because they are not
equally well established, and `warnings` says so when the release is outside
the range for its release type. Corrected jets were checked over 0.084–0.285
kg/s, 0.6–4.2 m/s wind and 0.79–6 m; the NASA pool evidence is one 9.1 m pond
at 9.2–10.3 kg/s, 1.55–6.3 m/s wind and a 33.8 m tower row.


A dense-gas model is being extended to a substance whose equilibrium
flammable mixture is entirely buoyant. The corrected cryogenic source can
nevertheless remain denser than air for several metres while it warms: in
Spadeadam test 6 it crosses near 3.2 m. At the lower flammable limit of 4 mol
%, at stoichiometric, and even at the *upper* limit of 75 mol %, the warmed
mixture is lighter than air.

Against the PRESLHY/HSL near-field array the fully conserved jet path passes
all three Hanna screening criteria for momentum-driven releases — **MG 1.047,
95 % CI [0.759, 1.402], VG 1.425, FAC2 0.84** over 62 arc maxima downstream of
the established source plane, from nine trials and three nozzle sizes. Seven
0.35/0.53 m arcs are upstream of that plane and are not extrapolated backwards.
The interval is bootstrapped over *trials* rather than points because sensors
within one trial share a release, wind and source estimate.

On the same momentum-driven population, the modelled/measured vertical
Gaussian standard-deviation ratio is **1.033** over 23 well-constrained fits
(median 1.042). An earlier 0.731 value mixed JETPLU's standard deviation with
the measurement fit's e-folding width; the fit, reduced data and validation
now use one explicit `exp(-0.5*((z-zc)/sigma_z)^2)` convention. See
[docs/gaussian-width-convention.md](docs/gaussian-width-convention.md).

The quiescent near-nozzle regime now has its own conserved-energy
axisymmetric model, checked against all nine Hecht--Panda simultaneous Raman
concentration/temperature releases. Relative errors in the four published
centreline/half-width slopes are -8.1%, -18.3%, -19.9% and +14.7%; all meet
the pre-registered 25% band, and the aggregate metrics reproduce a separate
official HyRAM 6.1 oracle within 0.75%. An end-to-end audit nevertheless found
the reproduced plug-to-Gaussian boundary loses 8.9--14.8% species flux. Two
first conservative alternatives missed the thermal centreline criterion. A
later component-enthalpy phase model now conserves every boundary invariant.
Against the final-journal fits its corrected-coverage errors are -24.86%,
-15.79%, -21.38% and +24.56%, provisionally 4/4 but with two metrics very near
the 25% boundary. It is the recommended dry-air axisymmetric research path,
but is not substituted into the atmospheric calculation. See
[docs/prereg-hecht-panda-hyram-closure.md](docs/prereg-hecht-panda-hyram-closure.md)
and the later
[journal benchmark correction](docs/hecht-panda-journal-benchmark-correction.md).

The qualification matters: the 2019 article revised the two mass-fit labels
from the 2017 conference version, and its aggregate plots list an otherwise
untabulated `4 bar, 45 K, 1.25 mm` series. Until the fit membership and
uncertainty are recovered, these pass counts are aggregate-fit checks rather
than raw-data validation.

The experimental case weighting has since been corrected from a uniform
549-point camera range to the 369 points implied by each Table-1 stitched-image
count. A new conservative H2/N2/O2/H2O equilibrium option includes ambient
water freezing and latent heat. At 100% RH it overpredicts the Raman
temperature-decay slope by 44.8%; at an unscored 40%-RH sensitivity all four
slopes lie within 25%. Because the experiment reports condensed moisture but
not laboratory RH, humidity is now an explicit input and no fitted default is
claimed. See [docs/prereg-humid-air-frost-upper-bound.md](docs/prereg-humid-air-frost-upper-bound.md).
The accepted dry-air thermodynamic correction is recorded separately in
[docs/prereg-phase-temperature-dependent-enthalpy.md](docs/prereg-phase-temperature-dependent-enthalpy.md).

The accepted configuration is now available without assembling its internal
flags by hand. The boundary below is the measured choked-gas plane; the helper
conservatively expands it to ambient pressure before running the axisymmetric
near field.

```python
from degali import (
    lh2_source_from_measured_throat,
    run_lh2_near_field_research,
)

source = lh2_source_from_measured_throat(
    throat_diameter=0.001,       # m
    throat_pressure=2.42e5,      # absolute Pa
    throat_temperature=37.4,     # K
    throat_density=1.65,         # kg/m3
    throat_velocity=498.2,       # m/s
)
near = run_lh2_near_field_research(source)
print(near.report())
```

For a horizontal release aligned with the wind, the same conserved near field
can now be handed to JETPLU through an audited phase-aware boundary:

```python
from degali import (
    lh2_source_from_measured_throat,
    run_lh2_crosswind_research,
)

horizontal_source = lh2_source_from_measured_throat(
    throat_diameter=0.001,
    throat_pressure=2.42e5,
    throat_temperature=37.4,
    throat_density=1.65,
    throat_velocity=498.2,
    theta=0.0,
    y=0.5,
)

coupled = run_lh2_crosswind_research(
    horizontal_source,
    wind=2.5,                       # m/s at wind_reference_height
    height=0.5,
    wind_reference_height=10.0,
)
print(coupled.report())
if coupled.accepted:
    trajectory = coupled.run(distmx=100.0)
```

The handoff preserves total mass, hydrogen and both momentum components,
transfers the accepted N2/O2/H2O phase/enthalpy profile, and separately screens
energy, H2 width and centre temperature. A failed screen blocks the downwind
run. In the strict-grid representative test, the accepted 0.08 m handoff
errors are 0.0585% energy, 1.57% H2 half-width and 0.94 K centre temperature; see
[the handoff results](docs/conservative-nearfield-crosswind-handoff-results.md).
Independent PRESLHY testing now covers seven momentum-dominated horizontal
releases. The initial conserved coupling improves concentration variance but
underpredicts vertical width; switching at the handoff to JETPLU's local
density-scaled shear entrainment repairs the mean width ratio from 0.710 to
0.994. That mechanism is supported, but a fully four-flux source boundary
still fails the unchanged temperature/profile compatibility screen in two
large-source trials. The coupled path therefore remains research-only and has
not replaced `assess()`; see
[the coupled PRESLHY results](docs/preslhy-coupled-crosswind-results.md).
The coupled runner's default near-field endpoint is the strictly confirmed
0.08 m boundary; longer or large-source domains require an explicit handoff
distance and a passing audit.

The function defaults to dry air, zero radiative absorption and the strict
81-point validation settings. Humidity, co-flow, argon phase change and grey
radiative absorptivity are explicit sensitivity inputs; enabling one adds an
applicability warning rather than borrowing the dry-air 4/4 claim. The legacy
DEGADIS and atmospheric `assess()` paths are unchanged.

Two additional prospective screens do not improve the default. Atmospheric
argon preserves a provisional 4/4 and all conservation checks but slightly worsens the
centreline-temperature error. A perfectly black `5B` plume envelope improves
that thermal error by about 0.83 percentage point but slightly worsens both
mass metrics; real absorption would be smaller. See
[docs/prereg-argon-phase-completeness.md](docs/prereg-argon-phase-completeness.md)
and [docs/prereg-radiation-upper-bound.md](docs/prereg-radiation-upper-bound.md).
The hydrogen spin-isomer caloric term is now explicit as well. A frozen
`ParaHydrogen` run improved the two mass slopes but worsened the primary
thermal-centreline error from -21.38% to -34.79%, so normal-hydrogen calorics
remain the recommended baseline; see
[docs/prereg-hydrogen-spin-isomer-enthalpy.md](docs/prereg-hydrogen-spin-isomer-enthalpy.md).
An independent thermal/species Gaussian profile with fixed published
`Pr_t/Sc_t` was then combined with the accepted phase calorics and a four-flux
boundary. It repairs both widths but worsens centreline temperature to
-29.14%, so it remains a research closure rather than the default; see
[docs/prereg-phase-two-scalar-four-flux.md](docs/prereg-phase-two-scalar-four-flux.md).

ELVHYS Tests 10/11 were screened from verified public raw files. They are not
yet scoreable: no H2 mass-flow channel is present and the official sensor
sheet and final D4.6 report disagree on nozzle elevation. See
[docs/elvhys-tcs-audit.md](docs/elvhys-tcs-audit.md).

Wind-steered releases, where the exit velocity is comparable with the wind, are
outside what a steady jet model can describe and the code says so. Getting
there meant discarding two plausible fixes that made things worse and one
comparison that was measuring the wrong quantity. See
[docs/lh2-model-improvements-2026-09-03.md](docs/lh2-model-improvements-2026-09-03.md).

## Documentation

| | |
|---|---|
| [docs/status.md](docs/status.md) | where the project stands |
| [docs/validation.md](docs/validation.md) | against the original Fortran |
| [docs/field-validation.md](docs/field-validation.md) | against LNG and ammonia field trials |
| [docs/data-inventory.md](docs/data-inventory.md) | every dataset, its format, and its traps |
| [docs/data-exhaustion.md](docs/data-exhaustion.md) | what each dataset was asked, and what it refused |
| [docs/lh2-plan.md](docs/lh2-plan.md) | liquid hydrogen: what can and cannot be attempted |
| **[docs/DATA_AND_REPRODUCTION.md](docs/DATA_AND_REPRODUCTION.md)** | **start here: the data, how to reproduce every number, and a program audit** |
| **[docs/lh2-model-improvements-2026-09-03.md](docs/lh2-model-improvements-2026-09-03.md)** | **liquid hydrogen: current conserved-source model, validation and remaining physics gaps** |
| [docs/lh2-recomputed.md](docs/lh2-recomputed.md) | historical liquid-hydrogen reconstruction, superseded by the update above |
| [docs/lh2-results.md](docs/lh2-results.md) | liquid hydrogen: what has been established |
| [docs/lh2-datasets.md](docs/lh2-datasets.md) | the other LH₂ datasets, and why one of them is better |
| [docs/elvhys-tcs-audit.md](docs/elvhys-tcs-audit.md) | verified ELVHYS Test-10/Test-11 reduction and why scoring stops |
| [docs/data-request-elvhys-source.md](docs/data-request-elvhys-source.md) | exact ELVHYS source, geometry and instrument fields to request |
| [docs/data-request-hecht-panda-humidity.md](docs/data-request-hecht-panda-humidity.md) | exact missing Raman humidity/raw-data fields and ready-to-send request |
| [docs/prereg-hydrogen-spin-isomer-enthalpy.md](docs/prereg-hydrogen-spin-isomer-enthalpy.md) | rejected para-hydrogen caloric sensitivity |
| [docs/prereg-preslhy-liquid-spin-source.md](docs/prereg-preslhy-liquid-spin-source.md) | spin-consistent measured-LH2 source and seven-trial para-hydrogen bound |
| [docs/prereg-phase-two-scalar-four-flux.md](docs/prereg-phase-two-scalar-four-flux.md) | conservative radial heat/species closure and its rejection |
| [docs/preslhy-coupled-crosswind-results.md](docs/preslhy-coupled-crosswind-results.md) | coupled LH2 field validation, accepted entrainment mechanism and remaining thermal-state gap |
| [docs/liftoff.md](docs/liftoff.md) | buoyant lift-off: the closure, and its first validation |
| [docs/integral-limits.md](docs/integral-limits.md) | the structural limits of an integral model for LH₂ |
| **[docs/HANDOVER.md](docs/HANDOVER.md)** | **state, method, data, mistakes, and what to do next — start here** |
| [docs/architecture.md](docs/architecture.md) | what is kept in Fortran shape, and why |
| [docs/release-readiness.md](docs/release-readiness.md) | the checks applied before publishing |
| [docs/references.md](docs/references.md) | every source, and whether it was read or quoted |
| [docs/claim-grading.md](docs/claim-grading.md) | which results are deterministic, which statistical, which only directional |
| [docs/prereg-entrainment.md](docs/prereg-entrainment.md) | a pre-registered test, and its falsification |
| [docs/prereg-trajectory.md](docs/prereg-trajectory.md) | the trajectory fault, isolated to the buoyancy balance |
| [docs/prereg-boussinesq.md](docs/prereg-boussinesq.md) | the correction the literature points at, and why it fails |
| [docs/prereg-expanded-source.md](docs/prereg-expanded-source.md) | historical three-correction stage and four rejected alternatives |
| [docs/local-validation.md](docs/local-validation.md) | measuring the sub-models separately: trajectory, spread, width |
| [docs/hydrogen-audit.md](docs/hydrogen-audit.md) | what liquid hydrogen breaks in a model built for LNG |

## Two backends

Every module exposes a `legacy` and a modern path, and the choice is one
argument:

```python
run_steady("B9.INP", backend="coolprop")
```
```
degali steady B9.INP --backend coolprop
```

**`legacy`** (the default) is bit-faithful to 1989, including approximations
that are not round-off:

- `SETDEN` builds the adiabatic mixing table on a grid computed with the
  Fortran `FLOAT()` intrinsic, i.e. in *single* precision. This puts a
  systematic 1e-7 relative error into every table node. Reproducing it is
  required to reproduce anything downstream.
- `ALPH` fits the wind-profile exponent by driving `ZBRENT` around an `RKGST`
  quadrature with `ERBNDZ = 0.005`. The quadrature error moves the fitted
  exponent by 2e-5 relative — more than the root-finder tolerance. `RKGST` is
  therefore ported exactly rather than replaced.
- `GSERIES` clamps its argument at 0.999 and stops at 7e-5 relative, which
  shifts `PHIHAT` by up to 2e-4.
- `ZBRENT` cannot converge below its hard-wired `EPS = 3e-8`.
- `GAMINC`, `SRC1O` and `CRFG` all thin their output against fixed relative
  criteria, so which points survive is part of the answer.

**`coolprop`** replaces the correlations with equations of state *and* fixes
the one place the original's numerics lose accuracy. It resolves the deck's
three-character label and molecular weight to a CoolProp fluid — `LNG` and
16.04 give methane — so the contaminant heat capacity and density come from an
EOS rather than the deck's two fitted constants. When nothing matches, water
and air still come from the EOS and the contaminant keeps its correlations.

Selecting it also switches `legacy_numerics` off: the enthalpy-to-temperature
inversion is solved to 1e-12 K instead of 1e-3 K, and heat capacities become
secants over accurately located intervals rather than quotients of two small
differences. That removes the residual documented below. The two choices are
independent — `make_thermo(legacy_numerics=False)` keeps the 1989
correlations with accurate numerics.

### What it changes

| Burro 9, steady | legacy | coolprop | change |
|---|---|---|---|
| source enthalpy | −409 207 J/kg | −434 442 J/kg | **−6.2 %** |
| distance to the UFL | 245.9 m | 243.8 m | −0.9 % |
| distance to the LFL | 497.1 m | 494.5 m | −0.5 % |
| mass above the LFL | 6732 kg | 6640 kg | −1.4 % |
| mass between UFL and LFL | 3510 kg | 3568 kg | +1.7 % |

| Burro 9, transient | legacy | coolprop | change |
|---|---|---|---|
| peak flammable mass | 6619 kg at 79 s | 6393 kg at 79 s | −3.4 % |
| peak at 400 m | 0.0725 | 0.0743 | +2.5 % |

The source enthalpy moves 6 %, the buoyancy excess 7 %, and the final
distances only half a per cent — the density change feeds gravity slumping as
a square root and entrainment suppression through Φ, and the two largely
cancel. Whether that cancellation holds against measurements is exactly the
kind of question the two backends exist to ask.

Jet cases are unaffected, and that is correct rather than a bug: a `.INO` deck
sets `ISOFL = 1` and supplies its own density table, so no mixing line is
built and there is nothing for the backend to change.

CoolProp costs about fifteen times the legacy runtime. Each property is
tabulated once on a 4001-point grid spanning 80–400 K and interpolated
thereafter, which keeps it a smooth function of temperature — the integrators
effectively differentiate these properties, so a rounded cache would put small
steps into them and the adaptive step control would chase them.

## Modernisation

| Original | Replacement |
|---|---|
| `RKGST` (Runge-Kutta-Gill, manual halving) | `scipy.integrate.solve_ivp` |
| `ZBRENT`, `LIMIT` | `scipy.optimize.brentq` |
| `GAMMA`, `INCGAMMA`, `ERF`, `GSERIES` | `scipy.special` |
| `AFGEN`, `AFGEN2` | `numpy.interp` |
| `WATVP` (Antoine-like fit) | CoolProp saturation line |
| `CPC`, `GASRHO`, `ADIABAT` | CoolProp `HEOS` |
| `SIMUL` (Gauss elimination) | `numpy.linalg.solve` |
| 27 `COMMON` blocks | dataclasses |
| Fortran unit 8/9/13 scratch files | in-memory results |

## Portability patches to the Fortran

Eight changes were needed to build the 2012 sources with gfortran on Linux.
All are portability fixes; none alters a computed value. They are enumerated in
`degali.validation.reference.PORTABILITY_PATCHES`. One is worth calling out:

`COMMON /ERROR/` interleaves eighteen `REAL*8`, one `INTEGER*4` (`NOBLpt`) and
two more `REAL*8`. `ESTRT1` overlays its read buffers onto that block with
`EQUIVALENCE`, which makes gfortran pad the block differently in that one
compilation unit than in every other. The result was that `CRFGER` and
`EPSILON` were silently read as zero — **disabling air entrainment into the
source blanket entirely**, with no diagnostic. Intel's `/align:dcommons` had
been masking it. Replacing the overlay with explicit assignment fixed it, and
`B9` matched the golden listing immediately.

## Three things that trip up a port

Confirmed against the Fortran by line number, and reproduced rather than
fixed. They are not all the same kind of finding, and calling them all "bugs"
would overstate two of them.

**`ADIABAT` does not compute `wa` on the `ifl=1` path.** `TPROP.for` lines
300–320 assign `wa` for `ifl` of 0, −1, −2 and 2; there is no such assignment
for 1, so on that path the routine uses whatever the caller passed, to form
the molecular weight and hence pick which panel of the mixing table to
interpolate on. `SZF.for` line 95 calls it that way and passes `walay` — a
name that appears nowhere else in that file, with no declaration and no
assignment. Under `/noauto` it is static, zero and never written, so `SZF`
runs its lookups as though the mixture contained no dry air, shifting the
layer density by up to 0.5 %. **This one is a defect.**

**`PSS` and `SSG` differ by one letter.** `ADDHEAT` writes `rholay`, `temlay`
and `cp`; the `ADIABAT` call after it puts its density in a different
variable (`rholam`), so the heated density survives while the temperature is
overwritten. Then `PSS.for` line 89 writes that temperature into `temlay` and
`SSG.for` line 92 into `temlam` — the two routines are the steady and
transient counterparts of each other and every other argument matches. So the
two stages compute the ground heat flux, and apply the cut-off that zeroes
it, from different temperatures. **This one is a defect too**, and the line
numbers pin it.

**`GAMINC` returns the unregularised incomplete gamma, deliberately.** It
computes the Numerical Recipes `gammp` and then undoes the normalisation, and
the comment above that line says why: "multiply the result by GAMMA(ALPHA) to
get the final value". In `PSS` the result is compared against a cap built from
Γ(1/(1+α)), which is consistent. **This is not a defect** — it is a trap for
anyone mapping it onto `scipy.special.gammainc`, which is regularised. Doing
so leaves the flammable-mass derivative low by exactly that factor while every
other quantity in the routine stays correct.

The distinction matters. A parallel reimplementation of SLAB reported six
defects in the original and, on obtaining Ermak's Fortran, found that only one
of them was: three were artefacts of a JavaScript transcription that had been
treated as the original, and two were deliberate sentinels misread as errors.
Attribution needs the source and a line number, not a plausible reading.

## Credits

DEGADIS was developed by Thomas O. Spicer and Jerry A. Havens, University of
Arkansas, under US Coast Guard contract DT-CG-23-80-C-20029, with subsequent
support from the Gas Research Institute, the American Petroleum Institute and
the US EPA. The original source is in the public domain via EPA SCRAM.
