# Data, reproduction, and a program audit

> **LH2 update (2026-09-03):** current source-state invariants, 62-arc
> corrected statistics, and Spadeadam detachment results are recorded in
> [`lh2-model-improvements-2026-09-03.md`](lh2-model-improvements-2026-09-03.md).

Everything needed to reproduce every number in this project from a clean
machine, plus an audit of what the code can be told to be and what is wrong
with it.

Written so that a session with no history can start here. Read §2 if you want
to run something today; read §1 if you want to know what the numbers rest on.

The suite is split by cost rather than documented with a count that will
drift. `pytest -m "not slow"` is the normal development check;
`pytest -m slow` runs whole-model, campaign and strict multiphase cases.

---

# Part 1 — Data

## 1.1 What ships with the repository

Four files, 512 kB, sufficient to recompute every statistic in the papers
without any archive.

| file | size | what |
|---|---|---|
| `reference/preslhy/e35_reduced.json` | 108 kB | the LH₂ campaign, reduced |
| `reference/rediphem_reduce_full.json` | 348 kB | the heavy-gas trials, reduced |
| `reference/preslhy/conditions.csv` | 8 kB | per-trial conditions, digitised |
| `reference/preslhy/farfield.csv` | 48 kB | far-field readings, 686 rows |

These are **derived products of this work**, not the datasets. Neither archive
is ours to redistribute. Everything that made them is in the repository, so
they can be regenerated rather than trusted.

### `preslhy/e35_reduced.json` — 24 trials

```
trials[24]
  trial, orientation, release_height_m, orifice_mm, tanker_barg
  flow_mean_gs, flow_peak_gs, window
  wind_ms, wind_ref_m, T_C, RH_pct
  sensors[520]        serial, x, y, z, z_axis, peak, mean
  vertical_fits[72]   x, centre, sigma_z, peak, r2, points, well_constrained
arc_fits[65]          trial, radius, height, centreline, offset, sigma_y,
                      r2, points, well_constrained
gaussian_width_definition
```

**`z` is above ground; `z_axis` is above the release axis.** Both are kept
deliberately — see §1.5.

**`flow_mean_gs` is the one to use.** It is the mean over the longest sustained
run above half the smoothed peak. `flow_peak_gs` is the peak of the same
window. On trial 10 they are 189.4 and 285.3 g/s, a factor of 1.5, and using
the wrong one moves MG by 0.2 and changes which trials pass the momentum
filter.

`sigma_y` and `sigma_z` are statistical standard deviations and the metadata
records the fitted profile explicitly:
`exp(-0.5*((coordinate-centre)/sigma)^2)`. Earlier reductions used the
e-folding width from `exp(-((coordinate-centre)/w)^2)` but named it sigma;
`w/sqrt(2)` is the standard deviation. The bundled reduction has been
migrated. `vertical()` converts legacy files with no convention metadata on
read.

`well_constrained` means `r² > 0.85`, at least four points, and the original
3 m e-folding-width rejection limit, now expressed as
`sigma_z < 3/sqrt(2) m`.

### `rediphem_reduce_full.json` — 313 trials, 221 usable

```
case         u0, z0, zr, rml, relhum, avtime, istab, stability, yclow,
             oodist, gmass0, tend, instantaneous, steady_state
source_term  rate, radius, temp, rho, wc, fracv, time, tend, enthalpy
gas          name, mw, temp, rho
ambient      tamb, pamb, humid, tsurf, ihtfl, iwtfl
arcs         maxima at z = 0.1, 1.0, 3.0, 8.0 m, averaging 18.4 s
case_error   present where `to_case` refused, with the reason
```

Series: HAMBURG 146, LATHEN 51, WSL 34, WSL_RED 34, FLADIS 14, TNO 13,
BURRO 8, EAGLE 4, TORTOISE 4, COYOTE 3, TI 2.

**It stores the built `Case`, not the specification.** That is the right
choice: `to_case` makes judgements — which source route a release type takes,
what to assume for a missing pool diameter — that should be made once and
recorded. Trials it refused carry `case_error`, so the skipped population
matches the one the original statistics skipped.

**Two fields are not in it and must not be left at their defaults:**

- **`ulc` and `llc`**, the levels of concern. With `llc` at its dataclass
  default the flammable-mass integral evaluates `log(cc / clow)` at 28 and
  `SERIES` raises its overflow guard on *every* trial.
  `cases_from_reduced` takes them from `trialcase.SUBSTANCES` and raises on an
  unknown substance rather than defaulting.
- **The dispersion coefficients**, which follow from roughness, stability
  class and averaging time and are recomputed. The Monin-Obukhov length *is*
  stored, because the archive measures it and it does not follow from the
  class.

### `preslhy/conditions.csv` and `farfield.csv`

Per-trial conditions as digitised from report D3.6, and 686 far-field readings
with `saturated` and `n_at_ceiling` flags.

**Ignore this file's flow columns.** `flow_meter_mean_gs` is a different
window, digitised from the report rather than computed, and `flow_used_gs`
takes `table4` nominal set-points for the 1 barg trials — 139.5 or 105.5 g/s
repeated across trials against measured peaks of 27–50. They were never used
and should be dropped from the file.

## 1.2 The source archives, and how to get them

None can be redistributed. Every reader takes an environment variable and the
tests skip cleanly without it.

### PRESLHY E3.5 — the primary LH₂ dataset

**DOI 10.35097/1481**, on RADAR (KIT). Resolve `https://doi.org/10.35097/1481`
in a browser. The page also holds ~10.5 GB of video — **take only the 24
`.xlsx`, about 19 MB**.

Twenty-four releases at HSL, September 2019. Five sheets per workbook; two
matter:

- **`Xensor`** — NREL thermal-conductivity sensors, **0 to 100 vol %**, about
  3 Hz. This is the reason to use the raw data: the far-field Dräger devices
  top out at 4 %, hydrogen's lower flammable limit, so every reading that
  matters there is censored.
- **`Flowmeter`** — Coriolis mass flow at 1 Hz, so the source is measured
  through the release rather than taken from a summary table.

Column headers carry the sensor serials
(`X2019_09_10_02EC42Output` → `02EC42`), matched by regex
`_(\d{2}[A-Z]{2}\d+)Output`. **The join to positions is the serial.** 25 of 31
columns have coordinates.

`PRESLHY_D3_6_Summary_of_Rainout_Experiments_V1_20.pdf` supplies Table A3
(pages 43–46, sensor positions) and Figure A4 (far-field stand layout). Public
deliverable, Fuel Cells and Hydrogen JU grant 779613.

Read with `openpyxl`, `read_only=True, data_only=True`.

### REDIPHEM — LNG and ammonia field trials

The JRC archive. `REDIPHEM_ROOT` points at the directory holding the series
folders.

**FLADIS ships no `CHANDEF.DAT`**, so its channel numbering cannot be
recovered. **Do not guess it.** Guessing selected channels reading 302 and
23.5 at 20 m downwind — not concentrations — and produced a clean-looking
MG 0.89 that meant nothing. Use SMEDIS instead.

### SMEDIS — the EU model-evaluation spreadsheets

30 trials, 1250 sensors, 79 arc reductions. Preferable to the raw archives:
every sensor is a row with its own position and value, **the wind direction is
recorded with its standard deviation**, and the arc reductions including
`sigma_y` are already computed.

The reader handles three format variants, including section headers where the
second column is empty — which is not the same as a column-header line, though
both start with `#`.

### NASA Witcofski

Witcofski and Chirivella (1984), *Int. J. Hydrogen Energy* **9**(5), 425–435.
Tables 1–4 used directly. Table 2 (time-resolved grab bottles, tower 5) is
scrambled by OCR; `pypdf` and OCR both interleave the columns. Someone with
the printed journal is the only route.

## 1.3 What each dataset settled

| dataset | settled |
|---|---|
| DEGADIS 2.1 Fortran + 5 EPA cases | the port, to 1e-12, six programs |
| REDIPHEM Burro (8) | MG 0.620 at 1 m, n=61 — published 0.811 |
| REDIPHEM Desert Tortoise (4) | MG 1.839; EPA's own score explained |
| SMEDIS FLADIS + DT (76 sensors) | the vertical defect on an independent substance |
| PRESLHY near field | historical MG 1.121 on all 69 arcs; fully conserved source MG 1.047 on 62 arcs downstream of its established plane |
| PRESLHY far field (15 points, 3–7 m) | corrected MG 0.904, VG 2.115; trial 20 remains the main outlier |
| PRESLHY vertical (23 filtered fits) | corrected model/measured `sigma_z` is **1.033**; mean/median signed centre error 0.036/0.002 m |
| NASA Witcofski (4) | buoyancy regime 4 of 4, RMS 2.9 m |
| EPA-450/4-90-018 | an external anchor and the height artefact |

## 1.4 What each dataset refused, and why

| dataset | refused because |
|---|---|
| PRESLHY far field, point-to-point | wind direction to 22.5° against stands 10–12° apart — solved by fitting the position out |
| PRESLHY near field, transients | travel time under a second against 0.3 s sampling |
| **PRESLHY Dräger sheet** | 90 channels at 1 Hz, values −275 to 1970 under a `%` header; **units unresolved, and this is the one place transients could be tested on hydrogen** |
| REDIPHEM FLADIS | no `CHANDEF.DAT`; recovered from SMEDIS |
| REDIPHEM Eagle | N₂O₄ dissociates to NO₂; molecular weight is temperature-dependent |
| REDIPHEM Thorney Island | the van Ulden momentum balance has no solution at H/D ≈ 1 — **the original Fortran stops in the same place**, and refining the grid a hundredfold does not move it |
| SMEDIS Thorney Island | peaks at 2060 under a `mean_C(%)` header |
| SMEDIS Prairie Grass | peaks at 235, same header |
| SMEDIS EMU | sensor heights read 244 553 m |
| BA-Propane (273 sensors) | MG 0.168 against the REDIPHEM reduction's 3.76 — **the same trials, two reductions, a factor of twenty apart** |
| PRESLHY E3.4 pool | concentrations at 35/45/55 cm only, no downwind distance |
| PRESLHY E3.1a | near-orifice discharge, 1.2 GB in nested archives |

## 1.5 Traps that have already cost a comparison

Each presented as a model defect and was a data-handling defect. All are now
covered by tests; they are listed because the next dataset will have its own.

**Table A3's `z` is above the release axis, not the ground.** The near-field
array was rebuilt for each release height — Figure 4 of D3.6 is captioned
"near-field array configured for releases at 1.5 m height". Reading it as a
ground height puts the 1.5 m releases' sensors a metre below their own axis,
where the Gaussian factor is `exp(-0.5*(1.0/0.15)²)` and the prediction is
about ten to the minus ten. Reported as the model failing catastrophically on ten
trials.

The data says the same independently: at 0.35 m downwind the 0.5 m sensor
reads 80.8 vol % for a 0.5 m release and 87.9 for a 1.5 m one, and a jet
cannot fall a metre in thirty-five centimetres.

**The flow window must be the longest sustained run.** Taking the outermost
crossings of half the peak gives, on trial 2, a window in which six per cent
of samples are flowing and a **mean mass flow of −0.5 g/s** — the meter
oscillates between −18.8 and +30.5 there. That number then drives the source.
Detect the longest run above 50 % of a 5-sample moving average.

**Over-range readings exist.** Four sensors exceed 100 vol %. Dropped.

**The far-field plume cannot be located from the wind record, and does not
need to be.** Stands sit on arcs at 10 and 14 m, 2 to 2.5 m apart across the
wind, so adjacent stands are 10–12° apart; the wind direction is a compass
point at 22.5°, once every five minutes. But five stands across an arc
constrain a Gaussian: fit it, and the centreline concentration and `sigma_y`
come out with the position as a nuisance parameter. 18 of 65 arcs are well
constrained.

**The source route must follow the release type.** Sending Thorney Island — a
14 m cylinder released at once — through the steady pool route produces a
source of 3970 kg/s over a 2050 m radius, which fails a degeneracy check and
was reported as "no usable measurements" while the trial carried eighteen
arcs.

**Thorney Island uses site grid coordinates** in SMEDIS, release at (400, 200),
so a sensor at x = 100 is 300 m *upwind*. The reader subtracts the origin.

**REDIPHEM arc maxima are keyed by raw sensor distance.** On Burro 9 the dict
holds 37, 40, 49, 55.6 and 57 m as separate entries; those are lateral
positions on two nominal arcs, and the off-centreline ones read 0.002 vol %
against neighbours reading 5. Dropping readings at or below 0.1 vol % is what
the published statistic did, and it is what takes VG from 1574 to 3.32.

## 1.6 The standing instruction

**When a validation result looks wrong, suspect the comparison before the
model.** Six times now a "model defect" has been a data-handling defect. The
seventh will be in whatever dataset arrives next.

---

# Part 2 — Reproduction

## 2.1 From nothing to a green suite

```bash
# 1. Toolchain. The reference implementation is built from source and run as
#    part of the suite: the port is checked against DEGADIS 2.1 itself, not
#    against recorded output.
apt-get install -y gfortran
pip install numpy scipy pytest CoolProp

# 2. The repository
tar xzf degali.tar.gz && cd degali
pip install -e ".[test]"

# 3. Build the oracle (also happens automatically on first test run)
cd reference/fortran && bash build.sh && cd ../..

# 4. Fast development regression
pytest -m "not slow" -q

# 5. Whole-model and high-cost research regression
pytest -m slow -q
```

The fast group should complete in well under the time required for the full
physics/campaign group. Field tests skip and say which variable would enable
them.

## 2.2 With the archives

```bash
export REDIPHEM_ROOT=/path/to/rediphem            # series folders
export DEGALI_E35_ROOT=/path/to/10.35097-1481/data/dataset
export DEGALI_E35_REPORT=/path/to/PRESLHY_D3_6_...pdf
export DEGALI_SMEDIS_ROOT=/path/to/smedis       # holding batch1,2,3
pytest -m "not slow" -q
pytest -m slow -q
```

## 2.3 Reproducing each headline number

Every one of these is computed by a test. None is transcribed.

```bash
# provenance: 24 workbooks -> 9 trials -> 69 arcs
pytest -k exclusion_chain_reproduces

# LH2 near-field concentration, both conventions
pytest -k "near_field_statistic_recomputed or convention_moves"

# the vertical sub-models
pytest -k "vertical_spread_correction_reproduces or trajectory"

# both LH2 configurations against their parameter dumps, 4 s.f.
pytest -k match_their_reference_dumps

# Burro
pytest -k burro_pairs

# the port itself
pytest -k "parity or epa_golden"
```

## 2.4 The LH₂ jet configuration

This did not exist in the package until recently; every LH₂ result had been
produced by assembling it by hand in a session, which is why none of them had
a test. It is now `validation.nearfield.hydrogen_jet`.

```python
from degali.validation.nearfield import hydrogen_jet, Trajectory, STEP, REACH

jp, y0 = hydrogen_jet(
    rate=0.189425, diameter=0.0254, wind=2.4667, height=0.5,
    ambient_temperature=288.58, relative_humidity=58.0,
    storage_pressure_barg=5.0, wind_reference_height=1.5,
    corrections=False,          # True switches all three on together
)
traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)
state = traj.at(6.05)           # centre height, sigma_z, sigma_y, centreline
value = traj.concentration_at(6.05, 0.0, 0.5)   # vol % at a sensor
```

### Five settings that are not obvious and are not defaults

| | value | why it matters |
|---|---|---|
| `distmx` | **0.2** | it is the integration **step**, not the limit |
| `smax` | 40.0 | the limit |
| averaging | **60 s** | not the 18.4 s used elsewhere |
| roughness | **0.001 m** | not 0.01 |
| `yclow` | **1e-5** | with 1e-4 the integration ends early and the far arcs drop |

**`distmx` is the one to remember.** Passing 40.0 there integrates the whole
plume in a single stride. The result is monotone, smooth, the right order of
magnitude, and wrong by a factor of two. Nothing in the output says it is
under-resolved.

### The two source planes

| | as shipped | corrected |
|---|---|---|
| `rho_exit` | 5.34689 — orifice, two-phase after flash | 2.54663 — expanded |
| `concentration` | 1.0 | 0.44141 |
| `diajet` | 0.02540 m — the orifice | 0.04602 m — equivalent bore |
| exit velocity | 69.9 m/s | — |
| `rhoe` | **1.33217 in both** — the saturated-vapour density, never the expanded one |

The baseline plane is worth stating because it is what the expanded-source
correction exists to fix: **the density is the two-phase value after the flash
while the concentration and the area are the orifice's.** One plane for one
quantity and a different plane for the other two. The correction makes the
three consistent, which is a coherence argument rather than a fitted
improvement.

## 2.5 The comparison conventions, which are decisions

Stated once here and applied everywhere, including to baselines.

**Concentration: at-sensor.** The prediction is evaluated at each sensor's own
position and maximised over the arc, the same operation the measurement gets.
By 6 m the model centreline sits above the topmost sensor, so it is a
concentration the array could not have measured.

For Burro this means DEGADIS's own vertical profile at the sensor height:

```python
pred = centre * math.exp(-((z / sigma_z) ** (alpha + 1))) * 100.0
```

At z = 1 m that is worth a factor of eight in the near field — far more than a
Gaussian, because the exponent exceeds one and `sigma_z` is small close in.
**The vertical profile being too steep, a defect documented separately, is what
makes the height matter so much.**

**The convention costs VG.** LH₂ near field goes from 1.44 to 19.1. That is the
tail sensitivity surfacing, not the model worsening: once the plume centre
climbs past the array the arc maximum is a tail value, and tails are where
small differences in `sigma_z` become large differences in concentration.

**The vertical fits use two different populations, and this was never written
down before.**

| statistic | basis | n | form |
|---|---|---|---|
| trajectory | well-constrained fits on horizontal trials | 42 | median |
| `sigma_z` ratio | the same, **momentum-filtered** | 23 | **mean of per-fit ratios** |

The momentum filter applies to concentration statistics and to the spread
ratio, and not to the trajectory fits, because a fitted plume centre is a
geometric measurement and the filter is a statement about where a steady jet
model applies. Mean of per-fit ratios and ratio of medians are different
statistics and are not substituted for one another.

The formerly reported 0.64 and 0.731 width ratios used unlike Gaussian width
definitions. Applying `sigma=w/sqrt(2)` to the measurement gives 0.901 and
1.033 on the same 23 fits. Nothing about the fitted curves or selected points
changes; the false factor-`sqrt(2)` deficit disappears.

**Bootstrap over trials, not points.** Readings within a trial share a
release, a wind and a source estimate. Point-level resampling gave
[0.68, 0.80] where trial-level gave [0.59, 0.92].

**State the aggregation with every statistic**: distance band, release
heights, filter, statistic form, n. Every discrepancy found in reconstructing
this work came from a number reported without its subset.

## 2.6 Rebuilding the reductions

```python
from degali.validation import rediphem
import json
json.dump(rediphem.reduce(), open("reference/rediphem_reduce_full.json", "w"))
```

```python
from degali.validation.nearfield import reduce
reduce(E35_ROOT, E35_REPORT, "reference/preslhy/conditions.csv",
       "reference/preslhy/nearfield.csv")
```

Deterministic given the same archive. **If a rebuild disagrees with the shipped
file, that is a finding** — record it rather than overwriting.

---

# Part 3 — Program audit

## 3.1 Can it still be the original DEGADIS?

**Yes, and this is now tested rather than asserted.**

`test_the_liquid_hydrogen_machinery_does_not_perturb_degadis` runs EPA's B9
deck, imports every LH₂ module, runs it again, and requires the arrays to be
identical. Both halves matter: an option defaulting to on would change the
ported model, and a module that patched a constant at import time would change
it only for programs that happened to import it — which is worse, because the
suite would still pass whenever it ran the two in isolation.

Every switch, and its default:

| switch | default | what it is |
|---|---|---|
| `JetCoefficients.alfa1` | **0.057** | DEGADIS 2.1's jet entrainment |
| `alfa2`, `sc`, `cd`, `delta` | 0.5, 1.42, 0.2, 2.15 | unchanged |
| `plume_transition` | False | implemented; measured pure-plume threshold is 0.716 |
| `density_scaled_entrainment` | False | Ricou–Spalding |
| `vertical_shear` | False | tested, under 5 %, not adopted |
| `rise_drag` | 0.0 | tested, no effect |
| `JetPlume.ground_effect` | False | |
| `liquid_fraction`, `evaporation_ratio` | 0.0 | |
| `bulk_air_phase_safe_source` | False | phase-safe N2/O2/Ar source bound; tested, not adopted |
| `non_boussinesq` | False | pre-registered and falsified |
| `spread_floor` | False | falsified — FLADIS never enters JETPLU |

The closure is pluggable: `DegadisClosure` is the default, with
`BuoyantClosure` and `UnifiedClosure` in `addons/`.

## 3.2 Can it be the LH₂ model?

Two entry points, both switchable, neither touching the default.

**High level.** `degali.lh2.assess` takes either a pool or a jet and routes
accordingly:

```python
assess(rate=9.2, pool_diameter=9.1, wind=1.6, at_distance=30.0)   # pool
assess(rate=0.26, orifice=0.0254, storage_pressure=6.0,
       wind=2.5, height=0.5, at_distance=4.0)                     # jet
```

It reports regime, lowest flammable height, distances to LFL and
stoichiometric, and the buoyancy crossover, each with its grade, and warns
when the release is outside the range that release type was checked over.

**Low level.** `validation.nearfield.hydrogen_jet(..., corrections=bool)` —
one flag switching all five adopted corrections together, so the as-shipped
and corrected runs come from one code path. Historical states can still be
requested explicitly through the two source-consistency switches.  The
off-by-default `bulk_air_phase_safe_source=True` path advances the source to
the first all-gas bulk-air state near 68 K.  It is an uncertainty bound, not
an adopted correction, because it worsens PRESLHY VG and centre-height MAE.

The preset registry has `DEGADIS 2.1` and `liquid hydrogen`.

## 3.3 Bugs found in this audit

**`check_range` silently ignored an unknown keyword.** It looked each name up
with `.get` and skipped a miss, so a typo turned a scope check into no check
at all with nothing saying so.

Raising instead immediately found a live instance: `lh2.assess` passes
`diameter=pool_diameter` on both paths, and the jet range has no diameter
limit, **so every jet was being checked on four criteria while the caller
believed it was five.** `None` still means "not supplied"; an unknown *name*
is now a `TypeError`. Pinned by
`test_a_scope_check_on_a_name_that_is_not_a_limit_is_an_error`.

### Previously found and fixed, listed so they are not reintroduced

| defect | was | is |
|---|---|---|
| scope warning checked the integration limit | every default call warned about 100 m while every reported number was in range | checks the furthest distance a reported answer relies on |
| one `VALIDATED` range for two campaigns | a jet asked about 30 m passed a limit set by a pool experiment | per release type, and a one-point range says so |
| stale statistics hard-coded in the report string | printed a superseded MG to users | read from `evidence.py` |
| absolute-path guard walked only `src/` | eleven machine paths in the tests; a clone skipped for the wrong reason | walks the tests too |
| `CITATION.cff` had no `authors` | CFF 1.2.0 requires it; Zenodo rejects | added, with a test blocking the placeholder at release |
| REDIPHEM cited as "Rise Data on Heavy Gas Dispersion" | not the title of anything | Nielsen & Ott, Risø-R-845(EN), 1995 |

## 3.4 Edge cases checked and found sound

`Trajectory.at` returns `None` outside the integrated range rather than
extrapolating, and callers skip on `None`. `concentration_at` returns 0.0
there and clamps to 100 % at the nozzle. `cases_from_reduced` raises on a
substance it has no levels of concern for. `Pair.visible` is true when the
centre equals the top sensor, which is the intended boundary.

## 3.5 Known-open, not bugs

**The `averaging` argument may not reach the REDIPHEM reduction.** The
reduction's `note` reads "peak values" although `averaging_time=18.4` was
passed. Worth confirming in the session that built it; if `averaging` falls to
zero internally it affects every REDIPHEM statistic.

**Burro MG is 0.620 against a published 0.811** with the same n=61 and
FAC2 within 0.02. Neither interpolation order (0.620 against 0.625) nor the
thermodynamic backend (CoolProp gives 0.576 and n=63) accounts for the
residual. Recorded, not tuned away.

**The trajectory fault.** At 5–7 m the measurement puts the plume 0.12 m
*below* the nozzle and the model puts it a metre above. The corrections take
0.877 m to 0.658 m at 6 m — about a quarter. The published `1.07 → 0.19 m` is
withdrawn at both ends: neither came from a reconstructible configuration and
they were taken over different subsets besides.

**`alfa1 = 0.0875` is derived from the primary measurement.** Papanicolaou and
List (1988), §3.1, pp. 353–354, give 0.0545 for a jet, 0.0875 for a plume and
0.716 for the pure-plume Richardson number. The commonly repeated 0.0533,
0.0833 and 0.557 are Fischer et al. (1979)'s proposed values quoted in that
paper; the earlier secondary-source attribution is superseded. Historical
parameter dumps still request 0.0833 explicitly. HyRAM's independent 0.082 is
corroboration rather than the derivation.

---

# Part 4 — Claim status

| claim | grade | computed |
|---|---|---|
| Fortran parity, 1e-12, six programs, five EPA cases | A | yes |
| three attribution findings, by line number | A | yes |
| both LH₂ configurations against parameter dumps, 4 s.f. | A | yes |
| LH₂ is buoyant across the whole flammable range | A | yes |
| provenance chain, 24 → 9 trials → 69 arcs | — | yes |
| historical near-field MG 1.121 at-sensor (0.729 centreline), n=69 | B | yes |
| mass-, enthalpy- and momentum-consistent source MG 1.047, VG 1.425, FAC2 0.84, n=62 (0.79–6 m) | B | yes |
| `sigma_z` 0.901 → 1.033, 23 fits, common standard-deviation convention, mean of per-fit ratios | B | yes |
| corrected mean / median signed centre-height error 0.036 / 0.002 m, 23 fits | B | yes |
| Burro MG 0.620, n=61 (published 0.811) | B | yes |
| lift-off RMS 2.9 m, regime 4 of 4 | C | yes |

Eight hypotheses were pre-registered and falsified: entrainment shortfall,
Boussinesq correction, notional nozzle, added mass, rainout source, FLADIS
channel guessing, Coyote pool radius, isotropy split misdiagnosis. They are
results, and `docs/` records each with its prediction and outcome.

At the start of this review **no grade-B claim was computed by anything.** All
of them are now.
