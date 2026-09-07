# Data inventory

Everything available to this project, what state it is in, and what it can and
cannot support. Written because most of the effort in a model-evaluation study
goes into finding out what the data actually is, and that finding-out is
usually thrown away.

---

## 1. DEGADIS 2.1 itself

Vendored in `reference/`. The 1989 EPA release with the 2012 SCRAM updates:
83 Fortran source files, 5 include files, 10 makefiles, and the five test
cases with their golden output.

| item | state |
|---|---|
| source | complete; builds with gfortran after eight portability patches |
| test cases | `B9`, `B9T`, `EX1`, `EX2`, `EX3` with `.INP`, `.ER1/2/3` and `.LIS` |
| user's guide | `degadis2.pdf`, OCR'd; equations are degraded, the input format is not |
| readme | `DEGADISreadme.pdf`, complete input-file specification |

**Usable without reservation.** This is the primary reference and everything
in `docs/validation.md` rests on it.

---

## 2. REDIPHEM — heavy gas field trials

Risø National Laboratory. 313 trials with measurements, of which the field
scale is what matters here.

| series | trials | substance | release | usable |
|---|---|---|---|---|
| BURRO | 8 | LNG | pool | yes |
| COYOTE | 3 | LNG | pool | yes |
| EAGLE | 4 | N₂O₄ | pool | measurements yes, properties no |
| TORTOISE | 4 | ammonia | jet | yes |
| FLADIS | 14 | ammonia | jet | **no — see below** |
| TI | 2 | R12 | puff | partly |
| LATHEN | 51 | propane | jet | 18, and the result is not meaningful |
| HAMBURG / WSL / TNO | 227 | SF₆ etc. | wind tunnel | outside DEGADIS's geometry |

### Format

Each trial is a directory of four files.

`SPECS.DAT` — conditions, one per line, `label : value`. Maps almost one to
one onto a DEGADIS deck.

`SETUP.DAT` — one line per channel: `id x y z type`.

`DATA.DBF` — **not a dBase file despite the extension.** A raw little-endian
`float32` stream: a channel count, then that many channel identifiers, then a
matrix with one row per sample whose first column is the time of day in
seconds. `-1234` marks a missing reading.

`CHANDEF.DAT` — one level up, in blocks of five lines: type number, quantity,
unit, instrument, notes.

### Traps found, all of which silently corrupt a result

**Channel numbers are per series.** Burro uses 21 and 22 for concentration,
Eagle 21 and 24. There is no global list. Hardcoding one drops Eagle entirely
without saying so.

**FLADIS ships no `CHANDEF.DAT`.** A plausible fallback list selects channels
reading 302 and 23.5 at 20 m downwind, which are not concentrations at all.
That produced a clean-looking FLADIS result — MG 0.89, VG 1.02, FAC2 1.00 —
that meant nothing. The reader now returns an empty set rather than guessing,
and FLADIS yields nothing until someone supplies the definitions.

**The release time is three lines and only the first is labelled.** Reading
the labelled line gives the hour and silently loses 37 minutes, which puts
every arrival time out by that much.

**Anemometer height is sometimes in the field name.** Thorney Island writes
`site average windspeed at 10m` rather than a separate height field.

**Sensor elevations differ between series.** Thorney Island samples at 0.4 m,
FLADIS at 1.5 m, Burro at 1, 3 and 8 m. Asking for "1 m" drops whole series
while reporting "no usable measurements".

**Crosswind sampling is thin.** Burro has at most two crosswind positions per
arc, and they are a symmetric pair — at 49 m on B7 one reads 15 mole per cent
and its mirror 0.09, so the plume was well off centre. Two points cannot
constrain a width, and the arc maximum is a *lower bound* on the true
centreline. Desert Tortoise has four to seven and is the only series that can
resolve lateral spread.

**The last line of some `SETUP.DAT` carries a DOS end-of-file byte.**

### A trap in the harness, not the data

Five times in this work a tool reported a true-sounding reason for a false
cause, and the last one was the harness itself.

`compare` chose the source route by default rather than by release type, so
Thorney Island -- a 14 m cylinder of Freon released all at once -- went
through the steady pool route and produced a source of **3970 kg/s over a
2050 m radius**. That failed a degeneracy check downstream. Worse, when the
run failed the series driver threw away the reason and rebuilt a *pool* deck
to ask what was wrong; a pool deck's verdict on a puff has nothing to do with
a van Ulden momentum balance that did not converge.

The trial was therefore reported as **"no usable measurements"** while
carrying eighteen arcs of them.

Both are fixed: the route follows `release_type`, and the failing case carries
its own reason. Thorney Island now reports the model limit it actually hits.
Where the message is still "no usable measurements" it is now accurate -- the
Lathen trials it names have every concentration channel at or upwind of the
release, which a test checks.

### Not redistributable

Tests that use it read `$REDIPHEM_ROOT` and skip cleanly without it.

---

## 3. Liquid hydrogen — PRESLHY and HSL

The new material, and a different problem from LNG. Three datasets, plus
literature.

### 3a. HSL E3.5 "rainout" trials — DOI 10.35097/1481

24 LH₂ release trials at the HSL facility, September 2019, on a 32 m concrete
pad. **This is the dataset with usable dispersion measurements.**

Already reduced into two CSVs in the project folder:

`lh2_e35_conditions_v2.csv` — 24 rows, 22 columns. Per trial: orientation
(horizontal / vertical up / vertical down), release height (0.5 or 1.5 m),
orifice (6, 12 or 25.4 mm), tanker pressure (1 or 5 barg), wind mean/max/gust
with reference height, air temperature, relative humidity, and **three
independent flow estimates** (meter peak, meter mean, and the report's Table 4
value) with a `flow_source` column recording which was used.

| condition | values |
|---|---|
| orifice | 6, 12, 25.4 mm |
| tanker pressure | 1 barg (9 trials), 5 barg (15) |
| release height | 0.5 m (13), 1.5 m (11) |
| orientation | horizontal (19), vertical (5) |
| wind speed | 0.57 to 4.17 m/s — **all very light** |
| air temperature | 12.8 to 18.2 °C |
| relative humidity | 51 to 81 % |
| flow rate | 84 to 285 g/s where recorded; 3 trials have none |

`lh2_e35_farfield_v2.csv` — 686 rows, one per (trial, sensor). Ten stands at
three heights (0.5, 1.5, 2.5 m), 30 Dräger devices per trial. Columns include
peak and mean H₂ concentration, an oxygen-depletion minimum, and — crucially —
a `saturated` flag and `n_at_ceiling` count.

**The measurement ceiling is the central problem.** The Dräger sensors
over-range above 4 % H₂ by volume. 57 of 686 rows are flagged saturated and
the maximum recorded value is exactly 4.00. Hydrogen's lower flammable limit
is 4 %. **So the instrument saturates exactly at the concentration a
flammability assessment cares about**, and every saturated reading is a lower
bound, not a measurement. Any comparison has to treat those as censored data
rather than as values.

**The far field is short.** The report states the furthest centreline sensor
is 14 m from the release point. Burro's arcs run to 800 m. Fourteen metres is
inside what DEGADIS would call the secondary source for most of these
releases, which is a serious structural obstacle rather than a detail.

**Stand coordinates are in a figure, not a table.** Table A4 of the report
gives serial number, stand number, height and thermocouple name — but not `x`
and `y`. Those are in Figure A4 as a sketch, and two layouts were used
depending on wind direction. So the downwind distance of each stand is not
machine-readable from what is here, and 534 of 686 rows carry a detection with
no distance attached to it.

Raw source: 24 `.xlsx` files (19 MB) plus 24 `.mp4` (10.5 GB), inventoried in
`lh2_inventory_1481.txt`. The videos are not needed.

### 3b. PRESLHY E3.4 pool trials — DOI 10.35097/1319

10 unignited LH₂ pool spills onto four substrates: gravel (3), sand (3),
concrete (3), water (1), March–April 2020, 383 MB of `.xlsx`.

Each file has two sheets: 79 thermocouple columns in Kelvin, and hydrogen
concentration at **three heights only — 35, 45 and 55 cm** — with sonic
anemometer wind in three components. Sampling is roughly 10 Hz over 90 000
rows.

Substrate is the independent variable, which is exactly the variable DEGADIS
represents through a single ground heat-transfer coefficient. That makes this
dataset interesting for the source term and close to useless for dispersion:
the concentration measurements span 20 cm of height and no downwind distance
at all.

### 3c. PRESLHY E3.1a discharge trials — DOI 10.35097/1187

Cold hydrogen discharge through 0.5, 1, 2 and 4 mm orifices at 80 K and
300 K. 1.2 GB in nine nested zips, **not extracted**. Release-rate and
near-orifice physics, not dispersion. Relevant to a source term, not to a
downwind comparison.

### 3d. Literature

| file | what it is |
|---|---|
| `PRESLHY_D3_6_...pdf` | the E3.5 report: conditions, layout, instrumentation, the 4 % ceiling |
| `rr985/986/987.pdf` | HSE research reports on LH₂ releases |
| `witcofski1984.pdf` | NASA LH₂ spill trials — the classic dataset |
| `aeat*.pdf` (6) | AEA Technology: HF thermodynamics, liftoff modelling, field-trial comparisons |
| `HGSYSTEM_docs.pdf` | HGSYSTEM, the model family that handles buoyant lift-off |
| `Mack*.pdf`, `MackExtension...` | buoyant plume rise, an extension to EFFECTS |
| `OlavR2022/2023.pdf`, `Andreas2023.pdf`, `applsci*.pdf` | recent LH₂ dispersion modelling |
| `paper_44`, `paper_161_0`, `115`, `160`, `xxiiipaper65`, `200822_1` | conference papers |

The presence of the liftoff and buoyant-plume-rise material is itself
informative: it says the people who worked on this concluded that a dense-gas
model alone does not cover LH₂.

---

## 4. SMEDIS

Three batches of spreadsheets from the EU model-evaluation exercise, plus
`equivsrc.txt`. Used for one thing and it mattered: the **equivalent source
terms** for Desert Tortoise and two Lathen trials — the plume state at the
point the flashing aerosol has fully evaporated. Those licensed the
first-principles flash calculation in `validation/flashing.py`, which
reproduces them to within 5 % and was then applied where none was published.

---

## 5. What the data will and will not support

**Supported.** Reproducing DEGADIS exactly. Evaluating it against LNG pool
spills and ammonia jets at the ground-level centreline. Testing its vertical
structure on Burro and its lateral structure on Desert Tortoise. Computing
flashing-jet source terms and checking them against published ones.

**Not supported.** FLADIS, without channel definitions. Eagle, without N₂O₄
properties — the substance dissociates to NO₂ so its effective molecular
weight is temperature-dependent, which is a real modelling problem and not an
oversight. Crosswind structure on Burro or Coyote. Anything above 4 % H₂ in
the LH₂ far field. Downwind distance for the LH₂ stands, unless the figure is
digitised or the coordinates found elsewhere.

**The honest summary of the LH₂ data**: the conditions are well characterised
and the concentrations are censored exactly at the threshold of interest, over
a fetch shorter than the model's own source region, with sensor positions that
are not in machine-readable form. That is a hard dataset to validate a
dispersion model against, and saying so up front is better than discovering it
halfway through.
