# Data that can address the trajectory weakness

The weakness is one thing: **the model lifts the plume and the measurement
does not.** Everything else — source, dilution, `sigma_y`, `sigma_z` after
correction — reproduces.

Two problems compound it. The fault itself is unfixed, and the validation
envelope is 0.35–6 m while a facility assessment needs tens to hundreds of
metres. The project folder turns out to hold material for both.

---

## 1. The decisive find: FFI 20/03101, the Spadeadam LH₂ tests

`21-03101-잠금_해제됨.pdf` is not a PDF — it is extracted text of **FFI-RAPPORT
20/03101**, *Large scale leakage of liquid hydrogen (LH2) — tests related to
bunkering and maritime use*, Aaneby, Gjesdal and Voie, Norwegian Defence
Research Establishment, 13 January 2021.

It contains the **per-sensor data tables**, not just the narrative. Extracted
to `reference/spadeadam/`. A later source audit recovered the two underlying
DNV GL reports from the official NPRA repository and found that the fifteen
tests are two different experiments, not one homogeneous release series:

| | |
|---|---|
| tests | **15**: 7 outdoor, 8 closed-room/ventilation-mast |
| sensor rows | **519** (357 at 30 m and beyond) |
| outdoor arcs | 30, 50, 100 m |
| closed-room field sensors | 2–50 m from the TCS-floor origin; up to 6.94 m high |
| per row | mean, peak, min, standard deviation, vol % |
| conditions | orifice, orientation, mass flow, wind speed high/low, direction, line pressures |

Release conditions:

| | range |
|---|---|
| mass flow | **0.162 – 0.833 kg/s** |
| orifice | 13 and 25.4 mm |
| wind | 2.3 – 6.7 m/s |
| atmospheric configuration | tests 1,2,3,5,7 outdoor downward; 4,6 outdoor horizontal; 8–14 ventilation mast; 15 no mast after the Test 14 explosion |

**Why the outdoor subset is the right dataset.** Our jet evidence spans
0.084–0.285 kg/s over 0.35–6 m. Tests 4 and 6 reach **three times the flow at
seventeen times the distance**, on the same substance, with the same 25.4 mm
orifice, and are horizontal releases at 0.50 m — the identical geometry to the
PRESLHY trials the model was fitted on. Tests 8–15 cannot extend that sample:
their tabulated nozzle conditions describe a release inside a 24 m3 container,
not the atmospheric source at the ventilation-mast outlet.

57 of 519 rows read above 100 vol % and are flagged `over_range`; they are
sensor faults and are excluded.

### The test it enables, run now

Spadeadam Test 6 is our exact configuration at five times the distance:
25.4 mm horizontal at 0.50 m, 0.833 kg/s, wind 2.5 m/s.

| x | model plume centre | model at 0.1 / 1.0 / 1.8 m | measured |
|---|---|---|---|
| 30 m | **11.0 m** as shipped, **7.1 m** corrected | 0.1 / 1.6 vol % | **21 %** |
| 50 m | 16.7 / 12.3 m | 0.04 / 0.42 | 2 % (plume missed the arc) |
| 100 m | 27.3 / 22.4 m | 0.02 / 0.10 | none recorded |

**At 30 m the model under-predicts ground-level concentration by a factor
between ten and two hundred, because it has lifted the plume seven to eleven
metres out of an array that reaches 1.8 m.** The measurement puts 21 vol % at
the ground there.

That is the same fault seen at 6 m in PRESLHY — a factor of two on the plume
height — appearing at 30 m as a factor of a hundred on concentration. It is
the clearest statement of the defect the project has, and it is at a distance
that matters for a facility.

Caveats, stated plainly: the wind reference height, the storage pressure and
the humidity are assumed, the arc has only five sensors 20 m apart so the
plume can miss them, and the reported 21 % is a peak rather than an average. A
factor of a hundred is far outside all of that, and the sign is unambiguous.

### What the high sensors do — and do not — add

Tests 9–15 carry sensors at radii 11.2–13.8 m and heights to **6.94 m**, but
those coordinates are relative to the release inside the TCS. The atmospheric
source is a 450 mm vertical mast whose outlet is at a different origin and
whose flow is cold gaseous hydrogen. Treating these points as a profile above
an outdoor ground release was a source-geometry error.

DNV GL Report 902696 provides inferred outlet source conditions for tests
8–14: 100 vol % hydrogen, −209 to −142 °C, 6.4–23.5 m/s and 0.183–0.673 kg/s.
They are now transcribed in `reference/spadeadam/mast_conditions.csv`.
`compare_mast` rotates the sensor coordinates into the wind frame, starts a
vertical cold-gas jet at the mast outlet and treats sub-0.5 vol % readings as
censored. The outdoor paths continue to refuse these tests.

The separate pre-registered comparison contains 111 downwind external sensor
rows. It produces zero false flammable predictions on nondetects, but also
predicts effectively zero at all 18 detected sensors. Ten detections lie
within 3.5 m crosswind of the modelled axis; the model centre is more than 8 m
above every detected sensor. The mast data therefore add an independent
source geometry showing the same plume over-rise, not an extra sample for the
outdoor performance statistic.

---

## 2. NASA White Sands, more than is currently used

`Andreas2023.pdf` and `MackExtension...Preprint.pdf` report vertical
concentration profiles from the sample bottles at **9.1, 18.3 and 33.8 m** for
test 6, with the plume centreline height at each. The project currently uses
Witcofski Tables 1–4 for the lift-off *height* only, at one distance.

Two things are available and unused:

- **Vertical profiles at three distances**, which constrain the trajectory
  rather than a single lift-off height. Mack reports the maximum at 9.1 m as
  60 vol %, at 18.3 m as 45 %, at 33.8 m as 19 %.
- **Two more tests.** Mack's Table 2 gives conditions for tests 2, 3, 4, 5, 6
  **and 7** — 4.23 kg/s at 4.5 m/s, and 1.66 kg/s at 3.1 m/s. The lift-off
  comparison uses four. Six would take the weakest claim in the package from
  n = 4 to n = 6.

The sample-bottle numbers are in figures rather than tables, so they need
digitising. The paper states the agreement in the text with specific values,
which is enough to check against.

---

## 3. Wind-tunnel data on exactly this failure mode

The four AEA Technology URAHFREP reports document a programme on buoyant plume
lift-off, and they describe our problem in almost the same words.

`aeatliftoffmodelling.pdf` concludes:

> the recent experiments of Hall and Walker indicate significant suppression of
> plume rise for wide sources, above that predicted by the simple model … models
> making similar plume rise assumptions to our simple model are likely to
> overestimate plume rise and hence underestimate ground level concentrations

And it records the trap we would walk into:

> to suppress plume rise, the vertical component of plume velocity must be
> reduced. However, reducing the vertical velocity component also has the
> undesirable effect of reducing dilution … suppression of entrainment acts
> against suppression of plume rise by virtue of the more concentrated plume

**That is a pre-registered prediction for free.** Any correction that suppresses
rise by damping vertical velocity should be expected to spoil the dilution,
which currently reproduces. If it does not, that is informative; if it does,
the mechanism is wrong and we know before running.

The same report offers the way out:

> Modifying the entrainment formula to depend upon a buoyant velocity scale
> rather than vertical component of velocity may circumvent such problems

The reports also give a lift-off **criterion** — Briggs' `Lp`, a bulk
Richardson number, critical value around 30; HGSYSTEM-MMES uses `Lp = 20` as
the transition from grounded to elevated. Our model has no such criterion: it
lifts continuously from the source.

The Hall and Walker data appear as figures against non-dimensional buoyancy
flux `F/u³L`, so they need digitising, but the correlations are quoted
numerically: first onset of rise at `F/u³W ≈ 0.01`, ground-level concentration
down to 10–20 % of maximum at 0.035, below 5 % above 0.3.

---

## 4. A mechanism, from EFFECTS

Mack's preprint is the one document that both diagnoses our fault and proposes
a specific fix:

> for strongly buoyant plumes, it was observed that vertical plume speeds are
> overpredicted by some integral models … In EFFECTS, not the 'added mass'
> concept is applied but an additional drag term is introduced … Applying the
> additional pressure drag then shows a good prediction of the plume trajectory
> during rise.

**This matters because added mass was pre-registered and falsified here**, and
Mack says explicitly that EFFECTS rejected added mass in favour of drag. Two
independent groups reaching the same conclusion about added mass is
corroboration; the drag is the untried half.

The `rise_drag` option exists in `JetCoefficients` and is a fixed coefficient
on `w|w|`. Mack's is **shape-dependent** — cylinder-in-crossflow drag with the
aspect ratio `AR = H/B` — so a wide flat plume and a compact one get different
drag at the same rise velocity. That is not what is implemented, and it is
precisely the dependence the AEA wind-tunnel work says is missing, since their
whole finding is that **wide** sources have their rise suppressed most.

Mack also extends the shear entrainment with the vertical velocity component,
which is our `vertical_shear` flag, and gives `C_g` = 0.025 for jets and heavy
gas against **0.035 for buoyant plumes** — a specific number for a coefficient
we currently leave at the heavy-gas value.

---

## 5. What this adds up to

**Three new constraints on the trajectory, all in hand:**

| source | distance | what it constrains |
|---|---|---|
| Spadeadam outdoor, tests 4 and 6 | 30–100 m | ground-level concentration under a horizontal jet |
| Spadeadam ventilation mast, tests 8–14 | 5–50 m from the mast coordinate frame | zero false LFL nondetects, but 18 detected signals missed through over-rise |
| NASA sample bottles, test 6 | 9.1 / 18.3 / 33.8 m | vertical profile and centreline height |

**Two mechanisms with a stated form, neither tried:**

- shape-dependent pressure drag (`AR = H/B`, cylinder in crossflow)
- entrainment on a buoyant velocity scale rather than the vertical component

**One criterion the model lacks entirely:** a lift-off threshold. Briggs `Lp`
around 20–30. The model currently lifts from the source with no transition,
which is the most likely reason it is at 11 m by 30 m downwind.

**And a pre-registered prediction, from AEA:** suppressing rise through
vertical velocity will damage the dilution that currently works. Write that
down before testing anything.

## 6. Where to be careful

The Spadeadam arcs are sparse — five sensors 20 m apart at 50 m, two 35 m
apart at 100 m — and Hansen documents the plume missing them repeatedly.
Reported maxima are **censored from below**. Use them as lower bounds on the
ground-level concentration, which is all that is needed here: the model
predicts 0.1 vol % where 21 % was measured, and a censored 21 % is still 21 %.

Five outdoor tests are **downward** releases that impinge on the ground. That
is a different source term from the horizontal jet path, and the ground-level
buoyant route is used for them. Tests 4 and 6 are the horizontal like-for-like
comparison. Tests 8–15 are neither group: they are closed-room tests and their
supply-nozzle orientation must not be used as the atmospheric source.

The validation suite now reads and classifies all fifteen tests, exercises the
two outdoor source paths, loads the mast outlet estimates for tests 8–14, and
raises if a closed-room trial is sent through either outdoor path.

## 7. Not in the folder, worth looking for

- **Hall and Walker's wind-tunnel tables.** The AEA reports plot them; the
  underlying data would let the shape dependence be fitted rather than assumed.
- **Witcofski Table 2**, time-resolved grab bottles, still scrambled by OCR.
- **The Spadeadam concentration time series.** The report gives mean, peak, min
  and standard deviation per sensor, not the traces. The traces would allow a
  transient comparison, which nothing in this project has yet achieved for
  hydrogen.
- **Selected raw DNV channels.** The official NPRA repository is reachable,
  and both DNV reports (853182 and 902696) are now stored in the project. The
  per-test raw archives range from hundreds to tens of thousands of megabytes;
  download only the channel groups needed for a pre-specified transient test.
