# Pre-registration: the buoyant trajectory at facility distances

> **Later width-normalisation audit:** width ratios in this historical test
> used a measured e-folding width under a sigma label. Multiply them by
> `sqrt(2)` for comparison with JETPLU's standard deviation; see
> `gaussian-width-convention.md`.

> **Current-model notice (2026-09-03).** This is a pre-registration and keeps
> the earlier baselines. The subsequently adopted mixed flashing-source state
> and ground-detachment result are reported in
> [`lh2-model-improvements-2026-09-03.md`](lh2-model-improvements-2026-09-03.md).

Written before any mechanism is implemented. Baseline measured, metrics fixed,
predictions stated. Falsification is the expected outcome for at least one of
the three.

> **Post-audit correction, 2026-09-03.** The mechanism predictions and the
> tests 4/6 baseline below remain as pre-registered. The later claim that all
> other non-horizontal tests were outdoor downward releases does not: DNV GL
> Report 902696 shows that tests 8--15 are a separate closed-room and
> ventilation-mast campaign. Test 15 was wrongly included in the downward
> statistic. The corrected outdoor result is n=5, MG 4.79, VG 12.32 and FAC2
> 0.00. Sections below are updated where this source-classification error
> affected a numerical claim.

---

## 1. The fault

The model lifts a buoyant plume that the measurement leaves on the ground.

It shows up wherever buoyancy dominates: **early** when the source momentum is
weak (PRESLHY trial 20, wind 0.57 m/s, exit 9.8 m/s), and **late** when the
momentum is spent (Spadeadam at 30 m, a release that is strongly
momentum-driven at the orifice). One defect, two regimes.

## 2. The baseline, measured now

Spadeadam tests 4 and 6: 25.4 mm horizontal at 0.50 m, 0.83 kg/s, winds 5.85
and 2.50 m/s. Sensors on 30, 50 and 100 m arcs at 0.1, 1.0 and 1.8 m.

### The primary metric: the vertical gradient

`c(1.8 m) / c(0.1 m)` at the bearing carrying the arc maximum. **This is the
metric because it survives everything else that is uncertain here**: it does
not depend on the absolute concentration, on whether the plume missed the arc
laterally, on the assumed storage pressure, humidity or wind reference height.
A grounded plume gives a ratio at or below one; a lifted plume gives more.

| | 30 m | 50 m | 100 m |
|---|---|---|---|
| measured, test 4 | **0.686** | 1.085 | 0.545 |
| measured, test 6 | **0.886** | 1.800 | — |
| model as shipped, test 4 | 1.628 | 1.259 | 1.061 |
| model corrected, test 4 | 1.034 | 1.048 | 1.018 |
| model as shipped, test 6 | 2.148 | 1.484 | 1.156 |
| model corrected, test 6 | 1.274 | 1.189 | 1.086 |

**Every modelled ratio exceeds one. The two strong measured arcs are below
it.** The 50 m ratios sit on 4.7 and 1.0 vol % signals and are noisy; the 30 m
arcs carry 17.2 and 21.0 vol % and are the ones to judge on.

### Secondary metrics

| | measured | as shipped | corrected |
|---|---|---|---|
| ground concentration at 30 m, test 4 | 17.2 % | 3.17 | 6.80 |
| ground concentration at 30 m, test 6 | 21.0 % | 0.107 | 1.70 |
| modelled plume centre at 30 m, test 4 | below 1.8 m | 3.86 m | 2.64 m |
| modelled plume centre at 30 m, test 6 | below 1.8 m | 10.95 m | 6.96 m |

The wind dependence is already visible and is the signature of the fault:
test 4 at 5.85 m/s is wrong by a factor of 2.5, test 6 at 2.50 m/s by a factor
of 12.

### The regression set, which must not degrade

| | value |
|---|---|
| PRESLHY near field, 66 arcs, at-sensor | MG 1.132 as shipped, 1.396 corrected |
| `sigma_z` ratio, 23 momentum-driven fits | 0.64 → 0.97 |
| trial 10 trajectory against the reference dump | 4 significant figures |
| DEGADIS 2.1 parity | 1e-12, unchanged |

## 3. The three mechanisms, and what each predicts

### M1 — a lift-off criterion

**The model has none.** It lifts continuously from the source. Every model
reviewed in the AEA survey has a threshold: Briggs' `Lp`, a bulk Richardson
number, critical value about 30; HGSYSTEM-MMES uses 20 as the transition from
a grounded to an elevated plume.

    Ri* = g H (rho_a - rho_m) / (rho_a u*^2)

**Predictions.**

1. This is the largest of the three effects, because it acts on *when* the
   plume starts to rise rather than on how fast, and the model currently
   starts at zero distance.
2. Test 6 at 30 m improves more than test 4, because a lower wind gives a
   smaller `u*` and therefore a larger `Ri*`, which delays lift-off further.
3. **The dilution is not touched.** M1 changes no entrainment term, so the
   near-field statistics and the `sigma_z` ratio should be unchanged to within
   rounding. If they move, M1 has been implemented wrongly.
4. It will not be sufficient alone: once lift-off occurs the rise rate is
   still the shipped one, so the 100 m arcs will remain too high.

### M2 — shape-dependent pressure drag

`rise_drag` exists and is a fixed coefficient on `w|w|`. Mack's is
**aspect-ratio dependent**: cylinder-in-crossflow drag with `AR = H / B`, so a
wide flat plume gets more drag than a compact one at the same rise velocity.
That dependence is the whole of the AEA finding — rise suppression is worst
for **wide** sources — and it is what the fixed coefficient cannot express.

**Predictions.**

1. Adding shape dependence changes the trajectory more than the fixed
   coefficient did, and most where the plume is widest, which is the far
   arcs.
2. **It should barely affect the near field**, where the plume is compact and
   `AR` is near one.
3. Added mass was pre-registered and falsified here, and Mack records EFFECTS
   rejecting it for the same reason. If drag also fails, the added-mass
   conclusion generalises: the vertical momentum balance is not where this
   fault lives.

### M3 — entrainment on a buoyant velocity scale

AEA states the trap explicitly and it is a prediction we did not have to
invent:

> to suppress plume rise, the vertical component of plume velocity must be
> reduced. However, reducing the vertical velocity component also has the
> undesirable effect of reducing dilution … suppression of entrainment acts
> against suppression of plume rise

and the escape:

> Modifying the entrainment formula to depend upon a buoyant velocity scale
> rather than vertical component of velocity may circumvent such problems

**Predictions.**

1. **M3 will degrade the near-field dilution** unless the buoyant velocity
   scale is used. That is AEA's prediction, not ours, and it is the sharpest
   test in this document: if a velocity-based suppression improves both the
   trajectory and the dilution, AEA is wrong about the coupling and that is a
   publishable result on its own.
2. On the buoyant scale, the trajectory improves and the dilution does not
   degrade by more than 5 % on MG.

## 4. Rules fixed in advance

- **All three are switchable and default to off.** DEGADIS 2.1 parity is not
  negotiable and is checked after every change.
- **The regression set is run before and after each mechanism**, and both are
  reported. A trajectory improvement bought with a near-field regression is a
  trade, not a fix, and must be presented as one.
- **The metric is the vertical gradient at 30 m on tests 4 and 6.** Chosen
  before implementing anything, for reasons given above. Absolute
  concentration is secondary because the arcs are sparse and censored.
- **Two tests is a thin sample and no combination of mechanisms will be
  claimed to be validated on it.** The most that can be claimed is that a
  mechanism moves the model in the measured direction without costing the near
  field.
- **A mechanism that works is not adopted until its coefficient has a source.**
  Mack gives `C_g` 0.035 for buoyant plumes against 0.025 for jets, and the
  cylinder drag data is cited; Briggs gives 20–30 for the threshold. A number
  fitted here is a fit and must be labelled one.

## 5. What would falsify the whole reframing

If a mechanism fixes the 30 m gradient on both Spadeadam tests **and** the
PRESLHY near field holds, the fault was the buoyant branch and this was the
right diagnosis.

If none of the three moves the gradient below one, the fault is not in the
vertical momentum balance or the entrainment, and the remaining candidates are
the source term at the transition to the buoyant regime, or the assumption of
a single Gaussian plume where the measurement sees a bifurcated one — Hansen
reports bifurcation on these very tests.

**That second outcome is at least as likely as the first.** Recording it now so
that it is a result rather than a disappointment.

---

# Outcomes

Recorded after running. Three predictions falsified, one mechanism found that
was not predicted.

## M1 — lift-off criterion: FALSIFIED

Implemented as `JetCoefficients.liftoff_richardson`, a smoothed threshold on
`Ri* = g H (rho_a - rho_m) / (rho_a u*^2)`.

**It fires nowhere.** `Ri*` runs 187 to 363 over the whole trajectory of
Spadeadam test 6, against Briggs' threshold of 20 to 30:

| x (m) | 2 | 5 | 10 | 20 | 30 |
|---|---|---|---|---|---|
| `Ri*` | 251 | 363 | 337 | 239 | 187 |

Setting the threshold anywhere in Briggs' range changes the 30 m plume centre
in the fourth decimal. **The criterion agrees the plume should lift.**

Prediction 1 said this would be the largest of the three. It is the smallest.
The reasoning was that the model had no criterion where every other model has
one, and that turned out to be true and irrelevant: the gap was real, and
filling it does nothing because the plume is an order of magnitude past the
threshold.

## M2 — shape-dependent pressure drag: FALSIFIED

Implemented as `JetCoefficients.shape_drag`, interpolating `Cd` between a
circular cylinder at 1.2 and a flat plate at 2.0 on the aspect ratio
`sigma_z / sigma_y`, acting on the plume width.

At a coefficient of 2.0 — larger than any classical drag coefficient justifies
— it moves the 30 m plume centre by ten per cent, 6.96 m to 6.29 m. Prediction
2 (little near-field effect) holds; prediction 1 (larger than the fixed
coefficient, concentrated at the far arcs) holds in direction and fails in
magnitude.

Prediction 3 now applies: **added mass was falsified here and drag is
falsified here, so the vertical momentum balance is not where this fault
lives.** EFFECTS reaching for drag after rejecting added mass is corroboration
that added mass is wrong, and not evidence that drag is right.

## M3 — vertical velocity in the shear entrainment: FALSIFIED

Two per cent on the 30 m plume centre.

## AEA's coupling prediction: NOT TESTED, and now moot

M3 was too small to move the dilution either way, so the sharpest test in this
document was not reached. It stands for whoever revisits the entrainment.

## The mechanism that works, which was not predicted

**`JetPlume.ground_effect` was in the package throughout and had never been
switched on for the liquid hydrogen path.** It reduces the entrainment
perimeter over the part of the cross-section in contact with the ground, and
scales the buoyancy by the fraction that has cleared.

By the model's own geometry the Spadeadam plume is in contact from 2 m to
about 48 m, so it applies over nearly the whole comparison.

| at 30 m, test 4 | off | on | measured |
|---|---|---|---|
| plume centre | 2.64 m | **1.59 m** | below 1.8 m |
| gradient `c(1.8)/c(0.1)` | 1.034 | **0.856** | 0.686 |
| c at 1 m | 6.92 % | **13.55 %** | 8.4 mean, 17.2 peak |

**The gradient goes below one for the first time** — the pre-registered
metric — at all three arcs, and the concentration lands between the measured
mean and peak instead of a factor of two under the mean.

**And it costs nothing.** The near field is unchanged: MG 1.132 against 1.130,
`sigma_z` ratio 0.64 → 0.97 against 0.64 → 0.95. The PRESLHY plume at 0.35 to
6 m is barely in contact, so the term is inactive there.

That it costs nothing is not luck and is worth stating: **ground contact does
not suppress rise through the vertical velocity.** It removes buoyancy the
ground is holding and entrainment area the ground is blocking. AEA's coupling
— suppressing rise costs dilution — is a statement about velocity-based
suppression, and this mechanism is outside it. That is a testable claim about
why the three predicted mechanisms failed and this one did not: **the fault
was a geometric constraint, not a force balance**, and all three predictions
were aimed at the force balance.

## What is still wrong

Test 6, the low-wind release, improves by a factor of three and remains
wrong: gradient 1.251 against a measured 0.886, concentration 5.43 % against a
measured mean of 15.4 %.

Same orifice, same height, the same 0.83 kg/s; the only difference from test 4
is wind, 2.50 m/s against 5.85. **The wind dependence of the rise is still too
strong.** Roughness is not the explanation — raising it from 0.001 to 0.03 m
makes both tests worse.

The open question is narrower than it was. Not "the buoyant branch is wrong"
but "what holds a plume down at 2.5 m/s that ground contact does not capture".

Four candidates were tried and none is it:

| tried | result |
|---|---|
| shape drag on top of ground contact | 2 % |
| vertical shear on top of ground contact | under 1 % |
| site roughness 0.001 to 0.03 m | **worse**, both tests |
| Pasquill class C through F | C helps and is not defensible at 2.5 m/s in December; E and F make it worse |

The aspect ratio rules out AEA's wide-source mechanism: both plumes are
`sigma_z / sigma_y` about 0.9, near circular, and AEA's unsolved case is the
wide flat one.

**And the release is not buoyancy conserving**, which is where the literature
stops. It leaves the orifice at 5.35 kg/m3, four times ambient; even the
saturated vapour at 1.33 is denser than the 1.27 air. It becomes buoyant
within half a metre, peaks near one, and has lost four fifths of its buoyancy
by thirty. There is no single buoyancy flux to put in a correlation.

AEA say so twice -- "lift-off parameters based on non-dimensional fluxes are
of limited use for non-buoyancy conserving flows", "the lift-off distance
correlations are unlikely to be valid for non-buoyancy conserving flows" --
and their own programme, which existed to solve exactly this for hydrogen
fluoride, closed reporting it unsolved:

> We have found it difficult to improve on the estimation of plume rise from
> wide sources ... without detrimentally affecting the dilution

So the remaining item is a known hard problem in the literature and not an
oversight here. That is worth stating in the paper rather than apologising
for.

## Discipline notes

- `ground_effect` was found by trying an option, not by predicting it. It is
  recorded as exploratory. The three that were predicted are recorded as
  falsified whether or not that is the flattering account.
- Every option added here defaults to off and DEGADIS 2.1 parity is unchanged.
- Two tests are not a validation and nothing here is presented as one.

---

# The aggregate, through the packaged reader

`validation/spadeadam.py` reads the fifteen releases and pairs the arcs. Six
arcs from the two horizontal ones can be modelled:

| `ground_effect` | n | MG | VG | FAC2 | within 2× |
|---|---|---|---|---|---|
| off | 6 | 3.224 | 6.00 | 0.17 | 1 of 6 |
| **on** | 6 | **1.513** | **1.60** | **0.67** | **4 of 6** |

The bias halves, the variance falls fourfold, and FAC2 goes from one arc in
six to four. **FAC2 and VG pass Hanna's criteria**, from a start where none of
the three was close. Every arc improves
individually, not only the aggregate.

Six arcs from two releases is thin and this is not a validation. It
establishes a direction and a size at a distance nothing in this package had
reached.

## The sample arithmetic, decided before any model ran

| | count |
|---|---|
| releases | 15 |
| closed-room/ventilation-mast tests — different atmospheric source | −8 |
| outdoor downward releases, which this path does not model | −5 |
| **modelled** | **2** |

Both exclusions are properties of the experiment. The reader returns and
classifies all fifteen so the arithmetic stays visible. `spadeadam.compare`
**raises** on a downward release, and both outdoor paths raise on a
closed-room trial, rather than substituting the supply-nozzle conditions for
the ventilation-mast outlet.

# The downward releases, and independent corroboration

Five outdoor downward releases carry real far-field signal. The other eight
non-horizontal rows belong to the closed-room/ventilation-mast experiment and
cannot be added to this comparison.

**There is no jet to integrate.** At 0.32 m with a 25.4 mm orifice the fluid
reaches the ground inside the jet's own development length, and
`initial_conditions_directed` refuses them for exactly that reason — a guard
that was already in the package and is right. The fluid arrives essentially
undiluted, so what is left is a ground-level source of unknown footprint.

**The footprint was measured before a value was chosen.** Sweeping it
eighty-fold, 0.05 m to 4 m, moves the 30 m concentration by seven per cent and
the plume centre by half a metre. The unknown does not matter, which is what
makes the five outdoor tests usable.

They run through `addons.LiftoffPlume`, the URAHFREP ground-truncated buoyant
plume, and not through `JetPlume`. That makes the two halves of the Spadeadam
comparison independent code paths:

| comparison | n | MG | VG | FAC2 |
|---|---|---|---|---|
| horizontal, `JetPlume`, ground contact off | 6 arcs | 3.22 | 6.00 | 0.17 |
| downward, `LiftoffPlume` | 5 tests | **4.79** | 12.32 | 0.00 |

**Two buoyant plume models, written by different people for different
purposes, over-predict dilution by more than a factor of three on the same
campaign.** The magnitudes are not identical after the source correction; the
common direction remains independent of either implementation.

And it is worth naming what the URAHFREP model is: **AEA Technology's own**,
built for this problem, and reported by them as over-predicting rise. It does
so here too, on liquid hydrogen, thirty years later.

The far-field outdoor sample consists of two horizontal and five downward
tests, evaluated through their respective source paths.

# NASA, and a third confirmation

**Tests 3 and 7 cannot extend the sample, and the reason is worth recording.**
Their release conditions are published — they sit in Mack's Table 2 beside the
others — but Witcofski's Table 1 records test 3 as yielding motion picture
only, and neither test appears in Table 3 or Table 4. **Conditions are not
measurements**, and the lift-off sample stays at four. That is now asserted in
the suite so the mistake cannot be made twice.

**Table 3 is a vertical profile and had never been used.** Maximum
concentration at the furthest tower row, 33.8 m, at 1, 9.4 and 18.6 m, for the
four tests that carry data. The package used only Table 4, one height per
test. Transcribed from the printed table, not digitised from a figure.

| test | wind (m/s) | measured peak at | model centre | ratio |
|---|---|---|---|---|
| 5 | 6.30 | 9.4 m | 6.2 m | 0.66 |
| 4 | 3.35 | 9.4 m | 14.8 m | 1.57 |
| 6 | 2.20 | ~14 m | 23.9 m | 1.71 |
| 2 | 1.55 | ≥ 18.6 m | 40.3 m | > 2.2 |

**The model over-predicts the height and monotonically worse as the wind
falls.** That is the Spadeadam result on a different campaign, four decades
and an ocean apart, at ten times the mass flow, from a 9.1 m pond instead of a
25 mm orifice, through the ground-level buoyant path rather than the jet.

**Three campaigns, two plume models, one fault.**

Test 5 is the only case the model puts too low, and it is the highest wind at
6.3 m/s, twice any other. The paper says the same from the data side: its
cloud travelled downstream fastest and gave the highest far-field
concentration, 29.2 % at 9.4 m.

# M4 — entrainment on a buoyant velocity scale

The one mechanism AEA recommended and nobody tried. `sqrt(g' H)` does not
scale with the wind and every other mechanism here does, which is why it is
the only one that changes the **shape** of the error rather than its size.

| coefficient | near MG | near VG | near FAC2 | `sigma_z` | far z6/z4 |
|---|---|---|---|---|---|
| 0 | 1.132 | **19.15** | 0.79 | **0.97** | 2.46 |
| 2 | 0.950 | **2.76** | 0.83 | 1.08 | 2.22 |
| 4 | 0.938 | **1.73** | 0.85 | 1.28 | 1.88 |
| 8 | 1.107 | 1.75 | 0.82 | 1.75 | 1.47 |

**AEA's coupling prediction is falsified on their own escape route.** The
pre-registration said the near-field dilution would degrade unless the buoyant
scale was used, and set five per cent on MG as the tolerance. It does not
degrade: MG moves *towards* unity and the variance collapses sevenfold, on
data the three corrections were never tuned against.

The far field improves in the way that matters: `z(2.3 m/s) / z(5.0 m/s)` falls
from 2.46 towards 1.5, and test 6's concentration at 1 m goes from 4.6 to
10.0 vol % against a measured mean of 15.4.

**It is not adopted.** The coefficient has no source, and this document says a
mechanism that works is not adopted until it does. The cost is visible in the
same table: `sigma_z` is the one quantity the three jet corrections were tuned
on, and this pushes it from 0.97 to 1.08 at a coefficient of 2 and 1.28 at 4.
That is a trade.

## Next

Two things, in order.

**Find or derive the coefficient.** The structure is right and the value is
fitted. Candidates: the plume-depth definition in `sqrt(g' H)` is a choice
(`2 delta sigma_z` here) and a different one absorbs part of the factor;
Mack's `C_g` of 0.035 for buoyant plumes against 0.025 for jets is a measured
ratio in the same place in the equations. Until then it is a fitted
coefficient with the right shape, and must be described that way.

**Calibrate on held-out data if no source exists.** Seven outdoor Spadeadam
tests and four NASA tests are now readable. Fitting on one campaign and
testing on the other is the minimum that would make a fitted coefficient
defensible. The ventilation-mast tests are a third source class, not extra
samples for either outdoor path.
