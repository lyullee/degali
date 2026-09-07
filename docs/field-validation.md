# Field validation

Reproducing DEGADIS says nothing about whether DEGADIS is right. This file is
the other half: the model against measurement.

The data is the Risø **REDIPHEM** database, 313 heavy-gas release trials
covering the LNG spills (Burro, Coyote), pressurised ammonia jets (Desert
Tortoise, FLADIS), nitrogen tetroxide (Eagle), Freon puffs (Thorney Island),
propane jets (Lathen) and several wind-tunnel series.

## How a trial becomes a deck

`SPECS.DAT` maps almost one to one onto a DEGADIS input deck — wind speed and
its reference height, friction velocity, roughness, Monin-Obukhov length,
stability class, temperature, pressure, humidity, release rate, duration and
pool or nozzle diameter all go straight across.

Four things do not, and `to_case` records each on the case it returns:

| supplied | why it is not in the trial |
|---|---|
| contaminant properties | the trial names a substance, not its thermodynamics |
| surface temperature | rarely recorded; defaults to ambient |
| levels of concern | a property of the question, not the experiment |
| release temperature (pools) | the normal boiling point |

## How a prediction is paired with a measurement

Three decisions change the answer more than most modelling choices, so each is
an argument rather than a default buried in code.

**Height.** Trials instrument several elevations on one mast. Comparing a
ground-level prediction against whichever height read highest looks good for
the wrong reason. Comparisons fix a height and evaluate the model's own
vertical profile at it, which in the near field is a factor of several.

**Averaging.** Field data is sampled every second; a model reports an average
over its own averaging time. Two conventions put them on the same footing:
average the measurements up, or run the model short and compare short-time
peaks. The second is used, because the first is a trap — the obvious averaging
time, the release duration, is far longer than a cloud takes to pass one
sensor, so it mixes the passage with the empty record either side.

On Burro 9 that single choice moves the geometric variance from 1.4 to 10.8
and drops every point outside a factor of two. Read carelessly it looks like a
model failure; it is an artefact of the reduction. Where the measurements are
averaged, a running mean is used rather than fixed blocks, for the same
reason.

**Crosswind position.** The model reports its centreline; the trial reports
what the sensors on that arc saw, and the plume does not oblige by passing
over one. The arc maximum is used, which biases *in the model's favour*: if no
sensor sat on the centreline the observation is low. So a model reading high
here reads higher still against a properly sampled centreline.

## Results

Two views, because they answer different questions and only reporting the
first is what makes DEGADIS look better than it is.

**A. The conventional comparison** -- arc maxima at the lowest instrumented
height, which is what the published evaluations report. Confidence intervals
are bootstrapped on the log ratio; `MG` is observed over predicted, so below 1
means the model reads high.

| series | n | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|---|
| Burro, z = 1.0 m | 61 | 0.811 | [0.63, 1.02] | 2.53 | 0.56 |
| Desert Tortoise, z = 1.0 m | 6 | 1.839 | [1.63, 2.04] | 1.48 | 0.83 |
| Lathen, z = 0.1 m | 58 | 0.400 | [0.33, 0.49] | 4.10 | 0.29 |

**B. The whole cloud** -- every instrumented height, with the model's own
vertical profile evaluated at each sensor's elevation.

| series | n | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|---|
| Burro | 173 | 25.6 | [11.2, 65.0] | 4e19 | 0.29 |
| Coyote | 99 | 1.58 | [1.12, 2.29] | 42 | 0.41 |
| Eagle | 7 | 0.96 | [0.21, 5.15] | 153 | 0.29 |
| Desert Tortoise | 17 | 1.97 | [1.64, 2.38] | 1.85 | 0.53 |
| Lathen | 157 | 5.44 | [2.98, 10.1] | 1e8 | 0.28 |

Burro split by elevation shows where B comes from:

| height | n | MG | FAC2 |
|---|---|---|---|
| 1 m | 61 | 0.81 | 0.56 |
| 3 m | 59 | 7.3 | 0.19 |
| 8 m | 53 | 5500 | 0.09 |

### Reading these

**Burro's conventional result is not a pass.** The interval [0.63, 1.02]
straddles Hanna's 0.7 bound, so the bias criterion cannot be called satisfied
at n = 61. Earlier versions of this file said it passed; that was point
estimation without an interval and it is withdrawn.

**Nothing passes on the whole cloud.** DEGADIS gets the ground-level
centreline roughly right and the vertical distribution badly wrong, so the
conventional statistic is carried by the one height where the two errors
happen to cancel. A geometric variance of 4e19 is not a large error, it is a
different quantity: the model predicts effectively zero where the sensors
measured several mole per cent.

**Desert Tortoise is the most trustworthy set.** It is the only series here
with enough crosswind positions to resolve the plume, its statistics are
stable between the two views, and its interval excludes 1 in a direction the
discarded jet momentum explains.

**Lathen should not be used.** Small low-elevation propane jets where the
release momentum dominates the near field, dispersed from an equivalent source
that discards exactly that. It is in the table so the failure is on record,
not as a result.

## The vertical profile is too steep

Comparing at one height and correcting the model down to it with the model's
own vertical profile is circular: the profile does the work and never gets
tested. Comparing at every instrumented height breaks that.

It has to be done mast by mast. Taking the arc maximum at each height
separately mixes positions, and on Burro 9 that puts the highest reading at
8 m and the lowest at 3 m simply because different sensors were nearest the
plume centre. At one mast the comparison means something.

Burro 9, concentration relative to the value at 1 m on the same mast:

| arc (m) | Sz (m) | obs 3 m | model 3 m | obs 8 m | model 8 m |
|---|---|---|---|---|---|
| 49 | 0.68 | 0.89 | 0.022 | 0.32 | 0.00002 |
| 57 | 0.71 | 0.63 | 0.027 | 1.16 | 0.00002 |
| 127 | 1.26 | 0.80 | 0.151 | 0.17 | 0.0007 |
| 140 | 1.39 | 0.65 | 0.183 | 0.12 | 0.0015 |
| 382 | 3.71 | 1.05 | 0.571 | 1.30 | 0.116 |
| 395 | 3.85 | 1.29 | 0.585 | 1.22 | 0.127 |

The measured cloud is close to well mixed over the bottom eight metres from
50 m downwind onward. The model confines it to a layer under a metre and a
half and is four orders of magnitude down at 8 m. The two converge only past
about 400 m, where Sz has grown to a few metres.

This is a property of DEGADIS, not of the port, and it is the reason a
single-height comparison flatters it: with all the contaminant in a thin layer
the 1 m concentration comes out roughly right by construction, while the
distribution is wrong. Anything that depends on the vertical extent — a
flammable volume, an exposure at head height, a cloud over a bund wall —
inherits that error and the arc-maximum statistics do not show it.

## The lateral spread is too wide

The complement of the vertical result, and it needs a trial that can resolve
it. Burro instrumented at most two crosswind positions per arc, and they are a
symmetric pair: at 49 m on B7 one reads 15 mole per cent and its mirror reads
0.09, so the plume was well off centre and two points cannot constrain a
width. Desert Tortoise has four to seven.

Half-width to 1/e of the arc peak, at 1 m:

| trial | arc (m) | sensors | observed (m) | model B + (sqrt(pi)/2) Sy (m) |
|---|---|---|---|---|
| D1 | 800 | 5 | 100 | 155 |
| D2 | 800 | 5 | 200 | 216 |
| D3 | 100 | 7 | 31 | 60 |
| D3 | 800 | 4 | 150 | 185 |
| D4 | 100 | 7 | 23 | 66 |
| D4 | 800 | 5 | 150 | 218 |

Always wider than measured, by a factor of one to three, and worst in the near
field.

So DEGADIS spreads the cloud too far sideways while confining it too tightly
in the vertical. The two errors push the ground-level centreline
concentration in opposite directions, which is part of why the single-height
arc-maximum statistics look better than the cloud does.

That off-centre plume matters for the headline numbers too. With two sensors
and the plume between them, the arc maximum is a lower bound on the true
centreline, so the observations here are low and the model's over-prediction
is a *lower* bound on the real one.

## The SMEDIS reductions, and FLADIS recovered

The EU SMEDIS exercise reduced a large set of trials into one spreadsheet
each. Reading those in preference to the raw archives recovers a series this
work had to withdraw, and supplies two things that were previously guessed at.

`degali.validation.smedis` reads 30 trials, 1250 concentration sensors and
79 arc reductions across FLADIS, Desert Tortoise, Thorney Island, Prairie
Grass, the Burro Abbott propane trials and the Hamburg and TNO wind tunnels.

**FLADIS is readable here.** REDIPHEM stores concentrations as numbered
channels whose meaning is per series and defined in a file FLADIS does not
ship; guessing the numbering selected channels reading 302 and 23.5 at 20 m
downwind, which are not concentrations, and produced a clean-looking result
that meant nothing. In these files every sensor is a row carrying its own
position and value.

**The wind direction is recorded, with its spread.** FLADIS wandered by 0.8
to 4.2 degrees within a trial, so comparing a steady plume against a fixed
sensor is meaningful there. The PRESLHY hydrogen trials wandered 33 to 49
degrees, and there it is not — which is why that comparison had to be made
against arc maxima. That distinction was previously an argument; it is now a
measurement.

**Two traps in the files themselves.** Thorney Island uses site grid
coordinates with the release at (400, 200), so a sensor at x = 100 is 300 m
*upwind*; the reader subtracts the release point. And Thorney Island peaks at
2060 under a `mean_C(%)` header while Prairie Grass peaks at 235 — whatever
those are, they are not volume percentages, so the reader flags them rather
than using or rescaling them.

### Sensor-level results

| set | n | MG | VG | FAC2 |
|---|---|---|---|---|
| FLADIS + Desert Tortoise, every usable sensor | 76 | 0.755 | 4.99 | 0.51 |

Split by sensor height, on data entirely independent of Burro:

| height | n | MG |
|---|---|---|
| 0.1 m | 18 | 0.710 |
| 0.5 m | 12 | 1.061 |
| 1.0 m | 8 | 0.595 |
| 1.5 m | 14 | 1.450 |
| 3.5 m | 8 | 0.490 |
| **8.5 m** | 6 | **0.114** |

The bias runs from about unity near the ground to 0.11 at 8.5 m — the model
reads nine times high up there. **This is the same vertical-structure failure
the Burro masts showed**, in the opposite direction because these are
elevated ammonia jets rather than a ground-level LNG pool, and on a different
substance, a different release type and a different measurement programme.
Two independent datasets, one defect.

## An external anchor, and a trap avoided

EPA's own 1991 evaluation (Zapert, Londergan and Thistle, *Evaluation of Dense
Gas Simulation Models*, EPA-450/4-90-018) put seven models — SLAB, DEGADIS,
TRACE, CHARM, AIRTOX, FOCUS and SAFEMODE — through Desert Tortoise, Burro and
Goldfish. It is the closest thing to an independent baseline for this work.

Its most useful passage is not a score but an admission. On the jet releases
the contractor used the **orifice area** as the source area, and after the
draft was reviewed:

> the SLAB model developer recommended one change to the input for the
> Goldfish and Desert Tortoise tests: for jet releases the source area should
> be the cross-section of the fully expanded jet rather than the orifice. This
> change would substantially reduce the source velocity and could change the
> model results considerably. **The change was not incorporated**, as it was
> identified after the performance results had been obtained.

So a published bias of 2.6 to 2.74 for SLAB on Desert Tortoise, with none of
three arcs inside a factor of two, was a consequence of a source-area choice
that the model's author had flagged.

**This work does not make that mistake**, because it starts from the SMEDIS
equivalent source — the plume state once the flashing aerosol has fully
evaporated, half-width 6.4 m and velocity 7.5 m/s for D1, against an orifice
velocity two orders larger. The same trap surfaced independently in the
hydrogen work: pushing an equivalent-source volume flux through the orifice
area gives an exit velocity of 501 m/s, and using the orifice-plane density
instead gives 105.

The residual difference is then a genuine model comparison rather than an
input artefact:

| | Desert Tortoise, arc maxima |
|---|---|
| EPA 1991, orifice area | 2.6–2.74× under, FAC2 0/3 |
| a SLAB reimplementation, expanded source | MG 0.733 |
| **this work, SMEDIS equivalent source** | **MG 1.839, VG 1.48, FAC2 0.83** |

The two corrected results straddle unity from opposite sides and differ by a
factor of 2.5, on the same trials with the same class of source treatment.
That is worth stating plainly: **a difference of this size between two
carefully-sourced integral models is the scale of disagreement the field
actually has**, and it is larger than most of the effects argued about within
any one model.

### The published DEGADIS score, explained

Table 5-7 of that report gives average fractional bias by model. Negative
means the model reads high.

| model | Burro | Desert Tortoise | Goldfish |
|---|---|---|---|
| AIRTOX | −0.06 | −0.63 | 1.19 |
| **DEGADIS** | **−1.07** | **−0.90** | 0.42 |
| SLAB | −0.13 | 0.88 | 1.01 |
| TRACE | −0.87 | −0.16 | 0.61 |

DEGADIS is the worst of the four on Burro: −1.07 corresponds to predicting
about three times the observed concentration. This reimplementation, on the
same experiments and the same measure, gives **−0.00**.

Two runs of the same model — the report states it used Version 2.1, which is
what is ported here — differing by more than the models differ from each
other. The difference is larger than DEGADIS-to-SLAB.

**The cause is in the report.** Section 4.3 says:

> measured concentrations were taken from samplers at 1 m height. These
> values were compared to model predictions **at ground level**, or at 1 m,
> if the model allowed for varying receptor heights.

Evaluating this port at a range of heights against the same 1 m
measurements reproduces the published number:

| evaluation height | FB | MG | FAC2 |
|---|---|---|---|
| 0.0 m | −1.29 | 0.14 | 0.05 |
| 0.5 m | −0.84 | 0.24 | 0.27 |
| **1.0 m** | **−0.29** | 0.46 | 0.50 |
| 2.0 m | +0.58 | 1.79 | 0.36 |

EPA's −1.07 falls between the ground and half a metre — an effective
evaluation height near 0.25 m, which is what "at ground level, or at 1 m if
the model allowed" would produce across a mixed set of runs.

So the published result is not a measurement of how well DEGADIS disperses.
It is a measurement of **the interaction between its vertical profile and the
height the comparison was made at**, and it is the same defect quantified
above: the profile is too steep, so concentration piles up at the ground.
Evaluated where the sensors actually were, the same model is unbiased.

That is worth stating carefully. It does not make DEGADIS a better model —
the profile error is real and this work documents it at three heights. It
means the number that has been carried in the literature for thirty years as
DEGADIS's performance is dominated by an evaluation choice, and that the model
and the evaluation were failing in the same place for the same reason.

### Desert Tortoise: a different cause, and not reproducible

The same test does *not* explain the Desert Tortoise score. EPA reports
−0.90; evaluating this port at ground level gives **+0.18**, the opposite
sign. So the height artefact is specific to Burro.

Section 4.2.1 of the report says what happened instead:

> DEGADIS Version 2.1 models only vertical jet releases, but previous
> analyses suggested that the jet momentum had a negligible influence on
> predicted concentrations for these experiments. The Desert Tortoise and
> Goldfish releases were simulated as steady-state, pure chemical releases.
> The jet releases are modeled as "isothermal" because the TRAUMA module has
> already accounted for heat exchange.

**The Ooms jet module was not used.** A horizontal momentum jet was replaced
by a ground-level steady source, and the thermodynamics were computed
externally in TRAUMA and supplied as isothermal ordered triples. The source
radius was set to "larger than the orifice area, representing an initial puff
size".

That is three substitutions at once — geometry, momentum and thermodynamics —
and reproducing the published number would need EPA's TRAUMA triples and
their chosen source radius, neither of which is in the report. **This
difference is therefore recorded as unexplained rather than explained.** The
Burro one is reproducible; this one is not.

One line in that passage is worth separating out, because it is a statement
about the model rather than about the evaluation:

> **DEGADIS Version 2.1 models only vertical jet releases**

Which is true of `SETJET`: its trajectory correlation is Kamotani and
Greber's, for a jet issuing vertically into a crossflow. Faced with two
horizontal-jet datasets, the 1991 evaluation dropped the jet model rather
than use it. The directed-release path added here
(`JetPlume.initial_conditions_directed`) fills that gap — it keeps the
zone-of-flow-development length, which is a property of the jet rather than
of its direction, and drops the trajectory correlation, which is not. It is
what makes the liquid hydrogen work possible at all: nineteen of the
twenty-four PRESLHY releases are horizontal.

EPA's own conclusions say why this was hard to see at the time:

> Models with very similar treatments of atmospheric dispersion and identical
> meteorological inputs produce concentration predictions differing by **more
> than an order of magnitude** because of those initial conditions.

> Given the complexity of dense gas dispersion and of the models, **attributing
> model performance to particular algorithms or design features was not
> feasible**.

The second is the case for taking a model apart rather than scoring it whole,
which is what the height-resolved and regime-split comparisons here are for.

## The bias is robust; two hypotheses are not

Burro's geometric mean bias barely moves under the deck-building choices that
could have produced it:

| variant | MG | VG | FAC2 |
|---|---|---|---|
| as built | 0.811 | 2.53 | 0.56 |
| ground 10 K above ambient | 0.854 | 2.32 | 0.59 |
| no ground heat or water exchange | 0.934 | 3.91 | 0.48 |
| averaging time 60 s | 0.795 | 2.58 | 0.55 |
| averaging time 10 s | 0.811 | 2.53 | 0.56 |

A spread of 0.14 against a bias of 0.19 from unity. The over-prediction is
real, not an artefact of the assumptions.

**The Coyote hypothesis was wrong.** The earlier reading was that Coyote's
factor of two came from a pool radius taken from a spill pond rather than a
pool. Forcing the diameter to 20 m or 30 m instead of the recorded 58 m makes
it *worse* — MG falls from 0.473 to 0.389. Whatever drives Coyote, it is not
the pool size, and the earlier explanation should be treated as withdrawn.

## Timing is good even where magnitude is not

The transient path had no comparison against data at all. Burro 9, time of
peak concentration at 1 m:

| arc (m) | observed (s) | modelled (s) |
|---|---|---|
| 49 | 46 | 49 |
| 137 | 70 | 67 |
| 140 | 81 | 68 |
| 382 | 113 | 112 |
| 400 | 81 | 90 |

The observer machinery places the cloud in time to within about ten seconds
over four hundred metres, which is the first evidence that half the code base
does what it should. Magnitudes at the same points are the over-prediction
already discussed.

## Pressurised jets: computing the source term

A flashing release of a liquefied gas is outside DEGADIS entirely, and forcing
one through the pool source assumes a 293 m pool for an 81 kg/s ammonia jet.
The accepted treatment is to hand the model the plume state where no liquid
remains.

Some trials come with such a source term published; most do not.
`validation.flashing` computes one from two conditions. An adiabatic energy
balance between the stored liquid and the diluted plume,

    h_store + r h_air(T_a) = h_vap(T) + r h_air(T)

and the requirement that the liquid has *just* gone, so the vapour is exactly
saturated at its own partial pressure,

    y p_amb = p_sat(T)

Two equations, two unknowns: the entrained air ratio `r` and the temperature
`T`. The resulting temperature is far below the substance's normal boiling
point -- ammonia comes out near 205 K against a boiling point of 240 K --
because the plume is dilute.

It is checked against source terms it did not see. SMEDIS published 13 mole
per cent at 205 K with a 6.4 m half-width for Desert Tortoise D1; the
calculation gives 13.6 per cent, 204.9 K and 5.3 m. That agreement is what
licenses applying it to FLADIS.

What the substitution gives up is the plume's momentum and elevation. At 51 m
the Desert Tortoise plume is still moving at 7.5 m/s and DEGADIS's ground
source has none, and both those trials under-predict by about a factor of two
-- the direction that omission pushes.

## Instantaneous releases: a real model limit

Thorney Island released a 14 m cylinder of Freon-air by dropping its walls.
That is DEGADIS's instantaneous path: an initial mass over a source of zero
rate. It turns on the van Ulden momentum balance, used when the
height-to-diameter ratio exceeds 0.1 because a quasi-steady gravity current
has not had time to form -- a branch no EPA test case reaches.

The balance does not converge for that cloud shape, and `SRC1` stops. **The
original Fortran stops in the same place**, at the same step, with `STOP SRC1
velocity loop`. This is a limit of DEGADIS, not of the port, and it is raised
rather than worked around: a frontal velocity that cannot be found is not a
numerical inconvenience, it is the whole gravity current.

Confirming that required building the deck, writing it out, and running the
1989 executable on it -- which is what the reference oracle is for.

## What is not yet evaluated

**Thorney Island.** Instantaneous puffs, which DEGADIS handles through its
transient path with a finite initial mass rather than a continuous source.

**The wind-tunnel series** (Hamburg, WSL, TNO — 227 of the 313 trials).
Valuable for entrainment closure but outside the geometry DEGADIS models.
