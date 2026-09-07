# Liquid hydrogen: what has been established

> **Superseded for the current corrected model (2026-09-03).** This document
> preserves the earlier investigation. Current source invariants and
> performance numbers are in `lh2-model-improvements-2026-09-03.md`.
> Its width-deficit discussion is additionally superseded by
> `gaussian-width-convention.md`: the fitted e-folding width was mislabeled as
> sigma, and the corrected width ratio is 1.033.

The DEGADIS work is in [validation.md](validation.md) and
[field-validation.md](field-validation.md). This is the hydrogen work: what
the model says, what the experiments say, and which of my own hypotheses the
experiments killed.

## The dense phase ends within a metre

On the adiabatic mixing line for a flashing LH2 release into ambient air at
15 °C:

| source | source rho/rho_a | dense above |
|---|---|---|
| pure vapour, no flash | 1.094 at 20 K | 99.9 mol % |
| flashing jet, 1 barg | 1.134 at 41 K | 85.2 mol % |
| flashing jet, 5 barg | 1.029 at 41 K | 90.2 mol % |

The flashing source matters: the liquid's latent heat has to come from the air
it entrains, so a real release is far colder at a given concentration than the
vapour line suggests, and the dense range is a hundred times wider than that
line implies.

It is still short. At the lower flammable limit of 4 mol % the cloud is
buoyant, at stoichiometric 29.5 mol % it is buoyant, and at the *upper*
flammable limit of 75 mol % it is buoyant. The dense phase exists only above
the UFL. Traced through the jet model, it is over within about half a metre of
the nozzle.

That is why the ground-level closure is swappable rather than fixed; see
[liftoff.md](liftoff.md).

## The jet model runs on all of them

Twenty-one of the 24 HSL trials — every one with a recorded flow rate — run
end to end with no failures, using the computed flashing source. For trial 10
(25.4 mm, 5 barg, 285 g/s, horizontal at 0.5 m):

| x (m) | plume z (m) | centreline (mol %) | rho/rho_a | Sz (m) |
|---|---|---|---|---|
| 0.5 | 0.50 | 80.7 | 0.81 | 0.06 |
| 3.1 | 0.52 | 36.0 | 0.85 | 0.40 |
| 7.1 | 0.71 | 21.3 | 0.91 | 0.78 |
| 14.0 | 1.49 | 9.5 | 0.96 | 1.37 |

Eighteen of the 21 give "no touchdown". That phrase is misleading and worth
stating precisely: the *centreline* leaves the ground, while the spread stays
at sensor height. At 14 m the centreline is 1.49 m up with Sz = 1.37 m, so a
sensor at 0.5 m sits 0.7 Sz below it and should read roughly half the
centreline — about 5 %, above the 4 % instrument ceiling. That is consistent
with the far-field measurements, which show 12 to 29 devices detecting and
saturating in nearly every trial.

## Two hypotheses I formed and the experiments killed

**Rainout as the source of ground-level detection.** The flash fraction is
only 6 % at 1 barg and 24 % at 5 barg, so most of the mass leaves as liquid,
and it seemed natural that droplets reaching the pad would make a secondary
ground-level source. The report is explicit:

> For un-impinged elevated releases, no evidence of rainout was observed
> during the releases.

Rainout appeared only *after* the valve closed and momentum collapsed, and in
low vertically-downward releases where pools formed. Nineteen of 24 trials are
un-impinged horizontal releases. So during the release there is no rainout —
which means the equivalent-source assumption that the flashed liquid fully
evaporates, the thing I had flagged as the weakest link, is supported by the
experiment rather than undermined by it.

**The near-field descent seen in the 2010 HSE trials.** Report RR986 records
the jet staying about a metre up and then moving toward the ground past 1.5 m.
That looked like evidence against the model's monotonic rise. But RR986 is a
60 l/min release — 71 g/s against E3.5's 84 to 285 g/s — and the same report
notes the cloud stayed on the ground above 5 m/s wind and became buoyant below
3 m/s. Lower momentum and different wind; not a comparison.

## What blocks a quantitative comparison, and what does not

The far-field Dräger devices range 0 to 4 vol %, and hydrogen's lower
flammable limit is 4 %. **The instrument saturates exactly at the
concentration a flammability assessment cares about.** Fifty-seven of 686
readings are at the ceiling and the maximum recorded value is exactly 4.00.
Those are censored observations, not measurements.

The far-field stand coordinates are in a figure, not a table, and two layouts
were used depending on wind direction. So 534 detections have no downwind
distance attached.

**The near-field array is different, and better.** It used Xensor XEN-5320
thermal conductivity sensors supplied by NREL, ranging 0 to 100 vol % — no
ceiling — and its coordinates *are* tabulated. Twenty-seven positions from
0.35 to 6.0 m downwind, at heights 0 to 1.0 m and offsets 0 and 1.0 m.

The concentrations themselves are plotted rather than tabulated, so they are
not machine-readable from the report. But the geometry is, which is enough to
state the prediction in advance.

## The prediction

`lh2_nearfield_prediction.csv` gives the model's concentration at each of the
27 near-field positions for all 14 horizontal trials with a recorded flow
rate: 378 points, with the plume centreline height and both dispersion
parameters alongside so the profile can be checked as well as the value.

Centreline, 0.5 m height, vol % hydrogen:

| trial | 0.35 m | 0.79 m | 1.78 m | 2.67 m | 4.0 m | 6.0 m |
|---|---|---|---|---|---|---|
| 10 (25.4 mm, 5 barg) | 86.7 | 71.1 | 47.6 | 38.6 | 30.6 | 23.3 |
| 11 (12 mm, 5 barg) | 72.1 | 46.6 | 27.5 | 21.2 | 17.5 | 14.5 |
| 12 (6 mm, 5 barg) | 50.9 | 29.5 | 15.9 | 11.7 | 9.3 | 7.5 |
| 16 (25.4 mm, 1 barg) | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 |

The last row is the sharpest test. Trials 16 and 17 release from 1.5 m rather
than 0.5 m, so the 0.5 m sensors sit a metre below a plume whose Sz is still
small, and the model predicts they see essentially nothing. If the
measurements show appreciable hydrogen there, the model's near-field vertical
structure is wrong in a way the other trials cannot reveal.

## Against the raw dataset

The report plots concentrations rather than tabulating them, and an earlier
version of this file compared against values read off three-dimensional
surface plots. The raw dataset (DOI 10.35097/1481) removes that: 24
workbooks, an uncensored `Xensor` sheet reading 0 to 100 vol %, and a
`Flowmeter` sheet at 1 Hz so the source is measured through the release. The
sensor serials in the workbook columns are the ones Table A3 gives positions
for, so the two join.

That comparison replaces the plot-read one entirely, and the numbers moved a
long way. Three things had to be got right first.

**The window.** Taking the first and last crossing of half the peak mass flow
gives, on trial 2, a window in which six per cent of the samples are flowing
and a mean mass flow of *minus* 0.5 g/s. The meter oscillates between −18.8
and +30.5 g/s there. The window is now the longest sustained run.

**The comparison quantity.** Point-to-point against a fixed sensor gives
MG 1.47 with **VG 79**. The wind direction wanders by 33 to 49 degrees within
every trial, so whether a given sensor was under the plume is what the
comparison measures. Taking the maximum over each arc instead — the reading
from whenever the plume crossed it — is what a steady model's centreline
corresponds to, and VG falls to 5.6.

**Applicability.** Even then the trials split into two populations:

| | trials | MG | between-trial spread (GSD) | range |
|---|---|---|---|---|
| momentum-driven | 9 | 0.738 | **1.16** | 0.61 – 0.95 |
| wind-steered | 8 | 0.387 | **2.94** | 0.019 – 1.38 |

The 1 barg releases leave the nozzle at 4 to 13 m/s into a 1.5 to 4 m/s wind,
so the wind decides where the plume goes. Nominally identical trials gave arc
maxima of 83 % and 4 %. **A steady jet model does not describe those**, and
`degali.lh2.assess` now warns below a velocity ratio of ten.

### Provenance

A reported statistic should be traceable to the published dataset without the
reader having to guess what happened to the rest. Every exclusion:

| | count |
|---|---|
| workbooks in the dataset (DOI 10.35097/1481) | 24 |
| excluded — not a horizontal release | −7 |
| **carried forward** | **17** |
| excluded — wind-steered, exit speed below ten times the wind | −8 |
| **trials in the reported statistic** | **9** |
| sensor readings on those trials | 197 |
| **reduced to arc maxima paired with a prediction** | **69** |

The wind-steered exclusion is the one that matters, and it is a statement
about where a steady jet model applies rather than a filter on the result: it
is applied by a criterion fixed in advance from the exit-velocity ratio, not
by looking at agreement. Those eight trials are reported separately below.

The reduction from 197 readings to 69 pairs is the arc-maximum step — one
value per distance per trial — not a rejection.

### The result

Momentum-driven releases, arc maxima against the model centreline:

**Superseded — recomputed from the reduction; the figures below were never computed by anything.** Recomputing gives MG 0.855, CI [0.672, 1.064], VG 1.30, FAC2 0.85 as shipped and MG 1.070 corrected, over 66 arcs, and the interval now includes 1. See `docs/lh2-recomputed.md`.
| n | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|
| **69** | **0.738** | **[0.591, 0.918]** | **1.41** | **0.83** |

The interval is bootstrapped over *trials* rather than over points, because
readings within one trial share a release, a wind and a source estimate and
are not independent. It excludes 1, so the model reads high by a real factor
of about 1.35, and it lies inside Hanna's bias band of 0.7 to 1.3.

Acceptance criteria throughout are Hanna, Chang and Strimaitis' — MG between
0.7 and 1.3, VG below 1.6, FAC2 above 0.5 — as used in the SMEDIS and LNG
model evaluation protocols. All three are met.

The residual has structure rather than being scatter:

| x (m) | observed | model | ratio |
|---|---|---|---|
| 0.35 | 84.2 | 96.7 | 1.15 |
| 1.19 | 78.9 | 82.8 | 1.05 |
| 2.67 | 44.8 | 64.1 | 1.43 |
| 6.00 | 22.1 | 32.5 | 1.47 |

That growth is not real. Regressing `ln(model/observed)` on `ln(x)` over all
69 pairs gives a slope of +0.018 with a bootstrap interval of [−0.168,
+0.171] — indistinguishable from zero. The table above is medians over a
seven-trial subset, and the trend in it does not survive the full set.

The residual is a **distance-independent factor of about 1.35**, and a
pre-registered set of tests (`prereg-entrainment.md`) rules out the source
rate, the stability class and the entrainment coefficient as causes. What
remains is the source normalisation or the measurement reduction, and the
array is not dense enough to distinguish them.

## The far field, out to 14 m

The far-field stands are fixed while the plume follows the wind, and the wind
direction is recorded only as a compass point every five minutes — 22.5
degree resolution against stands 10 to 12 degrees apart. So the plume cannot
be *placed* from the wind record.

It does not need to be. Five stands across an arc constrain a Gaussian, and
fitting one recovers the centreline concentration and the lateral spread from
the measurements themselves, with the plume's position as a nuisance
parameter. 65 arcs fit; 18 are well constrained, meaning the fitted centre
lands inside the stand pattern rather than on its edge.

| arc | fits | median centreline | median sigma_y |
|---|---|---|---|
| 10 m | 11 | 3.54 % | 3.82 m |
| 14 m | 7 | 1.77 % | 4.38 m |

**The bias reverses.** Against the model centreline at sensor height the far
field gives MG 3.6 — the model reads several times *low* — where the near
field gave 0.74, reading a third high. No single coefficient does both.

The trajectory does. For a 25.4 mm release at 5 barg the modelled plume climbs
from 1.4 m at 6 m downwind to 2.6 m at 10 m and 3.6 m at 14 m, while `Sz` is
only 1.0 to 1.4 m. Sensors at 0.5 to 2.5 m are inside the plume near the
source and underneath it further out: the predicted value at 1.5 m falls from
40 % at 6 m to 1.1 % at 14 m, while the fitted measurements are still 0.3 to
6.7 %. **The modelled plume rises too fast**, and that is the one defect that
explains both signs.

The fitted widths also give the meander directly. Measured `sigma_y` is 4.35 m
where the model has 1.14, and `sigma_y_total**2 = sigma_y_plume**2 + (x
sigma_theta)**2` puts `sigma_theta` near 17 degrees — consistent with the
compass record, and derived rather than assumed.

### Non-horizontal releases

One well-constrained arc fit each, so this is a direction rather than a
number. The ordering is what the physics requires:

| orientation | sigma_y | centreline |
|---|---|---|
| vertical down | 1.32 m | 6.69 % |
| horizontal | 4.38 m | 2.93 % |
| vertical up | 5.51 m | 1.56 % |

A downward jet impinges and stays compact; an upward one spreads and dilutes;
horizontal falls between.

Three of the five non-horizontal trials also leave at only three to five times
the wind speed, so the momentum criterion fixed from the horizontal trials
classifies them as wind-steered without being adjusted.

## What this campaign cannot test

**Transient behaviour.** See below.

**Transient behaviour.** The near-field array spans 0.35 to 6 m from a jet
leaving at 16 to 660 m/s, so the travel time is under a second everywhere
against a 0.3 s sampling interval. The measured arrivals have a median of a
few seconds, which is the right answer and an uninformative one — any model of
this release predicts the same. Their tails run to minus a hundred seconds,
which is hydrogen left over from an earlier release rather than physics.

The transient path is validated against Burro instead, where the cloud takes a
hundred seconds to cross four hundred metres and the modelled peak times land
within about ten seconds of the measured ones over that distance.

**The far field.** There are 686 readings from ten fixed stands at 9 to 14 m,
and they stay unused. Pairing a fixed stand with a model position needs the
wind direction through the release, and the dataset's only source for it is a
domestic weather station — the same sheet carries indoor temperature and
weekly rainfall — reporting a 16-point compass *string* every five minutes.

A release lasts two to four minutes, so every trial contains exactly one
direction sample, quantised to 22.5 degrees, and several read "---". At 14 m
that quantisation alone is ±2.7 m of lateral position against stands spaced
2.5 m apart. **The array cannot be placed relative to the plume even in
principle**, and no amount of work on the model changes that.

Earlier notes in this project recorded the far field as "wind direction not
recorded". It is recorded; it is just far too coarse, which is a different
statement and a checkable one.

**Non-horizontal releases.** Four of the five have usable flow. Three leave at
three to five times the wind speed and are wind-steered by the same criterion
that separated the horizontal trials — they fail against the model exactly as
that criterion predicts, so their failure is not new information. The fourth,
a downward jet at 276 times the wind, reaches the ground inside its own zone
of flow development and produces no comparable arcs.

The useful part is that the criterion was fixed from the horizontal trials and
then classified the vertical ones correctly without being adjusted.

## Known gaps

**Air condensation is not modelled.** The mixing line passes oxygen's boiling
point at 81 mol % hydrogen and nitrogen's at 84 mol %. Real LH2 clouds liquefy
air, oxygen preferentially — a hazard in its own right. It happens well inside
the buoyant range and does not change the buoyancy conclusion.

**The 1 barg releases are two-phase at the flow meter.** The report's own mass
flow analysis classifies all three 1 barg nozzle sizes as "gassy" at the
meter, while the 5 barg 6 mm and 12 mm releases are liquid. The
equivalent-source calculation starts from saturated liquid, which is right for
the 5 barg cases and wrong for the nine 1 barg ones.

**Ground heat transfer is extrapolated**, and the rest of the audit is in
[hydrogen-audit.md](hydrogen-audit.md).
