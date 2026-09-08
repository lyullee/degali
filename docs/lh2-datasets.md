# Other liquid hydrogen datasets in the project material

The PRESLHY/HSL E3.5 trials are not the only LH₂ dispersion data here, and
they are not the best-suited. This records what else is present, because the
E3.5 comparison is limited by measurement censoring and a fourteen-metre
fetch, and at least one of these is free of both.

## NASA Langley — Witcofski and Chirivella (1984)

*Int. J. Hydrogen Energy* **9**(5) 425–435. Seven unignited LH₂ spills, 1980.
The classic LH₂ dispersion dataset and still the most cited.

**Why it matters more than E3.5 for this work.**

*Measurement censoring is solved.* The paper is explicit that the catalytic
hydrogen sensors are "questionable" above 4 % and "meaningless" above 8 % —
the same lower-flammable-limit ceiling that caps the PRESLHY Dräger devices.
NASA got round it with **evacuated 500 ml grab bottles** analysed afterwards
on a gas chromatograph, and says the bottle data "were quite accurate and were
the preferred source". Table 2 records bottle concentrations up to 45.6 vol %.

*The vertical range is an order of magnitude larger.* Bottle clusters sit at
**1, 9.4 and 18.6 m**, sensors at 0.9, 9.1 and 18.3 m. PRESLHY samples 0.5 to
2.5 m. Testing a model's vertical structure — the thing the Burro work showed
that conventional statistics hide — needs the taller array.

*It reaches far enough to matter.* Nine towers, 19.5 m tall, on rings out to
36.6 m, with tower 5 at 18.3 m. PRESLHY's furthest centreline sensor is 14 m.

*It shows lift-off directly.* Table 3 gives maximum concentrations at the
furthest towers; Test 6 recorded **19.0 % at 18.6 m and 18.7 % at 9.4 m**
while the ground-level reading was 2.7 %. Table 4 tabulates the *minimum*
height at which a flammable mixture was found — 18.3 m for Test 2, 6.4 m for
Test 4, 3.4 m for Test 6, and 0.3 m for Test 5, the windiest.

That last table is close to an ideal test of a lift-off model: the cloud is
aloft in low wind and on the ground in high wind, and the transition is
tabulated. It is exactly the phenomenon `degali.addons.liftoff` implements
and cannot currently be checked against.

**Conditions**, from Table 1 (spill time s, wind m/s, air °C, RH %):

| test | date | spill time | wind | T | RH |
|---|---|---|---|---|---|
| 1 | 1 Aug | 60 | 2.7–3.1 | 30 | 18 |
| 2 | 25 Sept | 40 | 1.3–1.8 | 24 | 49 |
| 3 | 10 Oct | 85 | 4.5 | 26 | 27 |
| 4 | 22 Oct | 33 | 3.1–3.6 | 15 | 43 |
| 5 | 24 Nov | 24 | 6.3 | 12 | 43 |
| 6 | 18 Dec | 35 | 2.2 | 15 | 29 |
| 7 | 18 Dec | 240 | 3.1 | 17 | 29 |

Tests 2, 4, 5 and 6 are the ones the authors treat as primary: short spill
times and therefore the highest rates.

**What still has to be established** before it can be used: the spill
quantity, the pond geometry, and the surface roughness. The narrative gives
tower positions and instrument heights; the remaining conditions are in the
paper's own text and in the NASA reports it cites.

## HSE 2010 campaign — RR985, RR986, RR987

`rr986` covers the unignited releases. Sensor positions **are** tabulated: five
mounts at 1.5, 3.0, 4.5, 6.0 and 7.5 m, at six heights from 0.25 to 2.75 m,
on the same 32 m concrete pad HSL later used for PRESLHY. Release rate 60
l/min, about 71 g/s — a quarter of the E3.5 rates, so momentum is much
smaller.

Concentrations are in figures rather than tables, so the same digitising
problem applies. But the report's *narrative* observations are usable as
qualitative tests, and two of them bear directly on the model:

- the jet stays about a metre up near the release and moves toward the ground
  past 1.5 m, then rises again when the valve closes;
- the cloud stayed near the ground above 5 m/s wind and became buoyant below
  3 m/s.

The second is a wind-speed threshold for lift-off and can be compared against
the Richardson-number criterion in `addons.liftoff` without any digitising.

## PRESLHY E3.4 — pool spills, DOI 10.35097/1319

Ten unignited spills onto gravel, sand, concrete and water. Concentrations at
**35, 45 and 55 cm only** — a 20 cm vertical span and no downwind distance —
so this is a source-term dataset, not a dispersion one. Substrate is the
independent variable, which DEGADIS represents through a single ground
heat-transfer coefficient. Worth using to test that coefficient; useless for
dispersion.

## ELVHYS WP4.2 — transfer connection space, DOI 10.18710/JXJP0H

The official archive contains 48 small-leak dispersion, ignition and
pressure-peaking tests in a 1 m3 stainless-steel enclosure. Tests 10 and 11
are nominal repeats: horizontal 1.0 mm release, about 2 barg and active
ventilation near 500 L/min. Unlike the earlier field sets, it supplies raw
20 Hz local H2, temperature, nozzle-pressure, fan-flow and ambient time series.

This is valuable but not a free-jet dataset. Walls, floor, ventilation,
finite-volume accumulation and vent exchange act on the observations. The
pre-registered Test-10/Test-11 screen found no measured hydrogen mass flow and
materially different actual nozzle pressure/temperature histories. The
archive dictionary now resolves the apparent 200-versus-250 mm conflict: the
horizontal nozzle is at 250 mm and the bottom samplers are at 200 mm. D4.6's
200 mm nozzle entry is retained as a report discrepancy. No quantitative
free-jet model score is defensible without the hydrogen source flow. See
[`elvhys-tcs-audit.md`](elvhys-tcs-audit.md) and the exact request in
[`data-request-elvhys-source.md`](data-request-elvhys-source.md).

## What this changes

The E3.5 comparison in [lh2-results.md](lh2-results.md) rests on five trials,
one height, values read from surface plots, and a concentration range under
one decade. Its reliability is correspondingly limited, and no amount of work
on that dataset alone will fix it.

Witcofski is the dataset that can: uncensored measurements, a vertical range
from 1 to 18.6 m, distances to 36.6 m, and a tabulated lift-off transition
against wind speed. Bringing it in is the highest-value next step for the
hydrogen work, and it needs no data that is not already in this folder.
