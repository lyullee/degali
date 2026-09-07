# Buoyancy: swapping the ground-level closure

An add-on, not a port. Everything in `degali.core` is DEGADIS 2.1 and is
validated against it; this implements published theory the original lacks, so
nothing here can be checked against the Fortran and it lives in
`degali.addons` for that reason.

## Sub-models are parts, not fixtures

DEGADIS makes one particular choice about what a ground-level cloud does with
its density difference: slump sideways while denser than air, do nothing once
lighter. That choice is right for LNG and wrong for hydrogen, so it is made
**swappable** rather than fixed.

A closure supplies the downwind model with a lateral spreading rate and,
optionally, extra state of its own. `DegadisClosure` reproduces the original
exactly and is the default — every parity test exercises it, and nothing
changes unless a caller asks for something else. `BuoyantClosure` adds the
vertical momentum balance.

```python
from degali.addons import make_buoyant
dw.closure = make_buoyant(dw)      # dense behaviour unchanged; buoyant clouds rise
```

The interface is narrow on purpose: a closure sees the local cloud state and
returns rates. It does not get to rewrite the mass or energy balances, which
are common to every variant and are validated against the Fortran.

## What the original does

DEGADIS's ground-level model has **no vertical momentum equation**. The
density difference enters `PSS` and `SSG` in exactly two places:

```fortran
DERY(iBEFF) = 0.D0
IF(delrho .GT. delrmn) DERY(iBEFF) = PRMT(9)*sqrt(delrho/rhoa)
...
gamma = (rho-rhoa)/cc
```

Lateral gravity spreading, guarded so it acts only while the cloud is
*denser* than air, and the density-concentration slope. Nothing else. A cloud
that becomes lighter than ambient simply stops spreading and disperses
passively, pinned to the ground for ever.

`JETPLU` **does** have buoyancy — its vertical momentum balance carries
`-RK1*gg*gamma*ccsysz`, which is why EX1's plume rises and never touches
down. So the gap is specific: buoyancy exists for an airborne jet and is
absent once the cloud is on the ground. A test asserts both halves of that
against the Fortran source.

## Why it matters for hydrogen

Cold hydrogen vapour is dense at 20 K — about 1.3 kg/m³ against air at 1.2 —
so an LH₂ release starts as a dense cloud and DEGADIS's dense phase is the
right description. But hydrogen's molecular weight is 2, and as the cloud
warms and dilutes it becomes violently buoyant: a 4 % hydrogen-air mixture at
ambient temperature is lighter than air.

So an LH₂ cloud passes *through* the regime DEGADIS handles into one it
cannot represent at all. A model that cannot lift the cloud off will keep a
flammable layer at head height that in reality has risen away.

The same is true of the release the theory was developed for — anhydrous
hydrogen fluoride, heavier than air on release and buoyant once it has reacted
with atmospheric moisture.

## The model

The integral model of Slawson and co-workers, extended to a ground-truncated
cross-section by AEA Technology under the EC URAHFREP project (report
AEAT/NOIL/27328006/001, June 2001, in the project material). Implemented from
that description; no code was taken from anywhere.

The cross-section is a **lozenge**: a rectangle of length `Ls` with
semicircular ends of radius `R`, truncated at ground level. One shape covers a
round source (`Ls = 0`), a line source (`Ls >> R`), and the wide shallow cloud
a dense phase leaves behind when it turns buoyant — which is the shape a
lift-off model actually has to start from, and the reason a round plume will
not do.

Six quantities along the plume axis:

| | |
|---|---|
| `m_g` | contaminant mass flux, conserved |
| `m` | total mass flux, grows by entrainment |
| `M_x`, `M_z` | momentum fluxes |
| `x`, `z` | trajectory |

with `dM_z/ds = g(rho_a - rho)A` the equation DEGADIS lacks, and entrainment
through the **free** perimeter only — the part not touching the ground:

    E = rho_a L_free [ alpha |u - u_a cos(theta)|
                     + beta |u_a sin(theta)| + gamma u_a ]

shear, cross-flow and ambient turbulence, with `alpha = 0.1`, `beta = 0.6`,
`gamma = 0.1`.

The report calls its own model "very crude", particularly in representing
ambient turbulence by `gamma u_a` and in neglecting velocity shear and ground
effects. It is implemented as stated rather than improved, because its virtue
is that it was compared against wind-tunnel data in exactly that form.

## Lift-off thresholds

Hall and Walker's wind-tunnel measurements, as the report reduces them, in
terms of a Richardson number formed on the buoyancy flux:

| `Ri*` | behaviour |
|---|---|
| ~2 | the concentration maximum leaves the ground |
| ~10 | ground-level concentration down to 10–20 % of the maximum |
| >70 | ground level below 5 % of the maximum |

Briggs' earlier estimates agree. Note the sign convention is opposite to the
dense-phase `Ri*` in `degali.core.atmosphere`: this one is positive for a
buoyant plume.

## Why hydrogen needs it, quantitatively

On the adiabatic mixing line for liquid hydrogen into ambient air at 15 °C:

| yc (mol %) | rho/rho_a | T (K) | |
|---|---|---|---|
| 100 | 1.094 | 20.4 | dense |
| **99.88** | **1.000** | **22.5** | the only crossing |
| 71.6 | 0.734 | 132 | lightest |
| 29.5 (stoichiometric) | 0.873 | 242 | buoyant |
| 4.0 (LFL) | 0.985 | 282 | buoyant |

The cloud is denser than air only above 99.9 mole per cent. **Every
concentration a flammability assessment cares about is buoyant**, and the
strongest buoyancy sits inside the flammable range.

LNG behaves the opposite way: at molecular weight 16 against air's 29 the
temperature effect wins while the cloud is cold, and Burro 9 is still dense at
5 mole per cent. Hydrogen at molecular weight 2 overwhelms it almost at once.

Two further hydrogen-specific points came out of the same calculation.

**The mixing table outgrows its array bound.** `IGEN = 42` is a Fortran
dimension, not physics. The adaptive thinning keeps whatever nodes linear
interpolation needs, and a line spanning 20 to 289 K — a factor of fourteen,
against LNG's under three — needs 46. The bound is raised rather than the
tolerance loosened, because loosening it would degrade every lookup silently.

**Air itself condenses.** The line passes 90 K, oxygen's boiling point, at
81 mole per cent hydrogen, and 77 K, nitrogen's, at 84 per cent. Real LH2
clouds do liquefy air, with oxygen condensing preferentially — a known hazard
in its own right. `Thermo` models water condensation but not air condensation.
It does not change the conclusion above, since it happens well inside the
buoyant range, and if anything it makes the gas phase lighter still.

## Two stages, and where one hands to the other

The closure and the standalone plume model divide the work, and the division
is forced by the physics rather than chosen for convenience.

The entrainment the downwind model supplies is calibrated for a cloud hugging
the ground: air enters through the top of a thin layer. A cloud that has left
the ground entrains through its **whole perimeter**, several times faster, and
that is what bounds a real plume's rise. Continuing to integrate with
ground-layer entrainment gives rise velocities above 20 m/s — an order of
magnitude faster than any buoyant plume goes — because the momentum is not
being diluted at anything like the right rate.

So `BuoyantClosure` decides *whether and where* the cloud leaves the ground,
and `LiftoffPlume`, which has perimeter entrainment, carries it after that.
On a worked hydrogen case the closure detects lift-off at 55 m with a rise of
0.4 m/s, and the plume model takes it to 36 m over the next 300 m, which is
the right order for a buoyant plume.

Carrying *momentum* rather than velocity is what makes even the first stage
bounded: the entrained air has no vertical momentum, so it dilutes what the
cloud has. Integrating `dw/dx` directly omits that term entirely.

## Coupling it to the dense phase

The intended use is at the point `PSS` gives up — where `delrho` falls below
`delrmn` and DEGADIS stops spreading the cloud. The dense phase supplies the
half-width (which becomes `Ls`), the concentration, the effective height and
the wind at that height; the lift-off model carries it from there.

Passing `AdiabaticTable.from_mass_fraction` as `density_of` couples the two to
the same mixing line, so a cloud becomes buoyant for thermodynamic reasons
rather than by assumption.

## First validation: the NASA spills

`docs/lh2-datasets.md` records why Witcofski and Chirivella's 1980 NASA
Langley trials are the right data for this. Table 4 of that paper tabulates
the **minimum height at which a flammable mixture was found** at the furthest
tower row, 33.8 m from a 9.1 m pond taking 5.7 m³ of LH₂:

| test | wind (m/s) | measured lowest flammable | model lower edge | error |
|---|---|---|---|---|
| 2 | 1.55 | 18.3 m | 15.6 m | −2.7 |
| 6 | 2.20 | 3.4 m | 6.4 m | +3.0 |
| 4 | 3.35 | 6.4 m | 2.3 m | −4.1 |
| 5 | 6.30 | 0.3 m | 0.0 m | −0.3 |

Regenerate with `test_lh2_assess_reproduces_the_nasa_regimes`, which computes
these rather than reading them. An earlier version of this table was stale by
up to 4 m against what the code produced, and quoted a mean error of −0.1 m
where the four values give −1.0.

The cloud is aloft in light wind and on the ground in strong wind, and the
transition is tabulated. That is the phenomenon this module exists for.

**The metric matters, and the obvious one is wrong.** Geometric mean bias and
factor-of-two counts are built for concentrations spanning decades. Applied to
a *height*, which is bounded below at zero, they mislead: test 5's model value
of 0.0 m against a measured 0.3 m is a thirty-centimetre disagreement that
both describe as "the cloud is on the ground", and a ratio calls it a total
failure. On absolute error the four tests give a mean of −1.0 m and an RMS of
2.9 m. The instruments are nine metres apart, so a 2.9 m RMS is smaller than the
spacing — but that is a statement about what the measurement could resolve,
not a claim that the model is accurate to 2.8 m. Two of the four individual
errors are 3.5 and 4.3 m, which on clouds ten to fifteen metres deep is a
quarter to a third of the depth.

Sorted into where the cloud is — grounded below 1 m, low between 1 and 10,
aloft above 10 — the model agrees with all four.

**What each measurement actually is** matters too. Test 2's 18.3 m is a
*sensor height*: flammable at 18.3 m and not at 9.1 m, so the true value lies
in [9.1, 18.3] and the model's 15.6 m is inside it. The other three come from
an adiabatic mixing model applied to thermocouple data rather than from a
direct concentration measurement, and carry whatever error that inversion has.

**The regime is right.** The model lifts the cloud clear of the ground for
test 2 and keeps it grounded for test 5 — the two extremes — and the rank
correlation across all four is 0.80. The Richardson criterion alone orders
them the same way: `Ri* = 110, 44, 13, 2.7` against measured heights of 18.3,
3.4, 6.4 and 0.3 m.

**The magnitude is not.** MG 1.71, FAC2 0.25 on four points, and with n = 4
neither number carries an interval worth quoting. The model places the lower
edge within a factor of two for the two extreme tests and is out by two to
three for the middle pair. Three of the four measured values are not
measurements of height at all but inversions of thermocouple data through a
mixing model, which carries its own error.

Three things had to be got right to reach even this, and each was wrong first:

* **The dense phase cannot be left out.** A parameterised density curve that
  never exceeds ambient makes the plume buoyant from the source and it
  rockets to 97 m. Feeding the real adiabatic mixing line, which is *denser*
  than air above 85 mole per cent, halves that.
* **The comparison quantity is the lower edge, not the centroid.** Table 4
  reports the bottom of the flammable region. Comparing it against the plume
  centre gives MG 0.17; comparing it against the cross-section's lower edge
  gives 1.71.
* **The concentration floor must be the flammable limit.** Integrating past
  it lets the plume keep rising on buoyancy it no longer has anywhere that
  matters.

### Concentration, against the grab bottles

Table 3 of the same paper gives the maximum concentration at the 33.8 m tower
row from the evacuated sample bottles — the uncensored measurement, the one
the authors say was "quite accurate and the preferred source".

| test | 9.4 m | 18.6 m | model |
|---|---|---|---|
| 2 | 0 | 0 | 0 / 1.4 |
| 4 | 4.2 | 0.5 | 3.1 / 3.1 |
| 5 | 29.2 | 0 | 5.4 / 0 |
| 6 | 18.7 | 19.0 | 2.2 / 2.2 |

**Which heights are in the cloud is mostly right** — 5 of 8 points agree on
whether a flammable mixture is present, and the two clean structural cases are
correct: test 5 has the cloud at 9.4 m and nothing at 18.6, test 2 has nothing
at either.

**The concentration is 2.5x low** (MG 2.45, FAC2 0.20 on the five points where
both are non-zero), and the reason is structural rather than a coefficient.
The lozenge is a *top-hat*: uniform across a cross-section whose radius has
grown to 8-26 m by the 33.8 m tower row. A grab bottle catches the peak of a
real, peaked profile. Comparing a cross-sectional mean against a peak
under-predicts by roughly the profile's peak-to-mean ratio, which for a
Gaussian is a factor of two to three.

### Giving the section a profile

So the section was given one. Integrating a Gaussian over the lozenge — a
one-dimensional profile across the flat part, two-dimensional over the
semicircular ends, cut at the same 2.15 sigma DEGADIS uses — gives a
peak-to-mean ratio of **2.57** for a circular section and 2.2 to 2.3 for the
lozenge shapes these plumes actually take.

The measured discrepancy was 2.45. That agreement is what identifies the
shape as the cause rather than any coefficient, and the correction is
arithmetic: no parameter was fitted.

**It improved the bias and made the scatter worse**: MG from 2.45 to 2.05, VG
from 21.8 to 63.8, and the flammable/not agreement from 5 of 8 to 4 of 8.

That is not a failure of the profile; it is the profile doing its job. A
top-hat section is insensitive to where its centre is — a point either falls
inside or it does not. A peaked one is very sensitive, so it converts the
*height* error, which is a known factor of two to six, into a concentration
error. Before, the two errors partly cancelled.

Test 6 shows it plainly: the bottles read 18.7 % at 9.4 m and 19.0 % at
18.6 m, nearly uniform over nine metres, while the model puts the centre at
26.6 m and so suppresses the lower level. The real cloud was broader and
lower than the model has it.

So the structural fix is correct and cannot pay off until the trajectory is.
The profile is kept, because it is right and because the next comparison needs
it to be meaningful.

### Two structural improvements, and what they changed

The model is a 2001 integral formulation and there is no reason to keep it as
it was written. Two things were wrong with it as an approximation, independent
of any data:

**Buoyancy used the density at the mean concentration.** The force is
`g(rho_a - rho_bar)A` with `rho_bar` the *mean density*, and on a curve as
bent as hydrogen's mixing line that is not the density at the mean
concentration. `_section_density` now integrates the mixing line over the
Gaussian profile — one-dimensionally across the flat part of the lozenge,
two-dimensionally over the ends.

This matters in principle for the reason Giannissi and co-workers give:
condensation and freezing of atmospheric humidity dominates LH₂ cloud
buoyancy, and it happens in a *spatially extended* shell where cold hydrogen
meets moist air. Evaluating the thermodynamics once at the section mean cannot
see that shell; integrating over the profile can, because the curvature the
condensation puts into `rho(c)` is what the averaging picks up.

**There was no form drag.** Mack and co-workers report that integral models
over-predict the rise of strongly buoyant plumes and that adding pressure drag
fixes the trajectory. `JETPLU` already applies exactly that to its own
inclined cross-section, so the same term was added here.

**Neither changed the answer here.** Drag moves the rise by 2 %, because these
plumes reach only about a metre per second and the resistance is quadratic.
Profile integration moves it by 0.1 m, because by the far field the plume is
dilute and `rho(c)` is nearly linear there — the curvature is at high
concentration, near the source, where the cross-section is small.

Both are kept. They are correct where the previous form was an approximation,
and they will matter for a release that rises faster or is compared closer in.

### A known deficiency of this model, inherited

The URAHFREP report that this implementation follows says of itself, in
comparing against Hall and Walker's wind tunnel data:

> the model substantially under-predicts the suppression of plume rise for
> wide sources

That is the direction of the residual here. Witcofski's pond is 9.1 m across
against a plume that rises tens of metres, which is a wide source, and the
model puts the cloud higher than measured on three of four tests. **The
over-prediction is a documented property of the formulation rather than
something introduced in porting it.**

### Added mass: tested elsewhere, and it fails

The obvious remedy is an added-mass term, `k = C_A (rho_a/rho)(B/h)`, which
predicts exactly the observed behaviour — a wider source rises less at the
same buoyancy flux. It has been pre-registered and tested against the Hall
and Walker geometry (six source widths spanning a factor of 64, six values of
`C_A`) in a parallel reimplementation of SLAB, and the result was:

| `C_A` | spread across widths | note |
|---|---|---|
| 0 | 3.70 | and the *wrong way round* |
| 0.20 | 1.46 | |
| 0.6366 (2/π) | 1.46 | still grounded at the threshold |

The spread saturates at 1.46 and **no value of `C_A` fits**. Tickle, who
develops DRIFT, reports the same from the same data: including added mass
made agreement with the wind tunnel *worse* and suppressed the rise too much.

So this is a hypothesis that has been formulated, pre-registered, tested and
falsified. It is recorded here so that it is not tried again.

### What was ruled out

Three candidate causes for the residual were tested and eliminated, which is
worth more than n = 4 of curve-fitting would have been:

* **Source geometry.** Setting the initial cross-section from the 9.1 m pond
  and a one to three metre cloud depth, rather than imposing a pure-hydrogen
  source and letting continuity give the area, moves the bias from 1.57 to
  1.66. Not the cause.
* **Buoyancy generation.** Giannissi and co-workers' CFD study of the same
  HSL trials finds that condensation and freezing of *atmospheric humidity*
  is the main contributor to LH₂ cloud buoyancy, with air condensation
  confined near the release. That mechanism is already in the mixing line
  here: raising relative humidity from 0 to 80 % strengthens the buoyancy by
  11 %. It also pushes the wrong way — more buoyancy means a higher plume,
  and the plume was already high.
* **Input uncertainty.** Propagating the wind ranges Table 1 gives and a
  spill quantity between 4.5 and 5.7 m³ produces a band that contains the
  measured value for test 2 and not for the others.

Tests 4 and 6 err in opposite directions — and the reason they cannot be
reconciled by any of this is that **the measurements are not monotonic in wind
speed and the model necessarily is**. Test 6 at 2.2 m/s has its cloud at
3.4 m; test 4 at 3.35 m/s has it at 6.4 m. More wind, higher cloud. The two
differ in humidity, 29 % against 43 %, in the direction humidity would push —
but the mixing line only moves 1.3 % in density between those, which cannot
produce a factor of two in height.

Set against the depth of the clouds, the residual is small. Table 3 shows test
6 flammable at 9.4 m *and* 18.6 m with a lowest edge at 3.4 m: a cloud at
least fifteen metres deep. The model's offsets are −2.7, −4.1, −0.3 and
+2.9 m on clouds of that size, with instruments spaced nine metres apart, and
three of the four measured values obtained by inverting thermocouple data
through a mixing model rather than measured directly. An earlier version of this file proposed tying the ambient-turbulence
entrainment coefficient to the friction velocity; that proposal is withdrawn,
because it was an attempt to fix something the metric had misrepresented.

So: the regime prediction agrees on four of four and the section shape is
confirmed by an independent route, but the *height* is a grade C result -- a
direction, on four points, three of them inversions of thermocouple data
rather than measurements of height. It should be written as "consistent
with", never asserted. No coefficient was fitted; the three entrainment
constants are the report's own. See [claim-grading.md](claim-grading.md).

## What is not claimed

This has been checked for internal consistency — contaminant flux is
conserved, buoyant plumes rise, neutral ones do not, dense ones sink, and the
lozenge reduces correctly to a circle and to a line — and against the
thresholds the report states. It has **not** been validated against
measurements, because the wind-tunnel data it was fitted to is not in hand:
the report says Hall and Walker's data was planned for REDIPHEM, and the WSL
series there is Richardson-number-parameterised dense puffs, not the buoyant
plume work.

Until that data is found, this is a documented implementation of a published
model, not a validated one, and results from it should say so.
