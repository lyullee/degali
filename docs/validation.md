# Validation record

Every module of degali is compared against DEGADIS 2.1 running in the same
repository, at full double precision, through purpose-built probes. This file
records what was compared and how closely it agreed. The suite enforces it.

**Summary.** All five EPA test cases run end to end in Python. Every module
below the integration drivers matches the Fortran to round-off; the drivers,
where adaptive step control and loose internal tolerances make bit-exactness
unattainable, match to between 1e-6 and 3 per cent depending on the quantity,
with each residual traced to a specific cause.


Generated against DEGADIS 2.1 (dated 12260) as distributed by EPA SCRAM.

## Layer 1 -- the oracle reproduces EPA's golden listings

The Fortran in `reference/fortran` is patched only for portability. Its output
is compared against the `.lis` files shipped with the EPA test cases.

| case | description | lines | non-timestamp mismatches |
|---|---|---|---|
| `b9`  | steady state, ground-level plume (Burro 9) | 190 | 0 |
| `b9t` | transient, ground-level plume (Burro 9) | 1086 | 8 |
| `ex1` | vertical jet, no touchdown (MIC) | 133 | 0 |
| `ex2` | vertical jet with touchdown (ammonia) | 398 | 0 |
| `ex3` | vertical jet with touchdown, simplified density | 470 | 0 |

The eight `b9t` lines differ in the seventh significant figure only, e.g.
`9.803300E-05` against `9.803299E-05`. The EPA readme states this is expected
between compilers.

## Layer 2 -- the Python port reproduces the oracle

Full double-precision state, extracted by `reference/fortran/PROBE.for` and
`PROBE2.for`. Relative deviation, Burro 9 conditions.

### Scalars

| quantity | Fortran | max relative deviation |
|---|---|---|
| `humid` (absolute humidity) | 5.36121525139709836e-03 | 0 |
| `rhoa` (ambient density) | 1.07181180843550306e+00 | 0 |
| `ustar` (friction velocity) | 2.18777676015542516e-01 | 4.7e-13 |
| `alpha` (wind power law, legacy path) | 1.06323781023970457e-01 | 3.3e-13 |
| `gammaf` = Gamma(1/(1+alpha)) | 1.06550696789987609e+00 | 2.1e-16 |
| `hmrte` (source enthalpy) | -4.09207317488309636e+05 | 0 |
| `hwrte` (surface water enthalpy) | 2.70424999999997863e+03 | 0 |

### Adiabatic mixing table (`SETDEN`)

28 nodes, 5 columns. Node count and node positions both match.

| column | max relative deviation |
|---|---|
| `yc` mole fraction | 0 |
| `cc` concentration | 6.9e-16 |
| `rho` density | 6.7e-16 |
| `h` enthalpy | 0 |
| `T` temperature | 5.2e-16 |

### Closures

| routine | points | max relative deviation |
|---|---|---|
| `PHIF` (5 prescriptions x 9 Ri x 4 Ri_t) | 180 | 0 |
| `PHIHAT` (5 densities x 5 fetches) | 25 | 3.2e-15 |
| `SURFAC` (4 prescriptions x 4 T x 3 h, both fluxes) | 96 | 0 |

### Input deck (`IO.FOR`)

Burro 9 deck, including the columns `IO.FOR` derives rather than reads
(stability defaults, humidity reconciliation, source enthalpy and density).

| quantity | max relative deviation |
|---|---|
| 18 scalars | 0 |
| source table, 4 rows x 8 columns | 0 |

### Secondary source blanket (`SRC1`)

Derivative vector and all diagnostics, with the Fortran's own `alpha` and
`u*` injected so the blanket algebra is isolated from the fit chain.

| branch | states | max relative deviation |
|---|---|---|
| gravity slumping (6 times x 5 mass/radius scalings) | 30 | 2.5e-13 |
| van Ulden momentum balance (6 x 5 x 3 momenta) | 90 | 1.2e-13 |

Refitting `alpha` in Python instead of injecting it raises the worst
derivative deviation to 1.1e-10. That is not a porting error: `D_mc` is
`erate - qstrll`, a difference of two numbers near 130 kg/s whose result is
order 1, so the 3e-13 residual in `alpha` is amplified by about 10^3. The
residual itself comes from `ZBRENT`'s hard-wired `EPS = 3e-8` in `ALPH`.

The momentum branch is unreachable from any EPA test case: all five set
`GMASS0 = 0`, so the height-to-diameter ratio is zero and `VUFLAG` is never
set. It is nonetheless the most intricate part of `SRC1`, so the probe forces
it on and sweeps it.

### Numerical parameter file (`.ER1`)

All 35 fields, read positionally with `FORMAT(10X, G10.4)`. Nine are
cross-checked against the Fortran commons; deviation 0.

### DEG1 driver, end to end (`SRC1O`, `RKGST`, `CRFG`)

| check | size | max relative deviation |
|---|---|---|
| RKGST step sequence (x, IHLF, 6 states, 3 diagnostics) | 200 steps | 6e-10 |
| printed source parameter table, after both thinning stages | 8 x 9 | 4.4e-6 |
| `/GEN3/` handoff vectors in `.TR2` | 8 x 7 | 2.9e-6 |
| `RM`, `SZM`, `EMAX`, `RMAX`, `TSC1`, `ALEPH`, `TEND` | 7 | < 1e-6 |
| secondary source `outcc`, `outsz`, `outb`, `outl`, `swcl`, `swal`, `senl`, `srhl` | 8 | < 1e-6 |

The listing and `.TR2` figures are at their print-precision floors (six and
seven significant figures respectively).

The step-sequence comparison deserves a note. RKGST's controller decides where
output is recorded, so reproducing the trajectory is not sufficient -- the
sequence of steps has to match too. It does, exactly, through step 200, which
covers the whole of B9's blanket phase (it terminates at t = 7.66 s, around
step 170). Past that the two histories separate, because the controller
compares a computed error against a fixed bound and round-off eventually puts
the two sides of that comparison on opposite sides of it. The solutions stay
close -- the radius still agrees to 1e-5 at step 1000 -- but the recorded
points no longer line up, so a step-by-step comparison stops being meaningful
there. This is a property of adaptive step control, not of the port.

### SZF: the blanket-free source layer

| check | size | max relative deviation |
|---|---|---|
| sigma_z, cclay, wclay, rholay over 4 fluxes x 4 lengths | 64 | 4.4e-15 |

### Downwind dispersion (`PSS`, `SSG`)

Derivative vector and all diagnostics, configured from the `.TR2` handoff
exactly as `DEG2S` configures itself.

| stage | states | max relative deviation |
|---|---|---|
| `PSS` dense phase (4 distances x 4 mass fluxes x 3 heat levels) | 48 | 9.0e-7 |
| `SSG` Gaussian phase (same grid, virtual origin 25 m) | 48 | 5.0e-6 |
| mixing table `DEG2S` rebuilds on the secondary-source composition | 28 x 5 | 3.2e-10 |

Every base quantity — concentration, density, temperature, sigma_z, layer
density — agrees to 3e-10. The larger figures are confined to the
flammable-mass derivatives, which depend on `gamma = (rho - rho_a)/cc`; with
`rho` near 1.09 and `rho_a` at 1.0718 that subtraction amplifies round-off by
about 300x, and `gamma` then enters a logarithm and an incomplete gamma
function.

Getting there required using the handoff's `deltay` rather than recomputing
it. `DEG2` reads it from `.TR2`, which carries seven significant figures, and
sigma_y is directly proportional to it; recomputing at full precision moved
the concentration by 5e-7 — larger than everything else in the routine
combined.

## Two bugs the parity testing caught

**`ERBND` is 0.0025, not 0.005.** The driver was initially run on guessed
defaults instead of reading `.ER1`. The trajectory then diverged from the
Fortran at output step 3. Working backwards from the ratio of the two step
sizes gave 2^(-1/4), which pinned the discrepancy to the exponent in RKGST's
step expansion, `0.9 h (tol/delt)^0.25`, and from there to `tol` being wrong by
a factor of two. That prompted writing the `.ER1` reader instead of guessing.

**`PRMT(21)` is `ERTE`, the primary source rate.** The steady-state exit at
label 122 reuses the local `Qstar`, recomputing it as `PRMT(21)/(pi R^2)`, and
then falls through to the terminal write -- so the last row of the printed
table carries the recomputed pair rather than the integrated take-up. Reading
`PRMT(21)` as the take-up instead left that one row 3.3e-4 out while all 70
other values in the table were exact, which is what made it findable.

### Steady-state downwind driver (`DEG2S`)

| check | size | max relative deviation |
|---|---|---|
| initial conditions after the material balance closes | 6 | 1.7e-10 |
| `.ER2` parameter file | 26 | 0 |
| dense-phase step sequence: distance, bisection count | 2 steps | 1e-12 |
| dense-phase step sequence: rho·u·H, sy², B_eff | 2 steps | 1e-6 |
| profile rows recorded | 35 vs 35 | — |
| mole fraction, interpolated to reference distances | 35 | 1.1e-2 |
| concentration | 35 | 2.5e-2 |
| sigma_z | 35 | 5e-3 |
| mass above the LFL | 1 | 1.3e-3 |
| mass between UFL and LFL | 1 | 3.3e-3 |

### The one place accuracy is lost

`ADDHEAT` resolves the heated temperature with `ZBRENT` at `acrit = 1e-3` K,
then forms the heat capacity as `dh/(temp − amt)`. Early in the dense phase
`temp − amt` is only a few hundredths of a kelvin, so the quotient inherits a
percent-level uncertainty; it enters the forced-convection term of `SURFAC`,
which sets the ground heat flux, which integrates back into `dh`.

At the first output point of the dense phase the distance and the bisection
count agree to 1e-12 and the other three integrated states to 1e-6, while the
added-heat state differs by 1.6 %. That asymmetry is what localises it. The
consequence is that output points land at slightly different distances, so
profile comparisons interpolate rather than align rows.

Tightening `acrit` would remove it but would no longer be DEGADIS. The
non-legacy path should instead reformulate the heat capacity so it is not a
quotient of two small differences.

### Transient observers (B9T)

DEGADIS handles an unsteady release by switching to a Lagrangian frame: it
releases `NOBS` imaginary observers over the source at staggered times and
follows each one downwind, so every observer becomes its own quasi-steady
plume.

| check | size | max relative deviation |
|---|---|---|
| release window `T01`, `DTOB`, `XEND` | 3 | 0 |
| kinematics and edge crossings, 30 observers x 7 quantities | 210 | **0** |
| `OB` derivatives and layer diagnostics (5 times x 3 widths x 4 masses) | 660 | 1.5e-15 |

One thing had to be reproduced to make the probe run at all. An observer that
cannot reach an edge -- released too late to cross the downwind edge before
the source ends, or too early to meet the upwind edge -- is not given a root
to solve; `SSSUP` hands it the endpoint. Calling `TUPF` unconditionally makes
the Fortran abort inside `LIMIT`, which is how the guards were found.

### Transient supervisor (`SSSUP`)

For each of the 30 observers: integrate `OB` across the source, then convert
what it collected into a starting plume.

| check | size | max relative deviation |
|---|---|---|
| five integrated observer states | 150 | 7e-10 |
| layer composition preserved by `OBOUT` | 150 | 7e-10 |
| derived source strength, half-width, concentration, sigma_z0 | 120 | 7e-10 |

All 30 observers then run `PSS` and `SSG` to completion, producing 2291
profile points. Each carries an arrival time as well as a distance,
`TS(t0, x) = t0 + (x + Rmax)^(1/(1+alpha))/aleph`, which is what turns a set of
steady profiles into a transient field.

Comparing against B9T's own listing needs `DEG3`, which sorts those points by
distance and interpolates across observers to produce cloud footprints at
fixed times. That is the next port.

**One mistake worth recording.** `RHOE` lives in `/parm/` and is read once
from the handoff; it is not the current mixing table's last row. `SSSUP`
rebuilds that table for every observer, so reading the clamp from the table
would let it drift -- the second observer clamped against the first
observer's extrapolated centreline instead of the release density. The
symptom was a uniform 13 % error in sigma_z0 while every other seeded quantity
was exact, which localised it immediately.

### Cloud snapshots (`DEG3`)

| check | size | max relative deviation |
|---|---|---|
| snapshot instants chosen by `GETTS` | 19 | exact |
| points per snapshot | 19 | **exact** |
| cloud extent | 19 | 3e-3 |
| peak mole fraction, t <= 100 s | 8 | 5e-2 |
| mass above the LFL, near its peak | 7 | 3.5e-2 |

The point count matching exactly at every instant is the strongest single
check in the suite: it says the same observers were in range at the same
times, which requires the release window, both edge crossings, the whole
downwind integration and the arrival-time mapping to all be right
simultaneously.

**One mistake worth recording.** `TDIST0` is the observer's *first* recorded
distance -- `SORTS` sets `TABLE(21)` once before the loop and never updates it
-- so `x - x0` is how far the plume has travelled from the source. Reading it
as the spacing between consecutive output points makes the along-wind
correction never fire, because consecutive points are far closer together than
the 130 m cut-off. The symptom was a silent one: results that looked plausible
and were 3-10 % high on flammable mass. Fixing it halved the error.

### Receptor time histories (`DEG4`)

The same sorting machinery as `DEG3` with a different choice of times:
`GETTD` brackets the window during which the cloud passes one receptor rather
than spanning the whole cloud lifetime, and `DOSOUT` interpolates each
snapshot to the receptor's own distance. The two views are assembled from the
same observer profiles by different routes -- `DEG3` interpolates in time at
fixed instants, `DEG4` in distance at fixed receptors -- and agree on the peak
at 400 m to 2 %.

### Jet/plume (`JETPLU`)

Compared against the full-precision `.VEL` trajectory dump, not the listing.

| case | rows compared | max relative deviation | touchdown |
|---|---|---|---|
| EX1 (MIC, buoyant, never lands) | 44 | 5.4e-6 | correctly reports none |
| EX2 (ammonia pipeline) | 31 | 3e-2 near the orifice, 4e-4 beyond | 3.0e-5 |
| EX3 (as EX2, simplified density) | 30 | 2e-3 | 5.2e-6 |

EX2's first ten metres are extremely stiff -- the excess centreline velocity
falls from 430 to 15 m/s -- so the two step sequences separate there and
reconverge downstream. The landing point, which is the only thing that passes
to the ground-level model, agrees to 3e-5 in all three quantities.

**One mistake worth recording.** `MODOUT` stores the values it will
interpolate with *before* adding the ground image, then reports the imaged
concentration. Recording the bare plume instead left the centreline 35 % low
where the plume approached the ground -- and exactly right everywhere else,
since the image is negligible while the plume is high.

### Jet-to-ground bridge (`DEGBRIDG`)

The landed plume is replaced by a circular source of the plume's half-width,
releasing the same mass rate but *already diluted* to the concentration the
plume had when it landed. That dilution is the point of the bridge: restarting
from pure contaminant would over-predict everything downwind. The touchdown
distance becomes `OODIST`, an offset added to every reported distance so the
two stages share one coordinate.

Compared against the `.INP` deck the Fortran `DEGBRIDG` wrote for the same
touchdown:

| check | size | max relative deviation |
|---|---|---|
| 17 scalars (EX2 and EX3) | 34 | 1.2e-7 |
| source table, 4 rows x 8 columns | 64 | 3.5e-7 |
| density table | 22 x 5, 2 x 5 | 0 |

These are at the deck's print precision of seven significant figures.

### End to end, in Python alone

EX2 and EX3 run jet, bridge, source and downwind without touching the
Fortran. Concentration interpolated onto the reference profile's own
distances agrees to 1e-4 over most of the range and 3e-2 at worst.

## Three bugs found in DEGADIS

Distinct from the two above, which were mistakes in the port. These are
defects in the original, confirmed against the Fortran and reproduced.

**`ADIABAT(ifl=1)` reads its `wa` argument** rather than computing it. It
writes `wa` only in the out-of-range branches, so a normal call uses whatever
the caller passed — to form the molecular weight, hence the mole fraction,
hence which panel of the mixing table is interpolated on. `SZLOCAL` passes an
uninitialised local; under `/noauto` that is static, zero, and never written.
So `SZF` runs every lookup as though the mixture contained no dry air,
selecting panels one or two positions off and shifting the layer density by up
to 0.5 %. Reproducing it took `SZF` from 4.7e-3 to 4.4e-15.

**`SURFAC` gets a mix of heated and adiabatic layer properties.** `ADDHEAT`
writes `rholay`, `temlay` and `cp`; the `ADIABAT` call after it writes its
density into `rholam`, a different variable, so the heated density survives
while the temperature is replaced. `PSS` and `SSG` then differ by a single
letter — `PSS` writes into `temlay`, `SSG` into `temlam` — so the two stages
compute the ground heat flux, and apply the cut-off that zeroes it, from
different temperatures. Whether that is intentional is not recoverable from
the source.

**`GAMINC` returns the unregularised incomplete gamma.** It computes the
Numerical Recipes `gammp`, which is regularised, then undoes the
normalisation: `gaminc = exp(log(aa) + gln)`. In `PSS` the result is compared
against a cap built from Γ(1/(1+α)), so treating it as regularised leaves the
flammable-mass derivative low by exactly that factor — 6.5 % for this wind
profile — while every other quantity in the routine remains exact. The
asymmetry is what made it findable.

## The one residual, and how it was removed

`ADDHEAT` resolves the heated temperature with `ZBRENT` at `acrit = 1e-3` K,
then forms the heat capacity as `dh/(temp − amt)`. Early in the dense phase
`temp − amt` is only a few hundredths of a kelvin, so the quotient inherits a
percent-level uncertainty; it enters the forced-convection term of `SURFAC`,
which sets the ground heat flux, which integrates back into `dh`.

At the first output point of the dense phase the distance and the bisection
count agree to 1e-12 and the other three integrated states to 1e-6, while the
added-heat state differs by 1.6 %. That asymmetry is what localises it.

`legacy_numerics=False` removes it: the inversion is solved to 1e-12 K, and
the same secant is then well conditioned however small the interval, because
it is exactly the mean heat capacity over the heating the cloud received.

**A wrong first attempt, worth recording.** Replacing the quotient with an
analytic derivative — a ±0.01 K centred difference — made 26 of 30 observers
in the transient case terminate after a single step with `IHLF = 11`, the
integrator bisecting to a standstill. Mixture enthalpy has a slope
discontinuity where water begins to condense: the derivative jumps there, a
secant over a finite interval does not, and narrowing the interval towards a
derivative recovers the jump. The heat capacity is therefore a ±5 K centred
secant — the same width the original uses, centred to remove its bias.

## Deliberate fidelity choices

Three findings required reproducing an approximation rather than improving it.

**Single-precision dilution grid.** `SETDEN` writes
`zbda = (float(i)/float(iils))/(1+humid)`. The Fortran `FLOAT` intrinsic
returns `REAL*4`, so `199/200` is `0.99500000476837158` rather than `0.995`.
Every table node inherits about 1e-6 relative error. Using the intended
double-precision grid moves all 28 nodes by that amount, which then propagates
through every downstream interpolation. `exact_grid=True` selects the intended
behaviour; the default reproduces the original.

**`ALPH` quadrature error.** The residual whose root defines `alpha` is
integrated by `RKGST` at `ERBNDZ = 0.005`. Evaluating that residual accurately
at the Fortran's own answer gives 4.17e-4, not zero: the reported `alpha` is
the root of the *approximated* integral. The accurate root differs by 1.9e-5
relative. Since `alpha` appears as an exponent throughout, `RKGST` is ported
exactly and used by default.

**`GSERIES` clamp.** The hypergeometric series backing `PHIHAT` clamps its
argument at 0.999 and stops at 7e-5 relative. Against `scipy.special.hyp2f1`
the difference reaches 1.6e-4 in `PHIHAT`.

## Backend comparison, whole model

Burro 9, the steady LNG pool spill, and the same spill as a transient. The jet
cases are unaffected because a `.INO` deck supplies its own density table and
no mixing line is built.

| Burro 9, steady | legacy | coolprop | change |
|---|---|---|---|
| source enthalpy | -409 207 J/kg | -434 442 J/kg | **-6.2 %** |
| distance to the UFL | 245.9 m | 243.8 m | -0.9 % |
| distance to the LFL | 497.1 m | 494.5 m | -0.5 % |
| mass above the LFL | 6732 kg | 6640 kg | -1.4 % |
| mass between UFL and LFL | 3510 kg | 3568 kg | +1.7 % |

| Burro 9, transient | legacy | coolprop | change |
|---|---|---|---|
| snapshots produced | 19 | 19 | structure preserved |
| peak flammable mass | 6619 kg at 79 s | 6393 kg at 79 s | -3.4 % |
| peak mole fraction at 400 m | 0.0725 | 0.0743 | +2.5 % |

At the level of the mixing table alone, water properties are negligible for
LNG (0.001 % on density) while the contaminant heat capacity moves the
buoyancy excess by 7 %. The whole-model numbers confirm that: everything above
traces back to the contaminant.

But the buoyancy excess moves 7 % and the final distances only half a per
cent. The density change enters gravity slumping as a square root and
entrainment suppression through the Phi correction, and the two largely
cancel. Whether that cancellation survives comparison with measurements is
precisely what the two backends exist to test.

### Cost, and why the properties are tabulated

CoolProp costs about fifteen times the legacy runtime. Each property is a
function of temperature alone at fixed ambient pressure, and the root finders
evaluate them hundreds of thousands of times over a narrow range, so each is
tabulated once on a 4001-point grid spanning 80-400 K and interpolated after
that.

Two details were forced by the model rather than chosen:

* **Interpolation, not a rounded cache.** Rounding the temperature to a cache
  key is faster, but it puts small steps into functions the integrators
  effectively differentiate, and the adaptive step control chases them.
* **The grid is built once over the full span.** Growing it outward as queries
  arrive cost 361 rebuilds on a single steady run, because the model walks its
  temperature range gradually rather than jumping to the extremes.

A third was a real bug. Part of the span has no single-phase data -- a
cryogenic release sits below its own saturation line at ambient pressure --
and one failed point originally abandoned the whole grid, silently dropping
the run back to the 1989 correlation while still reporting success. The
symptom was that CoolProp results moved *towards* legacy instead of away. The
nearest valid value is now held across the gap.
