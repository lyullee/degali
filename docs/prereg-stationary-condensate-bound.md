# Pre-registration: stationary condensed-air momentum bound

Date fixed: 2026-09-04, before any PRESLHY or Spadeadam result from this
candidate was calculated.

## Why a bound comes before another fitted particle model

The current transported source assumes condensed N2/O2 and gas have one axial
velocity.  Li et al. (2026), assumptions 6--7 and equations 13--15, permit
different phase velocities and set the condensed-phase velocity to negligible.
These are opposite kinematic limits.

Source-only response diagnostics were calculated without looking at field
residuals.  For a 1 um solid-N2 particle, the Stokes relaxation time at the
20.37 K endpoint is about 57 us.  Relative to the post-endpoint transport
residence time, `St = tau_p/t_source = 0.0003--0.065`, supporting the carried
limit there.  Relative to the previously omitted evaporation-zone residence,
however, `St = 0.006--1.13`; finite axial slip can matter for the fastest small
orifices.  No measured particle-size distribution or formation history exists
for choosing an intermediate relaxation exposure.

Therefore candidate S0 implements Li's zero-condensate-velocity limit with the
present stable solid-N2/O2 thermodynamics.  It is a model-form bound, not a
calibration.  If even the opposite bound moves errors in the wrong direction,
an intermediate slip parameter must not be tuned to field data.

## Equations fixed in advance

At the H2 evaporation endpoint, let `r_g` and `r_c` be gaseous and condensed
air mass per unit H2, `q = u_p/u_g`, and `j = J/m_H2`.  Momentum and kinetic
energy per unit H2 are

`u_g = j / (1 + r_g + q*r_c)`

`k = 0.5*u_g**2 * (1 + r_g + q**2*r_c)`.

The existing total-energy root becomes `h_phase + k = h_store`.  The carried
model is `q=1`; S0 is `q=0`.

In the subsequent plug source, every equilibrium condensed increment is
removed with zero axial momentum, as required by the S0 boundary.  It carries
its own condensed-phase enthalpy out of the atmospheric stream.  Remaining
gas carries all conserved axial momentum, so its kinetic term is
`J**2/(2*m_gas)`.  Component mass, axial momentum and total energy including
the deposited phase must close below `1e-6`.  No settling speed or particle
diameter enters S0.

## Comparator and predictions

S0 uses the energy-consistent no-pressure-thrust source.  The comparator is
the carried 1 um E source; the validated corrected default is reported
separately.

- S0 must deposit at least the initial condensed-air mass and hand off a
  colder, more H2-rich gas than E.
- Removing cold solid mass can make the remaining stream less negative in
  enthalpy, while retaining low H2 temperature can make it denser.  The net
  buoyancy and concentration directions are not pre-assigned.
- S0 should bound any later finite-relaxation result; an intermediate model
  outside the S0--E source-state range is an implementation error.

## Fixed validation and adoption rule

- PRESLHY horizontal momentum-driven window-mean sensor-arc maxima on common
  arcs.
- PRESLHY common well-constrained vertical fits with the momentum filter.
- Spadeadam horizontal tests 4 and 6 at low-mast wind.

Adoption requires PRESLHY MG 0.7--1.3 and no farther from one than the
corrected default, VG no greater than 1.455 and FAC2 at least 0.82; centre MAE
no greater than 0.160 m and absolute mean width-ratio error no greater than
0.033; Spadeadam MG closer to one, VG no greater than 1.506 and FAC2 at least
0.65.  Both campaigns and all conservation checks must pass.

If rejected, report handoff distance/temperature/H2 fraction, deposited mass,
all invariant residuals, PRESLHY MG/VG/FAC2 and geometry, and Spadeadam
MG/VG/FAC2 and centre heights.  S0 remains a diagnostic bound only; it is not
made selectable in the validated preset.

## Results after pre-registration

All eleven sources close and reach the S0 handoff in `0.0008--0.0045 m` at
`21.171--21.200 K`.  The retained gas is H2 to more than eight nines by mass;
the deposited N2/O2 is `0.680--0.894 kg/kg-H2`.  Maximum relative mass,
momentum and total-energy residuals are `2.1e-16`, zero and `6.2e-12`.
Particle diameter has no numerical effect, as required for the declared
zero-velocity bound.

S0 reaches 66 PRESLHY arcs and gives full-sample MG/VG/FAC2
`1.259/2.196/0.879`.  On the 62 arcs shared with the corrected default:

| model | MG | VG | FAC2 |
|---|---:|---:|---:|
| corrected default | 1.047 | 1.425 | 0.839 |
| S0 | **1.285** | **2.308** | 0.871 |

On the narrower 38-arc sample shared with E, E gives
`0.965/2.094/0.842` and S0 gives `1.488/3.753/0.816`.  Thus the opposite
kinematic limit crosses the bias but increases rather than bounds down the
residual variance.

The 23 common vertical fits are also worse:

| model | centre-height MAE | mean width ratio |
|---|---:|---:|
| corrected default | 0.145 m | 1.033 |
| S0 | **0.321 m** | **1.269** |

On the 17 fits shared with E, E is `0.284 m / 1.014` and S0 is
`0.422 m / 1.287`.  Removing stationary cold solids leaves a nearly pure H2
stream that becomes more buoyant downstream; its low initial temperature does
not hold the trajectory down.

The independent six-arc Spadeadam comparison is:

| model | MG | VG | FAC2 |
|---|---:|---:|---:|
| corrected default | 0.869 | 1.476 | 0.667 |
| E, carried 1 um | 1.336 | 1.583 | 0.833 |
| S0 | **1.563** | **1.827** | 0.833 |

Test-6 centre height at 30/50/100 m rises from the default
`3.35/6.78/13.65 m` to `5.04/9.03/16.99 m` under S0.

## Decision

S0 is rejected.  It remains an explicit diagnostic option, not a preset.
Both source-kinematic limits now fail the same independent direction test:
the carried E source and stationary S0 source over-rise relative to the
corrected default.  Without a measured particle formation history, selecting
an intermediate relaxation exposure from field residuals would be tuning.
The next defensible source question is finite-rate heat/mass transfer, which
can delay condensed-air sublimation, rather than another axial-slip scalar.
