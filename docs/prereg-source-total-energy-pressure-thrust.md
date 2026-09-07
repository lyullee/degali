# Pre-registration: storage-to-source total energy and pressure thrust

Date fixed: 2026-09-04, before running PRESLHY or Spadeadam with either
candidate below.

## Fault found without looking at a new residual

The transported condensed-air source closes mass, momentum and energy from
its first multiphase station onward.  However, that first station was built in
two independent steps:

1. the storage-liquid enthalpy was balanced against H2 plus condensed-air
   phase enthalpy; and
2. a velocity was then assigned from source momentum and its kinetic energy
   was added to the transport ledger.

Consequently the storage-to-source control volume creates
`0.5 * m_total * u_source**2`.  This is small for a slow release but material
for the large Spadeadam jets.  The defect is a conservation error, not a
coefficient choice.

Separately, the adopted source momentum contains only `m_dot * u_orifice`.
For a choked release the throat pressure is above ambient and the physical
momentum flux also contains `(P_throat - P_ambient) * A_orifice`.  Existing
comments correctly reserve that term for the preceding notional-nozzle zone,
but no production path calculates it.

## Primary equations fixed in advance

HyRAM+ 6.0 section 3.2 supplies the complete sequence.

1. At each candidate throat pressure, keep storage entropy fixed and solve
   `u = sqrt(2 * (h0 - h(P,s0)))` (equation 44).
2. Select the pressure that maximises `rho * u` between ambient and storage
   pressure (equation 45).  An interior maximum is choked.
3. Use the measured mass flow to infer the discharge coefficient
   `Cd = m_dot / (rho_t * u_t * A)` rather than replacing the measured flow
   (equation 46).  Report every inferred value; `Cd > 1` is physically
   inconsistent and bars adoption for that trial/population.
4. Conserve notional-nozzle mass and momentum (equations 47--50):

   `J = m_dot * u_t * Cd + (P_t - P_a) * A`

   `u_eff = J / m_dot` and `A_eff = m_dot / (rho_eff * u_eff)`.

5. Close notional-nozzle energy with
   `h_eff + 0.5*u_eff**2 = h_t + 0.5*u_t**2 = h0`
   (equation 51).  The ambient-pressure density is solved from that enthalpy,
   not imposed from storage temperature.
6. At the first fully evaporated-H2/solid-air station, solve phase enthalpy
   and kinetic energy together:

   `h_store = h_phase(T,r) + 0.5*J_specific**2/(1+r)`.

   The H2/N2/O2 partial-pressure closure continues to determine `T`; the air
   ratio `r` is now the positive root of the total-energy equation.

No pressure, velocity, discharge coefficient, temperature, particle size or
entrainment coefficient will be fitted to concentration data.

## Two candidates, evaluated separately

### E: energy-only correction

Keep the currently adopted no-pressure-thrust specific momentum, but include
its kinetic energy in the evaporation endpoint.  This isolates the existing
energy defect without changing the momentum boundary.

### P: full isentropic-throat and pressure-thrust correction

Use the maximum-mass-flux throat and the energy-conserving Yuceil--Otugen
notional nozzle above.  Feed its conserved momentum and total energy into the
same condensed-air transport model.  Candidate P is evaluated only on trials
whose inferred `Cd` is physically admissible.

The fixed particle case for this source test is 1 um.  The preceding
pre-registration showed that 1 and 10 um are nearly carried and give the same
answer, while 100 um requires a separate axial-slip model and fails badly.
Selecting 1 um here avoids adding a second model-form variable; it is not
reselected after viewing these results.

## Predictions

- E removes positive artificial energy, requires at least as much condensed
  air at the H2 evaporation endpoint, and cannot increase the handoff H2 mass
  fraction.
- P has greater axial momentum than its own discharge-coefficient-weighted
  throat advection whenever the throat is choked.  Its momentum relative to E
  is not fixed because E obtains velocity from a different, post-flash
  homogeneous density.  If P is faster, it should reduce residence time for
  buoyant rise; if slower, the direction reverses.  The concentration-bias
  direction is not fixed.
- Both candidates must close storage-to-handoff component mass, axial
  momentum and total energy below `1e-6` relative at every source.
- If P infers `Cd > 1`, produces no ambient-pressure enthalpy state, or places
  the source downstream of the comparison range, that trial is reported as a
  source-model incompatibility rather than silently discarded.

Protocol correction recorded before any concentration or geometry metric was
run: the first version of the preceding momentum prediction incorrectly said
P must exceed E for every choked source.  The statement is true only relative
to P's own throat-advection term.  Source-only diagnostics exposed the unequal
density definitions; no acceptance threshold or candidate definition changed.

## Fixed comparisons and adoption rule

- PRESLHY E3.5: the same horizontal momentum-driven, window-mean,
  sensor-position arc maxima, on common arcs.
- PRESLHY vertical geometry: common well-constrained fits with the same
  momentum filter.
- Spadeadam horizontal tests 4 and 6 at low-mast wind.

The corrected default is the comparator.  A candidate is adopted only if:

- PRESLHY MG remains 0.7--1.3 and moves no farther from one, VG does not
  exceed 1.455 and FAC2 does not fall below 0.82;
- common-fit centre-height MAE does not exceed 0.160 m and absolute mean
  width-ratio error does not exceed 0.033;
- Spadeadam MG moves closer to one, VG does not exceed 1.506 and FAC2 does not
  fall below 0.65;
- the source invariants and admissibility checks above pass; and
- the same candidate passes both campaigns.  A pressure-thrust option that
  helps only one pressure range remains a sensitivity, not a correction.

## Required reporting if rejected

Report the throat-pressure ratio, inferred discharge-coefficient range,
pressure-thrust / advective-momentum ratio, endpoint air-ratio change,
handoff range, all three concentration metrics, centre MAE, width ratio and
the source invariant residuals.  A conservation correction remains in the
research source even if its dispersion score worsens, but it is not promoted
to the validated default without the campaign criteria.

## Results after pre-registration

The implementation follows the current HyRAM+ source code's discharge-
coefficient convention.  The HEM throat lies at `P_t/P_0 = 0.6675--0.7316`.
Across the nine eligible PRESLHY sources and Spadeadam tests 4 and 6, inferred
`Cd = 0.101--0.819` and pressure thrust is `1.50--93.95` times the
coefficient-weighted advective momentum.  No `Cd` exceeds one.  Nevertheless,
PRESLHY trials 10, 22 and 25 have no atmospheric-pressure thermodynamic state:
with their low measured flow and large geometric area, the pressure-work term
would require more kinetic energy than the stored stream carries.  They are
reported as incompatible and are not clipped or silently assigned a state.

E increases the endpoint air/H2 mass ratio in every source.  The static range
`0.656699--0.807192` becomes `0.659088--0.807243`; the smallest change occurs
in the low-momentum 1-barg trials and the largest in the fast small-orifice
trials.  P's eight admissible endpoints span `0.678683--0.849530`.  The
notional-nozzle mass, momentum and energy residuals are at or below
`1.4e-16` in the admissible cases.

All eleven E source marches reach the single-phase handoff at
`0.125--0.797 m`; the eight admissible P marches reach it at
`0.073--0.925 m`.  Their maximum relative mass / momentum / energy residuals
are respectively `7.4e-15 / 1.4e-15 / 8.4e-11` for E and
`7.6e-15 / 1.4e-15 / 9.9e-11` for P.

The PRESLHY concentration comparison, restricted separately to arcs reached
by the corrected default and each candidate, is:

| candidate | n | default MG / VG / FAC2 | candidate MG / VG / FAC2 |
|---|---:|---:|---:|
| E | 38 | 1.056 / 1.630 / 0.789 | **0.965 / 2.094 / 0.842** |
| P | 33 | 1.182 / 1.636 / 0.758 | **1.201 / 1.505 / 0.758** |

E moves bias toward one but exceeds the fixed VG limit.  P reduces variance
on its restricted, admissible sample but moves MG away from one and cannot
represent three of the nine eligible sources.

The matching vertical-geometry results are:

| candidate | n | default centre MAE / width ratio | candidate centre MAE / width ratio |
|---|---:|---:|---:|
| E | 17 | 0.187 m / 1.113 | **0.284 m / 1.014** |
| P | 15 | 0.186 m / 1.083 | **0.145 m / 1.083** |

E corrects mean width but moves the plume centre too high.  P improves centre
MAE, but leaves the width error outside the fixed `0.033` tolerance.

The independent six-arc Spadeadam result is decisive:

| model | MG | VG | FAC2 |
|---|---:|---:|---:|
| corrected default | 0.869 | 1.476 | 0.667 |
| E | **1.336** | **1.583** | 0.833 |
| P | **1.466** | **1.688** | 0.833 |

Both candidates move MG farther from one and exceed their pre-registered VG
ceiling.  Test-6 centre height at 30/50/100 m changes from
`3.35/6.78/13.65 m` to `4.60/8.60/16.11 m` for E and
`4.97/9.13/16.81 m` for P.

## Decision

Neither E nor P is promoted to the corrected default.  E remains available as
the explicitly selected, conservation-correct research source because it
removes a real storage-to-source energy defect; its worse field score shows
that the current common-velocity equilibrium handoff is not a better
dispersion source.  P remains a further research sensitivity only.  Its one
local improvement in centre height does not survive concentration, width,
independent-campaign or source-applicability checks.
