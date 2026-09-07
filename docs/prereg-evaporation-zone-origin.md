# Pre-registration: evaporation-zone axial origin

Date fixed: 2026-09-04, before running any concentration or geometry
comparison with this candidate.

## Fault found from the source equations

The multiphase source plane already contains the ambient air required to
finish evaporating H2, but it is currently assigned zero axial distance.  The
subsequent condensed-air transport distance is preserved when JetPlume starts;
the preceding evaporation/entrainment distance is not.  Thus the source state
and the source coordinate describe different control volumes.

Li et al. (2026) equation 22 gives the Zone-III length

`S_evap = (1 - Y_H2) * m_air / (E * rho_ambient)`

and equation 24 gives the momentum-driven local entrainment scale

`E = beta_A * sqrt(J / rho_ambient)`, with `beta_A = 0.281`.

The present source has stable solid N2/O2 rather than the paper's metastable
liquid-N2-only state, so its own solved `Y_H2`, total entrained-air flow and
conserved momentum `J` enter the same geometric relation.  No temperature,
distance, coefficient or particle diameter is fitted.

Source-only diagnostics, inspected before outcome metrics, give
`S_evap = 0.012--0.050 m` over the nine eligible PRESLHY sources and
`0.072--0.078 m` for Spadeadam tests 4 and 6.

## Candidate and comparator

Candidate Z adds `S_evap` to the explicit 1 um transported-source handoff
coordinate.  It uses the energy-consistent E source from
`prereg-source-total-energy-pressure-thrust.md`; pressure thrust stays off
because it failed campaign and source-applicability checks.  Thermodynamic
state, mass, momentum, energy, entrainment and particle dropout are unchanged.

The comparator is the identical E source without `S_evap`.  The validated
corrected default is also reported so that a local improvement over an already
rejected source cannot be mistaken for adoption.

## Predictions fixed in advance

- Candidate Z starts JetPlume farther downstream by exactly `S_evap` and
  changes no handoff state invariant.
- At a fixed sensor x, Z gives JetPlume less development distance.  It should
  not increase plume-centre height and should generally increase predicted
  concentration.
- Effects should be largest on the closest PRESLHY arcs and smaller at the
  30--100 m Spadeadam arcs.

## Fixed comparisons and adoption rule

- PRESLHY E3.5 horizontal momentum-driven, window-mean, sensor-position arc
  maxima on common arcs.
- PRESLHY well-constrained vertical fits with the same momentum filter.
- Spadeadam horizontal tests 4 and 6 at low-mast wind.

Z can replace the research E source only if every source coordinate advances
by its equation-22 value, all state variables and conservation residuals are
unchanged to numerical tolerance, and its common-data centre MAE and MG move
toward the corrected default without worse VG or FAC2.

It can enter the validated default only if the stronger campaign criteria
remain satisfied: PRESLHY MG 0.7--1.3 and no farther from one, VG no greater
than 1.455, FAC2 at least 0.82; centre MAE no greater than 0.160 m and absolute
width-ratio error no greater than 0.033; Spadeadam MG closer to one, VG no
greater than 1.506 and FAC2 at least 0.65.  The same candidate must pass both
campaigns.

## Required rejection report

Report the Zone-III distance range, common-arc MG/VG/FAC2, centre MAE and
width ratio, Spadeadam MG/VG/FAC2 and test-6 centre heights.  A coordinate
bookkeeping correction may remain available as a research option if it is
conserved but fails the field criteria; it does not become default by being
dimensionally correct.

## Results after pre-registration

The implementation adds `0.012--0.050 m` to the nine PRESLHY E-source
coordinates and `0.072--0.078 m` to the two Spadeadam coordinates.  A direct
integration test requires the complete transported source object to be
identical with the option on and off; only the initial axial coordinate may
change.

On all 38 common PRESLHY arcs, Z changes E's MG/VG/FAC2 from
`0.965/2.094/0.842` to `0.954/2.076/0.842`.  Scatter improves slightly, but
MG moves farther from one and remains much worse in variance than the
corrected default's `1.056/1.630/0.789` on the same subset.

On the 17 common vertical fits, E centre-height MAE / mean width ratio changes
from `0.284 m / 1.014` to `0.280 m / 1.005`.  The coordinate correction moves
both quantities in the predicted direction, but the centre improvement is
only 4 mm and remains far outside the `0.160 m` adoption limit.

The independent Spadeadam MG/VG/FAC2 changes from E's
`1.336/1.583/0.833` to `1.331/1.577/0.833`, still worse than the corrected
default's `0.869/1.476/0.667`.  Test-6 centre heights at 30/50/100 m change
only from `4.60/8.60/16.11 m` to `4.58/8.58/16.10 m`.

## Decision

Z is not adopted and does not replace E as the preferred research comparator,
because its PRESLHY MG moves in the wrong direction and neither campaign
approaches the default.  The explicit option remains useful for coordinate
audits: it fixes a real source-origin omission without pretending that a
12--78 mm shift repairs the much larger equilibrium-source buoyancy error.
