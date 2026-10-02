# Sandia steady-crosswind LH2 pool contour screen

Date: 2026-09-17

## Scope

This protocol reduces selected contours from Hecht, *Liquid Hydrogen Pooling
and Vaporization Experiments*, SAND2025-08847, DOI
10.2172/2585576. The source PDF and any digitised coordinates remain under
the ignored `reference/` tree. They are not package data.

The report contains 16 continuous LH2 spills in a large tunnel with controlled
cross-wind. For low-flow, contained, instrumented-concrete tests, its measured
substrate heat transfer is reported to agree with the latent power required to
vaporize the inflowing LH2. Only those cases may test the current pure-vapour
pool source using the measured liquid inflow as its vapour mass rate. A spill
that runs off the instrumented substrate is excluded: its pool-vapour source
cannot be inferred from the nominal liquid feed alone.

This is a deliberately new regime screen, not a replication of the historical
NASA pool comparison: the latter used one 9.1 m pond at roughly 9--10 kg/s.
It provides a lift-off-height check, but not a small-pool concentration claim.

## Observation reduction

For each selected figure panel, record locally:

1. report figure and panel;
2. contour mole fraction;
3. the furthest down-wind plane at which that contour is visibly resolved;
4. one grid/line-width spatial resolution; and
5. whether the heat-transfer evidence supports complete inflow vaporization.

The resulting value is a **lower bound** on the real centreline reach. An
off-centre contour at a given plane proves that the true centreline reaches at
least that mole fraction there. It does not identify the centreline, plume
width, instantaneous maximum, source heat flux, or an exact contour endpoint.

## Decision rule

Compare the model centreline reach to the same mole fraction. If it terminates
before the resolved plane, allowing only the stated plotting resolution, the
model is falsified for that case. Reaching farther is merely `not_falsified`;
it is not an accuracy score, a fitted tolerance, or a default-promotion gate.

`tools/audit_sandia_pool_contours.py` reads two *local* JSON files: the
reduction and declared model reaches. The library deliberately provides no
embedded observation table or default source rate.

## Exclusions

- No high-flow spill-off case is used to calibrate a pool-vapour fraction.
- No concentration contour is compared as though it were a centreline sensor.
- No transient video plume angle is used to set a steady trajectory; the
  report itself identifies projective-transform disagreement.
- No temperature-to-H2 conversion is silently rescaled to match a model.

The next stage, if the pure-vapour pool path is falsified by contained cases,
is an explicit pool energy/evaporation boundary using independently declared
substrate, gas-convection and radiation terms. It is not an entrainment
coefficient adjustment.
