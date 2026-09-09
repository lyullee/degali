# Pre-registration: stored wind-vector fields through the five-sensor operator

Date frozen: 2026-09-08, before projecting the stored yawed trajectories at
the paired PRESLHY vertical crosses.

## Question

Can the independently measured Trial 10 wind direction and/or magnitude
explain the excessive scalar amplitude at the fixed y=0 instrument line without
hiding the thermal/species variance error?  This is a stored-field observation
audit.  It does not fit wind direction, plume displacement, diffusivity or a
thermal source to temperature or hydrogen measurements.

## Frozen scope

Use Trial 10 only.  Trial 23's far-field wind direction was reported faulty and
its sparse local compass records are not simultaneous plume-direction inputs.
Do not select one of them from its temperature response.

Compare these already completed fields:

1. planar density-profile control;
2. direction-only yaw, step 0.02 m;
3. observed-magnitude-only, step 0.02 m;
4. observed magnitude plus direction, step 0.02 m;
5. observed magnitude plus direction, step 0.01 m.

The observed direction is the previously frozen 49-sample vector mean for the
same release window, 21.304311 degrees relative to the nozzle axis.  The
observed speed is the magnitude of that mean vector, 1.922264 m/s at 3 m.
Neither value may be changed in this audit.  No trajectory is re-integrated.

## Stored-field integrity

- Hash every field, completion record, input record, observation input and
  executable source used for reconstruction.
- Require each stored field's completed flag, accepted handoff, full coverage
  and global six-flux balance residual below 1e-5.
- At nine fixed state indices, recompute the six fluxes and sources with the
  current reconstructed model.  The maximum component-scaled difference from
  the stored arrays must be at most 1e-8.
- The two measured-vector steps must use identical physical wind inputs.

## Identical observation operator

At x=1.78 and 4.00 m, y=0, evaluate the yawed trajectory at the uniquely
bracketed horizontal-normal cut.  Sample the five absolute receptor heights
formed by adding `[-0.50, -0.25, 0, +0.25, +0.50] m` to the 0.50 m release
height.  Convert the yawed receptor's H2 mole fraction from 0--1 to vol-% before
passing it to the existing operator.

Use the same Trial 10 raw workbook, window, nonnegative observation convention,
signal gates and trapezoidal finite-line moments as
`prereg-model-thermal-profile-observation.md`.  Compare model values with the
observed median and IQR.  Do not fit a continuous Gaussian to either profile.

## Frozen diagnostics and gates

Keep centre amplitude, finite-line zeroth moment and truncated variance
separate for thermal deficit and H2.  For each group, report the median absolute
log ratio across the two stations and the station-resolved model/observed ratio.

A wind cell supports a lateral-displacement explanation only if, relative to
the planar control, all four amplitude/load summaries improve:

- thermal centre amplitude;
- thermal finite-line zeroth moment;
- H2 centre amplitude;
- H2 finite-line zeroth moment;

and neither the thermal nor H2 variance summary worsens.  At least one of the
six inequalities must be strict.  A lower combined score cannot override a
failed individual group.

For measured-vector step convergence, compare the exact ten projected values
(five thermal deficits plus five H2 values) at both stations.  The maximum
relative difference, scaled by `max(abs(step0.01), abs(step0.02), 1e-6)`, must
be at most 0.005.  The same limit applies separately to each reported moment.

Passing the displacement gate would identify a mechanism, not promote the
yawed model.  It must still retain the earlier concentration/geometry and
conservation limitations.  Failure directs the next implementation to a
separate thermal second-moment closure plus an independently constrained
species/trajectory correction.  No default changes in either case.
