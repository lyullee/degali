# Pre-registration: thermal-moment field observation comparison

Date frozen: 2026-09-10, before extending the direct-flux candidate to 4 m or
calculating its sensor predictions.

## Question

Does the coefficient-free, independently transported enthalpy second moment
improve both the temperature-deficit amplitude and its truncated vertical
variance when the candidate and control are passed through the same finite
PRESLHY sensor operator?

This is an observation comparison, not a parameter fit.  The numerical
six-trial extension is already complete and is not rescored here.

## Frozen population and observations

- PRESLHY E3.5 Trials 10 and 23 only.  These are the two trials for which the
  paired temperature and H2 profiles have already passed the frozen timing and
  signal gates.
- Stations x = 1.78 and 4.00 m.
- Five vertical offsets from the release axis: -0.50, -0.25, 0, +0.25 and
  +0.50 m.  Add the trial release height before querying the model.
- Reuse the observed medians and interquartile ranges in the immutable
  `model_thermal_profile_observation_2026-09-08.json` evidence.  Do not select
  a time, sensor or station after seeing the new prediction.
- Temperature deficit, H2 concentration, finite-line zeroth moment, centroid,
  variance and width use the exact operator frozen in
  `prereg-model-thermal-profile-observation.md`.

## Frozen model paths

The control is the stored density-profile trajectory
`phase_ambient_consistency_field_complete_2026-09-05.json`, replayed without
modification.

The candidate starts from the accepted six-moment boundary for each trial and
marches mass, H2 mass, horizontal and vertical momentum, total energy and the
advective-enthalpy second moment directly.  Retain:

- `thermal_species_ratio = 1`;
- `mechanical_work = reduced_buoyancy_work`;
- the default linear phase lookup;
- free planar transport, with zero yaw and no ground-contact dynamics;
- the exact positive analytic width split;
- the established inverse-only step rejection and eight-success growth rule.

The receptor calculation retains the ground image as an observation operator.
Species and enthalpy images are superposed separately: for exponent `q` and
`p = 1 / beta_H^2`, `C` is proportional to the sum of `exp(-q)` images while
`H` is proportional to the sum of `exp(-p q)` images.  Raising the already
summed species image to `p` is forbidden because it creates a nonphysical
cross term.

## Frozen numerics and failure rules

- Coarse nominal arc step: 0.005 m.
- Refined nominal arc step: 0.0025 m.
- Target x: 4.00 m; maximum arc length: 5.0 m; accepted-step ceiling: 12000.
- Stage inverse tolerance: 2e-5; accepted endpoint tolerance: 1e-8.
- Minimum adaptive step: 7.8125e-5 m.
- Stop on the first unhandled inverse, phase, weak-budget, inward-flow,
  diffusivity, curvature or receptor-domain failure.  Do not clip or
  extrapolate.
- Independently reconstructed terminal balance must be at most 1e-5.
- Coarse/refined terminal physical-parameter and six-flux differences must
  each be at most 0.005.
- Across all 40 paired scalar predictions (2 trials x 2 stations x 5 sensors x
  temperature/H2), the maximum coarse/refined scaled difference must be at
  most 0.005.  Temperature is scaled by the larger ambient-relative deficit
  with a 1 K floor; H2 vol-% by the larger value with a 0.1 vol-% floor.
- Every executable input is hashed before and after the run.  Any change
  invalidates the result.

## Frozen decision

As before, calculate the median absolute log ratio separately for the four
thermal centre amplitudes and the four thermal variances.  The candidate is
directionally better only if both medians are no larger than the control and
at least one is strictly smaller, in both coarse and refined runs.  It may
advance as a research candidate only if that observation decision and every
numerical gate above pass.

H2 amplitude and variance, thermal/H2 zeroth moments, centroids, individual
station ratios and IQR membership remain visible diagnostics; they may not be
merged into the primary decision or used to tune a coefficient.  Passing does
not validate safety distances or promote a package default.  Failure is
retained and reported without retuning.
