# Pre-registration: model projection through the thermal-profile observation operator

Date frozen: 2026-09-08, before evaluating either stored trajectory at the
paired five-sensor crosses.

## Question

When the current stored model trajectories are sampled at exactly the same
finite vertical coordinates as the PRESLHY Trial 10 and 23 instruments, do
they reproduce both the centre thermal-deficit amplitude and the thermal
second moment?  This audit distinguishes a profile-width error from an energy
amplitude error.  It does not fit a diffusivity, width, heat source or lag.

## Frozen model paths

Compare two already stored, non-reintegrated trajectories:

1. `phase_ambient_consistency_field_complete_2026-09-05.json`: the
   density-profile control with the consistent ambient phase closure.
2. `gaussian_enthalpy_profile_allflux_field_complete_2026-09-05.json`: the
   Gaussian enthalpy-profile path with the same thermal/species Gaussian width.

The second path tests a different thermodynamic profile representation, not an
independently transported thermal width.  Neither path may be modified or
re-integrated in this audit.  Record the SHA-256 hash of both files and of every
executable input.

## Frozen observations and coordinates

- Use the same official PRESLHY E3.5 Trial 10 and 23 workbooks and release
  windows as `prereg-observed-thermal-profile-moments.md`.
- Use x=1.78 and 4.00 m, y=0, with offsets
  `[-0.50, -0.25, 0, +0.25, +0.50] m` from the release axis.
- Convert each offset to an absolute model receptor height by adding the
  trial-specific release height.  Do not compare the relative offset directly
  with the model's ground-referenced height.
- Evaluate the model's ground-image receptor operator at those exact points.
  Do not evaluate a continuous-model Gaussian width and compare it directly
  with a five-point observed width.

## Frozen reductions

For each observed time and model profile, retain the definitions

    theta_i = max(T_ambient - T_i, 0)
    c_i     = clip(X_H2,i, 0, 100)

and calculate `M0`, centroid and truncated variance by the same trapezoidal
five-point operator.  The model profile must pass the same centre signal gates:
5 K for temperature deficit and 1 vol-% for hydrogen.  Failure is reported,
not replaced by a small positive value.

For each observed station/scalar report the median, 25th and 75th percentiles
of:

- centre amplitude;
- zeroth moment;
- centroid;
- truncated variance and its square-root width.

For each model report its single steady value and its ratio to the observed
median.  Also report whether that value lies inside the observed interquartile
range.  A zero or sign-changing denominator is not converted to a ratio.

## Frozen comparison and decision

The primary thermal comparison contains eight equally visible quantities:
centre temperature-deficit amplitude and truncated thermal variance at two
stations in each of two trials.  Summarise the absolute log ratio separately
for the four amplitudes and four variances; do not merge them into one score
that can hide a width/amplitude trade-off.  Hydrogen centre amplitude and
variance are a separate geometry/species diagnostic.

A stored path is only directionally better if both its median thermal-amplitude
log error and median thermal-variance log error are no larger than the control,
with at least one strictly smaller.  Regardless of that result, neither stored
path is a promotable independent thermal closure because neither evolves a
separate thermal second-moment state.

The next transport candidate may advance only if it:

1. improves or preserves both thermal error groups under this identical
   observation operator;
2. preserves mass, species, momentum and total-energy balances;
3. passes step refinement without changing the observational conclusion; and
4. does not select its coefficient from these same four profiles.

Observed maxima and minima may be retained as diagnostics, but the primary
steady comparison remains the pre-specified median and IQR.  No default model
setting changes in this audit.
