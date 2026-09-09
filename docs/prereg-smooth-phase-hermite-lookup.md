# Pre-registration: C1 phase-property lookup confirmation

Date frozen: 2026-09-10, after an exploratory 4,000-point screen and before
the independent confirmation run or production implementation.

## Retained exploratory result

The existing 201 by 161 bilinear lookup is continuous but has derivative
jumps at arbitrary grid lines. Trial 24 approached the nearly coincident
57.325 K and H2 mass-fraction 0.30625 cell edges and could not continue in
flux space. A conventional bicubic spline reduced lookup error but produced
21 nonnegative `dH/drho` points in a 513,441-point screen and is rejected.

A tensor PCHIP retained negative `dH/drho` in 30,000 random points and reduced
error against the existing exact equilibrium oracle. A faster exploratory C1
bicubic Hermite surface, using PCHIP nodal slopes and a symmetric mixed nodal
derivative, matched the tensor-PCHIP errors and retained negative `dH/drho`
at all 513,441 dense points. These exploratory values are not confirmation.

## Frozen confirmation

- Rebuild the same trial 24 thermodynamics without fitting observations.
- Use random seed 240911 and 8,000 uniformly sampled interior lookup points,
  distinct from the exploratory seed and sample.
- Compare bilinear, tensor PCHIP and C1 Hermite temperature and volumetric
  enthalpy directly with `_condensed_air_state_exact` at identical states.
- Require Hermite RMS and p99 errors to be lower than bilinear for both
  temperature and enthalpy. Maximum relative enthalpy error is diagnostic
  because exact enthalpy crosses zero.
- Require Hermite reproduction at every source lookup node within `1e-10 K`
  and `1e-6` enthalpy-density units.
- Require finite and strictly negative fixed-H2-density `dH/drho` at the same
  801 by 641 grid used exploratorily: 513,441 points total.
- Retain the conventional bicubic negative control and the tensor-PCHIP
  reference in the output. No observational scores are calculated.
- Hash the screening code, this registration, and the existing exact phase
  oracle before and after the calculation; write a new immutable JSON file.

Passing this screen permits an opt-in implementation and unit tests only. It
does not permit replacing the package default, claiming field validation, or
discarding the retained trial 24 integration failures. Trial 24 continuation,
full phase-domain comparison and existing regression tests remain separate
gates.
