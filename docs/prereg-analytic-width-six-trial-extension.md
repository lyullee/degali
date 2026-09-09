# Supplemental preregistration: analytic-width six-trial extension

Date frozen: 2026-09-10 (Asia/Seoul), after the Trial 24 x=0.43 m pilot and
before any analytic-width march to the six frozen observation targets.

## Reason for the supplement

The original six-trial extension preregistration fixes the six stored boundary
moments and requires their prior numerical state to replay to `1e-9`.  Replacing
only the independent-enthalpy add-on's quantized width root with the exact
positive algebraic root changes the coordinate representation of those same
moments.  Direct replay errors of the old states are between `3.81e-10` and
`5.78e-8` across the six trials.  Reusing those states would therefore mix old
geometry coordinates with the corrected geometry.

## Frozen boundary transfer

- Keep every stored boundary moment, source, ambient state, target station,
  transport choice, coefficient, and acceptance tolerance unchanged.
- For each selected trial, project its stored six moments onto the linear
  phase-interpolation section using the stored state and thermal-width ratio as
  the initial guess.
- Require projection success, maximum scaled moment residual `1e-8`, maximum
  coarse/fine quadrature residual `1e-5`, and thermal width strictly inside
  0.5--2.0.
- Record both the old-state replay error and the complete corrected projection.
- The package default phase interpolation remains `linear`; the C1 result is a
  separate pilot and is not substituted into this six-trial extension.

## Frozen march and decisions

All targets, coarse/refined steps, adaptive minimum, inverse tolerances, local
physical gates, independent order-16 balances, convergence limit, hashing, and
no-observation-scoring decision remain exactly those in
`prereg-thermal-moment-flux-space-extension.md` and its registered adaptive
supplements.

Trials 12 and 24 shall be run first because they were the retained failures.
The other four may be run only after that pair is recorded.  No full six-trial
claim is permitted until all six rows pass in one or combined immutable files.

This projection is a change of numerical coordinates for fixed boundary
moments.  It is not a fit to downstream observations.
