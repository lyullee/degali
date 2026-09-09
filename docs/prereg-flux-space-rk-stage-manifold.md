# Pre-registration: constrained RK stage reconstruction on the flux manifold

Date frozen: 2026-09-09, after retaining two failed trial 11/23 extension
runs and before adding a distinct intermediate-stage tolerance.

## Diagnosis retained

At the second 0.005 m step, trial 11 and trial 23 stop only at the third RK4
stage. Their phase-partitioned six-flux residuals are `7.55e-6` and `1.04e-5`.
A numerical-Jacobian retry gives the same minima, ruling out the initially
suspected optimizer-Jacobian stagnation. When these diagnostic intermediate
sections are allowed, the following RK stage and accepted endpoint again
reconstruct near `3e-11`.

The finite-step RK affine combination therefore leaves the nonlinear
four-parameter section-flux manifold at an internal stage. It is not an
accepted physical state and cannot always satisfy the endpoint inverse limit.

## Frozen constrained-stage rule

- Keep the six physical fluxes as the unmodified primary RK state.
- Keep the accepted endpoint inverse limit at `1e-8`.
- Keep the same parameter bounds, phase fluxes, source closure and analytic
  Jacobian; do not add a numerical-Jacobian retry.
- Permit a maximum scaled residual of `2e-5` only when reconstructing the four
  internal RK source-evaluation stages.
- Do not overwrite or project the primary six-flux vector to the approximate
  stage section. The approximation is used only to evaluate its local source.
- Record the largest internal-stage residual.
- Any internal residual above `2e-5`, any accepted endpoint residual above
  `1e-8`, or any existing physical gate failure stops the trial.

This rule may lower formal RK order on the constrained manifold. It is
accepted only if the already frozen 0.005/0.0025 m terminal state and flux
differences remain below `0.005`, and independent order-16 balances remain
below `1e-5`. Rerun trials 11 and 23 into a third immutable file; retain both
earlier failed files. No receptor score or model promotion is allowed here.
