# Pre-registration: endpoint-manifold step rejection

Date frozen: 2026-09-09, after retaining the trial 23 refined endpoint failure
and before adding inverse-specific step rejection.

## Retained failure

With the separately registered internal-stage rule, trial 11 passes. Trial 23
reaches x=0.79 m at the 0.005 m step, but its 0.0025 m march stops at the fourth
endpoint: maximum exact six-flux inverse residual `2.13e-8` versus the unchanged
`1e-8` limit. The local weak residual, inward flow and positive diffusivity
gates pass. The residual is distributed across mass, momentum, energy and the
thermal second moment, consistent with finite-step departure from the
nonlinear section manifold.

## Frozen rejection rule

- Introduce distinct inverse and local-closure exceptions.
- Only a six-flux inverse failure may reject and halve the entire RK step.
- A phase, weak-budget, inward-flow, diffusivity, curvature or other local
  physical failure still stops immediately.
- Keep the accepted endpoint inverse limit at `1e-8` and the internal-stage
  limit at `2e-5`; do not project or alter the primary flux vector.
- Start from the already frozen coarse/refined maximum steps of 0.005 and
  0.0025 m. After a rejection, retain the halved step for the rest of that
  march; do not regrow it.
- Minimum ordinary step: 0.0003125 m. A target-shortened final step may be
  smaller only when it succeeds; a failed step at or below the minimum stops.
- Record rejected-step count and minimum accepted step.

Rerun trials 11 and 23 in a new immutable file. Both independent order-16
balance and the existing coarse/refined `0.005` convergence gate remain
mandatory. No observational scoring or physical-model promotion follows from
this numerical repair alone.
