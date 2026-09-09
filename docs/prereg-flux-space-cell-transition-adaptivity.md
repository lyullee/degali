# Pre-registration: phase-cell transition continuation

Date frozen: 2026-09-10, after retaining the trial 12/24 extension failures
and transition diagnostics, and before changing the marcher.

This registration supersedes only the earlier 0.0003125 m minimum and
no-regrowth clauses. All earlier physical gates and tolerances remain active.

## Retained diagnosis

The first trial 12/24 extension file is immutable and failed its gates.

- Trial 12 coarse reached x=1.78 m without rejection. Its refined path stopped
  after 23 accepted steps and three halvings. At the failed endpoint, the RK4
  last-stage warm start converged to scaled residual `4.364e-3`. The same
  physical target reconstructed to approximately `2.53e-11` from the previous
  accepted section with the unchanged analytic Jacobian, and also reconstructed
  with a numerical Jacobian. This is a warm-start branch-selection failure,
  not loss of a physical endpoint.
- Trial 24 stopped after 0.0071875 m. At the minimum 0.0003125 m step its third
  RK stage had residual `1.227e-3`. Analytic and numerical Jacobians, finite
  difference steps from `1e-4` through `1e-2`, and thermal-width seeds from
  0.8 through 1.2 did not recover a root. Euler, midpoint and Heun trial steps
  also did not meet the endpoint tolerance.
- From the same last accepted trial 24 section, RK4 at 0.00015625 m passed all
  internal stages but not the endpoint. RK4 at 0.000078125 m passed all stages
  and the endpoint at approximately `2.29e-11`. Thus the physical section has
  not disappeared; the local flux-to-section map has a narrow high-curvature
  transition associated with the piecewise-bilinear phase table.

## Frozen continuation rule

- Keep the six transported physical fluxes unmodified.
- Keep the internal-stage inverse tolerance at `2e-5` and accepted endpoint
  tolerance at `1e-8`.
- At an accepted endpoint only, first use the last RK-stage section as before.
  If that inverse fails, retry the identical target once from the previous
  accepted section with the same analytic Jacobian and bounds. Do not use a
  numerical-Jacobian fallback or alter the target. Record attempts and
  successful recoveries.
- If an internal inverse fails, or both endpoint starts fail, reject and halve
  the whole step. Physical closure, phase-domain, weak-budget, inward-flow,
  diffusion and curvature failures still stop immediately.
- Lower the ordinary minimum step to 0.000078125 m. Do not go below it except
  for a successful target-shortened final step.
- After eight consecutive accepted steps at a reduced step, double the next
  ordinary step, capped by that march's original 0.005 or 0.0025 m maximum.
  A failed growth attempt is an ordinary rejected step and halves again. Reset
  the success counter on every rejection and growth.
- Keep the accepted-step ceiling at 4000 and record rejected steps, minimum
  accepted step, step growths, endpoint fallback attempts and endpoint fallback
  successes.

## Frozen decision

Rerun trials 12 and 24 into a new immutable evidence file. Both must reach the
frozen target, retain independent order-16 balance below `1e-5`, and retain
coarse/refined terminal parameter and flux differences below `0.005`. The
earlier failed evidence remains retained. No observation score, coefficient
fit, default change or physical-model promotion follows from this numerical
continuation repair.
