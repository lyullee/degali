# Pre-registration: flux-space march of the conservative thermal-moment state

Date frozen: 2026-09-09, before any Trial 10 flux-space reach calculation.

## Question

Does direct numerical integration of the six transported section quantities
remove the accumulated conservation drift seen when the eight profile
parameters are integrated and the fluxes are only checked afterwards?

This is a numerical reformulation of the already selected local closure.  It
does not add or fit a diffusivity, heat-transfer coefficient, source term, or
profile parameter.

## Frozen model and state

- Trial 10 only for the first reach test.
- Start from the accepted six-moment boundary in
  `buoyancy_constrained_enthalpy_width_interface_2026-09-05.json`.
- Reconstruct the same measured pipe/HEM source and ambient state.
- Retain `thermal_species_ratio=1` and
  `mechanical_work=reduced_buoyancy_work`.
- Retain the free planar section.  Ground interaction and measured yaw are not
  combined with this stage.
- The six primary transported values are total mass flux, hydrogen mass flux,
  the two laboratory momentum components, total-energy flux, and the physical
  enthalpy second moment.

## State inversion and march

At every Runge--Kutta stage, recover a physical section from the six primary
values.  Momentum direction is `atan2(Pz, Px)`.  Hydrogen centre density is
eliminated using the analytic hydrogen-flux integral.  Solve the remaining
four bounded logarithmic variables (centre mixture density, section area,
excess velocity and thermal-width ratio) against total mass, momentum
magnitude, total energy and enthalpy second moment.  The permitted thermal
width remains `0.5 < beta_H < 2`; no state is clipped into the domain.

Use the phase-partitioned section fluxes and their established analytic
Jacobian in the inverse.  Require the maximum component-scaled mismatch of all
six reconstructed values to be at most `1e-8`.  A failed inverse or failed
local closure rejects the step rather than extrapolating it.

Advance the six values with classical fourth-order Runge--Kutta and advance
the centreline coordinates using `dx/ds=cos(theta)` and
`dz/ds=sin(theta)`.  The source vector is exactly the existing five-source
ledger plus its existing thermal-second-moment rate.

## Frozen numerical decisions

- Reach the first Trial 10 profile station at x = 1.78 m.
- Coarse maximum arc-length step: 0.005 m.
- Refined maximum arc-length step: 0.0025 m.
- Shorten the final step to approach the target; do not use a receptor
  interpolation to conceal an overshoot.
- Stop at the first inversion, phase, positive-diffusion, inflow, curvature or
  weak-budget failure.
- Maximum accepted-step budget: 4000.

## Acceptance and decisions

- Replay the frozen initial boundary and source before marching.
- The maximum component-scaled terminal balance error for the six directly
  integrated quantities must be at most `1e-5` under an independent order-16
  flux reconstruction.
- Coarse and refined terminal physical parameters and six independently
  reconstructed quantities must each agree within `0.005`.
- Record inverse residuals, weak residuals, sampled diffusivity minima,
  pointwise edge heat defect and every stop reason.
- Hash the inputs, code, tests and pre-registration before and after the run.
- Do not score measured receptors, promote a model default, or release a new
  package unless all preceding gates pass.  A failure remains a reported
  domain/numerical result and is not repaired by fitted physics.
