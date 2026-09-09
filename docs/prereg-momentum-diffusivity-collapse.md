# Pre-registration: Trial 10 momentum-diffusivity collapse diagnosis

Date frozen: 2026-09-10, after the stopped 4 m observation attempt and before
sampling the radial transport fields at its limiting section.

## Question

The Trial 10 direct-flux path stops near x=3.346 m because the inferred
momentum diffusivity becomes negative at an RK stage.  Is this a numerical
flux-inverse failure, a global loss of positive momentum transport, or a
localized radial incompatibility between the transported mean profiles and
the algebraic shear-stress closure?

## Frozen calculation

- Use the last accepted Trial 10 coarse section from
  `thermal_moment_observation_validation_first_attempt_2026-09-10.json`.
- Reconstruct the exact accepted closure and the first half-stage of the
  rejected 0.005 m step.  Do not alter a coefficient, flux or state.
- Retain `D_h/D_C=1`, reduced buoyancy work, order-8 integration, 257 phase
  probes and the default linear phase table.
- Sample 2,049 uniform radial-q points plus every phase-cell midpoint.
- At the minimum of species and momentum diffusivity record q, normalized
  radius `sqrt(2q)`, physical lateral and normal distances, affine origin and
  thermal-width-rate response, stress, velocity gradient and density.
- Record inverse residual, weak residual, thermal-width rate and every local
  validity gate at the accepted and rejected-stage sections.

This diagnostic cannot change the failed observation decision.  A localized
negative eddy diffusivity indicates a model-form incompatibility and motivates
an independently transported stress/TKE closure; it is not permission to clip
the diffusivity.  A uniform sign loss instead points to the bulk momentum
source or mean-profile evolution.  No default change or field score is made.
