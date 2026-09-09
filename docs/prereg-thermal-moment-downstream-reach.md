# Pre-registration: guarded downstream reach of the conservative thermal moment state

Date frozen: 2026-09-08, before the reach calculation.

## Question

Can the already implemented six-moment ambient-reservoir closure leave the
millimetre verification interval and reach the first Trial 10 profile station
without violating its own phase, positive-diffusion, inflow, curvature or
conservation conditions?

This stage does not select a new diffusivity and does not fit temperature,
hydrogen, wind direction or a thermal-width rate. It adds a guarded adaptive
marching interface around the existing local conservative operator. A failed
reach is retained as a model-domain result, not bypassed by extrapolation.

## Frozen model and initial condition

- Trial 10 only.
- Start from the accepted six-moment boundary in
  `buoyancy_constrained_enthalpy_width_interface_2026-09-05.json`.
- Reconstruct the original measured pipe/HEM source and ambient state.
- Use the existing natural-flux ambient-reservoir weak closure with
  `thermal_species_ratio=1` and
  `mechanical_work=reduced_buoyancy_work`.
- The ratio one is the previously tested equal eddy-diffusion reference, not a
  value inferred from the new five-sensor profiles.
- Retain the free planar section. Yaw and ground interaction are excluded
  because the existing thermal-moment weak geometry has no validated closure
  for them.

## Guarded adaptive march

Integrate the eight log/trajectory/thermal-width parameters together with the
five cumulative conserved sources and the cumulative physical thermal second
moment. Use explicit midpoint step doubling: compare one full step with two
half steps and accept the two-half result only when their maximum scaled
parameter difference is at most `2e-5`.

- initial trial step: 0.002 m;
- maximum step: 0.005 m;
- minimum step: 1e-6 m;
- maximum accepted steps: 4000;
- target: the first observed vertical cross at x=1.78 m, with an arc-length
  ceiling of 1.25 m from the boundary;
- accepted steps may grow by at most 1.5 and rejected steps are halved;
- every right-hand-side evaluation must pass the existing weak residual,
  positive sampled species/momentum diffusivity, inward mass flux and curvature
  gates. No negative value is clipped.

Stop at the last accepted state if a minimum step still cannot produce a valid
or accurate update. Record the exception, attempted step and position.

## Integrity and decisions

- Hash all raw inputs, stored boundary, executable code and tests before and
  after the calculation.
- Replay the frozen boundary moments to `1e-9` and the original stored source
  to `1e-8` before marching.
- At the final accepted state, compare the independently recomputed five
  fluxes plus physical thermal second moment with initial values plus cumulative
  sources. The maximum component-scaled balance residual must be <= `5e-4`.
  This relaxed reach-screen limit does not supersede the earlier short-segment
  `1e-5` verification.
- If x=1.78 m is reached, run a second march with half the initial and maximum
  step and require terminal parameters and all six recomputed moment/flux values
  to agree within 0.005 before a separate receptor implementation is allowed.
- If it is not reached, do not calculate a sensor score or promote the model.
  The first failed internal gate becomes the next closure target.
- No existing default, stored calculation or public release is changed.
