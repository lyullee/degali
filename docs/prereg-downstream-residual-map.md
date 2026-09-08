# Pre-registration: downstream thermal and energy residual map

Date: 2026-09-08

## Purpose

Combine already sealed Trial-10/Trial-23 diagnostics into one distance-resolved
map before selecting another physical closure. No trajectory is reintegrated
and no coefficient is selected in this step.

## Frozen inputs

- `downstream_cold_envelope_2026-09-06.json`: pressure-loss-density sensor and
  section temperatures;
- `mechanical_thermal_scale_2026-09-06.json`: thermal flux, mean kinetic flux,
  accounted mechanical input and heat-capacity flux;
- `trial10_vertical_momentum_budget_2026-09-05.json`: Trial-10 centre density,
  vertical momentum and cumulative buoyancy;
- `joint_tke_normal_budget_2026-09-06.json`: independent initial TKE/normal
  work feasibility bounds.

Use only the `pressure_loss_density` variant and Trials 10 and 23. Match rows
by trial and requested distance, with a maximum saved-section distance error
of 0.011 m. Preserve minimum, p05 and median temperature observations as
separate statistics.

## Derived quantities

For each centreline thermocouple distance, record:

- predicted minus observed temperature;
- conditional specific-enthalpy gap from the sealed audit;
- that positive gap divided by local mean kinetic energy;
- cumulative accounted mechanical input divided by the magnitude of the
  initial thermal flux;
- the uniform sensible-temperature scale of the accounted mechanical energy.

For Trial 10, attach the nearest sealed cumulative-buoyancy value only when
its distance differs by at most 0.011 m. Do not interpolate a new force
history. Report the order-16 TKE/normal-work feasibility bounds separately;
they are necessary initial constraints, not measured initial conditions.

## Interpretation gates

1. If the positive enthalpy gap exceeds local mean kinetic energy by more than
   a factor of 10, direct mean-kinetic-energy thermalisation is rejected as the
   primary explanation at that point.
2. If the accounted mechanical-energy temperature scale is below 1 K, it is
   classified as thermally small relative to the observed multi-kelvin gap.
3. A change of sign between minimum and median temperature residuals is marked
   as an observation-statistic/time-response ambiguity, not averaged away.
4. No TKE value is adopted without an independent measured boundary.

