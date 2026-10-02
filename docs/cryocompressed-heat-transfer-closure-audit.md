# Cryo-compressed heat-transfer closure audit

## Decision

DEGALI does **not** silently select an internal tank or discharge-pipe heat-
transfer coefficient for cryo-compressed hydrogen blowdown.  The available
`TransientTankWall` and `TransientPipeWall` models therefore require explicit
inner and outer coefficients.  The simpler prescribed-`UA` paths remain
explicit sensitivities only.

This is a scope guard, not an omission of heat transfer.  A wrong universal
coefficient can be more misleading than a declared sensitivity because it
changes both the blowdown mass flux and the onset of the two-phase tank state.

## What the public literature establishes

Cirrone et al. (2023) model E3 cryo-compressed releases with real-gas
properties, transient wall conduction, and a separate discharge-line energy
balance.  Its pipe equation uses a convection coefficient, pipe inner area,
and the inner wall temperature; the latter is either resolved by one-
dimensional conduction or deliberately prescribed in the simplified case.
That supports DEGALI's separation of the tank and pipe thermal boundaries.

Molkov, Dadashzadeh and Makarov (2019) derive a changing tank coefficient
from Nusselt correlations and a characteristic in-tank velocity, rather than
from inlet velocity alone.  However, their published property domain is
270--350 K and 0.1--77 MPa, with validation for 29--74 L high-pressure gas
filling tanks.  The paper explicitly warns that a coefficient depends on tank
thermal properties, dimensions, and flow conditions.  It cannot identify a
coefficient for the 2.815 L, LN2-immersed, rapidly discharging 80 K E3.1
vessel.

The same paper reports that constant-coefficient models can overpredict tank
temperature by 8--12 C or more in other filling cases.  That is direct reason
not to turn its correlation into a hidden 80 K blowdown default.

## Consequence for E3.1

The published E3.1 geometry and bath condition can define wall material,
external temperature, external coefficient, and pipe dimensions.  They do
not provide the evolving inner tank-wall or pipe-wall temperature, nor a
measured in-tank circulation field.  The E3.1 T4 signal is a welded flow
thermocouple, not a wall boundary.  It must remain a diagnostic.

An admissible study may therefore run a declared envelope of inner
coefficients and initial wall temperatures, report the resulting source
history spread, and retain the no-heat path as a lower-layer comparison.  It
must not tune an inner coefficient to the pressure trace and then call the
result a predictive validation.

## Inputs needed to promote a closure

1. A cryogenic blowdown test with time-resolved tank inner-wall or liner
   temperature, pipe-wall temperature, and storage pressure.
2. Vessel geometry and orientation sufficient to define the characteristic
   natural/forced/combined convection length scale.
3. A stated source of in-tank circulation or a separately validated
   characteristic-velocity model for the release regime.
4. An independent hold-out pressure/temperature trace to evaluate the closure
   without reselecting its coefficient.

## Sources

- Cirrone, D., Makarov, D., Kashkarov, S., Friedrich, A., & Molkov, V.
  (2023). *Physical model of non-adiabatic blowdown of cryo-compressed
  hydrogen storage tanks*. International Journal of Hydrogen Energy, 48(90),
  35387--35406. https://doi.org/10.1016/j.ijhydene.2023.05.182
- Molkov, V., Dadashzadeh, M., & Makarov, D. (2019). *Physical model of
  onboard hydrogen storage tank thermal behaviour during fuelling*.
  International Journal of Hydrogen Energy, 44(8), 4374--4384.
  https://doi.org/10.1016/j.ijhydene.2018.12.115
