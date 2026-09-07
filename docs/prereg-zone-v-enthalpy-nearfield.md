# Pre-registration: Li equation-35 enthalpy balance before 10D

Date: 2026-09-05

## Defect isolated before calculation

Li et al. (2026) define the Zone-V energy equation as

`d/ds integral rho*v*(h-h_amb) dA = 0`  (equation 35).

The present independent-energy comparison selected Li's enthalpy transport
only at the 10D crosswind interface and downstream.  The upstream conserved
Gaussian march from the Station-4 boundary to 10D still transports

`rho*v*(h-h_amb) + 0.5*rho*v^3`.

It therefore follows the HyRAM+ total-energy convention during the first 3.8D
of Zone V, then switches to Li's equation only after the fixed handoff.  That
is a physically inconsistent hybrid.  The missing comparison is to apply
Li equation 35 from the beginning of established flow, where the high velocity
makes the distinction largest.  It tests whether resolved kinetic-energy
loss has been converted into internal heating too early, contributing to the
warm/light 10D state and excessive downstream buoyancy.

## Frozen implementation

- Add `energy_transport="total"|"enthalpy"` to the conserved Gaussian jet.
  `total` exactly preserves every existing result and remains the default.
- For `enthalpy`, use the already implemented component/phase enthalpy
  integral as the fifth establishment target and ODE flux.  Entrained ambient
  air has zero ambient-relative enthalpy; radiation, if selected, remains an
  external enthalpy source.  Do not add or fit a turbulent-dissipation term.
- Use `source_flux` Zone-IV mapping, equilibrium N2/O2/H2O thermodynamics,
  normal hydrogen, local-shear crosswind entrainment, corrected Houf velocity
  width and Li enthalpy transport downstream.  All coefficients, the 10D
  handoff and numerical thresholds remain fixed.

## Test order and gates

1. Unit tests must show that `total` is unchanged and both forms conserve
   their selected energy flux below `2e-4`.
2. Run trial 10 through the 10D interface.  Require four establishment and
   five interface residuals below `1e-8`, energy quadrature below `1e-5`, H2
   width mismatch below 5% and centre-temperature mismatch below 2 K.
3. Before downstream scoring, report Station-4 and 10D temperature, density,
   total mass, enthalpy flux and resolved kinetic-energy flux for both energy
   forms.  The hypothesis direction is a colder/denser and less-buoyant 10D
   state under equation 35.
4. Only if the interface passes and trial-10 centre-height absolute error falls
   at 1.78, 4 and 6 m without a greater-than-10-percentage-point width loss at
   6 m may the seven-trial field run proceed.
5. Failure retains the option only for equation-level reproduction.  No
   partial field score or post-result dissipation coefficient may be used.

## Result

Rejected at the first boundary gate; no trial-10 or field score was run.

On the existing representative measured-throat regression source, the
`source_flux` plug-to-Gaussian solve cannot simultaneously conserve total
mass, H2 species, axial momentum and Li's enthalpy flux.  Its final normalized
residuals are `+33.19%` mass, `-58.41%` species, `-2.69e-8` momentum and
`-27.46%` enthalpy, far outside the fixed `1e-8` limit.  This is not an ODE
or quadrature error: the four target integrals are incompatible with the
four-parameter Gaussian state under this boundary construction.

The option is retained solely so equation 35 can be reproduced explicitly,
and the regression requires this incompatible conservative boundary to fail
closed.  Applying it only to the lower-velocity PRESLHY source after observing
its result would violate the applicability and test order.  A viable
replacement needs an additional turbulent-kinetic-energy/profile state across
Zone IV, not a switch that discards resolved kinetic energy while preserving
the same four Gaussian unknowns.
