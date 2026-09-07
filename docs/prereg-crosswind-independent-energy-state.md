# Pre-registration: independent-energy crosswind Gaussian plume

Date frozen: 2026-09-05, after the four-flux/local-shear diagnosis and before
calculating any result from the new state system.

## Why another state is required

The current JETPLU state has four physical unknowns: centre contaminant mass
concentration, area-width product, trajectory angle and excess axial velocity.
They are fixed by total mass, H2 mass and the two momentum balances. Mixture
density, temperature and enthalpy are then all read from one adiabatic table as
functions of that single concentration. Total energy can be audited at a
handoff, but it is not transported downstream.

That restriction is now observable rather than hypothetical. The conservative
four-flux near-field boundary closes source total mass, H2 mass, momentum and
energy below `1.4e-15`, yet two of seven PRESLHY sources cannot pass the existing
JETPLU temperature/profile screen at any tested boundary from 10D through 20D.
Moving the boundary makes both mismatches worse. Width tuning and handoff
location are therefore excluded as remedies.

HyRAM+ provides the coefficient-free structural precedent. Its unignited
Gaussian model transports centre velocity, width, centre density, centre fuel
mass fraction and trajectory angle, and solves continuity, species, two
momentum components and total energy. The present candidate adopts those two
independent thermodynamic/species centre states while retaining the crosswind
and atmospheric terms already validated in JETPLU.

## Frozen state and profiles

The seven integrated variables are

`q = (rho_c, Y_c, sysz, theta, u_c, x, z)`.

The first five are determined by five conservation equations; `x` and `z`
follow `dx/ds = cos(theta)` and `dz/ds = sin(theta)`. `sysz = sigma_y sigma_z`
retains JETPLU's elliptic cross-section and its existing atmospheric split.
No new adjustable state or field-fitted coefficient is introduced.

On JETPLU's existing finite Gaussian domain, let `G` be its scalar Gaussian
and let `G_u = G**lambda**2`, with the already published `lambda=1.16`. The
profiles are fixed as

```
u_s   = u_a cos(theta) + u_c G_u
rho   = rho_a + (rho_c - rho_a) G
rho Y = rho_c Y_c G
```

Thus density and H2 mass density follow the HyRAM+ profile form independently;
the local mass fraction is `Y=(rho Y)/rho`. Temperature, phase inventories and
enthalpy are obtained from the existing density-plus-composition equilibrium
property closure, including N2/O2/H2O phase change and component `h(T)`. They
are not taken from a one-dimensional concentration table.

## Frozen balances

The five section integrals are total mass, H2 mass, horizontal momentum,
vertical momentum and total energy relative to ambient enthalpy:

```
M    = integral rho u_s dA
M_H2 = integral rho Y u_s dA
P_x  = cos(theta) integral rho u_s**2 dA
P_z  = sin(theta) integral rho u_s**2 dA
E    = integral [rho u_s (h-h_a) + 0.5 rho u_s**3] dA
```

The right-hand sides retain JETPLU's local density-scaled shear entrainment,
cross-flow entrainment, ambient-spread growth, buoyancy and perpendicular form
drag. Source-momentum entrainment ends at the near-field boundary, as supported
by the earlier seven-trial mechanism test. Entrained air supplies total mass,
horizontal momentum and kinetic energy at the local ambient-wind state. It
supplies no H2. Form drag only rotates the axial momentum vector and therefore
does no axial work; no empirical drag-energy term is added. Radiation remains
off for this comparison.

The five state derivatives are obtained from the numerical Jacobian of these
five physical fluxes, following the independently tested conserved
axisymmetric implementation. Existing JETPLU equations and six-state behavior
remain unchanged; this is a separate research path.

## Handoff and numerical gates

- Use the same seven condition-selected PRESLHY trials: 10, 11, 12, 22, 23,
  24 and 25.
- Use four-flux `entrained_mass` near-field establishment and the fixed 10D
  handoff. Do not search handoff position in the primary test.
- At the boundary, solve all five physical states against the five near-field
  fluxes. Every relative flux residual must be below `1e-8`.
- Retain the existing independent profile screens: H2 half-width within 5%
  and centre temperature within 2 K. All seven interfaces must pass.
- Retain adaptive 32/64/128/256/512-point energy quadrature and the unchanged
  `1e-5` quadrature-convergence gate.
- Reject non-positive density, width or excess velocity, `Y_c` outside
  `(0,1]`, angles outside `(-pi/2,pi/2)`, or a singular/ill-conditioned flux
  Jacobian.
- Unit tests must first prove source-plane five-flux closure, finite
  derivatives, monotone H2 dilution and step-refined downstream balance
  closure before any field score is used.

## Frozen field comparison and decision

Use the same reduced JSON, atmospheric inputs, sensor coordinates, downstream
arc restriction, vertical fits, solver tolerances and corrected JETPLU
baseline as `prereg-preslhy-coupled-crosswind-validation.md`. No observation
may participate in the interface solve or set a coefficient.

Promotion requires all seven interfaces and, on exactly common keys:

1. lower `abs(log(MG))` than the corrected baseline;
2. lower VG;
3. FAC2 no lower than baseline;
4. mean model/measured `sigma_z` closer to one; and
5. centre-height MAE no larger.

The result will be reported even if rejected. A successful interface test alone
establishes that the missing energy degree of freedom was diagnosed correctly;
it does not by itself validate downstream dispersion.

## Interface result

The boundary hypothesis is supported: all seven fixed 10D interfaces pass.
The maximum relative five-flux residual is `6.45e-16`; maximum energy-
quadrature change is `8.37e-6`; maximum H2-width mismatch is 4.761%; and maximum
centre-temperature mismatch is 1.990 K. Trials 10 and 25, which the single-
scalar JETPLU interface rejected, now pass without moving the boundary or
relaxing a gate. See `preslhy-independent-energy-interface-results.md`.

Downstream promotion remains untested. Per the frozen sequence, the next work
is the energy ODE and its unit-level balance/convergence tests, followed by the
seven-trial field comparison.

## Downstream numerical addendum

Frozen before calculating any downstream sensor prediction:

- Integrate the transformed positive state with deterministic fourth-order
  Runge--Kutta and a primary maximum arc step of 0.02 m. Repeat at 0.01 m as a
  numerical sensitivity; do not choose between them from agreement with data.
- Use the adaptive interface quadrature order downstream: 64 points for trials
  10--24 and 128 for trial 25 in the present boundary result.
- Integrate only far enough to cover the last selected concentration sensor or
  vertical fit, plus a fixed 0.25 m margin. First-order downstream equations are
  causal, so calculating the unused remainder to 40 m cannot change a sensor
  prediction.
- Evaluate receptors with the same direct-plus-ground-image scalar Gaussian as
  the existing validation. Apply that shape to both density excess and H2 mass
  density before converting local mass fraction to mole fraction. Fail closed
  if the superposition produces a non-physical mixture.
- Report the direct section-flux minus integrated-source balance. The field
  result is not eligible for promotion if state positivity fails, if any
  integration fails, or if halving the step materially changes the aggregate
  decision metrics.

### Numerical audit correction before aggregate scoring

The first trial-10 pilot integrated transformed physical states directly. It
remained positive but accumulated a 4.431% flux-minus-source residual at the
frozen 0.02 m step. No aggregate statistic was calculated or inspected. That
integrator is rejected on its internal conservation result.

The primary RK4 step and every physical source term remain unchanged. The
accepted numerical form advances the five conserved section fluxes directly
and inverts them back to positive physical states at each RK stage. This makes
the discretised mass, species, vector-momentum and energy equations the
integrated variables rather than relying on a noisy finite-difference flux
Jacobian over piecewise phase tables. The state-Jacobian implementation remains
available as a diagnostic. Field scoring proceeds only with the conservative-
flux form and still requires the frozen 0.02/0.01 m comparison.

## Downstream result

All seven interfaces and integrations pass. Primary 0.02 m results are MG
1.118, VG 1.176, FAC2 0.976, mean model/measured `sigma_z` 1.066 and centre-
height MAE 0.078 m. The corrected baseline on identical keys is 1.074, 1.213,
0.952, 1.091 and 0.054 m, respectively. The 0.01 m repeat gives the same
reported concentration metrics and geometry 1.066/0.078 m; the precise width
ratio changes by `6.9e-5` and centre MAE by `2.7e-6 m`.

The candidate is not promoted because `abs(log(MG))` and centre-height MAE
worsen. VG, FAC2 and width improve, and the maximum long-range balance residual
is only `1.04e-7`, so this rejection is physical rather than numerical. The
independent-energy state remains the structural base for the next buoyant-rise
audit.
