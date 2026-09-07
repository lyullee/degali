# Pre-registration: independent temperature/species profiles

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The rejected density-width split demonstrated that density cannot stand in
for the thermal scalar in a light, variable-composition gas. The RANS energy
equation used for cryogenic hydrogen by Cirrone, Makarov and Molkov contains
the turbulent heat flux `(mu_t cp / Pr_t) grad(T)`, while the species equation
contains `(mu_t / Sc_t) grad(Y)`. The integral model must therefore allow the
measured temperature and hydrogen-species profiles to have distinct widths,
while conserving the integrated energy rather than deriving temperature from
a prescribed density Gaussian.

## Candidate fixed before running

Use the published velocity and hydrogen mass-density shapes and an independent
temperature shape:

```
v       = v_c exp[-(r/B)^2]
rho Y   = rho_c Y_c exp[-(r/(lambda_Y B))^2]
T - T_a = (T_c - T_a) exp[-(r/(lambda_T B))^2]
```

The independently sourced values remain `lambda_Y=1.16`, `Pr_t=0.85` and
`Sc_t=0.70`, giving

```
lambda_T = lambda_Y sqrt(Sc_t / Pr_t) = 1.0526828470.
```

At every radius, hydrogen mass fraction and density are obtained from the
specified `rho Y`, temperature, ambient pressure, ideal-mixture molecular
weight and the equation of state. Mixture enthalpy is then `cp(Y) T`. The ODE
still integrates mass, horizontal/vertical momentum, hydrogen species and
total enthalpy plus kinetic energy fluxes. The four-flux conservative
`entrained_mass` establishment is mandatory.

No Raman slope or radial coefficient sets either width. Equilibrium air
condensation is disabled to isolate this turbulent-transport closure. The
existing density-profile model and production JETPLU defaults are unchanged.

## Decision rule

1. Boundary mass, species, momentum and energy residuals below `1e-8` for all
   nine sources; downstream species and energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients within 33--64 (mass) and 21--49 (temperature).
4. Every sampled density, mass fraction and temperature is physical.
5. A 161-point radial rerun must change every slope by less than 0.2%.

Failure of any criterion rejects the candidate. `Pr_t`, `Sc_t`, Gaussian
exponents, source state and virtual origin will not be fitted.

---

## RESULTS

Rejected before regression by criterion 4. The first source's published
establishment guess has `Y_c=0.872`; because `lambda_T < lambda_Y`, the
specified hydrogen mass density remains high at radii where temperature has
already recovered. At those points it exceeds the maximum pure-H2 density
allowed by pressure and temperature, so the equation-of-state denominator is
non-positive. The four-flux optimizer cannot even evaluate its starting
profile.

This is a structural incompatibility, not a convergence tolerance to relax.
The candidate is removed from executable options. Hecht--Panda fit *mass
fraction* and temperature as Gaussian quantities, whereas the inherited
integral model prescribes hydrogen mass density. A separate pre-registration
therefore tests independent `Y` and `T` profiles without changing this failed
candidate after seeing the result.
