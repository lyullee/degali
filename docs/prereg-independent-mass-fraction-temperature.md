# Pre-registration: independent mass-fraction/temperature profiles

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

Hecht--Panda fit both hydrogen mass fraction and temperature with radial
Gaussians. Their measured temperature-to-mass half-width ratio is 0.9516,
whereas a single thermodynamic/species profile in the conserved integral model
gives about 1.25. A direct temperature Gaussian is also the self-similar
counterpart of the `(mu_t cp / Pr_t) grad(T)` term in the independently sourced
cryogenic-hydrogen RANS energy equation.

## Candidate fixed before running

Use independent, everywhere-admissible mass-fraction and temperature shapes:

```
v       = v_c exp[-(r/B)^2]
Y       = Y_c exp[-(r/(lambda_Y B))^2]
T - T_a = (T_c - T_a) exp[-(r/(lambda_T B))^2]
```

Set `lambda_Y=1.16` and, independently of the Raman result, use Cirrone et
al.'s `Pr_t=0.85`, `Sc_t=0.70`:

```
lambda_T = lambda_Y sqrt(Sc_t / Pr_t) = 1.0526828470.
```

Local mixture molecular weight, density and enthalpy follow from `Y`, `T`,
ambient pressure, the ideal-mixture equation of state and `cp(Y)`. The ODE
integrates total mass, two momentum components, hydrogen species and total
enthalpy plus kinetic energy. The plug boundary must use the four-flux
`entrained_mass` solve. No Raman result sets a coefficient.

Equilibrium air condensation is disabled to isolate differential turbulent
transport. Existing density-based and production JETPLU defaults are
unchanged.

## Decision rule

1. Boundary mass, species, momentum and energy residuals below `1e-8` for all
   nine sources; downstream species and energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients within 33--64 (mass) and 21--49 (temperature).
4. Every sampled density, mass fraction and temperature is physical.
5. A 161-point radial rerun must change every slope by less than 0.2%.

Failure of any criterion rejects the candidate. `Pr_t`, `Sc_t`, profile
exponents, source state and virtual origin will not be fitted.

---

## RESULTS

Rejected. The 121-point result over all 549 frozen samples is:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21645 | 0.06472 | 0.01727 | 0.05873 | 50.50 | 61.32 | 3/4 |

The two independently fixed widths repair both width slopes: relative errors
are -0.5% for mass and -5.1% for temperature. Mass centreline decay remains
within the frozen band at -17.6%. Centreline warming, however, worsens to
**-38.9%**, and the temperature radial coefficient 61.32 exceeds the printed
21--49 case range. Criteria 2 and 3 fail, so this constant-heat-capacity
candidate is not adopted.

All nine boundaries close below `2.8e-14`; maximum downstream species and
energy drift are `1.7e-7` and `1.6e-7`. Sampled mass fraction stays below
0.133, density remains positive and temperature lies within 75.9--295 K. A
161-point rerun gives 0.21647, 0.06472, 0.01727 and 0.05873; every slope changes
by less than 0.01%.

This result separates two errors that the original one-scalar model conflated:
independent turbulent diffusivities explain the observed radial widths, but
the axial energy-temperature conversion is still deficient. The next frozen
test retains these profiles and replaces ambient-temperature constant `cp*T`
with reference-quality temperature-dependent ideal-gas component enthalpies.
