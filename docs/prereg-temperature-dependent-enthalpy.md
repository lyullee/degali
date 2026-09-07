# Pre-registration: temperature-dependent ideal-gas enthalpy

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The independent mass-fraction/temperature profiles reproduce the two Raman
half-width slopes within 5.1%, but their centreline warming is 38.9% too slow.
The remaining energy model evaluates both component enthalpies as `cp(295 K) T`
over roughly 43--295 K. For normal hydrogen, rotational heat capacity changes
substantially across this interval. A single ambient-temperature heat capacity
overstates the energy required to warm cold hydrogen.

## Candidate fixed before running

Retain exactly the pre-registered profiles and transport ratios
`lambda_Y=1.16`, `Pr_t=0.85`, `Sc_t=0.70` and
`lambda_T=1.0526828470`. Replace only the caloric equation:

```
h_mix(T,Y) = Y h_H2^0(T) + (1-Y) h_air^0(T)
```

where each ideal-gas component enthalpy is evaluated from the installed
CoolProp Helmholtz ideal term. CoolProp's hydrogen implementation cites the
reference equation of state of Leachman et al., J. Phys. Chem. Ref. Data 38
(2009) 721--748. Tables cover 14.1--400 K and are interpolated monotonically.
CoolProp does not expose the Air ideal term below its first valid point near
25 K; that unused leading table segment is linearly extended with the first
two valid points solely to keep nonlinear trial states evaluable. Acceptance
requires every final profile temperature to remain above the extension.
Only enthalpy differences matter; the same component reference states are
used at the source, every radial point and the ambient subtraction, while
hydrogen species mass is conserved. For numerical conditioning, each
component is shifted to zero at ambient temperature; this cannot alter either
balance because total and hydrogen mass fluxes are separately conserved.

Density remains the same ideal-mixture equation of state. Equilibrium air
condensation remains disabled so the test changes one assumption only. The
four-flux `entrained_mass` establishment is mandatory. No coefficient or
property value is inferred from Raman data, and all existing defaults remain
unchanged.

## Decision rule

1. Boundary mass, species, momentum and energy residuals below `1e-8` for all
   nine sources; downstream species and energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients within 33--64 (mass) and 21--49 (temperature).
4. Every sampled density, mass fraction and temperature is physical.
5. A 161-point radial rerun must change every slope by less than 0.2%.

Failure of any criterion rejects the candidate. No heat-capacity multiplier,
enthalpy offset, profile width or virtual origin will be fitted.

---

## RESULTS

Rejected before regression by criteria 4 and 5. The enthalpy reference-offset
invariance test passed, but the four-flux boundary selected centreline
temperatures of 9.8--22.8 K for eight solvable cases; case 1 produced a
non-finite boundary Jacobian. These states are below the approximately 25 K
first valid Air ideal-gas table point and several are below hydrogen's 13.96 K
triple point. Six subsequent integrations also reached a singular flux
Jacobian before 0.13 m.

Starting from the constant-`cp` conservative solution and continuing to the
temperature-dependent property did not find another branch. It returned the
same low-temperature roots, with boundary Jacobian condition numbers of about
`0.76e12`--`1.40e12`. Thus the failure is not repaired by an enthalpy offset,
optimizer tolerance or initial guess.

The diagnostic also exposes the earlier assumption that matters first: even
the constant-`cp` conservative establishment has centre temperatures of only
8.7--17.9 K, despite 27--46 K ambient-plane sources. A physically admissible
initial entrainment/heating zone is required before refining gas heat
capacity. The off-default property-table implementation remains tested, but
this failed variant is excluded from the aggregate `compare()` path.
