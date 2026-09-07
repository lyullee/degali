# Pre-registration: differential thermal/species spreading

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Structural diagnosis

Hecht--Panda report a temperature half-width divided by hydrogen-mass
half-width of `0.06188 / 0.06503 = 0.9516`. The published integral model uses
one Gaussian spreading ratio for both thermodynamic and species scalars. Its
corresponding ratio is about 1.23--1.26, and remains above one after a fully
conservative establishment and after distributed equilibrium N2/O2
condensation. This sign-stable discrepancy cannot be removed by quadrature
refinement or a latent-heat multiplier.

Gradient-diffusion closures transport heat with `nu_t / Pr_t` and species with
`nu_t / Sc_t`. Cirrone, Makarov and Molkov, *Thermal radiation from cryogenic
hydrogen jet fires*, IJHE 44 (2019) 8874--8885, equations (3)--(4), use
`Pr_t=0.85` and `Sc_t=0.70` for cryogenic hydrogen. These values and the
following conversion are fixed independently of the Raman result.

## Candidate fixed before running

Retain the published velocity and hydrogen-species Gaussian profiles and its
species spreading ratio `lambda_Y=1.16`. Give the thermodynamic density
contrast its own Gaussian width

```
lambda_h = lambda_Y * sqrt(Sc_t / Pr_t)
         = 1.16 * sqrt(0.70 / 0.85)
         = 1.0528.
```

The square-root mapping follows Gaussian diffusion width proportional to the
square root of scalar diffusivity. Temperature remains derived from local
pressure, mixture molecular weight and density; no temperature data enter the
profile. The plug-to-Gaussian boundary uses the four-flux `entrained_mass`
solve, conserving total mass, hydrogen species, momentum and total enthalpy
plus kinetic energy. Distributed condensation is disabled so this test
isolates differential turbulent transport.

The existing one-ratio model is recovered exactly when the optional thermal
ratio is omitted or equals `lambda_Y`. The production JETPLU path and all
default parameters remain unchanged.

## Decision rule

1. Boundary mass, species, momentum and energy residuals below `1e-8` for all
   nine sources; downstream species and energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients within the printed case ranges 33--64 (mass)
   and 21--49 (temperature).
4. Every sampled density, mass fraction and temperature is physical.
5. A 161-point radial rerun must change every slope by less than 0.2%.

Failure of any criterion rejects this closure. Neither `Pr_t`, `Sc_t`, their
ratio, the Gaussian exponent nor a virtual origin will be fitted.

---

## RESULTS

Rejected. The 121-point result over all 549 frozen samples is:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21601 | 0.06079 | 0.01977 | 0.07719 | 57.73 | 34.86 | 3/4 |

Relative errors are -17.7%, -6.5%, **-30.1%**, and +24.7%. The centreline
temperature criterion fails. More importantly, narrowing the density-contrast
profile makes the physical temperature half-width slightly *wider* than the
single-ratio conservative model (0.07719 versus 0.07641). Density is a coupled
result of composition and enthalpy, so assigning it the heat-diffusion width is
not equivalent to transporting enthalpy.

All nine boundaries close below `1.7e-15`; maximum downstream species and
energy drift are `3.3e-8` and `1.2e-7`. Sampled mass fractions remain in
`[0, 0.155]`, density stays positive, and temperature remains within
116.9--294.9 K. A 161-point rerun gives 0.21604, 0.06079, 0.01977 and 0.07719;
each slope changes by less than 0.01%. The rejection is physical rather than a
solver or quadrature failure.

The next admissible structural test must transport mixture enthalpy itself as
an independent Gaussian scalar and obtain density from enthalpy, composition,
pressure and the equation of state. This rejected density-width shortcut will
not be used in the production path.

### Final-journal source correction

The final article revises the mass half-width to 0.07069 while leaving the
temperature half-width at 0.06188. The observed width ratio is therefore
0.8754, not the 0.9516 stated in the frozen rationale. This strengthens the
case for distinct scalar/thermal transport, but does not rescue this rejected
density-width shortcut or its failed centreline-temperature metric. See
`hecht-panda-journal-benchmark-correction.md`.
