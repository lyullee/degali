# Pre-registration: distributed equilibrium air condensation

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

Both coefficient-free conservative Gaussian boundaries close mass, species,
momentum and total energy but underpredict the Hecht--Panda centreline
temperature slope by about 31%. Their common energy model treats all entrained
air as gas. At 40--80 K that is not an admissible equilibrium state: N2 and O2
condense/freeze, releasing latent heat. Li et al. (2026) identify the same
opposing density/latent-heat effects for cryogenic hydrogen jets.

## Candidate fixed before running

Use the scalar-peak conservative boundary, which already closes the three
plug invariants without assuming an integrated development-zone air mass.
Within every Gaussian radial quadrature point:

1. retain the prescribed total H2 and dry-air mass densities;
2. split dry air into N2/O2 at atmospheric composition;
3. solve temperature from bulk specific volume, including gas volume and
   condensed-phase volume;
4. enforce component partial-pressure equilibrium with the stable
   solid-vapour relation below each triple point and liquid-vapour relation
   above it;
5. include gas and condensed enthalpy in the single total-energy flux, with
   phase-change enthalpy exactly once.

Condensate shares the gas axial velocity. Argon, water/CO2 frost and finite
phase-change kinetics are omitted. Therefore this is explicitly the
instantaneous-equilibrium/no-slip bound, not a fitted production closure.
No Raman quantity sets a coefficient or phase parameter.

## Decision rule

1. Boundary species, momentum and phase-aware energy residuals below `1e-8`
   for all nine sources; downstream species and energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients within 33--64 (mass) and 21--49 (temperature).
4. Temperature and condensed fractions must remain physical at every sampled
   point.
5. Failure of any criterion rejects the candidate. No latent multiplier,
   particle diameter, virtual origin or profile ratio will be tuned.

---

## RESULTS

Rejected. The 81-point radial quadrature result over all 549 samples is:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21490 | 0.05976 | 0.01951 | 0.07515 | 59.9 | 36.9 | 3/4 |

Relative errors are -18.2%, -8.1%, **-31.0%**, and +21.4%. The distributed
phase balance moves the failed centreline-temperature error by only 0.3
percentage point from the scalar-constrained all-gas result. Thus criterion 2
fails and the option is not adopted.

The sign is physically informative: latent heat raises the just-established
temperature slightly, but subsequent re-evaporation consumes that heat. An
instantaneous equilibrium/no-slip closure cannot reproduce the missing
warming. A finite-rate two-temperature gas/particle model would require
particle number/size or interfacial area, already identified as unavailable
in `finite-rate-sublimation-audit.md`; no such value is fitted here.

A 121-point radial rerun gives 0.21496, 0.05976, 0.01951 and 0.07514; every
slope changes by less than 0.03%. The rejection is therefore not a coarse
quadrature artefact.
