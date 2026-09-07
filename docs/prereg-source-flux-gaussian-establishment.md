# Pre-registration: source-flux Gaussian establishment bound

Date: 2026-09-05

## Physical defect isolated before calculation

The PRESLHY axisymmetric source is already the air-loaded Station-3 endpoint
of the flash/evaporation and initial-heating calculation. Li et al. (2026)
place momentum entrainment `E_mom` in Zone III (their equations 22 and 24),
then define Zone IV separately by its length and plug-to-Gaussian profile
relations (equations 25--28). The established-flow mass source begins in Zone
V (equations 31 and 43--48).

The current `entrained_mass` boundary instead adds

`rho_ambient * E_mom * (S4-S3)`

once more while mapping the already air-loaded source through Zone IV. That
extension is not stated in Li or the HyRAM source. It was introduced locally
to close four Gaussian fluxes, and the earlier Raman audit independently found
that it admits excessive air over the same empirical development length.

## Frozen coefficient-free bound

- Add a `source_flux` Gaussian-establishment option. Keep the published Zone-IV
  length and boundary location, but solve the four Gaussian unknowns against
  the Station-3 source-plane total mass, H2 species, axial momentum and total
  energy exactly. No ambient mass, co-flow momentum or co-flow kinetic energy
  is added inside Zone IV.
- Begin the existing `beta_A=0.28` momentum entrainment only after the Gaussian
  established-flow boundary. Preserve phase thermodynamics, component
  enthalpy, `lambda=1.16`, the corrected Houf velocity width, local-shear
  crosswind entrainment, the fixed 10D handoff and every numerical threshold.
- This is a zero-net-Zone-IV-entrainment lower bound, not permission to fit
  `beta_A` or a development-length multiplier. The existing `entrained_mass`
  path remains available for reproduction.

## Test order and unchanged gates

1. Run PRESLHY trial 10 only. Require four establishment residuals and five
   interface residuals below `1e-8`, energy quadrature below `1e-5`, H2-width
   mismatch below 5%, and centre-temperature mismatch below 2 K.
2. Report source, established-boundary and 10D total mass, temperature,
   density and section buoyancy before looking at downstream height.
3. If trial 10 passes, compare its 1.78/4/6 m centre-height and force budget.
   Then run the frozen seven interfaces. Only 7/7 permits the 42-arc and
   17-section field comparison.
4. A useful direction must reduce the corrected equilibrium candidate's
   trial-10 rise and seven-trial centre-height MAE without worsening
   `abs(log MG)`, VG, FAC2 or vertical-width error relative to the corrected
   JETPLU baseline. No partial field score may be used for promotion.
5. If the lower bound becomes excessively cold/dense or fails an interface,
   retain it only as a diagnostic. A finite Zone-IV entrainment closure then
   requires independent development-zone data; the coefficient must not be
   inferred from the PRESLHY trajectory.

## Result

The source-flux mapping is feasible and conservative. For trial 10 it retains
the Station-3 total mass `0.430660 kg/s` through the 6.2D profile transition,
with all four establishment residuals below `3e-16`. The previous
`entrained_mass` boundary added ambient mass there. After the remaining 3.8D
of established-flow integration, trial-10 mass at 10D falls from `1.265796`
to `0.748013 kg/s`, centre temperature from `64.315` to `46.967 K`, and
section buoyancy from `0.01723` to `0.01512 N/m`.

All seven fixed interfaces pass without relaxing a gate. The largest
five-flux residual is `1.58e-15`, quadrature change `9.13e-6`, H2-width
mismatch 4.019%, and centre-temperature mismatch 1.642 K. Trial 10's 6 m
centre height falls from 0.7649 to 0.7375 m and cumulative buoyancy from 5.695
to 5.239 N. This confirms that the extra Zone-IV mass was a real contributor,
but not the dominant downstream rise mechanism.

The full frozen field result is:

| model | MG | VG | FAC2 | mean sigma ratio | centre MAE, m |
|---|---:|---:|---:|---:|---:|
| corrected JETPLU baseline | 1.074 | 1.213 | 0.952 | 1.091 | 0.0542 |
| previous `entrained_mass` candidate | 1.062 | 1.163 | 0.976 | 0.976 | 0.0758 |
| `source_flux` bound | **1.012** | **1.162** | **0.976** | 0.899 | 0.0707 |

The source-flux bound materially improves concentration and reduces the
candidate centre-height error, but it is **not promoted**: its width deviation
from one is slightly larger than baseline and its centre-height MAE remains
31% higher than baseline. The result also shows why fitting `beta_A` in Zone
IV would be unsafe. By 6 m the local-shear established-flow model has largely
replaced the removed mass and warming, so another Zone-IV multiplier cannot
repair the remaining trajectory without trading away width.

`source_flux` remains an explicit research bound and `entrained_mass` remains
available for reproduction. Neither is a production default. The next audit
must address the downstream vertical-force/ground-interaction structure, not
fit a development-zone entrainment coefficient.

Machine-readable values are in
`reference/preslhy/source_flux_establishment_2026-09-05.json`.

## Regression verification

The boundary-option unit tests pass. The final complete non-slow regression
passes `114` tests (`127` dependency/data skips and `58` slow tests
deselected), and the representative phase-profile crosswind handoff passes
separately. No compatibility threshold or existing default was changed.
