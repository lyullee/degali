# Pre-registration: entraining conservative Gaussian establishment

**Frozen before the first result from this candidate.** Do not edit above the
`RESULTS` line after a candidate result is known.

## Candidate fixed before running

The fixed-centreline-velocity candidate has been falsified. For the same
published development length, calculate the air added by the already frozen
source-momentum entrainment flux,

`m_air = rho_ambient E_mom S_development`.

The established total mass target is the plug mass plus this air. Solve the
four Gaussian unknowns -- centreline velocity, width, centreline density and
centreline fuel mass fraction -- against four direct cross-sectional fluxes:

1. total mass, including the predicted entrained air,
2. hydrogen species,
3. axial momentum, and
4. total enthalpy plus kinetic energy relative to ambient.

Position, angle, development length, profile equations and all entrainment
constants remain unchanged. No observation or adjustable coefficient enters
the boundary solve.

## Predictions and decision rule

1. All four normalized establishment residuals must be below `1e-8` for all
   nine sources, with positive velocity, width, density and mass fraction.
2. All four Raman slopes must remain within the frozen 25% band.
3. Median radial coefficients must be inside 33--64 for mass and 21--49 for
   temperature. Boundary rounding is reported, not exploited.
4. Species and energy flux drift through the ODE must remain below `2e-4` in
   the regression test.
5. If the nonlinear system has no positive solution or any rule fails, retain
   the published establishment model and report its boundary loss explicitly.
6. The DEGADIS/Fortran and atmospheric cross-wind defaults remain unchanged.

---

## RESULTS

Rejected. All nine boundary systems found positive solutions and closed all
four normalized flux residuals below `5e-12`. The Raman result was:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21663 | 0.06136 | 0.01947 | 0.07641 | 56.6 | 35.6 | 3/4 |

The relative errors are -17.5%, -5.6%, **-31.1%**, and +23.5%. Thus rule 2
fails on centreline temperature despite exact boundary and downstream
conservation. The source-momentum flux integrated as constant over the whole
development length admits too much air: in the first case it is about 2.8
times the plug source mass. No length or multiplier is adjusted to repair it.

The published scalar-peak relationship supplies a non-fitted alternative
fourth constraint without prescribing the incompatible plug centreline
velocity or an assumed development-zone mass integral. It is tested under a
new freeze in `prereg-scalar-constrained-establishment.md`.
