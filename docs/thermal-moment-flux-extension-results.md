# Direct-flux thermal-moment extension results

Date: 2026-09-10 (Asia/Seoul)

## Scope

This stage tests whether the research-only independent-enthalpy-width model can
transport its six conserved section quantities from the stored PRESLHY
near-field interfaces to the first frozen downstream profile stations. It is a
numerical applicability test, not an observation fit.

The transported primary variables are total mass, H2 mass, horizontal and
vertical momentum, total energy, and the physical second transverse moment of
advective enthalpy. A bounded physical Gaussian section is reconstructed at
every Runge--Kutta stage. The calculation fixes `D_h/D_C=1` and the registered
reduced-buoyancy-work ledger; neither value was fitted to these trials.

## Resolved numerical defect

The inherited JETPLU width split enforces

```text
sigma_y sigma_z = A
sigma_y^2 - sigma_z^2 = sigma_ya^2 - sigma_za^2 = D
```

through a Brent solve with absolute tolerance `1e-4`. Near Trial 24, the
six-flux inverse requires much smaller width changes. The computed widths then
remain on a numerical plateau, so the actual M2 area derivative disagrees with
the registered Jacobian and the inverse cannot converge by reducing arc step.

For the independent Gaussian-enthalpy add-on only, DEGALI now uses the
cancellation-safe positive algebraic solution

```text
t = sigma_y^2
t^2 - D t - A^2 = 0
```

while preserving the existing spread-floor rule. The original/core DEGADIS
splitter is untouched. Randomized tests over six decades satisfy both geometry
identities to floating-point precision, and the retained Trial 24 M2 Jacobian
agrees across three independent difference sizes.

## Preregistered result

The stored six boundary moments were held fixed and reprojected onto the exact
width coordinates. Both 0.005 m and 0.0025 m nominal arc steps were marched to
the preregistered station, with inverse-only step rejection, strict endpoint
inversion, independent order-16 reconstruction, and all local physical gates.

| Trial | target x (m) | max parameter difference | max flux difference | min sampled diffusivity | max independent balance |
|---:|---:|---:|---:|---:|---:|
| 11 | 0.79 | 1.130e-3 | 7.959e-5 | 4.432 | 2.965e-11 |
| 12 | 1.78 | 8.639e-5 | 1.571e-4 | 0.958 | 5.046e-11 |
| 22 | 1.78 | 1.781e-5 | 3.377e-5 | 0.741 | 3.912e-11 |
| 23 | 0.79 | 4.065e-4 | 1.147e-4 | 3.586 | 2.578e-11 |
| 24 | 1.78 | 2.404e-5 | 4.423e-5 | 0.731 | 4.514e-11 |
| 25 | 1.78 | 1.887e-5 | 3.567e-5 | 0.772 | 4.040e-11 |

All six trials passed the frozen `0.005` convergence limit, reached their
targets at both step levels, retained positive sampled species/momentum mixing,
and passed independent balance. Trial 24 required guarded initial step
reduction; the reduced step later grew back under the preregistered eight-
success rule.

## Phase interpolation sensitivity

A separate nodally exact C1 bicubic-Hermite lookup, based on PCHIP nodal slopes,
was screened against 8,000 independently sampled exact property states. It
reduced temperature RMS lookup error from `1.1486e-4` to `4.1887e-5 K`, retained
negative `dH/drho` over 513,441 dense points, and reproduced table nodes. Its
Trial 24 pilot also reached x=0.43 m after the analytic width correction.

The six-trial result above deliberately uses the default linear lookup. C1 is
available only as an explicit research sensitivity and was not promoted.

## Interpretation and limitations

This closes the earlier numerical conservation/solvability gap for the six
selected boundaries. It does not show that predicted concentration,
temperature, trajectory or spread agrees with downstream measurements. No
sensor was scored and no parameter was calibrated in this stage.

The path remains outside a validated design domain for obstacles, ground-cut
sections, indoor releases, pool spreading, arbitrary jet directions, ignition,
combustion and explosion. Its unity thermal/species diffusivity and reduced
mechanical-work allocation remain research assumptions requiring independent
physical validation.

The immutable audit JSON and external experiment files remain local and are
excluded from distribution. Public tools record the calculation, provenance
requirements and acceptance logic without redistributing third-party data.
