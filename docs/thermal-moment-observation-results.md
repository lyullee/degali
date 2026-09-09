# Thermal-moment observation and constitutive-domain results

Date: 2026-09-10 (Asia/Seoul)

## Scope

This stage subjects the coefficient-free, direct-flux thermal-moment path to
the same finite sensor operator as the PRESLHY observations.  Trials 10 and 23
are evaluated at x = 1.78 and 4.00 m on the five registered vertical sensors.
The control, candidate population, statistics and acceptance rule were frozen
before the new predictions were calculated.  No observation was used to fit a
coefficient.

The candidate transports total mass, H2 mass, horizontal and vertical
momentum, total energy and the advective-enthalpy second moment.  It retains
`D_h/D_C = 1`, reduced buoyancy work and the default linear phase lookup.

## Corrected ground-image observation operator

For an independent enthalpy width, the direct and ground-reflected images must
be superposed separately for species and enthalpy:

```text
C is proportional to exp(-q_y) [exp(-q_d) + exp(-q_i)]
H is proportional to exp(-p q_y) [exp(-p q_d) + exp(-p q_i)]
p = 1 / beta_H^2
```

Raising the already summed species image to `p` would introduce a cross term
that is not present in either physical image.  The implementation now follows
the separate-image expression and exactly recovers the former receptor when
`beta_H = 1`.

## Why the original full-square gate stopped

With the original `full_square` positivity check, Trial 10 stopped at
x = 3.346286 m.  The weak balance was still `6.13e-11`; the failure was a
single negative sampled momentum diffusivity at the outer square corner.

The radial constitutive identity is derived by integrating around a complete
constant-q contour.  Such a contour is complete only inside the square's
inscribed circle, `q <= q0`.  In the clipped corner, `q0 < q <= qmax`, side
flux terms also enter and the complete-contour identity is not a local eddy
diffusivity formula.  At the rejected half stage, the diffusivity remained
positive through `q/qmax = 0.9999`; the sign loss appeared only between
0.9999 and 0.99999 as the surviving angular measure fell below
`4.0e-4 rad`.  This locates a coordinate-domain inconsistency rather than a
bulk negative-transport region.

An opt-in `radial_core` validity domain was therefore preregistered.  It checks
local species and momentum diffusivities only on complete radial contours
`0 <= q <= q0`.  It does **not** clip or floor a diffusivity.  Full-square
mass, momentum, energy and thermal-moment integrals, reservoir side fluxes,
edge defects and source terms are unchanged.  `full_square` remains the
package default.

## Numerical result with the radial-core domain

Both nominal arc steps reached 4.00 m for both trials and passed every frozen
numerical gate.

| Trial | max independent balance, coarse/refined | max terminal parameter difference | max terminal flux difference | initial/final refined thermal-width ratio |
|---:|---:|---:|---:|---:|
| 10 | 6.788e-11 / 6.788e-11 | 2.886e-5 | 6.079e-5 | 1.04810 / 0.99688 |
| 23 | 5.797e-11 / 5.797e-11 | 5.762e-6 | 6.914e-5 | 1.03304 / 0.99740 |

The maximum scaled coarse/refined difference across all 40 sensor predictions
was `3.534e-4`, below the frozen `0.005` limit.  This demonstrates numerical
reach, conservation and step convergence; it is not an accuracy score.

## Observation result

The registered primary metrics are median absolute log ratios over the four
matched profiles.  Lower is better.

| Primary metric | control | candidate, coarse | candidate, refined | decision |
|---|---:|---:|---:|---|
| Thermal centre amplitude | 0.179243 | 0.190701 | 0.190702 | worse |
| Thermal variance | 0.268331 | 0.251458 | 0.251463 | better |

The candidate improves thermal variance by about 6.3% but worsens centre
amplitude by about 6.4%.  Because the preregistered decision requires both
components to be no worse, the observation gate fails at both step sizes.
The candidate may not advance and no default is changed.

The refined model/observation ratios show where the remaining mismatch lies:

| Trial | x (m) | thermal centre | thermal variance | thermal M0 | H2 centre | H2 variance | H2 M0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1.78 | 1.012 | 1.494 | 1.269 | 1.736 | 1.640 | 2.336 |
| 10 | 4.00 | 1.862 | 1.248 | 2.397 | 5.407 | 1.708 | 7.673 |
| 23 | 1.78 | 0.924 | 1.325 | 1.067 | 0.968 | 1.374 | 1.189 |
| 23 | 4.00 | 1.353 | 1.051 | 1.386 | 1.430 | 1.312 | 1.688 |

These are correlated profiles from two experiments, not eight independent
validation experiments.  M0 is a five-sensor finite-line diagnostic rather
than a complete plume cross-section integral.

## Physical interpretation and next boundary

The transported thermal-width ratio rapidly relaxes from 1.048/1.033 to about
one.  Under the present instantaneous closure, shear production is deposited
directly into mean enthalpy and `D_h/D_C` is fixed at one.  The new result says
that independently carrying the thermal second moment is numerically viable,
but this closure cannot simultaneously correct amplitude and spread.

DEGALI already contains a research-only finite-TKE operator that separates
mean-shear production into `Q = rho k` and transfers only positive dissipation
from Q to enthalpy, with a verified total-energy identity.  It must not be
activated with arbitrary coefficients: the available PRESLHY workbooks have
no velocity RMS, Reynolds stress, k or epsilon with which to identify its
initial state and dissipation closure.  The next scientifically defensible
step is an independently measured cryogenic-jet mean-velocity and turbulence
profile, or a clearly labelled end-member sensitivity that cannot be promoted
as validation.  Retuning a thermal diffusivity from the four temperature
profiles is not justified.

The raw experimental files and immutable audit JSON remain local and excluded
from distribution.  Public preregistration and provenance documents identify
the source and reproduce the acceptance logic without redistributing third-
party data.
