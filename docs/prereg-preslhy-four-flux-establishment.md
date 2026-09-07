# Pre-registration: four-flux source establishment for PRESLHY coupling

Date frozen: 2026-09-05, after rejecting the scalar-peak coupled candidate on
vertical width and before calculating any four-flux coupled field prediction.

## Structural finding

The `scalar_peak` plug-to-Gaussian conversion fixes the published centreline
fuel fraction and solves only three equations: H2 mass, momentum and total
energy. It does not solve the total-mass equation. On the independent 1 mm
regression source, its Gaussian mass flux is 20.58% below the source plus the
air explicitly entrained over the formation length, even though the three
reported residuals are at machine precision.

That omission is unacceptable for a path advertised as conserved from the
atmospheric source through the crosswind handoff. It is also a physically
plausible cause of the PRESLHY scalar-peak candidate's vertical-width ratio
of 0.710. No measured width is used to select a coefficient or new profile.

## Fixed candidate

- Change only Gaussian establishment from `scalar_peak` to the existing
  `entrained_mass` four-flux solve.
- Solve centre velocity, width, density and H2 fraction against total mass,
  H2 mass, vector-momentum magnitude and total energy at the same published
  formation distance.
- Retain the phase-equilibrium density profile, normal-H2 component enthalpy,
  source-momentum entrainment beta 0.28, spreading ratio 1.16, humidity,
  `10D` handoff, phase-manifold mapping, adaptive energy quadrature and every
  numerical setting in
  `prereg-preslhy-coupled-crosswind-validation.md`.
- Retain the same condition-only trial set 10, 11, 12, 22, 23, 24 and 25 and
  the same common-arc/sensor comparison.

The scalar-peak option remains available only to reproduce the accepted Raman
temperature study. It is not relabelled as four-flux conservative.

## Additional gate

At all seven sources, each plug-to-Gaussian residual for total mass, H2 mass,
momentum and total energy must be at most `1e-8`. All previously frozen
near-field and JETPLU handoff gates also remain.

## Decision

Compare against the same existing corrected JETPLU baseline. Promotion of
the coupled field path still requires all six criteria from the original
PRESLHY pre-registration: 7/7 interfaces, improved `abs(log(MG))`, lower VG,
non-lower FAC2, width ratio closer to one and non-larger centre-height MAE.

If the four-flux source boundary fails or the field criteria fail, do not tune
the formation distance, spreading ratio or entrainment coefficient. Retain
the model as a physically conservative research closure and identify the next
missing mechanism explicitly.

## Results

Rejected. All seven source boundaries closed total mass, H2 mass, momentum
and energy; the largest absolute residual was `1.31e-15`. Five handoffs
passed, but trial 10 missed the 2 K centre-temperature screen at 2.069 K and
trial 25 missed it at 2.106 K as well as the 5% H2-width screen at 5.036%.
The seven-of-seven interface condition therefore failed.

The surviving five-trial diagnostic again points away from source-boundary
mass as the width cause. On 33 common concentration arcs, baseline versus
candidate was MG 1.146 versus 0.970, VG 1.240 versus 1.197 and FAC2 0.970 for
both. On 12 vertical fits, centre MAE was 0.040 m for both, while the width
ratio was 1.149 for the baseline and **0.713** for the candidate. Closing the
previously omitted total-mass equation did not repair downstream vertical
spread.

The next isolated mechanism is therefore the entrainment representation after
the handoff. The current coupled path keeps HyRAM's constant source-momentum
entrainment and disables JETPLU's local density-scaled shear term; the latter
is part of the independently validated corrected field model.
