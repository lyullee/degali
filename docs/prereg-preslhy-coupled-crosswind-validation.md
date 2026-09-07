# Pre-registration: PRESLHY coupled near-field/crosswind validation

Date frozen: 2026-09-05, before calculating any coupled-model sensor
prediction.

## Question

Does the conserved axisymmetric near field followed by the conservative
JETPLU handoff reduce the existing corrected model's error on independent
PRESLHY E3.5 horizontal releases, without degrading plume width or centre
height?

This is a validation of the new coupling, not a refit.  No coefficient or
screen below may be changed after looking at the coupled predictions.

## Fixed data and population

- Input: `reference/preslhy/e35_reduced.json`.
- Source flow: `flow_mean_gs`, matching the current reported validation.
- Only unobstructed horizontal releases are eligible.
- The existing orifice velocity/wind screen remains 10.
- The coupled model's stricter applicability screen is also enforced at the
  already air-loaded atmospheric source plane: source velocity divided by
  local logarithmic-profile wind must be at least 10.
- Those condition-only rules select trials **10, 11, 12, 22, 23, 24 and 25**.
  Trials 20 and 21 pass the old orifice screen but fail the atmospheric-source
  screen; they are excluded before inspecting coupled predictions.
- An arc or vertical fit is used only when its downwind coordinate is at or
  beyond the handoff coordinate and both models reach it.  The existing model
  is then restricted to exactly the same trial/arc set.
- Observed concentration is the maximum reported peak across the sensors at
  one downwind arc.  Prediction is evaluated at every corresponding sensor
  and maximised over the same sensors.  Concentrations are never compared to
  an unobserved model centreline.
- Vertical geometry uses only the reduced file's `well_constrained` Gaussian
  fits and its declared standard-deviation convention.

## Fixed source and solver

The candidate begins from the exact corrected atmospheric source already
constructed by `validation.nearfield.hydrogen_jet`: expanded diameter,
momentum-consistent total velocity, mixed density, flash temperature and H2
mass fraction.  The H2 mass flow remains the measured window mean; the total
flow is H2 flow divided by source mass fraction.  No second flash or empirical
source adjustment is allowed.

For every selected trial:

- release angle: horizontal and aligned with the mean wind;
- atmosphere, RH, reference wind height, roughness 0.001 m, stability D and
  60 s averaging: the same inputs as the existing corrected calculation;
- near-field handoff arclength: **10 atmospheric-source diameters**;
- near-field radial points: **41**;
- maximum near-field step: `min(0.001 m, source diameter / 20)`;
- relative tolerance: **2e-6**;
- minimum H2 mass fraction: **7e-4**;
- thermodynamic projection: **phase_manifold**;
- JETPLU output step: **0.02 m**, tolerance **1e-5**, reach **40 m**.

The phase-manifold mapping is the previously registered domain extension; it
is not chosen from these field results.  The crosswind coefficients remain
those already fixed in `run_lh2_crosswind_research` (`Sc=1.16^2`, continuous
momentum entrainment beta 0.28 and Houf--Schefer buoyancy entrainment).  There
is no trial-specific coefficient.

## Interface gates

Every one of the seven trials must pass the existing frozen handoff gates:

- native total-mass, H2-mass and vector-momentum relative residuals at most
  `1e-8`;
- total-energy flux mismatch at most 2%;
- H2 half-width mismatch at most 5%;
- centre-temperature mismatch at most 2 K;
- energy quadrature convergence residual at most `1e-5`;
- conservative near-field audit and source velocity/local-wind ratio at least
  10.

Failure of any trial rejects the coupled path as a general replacement on
this population.  Results from a surviving subset may be reported only as a
diagnostic.

## Metrics and decision

On the exact common concentration arcs, calculate MG, VG, FAC2 and a
trial-bootstrap 95% interval for MG.  On the exact common vertical fits,
calculate mean model/measured `sigma_z`, mean signed centre-height error and
mean absolute centre-height error.

Promotion over the current corrected model requires all of the following:

1. all seven interface gates pass;
2. `abs(log(MG))` is smaller;
3. VG is smaller;
4. FAC2 is not smaller;
5. the mean `sigma_z` ratio is closer to one;
6. centre-height mean absolute error is not larger.

If any condition fails, the current measured-source JETPLU path remains the
validated default and the coupled route remains research-only.  Failed
criteria will not be repaired by changing handoff distance, entrainment,
Schmidt number, source flow choice or trial membership.

## Audit correction recorded after the first run

The first version of this file accidentally transcribed the already existing
handoff limits as 1% energy, 5 K temperature and `1e-4` quadrature.  The code
and the earlier phase-profile pre-registration, both predating this test,
actually enforce 2%, 2 K and `1e-5`; the bullets above have been corrected to
match that immutable implementation.  This correction cannot improve the
first field result: every trial was already inside both the mistyped and real
energy/temperature limits.  Five trials failed only because the fixed 32/64
quadrature pair had not converged to the real `1e-5` threshold.

The numerical remediation is recorded separately in
`preslhy-handoff-quadrature-remediation.md`.  It retains the `1e-5` threshold
and only increases integration order; no physical parameter, population,
handoff station or promotion criterion changes.

## Results

After adaptive quadrature, all seven scalar-peak interfaces passed. Orders 64
or 128 were sufficient. Maximum native flux residual was `1.10e-15`; energy,
H2 half-width and centre-temperature mismatches were at most 0.451%, 4.484%
and 1.245 K.

On 42 common arcs, the corrected baseline and coupled candidate gave:

| model | n | MG | trial-bootstrap 95% interval | VG | FAC2 |
|---|---:|---:|---:|---:|---:|
| corrected JETPLU baseline | 42 | 1.074 | 0.863--1.278 | 1.213 | 0.952 |
| scalar-peak coupled | 42 | **0.925** | 0.765--1.081 | **1.183** | 0.952 |

On 17 common vertical fits, centre-height MAE improved from 0.054 to 0.049 m,
but mean model/measured `sigma_z` moved from 1.091 to **0.710**. The candidate
therefore fails the frozen width criterion and is not promoted despite its
better concentration and centre-height statistics.

Code audit then found that `scalar_peak` solves H2 mass, momentum and energy
but not total mass at Gaussian establishment. The four-flux follow-up is
recorded in `prereg-preslhy-four-flux-establishment.md`.
