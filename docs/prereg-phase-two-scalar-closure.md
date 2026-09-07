# Pre-registration: phase-equilibrium two-scalar radial closure

Date frozen: 2026-09-04, before the first combined model run.

Do not edit above the `RESULTS` line after a candidate result is known.

## Physical question

The recommended dry phase-equilibrium model prescribes Gaussian mixture
density and hydrogen mass-density profiles with the same spreading ratio.
Temperature is then inferred from those two profiles. This produces a
temperature half-width wider than the mass-fraction half-width, whereas the
Hecht--Panda fits give `0.06188` and `0.06503 mm/mm`, respectively.

The existing two-scalar experiment prescribed separate fuel-fraction and
temperature profiles using the published gradient-diffusion values
`Pr_t = 0.85` and `Sc_t = 0.70`, but it was tested before the physically
admissible N2/O2 phase equilibrium and component `h(T)` model existed. Its
gas-only caloric branch failed or under-warmed the centre. The missing
combination is therefore tested without fitting either diffusivity.

## Candidate fixed before running

1. Retain `scalar_peak` conservative establishment, dry N2/O2 equilibrium,
   normal-hydrogen component enthalpy, the accepted entrainment closure and
   the species spreading ratio `lambda_Y = 1.16`.
2. Prescribe independent Gaussian fuel mass-fraction and temperature
   departures. Fix the thermal ratio from the external transport constants:
   `lambda_T = 1.16 * sqrt(0.70 / 0.85) = 1.052682847`.
3. At each radius, use that ideal temperature/composition pair only to define
   total density, then solve the existing local N2/O2 phase equilibrium and
   component enthalpy. No phase, latent heat or density correction is bypassed.
4. Use 81 radial points, 0.25 mm maximum step and relative tolerance `5e-8`.
   Run all nine conditions and score both 549- and corrected 369-point
   protocols. Fit no coefficient.

## Decision rule

Relative to the recommended dry density-profile baseline:

1. all four corrected 369-point slope errors must remain within 25%;
2. the sum of their absolute errors must decrease;
3. each thermal absolute error must decrease, because radial thermal shape is
   the only targeted defect;
4. neither mass error may worsen by more than two percentage points;
5. both median radial coefficients must remain in the reported case ranges;
6. maximum boundary residual must be below `1e-8`, downstream species and
   energy drift below `2e-4`, all state minima positive;
7. if the primary criteria pass, a 121-point rerun must move every slope by
   less than 0.2%.

Failure of any accuracy item rejects the combination while retaining it as an
explicit research closure. The existing dry density-profile model remains the
default during the test.

---

## RESULTS

Rejected before the nine-case downstream regression.

The combination was first exercised on the representative 1 mm regression
source used by the axisymmetric unit tests. With the frozen
`scalar_peak` boundary, the three adjustable Gaussian variables could not
simultaneously close species, momentum and energy after the independent
temperature profile was introduced. The final relative residuals were
`+3.927%` species, `-0.295%` momentum and `+3.939%` energy, far above the
`1e-8` boundary criterion.

This is a structural incompatibility at the plug-to-Gaussian handoff, not a
downstream numerical-integration failure. Running all nine cases or tuning the
thermal ratio after that failure would violate the frozen decision order. The
previous constructor guard requiring the density profile with equilibrium air
condensation was therefore restored, and a regression test now keeps this
invalid combination from being selected accidentally.

A future enthalpy-flux formulation could transport an independent thermal
shape by adding a separately evolved energy-profile state and solving a new
four-flux establishment problem. That is a new model, not a safe parameter
switch in the present five-balance closure.

### Final-journal source correction

The final observed temperature-to-mass half-width ratio is
`0.06188 / 0.07069 = 0.8754`, rather than 0.9516. The later conservative
four-flux test already uses an independently published Prandtl/Schmidt ratio
and remains the applicable structural audit; this failed three-variable
boundary is not retried. See `hecht-panda-journal-benchmark-correction.md`.
