# Pre-registration: temperature-dependent enthalpy with air phase equilibrium

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The dry equilibrium N2/O2 model uses constant heat capacities evaluated near
ambient temperature over a 45--295 K domain. For hydrogen this overstates the
50--295 K ideal-gas enthalpy increase by about 11%, requiring too much
entrainment to warm the jet and potentially causing the observed
centreline-temperature decay slope to be too small. The earlier gas-only
two-scalar test was thermodynamically inadmissible below the air property
range; simultaneous condensed-air equilibrium removes that specific defect.

## Candidate fixed before running

1. Use low-density Helmholtz ideal-gas enthalpy from CoolProp for H2, N2, O2
   and H2O on 14.1--400 K tables. Shift each component reference to zero at
   the reported 295 K ambient temperature.
2. Within the existing N2/O2 equilibrium model, calculate gas sensible
   enthalpy from those component tables and condensed enthalpy as saturated
   vapor enthalpy minus the existing phase-change latent heat.
3. Keep ambient humidity zero to isolate caloric properties. Keep the frozen
   `scalar_peak` conservative establishment, density/species Gaussian
   profiles, entrainment and spreading coefficients unchanged.
4. Use 81 radial points and the previously converged phase-integration
   settings (0.25 mm maximum step, relative tolerance `5e-8`). Fit no heat
   capacity or temperature parameter.

## Decision rule

1. Boundary residuals below `1e-8`; downstream species and total-energy drift
   below `2e-4` in every case.
2. All four printed Raman slopes within 25%, and both median radial
   coefficients within their reported case ranges.
3. All component states remain inside the tabulated temperature range and all
   physical state variables remain positive.
4. A 121-point rerun changes every slope by less than 0.2%.

Failure of any item rejects the candidate. The existing constant-property
default remains available and unchanged.

---

## RESULTS

Accepted for the conserved axisymmetric Raman research path. At 81 radial
points and the strict integration settings, the original 549-point audit
predicts slopes `0.21505534`, `0.05964072`, `0.02245221`, and `0.07579373`.
Their relative errors are -18.11%, -8.29%, -20.55%, and +22.49%; all four
meet the 25% rule. Median radial coefficients 60.21 and 36.22 are inside the
reported mass and temperature ranges.

Maximum nine-case boundary residual is `1.50e-14`, species drift `1.53e-5`,
and total-energy drift `3.77e-5`. Integrated temperatures remain
51.19--231.52 K and all state variables are positive. A 121-point rerun
changes every slope by less than 0.05%.

Under the adopted unequal 369-point experimental coverage, the slopes are
`0.20820525`, `0.05952579`, `0.02221798`, and `0.07707466`, with errors
-20.71%, -8.46%, -21.38%, and +24.56%; the candidate therefore remains 4/4
under the primary protocol. This is the first tested establishment that both
conserves the source-to-Gaussian invariants and passes all four printed Raman
slopes. It supersedes constant-`cp` as the recommended dry-air phase research
model, but does not become a production atmospheric default because the
validation experiment's humidity remains unknown.

### Final-journal source correction

Using the final 2019 mass-fit labels, the 549-point errors are -22.39%,
-15.63%, -20.55% and +22.49%. The primary 369-point errors are -24.86%,
-15.79%, -21.38% and +24.56%. The candidate stays 4/4 and remains the
recommended conservative dry research closure, but two primary metrics are
now very close to the frozen boundary. Because the aggregate plots contain an
untabulated tenth series and no fit uncertainty, the pass is provisional; see
`hecht-panda-journal-benchmark-correction.md`.
