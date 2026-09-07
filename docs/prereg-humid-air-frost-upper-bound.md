# Pre-registration: humid-air frost upper bound

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The Hecht--Panda Raman experiment explicitly reports condensed ambient
moisture in the cold jet, while the current equilibrium phase bound contains
only dry-air nitrogen and oxygen. Freezing entrained water releases a large
latent heat per unit mass and can therefore reduce the conservative model's
underprediction of centreline warming.

## Candidate fixed before running

1. Use 100% relative humidity at the reported 295 K and 1 atm ambient state.
   Calculate absolute humidity from the water saturation pressure and ideal
   dry-air/water molecular weights; fit no humidity value to Raman results.
2. Carry that water mass with every unit of entrained ambient gas. At each
   radial quadrature point solve simultaneous H2/N2/O2/H2O ideal-gas partial
   pressures and equilibrium condensed fractions at the local bulk density.
3. Below the water triple point use the already implemented DEGALI
   ice-sublimation pressure continuation. Include vapor sensible enthalpy,
   vaporisation plus fusion latent heat, and condensed-water volume with the
   existing core constants.
4. Retain the frozen `scalar_peak` Gaussian establishment and all dry-air
   N2/O2 equilibrium, entrainment, spreading and source parameters. Treat
   condensate as locally equilibrated and no-slip. This deliberately
   maximises thermal feedback and is an upper bound, not a calibrated model.

## Decision rule

1. Boundary invariant residuals below `1e-8`; downstream species and total
   energy drift below `2e-4` in every case.
2. All four printed Raman slopes within 25%, with median radial coefficients
   inside 33--64 (mass) and 21--49 (temperature).
3. Every sampled phase fraction, density and temperature is physical.
4. A 121-point radial rerun changes every slope by less than 0.2%.

If the 100%-RH upper bound still fails centreline temperature, unknown
laboratory humidity cannot rescue this closure and no humidity measurement
will be requested. If it passes, the experimental RH becomes required before
claiming predictive validation at realistic humidity.

---

## RESULTS

The phase implementation is accepted, but the 100%-RH validation candidate is
rejected. The independently calculated saturation humidity is
`0.01651706 kg H2O/kg dry air`. At 81 radial points the four slopes are
`0.20101190`, `0.05797354`, `0.04092011`, and `0.07148364`, with relative
errors -23.45%, -10.85%, +44.80%, and +15.52%. Three of four pass; both
median radial coefficients (63.78 and 40.44) are inside their source ranges.

The nine-case maximum boundary residual is `1.08e-15`, species drift
`4.80e-6`, and total-energy drift `9.71e-5`. Integrated temperatures remain
45.73--255.99 K, densities and mass fractions remain positive, and a
121-point rerun changes every slope by less than 0.02%. The candidate fails
only because maximum frost feedback overpredicts centreline warming.

This is nevertheless a material model correction. The dry equilibrium bound
underpredicts the temperature-decay slope by about 31%, whereas the saturated
bound overpredicts it by about 45%. An explicitly labelled 40%-RH sensitivity
case gives slopes `0.20904988`, `0.05906529`, `0.02567265`, and `0.07390881`,
all within 25%, without changing any jet coefficient. This is not accepted as
a calibrated validation because the paper reports condensed moisture but not
laboratory relative humidity. The actual experimental RH is required to score
the humid model predictively.

### Final-journal source correction

The revised 100%-RH errors are -27.46%, -17.99%, +44.80% and +15.52%, so the
upper bound now passes only 2/4. The unscored 40%-RH sensitivity remains 4/4
at -24.56%, -16.45%, -9.16% and +19.44%. Its interpretation is unchanged:
unknown measured humidity and aggregate-fit provenance forbid calibration.
See `hecht-panda-journal-benchmark-correction.md`.
