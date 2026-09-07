# Pre-registration: Hecht--Panda 0.3 m/s co-flow

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The Raman experiment surrounded each vertical hydrogen jet with a measured
0.3 m/s air co-flow through a 19 cm honeycomb. The present axisymmetric model
assumes stationary ambient air. Carrying the co-flow's mass, axial momentum
and kinetic energy may slightly reduce relative shear and alter the decay
rates without changing an empirical jet coefficient.

## Candidate fixed before running

1. Set the ambient co-flow to 0.3 m/s, aligned with the source axis.
2. Evaluate source and local entrainment using jet/co-flow relative velocity.
3. Add entrained-air momentum `rho_a E U_a` in the fixed co-flow direction and
   kinetic-energy influx `rho_a E U_a^2/2` to the five flux balances.
4. Include the same development-zone contributions in the conservative
   plug-to-Gaussian target fluxes.
5. Use the dry, constant-`cp`, four-flux `entrained_mass` establishment at 121
   radial points. Do not combine humidity, phase equilibrium, initial heating
   or differential scalar spreading in this isolated test.

## Decision rule

1. Boundary residuals below `1e-8`; downstream species and co-flow-corrected
   total-energy drift below `2e-4` in every case.
2. All four printed Raman slopes within 25%, and both median radial
   coefficients inside the reported case ranges.
3. A 161-point rerun changes every slope by less than 0.2%.
4. All states remain physical.

Failure of any item rejects co-flow as a sufficient model correction. The
measured 0.3 m/s is fixed and is not fitted.

---

## RESULTS

Rejected as a sufficient correction. At 121 radial points the four slopes are
`0.21646845`, `0.06091483`, `0.01945153`, and `0.07586941`, with errors
-17.57%, -6.33%, -31.17%, and +22.61%. Three of four pass, and the median
radial coefficients 57.35 and 36.02 are inside their experimental ranges, but
the centreline-temperature criterion still fails.

The implementation is numerically conservative: maximum nine-case boundary
residual `1.52e-15`, species drift `3.27e-8`, and co-flow-corrected
total-energy drift `1.03e-7`. All sampled state variables are physical. A
161-point rerun changes every slope by less than 0.01%. Relative to the
stationary-air conservative baseline, every slope changes by less than 0.8%.
The measured co-flow is retained as an explicit boundary condition, but it
cannot explain the missing thermal feedback and is not enabled by default.
