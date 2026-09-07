# Pre-registration: dew-point-heated two-scalar jet

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The dry-air-dew-point plug zone makes the conservative source physically
admissible and predicts both mass metrics, but retains the one-scalar thermal
width error. Independent mass-fraction/temperature profiles predict both
half-width slopes, but without the plug heating their Gaussian boundary enters
an unphysical cryogenic range. These closures act on successive, distinct
regions and should therefore be combined without changing either parameter.

## Candidate fixed before running

1. Apply HyRAM+ equations 58--64 from the ambient-pressure source to
   `T_min=82.1508333771 K`, the independently calculated 1 atm dry-air dew
   point. Conserve fuel mass, advective momentum and total energy.
2. Convert that heated plug with the four-flux `entrained_mass` boundary.
3. Use independent Gaussian `Y` and `T` profiles with `lambda_Y=1.16` and
   `lambda_T=1.16 sqrt(0.70/0.85)=1.0526828470`, based on the separately
   sourced turbulent Schmidt and Prandtl numbers.
4. Use the same constant-`cp*T` caloric relation throughout this isolated
   combination. Equilibrium condensation and fitted parameters remain off.

The physical orifice remains coordinate zero; the initial heating length and
Gaussian development length are both included before comparison at 40--100
mm. Existing defaults and production JETPLU are unchanged.

## Decision rule

1. Initial-zone and Gaussian-boundary invariant residuals below `1e-8` for all
   nine cases; downstream species/energy drift below `2e-4`.
2. All four printed Raman slopes within 25%.
3. Median radial coefficients inside 33--64 (mass) and 21--49 (temperature).
4. Every source and sampled profile is physical.
5. A 161-point rerun changes every slope by less than 0.2%.

Failure of any criterion rejects the candidate. No temperature, width,
entrainment, source state or virtual origin will be fitted.

---

## RESULTS

Rejected. At 121 radial points the four predicted slopes were
`0.21612722`, `0.06637858`, `0.01728037`, and `0.06023758`. Their relative
errors were -17.70%, +2.07%, -38.85%, and -2.65%, so only three of four met
the 25% rule. The mass radial coefficient, 47.86, was inside its 33--64
range, but the temperature coefficient, 58.11, exceeded its 21--49 range.

All conservation and physical-state checks passed: maximum initial-heating
residual `3.16e-16`, Gaussian-boundary residual `4.71e-15`, downstream
species drift `1.34e-10`, and total-energy drift `2.61e-10`. Temperatures
over the integrated solutions remained 66.84--214.09 K. A 161-point rerun
changed every slope by at most 0.014%. The failure is physical rather than
numerical: combining the two individually motivated closures fixes both
width slopes but further worsens centreline warming.
