# Pre-registration: compatible handoff search for the four-flux field path

Date frozen: 2026-09-05, before calculating the combined candidate.

## Rationale

Four-flux establishment closes source total mass, H2 mass, momentum and
energy, but at a fixed 10D handoff trials 10 and 25 narrowly miss the existing
profile-compatibility screen. Trial 10 has a 2.069 K temperature mismatch;
trial 25 has 2.106 K temperature and 5.036% H2-width mismatches. The limits
remain 2 K and 5% and will not be relaxed.

The separately isolated downstream local-shear transition repairs the field
width deficit. A fully conservative combined path therefore needs a handoff
station chosen by interface compatibility, not by a concentration residual.

## Fixed search

- Use four-flux `entrained_mass` establishment and the accepted post-handoff
  local density-scaled shear closure.
- Integrate the conserved axisymmetric near field through 20 atmospheric-
  source diameters.
- Test handoffs at 10D, 11D, ..., 20D, in that order.
- Select the first station that passes every already frozen native-flux,
  energy, H2-width, centre-temperature, quadrature and state-domain screen.
- If no station passes by 20D, reject that trial and the combined candidate.
- A failed trial/station may not be rescued by interpolation, a fractional-D
  search, threshold relaxation or inspection of measured concentration.
- Retain the same seven trials, reduced data, sensor/arc rules, coefficients,
  atmosphere, solver tolerances and quadrature remediation as the earlier
  PRESLHY protocol.

## Decision

All seven trials must find an accepted station. Against the corrected JETPLU
baseline on exactly common downstream arcs, the combined candidate must:

1. improve `abs(log(MG))`;
2. lower VG;
3. not lower FAC2;
4. move mean model/measured `sigma_z` closer to one; and
5. not increase centre-height MAE.

All four source-establishment flux residuals must remain below `1e-8`. No
field observation participates in the handoff search.

## Results

Rejected. Trials 11, 12, 22, 23 and 24 passed at the first 10D station.
Trials 10 and 25 found no accepted station through 20D. At 20D their
centre-temperature mismatches had grown to 6.415 and 6.553 K and their H2
half-width mismatches to 8.179% and 8.347%, respectively. Moving the boundary
downstream makes the incompatibility worse rather than resolving a short
formation-zone transient.

The five accepted trials still confirm the local-shear diagnosis. Across 12
common vertical fits their candidate width ratio is 1.019, versus 1.149 for
the corrected baseline, and centre-height MAE is unchanged at 0.040 m. Across
33 concentration arcs, baseline versus candidate MG/VG/FAC2 is 1.146/1.240/0.970
versus 1.165/1.231/0.970.

The combined candidate is not promoted because the seven-of-seven interface
gate fails. The next boundary model must carry an independent thermal/energy
profile state; adjusting handoff position or a width coefficient is ruled
out by this result.
