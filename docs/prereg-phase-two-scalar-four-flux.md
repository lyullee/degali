# Pre-registration: four-flux phase two-scalar model

Date frozen: 2026-09-04, before enabling or running this combination.

Do not edit above the `RESULTS` line after a candidate result is known.

## Rationale

The first phase/two-scalar candidate failed because `scalar_peak` fixes
centreline fuel fraction and leaves only three Gaussian variables to close
species, momentum and energy. The existing `entrained_mass` establishment
instead solves velocity, width, density and centreline fraction against total
mass, fuel species, momentum and total energy. It therefore has the four
degrees of freedom required by the changed radial shapes.

This is a distinct boundary closure, not a retry with relaxed tolerances.

## Candidate fixed before running

1. Use four-flux `entrained_mass` establishment, independent Gaussian fuel
   mass-fraction and temperature departures, dry N2/O2 equilibrium and normal
   H2/N2/O2/H2O component enthalpy.
2. Keep `lambda_Y = 1.16` and fix
   `lambda_T = 1.16 * sqrt(0.70 / 0.85) = 1.052682847` from the published
   turbulent Schmidt and Prandtl values. Fit nothing.
3. First require the representative 1 mm boundary regression to close below
   `1e-8`. Only then run all nine Raman conditions at 81 radial points,
   0.25 mm maximum step and `5e-8` relative tolerance.
4. Score the corrected 369-point protocol first and retain the 549-point
   audit.

## Decision rule

1. all four 369-point slope errors must be within 25%;
2. both thermal absolute errors and the sum of all four absolute errors must
   improve over the recommended dry density-profile model;
3. neither mass error may worsen by more than two percentage points;
4. both radial coefficients must remain within the reported ranges;
5. maximum boundary residual below `1e-8`, species and energy drift below
   `2e-4`, positive state minima;
6. if those pass, a 121-point rerun must change each slope by less than 0.2%.

If the representative boundary fails, the candidate is rejected without a
nine-case run. If the boundary closes but an accuracy item fails, the option
is retained only for research and does not replace `scalar_peak`.

---

## RESULTS

Rejected as the default; retained as an explicit fully conservative research
closure.

The representative boundary passed, so the frozen nine-case run was
completed. On the corrected 369-point protocol the slopes are `0.21000150`,
`0.06636448`, `0.02002639` and `0.06047750`; relative errors are -20.03%,
+2.05%, -29.14% and -2.27%. The candidate corrects the radial-shape defect:
the temperature half-width error falls by 22.29 percentage points and mass
half-width absolute error also falls. But centreline-temperature error worsens
by 7.76 percentage points and crosses the 25% rejection boundary. It passes
only three of four slopes.

The four-error absolute sum falls from 75.11% to 53.48%, but the targeted rule
required both thermal errors to improve, not compensation of a failed centre
metric by a very good width metric. Median radial coefficients are 47.60 for
mass and 57.37 for temperature; the latter also lies above the reported
21--49 case range.

The result is numerically conservative: maximum boundary residual is
`1.42e-14`, species drift `8.88e-6`, and energy drift `2.88e-5`. State minima
remain positive, though the minimum radial temperature reaches 16.00 K. The
historical 549-point audit gives slopes `0.21616477`, `0.06578344`,
`0.02027905` and `0.05997515`; its centre-temperature error is likewise
-28.24%. A 121-point rerun was not made after the primary accuracy criterion
failed.

Conclusion: independently prescribed thermal spreading repairs the width but
does not supply the missing centreline warming. The next defensible model
change must alter a measured heat/source term or add an independently
constrained energy-flux closure; fitting `lambda_T` cannot solve both thermal
observables.

### Final-journal source correction

The revised primary errors are -24.21%, -6.12%, -29.14% and -2.27%. The
candidate still passes only 3/4 and remains rejected by the same thermal-centre
and radial-coefficient criteria. See
`hecht-panda-journal-benchmark-correction.md` for fit provenance.
