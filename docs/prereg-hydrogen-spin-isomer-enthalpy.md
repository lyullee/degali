# Pre-registration: hydrogen spin-isomer caloric sensitivity

Date frozen: 2026-09-04, before the first para-hydrogen model run.

Do not edit above the `RESULTS` line after a candidate result is known.

## Physical question

The accepted component-enthalpy model currently evaluates the hydrogen
caloric term with CoolProp `Hydrogen`, i.e. normal hydrogen.  That is a
potentially important assumption for a stream derived from stored liquid
hydrogen.  NIST reports that ortho/para composition changes gas enthalpy and
heat capacity much more strongly than pressure-volume-temperature behaviour,
and provides separate reference equations for para-, normal- and
ortho-hydrogen.  CoolProp exposes the corresponding `ParaHydrogen` equation
of state from Leachman et al.

No in-flight equilibrium conversion is added.  NBS Monograph 168 reports that
uncatalysed self-conversion is very slow (a half-life longer than one year at
liquid-air temperatures), whereas the Raman jet residence time is milliseconds.
The candidate therefore preserves a fixed spin composition along the jet.

Sources:

- Leachman et al., *Fundamental Equations of State for Parahydrogen, Normal
  Hydrogen, and Orthohydrogen*, DOI 10.1063/1.3160306.
- McCarty, Hord and Roder, *Selected Properties of Hydrogen*, NBS Monograph
  168, section 4, ortho-para conversion.
- CoolProp `ParaHydrogen` fluid documentation.

## Candidate fixed before running

1. Replace only the fuel ideal-gas component enthalpy table from
   `Hydrogen` to `ParaHydrogen`; keep the public species identity, molecular
   weight, source expansion, density relation and N2/O2 equilibrium unchanged.
2. Retain dry air, `scalar_peak` conservative establishment, component
   temperature-dependent enthalpy, 81 radial points, a 0.25 mm maximum step
   and relative tolerance `5e-8`.
3. Run all nine Hecht--Panda conditions once and score both the historical
   549-point audit and the corrected unequal 369-point camera coverage.
4. Fit no spin fraction, heat capacity, entrainment coefficient, spreading
   coefficient or temperature offset.

## Decision rule

Relative to the accepted normal-hydrogen dry baseline on the primary
369-point protocol:

1. all four printed slopes must remain within the frozen 25% band;
2. the sum of the four absolute relative errors must decrease;
3. neither mass metric may worsen by more than one percentage point;
4. the centreline-temperature error must decrease without the temperature
   half-width crossing the 25% limit;
5. maximum boundary residual must remain below `1e-8`, and downstream species
   and externally corrected total-energy drift below `2e-4`.

The accepted normal-hydrogen baseline slopes are `0.20820525`, `0.05952579`,
`0.02221798` and `0.07707466`; its relative errors are -20.71%, -8.46%,
-21.38% and +24.56%.

Even if all numerical criteria pass, the result is labelled a
`ParaHydrogen` sensitivity rather than a default change until the released
hydrogen's spin composition is documented.  A failed candidate remains an
explicit, auditable off-default option.

---

## RESULTS

Rejected as the default; retained as an explicit caloric sensitivity.

The corrected 369-point slopes are `0.21150182`, `0.06019488`,
`0.01842779` and `0.07709132`.  Their relative errors are -19.46%, -7.44%,
-34.79% and +24.58%, so only three of four slopes meet the 25% rule.  The
mass-centre and mass-width errors improve by 1.26 and 1.03 percentage points,
but the centreline-temperature error worsens by 13.41 percentage points.  The
sum of absolute errors rises from 75.11% to 86.27%.

The historical 549-point slopes are `0.21802688`, `0.06027028`,
`0.01871118` and `0.07584022`, with errors -16.97%, -7.32%, -33.79% and
+22.56%.  The failed thermal-centre criterion is therefore not caused by the
camera-coverage correction.

Numerically, the candidate is sound: the maximum boundary residual is
`4.64e-15`, species drift `1.61e-5`, and energy drift `3.39e-5`; integrated
temperatures remain 49.70--222.98 K and all reported state minima are
positive.  A 121-point rerun was not performed because the primary accuracy
criterion already failed decisively.

This result also fixes the sign of the diagnosis.  Preserving para-hydrogen
calorics makes the computed plume warm more slowly, whereas the Raman centre
already warms too slowly.  Missing spin-isomer composition remains a source
uncertainty, but it cannot be invoked as the unmodelled heat source needed to
close the present residual.  The recommended dry research configuration
therefore continues to use normal-hydrogen calorics.

### Final-journal source correction

The revised primary errors are -23.67%, -14.85%, -34.79% and +24.58%.
Para-hydrogen still passes only 3/4, so the rejection is unchanged. See
`hecht-panda-journal-benchmark-correction.md` for the provisional aggregate-fit
qualification.
