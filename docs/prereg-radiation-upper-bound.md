# Pre-registration: blackbody radiation upper bound

**Frozen before the first radiation-bound result.** Do not edit above the
`RESULTS` line after a candidate result is known.

## Question

Could omitted room-temperature thermal radiation explain a meaningful part
of the remaining Hecht--Panda near-field warming error?

## Bound fixed before running

No gas or frost emissivity is fitted. At every integrated station, assume the
entire cylinder of radius `5B` is a perfectly absorbing black surface held at
the centreline temperature. Add

`dQ/dS = 2 pi (5B) sigma (T_ambient^4 - T_center^4)`

along the trajectory, using the exact 2022 CODATA Stefan--Boltzmann constant
`5.670374419e-8 W m-2 K-4`. This deliberately overstates real absorption:
the outer Gaussian plume is warmer than the centreline, clear H2/air is not a
black surface, and the same radiative field cannot be fully absorbed at every
radius.

Run the unchanged accepted 81-point dry-air model for all nine sources. For
each case compare cumulative blackbody heat through 0.13 m with (a) the
absolute initial cross-sectional enthalpy-deficit flux and (b) the source H2
mass flow multiplied by the absolute centreline mixture-enthalpy change
between the first and last state. The latter is only a common-scale warming
denominator, not an additional conservation equation.

## Decision rule

If the largest blackbody energy ratio is below 1%, radiation is screened out
of the near-field ODE and remains a documented negligible upper bound. If it
is 1% or larger, implement an explicit emissivity/optical-depth sensitivity,
without fitting to Raman residuals, and repeat the conservation and four-slope
tests before any adoption decision.

### Follow-up fixed after the bound and before any radiative ODE result

The blackbody bound exceeded 1%, so the triggered follow-up is:

1. add a dimensionless effective absorptivity `a_eff` to the energy source,
   with `a_eff = 1 - exp(-tau)` as the equivalent grey optical-depth mapping;
2. integrate absorbed radiation in a separate cumulative energy ledger and
   require `(energy flux - absorbed radiation)` to conserve within `2e-4`;
3. run the full-black limit `a_eff=1` through all nine cases at the accepted
   81-point/0.25-mm/`5e-8` settings and corrected 369-point coverage;
4. infer the monotone lower-absorptivity sensitivity from the explicit
   blackbody ratios and report the optical depth at which each energy ratio
   reaches 1%. No optical depth is selected from the Raman residuals.

Full-black radiation is adopted only if it preserves all existing acceptance
checks and improves the absolute error of every affected thermal metric
without degrading either mass metric. Otherwise the option remains a bounded
sensitivity and the recommended configuration remains radiatively adiabatic.

---

## RESULTS

The deliberately excessive perfect-black cylinder adds at most `11.8899 W`
over 0.13 m. Its largest ratio is `1.3820%` of the initial cross-sectional
enthalpy deficit and `1.6905%` of the centreline warming scale. The initial
1% screen is therefore not sufficient to discard radiation, and the triggered
ODE sensitivity above is required.

The full-black ODE then gives corrected-369-point slopes `0.20788598`,
`0.05950367`, `0.02245214`, and `0.07706592`. All four remain inside 25%,
and the external-heat-corrected energy drift is `3.44e-5`. Compared with the
accepted adiabatic model, the centreline-temperature error improves by about
0.83 percentage point and the temperature-width error improves by only 0.014
point, while both mass errors worsen (about 0.12 and 0.034 point).

The blackbody energy ratios reach 1% only at effective absorptivities 0.724
and 0.592, corresponding to grey optical depths 1.286 and 0.895. At optical
depths 0.001/0.01/0.1/0.5/1.0, the larger warming-scale ratios are
0.0017%/0.0168%/0.161%/0.665%/1.069%.

**Decision:** reject radiation as a recommended correction. Even the
physically excessive perfect-black limit fails the pre-registered requirement
not to degrade the mass metrics. Keep `radiative_absorptivity` as an explicit
off-default external-energy sensitivity; do not infer it from the Raman
residual. The recommended dry-air model remains radiatively adiabatic.

### Final-journal source correction

The revised full-black errors are -24.978%, -15.824%, -20.552% and +24.541%.
The first value is only 0.022 percentage point inside the 25% band; without
fit uncertainty this 4/4 label has no meaningful selection margin. The
radiation rejection is unchanged. See
`hecht-panda-journal-benchmark-correction.md`.
