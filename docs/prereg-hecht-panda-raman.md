# Pre-registration: Hecht--Panda cryogenic-H2 Raman jet

**Frozen before running DEGALI on these conditions.**  Do not edit above the
`RESULTS` line after the first model run.  Corrections and interpretation go
below that line.

## Purpose

Hecht and Panda measured simultaneous hydrogen concentration and temperature
fields from nine 50--61 K, 2--5 bar absolute jets through 1.00 and 1.25 mm
orifices.  The paper reports aggregate centreline-decay and radial-width fits,
so it can independently test the near-nozzle mixing and warming closures that
the outdoor PRESLHY arrays cannot resolve.

This is a cold **gaseous, momentum-dominated** source test.  It does not test
the liquid flash, condensed-air particle transport, ground interaction, wind
bending, or the far-field DEGADIS handoff.

## Frozen observations

The machine-readable transcription is
`reference/lh2/hecht-panda-raman-2017.json`.  The comparison uses only values
printed in Table 1 and Figures 5--8:

| response | reported aggregate relation |
|---|---:|
| centreline H2 mass fraction | `1/Y_cl = 0.2626 z/(d sqrt(rho0/rhoa))` |
| H2 half-width | `r_1/2(mm) = 0.06503 z/r_noz` |
| radial H2 profile | `Y/Y_cl = exp(-49 (r/z)^2)`; case fits 33--64 |
| centreline temperature | `Tinf/(Tinf-T_cl) = 0.02826 z/(d sqrt(rho0/rhoa))` |
| temperature half-width | `r_1/2(mm) = 0.06188 z/r_noz` |
| radial temperature profile | `Theta = exp(-42 (r/z)^2)`; case fits 21--49 |

The aggregate slope fits are digitisation-free.  Figure contours are retained
as qualitative corroboration only because no numerical image data or
uncertainty table accompanies the paper.  The authors explicitly note ice
deflection, possible out-of-plane offsets, and appreciable temperature noise.

## Frozen source reconstruction

1. The measured/calculated throat values in Table 1 (`T_throat`, `P_throat`,
   `rho_throat`, `v_throat`) define one internally consistent choked plane.
   Mass flow is `rho_throat v_throat pi d^2/4`.  Before the ambient-pressure
   plume calculation, that given plane is expanded with the Yuceil--Otugen
   mass, pressure-thrust and total-energy balances (HyRAM+ equations 47--51).
   The printed throat state is retained; it is not replaced by a new choked-
   flow calculation.
2. Ambient air is 295 K and 1 atm.  The 0.3 m/s coflow is less than 0.07% of
   the 498--559 m/s throat speed and is neglected in the momentum balance.
3. The model starts after its existing 7.7-effective-diameter Gaussian-
   development length.  The expansion plane has negligible assigned axial
   length, and model coordinates are converted back to physical distance from
   the physical orifice before comparison.
4. Only the measured range `z = 40--100 mm` is scored, sampled every 1 mm.
5. Stagnation density `rho0` is evaluated from the printed nozzle pressure and
   temperature with CoolProp.  The table's throat properties are not replaced
   by a new choked-flow calculation.
6. A model mass-fraction half-width is the radius at `Y/Y_cl = 0.5`.  Since
   JETPLU carries `exp[-r^2/(2 sigma^2)]`, this is
   `sqrt(2 ln 2) sigma`.  Temperature half-width is found from the adiabatic
   table where `(T-Tinf)/(T_cl-Tinf) = 0.5`; it is not assumed equal to the
   concentration width.

## Frozen model variants

- `legacy`: DEGADIS JETPLU constants, `alpha1=0.057`, `Sc=1.42`, no density
  scaling.
- `current_lh2`: the present corrected LH2 choice, `alpha1=0.0875`, `Sc=1.42`,
  Ricou--Spalding density scaling.
- `alpha_only`: `alpha1=0.0875`, `Sc=1.42`, without density scaling.  This is
  a mechanism isolation, not a fitted candidate.

No coefficient will be optimized against these nine cases in this phase.

## Predictions and decision rules

1. A credible variant must put both aggregate centreline slopes and both
   aggregate half-width slopes within 25% of the printed values.
2. Its equivalent radial Gaussian coefficients must fall inside the printed
   case ranges (33--64 for H2 and 21--49 for temperature).
3. Density scaling is supported only if `current_lh2` lowers the cross-case
   normalized-trend RMSE for centreline mass fraction relative to
   `alpha_only`; matching one aggregate slope is insufficient.
4. A change is eligible for the corrected production path only if it improves
   at least three of the four slope errors, worsens none by more than five
   percentage points, leaves the Fortran-parity path untouched, and does not
   degrade the existing PRESLHY concentration statistics.
5. If no non-fitted published mechanism satisfies rule 4, the result is an
   identified model-form shortfall.  The production constants stay unchanged.

---

## RESULTS

The first run used 61 points per case over 40--100 mm (549 points per
variant):

| variant | mass decay | mass width | T decay | T width | `A_Y` | `A_T` | slopes within 25% |
|---|---:|---:|---:|---:|---:|---:|---:|
| observed | 0.26260 | 0.06503 | 0.02826 | 0.06188 | 49 | 42 | - |
| legacy | 0.34296 | 0.11519 | 0.05813 | 0.13506 | 16.3 | 11.6 | 0/4 |
| current LH2 | 0.58235 | 0.19238 | 0.08448 | 0.21382 | 5.8 | 4.6 | 0/4 |
| alpha only | 0.52407 | 0.17377 | 0.07809 | 0.19496 | 7.2 | 5.6 | 0/4 |

The current LH2 path fails all four slope criteria.  Its centreline mass-trend
relative RMSE is 1.220, compared with 0.987 without density scaling, so frozen
rule 3 rejects the Ricou--Spalding local-density scaling for this cold
free-jet application.  Even the legacy constants overpredict mass width by
77% and temperature width by 118%; the deficiency is therefore not repaired
by reverting only the entrainment coefficient.

This combination - near-er mass decay but much too-large width - is the same
structural distinction identified by Hecht and Panda.  A single mixing-rate
coefficient cannot independently set the velocity-profile width relative to
the scalar profile.  The next phase must therefore test the published/profile
width-ratio closure and a conserved total-energy equation separately.  No
production constant is changed by this failed first phase.

### Regression-definition correction

The table above is retained as the first-run audit trail, but its two
centreline numbers are not valid comparisons.  The plot labels print a slope
(`fit proportional to ...`), not a zero-intercept equation.  Figure 7 makes
the distinction decisive: `Tinf/(Tinf-T_cl)` approaches one at a pure-source
limit and cannot pass through zero.  The implementation therefore now uses an
ordinary straight-line slope with a fitted intercept for the two centreline
decays.  Half-widths remain forced through zero.  The cross-case trend RMSE is
the relative residual about each variant's own fitted line, because the paper
does not print the fitted intercept.  Corrected results follow in the next
pre-registration document; the frozen candidate rules are unchanged.

The completed corrected comparison is:

| variant | mass decay | mass width | T decay | T width | `A_Y` | `A_T` | slopes within 25% |
|---|---:|---:|---:|---:|---:|---:|---:|
| observed | 0.26260 | 0.06503 | 0.02826 | 0.06188 | 49 | 42 | - |
| legacy | 0.35457 | 0.11519 | 0.04304 | 0.13506 | 16.3 | 11.6 | 0/4 |
| current LH2 | 0.57855 | 0.19238 | 0.06819 | 0.21382 | 5.8 | 4.6 | 0/4 |
| alpha only | 0.53746 | 0.17377 | 0.06350 | 0.19496 | 7.2 | 5.6 | 0/4 |
| HyRAM momentum only | 0.20827 | 0.06533 | 0.02748 | 0.08208 | 52.0 | 31.8 | 3/4 |
| conserved energy | 0.25456 | 0.05774 | 0.02264 | 0.07098 | 64.2 | 41.2 | **4/4** |

The baseline/current/alpha variants remain rejected after the correction. The
published source-momentum closure repairs the mass profile but not the thermal
width; the separately pre-registered total-energy model is the first to pass
all four slopes. See `prereg-hecht-panda-hyram-closure.md` for its construction,
official-oracle check and acceptance boundary.

### Final-journal source correction

The frozen observations and historical decisions above used the ICHS 2017
fit labels. The final 2019 article prints revised mass slopes 0.2771 and
0.07069; it also shows an otherwise untabulated tenth series in the aggregate
legends. Active scores and the resulting evidence limitation are superseded by
`hecht-panda-journal-benchmark-correction.md`.
