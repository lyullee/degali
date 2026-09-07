# Hecht--Panda final-journal benchmark correction

Date: 2026-09-04

## Outcome

The active Raman benchmark now uses the numerical fit labels in the final
2019 journal manuscript, not the earlier ICHS 2017 conference manuscript.
This is a source correction, not a model calibration: no predicted trajectory,
closure coefficient or thermodynamic result was changed.

The correction preserves the recommended dry conserved-axisymmetric research
model, but materially weakens how confidently its binary 4/4 result can be
stated. Two of its four errors are now less than half a percentage point from
the pre-registered 25% boundary, and the paper supplies neither uncertainty
bars for the aggregate slopes nor unambiguous membership of the fitted cases.
The result must therefore be described as a **provisional aggregate-fit
check**, not a raw-data validation.

## Versioned observations

| response | ICHS 2017 | final journal 2019 | active value |
|---|---:|---:|---:|
| inverse centreline H2 mass-fraction slope | 0.2626 | 0.2771 | **0.2771** |
| H2 half-width slope | 0.06503 | 0.07069 | **0.07069** |
| inverse centreline-temperature slope | 0.02826 | 0.02826 | **0.02826** |
| temperature half-width slope | 0.06188 | 0.06188 | **0.06188** |

The final values are printed directly in Figures 6 and 8 of the Sandia/OSTI
author manuscript for the journal article. The manuscript's surrounding text
rounds the two revised mass values to 0.277 and 0.071. No curve digitisation
was used.

## Unresolved condition-count provenance

Table 1 contains nine source conditions and the radial-profile figures contain
nine condition panels. In contrast, both aggregate centreline/half-width
figures list a tenth series, `4 bar, 45 K, 1.25 mm`. The prose still describes
fits to nine conditions and the conclusion gives a 50--61 K nozzle range,
which excludes 45 K.

The extra condition cannot be reconstructed safely: Table 1 does not give its
number of stitched heights, throat temperature, throat pressure, throat
density or throat velocity. It has therefore **not** been invented or copied
from the 4 bar, 54 K row. The model predictions cover the nine fully tabulated
conditions. Whether each published aggregate fit used nine or ten series is
unresolved until the fit inputs or reduced arrays are obtained from the
authors.

## Density-normalization audit

Both centreline regressions normalize distance with the nozzle stagnation
density `rho0`, but the manuscript does not state the property routine used
for that density. Replacing the real-gas convention with an ideal-gas value
would be an attractive but unsupported way to move a borderline slope.

The independent check is Table 1 itself. At the nine printed throat
temperature/pressure pairs, CoolProp hydrogen densities reproduce the printed
throat densities with 0.317% RMS error and 0.550% maximum error. Ideal-gas
densities have 2.73% RMS error and 4.84% maximum error. The table therefore
supports retaining the existing real-gas normalization. Across the stagnation
states the real/ideal density ratio is 1.008--1.042; this ambiguity is too small
to explain the thermal residual and is not used as a tuning parameter.

## Effect on the corrected 369-point comparison

| candidate | mass centre | mass width | T centre | T width | within 25% |
|---|---:|---:|---:|---:|---:|
| recommended dry component-enthalpy/phase model | -24.86% | -15.79% | -21.38% | +24.56% | **4/4** |
| argon phase sensitivity | -24.64% | -15.72% | -22.49% | +23.90% | 4/4 |
| perfect-black radiation bound | -24.98% | -15.82% | -20.55% | +24.54% | 4/4 |
| para-H2 enthalpy sensitivity | -23.67% | -14.85% | -34.79% | +24.58% | 3/4 |
| four-flux phase/two-scalar closure | -24.21% | -6.12% | -29.14% | -2.27% | 3/4 |
| independent HyRAM 6.1 oracle | -8.82% | -18.21% | -20.21% | +14.78% | 4/4 |

The model-selection conclusions do not reverse. The dry component-enthalpy
candidate remains the recommended conservative research closure; para-H2 and
the four-flux candidate still fail the centreline-temperature criterion.
However, the perfect-black case is only 0.022 percentage point inside the
mass-centre threshold, so its 4/4 label carries essentially no discriminating
weight without measurement/regression uncertainty.

## Required resolution

Request the exact row or run record for `4 bar, 45 K, 1.25 mm` and, separately
for each of the four aggregate fits, the included run IDs, fitted intercept,
weighting and uncertainty/covariance. Reduced median fields remain preferred
because they permit a like-for-like re-fit of the nine tabulated conditions.

Until those records arrive:

1. retain both publication versions in the machine-readable transcription;
2. use the final journal values as the active benchmark;
3. do not add a guessed tenth model condition;
4. report 4/4 only with the provisional aggregate-fit qualifier; and
5. do not tune a physical coefficient to the remaining sub-threshold margins.
