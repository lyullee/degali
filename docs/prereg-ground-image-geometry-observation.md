# Pre-registration: ground-image geometry observation operator

Date: 2026-09-05

## Why this audit is needed

The PRESLHY vertical centres and widths are not direct flow variables. They are
parameters of a single Gaussian fitted to four or five fixed sensor heights.
The validation concentration operator, however, evaluates a direct Gaussian
plus its mirror image below the ground. The current geometry comparison then
compares the fitted measurement with the unreflected internal trajectory
centre and width. Those are identical only when the mirror image is
negligible.

This distinction is largest in the 0.5 m releases. For example, trial 10 at
6 m has a measured fitted standard deviation of 0.594 m, comparable with and
larger than its fitted centre height of 0.227 m. It is therefore not valid to
assume before calculation that the ground image is negligible there.

This is an observation-operator audit. It changes no plume equation, source
state, coefficient, integration step, concentration prediction or selected
measurement.

## Frozen calculation

For every one of the same 17 well-constrained vertical fits used in the
independent-energy field result:

1. Retain the already selected trial and distance.
2. Take the exact on-axis sensor heights (`abs(y) < 0.01`) at that distance
   whose measured peak exceeded the original 0.05 vol-% fit-input floor.
3. Evaluate the model at those heights with its existing direct-plus-image
   concentration operator.
4. Fit those model values with the same single-Gaussian equation, initial
   width and parameter bounds used for the measured profile:

   ```text
   C(z) = C0 exp[-0.5 ((z-zc)/sigma_z)^2]
   C0 in [0,110], zc in [-1,6] m,
   sigma_z in [0.03/sqrt(2), 10/sqrt(2)] m.
   ```

5. Keep every original fit. Do not discard a candidate profile from its fit
   quality; report the minimum model-profile R-squared instead.
6. Apply this operator to both the corrected JETPLU baseline and the
   independent-energy candidate on identical `(trial, x)` keys.

The primary independent-energy integration remains the already frozen 0.02 m
conserved-flux RK4 calculation. No 0.01 m repeat is required because the
underlying state-convergence result is already recorded and this audit is a
deterministic post-processing operation.

## Predictions fixed before calculation

- The fitted model centre will differ materially from the internal trajectory
  centre only for sections whose lower Gaussian envelope overlaps the ground.
- The largest change will be in the 0.5 m releases, especially trial 10 at
  6 m. The 1.5 m releases should change much less.
- Re-fitting the direct-plus-image field should lower the apparent model
  centre for ground-overlapping profiles. It may also change the apparent
  width; its direction is not fixed in advance.
- Concentration MG, VG and FAC2 must be bit-for-bit unchanged because no
  concentration field is changed.

## Decision rule

The sensor-fit geometry operator replaces the internal-state comparison for
claims about what the PRESLHY sensor array observed if all 17 model profiles
fit successfully and their minimum R-squared is at least 0.85. If that gate
passes, promotion of the independent-energy model is re-evaluated using the
unchanged all-metric rule and the re-fitted centre MAE and mean width ratio.

If the fit gate fails, or if candidate centre MAE remains worse than the
baseline after applying the common operator, the earlier `do not promote`
decision stands. In either outcome, internal trajectory centre and width
remain available as model-state diagnostics and must not be relabelled as
sensor-fitted quantities.

## Result

All 17 baseline and candidate profiles fitted successfully. The minimum
model-profile R-squared is 0.984 for the baseline and 0.971 for the
independent-energy candidate, so the frozen observation-operator gate passes.
Concentration remains exactly MG 1.074/1.118, VG 1.213/1.176 and FAC2
0.952/0.976 for baseline/candidate, as required.

| common sensor-fit geometry | baseline | independent energy |
|---|---:|---:|
| sections | 17 | 17 |
| mean model/measured `sigma_z` | 1.360 | 1.420 |
| centre bias, m | -0.033 | +0.087 |
| centre MAE, m | 0.070 | 0.137 |
| minimum model-profile R-squared | 0.984 | 0.971 |

The prediction that the ground image matters most for the 0.5 m releases is
confirmed. For baseline trial 10 at 6 m, the internal centre of 0.541 m maps
to a sensor-fit centre of 0.150 m. The candidate profile at the same sparse
heights is broad enough that its fitted centre and width become 1.585 m and
2.550 m, despite an R-squared of 0.971; this also exposes weak parameter
identifiability when the instrument span is shorter than the predicted
profile width.

The candidate remains **not promoted**. Applying the common observation
operator does not remove its centre-height defect and makes its excessive
sensor-visible width explicit. The original internal-state result remains a
useful state diagnostic, while this result is the appropriate claim about
what the finite PRESLHY vertical arrays would observe.
