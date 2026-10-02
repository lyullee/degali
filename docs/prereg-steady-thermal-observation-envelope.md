# Pre-registration: steady thermal prediction against the observed time envelope

Date frozen: 2026-09-17, before calculating the comparison.

## Question

The downstream PRESLHY thermal residual changes substantially when a steady
model is compared with the observed minimum rather than its median. Is that
difference large enough in the existing paired raw records that it must be
separated from a mean thermal-transport error?

This is a validation-statistic audit. It changes no state equation, source
term, diffusivity, or default option.

## Frozen inputs

- Source: Lyons, Coldrick and Atkinson (2023), PRESLHY E3.5,
  [DOI 10.35097/1481](https://doi.org/10.35097/1481), CC BY-SA 4.0.
- Read the Trial 10 and 23 workbooks in place, using the already frozen
  Flexlogger windows 30--78 and 28--125 respectively.
- Use the co-located release-axis temperature-deficit trace at x=1.78 and
  4.00 m extracted by `read_paired_profiles`. Deficit is
  `max(T_ambient - T, 0)`.
- Retain only recorded samples with centre deficit at least 5 K, the same
  signal gate used for the five-point profile comparison.
- Compare the two already stored, replayed predictions without reintegrating:
  `density_profile_control` and `same_width_enthalpy_profile`.

## Frozen calculation

For every trial, station and model variant, report the selected sample count,
minimum, 5th/25th/50th/75th/95th percentiles, mean, sample standard deviation,
and model location within the empirical distribution. The latter is the
mid-rank:

```
(number strictly below + 0.5 * number equal) / number selected.
```

Also report prediction-minus-mean and prediction-minus-median and whether the
prediction lies inside the recorded 5--95% envelope.

## Interpretation fixed before calculation

- These statistics are empirical descriptions of a finite, serially
  correlated release window. They are not confidence intervals and no
  normality, independent-sample count, or sensor-time-constant correction is
  introduced.
- A value inside the 5--95% envelope does **not** validate a steady mean
  closure. A value outside it identifies an error that cannot be dismissed
  merely by choosing a different recorded instantaneous value.
- A positive prediction-minus-mean means that the model retains too large a
  temperature deficit, i.e. predicts a plume that is too cold. It must not be
  corrected by fitting a heat-transfer or diffusion coefficient from these
  four comparisons.
- The result is used only to decide whether the remaining downstream defect
  belongs to mean heat/density transport or unresolved temporal sampling.
