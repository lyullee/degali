# Steady thermal prediction versus observed time envelope

Date: 2026-09-17

## Result

The paired PRESLHY raw records confirm a directional downstream thermal error:
the density-profile control retains a **too-large temperature deficit** at
4.00 m. In physical temperature terms, it remains too cold; this is not a
case in which the model has warmed the plume too rapidly.

The conclusion comes from the recorded release-window distribution rather
than from a selected instantaneous minimum. It does not change a coefficient
or package default.

## Method and provenance

- Source data: Lyons, Coldrick and Atkinson (2023), PRESLHY E3.5,
  [DOI 10.35097/1481](https://doi.org/10.35097/1481), CC BY-SA 4.0.
- Trial 10 uses 49 one-second values and Trial 23 uses 98 values after the
  previously fixed 5 K release-axis signal gate.
- The raw workbooks were read in place and are not distributed. Their earlier
  recorded SHA-256 values remain unchanged.
- A steady value is placed in the observed record using a tie-aware empirical
  mid-rank; this is not a normal-distribution percentile or confidence
  interval. The frozen calculation is in
  [prereg-steady-thermal-observation-envelope.md](prereg-steady-thermal-observation-envelope.md).

## Direct comparison

Temperature deficit is `T_ambient - T`; a positive model-minus-observation
therefore means that the predicted temperature is too low.

| Trial | x (m) | observed mean / median deficit (K) | control deficit (K) | control percentile rank | control − mean (K) | result |
|---:|---:|---:|---:|---:|---:|---|
| 10 | 1.78 | 169.67 / 163.62 | 161.48 | 46.9% | -8.19 | central part of record |
| 10 | 4.00 | 56.60 / 53.76 | 94.93 | 79.6% | +38.33 | too cold downstream |
| 23 | 1.78 | 164.81 / 167.82 | 150.91 | 16.3% | -13.90 | too warm upstream |
| 23 | 4.00 | 69.30 / 64.69 | 83.26 | 78.6% | +13.96 | too cold downstream |

The independent same-width enthalpy-profile candidate has the same pattern
and increases the downstream cold bias: +43.31 K in Trial 10 and +18.06 K in
Trial 23 relative to their means. It remains rejected by the earlier joint
amplitude/variance criterion.

All four steady values happen to fall inside their wide finite-record 5--95%
envelopes. That fact does not erase the pattern: both controls move from a
low-to-central upstream rank to approximately the 79th percentile downstream.
For Trial 10 the observed 5--95% range at 4 m is 15.83--107.29 K, showing why
comparison with the single minimum alone would falsely exaggerate the mean
error while comparison with the median alone would conceal its time-domain
uncertainty.

## Source-variation-independent check

Using the pre-existing maximum-correlation thermal delay (1 s for Trial 10,
0 s for Trial 23), the simultaneous downstream/upstream centre-deficit ratio
is:

| Trial | valid paired seconds | observed ratio: p05 / median / p95 | control ratio | same-width candidate ratio |
|---:|---:|---:|---:|---:|
| 10 | 48 | 0.112 / 0.345 / 0.502 | 0.588 | 0.603 |
| 23 | 98 | 0.289 / 0.396 / 0.563 | 0.552 | 0.564 |

Trial 10’s control is above the observed 95th percentile even after the
upstream release amplitude is divided out. Trial 23 is closer but still above
the observed median. Thus a source-amplitude adjustment cannot explain the
common direction, and the same-width enthalpy candidate worsens both ratios.

## Physical disposition

The current thermodynamic path already contains humid-air H2O phase
equilibrium, including water ice/liquid saturation, latent heat, N2/O2/Ar
condensation, and a self-consistent humid ambient endpoint. The remaining
deficit is therefore **not** evidence that water condensation or its latent
heat was omitted. A zero-cost mean kinetic-energy-to-heat explanation was
also previously rejected: available mean kinetic energy accounts for only
0.015--0.123 K, versus multi-kelvin to tens-of-kelvin residuals.

The remaining physically admissible route is an explicitly parameterized
turbulent-energy/thermal transport closure: shear production must first enter
resolved `rho*k`, then only measured or independently bounded dissipation can
enter enthalpy. DEGALI has this conservation-preserving research operator,
but the PRESLHY records contain no velocity RMS, Reynolds stresses, integral
length scale, or dissipation measurement needed to set it without fitting.
It remains off by default.

This completes the present downstream mean-error diagnosis. No unmeasured
heat-exchange multiplier has been enabled merely to force the four profiles
into agreement.
