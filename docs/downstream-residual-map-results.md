# Downstream thermal and energy residual map

Date: 2026-09-08

## Outcome

The dominant cold-core discrepancy cannot be supplied by direct conversion of
the resolved mean kinetic energy. Across the positive-gap centreline points,
the conditional enthalpy requirement is approximately 150 to 1,153 times the
local mean kinetic energy. The already accounted mechanical-energy input has
a uniform sensible-temperature scale of only 0.015 to 0.123 K.

No trajectory was reintegrated. The result combines sealed pressure-loss
source trajectories and independent diagnostics under the rules fixed in
[`prereg-downstream-residual-map.md`](prereg-downstream-residual-map.md).

## Distance-resolved result

Positive temperature residual means that the model is warmer than the
observation.

| trial | x (m) | model - observed minimum (K) | model - observed median (K) | minimum-basis enthalpy gap / mean KE | accounted mechanical scale (K) |
|---:|---:|---:|---:|---:|---:|
| 10 | 1.19 | 36.75 | 20.95 | 1065.75 | 0.0571 |
| 10 | 1.78 | 57.96 | 3.19 | 1152.61 | 0.0435 |
| 10 | 2.67 | 24.21 | -22.27 | 771.41 | 0.0319 |
| 10 | 4.00 | 19.64 | -40.16 | 1074.12 | 0.0227 |
| 10 | 6.00 | -9.53 | -48.45 | 0 | 0.0152 |
| 23 | 0.79 | 7.20 | 1.86 | 150.40 | 0.1227 |
| 23 | 1.19 | 33.24 | 29.94 | 296.77 | 0.0941 |
| 23 | 1.78 | 41.68 | 18.19 | 631.42 | 0.0695 |
| 23 | 2.67 | 26.76 | 6.54 | 733.53 | 0.0503 |
| 23 | 4.00 | 19.51 | -17.44 | 938.78 | 0.0349 |
| 23 | 6.00 | 13.28 | -29.69 | 1134.48 | 0.0230 |

The zero ratio at Trial 10, 6 m means that the minimum-basis model is already
colder than the observation; there is no positive missing-enthalpy requirement
at that point.

## What the sign changes mean

Trial 10 at 2.67 and 4 m, and Trial 23 at 4 and 6 m, change residual sign when
the comparison moves from the time-window minimum to the median. A single
steady prediction is warmer than the cold transient but colder than the
typical time-window state. Averaging these values into one temperature target
would hide a first-order measurement/model mismatch.

This pattern supports separating:

1. a transported mean thermal/species profile, which controls the steady
   cross-section shape;
2. unresolved intermittency and sensor response, which control minimum versus
   median observations;
3. trajectory buoyancy, which must not be tuned with the same temperature
   statistic.

## Buoyancy relation in Trial 10

The sealed Trial-10 momentum audit gives cumulative candidate buoyancy of
0.178 N at 1.78 m, 2.110 N at 4 m and 5.695 N at 6 m. The thermal residual
does not remain one-signed while this upward impulse grows. A single global
heat offset or buoyancy multiplier therefore cannot repair the full trajectory
without changing the wrong parts of the record.

## Initial TKE and normal-stress bounds

The existing order-16 feasibility calculation is preserved, not converted
into an initial condition:

| trial | mean kinetic flux (W) | shear-bound flux (W) | minimum TKE fraction when normal-work fraction <= 0.05 | minimum normal-work fraction when TKE fraction <= 0.12 |
|---:|---:|---:|---:|---:|
| 10 | 347.28 | 36.67 | 0.1365 | 0.0630 |
| 23 | 521.89 | 59.25 | 0.1539 | 0.0811 |

These are compatibility bounds. They show that an arbitrarily tiny transported
TKE state is not consistent with the sealed joint momentum/TKE construction,
but they do not identify the actual experimental value.

## Model decision

- Reject resolved mean-KE thermalisation as the main cold-core correction.
- Do not introduce a single fitted heat offset or buoyancy multiplier.
- Keep minimum, p05 and median observations separate.
- The next bounded implementation should transport an independent thermal
  profile moment and retain TKE/normal stress as an uncertainty interval until
  a measured velocity boundary is available.

The machine-readable local audit is
`reference/preslhy/downstream_residual_map_2026-09-08.json`. It includes input
hashes and is excluded from the distributed package.

