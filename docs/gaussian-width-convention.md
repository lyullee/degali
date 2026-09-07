# Gaussian width convention correction

## Finding

The PRESLHY fitting code and `JETPLU` used the same Gaussian family but named
different width parameters `sigma`.

The former measurement fit was

```text
C(z) = C0 exp(-((z-zc)/w)^2),
```

while `JETPLU.FOR` and `Trajectory.concentration_at` use

```text
C(z) = C0 exp(-0.5*((z-zc)/sigma)^2).
```

They are identical only when

```text
sigma = w/sqrt(2).
```

The fit code stored `w` under the field names `sigma_z` and `sigma_y`, then
compared it directly with the model's true standard deviation.  This created
an exact factor-`sqrt(2)` apparent width deficit.

## Correction

- `fit_vertical` and `fit_arcs` now fit the standard Gaussian form with the
  factor `0.5` in the exponent.
- Their initial guesses, bounds and the broad-fit rejection limit were all
  divided by `sqrt(2)`, so the fitted curves and accepted observations are
  unchanged; only the reported width coordinate changes.
- The 72 vertical and 65 lateral fitted widths in
  `reference/preslhy/e35_reduced.json` were divided by `sqrt(2)` and the file
  now declares its formula in `gaussian_width_definition`.
- `vertical` remains backward-compatible with legacy reductions that lack the
  metadata: it treats their stored value as the old e-folding width and
  converts it once on read.

No plume equation, source term, coefficient, concentration value, fitted
centre, fit quality or selected sample changed.

## Corrected result

On the same 23 well-constrained, momentum-driven PRESLHY fits:

| model | legacy comparison | common-standard-deviation comparison |
|---|---:|---:|
| historical reconstruction | 0.637 | **0.901** |
| mass/enthalpy/momentum-consistent source | 0.731 | **1.033** |

The corrected model therefore has a mean vertical-width bias of +3.3%, not a
26.9% deficit.  The median per-fit ratio is 1.042.  The earlier conclusion
that missing entrainment or cross-section anisotropy was required to repair
`sigma_z` is withdrawn.

The previously registered circular-section upper-bound diagnostic changes
from 0.776 on the legacy scale to 1.097 on the correct scale.  It still shows
that forcing a circular section is unnecessary and would now over-widen the
plume.

## Remaining scope

This correction does not alter the already computed concentration statistics
or centre-height errors.  It also does not prove the lateral plume width,
because the far-field lateral fits contain wind-direction meander as well as
turbulent spread.  It removes one false target from the physics backlog; air
condensation/slip and independent far-field trajectory validation remain open.
