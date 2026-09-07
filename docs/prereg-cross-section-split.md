# Pre-registration: near-field cross-section split diagnostic

> **Superseded measurement scale.** This diagnostic triggered a subsequent
> validation audit: the measurement width below was the Gaussian e-folding
> width mislabeled as `sigma`, whereas JETPLU reports a standard deviation.
> The original registration and run are retained as an audit trail; the
> corrected interpretation is appended at the end.

Written after inspecting the current 23 PRESLHY vertical fits, but before
rerunning the coupled plume equations with the diagnostic below.  Do not edit
above the `RESULTS` line after the coupled comparison has run.

## Fault under test

`JETPLU` integrates only the product `sigma_y sigma_z`.  It then recovers the
two widths from

```text
sigma_y^2 - sigma_z^2 = sigma_ya^2 - sigma_za^2,
```

which is equivalent to adding the same jet-generated variance to the passive
horizontal and vertical variances.  The PRESLHY data resolve `sigma_z`, but do
not have enough lateral sensors at the same arcs to resolve `sigma_y`.

A direct audit of the corrected model gives a median
`sigma_y/sigma_z = 1.13` over the 23 fixed fits.  Merely replacing both widths
by `sqrt(sigma_y sigma_z)`, without reintegration, increases the median
vertical width only from 0.287 to 0.299 m.  That preliminary algebra says the
split can explain only a small part of the 27% mean vertical-width deficit.
It does not include the feedback of shape on perimeter entrainment, wind
averaging, trajectory and concentration, so the coupled upper-bound test is
still required.

## Diagnostic and predictions fixed in advance

The diagnostic forces

```text
sigma_y = sigma_z = sqrt(sigma_y sigma_z)
```

at every derivative evaluation.  It preserves the integrated area exactly,
adds no coefficient, and is deliberately an upper bound: it suppresses real
far-field atmospheric anisotropy and therefore is not a candidate production
closure by itself.

- **S1, split-primary hypothesis:** if the split is the main cause, the mean
  `modelled_sigma_z/measured_sigma_z` on the same 23 momentum-filtered fits
  must rise from 0.731 to at least **0.90**.
- **S2, area-primary diagnosis:** a result below **0.80** rejects S1 and points
  to insufficient total `sigma_y sigma_z` growth (or a measurement-model
  width-definition mismatch), rather than horizontal/vertical allocation.
- The 62 common PRESLHY concentration arcs will also be reported as
  MG/VG/FAC2.  They are a diagnostic side effect, not an adoption criterion.
- The historical/default path is not changed, and no coefficient is fitted.

## Decision rule

The circular split will not be adopted even if S1 passes, because it has the
wrong passive-plume asymptote.  A pass would justify implementing a published
near-to-far aspect-ratio transition.  Failure directs the next audit to total
entrainment, source-width/profile conversion, and the crossflow vortex-pair
entrainment term.

---

## RESULTS

The fully coupled circular-section run gives:

| configuration | vertical fits | mean `sigma_z` ratio | centre MAE (m) | concentration arcs | MG | VG | FAC2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| corrected default | 23 | 0.731 | 0.145 | 62 | 1.047 | 1.425 | 0.84 |
| circular diagnostic | 23 | **0.776** | 0.145 | 62 | 1.023 | 1.405 | 0.82 |

S1 fails and S2 passes.  Even the deliberately extreme circular split removes
only 4.5 percentage points of the 26.9-point vertical-width deficit.  The
trajectory is unchanged to the reported centre-error precision.  MG and VG
move slightly towards one while FAC2 loses one of 62 arcs; these small changes
do not rescue the geometric hypothesis.

## Decision

The circular split is rejected and is not added to the production code.  The
existing variance-addition split is not established as exact, but it is not
the main source of the PRESLHY vertical-width error.  The next audit is the
growth of total `sigma_y sigma_z`, including source-profile conversion and the
crossflow/vortex entrainment missing from a straight-jet shear closure.

## Subsequent Gaussian-convention audit

The next audit found that the apparent deficit was not a total-area error.
The PRESLHY fit used `exp(-(z/w)^2)` and stored `w` as `sigma_z`; JETPLU uses
`exp(-0.5*(z/sigma_z)^2)`. Thus the measurement standard deviation is
`w/sqrt(2)` and every model/measurement ratio above must be multiplied by
`sqrt(2)`:

| configuration | common-definition ratio |
|---|---:|
| corrected default | **1.033** |
| circular diagnostic | **1.097** |

There is no 27% width deficit to repair. The circular diagnostic remains
rejected, now because it increases a small +3.3% mean bias to +9.7%. See
`gaussian-width-convention.md`.
