# Pre-registration: Raman case-coverage correction

**Frozen before the first corrected-coverage model result.** Do not edit above
the `RESULTS` line after a result is known.

## Source correction

The first Raman protocol sampled every one of the nine cases uniformly from
40 through 100 mm, producing 549 synthetic model points. That is not the
experimental weighting reported by Hecht--Panda. Table 1 gives `n_heights`
from 2 to 6; the methods state that the nozzle was moved in 10 mm increments
and images at each height were stitched. Figure 2 labels a six-height stitched
image from 40 to 100 mm. Thus a case with `n_heights=n` covers

```
z = 40, 41, ..., 40 + 10 n  mm.
```

The nine coverages contain 369 one-millimetre model samples, preserving the
reported unequal experimental support. This correction is determined entirely
from Table 1, the methods text and the labelled six-height image; no model
residual enters it.

## Frozen comparison

Reapply the same fitted-intercept centreline regressions and forced-origin
half-width regressions to:

1. the independent reproduction of the published Gaussian establishment;
2. the four-flux conservative `entrained_mass` establishment;
3. the dew-point-heated, independent `Y`--`T` two-scalar candidate.

All source equations, coefficients and integration tolerances remain
unchanged. Observations remain the numerical fits printed in Figures 5--8.

## Decision rule

The corrected coverage becomes the primary aggregate protocol if and only if
the generated endpoints agree with every Table-1 `n_heights` and the Figure-2
six-height 40--100 mm example. Adoption does not depend on improving any
model. Report both old and corrected results; do not erase the earlier frozen
549-point audit.

---

## RESULTS

Adopted. The generated endpoints agree with all nine Table-1 `n_heights`
values and reproduce the Figure-2 40--100 mm extent for six heights. The
corrected protocol contains 369 samples rather than 549.

The independently reproduced published establishment predicted slopes
`0.24916660`, `0.05738458`, `0.02275205`, and `0.07174457`: all four are
within 25%, although its median mass radial coefficient (64.91) is just
outside the reported case range ending at 64. The conservative establishment
predicted `0.20978143`, `0.06142980`, `0.01936885`, and `0.07780172` (two of
four within 25%). The heated two-scalar candidate predicted `0.20974082`,
`0.06705702`, `0.01711375`, and `0.06085325` (three of four). The correction
does not rescue either rejected candidate and the original 549-point audit is
retained for traceability.

### Final-journal source correction

The unequal 369-point sampling remains the primary model protocol, but the
active observations are now the final-journal fits. The source also leaves it
unclear whether an untabulated tenth series contributed to either aggregate
figure. See `hecht-panda-journal-benchmark-correction.md`; do not interpret
369 as a recovered experimental point count until the original fit inputs are
available.
