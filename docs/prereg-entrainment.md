# Pre-registration: the entrainment shortfall in the cryogenic jet

**Registered before any of the tests below were run.**
Do not edit above the `RESULTS` line once results are known.

---

## Background

**Superseded — recomputed from the reduction; the figures below were never computed by anything.** Recomputing gives MG 0.855, CI [0.672, 1.064], VG 1.30, FAC2 0.85 as shipped and MG 1.070 corrected, over 66 arcs, and the interval now includes 1. See `docs/lh2-recomputed.md`.
Against the PRESLHY E3.5 raw data, momentum-driven releases give MG 0.738 with
a 95 % interval of [0.591, 0.918] over 69 arc maxima. The model reads high by
about a third, and the residual grows with distance:

| x (m) | observed | model | ratio |
|---|---|---|---|
| 0.35 | 84.2 | 96.7 | 1.15 |
| 1.19 | 78.9 | 82.8 | 1.05 |
| 2.67 | 44.8 | 64.1 | 1.43 |
| 6.00 | 22.1 | 32.5 | 1.47 |

A bias that is small at the source and grows downstream is what too little
entrainment looks like: the concentration is right where nothing has been
entrained yet, and increasingly too high as the deficit accumulates.

`JETPLU` entrains through

    E = alfa1 |u_c| + alfa2 |u_a sin(theta)| + (curvature and Richardson terms)

with `alfa1 = 0.057` and `alfa2 = 0.5` from the routine's own `DATA`
statement. Those are Ooms' values for a buoyant jet in a crossflow, fitted on
ordinary-temperature gases.

**Four sessions of this work have produced four "model defects" that turned
out to be errors in how the data was being compared.** So the first question
is not which coefficient to move; it is whether the shortfall survives the
comparisons that have already caught four such errors.

## What is being tested

**The diagnosis, not a fitted correction.** No coefficient will be adjusted to
improve agreement. The tests ask whether the residual behaves the way an
entrainment deficit must behave, and whether any *published* alternative
closure removes it.

### Predictions

| # | Prediction | Criterion |
|---|---|---|
| **P-E1** | The residual is a distance effect, not a per-trial offset | fitting `ln(model/obs) = a + b ln(x)` over the 69 pairs gives **b > 0** with a bootstrap interval excluding zero |
| **P-E2** | It is not the source term | recomputing with the flow-meter peak instead of the window mean, and with the report's Table 4 rate, moves MG by **< 0.10** |
| **P-E3** | It is not the ambient dispersion parameters | switching stability class D to C and to E moves MG by **< 0.10** |
| **P-E4** | It is not specific to hydrogen | the same distance-growing residual appears on the SMEDIS ammonia jets (FLADIS, Desert Tortoise), **b > 0** there too |
| **P-E5** | A published cross-flow coefficient in the accepted range removes most of it | with `alfa2` at Ooms' upper published value, MG moves **toward 1 by at least half the gap**, without VG worsening by more than 20 % |
| **NC1** | Negative control: the EPA test cases are unaffected by anything not enabled | B9, B9T, EX1, EX2, EX3 reproduce to 1e-12 throughout |
| **NC2** | Negative control: a coefficient change must move the answer | setting `alfa1` to zero must change the LH₂ result substantially, or the sensitivity test is not testing what it claims |

### What would falsify the diagnosis

- **P-E1 fails** — the residual is a constant offset, so it is a source or
  normalisation problem and the entrainment reading is wrong.
- **P-E2 or P-E3 fails** — the shortfall is an artefact of an input choice.
- **P-E4 fails** — it is specific to cryogenic hydrogen, and the mechanism is
  more likely the droplet phase than the entrainment coefficient.
- **P-E5 fails** — the deficit cannot be expressed as a coefficient within the
  published range, and any fix would be a refit rather than a substitution.

### What will not be done

No value of `alfa1` or `alfa2` outside the range those coefficients are
published with. No optimisation against the PRESLHY data. If P-E5 fails, the
result is recorded as an open shortfall and the coefficients stay as DEGADIS
has them.

---

## RESULTS

**P-E1 fails. The entrainment diagnosis is wrong.**

| test | result | verdict |
|---|---|---|
| **P-E1** distance slope | b = **+0.018**, CI [−0.168, +0.171] | **FAIL** |
| P-E2 flow-meter peak instead of window mean | MG 0.702, Δ −0.036 | pass |
| P-E2 report Table 4 rate | MG 0.692, Δ −0.046 | pass |
| P-E3 stability C | MG 0.824, Δ +0.086 | pass |
| P-E3 stability E | MG 0.706, Δ −0.032 | pass |
| NC2 `alfa1 = 0` | MG 0.614, Δ −0.124 | pass |

Regressing `ln(model/observed)` on `ln(x)` over all 69 pairs gives a slope
indistinguishable from zero. **The residual does not grow with distance.**

The apparent growth — 1.15 near the source to 1.47 at 6 m — came from a table
of medians over a seven-trial subset. It does not survive a regression on the
full set. That is the fifth time in this work that a pattern seen in a
convenient summary has failed on the underlying pairs, and it is exactly what
the pre-registration was for: the criterion was fixed before the regression
was run, so there was nothing to negotiate afterwards.

P-E4 and P-E5 were not run. Both are conditional on P-E1: there is no
distance-growing deficit to look for in the ammonia jets, and no reason to
try a cross-flow coefficient against a residual that does not behave like an
entrainment shortfall.

### What the residual is instead

A constant factor of about 1.35, independent of distance, and not explained by

- the source rate (P-E2: three different rates, all within 0.05 of each other),
- the ambient dispersion parameters (P-E3: stability C to E spans 0.12),
- nor by the entrainment coefficient the model already has (NC2 shows the
  machinery responds, so a null result elsewhere is informative).

A distance-independent factor points at the source or at normalisation rather
than at anything that accumulates downstream. Candidates, none tested:

- the profile constants `rk1`, `rk2`, which set the ratio between the
  contaminant flux and the centreline value for the assumed Gaussian shape;
- the reduction on the measurement side, where an arc maximum over a sparse
  array is a lower bound on the true centreline whenever the plume is not
  exactly over a sensor;
- the flow meter, which the PRESLHY report classifies as reading a "gassy"
  two-phase stream for several of these nozzles.

The second of those would mean the model is not biased at all and the
observations are low. That cannot be resolved with an array of this density,
and it should not be argued either way without a test that can distinguish
them.

**No coefficient was changed.** The pre-registration said that if the
diagnosis failed the coefficients stay as DEGADIS has them, and they do.
