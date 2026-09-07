# Measuring the sub-models, not just the answer

> **Historical analysis.** Current results after flash-state and source-
> momentum conservation are in `lh2-model-improvements-2026-09-03.md`.
> The width fits in this analysis also used an e-folding parameter under a
> sigma label. The factor-`sqrt(2)` correction and current 1.033 ratio are in
> `gaussian-width-convention.md`; width-deficit conclusions below are retired.

A concentration comparison gives one number per point, and everything the
model does wrong arrives in it mixed together. On the PRESLHY hydrogen data
that mixing produced a signal no single coefficient can make: **MG 0.74 near
the source and 3.6 in the far field** — a third high, then several times low.

The near-field array carries four or five heights at each distance and the far
field five stands across each arc. That is enough to fit the quantities the
model actually computes — the plume's centre, its vertical spread, its lateral
spread — and to test each one against measurement on its own.

## Three sub-models, three verdicts

### Lateral spread: correct

| | |
|---|---|
| measured `sigma_y` at 14 m | 4.35 m |
| modelled | 1.14 m |

A factor of nearly four, and none of it is a model error. The measurement is a
maximum over a release during which the wind direction wandered; the fitted
width therefore contains the meander. Removing it,

    sigma_y(total)^2 = sigma_y(plume)^2 + (x sigma_theta)^2

a `sigma_theta` of 15 to 17 degrees leaves a plume width of 0.8 to 2.2 m, and
the model's 1.14 m sits inside that. The 17 degrees is itself derived from
these fits rather than assumed, and is consistent with the compass record.

### Vertical spread: 1.4 times too narrow

The vertical carries no meander term — a plume does not wander up and down
with the wind direction — so the fitted `sigma_z` is directly comparable.

| x (m) | measured | modelled | ratio |
|---|---|---|---|
| 0.6–1.5 | 0.14 | 0.09 | 1.4 |
| 1.5–3.0 | 0.34 | 0.18 | 1.9 |
| 3.0–5.0 | 0.53 | 0.39 | 1.4 |
| 5.0–7.0 | 0.91 | 0.63 | 1.4 |

Consistent across the range, so it is not an accumulating error.

**It cannot be fixed with the entrainment coefficient alone.** Papanicolaou
and List (1988) measured 0.0545 for a pure jet and 0.0875 for a pure plume;
`JETPLU` uses 0.057, at the jet end and inside the measured range. Sweeping
`alfa1` from 0.052 to 0.100 — past the plume value — moves the ratio only from
0.59 to 0.78. `JETPLU` sets `sigma_z` by splitting the product
`sigma_y sigma_z` in the ratio of the *ambient* dispersion parameters, and at
0.35 to 6 m those are tiny numbers whose ratio nonetheless governs the split.

### Trajectory: the dominant fault

| x (m) | measured rise | modelled rise |
|---|---|---|
| 0.6–1.5 | +0.00 | 0.00 |
| 1.5–3.0 | −0.01 | +0.03 |
| 3.0–5.0 | +0.06 | +0.40 |
| **5.0–7.0** | **−0.12** | **+1.07** |

**The measured plume does not rise.** Over six metres its centre stays within
a tenth of a metre of the release height. The model lifts it by more than a
metre.

That gap is about 1.5 times the measured `sigma_z`, so the predicted
concentration at sensor height is cut by `exp(-1.5**2)` — a factor of ten.
That is the far-field MG of 3.6. Near the source the sensors are still inside
the plume, the trajectory error has not yet developed, and the narrow spread
dominates instead, giving MG 0.74. **One model, two faults, opposite signs,
and they cancel where the array is densest.**

## Cross-checks in the literature

**Papanicolaou and List (1988)**, read at source, derive `alpha = 0.0545` for
the jet, `alpha = 0.0875` for the plume and `Ri_p = 0.716` from their own
measurements. The frequently repeated 0.0533, 0.0833 and `Ri_p = 0.557` are
Fischer et al. (1979)'s proposed values, quoted for comparison in the paper;
the earlier secondary-source attribution is superseded. DEGADIS's 0.057 is
still at the jet end, so the sweep above remains a check rather than a fit.

**Xiao and co-workers (2009)**, a non-Boussinesq integral model built for
hydrogen safety, report that the normalised trajectory "will not collapse when
the Froude number is small, which means the Boussinesq approximation is
invalid when the buoyancy effect is comparable with the momentum effect".
A flashing LH₂ jet leaves the orifice at 4.4 times ambient density and is
buoyant within a metre — squarely in that regime. Their remedy is an
entrainment law in the Richardson number and the trajectory angle.

**Winters and Houf**, whose model is in HyRAM, divide a cryogenic release into
four zones — underexpanded flow, initial entrainment and heating, flow
establishment, established flow. `JETPLU` has no such division.

So the trajectory fault is the one the hydrogen literature predicts: a
Boussinesq buoyancy term applied where the density ratio is four.

## What was tried, and what is left

Each candidate was registered before it was run
(`prereg-trajectory.md`, `prereg-boussinesq.md`) and tested against the rise
rather than against a concentration that three faults contribute to.

| candidate | effect on the rise at 5–6 m |
|---|---|
| source state — mixing line from the flashing source | −0.22 m |
| form drag | 0.00 m |
| **non-Boussinesq buoyancy** | **+0.27 m, the wrong way** |
| ground effect | 0.00 m |
| droplet mass | −0.04 m |
| **needed** | **−0.93 m** |

The non-Boussinesq result is the informative one. It is the correction the
hydrogen literature points at, the approximation it targets is genuinely
present, and applying it makes the model worse — because dividing by the
parcel density rather than the ambient makes a *light* plume rise faster. That
was predicted from the algebra before the run, which is why a 1.34× worse
answer could not be read as needing tuning.

What is left is the slip between the droplets and the gas. An integral model
of this kind carries one velocity for the cross-section, so it can represent
the mass the droplets add — which is why `liquid_fraction` helps a little —
and not their falling. **The trajectory error is probably not fixable within
the formulation**, and closing it needs a second velocity field for the
dispersed phase.

## Why this is worth doing

Every fault above is now a *measured target* rather than a residual to be
argued about: keep `sigma_y`, widen `sigma_z` by 1.4, and suppress the rise.
Any change can be checked against the quantity it was meant to affect instead
of against a concentration that three faults contribute to.

It also explains a pattern that recurred through this work. Four times a
"model defect" turned out to be an error in how the data was being compared,
and a fifth — an entrainment shortfall that appeared to grow with distance —
was pre-registered, tested and falsified. Aggregate statistics on a model with
more than one fault are very good at producing plausible wrong diagnoses.


## The fit basis, written down

The published liquid hydrogen sub-model figures rest on **42 vertical fits**,
and which 42 was never recorded as a criterion — only as a number. Two
readings were possible and they give different answers, so it is written here.

**Every well-constrained fit on a horizontal trial, with no momentum filter.**
Well-constrained means `r² > 0.85`, at least four points, and `sigma_z < 3 m`.
For the far field the same convention gives 18 well-constrained arc fits of
65.

The momentum filter — exit speed at least ten times the wind — applies to the
**concentration** statistics only, never to the fits. The reason is not
convenience: a fitted plume centre is a geometric measurement and does not
care why the plume is where it is, while the filter is a statement about where
a steady jet model applies. Filtering the fits as well gives 11 and a
different spread ratio.

The spread ratio is a different population again: **23 fits, momentum-filtered
after all**, and the statistic is the **mean of per-fit ratios**, not the ratio
of medians — on this data those differ by 0.11, which is larger than the
correction being measured on some subsets.

Two populations and two statistic forms, in adjacent rows of the same table,
distinguished nowhere. That is the whole lesson.

**State the aggregation with every statistic**: distance band, release
heights, filter, statistic form, n. Every discrepancy found in reconstructing this work came
from a number reported without its subset — a rise quoted as one figure that
was three medians over three bands, a spread ratio over 42 fits described
alongside a concentration statistic over 9 trials. None of it was wrong; all
of it was unlabelled.
