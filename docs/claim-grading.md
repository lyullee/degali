# Grading the claims

Not every number here is supported the same way, and a paper that presents
them as though they were invites a reviewer to discount the strong ones along
with the weak. The grades below separate them by *kind of evidence* rather
than by topic.

| grade | meaning | how to write it |
|---|---|---|
| **A** | deterministic; re-running gives the same number, and sample size is irrelevant | "is", "reproduces" — assert |
| **B** | campaign statistic with sampling uncertainty quantified | "estimates", "supports"; state the interval and n |
| **C** | direction only; the interval includes the null | "suggests", "is consistent with" — **do not assert** |

## Grade A — deterministic

| claim | value | basis |
|---|---|---|
| Reproduces DEGADIS 2.1 | 1e-12 across all six programs | probes against the Fortran, enforced by the suite |
| The five EPA cases reproduce | exactly | golden `.LIS` comparison |
| `ADIABAT` assigns `wa` on the branches for `ifl` 0, −1, −2, 2 and not for 1 | — | `TPROP.for`; `SZF.for` 95 passes `walay`, which appears once in the file |
| `PSS` and `SSG` write the layer temperature to different variables | `temlay` vs `temlam` | `PSS.for` 89, `SSG.for` 92; every other argument matches |
| `GAMINC` is unregularised **by design** | factor Γ(1/(1+α)) | `INCGAMMA.for`: the comment above the line says so — a porting trap, not a defect |
| DEGADIS has no vertical momentum at ground level | — | `PSS`/`SSG` carry gravity only in the common block |
| `JETPLU` *does* carry buoyancy | `-RK1*gg*gamma*ccsysz` | source; EX1 never touches down |
| LH₂ is buoyant across the whole flammable range | rho/rho_a is 0.98 at LFL, 0.87 at stoichiometric, 0.69 at UFL | CoolProp mixing line, no fitting |
| The adopted Spadeadam test-6 source changes from dense to buoyant near 3.2 m | — | conserved flash state and CoolProp mixing table, no fitted coefficient |
| The mixing table outgrows `IGEN = 42` | needs 46 nodes | adaptive thinning, deterministic |
| Hydrogen's critical temperature is below ambient | 33.1 K | property data; no saturation line above it |
| Two readings of Hall & Walker agree to 30 % | Ri* 7 vs 10 | independent reductions of the same experiments |
| EPA's published DEGADIS bias on Burro is an evaluation-height artefact | FB −1.29 at ground vs −0.29 at 1 m, against −1.07 published | reproduced with this port |
| The PRESLHY far field cannot locate the plume | wind resolution 22.5°, stand spacing 10–12° | measurement programme |
| The near field cannot test transient behaviour | travel time < 1 s, sampling 0.3 s | geometry |

## Grade B — statistically supported

| claim | value | n |
|---|---|---|
| LH₂ near-field concentration, fully conserved source | MG 1.047, CI [0.759, 1.402], VG 1.425, FAC2 0.84 | 62 arc maxima at 0.79–6 m, 9 trials |
| The bias is the same at both release heights | 0.722 vs 0.746 | 23 and 46 |
| Wind-steered releases are a different population | VG 3.5 vs 1.4; GSD 2.94 vs 1.16 | 13 vs 7 trials |
| Burro reads high at the lowest height | MG 0.811, CI [0.63, 1.02] | 61 |
| DEGADIS's vertical profile is too steep | MG 7.3 at 3 m, 5500 at 8 m | 59, 53 |
| The lateral spread is too wide | 1.2–2.9× | 6 arcs, Desert Tortoise |
| The vertical defect appears on an independent dataset | MG 0.71 at 0.1 m falling to 0.11 at 8.5 m | 76 SMEDIS sensors, ammonia jets |
| Corrected LH₂ plume-centre bias is small on the filtered geometry set | mean/median signed error +0.036/+0.002 m; MAE 0.145 m | 23 vertical fits |
| Corrected vertical spread agrees on average | model/measured Gaussian standard deviation = **1.033** (median 1.042) | same 23 fits |
| The modelled lateral spread is correct once meander is removed | 1.14 m against 0.8–2.2 m | 18 arc fits |

Passing a screening bound is not proof of zero bias. In particular, the
corrected LH₂ MG interval includes unity and also extends beyond the nominal
upper screening bound; quote the interval with the point estimate.

## Falsified, and recorded so

| hypothesis | how it died |
|---|---|
| FLADIS channel numbering could be guessed | the guessed channels read 302 and 23.5 at 20 m and are not concentrations |
| Coyote's bias is the pool radius | forcing smaller diameters makes it worse |
| LH₂ rainout feeds a ground-level source | the report: no rainout during un-impinged elevated releases |
| the notional nozzle fixes the near field | FAC2 0.69 → 0.48; a wider source is slower and dilutes less |
| the section is too narrow by a factor of two | measured on mixed release heights; the real factor is 1.4 and widening made it worse |
| corrected `sigma_z` remains 27% too narrow | measurement fit stored the e-folding width `w` as sigma; after `sigma=w/sqrt(2)`, the ratio is 1.033 |
| an all-gas 68 K handoff fixes the unmodelled air-condensation zone | phase state becomes valid and FAC2 improves, but PRESLHY VG worsens 1.611 → 2.027 and centre-height MAE 0.160 → 0.274 m on common data |
| the residual is an entrainment shortfall growing with distance | pre-registered; slope +0.018, CI [−0.168, +0.171] |
| the trajectory fault is the Boussinesq approximation | pre-registered; the correction increases the rise 0.81 → 1.08 m |
| added mass suppresses rise for wide sources | tested elsewhere on wind-tunnel data across six widths; no coefficient fits |

## Grade C — direction only

| claim | why it is only C |
|---|---|
| Lift-off height, RMS 2.9 m | n = 4, and three of the four measurements are inversions of thermocouple data through a mixing model |
| The buoyancy regime is right | 4 of 4 is 4 of 4; the categorical agreement is strong but the sample is tiny |
| The residual is an entrainment shortfall | **falsified.** Pre-registered: the distance slope is +0.018, CI [−0.168, +0.171]. The apparent growth was medians over a seven-trial subset and does not survive the full 69 pairs. Source rate, stability and the entrainment coefficient are all ruled out; a distance-independent factor of 1.35 remains, unexplained |

## A note on intervals

Bootstrap intervals here resample **trials**, not points. Readings within one
trial share a release rate, a wind and a source estimate, so treating them as
independent narrows the interval by roughly the square root of the number of
sensors and overstates the evidence. On the LH₂ set, point-level resampling
gives [0.68, 0.80] where trial-level gives [0.59, 0.92].

## Audit findings

Two error-concealment patterns were found by static inspection rather than by
a failing test, which is why they are recorded here.

**Silent equation-of-state fallbacks.** Three lookups in `CoolPropBackend`
fell back to a 1989 correlation or to `None` when CoolProp raised, and said
nothing — so a run could use the legacy correlation for one property over one
temperature range while reporting itself as `coolprop`. Behaviour unchanged;
they are now counted and `fallback_report()` says what happened. Nothing in
the published results used a fallback.

**A two-phase gap filled with the nearest valid value.** Found earlier in the
same class: the interpolation grid was filling the liquid-vapour gap by
carrying the last good number forward, which is a plausible-looking answer
with no physics in it.

Neither changed a published number. Both are the kind of defect that hides a
wrong answer rather than producing an obviously wrong one, and the second was
only found because the first prompted a look.

## What this is for

The LH₂ concentration result (B) and the structural results about DEGADIS (A)
do not depend on the lift-off sample size. Writing them as though they did —
or writing the lift-off result as though it were as firm as they are — is what
the grading is meant to prevent.
