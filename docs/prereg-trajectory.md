# Pre-registration: the trajectory and the vertical spread

> **Later width-normalisation audit:** the vertical-width premise below was
> false. The fit stored an e-folding width as sigma; on JETPLU's
> standard-deviation convention the corrected ratio is 1.033. The trajectory
> findings remain valid. See `gaussian-width-convention.md`.

**Registered before any of the tests below were run.**
Do not edit above the `RESULTS` line once results are known.

---

## Background

Fitting a Gaussian in height at each near-field arc gives the plume's centre
and its vertical spread directly, and against those the model has two separate
faults (`local-validation.md`):

| quantity | measured | modelled |
|---|---|---|
| rise at 5–7 m | −0.12 m | +1.07 m |
| `sigma_z` at 5–7 m | 0.91 m | 0.63 m |
| `sigma_y` at 14 m, meander removed | 0.8–2.2 m | 1.14 m |

The lateral spread is right. The vertical spread is 1.4 times narrow, flat
with distance. The trajectory climbs when the measurement says it does not,
and the gap is about 1.5 `sigma_z`, which removes the sensors from the plume.

Two prior attempts at a single-cause explanation failed: added mass, tested
and falsified elsewhere on wind-tunnel data, and an entrainment shortfall,
pre-registered here and falsified. So this registration treats the two faults
as separate and does not assume a common cause.

## Hypothesis

`JETPLU`'s vertical momentum balance is Boussinesq: buoyancy enters as
`-RK1 g gamma ccsysz` with `gamma = (rho - rho_a)/cc`, linear in the density
difference. Xiao and co-workers report that for hydrogen "the normalised
trajectory will not collapse when the Froude number is small, which means the
Boussinesq approximation is invalid when the buoyancy effect is comparable
with the momentum effect". A flashing LH₂ jet leaves at 4.4 times ambient
density and is buoyant within a metre, so it is in that regime throughout.

## Predictions

| # | Prediction | Criterion |
|---|---|---|
| **P-T1** | The rise error is not a source-state error | recomputing the mixing line from the flashing source rather than saturated vapour changes the modelled rise at 5–7 m by **< 0.2 m** |
| **P-T2** | The rise error is not the drag term | enabling the `JETPLU` form drag changes the modelled rise by **< 0.2 m** |
| **P-T3** | The rise error scales with the buoyancy term | halving the buoyancy coefficient roughly halves the rise, confirming the term is what drives it and not the initial angle or the entrainment |
| **P-T4** | The measured rise is genuinely near zero, not an artefact of the sensor span | restricting the fits to arcs where the plume centre sits at least 0.2 m inside the top and bottom sensors leaves the median rise **within 0.3 m of zero** |
| **P-T5** | The spread deficit is in the split, not the total | the modelled product `sigma_y sigma_z` at 5–7 m is within **30 %** of `sigma_y(plume) sigma_z(measured)`, while `sigma_z` alone is out by 40 % |
| **NC3** | Negative control | the five EPA cases reproduce to 1e-12 with every option off |

## What would falsify each reading

- **P-T1 or P-T2 fails** — the rise is a source or drag problem, not the
  buoyancy formulation, and the Boussinesq reading is wrong.
- **P-T3 fails** — something other than the buoyancy term is lifting the
  plume, and the diagnosis is misplaced.
- **P-T4 fails** — the measured "no rise" is an artefact of an array that
  cannot see a plume once it leaves the top sensor, and the comparison is
  invalid rather than the model.
- **P-T5 fails** — the vertical spread deficit is a real shortfall in total
  entrainment after all, not a mis-split, and `sigma_y` agreeing was luck.

## What will not be done

No coefficient will be fitted to the PRESLHY data. If the Boussinesq reading
survives, the correction to try is the published non-Boussinesq form, with its
own coefficients, and it will be registered separately before being run.

---

## RESULTS

| test | result | verdict |
|---|---|---|
| **P-T1** mixing line from the flashing source | rise 0.81 → 0.60 m, Δ **−0.22 m** | **marginal fail** (threshold 0.20) |
| **P-T2** form drag enabled | Δ **0.00 m** | pass |
| **P-T3** buoyancy × 0.5 | rise 0.81 → 0.46 m | pass |
| **P-T3** buoyancy × 0.0 | rise 0.81 → **0.00 m** | pass |
| **P-T4** fits with the centre inside the sensor span | median rise **+0.04 m**, n = 15 | pass |
| **P-T5** spread deficit in the split | **not evaluable** — see below | — |

### The rise is the buoyancy term, and nothing else

Halving the buoyancy coefficient halves the rise; zeroing it removes the rise
entirely. Form drag makes no difference at all. So the climb is not the
initial angle, not the entrainment and not the drag: it is the buoyancy
balance, exactly as read.

**P-T1 fails marginally.** Rebuilding the mixing line from the flashing
source rather than saturated vapour lowers the rise by 0.22 m against a 0.20 m
threshold. That is a real contribution and it is in the right direction — the
flashing line is colder and denser at a given concentration — but it accounts
for a quarter of the 0.9 m discrepancy, not for it. The threshold was set at
0.2 m in advance and 0.22 exceeds it, so this is recorded as a fail rather
than argued down.

**P-T4 passes, and it mattered.** The near-field sensors span only 1.25 to
2.00 m above ground on the axis, so a plume that climbed past the top sensor
would be invisible and could masquerade as one that never climbed. Restricting
the fits to those whose centre sits at least 0.2 m inside that span leaves 15
of 18, with a median rise of +0.04 m. The measured plume genuinely stays put.

### P-T5 could not be evaluated as written

It compares the modelled product `sigma_y sigma_z` at 5–7 m against the
measured one. But `sigma_z` is measured in the near field, at 0.35 to 6 m,
while `sigma_y` is measured in the far field, at 10 and 14 m: the near-field
array carries only `y = 0` and `±1 m`, too sparse to fit a lateral profile.
Comparing the model's 5–6 m lateral spread against a measurement made at 14 m
tests nothing.

The prediction was written without checking that the two measurements exist at
the same distance. It is recorded as not evaluable rather than answered with
the mismatched comparison, and the question — whether the vertical deficit is
a mis-split or a shortfall in the total — remains open.

### Where this leaves the diagnosis

The trajectory fault is isolated to the buoyancy balance, with a quarter of it
attributable to the source state. That is consistent with the Boussinesq
reading but does not establish it: a Boussinesq term and a correct term differ
in *how* they scale with density ratio, and testing that needs releases
spanning a range of ratios rather than one campaign at 4.4.

**No coefficient was changed.**
