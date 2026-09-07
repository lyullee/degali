# Pre-registration: a non-Boussinesq buoyancy term

**Registered before the correction was implemented or run.**
Do not edit above the `RESULTS` line once results are known.

---

## Background

The trajectory fault is isolated (`prereg-trajectory.md`): halving the
buoyancy coefficient halves the modelled rise and zeroing it removes the rise
entirely, while form drag makes no difference. The measured plume does not
rise at all over six metres — median +0.04 m across 15 fits whose centre sits
inside the sensor span — where the model lifts it.

The magnitude quoted here when this was written, 1.07 m, is superseded: with
the configuration validated against a full parameter dump the model lifts the
0.5 m releases 0.24 m beyond 3 m and about 0.9 m beyond 5 m. The *fault* is
unchanged and so is everything this pre-registration concluded from it —
buoyancy carries the rise, form drag does not — but the number was a median
over an unstated subset. See `lh2-recomputed.md`.

`JETPLU` writes the vertical momentum source as

    -RK1 * g * gamma * cc * sy * sz,    gamma = (rho - rho_a) / cc

which is linear in the density difference: a Boussinesq form. Xiao and
co-workers, building a non-Boussinesq integral model for hydrogen safety,
report that "the normalised trajectory will not collapse when the Froude
number is small, which means the Boussinesq approximation is invalid when the
buoyancy effect is comparable with the momentum effect".

A flashing LH₂ jet leaves the orifice at 4.4 times ambient density and is
lighter than air within a metre, so it spans the regime where that applies.

## The change to be tested

Boussinesq buoyancy accelerates a parcel at `g (rho_a - rho)/rho_a`. The
non-Boussinesq form divides by the *parcel's own* density instead:

    a = g (rho_a - rho) / rho

The two agree when `rho` is close to `rho_a` and differ by the density ratio
when it is not. For a plume lighter than ambient — `rho < rho_a` — the
non-Boussinesq form gives a *larger* acceleration, not a smaller one.

**That is the wrong direction for this fault**, and saying so before running is
the point of registering it. The prediction below is therefore that the
correction fails, and the test is whether it fails for the reason expected.

## Predictions

| # | Prediction | Criterion |
|---|---|---|
| **P-B1** | The non-Boussinesq form increases the modelled rise | rise at 5–7 m goes **up**, not down |
| **P-B2** | The increase is of order the density ratio at the point where buoyancy is strongest | the change is between 1.1× and 3× the Boussinesq rise |
| **P-B3** | It does not fix the concentration bias either | MG moves **away** from 1, or by less than 0.05 |
| **NC4** | With the correction disabled, every result is bit-identical | the five EPA cases and the LH₂ statistics reproduce exactly |

## What each outcome means

- **P-B1 confirmed** — the Boussinesq reading of the fault was wrong. The term
  is not the issue; something else is holding the real plume down, and the
  candidates are the droplet phase (mass without volume, already implemented
  and too weak) and ground effect (a jet released at 0.5 m with `Sz` growing
  past that is in contact with the surface).
- **P-B1 refuted** — the non-Boussinesq form *reduces* the rise, contrary to
  the algebra above, which would mean the implementation is wrong rather than
  the physics.

Either way the correction is not adopted on this evidence: it would need to
both reduce the rise and be published for this regime, and the algebra says it
cannot do the first.

## What will not be done

No coefficient will be fitted. If P-B1 confirms, the registration closes with
the Boussinesq reading falsified and the ground-effect hypothesis registered
separately.

---

## RESULTS

All three predictions confirmed, which means **the correction fails**, as
registered.

| | rise at 5–6 m | MG | VG | FAC2 |
|---|---|---|---|---|
| Boussinesq (as shipped) | 0.81 m | 0.738 | 1.41 | 0.83 |
| non-Boussinesq | **1.08 m** | 0.754 | 1.39 | 0.84 |
| measured | **−0.12 m** | | | |

| test | result | verdict |
|---|---|---|
| P-B1 rise increases | 0.81 → 1.08 m | confirmed |
| P-B2 increase between 1.1× and 3× | ratio 1.34 | confirmed |
| P-B3 MG not improved by more than 0.05 | 0.738 → 0.754 | confirmed |

**The Boussinesq reading of the trajectory fault is falsified.** The
approximation is genuinely present in `JETPLU` and Xiao and co-workers are
right that it is invalid at this density ratio, but correcting it moves the
plume the wrong way: dividing by the parcel density rather than the ambient
makes a light plume rise *faster*, and the model already rises too fast.

That was predictable from the algebra and was predicted here before the run.
Registering it was worth the effort precisely because the correction is the
one the literature points at, and running it without a stated expectation
would have invited reading a 1.34× worse answer as "the term needed tuning".

The option is kept, defaulted off. It is the more correct form for a release
*denser* than air, where it increases the sink rate, and this campaign has
none of those to test it on.

### What this leaves

Something holds the real plume down that the model does not have. Two
candidates remain, and neither is buoyancy:

**The droplet phase.** A flashing LH₂ jet leaves as a two-phase fluid; at
5 barg only a quarter has flashed. Droplets carry mass without volume and fall
relative to the gas. `JetPlume.liquid_fraction` implements the mass part and
moves the rise by about 5 per cent — right direction, an order too small. The
falling part, the "slip" that Giannissi and co-workers describe, is not
implemented and cannot be within a single-velocity formulation.

**Ground effect.** The releases are 0.5 or 1.5 m up with `Sz` reaching 0.9 m
by six metres, so the plume is in contact with the surface over most of the
measured range. `JETPLU` models a free plume until touchdown: it entrains
through its whole perimeter and rises unrestrained. `ground_effect` exists and
changed the near-field answer by 0.01 in MG, but it was tested before the
trajectory was measured and against concentrations rather than against the
rise. It deserves a test against the quantity it should affect.

### Both remaining candidates were then tested against the rise

They had been tried before, but against *concentration*, which three faults
contribute to. Against the trajectory itself:

| variant | rise at 5–6 m | MG |
|---|---|---|
| baseline | 0.81 m | 0.738 |
| ground effect | **0.81 m** | 0.737 |
| droplet mass | 0.77 m | 0.744 |
| both | 0.77 m | 0.742 |
| **measured** | **−0.12 m** | |

Ground effect moves the rise by 0.00 m and droplet mass by 0.04 m. The gap is
0.93 m.

### The trajectory fault is unexplained, and probably not fixable here

Five candidates, all tested, all insufficient:

| candidate | effect on the rise |
|---|---|
| source state (flashing mixing line) | −0.22 m — a quarter of the gap |
| form drag | 0.00 m |
| non-Boussinesq buoyancy | **+0.27 m — the wrong way** |
| ground effect | 0.00 m |
| droplet mass | −0.04 m |

What remains is the slip Giannissi and co-workers describe: liquid droplets
falling relative to the gas, leaving a vapour-enriched cloud behind and
carrying momentum downward. `liquid_fraction` implements the *mass* those
droplets add, which is why it helps a little. It cannot implement their
falling, because an integral model of this kind carries one velocity for the
whole cross-section.

So the honest position is that this trajectory error is unlikely to be
fixable within the formulation. Closing it needs a second velocity field for
the dispersed phase, and that is a different class of model rather than a
coefficient or a term.

That is a bounded, stated limitation rather than an open question, which is
what the sequence of registrations was for.
