# Pre-registration: starting the jet from the expanded source

> **Later width-normalisation audit:** width ratios in this historical test
> used a measured e-folding width under a sigma label. Multiply them by
> `sqrt(2)` for comparison with JETPLU's standard deviation; see
> `gaussian-width-convention.md`.

**Registered before the change was adopted.**
Do not edit above the `RESULTS` line once results are known.

---

## Background

Five candidates for the trajectory fault were registered and tested
(`prereg-trajectory.md`, `prereg-boussinesq.md`), and none closed it. All five
were run on the same starting state: the orifice-plane two-phase density with
a contaminant mass fraction of one, through the orifice area.

That starting state is wrong, and it is the mistake EPA's 1991 evaluation
records the SLAB developer objecting to:

> for jet releases the source area should be the cross-section of the fully
> expanded jet rather than the orifice. This change would substantially
> reduce the source velocity and could change the model results considerably.

This work took the density from after the flash but the concentration and the
area from the orifice. Sweeping the starting state along the expansion, with
`rho u^2 A` conserved so the momentum is carried correctly, moves all three
measured quantities monotonically:

| start | `sigma_z` ratio | rise at 5–7 m | concentration MG | start velocity |
|---|---|---|---|---|
| orifice | 0.64 | 1.07 m | 0.737 | 70 m/s |
| half way | 0.65 | 0.84 m | 0.820 | 81 m/s |
| **expanded** | **0.70** | **0.60 m** | **1.000** | 101 m/s |
| measured | 1.00 | **−0.12 m** | 1.00 | |

## What is being tested

Whether to adopt the expanded source as the jet starting state. The three
quantities are judged **separately**, because a concentration comparison is
where the faults cancel — that cancellation is what made the first LH₂ result
look better than it was, and MG reaching exactly 1.000 here is more likely to
be two errors meeting than one error gone.

## Predictions

| # | Prediction | Criterion |
|---|---|---|
| **P-X1** | The rise improves and does not close | rise falls by **more than a third** from 1.07 m and remains **above 0.3 m** — better, not fixed |
| **P-X2** | The vertical spread barely responds | `sigma_z` ratio moves by **less than 0.15**, confirming it is a separate fault rather than the same one |
| **P-X3** | The start velocity stays physical | between 50 and 200 m/s for the 25.4 mm releases, not the 457 m/s that squeezing the expanded source through the orifice gives |
| **P-X4** | Concentration MG near 1 is not evidence | with the rise still wrong by 0.7 m and `sigma_z` by 30 %, MG within 0.05 of unity must be **reported as coincidence**, not as agreement |
| **NC6** | The change touches nothing else | the five EPA cases reproduce to 1e-12; the REDIPHEM statistics are unchanged, since neither uses this path |

## What each outcome means

- **P-X1 and P-X2 both confirmed** — the starting state was one of the five
  candidates and the largest single contributor found so far, and the vertical
  spread is genuinely a second fault.
- **P-X1 refuted upward** (rise closes to under 0.3 m) — then the trajectory
  fault was *entirely* the starting state, and the earlier conclusion that it
  is unfixable within the formulation was wrong by more than it already is.
- **P-X2 refuted** (`sigma_z` moves by more than 0.15) — then the two faults
  share a cause and should be treated together.

## What will not be done

The starting state is fixed by the flash calculation and by conserving
`rho u^2 A`; nothing in it is fitted. The expanded source will not be pushed
through the orifice area to improve the statistics, whatever they do:
that gives 457 m/s for a 25.4 mm release, and it is the error the 1991
evaluation records.

---

## RESULTS

| test | result | verdict |
|---|---|---|
| **P-X1** rise falls by more than a third, stays above 0.3 m | 1.07 → 0.60 m | confirmed |
| **P-X2** `sigma_z` ratio moves less than 0.15 | 0.64 → 0.70 | confirmed |
| **P-X3** start velocity 50–200 m/s | **24–691 m/s** | **refuted** |
| **P-X4** MG near 1 is coincidence | MG 1.000 with the rise still 0.72 m out | held |

**P-X3 failed because the criterion was written from one nozzle.** The 25.4 mm
releases start at 139 m/s, inside the band. The 6 mm at 5 barg starts at
691 m/s, and the 12 mm at 677 — Mach 2.7 and 1.9 against the speed of sound in
saturated hydrogen vapour. Those are under-expanded releases and a subsonic
plume model does not describe their first metre.

Capping the start at the sonic plane by continuity changes nothing measurable
(`sigma_z` ratio 0.70 → 0.69, rise unchanged): the initial velocity is
forgotten within a metre, and what survives is the cross-section, which
continuity preserves either way.

## What followed

Two further corrections were found by the same route — reading what the
published relation actually says rather than what the model assumes.

**Density-scaled entrainment.** Ricou and Spalding's measurement is that the
entrained mass flux scales as `sqrt(rho_ambient / rho_jet)`, restated by Panda
and Hecht for cryogenic hydrogen as "scaling inversely by the square root of
the density of the jet at the nozzle". `JETPLU` has no such scaling. Applied
with the sign inverted it moved `sigma_z` the wrong way, 0.70 → 0.65; applied
correctly, 0.70 → 0.76.

**The plume entrainment coefficient.** This preregistered run used
`alfa1 = 0.0833`, a historical Fischer et al. value whose attribution through
secondary sources is now superseded. Papanicolaou and List (1988) actually
measured 0.0875 for a pure plume and 0.0545 for a pure jet; `JETPLU` ships
0.057, at the jet end. The default has therefore moved to 0.0875 while the
historical dump-reproduction tests request 0.0833 explicitly. With the
measured coefficient the corrected `sigma_z` ratio is about **0.995**.

### Where that leaves the model

**Superseded — recomputed from the reduction; the figures below were never computed by anything.** Recomputing gives MG 0.855, CI [0.672, 1.064], VG 1.30, FAC2 0.85 as shipped and MG 1.070 corrected, over 66 arcs, and the interval now includes 1. See `docs/lh2-recomputed.md`.
| | as shipped | corrected | measured |
|---|---|---|---|
| `sigma_z` ratio | 0.64 | **0.94** | 1.00 |
| rise, 0.5 m releases | — | **+0.19 m** | +0.00 m |
| rise, 1.5 m releases | — | +0.51 m | +0.03 m |
| concentration MG | 0.737 | 1.228 | 1.00 |

The vertical spread deficit is closed to six per cent and the trajectory error
on the low releases to 0.19 m.

### The array can see a risen plume, and does not

Worth checking, because the near-field sensors span only 0 to 1.0 m above
ground for a 0.5 m release. If the plume rose past the top sensor a Gaussian
fit would place its centre near that sensor and "no rise" would be an artefact
of the array rather than a measurement.

With the corrections applied the modelled centre for those releases sits at
**0.69 m, inside the array**. A plume there would have been fitted. It was not:
the measured rise is +0.00 m against a modelled +0.19.

So the residual trajectory error is real and not a censoring artefact. An
earlier version of this work asserted that from fits that did not distinguish
the two, which was luck rather than evidence.

### What did not work

| correction | source | effect on the rise |
|---|---|---|
| additional pressure drag | Mack et al., EFFECTS | 0.02 m at their `C_d = 0.39` |
| vertical component in the shear velocity | Mack et al., EFFECTS | under 5 % |
| sonic cap on the start | — | none |

Both EFFECTS corrections are aimed at plumes rising fast enough for a
quadratic drag and a vector shear velocity to matter — Witcofski's pool rises
twenty metres. These jets rise 0.2 to 0.5 m over six, so the same terms are
negligible. **The corrections are right for their regime and do not transfer
to this one**, which is worth recording because the papers do not say so.
