# Pre-registration: the momentum filter is a ratio and should be a Richardson number

> **Primary-source audit, 2026-09-03 — replacement not adopted.** Schefer et
> al. (2008) characterise gaseous hydrogen jets with the densimetric Froude
> number `Fr = u/sqrt(g d |rho_a-rho_j|/rho_j)`. The Richardson number proposed
> below instead divides the density contrast by `rho_a`; it is therefore not
> simply `1/Fr²`, and its `Ri < 1e-4` threshold cannot be borrowed from the
> Schefer regimes. More importantly, the paper's Fr=58–268 experiments are
> small vertical gaseous jets in a quiescent ambient, whereas a flashing LH2
> source changes density sign as it warms and is then bent by crosswind. The
> source plane and density must be pre-specified before a Froude filter can be
> scored. The proposal below is retained as research history, not applied to
> any reported statistic. See `docs/source-froude-applicability.md`.

**Written before testing the replacement.** The defect below was found while
checking a variance statistic, so the ordering matters and is stated plainly:
the failing trial was identified first, and the criterion is being questioned
second. That is the wrong order and it is why this is pre-registered rather
than applied.

## The defect, which is visible without the outcome

The applicability filter keeps a trial when the exit velocity is at least ten
times the wind. It is a ratio with the wind in the denominator, so **a small
wind makes a weak jet look momentum-dominated**, which is the opposite of what
the filter is for.

| trial | wind (m/s) | exit (m/s) | ratio | kept? |
|---|---|---|---|---|
| 11 | 2.70 | 421.7 | 156.2 | yes |
| 24 | 1.90 | 343.7 | 180.9 | yes |
| 10 | 2.47 | 69.9 | 28.3 | yes |
| **20** | **0.57** | **9.8** | **17.3** | **yes** |
| 21 | 0.90 | 16.4 | 18.2 | yes |

Trial 20 has the **lowest exit velocity in the campaign**, a factor of forty
below trial 11, and is kept on the same criterion. Nothing about that argument
requires knowing how trial 20 scores.

## What it costs

Trial 20 carries the entire far-arc variance of the near-field statistic.

| 3–7 m | n | MG | VG | FAC2 |
|---|---|---|---|---|
| all trials | 15 | 2.171 | **50.72** | 0.67 |
| without trial 20 | 13 | 1.107 | **1.35** | 0.77 |

Both are reported. Neither is "the" answer until the criterion is settled.

## The proposed replacement

A source Richardson number, buoyancy over momentum at the orifice:

    Ri = g d (rho_a - rho_j) / (rho_a u_exit^2)

which does not divide by the wind:

| trial | 11 | 24 | 12 | 23 | 10 | 22 | 25 | 21 | **20** |
|---|---|---|---|---|---|---|---|---|---|
| Ri | 6e-8 | 4e-8 | 2e-8 | 3e-7 | 5e-6 | 6e-6 | 6e-6 | 2e-5 | **1e-4** |

Trial 20 is two to four orders above the strongly momentum-driven trials and a
factor of five above the next weakest. **The separation exists; where to cut it
does not follow from the data and must be fixed in advance.**

## Predictions, fixed now

1. A threshold at `Ri < 1e-4` excludes trial 20 and keeps the other eight.
   Under it the 3–7 m band gives MG ≈ 1.11, VG ≈ 1.35, FAC2 ≈ 0.77, and the
   full 66-arc statistic improves.
2. **This is not evidence the model is better.** It is a narrower claim over a
   narrower population, and the paper must say so: the statistic would then
   describe momentum-driven releases in winds above about 1 m/s.
3. The trials the new filter removes should fail *for the stated physical
   reason* — buoyancy-dominated at the source — and not merely score badly. If
   a trial with low `Ri` also scores badly, the filter is not the explanation.
4. If instead the threshold is set anywhere that keeps trial 20, the reported
   variance must carry it. Excluding a trial because it scores badly, under a
   criterion chosen after seeing the score, is not a filter.

## What would falsify the reframing

Trial 20 fails because the model lifts the plume 4.9 m at 6 m downwind against
a top sensor at 1.75 m. If the same lift appears in trials the new filter
*keeps*, then low `Ri` is not what distinguishes them and the buoyant branch is
simply wrong across the board — which is the conclusion the Spadeadam
comparison at 30 m already points to, in a release that is strongly
momentum-driven at the source.

**That second possibility is the more likely one** and it is why this
pre-registration exists rather than a patch.
