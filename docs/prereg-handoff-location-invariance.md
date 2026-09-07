# Pre-registration: handoff-location invariance

Date frozen: 2026-09-04

## Question

Does the coupled prediction depend materially on choosing 0.08, 0.10 or
0.13 m as the boundary between two models, even when every individual
interface passes conservation and thermodynamic screens?

## Fixed run

Use the same representative 5 bar-class measured-throat case and atmosphere
as the accepted interface test. Compute one near-field solution to 0.13 m at
the 41-point pilot resolution, then create independent phase-profile JETPLU
models at `S = 0.08`, `0.10`, and `0.13 m`. Do not reuse a thermodynamic table
between handoffs.

Run each accepted state to at least 2 m. At common downwind coordinates 0.5,
1.0 and 2.0 m compare centreline mole fraction, centre height, lateral sigma
and vertical sigma.

## Frozen screen

All three interfaces must pass. Relative to the three-boundary mean at each
coordinate:

- centreline mole-fraction range at most 5%,
- lateral- and vertical-sigma range at most 5%,
- centre-height absolute range at most 0.02 m.

Report the maximum range for each quantity. Failure makes handoff position an
explicit uncertainty input; it must not be hidden by selecting the best
station.

## Phase-profile result

The 0.08 and 0.10 m interfaces passed, but the 0.13 m projection reached the
instantaneous profile table's centre-concentration ceiling and left
mass/species/momentum residuals of 0.160%/0.169%/0.0184%. The three-boundary
screen therefore fails before downstream invariance is scored.

The phase-manifold domain extension had already been pre-registered in
`prereg-crosswind-handoff-phase-manifold.md` before this run. Apply that
unchanged table construction to all three positions as the only follow-up;
the downstream comparison and limits above remain unchanged.

## Phase-manifold result

All three phase-manifold interfaces pass. The largest downstream boundary
ranges over 0.5, 1.0 and 2.0 m are:

| quantity | maximum range | limit | result |
|---|---:|---:|---|
| centreline mass concentration | 0.305% of mean | 5% | pass |
| centreline mole fraction | 0.198% of mean | 5% | pass |
| lateral sigma | 0.082% of mean | 5% | pass |
| vertical sigma | 0.347% of mean | 5% | pass |
| centre height | 0.000307 m | 0.02 m | pass |

The interface diagnostics also remain inside their original limits at every
station. Phase-manifold is therefore position-invariant over the tested
interval. Per its separate pre-registration, it remains an explicit extended
domain rather than replacing the measured-single-phase `phase_profile`
default, because the reconstructed Spadeadam test-4 source failed its physical
interface/applicability screens.

The public coupled runner defaults to the strictly confirmed 0.08 m
near-field endpoint. A caller requesting a longer near-field domain must
select and audit its handoff explicitly rather than silently inheriting the
final numerical station.
