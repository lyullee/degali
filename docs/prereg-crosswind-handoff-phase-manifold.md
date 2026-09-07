# Pre-registration addendum: phase-manifold table domain

Date frozen: 2026-09-04

## Trigger

The phase-profile interface passed the measured-throat representative case,
but the pre-registered Spadeadam pilot stopped before concentration scoring.
At `S = 3.0 m`, tests 4/6 left 3.31%/0.575% total-mass and 3.94%/0.622%
hydrogen residuals. In both cases the projection reached the transferred
table's maximum concentration: the instantaneous radial profile ended at its
current centre and could not represent the slightly more concentrated state
required by JETPLU's different finite-section convention.

No Spadeadam downwind prediction was produced or inspected.

## Candidate

Retain the full handoff radial profile up to its centre concentration. Above
that concentration, append the centre states already computed by the same
conserved near-field solution between Gaussian establishment and the handoff.
Each appended node carries the model's own density, phase-equilibrium
temperature and component specific enthalpy. This forms a phase
**manifold** over the actually visited near-field states; it does not
extrapolate a property or add a fitted coefficient.

Only upstream centre nodes whose concentration is strictly above the handoff
centre are appended, avoiding a multivalued overlap with the radial profile.
The ambient row and all lower concentrations remain the handoff radial state.

## Frozen tests

1. Re-run the representative 0.08 m measured-throat interface. It must retain
   all original interface passes.
2. Re-run Spadeadam source reconstructions 4 and 6 at the already fixed 3.0 m
   boundary. All interface criteria must pass before any sensor prediction is
   produced.
3. Only if both pass, score the six arcs with the unchanged metrics and
   decision limits in `prereg-coupled-crosswind-spadeadam.md`.

Failure leaves `phase_profile` as the measured-single-phase research default
and records the reconstructed field sources as outside the coupled model's
representable domain.
