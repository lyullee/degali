# Pre-registration addendum: flux-averaged handoff energy closure

Date frozen: 2026-09-04

## Reason for the second candidate

The centre-state mixing-line candidate in
`prereg-conservative-nearfield-crosswind-handoff.md` was run first and failed
unchanged. At 0.025 m its four native JETPLU balance fluxes closed below
`1e-15`, while hydrogen half-width differed by 5.05% and centre temperature
by 3.97 K. The preliminary numerical energy integral gave 6.05%; the later
exact Gaussian energy integral corrected that diagnostic to 13.46% without
changing the already frozen rejection. No experimental downwind value was
inspected or used in reaching this result.

The failure is consistent with loss of the near-field radial enthalpy profile
when it is collapsed onto JETPLU's one-parameter adiabatic mixing line.

## Candidate

Keep the endpoint composition fixed at the near-field centre composition, but
replace the endpoint *centre temperature* anchor with a virtual source
temperature determined by the target integral total-energy flux. At every
trial virtual temperature:

1. rebuild the mixed-source adiabatic table,
2. project total mass, hydrogen, and both momentum components exactly, and
3. calculate total energy independently over the JETPLU Gaussian section.

A bracketed scalar root then selects the virtual endpoint temperature whose
integral energy equals the near-field energy. This temperature is a sub-grid,
flux-averaged closure; it is not asserted to be the physical centre
temperature. The projected physical centre temperature remains an independent
screen.

The root search is restricted to 40 K below or above the measured near-field
centre temperature, clipped to 14.1 K and the ambient temperature. It may not
extrapolate thermodynamic properties outside that bracket. Failure to bracket
the energy target rejects the candidate.

## Frozen acceptance conditions

The original conditions remain, except that the total-energy residual is now
required to be at most `1e-8`, rather than 2%, because energy is an explicit
constraint in this candidate. The independent physical screens remain:

- hydrogen mass-density half-width mismatch at most 5%,
- projected centre-temperature mismatch at most 2 K,
- native mass/species/momentum residuals at most `1e-8`, and
- finite, positive state inside the rebuilt table.

The source-momentum and Houf--Schefer entrainment terms remain copied from the
near-field model. No entrainment or dispersion coefficient is fitted.

## Interpretation

A pass establishes only a conservative and thermodynamically compatible
interface. It does not validate downwind concentration or ground interaction;
those require independent PRESLHY/FLADIS comparisons after this screen.
