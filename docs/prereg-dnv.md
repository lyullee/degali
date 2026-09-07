# Pre-registration: the DNV corrections

> **Scope notice:** the FLADIS aspect-ratio premise below is independent of
> the later PRESLHY Gaussian-width naming error. No same-arc PRESLHY lateral
> width exists. Current PRESLHY vertical-width results are in
> `gaussian-width-convention.md`.

**Registered before any of the three changes were implemented or run.**
Do not edit above the `RESULTS` line once results are known.

---

## Background

Hart and Harper (DNV, *Hazards 32*, 2022) describe three changes made to the
Unified Dispersion Model for Phast 8.61. Two of them address faults this work
has measured independently, and the third addresses the trajectory.

**The spread floor.** In the UDM the heavy-gas spread rate goes as
`sqrt(max(0, rho_cld - rho_a))`, which is zero once the cloud is buoyant, so
"any entrainment of air into the cloud therefore is constrained to make the
cloud taller rather than wider" and "an unrealistic cloud geometry can have
developed". Their fix is a floor:

    dW/dt = max( dW/dt|near-field , dW/dt|heavy , dW/dt|passive )

DEGADIS does the same thing by a different route. `JETPLU` splits the product
`sigma_y sigma_z` by requiring the plume's excess spread to be isotropic,
`sy^2 - sz^2 = sya^2 - sza^2`. Measured against FLADIS, where both spreads are
resolved on the same arc, the model gives `sigma_y/sigma_z = 3.03` where the
measurement gives 1.66 to 1.91: **too flat**, in the same direction.

**The lift-off dwell.** Their criterion is Briggs',
`Ri* = g[rho_cld - rho_a] H_eff / (rho_a u*^2) < -20`, with the addition that
it must hold continuously for

    t_lo = sqrt( 2L / [g (rho_a - rho_cld)/rho_a] ),  L = sqrt(W_eff H_eff)/2

because "clouds that are only transiently buoyant may lift off" otherwise.
`addons.liftoff` has the threshold but no dwell.

**Crosswind entrainment suppression.** Their extended Morton model zeroes
crosswind entrainment inside the potential core and phases it in over a
suppression length, with a decaying drag:

    L_core = 6.4 D / (1 + 4.6/R),   L_supp = (1 + sqrt(rho_0/rho_a)) L_core
    C_d = 0.39 decaying linearly to zero over L_drag = 3 L_supp

## Predictions

| # | Prediction | Criterion |
|---|---|---|
| **P-D1** | A spread floor widens the modelled cross-section | modelled `sigma_y/sigma_z` at the FLADIS arcs moves from 3.03 **towards** the measured 1.66–1.91 |
| **P-D2** | It does not fix the total | the modelled product `sigma_y sigma_z` stays below half the measured value, because a floor redistributes rather than adds |
| **P-D3** | The lift-off dwell changes nothing on the NASA spills | those clouds are buoyant continuously, so a dwell requirement cannot reclassify them; the regime stays 4 of 4 |
| **P-D4** | The crosswind suppression does **not** help the LH₂ jets | DNV state the extended model is "deliberately formulated to not significantly change predictions where `R >> 20`", and these releases have `R` of 10 to 660; MG moves by less than 0.05 |
| **NC5** | With all three off, everything is bit-identical | the five EPA cases and the LH₂ statistics reproduce exactly |

## What each outcome means

- **P-D1 confirmed and P-D2 confirmed** — the isotropy split is a real fault
  and the floor is the right shape of fix, but the vertical deficit is
  *also* a shortfall in total entrainment. Two faults, not one.
- **P-D1 refuted** — DEGADIS's split does not respond to a floor the way the
  UDM's spread rate does, and the analogy between the two models is wrong.
- **P-D4 refuted** — the suppression helps despite DNV's statement about
  `R >> 20`, which would be worth reporting to them as much as using.

## What will not be done

No coefficient will be fitted. `C_d = 0.39`, the `L_core` and `L_supp`
relations and the Briggs threshold of −20 are taken as published. If a
prediction fails, it is recorded; the published values are not adjusted to
rescue it.

---

## RESULTS

*(appended after running; nothing above this line changes)*
