# Pre-registration: conservative plug-to-Gaussian establishment

**Frozen before the first result from the candidate.** Do not edit above the
`RESULTS` line after a candidate result is known.

## Defect isolated

The published Winters/HyRAM algebra supplies a smooth established-flow state,
but a direct flux audit found that its first Hecht--Panda condition carries
about 15% less hydrogen species flux after the plug-to-Gaussian conversion
than at the plug source. The downstream ODE conserves the flux it receives;
it cannot restore a boundary-condition loss.

This is separate from the notional-nozzle calculation, whose measured-throat
mass, momentum and total-energy residuals are already below `1e-8`.

## Candidate fixed before running

Keep the published centreline-velocity assumption and development length.
At the Gaussian start, solve only for width, centre density and centre fuel
mass fraction such that the numerically integrated Gaussian profiles exactly
match the plug source's:

1. hydrogen species flux,
2. axial momentum flux, and
3. total enthalpy-plus-kinetic-energy flux relative to ambient.

The total mass flux is allowed to grow because ambient entrainment is the
physical purpose of the development zone. No empirical coefficient or Raman
observation enters the solve. The published algebraic state remains available
for official-oracle comparison.

## Predictions and decision rule

1. All three normalized establishment residuals must be below `1e-8` for all
   nine sources.
2. All four Raman slopes must remain inside the already frozen 25% band.
3. Both median radial coefficients should lie inside the paper's printed case
   ranges, 33--64 for mass and 21--49 for temperature. A rounded boundary
   miss must be reported rather than hidden.
4. The candidate is accepted only if it meets 1--3 without changing
   `lambda=1.16`, `beta_A=0.28`, `alpha=0.082`, the camera range or the
   regression definitions.
5. Existing DEGADIS/Fortran paths must remain unchanged.

---

## RESULTS

Rejected. On the first source the best positive-state solution leaves
normalized residuals of **-7.79% species, +9.61% momentum and -3.47% total
energy**. The three flux constraints have no simultaneous physical solution
when the centreline velocity is fixed to the plug velocity. This candidate is
not used and the acceptance tests were not run on the remaining cases.

The failure identifies the missing degree of freedom: an established
Gaussian profile has already entrained air, so its centreline velocity cannot
also be prescribed independently. A second, separately frozen candidate adds
the total mass acquired over the development length and solves centreline
velocity as the fourth unknown; see
`prereg-entraining-gaussian-establishment.md`.
