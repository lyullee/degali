# DEGALI hybrid-AMR session target

## Objective

Replace the whole-domain uniform-grid transient path with a multifidelity
obstacle path. Use the existing integral transport model away from obstacles,
and activate a locally refined, geometry-resolved patch only where an obstacle
creates stagnation, bypass shear, or a wake. Couple the two paths by conserved
species, momentum, and sensible-enthalpy face fluxes.

## Priority order

1. Local AMR/overset patch around obstacle faces, edges, and the resolved wake.
2. Conservative coarse/patch exchange with a shared event clock and residual audit.
3. Wake-only unsteady closure (URANS/LES-like scalar mixing), not whole-domain LES.
4. ROM/cache for repeated source, wind, and sensor queries after the physical path
   is validated.

## Acceptance gates

- No obstacle: hybrid and integral paths agree within 1% for mass, peak, and arrival.
- Obstacle case: interface species and enthalpy residuals <= 1e-8 relative.
- Refinement: obstacle peak and 4 vol% arrival changes <= 10% between final two
  local patch levels.
- Runtime: 600 s event <= 30 s for the fast path on the reference workstation.
- Accuracy: report peak, dose, arrival, and low-concentration metrics separately;
  never use aggregate MAPE alone as a promotion gate.
- Any missing wake or sensor evidence remains `conditional`/`withheld`.

The existing full-domain 2.5 m result (about 101 minutes) remains a baseline,
not a production target. No speed or accuracy claim is promoted until these
gates are measured on the hybrid execution path.
