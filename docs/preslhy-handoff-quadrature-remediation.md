# PRESLHY handoff energy-quadrature remediation

Date: 2026-09-05, after the first pre-registered PRESLHY coupled run and before
its confirmatory rerun.

## Observed numerical failure

The initial run accepted trials 12 and 23 but rejected trials 10, 11, 22, 24
and 25.  All seven closed total mass, hydrogen mass and vector momentum to
about machine precision.  All seven were also within the pre-existing 2%
energy, 5% H2 half-width and 2 K centre-temperature screens.  The five
rejections came solely from the difference between 32- and 64-point Gaussian
energy quadrature, which ranged from `1.51e-5` to `4.90e-5` against the fixed
`1e-5` convergence limit.

This is a resolution failure in the diagnostic integral, not evidence for a
new physical term and not permission to relax an acceptance threshold.

## Fixed remediation

For phase-profile and phase-manifold handoffs, start at the caller's requested
quadrature order and repeatedly double it.  Compare each new energy flux to
the preceding order and stop at the first relative change at or below
`1e-5`.  Permit at most four doublings (the default sequence is
32/64/128/256/512).  If 512 points still fail, the interface remains rejected.

The final, finer energy flux is used for the unchanged 2% physical-energy
screen.  The actual final point count and last relative change are reported.
Nothing else in the source, near-field ODE, thermodynamic manifold, projection
or downstream coefficients changes.

## Decision rule

Rerun the original seven-trial protocol exactly.  The coupled candidate can
only proceed to the already frozen field-performance decision if all seven
interfaces meet the same original screens.  Subset statistics from the first
run remain diagnostic and cannot be used for promotion.
