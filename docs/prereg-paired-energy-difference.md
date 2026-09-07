# Common-grid actual-energy difference diagnosis

2026-09-06. After the phase-consistent witness has stable rates and passed
retained volume equations, but its separate-integral finite difference fails.
The old result remains failed under its original complete screen.

In the earlier fixed candidate, the total energy derivative is only about
5.38 while the enthalpy/modal fluxes are large. Only the global energy row
dominates the separate-integral FD failure: relative4.81e-4 at h and9.62e-5
at h/2, versus<=1.7e-7 for all other rows. Cancellation of nearly equal
large integrals is a plausible numerical cause, not yet assumed proven.

Use NEW verification-only paired_moment_difference on the fine
phase-consistent witness. Independently perturb the actual encoded state
and scalar shapes by±h*r as before. Union BOTH states' angular and radial
phase cuts, evaluate their actual EOS fields at common points, form product
differences locally, then accurately sum with math.fsum. Do NOT use EOS
derivatives to construct the differences or alter the witness. Field values
remain binary64; do not claim arbitrary precision.

The conserved products include A*rho, A*C, A*H and actual velocity/theta.
Factor u²/u³ product differences instead of subtracting huge total energy
integrals. Use batches of8, radial8/angular8. Save all5+2m derivatives,
energy-flux absolute term scale, peak batch nodes and numerical errors.

Evaluate four preselected distances relative to the original h:0.5,1,2,4.
The original h and h/2 MUST both meet1e-5 relative derivative error and
their mutual refinement<=1e-5 to count as a same-step successful diagnostic.
The larger two steps are sensitivity evidence only, not post-hoc replacement
of a failed original pair. Stop if a perturbation leaves the physical/logshape
domain; do not expand the shape trust bound.

If the paired method passes, it supports a numerical cancellation diagnosis,
not physical stress closure/downstream/observed validation. If it fails,
report the row/step dependence and pursue a separate conditioning/error
study. Freeze code/tests/spec/input and refuse overwrite.
