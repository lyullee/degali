# Phase-consistent reconstruction of the local flux witness

2026-09-06. Prepared after the old unsplit-angular witness failed the first
independent weak-equation check, before this reconstructed candidate. Its
failure and the memory-stopped audit remain unchanged.

The memory-bounded independent audit supplies current-phase4/8 and8/16
retained matrices/RHS. Recompute original baseline rates and their dependence
on the8 solenoidal amplitudes from EACH matrix. Use the SAME16/96 constraint
mesh, K4 basis,4% edge/production/inflow constraints and minimax amplitude
metric for both candidates. This isolates the volume operator refinement
from changing LP witness sampling. No new physical coefficient or data fit.

Only run after the streaming audit is complete, its dense coarse parity
<=1e-8 and phase matrix/RHS refinement<=1e-5. That audit need not pass its
old frozen witness: this separately registered reconstruction specifically
addresses that candidate's inaccurate volume equation input.

Record the helper's unsplit constraint-grid heat/mass diagnostics without
calling them independent volume checks. Decide volume validity with the
phase-split retained matrices and a fresh streamed fine full-vector heat/mass
audit of the NEW fixed candidate. Cross-evaluate the fine candidate in both
phase matrices, requiring normalized retained residuals<=1e-8. Require
candidate rate/amplitude refinement<=1e-5, actual moment central-difference
tests and their halved-step refinement<=1e-5; independent fixed65 rays/
radial16 must satisfy5% C/H/P, inward mass and nonnegative production/radial
diffusivity. Reapply curvature and fine mass/heat numerical gates1e-8.

Preserve LP infeasibility, numerical/physical failures, nonradial stress and
selector sensitivity. Even a passed phase-consistent local witness is NOT
a physical law determining transverse circulation, a scalar-chi stress
closure, downstream validation or measured accuracy. Do not integrate it
or promote defaults. Freeze all inputs/code/tests/spec and refuse overwrite.
