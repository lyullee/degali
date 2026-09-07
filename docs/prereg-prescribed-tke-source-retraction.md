# Prescribed-Q source energy retraction (boundary-only)

2026-09-06. Follow-up to the finite-TKE operator; no real Q is selected here.

Adding a positive Q flux to an unchanged, already energy-matched initial
section would create energy. For any explicitly supplied Q profile, its flux
must be included in the existing source budget, not added to the target.

At fixed geometry, mean velocity and Q density profile, Q advection is
independent of the C/H shape coefficients. Therefore subtract its actual flux
from the *mean-plus-thermal* energy target, retract the C/H shapes under the
existing six moment constraints, then verify H+meanKE+Q against the original
total target. No energy is added to the EOS; density is recalculated from H/C.
The existing constrained, phase-split retraction can be reused on a detached
projection view. Never modify the original projection or frozen target.

Here the six targets are the enriched projection's [mass, H2, axial momentum,
total energy, buoyancy force, thermal radial second moment]. This is NOT the
different [mass,H2,Px,Pz,energy,force] source array. Direction and centers are
fixed. Buoyancy and the thermal second moment are matching constraints, not
new conserved fluxes. They do not include Q.

The caller supplies the total target explicitly, with Q already counted in
its energy convention. Q0 and its square modes are explicit input, not fitted
to scores, inferred from an assumed turbulence intensity or clipped to a
covariance bound. No dissipation or diffusion constant enters this source
energy calculation. Natural/gradient boundary agreement, Q realizability,
updated weak rates and a source-origin turbulence law remain separate gates.

Use independent rectangular orders32/64 for Q, existing 8/8 phase retraction,
and independent 8/16 final moment checks. Require final scaled six-moment
errors<=1e-8, moment refinement<=1e-5, Q flux refinement<=1e-10. Report any
shape/phase/trust-region failure; do not relax the original .1 shape bound.
Only synthetic known-target and energy-reallocation tests initially. Passing
this helper does not make a real initial Q or its boundary conditions valid.
