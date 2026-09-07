# Inequality-constrained flux feasibility witnesses

2026-09-06. After both fixed-field solenoidal hard-boundary-moment pilots
failed: they retained conservation but generated very large oscillatory
stress, negative production and failed refinement. Before this new test.

Question: is there ANY small K4 solenoidal correction preserving all global
and original scalar weak equations while meeting the existing inequalities?
This is an existence/feasibility diagnostic, NOT a new turbulence closure.
The minimax amplitude objective is a numerical witness selector, not an
asserted physical minimum-energy or minimum-dissipation principle.

Keep trial10's initial and rejected8-step endpoint fields, source and chi_C.
Use the saved immutable K4 matrices from the preceding pilot. Eliminate the
original state rates with their original nonsingular weak matrix:

    r(alpha)=r0-M_active^-1 M_flux alpha, with kinematic x/z unchanged.

All local boundary residuals, full-vector shear production and mass flux
are then affine in the8 mass/momentum redistribution amplitudes alpha.
Require on the declared meshes:

- |normalized C/H/P edge residual|<=.04, the existing initialization target;
- nonnegative shear production everywhere sampled, not clipped;
- inward face mass<=-1e-6*|old reference face mass| (numerical strictness margin);
- |theta_dot|*physical half-width<=.1.

Minimize max(|alpha_i|/s_i) using a linear program. Mass amplitude scale is
source_M/(16q0); momentum scale multiplies this by max(center axial speed,
ambient speed,1 m/s). These scales define the numerical witness metric,
not new physical constants. Do not silently use returned optimal amplitudes
as a downstream model. A successful LP must still pass independently
recomputed actual residuals/production, not just solver status.

Evaluate radial8/angular48 and16/96 K4 meshes with65 extra cosine face
points, same current phase cuts as the saved pilot. Verify each successful
fine witness with65 uniform-angle/radial16 rays and original5% edge gate,
nonnegative production/radial diffusivity, inward mass and curvature. Require
linear/heat/mass and LP feasibility residuals<=1e-8. Record coefficient/rate
changes across resolutions; do not treat this LP selector as a smooth ODE.

An initial zero-amplitude witness is a useful consistency check. Failure or
infeasibility is retained. Success would establish a sampled conservative
alternative exists, not isotropic stress, Reynolds-tensor realizability,
TKE closure, stable time integration, or improved observed accuracy.
Freeze hashes; no overwritten failed evidence/default promotion/field score.
