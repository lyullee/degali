# Post-screen diagnosis of trial10's downstream boundary exit

2026-09-06. Written after the primary2/4/8-step trial10 failures, before this
diagnostic run. This is not a new acceptance screen and does not relabel the
original rejected segments. Trial selection is explicitly post-screen.

Input: immutable completed `enriched_segments_a_2026-09-06.json` plus its
90 frozen dependencies. Reconstruct the actual source and accepted initial
shape. For each failed4/8-division primary record, start at its last ACCEPTED
state and reconstruct the identical rejected stage: midpoint for4 divisions,
endpoint for8 divisions. Disable only the
boundary/constitutive acceptance switch for the diagnostic candidate; retain
EOS and logshape evaluation/rejection. Store start, midpoint, rejected
endpoint, rates and distances, not just the last accepted point.

At each reconstructed endpoint compare the usual fresh-mesh8/48 operator
with independent moving-phase4/8 and8/16 operators. Recompute65 fixed-angle,
radial16 boundary witnesses with the fine rates, saving physical terms,
normalization scales, coordinates, inferred chi_P and mass inflow. Compare
the two endpoints at their common requested distance, including encoded
state difference and peak residual difference. The4-division midpoint is an
explicit internal stage, not a completed second-order endpoint; this is NOT
a Richardson convergence estimate. No temperature, state, source
or diffusivity refit, stress clipping, altered acceptance threshold, adaptive
shortening, observed scoring or default promotion.

A matched endpoint and refined edge failure supports a real exit of this
finite-mode conditional ODE from the imposed boundary tolerance, not failure
of every possible physical model. It alone does not distinguish missing
transverse redistribution, inadequate modal boundary enforcement, velocity
shape, or turbulence closure. As a local sensitivity, evaluate the SAME
8-step failed endpoint with constant-geometric-diffusivity mixing. This is
not an alternative trajectory or an adoption option. Record both results.

Compare refinement with existing1e-5 numerical and5% edge reference limits,
but do not redefine passing status for the old run. Hash all dependencies
before/after, preserve diagnostics even when they fail, and refuse overwrite.
