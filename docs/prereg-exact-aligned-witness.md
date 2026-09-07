# Re-solve conservative aligned flow on coherent exact geometry

Registered2026-09-06 before calculation. Run only after the full directional
energy gate in exact_geometry_energy passes. Keep all earlier failures and
the old-direction/new-source mismatch. That mismatch is not repaired by
reinterpreting a source: solve all retained equations again on the new geometry.

Use the same frozen trial10 failed physical state/scalar coefficients/source
reconstruction, prescribed mixing and four aligned circulation modes. This is
a LOCAL conditional flow witness, not an initial-state refit, new trajectory,
physical circulation closure or observed validation. Do not claim preservation
of old six source-plane moments after the small numerical geometry change.

Use the freshly verified exact-geometry state matrices4/8 and8/16. Independently
reintegrate aligned weak columns at both resolutions. Same LP4% edge target,
nonnegative production, strict inward mass and curvature bound. Same amplitude
normalization and objective as the preceding aligned pilot. Retain ALL original
global/weak equations; no row replacement or residual clipping.

Require unit-amplitude mass/momentum/global-heat correction<=1e-7 absolute;
rate/amplitude/column refinement<=1e-5; both cross-volume residuals<=1e-8;
65 independent rays/order16: edges<=5%, production/diffusivity>=0, inward mass,
nonradial fraction<=1e-10. Independent actual H+K differences for the NEW fine
direction must also pass at original h/2,h on4/8 and8/16. Compare to freshly
integrated exact-geometry Jacobian; max(abs(expected E'),1) scale, error and all
mesh/step changes<=1e-5. This time E' must also match the new energy source
within the retained-equation1e-8 gate. No unverified direction substitution.

Do not advance downstream or fit observations on success. The remaining
physical circulation and scalar/TKE mixing laws must still be determined.
