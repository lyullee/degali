# Coherent exact geometry: full energy directional verification

Registered2026-09-06 before real-case execution. Trial10 failed physical
state and the frozen phase-consistent witness direction are kept. New geometry
values AND their derivatives must use the SAME exact width constraint and
log-profile wind average. Preserve original JETPLU and all earlier audit hashes.

Use a read-only geometry wrapper: solve sy*sn=A and sy²-sn²=sya²-sna² by the
positive algebraic root; retain spread floor, ambient spreading, wind constants,
stability branches and averaging. Use expm1 for the difference of powers.
Compute logA/theta/x/z tangents with complex step1e-24. Independently test these
against60-digit mpmath values/derivatives on neutral/stable/unstable cases.
Quiescent wind does not imply equal widths. No derivatives are claimed exactly
at the x=0 or center roughness-height branch switches.

Construct a coherent new field view. On phase-split4/8 and8/16 streaming
quadrature recompute the actual advective Jacobian, retained operator and
source. The OLD direction is deliberately not re-solved: it need not conserve
the NEW source. Report that residual; do not call it a new conservative witness.
This audit isolates differential consistency before any new flow solve.

Compare the full energy derivative along that fixed direction to actual
H-flux plus kinetic-flux central differences, at original h/2,h. H uses the
registered smooth64-point70-digit scalar-prefactor method with exact-constraint
wind. K uses actual plus/minus binary64 density with same-table root polishing,
common union phase grids and high-precision scalar area/speed differences.
No analytic EOS derivative forms the actual difference. Density is NOT
arbitrary precision, and the original frozen full-energy failures remain failed.

Run kinetic differences4/8 and8/16 at each distance. Require relative error
and all mesh/step refinements<=1e-5, denominator max(abs(new expected E'),1).
Require complex-step1e-24/1e-28 geometry tangent difference<=1e-10. Report
original-versus-exact width/wind changes, new-versus-old E' and the source
residual. No relaxed gate, fitted parameters, observation score, default change,
downstream integration or physical circulation-closure claim.
