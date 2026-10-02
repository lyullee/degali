# Sandia LH2 pool contour screen: expanded result

Date: 2026-09-17

The public SAND2025-08847 report has been read and its contained-pool
criterion has been separated from spill-off cases. The original Test 3 trace
was removed: despite its low nominal flow, the report explicitly says that
its pool flowed off the substrate, so inflow cannot be equated to a local
vapour source. A contained Test 5 trace was reduced instead, alongside Test
8. Both are local, figure-traceable 4 mol% lower-bound checks against the
existing pure-vapour `LiftoffPlume` pool path.

The first audit appeared to reject the path, but it compared the report's
spatial contour with the conserved **cross-sectional mean** in the integral
plume. That is not the model quantity at a contour point. `LiftoffPlume`
already represents a Gaussian section and exposes its axis peak; the public
trajectory and LFL distance now convert that peak to mole fraction. The
integral mass, momentum, buoyancy and entrainment balances remain on the mean.
With this correction, Test 8 reaches its figure-derived 4 mol% lower bound
after the declared one-grid allowance and remains `not_falsified`. The added,
contained Test 5 lower bound is not reached, even when the model is given the
more favourable of the report's two pool-size assessments. It is therefore
`falsified_below_resolved_reach`. The one-grid allowance is a plotting
resolution, not a fitted tolerance.

This is not an MG/VG/FAC2 result and does not identify a numerical correction.
It does, however, reject the present use of inflow-as-complete-vapour-source
for this contained, cross-wind pool regime. The result applies to the existing
path **as extrapolated** from the historical NASA regime test. That earlier
evidence concerns one 9.1 m pond at about 9--10 kg/s; it does not establish a
small-pool source boundary. No pool source rate, entrainment coefficient,
trajectory parameter or default was changed.

The local reduction uses report figures/panels and one-grid spatial allowance;
neither digitised coordinates nor the report PDF is distributed. Reproduce the
same decision only with authorized local inputs using
`audit_sandia_pool_contours.py` and the protocol in
`prereg-sandia-lh2-pool-contour-screen.md`.
