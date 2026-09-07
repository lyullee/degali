# Pre-registration: independent-energy ground-contact bounds

Date: 2026-09-05

## Physical question isolated before calculation

The accepted seven-state crosswind implementation transports five conserved
fluxes through a full elliptical Gaussian section.  Its observation operator
uses the ground image, but its dynamics still entrain through the full curved
perimeter and apply the full unbounded-fluid buoyancy source even when the
lower finite-domain edge intersects the ground.  In trial 10 the
`source_flux` candidate is already in this regime by 6 m: its centre is about
0.737 m and its vertical sigma about 0.554 m, while the finite-domain lower
edge lies `delta*sigma_z*cos(theta)` below the centre.

This is a mismatch between the dynamic and observation geometries.  It is a
plausible explanation for part of the remaining excessive rise, but the
existing JETPLU ground option previously prevented a low-wind Spadeadam plume
from detaching.  Ground contact is therefore tested only as two explicit
coefficient-free bounds.  Neither is presumed to be a production correction.

The PRESLHY apparatus placed horizontal releases at 0.5 or 1.5 m above the
ground and photographed a low visible cloud.  The experiment therefore
establishes that surface contact is relevant, but does not directly measure a
ground reaction or contact-force law.  No coefficient will be inferred from
the measured trajectory.

## Frozen implementations

For ellipse semi-axes `a=delta*sigma_y`, `b=delta*sigma_z` and trajectory
angle `theta`, contact occurs when

`z_c < b*cos(theta)`.

The cut in the normal cross-section is

`z_cut = clip(z_c/cos(theta), -b, b)`.

Both bounds use the existing analytic elliptical-segment routine for the
free curved perimeter and the existing smooth cleared fraction

`f_clear = clip(z_c/(b*cos(theta)), 0, 1)`.

They retain the full Gaussian conserved-flux definition; only downstream
source terms change, exactly as in the existing JETPLU sensitivity.

1. `geometry`: use the exposed curved perimeter in shear, cross-flow and drag,
   and multiply a positive net vertical-momentum source by `f_clear` while in
   contact.
2. `surface_layer`: apply `geometry` and add the already implemented DEGADIS
   top-surface entrainment over the contact chord,
   `w_e = kappa*u_star*(1+alpha_wind)/Phi(Ri*)`.  This adds no new constant.

`free` reproduces the current independent-energy path exactly.  Interface
projection is unchanged because contact modifies downstream source terms, not
the five boundary fluxes.

## Test order and decision gates

1. Add unit tests requiring exact `free` reproduction, a smaller exposed
   perimeter and positive top entrainment in contact, and exact convergence to
   the free form after clearance.
2. Run trial 10 only with the frozen `source_flux`, enthalpy transport,
   corrected Houf velocity width, local-shear entrainment and 0.02 m maximum
   step.  Require the existing interface pass and maximum downstream balance
   residual below `1e-7`.
3. Report centre height and sigma at 1.78, 4 and 6 m, cumulative buoyancy and
   drag at 6 m, and the first/last contact locations.  A viable direction must
   reduce absolute centre-height error at all three sections.  At 6 m it must
   not worsen absolute vertical-width error by more than 10 percentage points.
4. Only a viable trial-10 direction may be run on all seven interfaces and the
   frozen 42-arc/17-section field comparison.  Promotion then requires all
   previously frozen concentration, width and centre-height metrics to improve
   together relative to the corrected JETPLU baseline.
5. If either bound holds the plume down discontinuously, prevents physical
   clearance, or fails conservation/numerical gates, retain it only as a
   diagnostic.  A calibrated ground-reaction law is outside the available
   measurements and must not be fitted here.

## Result

Both bounds are numerically conservative: their maximum downstream balance
residuals are `3.04e-8`, below the fixed `1e-7` gate.  The free mode exactly
reproduces the prior trial-10 states.  Contact begins at about 2.597 m, so no
ground treatment can change the 1.78 m section in this geometry.

| mode | centre z at 1.78 m | centre z at 4 m | centre z at 6 m | sigma_z at 6 m |
|---|---:|---:|---:|---:|
| measured | 0.4816 | 0.4963 | 0.2266 | 0.5936 |
| `free` | 0.5015 | 0.5672 | 0.7375 | 0.5542 |
| `geometry` | 0.5015 | 0.5631 | 0.6969 | 0.5369 |
| `surface_layer` | 0.5015 | 0.5630 | 0.6955 | 0.5438 |

At 6 m the raw cumulative buoyancy changes little (`5.256 N` free and
`5.188 N` with the surface-layer bound), while the cleared-fraction reaction
reduces cumulative net vertical forcing from `5.189 N` to `3.679 N`.  This
only lowers the predicted centre by 0.042 m; about 0.469 m of the measured
height error remains.  The surface-layer term recovers some of the width lost
when the buried perimeter is removed, but does not materially change rise.

The pre-registered all-three-section improvement gate fails because the
1.78 m section is unchanged.  More importantly, the physical effect is much
too small to resolve the 6 m residual and the section remains in contact at
the end of the 8.28 m simulated range.  Per protocol no seven-trial or field
score is run.  Both options remain explicit off-default diagnostic bounds;
neither is promoted.  The next model change must address the distributed
buoyancy/vertical-momentum structure itself rather than strengthening a
ground-contact multiplier.

Machine-readable values are in
`reference/preslhy/independent_energy_ground_contact_2026-09-05.json`.
