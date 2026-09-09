# Pre-registration: radial-core constitutive domain

Date frozen: 2026-09-10, after locating the Trial 10 sign loss and before
implementing or marching the radial-core candidate.

## Physical issue

The reduced radial stress/diffusion identity is derived by integrating over a
complete constant-q contour.  In the finite square Gaussian quadrature domain,
that contour is complete only through the inscribed circle `q <= q0`.  For
`q0 < q <= qmax`, the contour is cut by the square sides and the balance also
contains flux through those sides.  Applying the complete-circle cumulative
identity there does not define a local eddy diffusivity.

The stopped Trial 10 calculation first reports negative momentum diffusivity
only in this clipped corner: at the rejected half-stage, 1 of 2,136 samples is
negative, and the sign change lies between 0.9999 and 0.99999 q/qmax where the
remaining angular measure falls below 4.0e-4 rad.

## Candidate

Add an explicit `radial_core` constitutive-validity option.  It checks species
and momentum diffusivity on `[0, q0]`, including every phase-cell midpoint in
that interval.  It does not clip, floor or change either diffusivity.

All square-domain conserved-flux integrals, moment integrals, phase cells,
reservoir side fluxes, edge heat defects, incoming-mass gate and source terms
remain unchanged.  The existing `full_square` option remains the default until
this experiment is complete.  The candidate therefore changes only where the
complete-circle constitutive identity is declared meaningful; it does not
discard the square-corner mass or energy.

## Frozen test

Repeat the pre-registered Trial 10/23 observation calculation to x=4.00 m with
the same boundaries, linear phase lookup, `D_h/D_C=1`, reduced buoyancy work,
0.005/0.0025 m nominal steps, inverse tolerances, sensor coordinates,
observation reductions and numerical gates in
`prereg-thermal-moment-observation-validation.md`.

Additional requirements:

- record the selected constitutive domain in every result;
- require nonnegative species and momentum diffusivity throughout the sampled
  radial core at every RK stage;
- retain the full-square corner diffusivity as a diagnostic and do not call it
  a physical transport coefficient;
- require both trials and both step sizes to reach 4 m, independent balance at
  most 1e-5, terminal and 40-sensor refinement differences at most 0.005;
- require the same componentwise improvement of thermal centre-amplitude and
  thermal-variance log errors in both step runs.

If any gate fails, retain `full_square` as the default and do not tune a
coefficient.  If all pass, the option may advance as a research candidate but
is not a validated safety-design default.  Raw experimental files remain
local and undistributed.
