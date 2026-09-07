# Recompute source centers/width with prescribed Q and normal stress

2026-09-06. Registered after preserving the failed fixed-center Q0=2 normal
source retraction. That failure was a bounded QP failure, NOT proof that no
physical source can exist. Neither the old bounds nor its result are edited.

## Explicit change of numerical source degrees of freedom

Keep the original six enriched source targets [M,H2,Pax,Etotal,B,Hsecond].
Keep x,z and direction fixed. Q0, Q shape and R_ss/k are explicit fixed inputs;
they are not optimized. Assume ambient pressure without an invented pressure
compensation field, as in the immediately preceding boundary-only test.

Allow six previously fixed quantities to adjust simultaneously:

    log rho_c, log C_c, log area, log excess u_c,
    C log-shape coefficient multiplying q,
    H log-shape coefficient multiplying q.

The latter two directions are represented exactly in the existing center-zero
square basis; no new velocity/scalar reference-width gauge is introduced.
Changing area changes the velocity width as well as geometry, whereas changing
C/H q coefficients changes their width relative to velocity. These are distinct.
This is matching the source moments, not fitting any observation or turbulence
input. The centers and width are NOT claimed unchanged in the new method.

Use coherent exact width/wind values and derivatives at each candidate. The
six actual moments include normal momentum A int(r Q), Q energy A int(Q u)
and normal stress work A int(r Q u). H and its second moment remain thermal;
the EOS uses H/C only. Differentiate buoyancy and physical H-second-moment
weights consistently with density and the changing widths.

## Numerical trust region and required checks

The four center/area/velocity log changes are limited to +/-0.1 from the
provided start. Retain the original sampled .0999 C/H log-shape bounds and
strict .1 field bounds. These are numerical search limits, not empirical
physical uncertainty estimates; a miss is retained as a bounded search miss.
Use at most20 Newton iterations, rank guard1e-12, and backtracking factors
1,.5,.25,.125,.0625,.03125,.015625. No state clipping or relaxed conservation.

Use phase-split 4/8 for the first solve, then refine at8/16. At the final state
require six-moment target error<=1e-8 and4/8-to8/16 quadrature difference<=1e-5.
Independent actual +/- source moment differences check the 6x6 Jacobian in
manufactured cases, including the changing physical second-moment weights.
The exact-geometry Q0=.02 and2, r=2/3 cases are both retained and reported;
neither is selected as a real LH2 initial condition.

Do not return or reuse old transport rates after moving the source. New
scalar mixing, Q/epsilon/circulation laws, full normal-stress/pressure transport
and all boundary/PSD gates must be recomputed before any real trajectory.
