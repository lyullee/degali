# Explicit axial normal-stress source matching

2026-09-06. Boundary-only extension of prescribed_tke_source, before its tests.
This is not a downstream stress/pressure or epsilon transport closure.

At a supplied source plane let Q=rho*k and let the caller explicitly prescribe
r=R_ss/k, a finite scalar in [0,2]. This interval follows from the nonnegative
trace but is not sufficient for a shear-bearing covariance. Do not select r
from a desired observation score. The equal-normal-variance value r=2/3 is
used only in manufactured tests, not adopted as LH2 anisotropy data.

Require an explicit fixed-ambient-pressure/no-pressure-compensation assumption.
The current EOS is not a variable-pressure mean-field solver, so a compensating
pressure cannot be silently invented. Retain both contributions:

    delta axial momentum = A integral(r*Q)
    delta axial stress work = A integral(r*Q*u).

The total source energy now includes H+meanKE+Q+r*Q*u integrated as appropriate.
Subtract the supplied normal momentum and normal work from the original target,
then reuse the existing explicit-Q energy retraction. Verify the original six
total target moments after adding all contributions back. Do not create energy,
adjust a source flow rate or mutate the original projection. A failed bounded
retraction is reported, not repaired by relaxing the old .1 log-shape bound.

Use the same independent32/64 smooth integrals and8/8-to8/16 phase checks as
the Q-only helper, with final six-moment error<=1e-8 and refinement<=1e-5.
Check known manufactured totals and reallocation inside unchanged totals.
Report numerical success separately from physical initial-boundary acceptance.

This does not establish r, Q, the remaining covariance, pressure transport,
normal production or their axial evolution. All must be consistent before
using the result in a physical trajectory. The previous Q-only audit remains
valid for its explicitly reduced approximation and is preserved unchanged.

## Exact-geometry manufactured audit cases

Registered before this separate audit: reconstruct the same exact-geometry
manufactured section as finite_tke_manufactured, keep its original mean-only
six target moments, and evaluate Q0=.02 and Q0=2 with zero Q modes and r=2/3.
Use both predetermined levels, including any infeasible larger case; do not
select the smaller level as a real initialization. Run the same bounded
retraction (maximum8 iterations) and independent final moment checks. Preserve
per-case failures, verify hashes and prohibit overwriting an existing result.
