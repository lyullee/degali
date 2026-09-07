# Isolate actual enthalpy-flux cancellation

Registered2026-09-06 before the trial10 calculation. The paired/polished full
energy difference FAILURES remain unchanged. This is a component diagnosis,
not an alternative way to label those full-energy gates as passed.

The enthalpy profile is smooth: H=Hc exp(-q/beta²+psi.pH). Consequently its
advective flux has two separable terms AHc Wparallel and AHc uc exp(-kappa q).
Evaluate the actual two prefactors in50-digit mpmath1.3.0, using the SAME stored
binary64 bilinear phase table and coefficients, and exact mathematical shifts
of the stored binary64 encoded origin/direction. Difference these scalars before
conversion to binary64. At each fixed quadrature point use the exact identity
Bplus exp(delta)-Bminus exp(-delta)=DeltaB cosh(delta)+2 meanB sinh(delta).
No EOS or enthalpy-flux derivative is used to form these actual-value differences.
The square basis/quadrature and vector evaluation remain binary64; this is NOT
an all-field arbitrary-precision run. Optional mpmath is research-only.

Two explicitly separate wind policies: binary_legacy retains actual existing
wind values at each shifted state; exact_constraint evaluates the same width
equation by its positive algebraic root and the same log-profile averaging
formula in arbitrary precision, preserving constants/floors/stability branches.
Measure the change in wind and widths. The latter is a numerical-reference
variant, not a silent change to the legacy Fortran oracle or adopted model.

Use original h/2,h. Compare tensor Gauss32/64, and50/70 digits. Independently
integrate the existing analytic BH@rates on64/96. Report absolute errors and
errors scaled both by this enthalpy derivative and by the original total-energy
derivative, whose cancellation scale matters. A component diagnosis is adequate
at1e-5 total-energy-scaled error and1e-5 refinement, but can never establish a
complete energy verification without the independent kinetic contribution.
Record scalar center/wind differences and differences from the existing tangent.
No observations, defaults, targets, fitted parameters or downstream integration.
