# Joint TKE and axial Reynolds-stress budget

2026-09-06. Registered before evaluating the seven-case tradeoff. This is a
coefficient-free consistency requirement, not a chosen turbulence model.

## Physical terms that cannot be silently conflated

The Favre mean momentum flux includes rho*R_ss, where R is the velocity
covariance. Its contribution to axial stress work is rho*u*R_ss. The total
energy flux already includes rho*u*k, but that does not replace stress work.
The sign conversion is important: NASA TMR writes tau=-rho*R and places the
stress divergence on the momentum RHS. Moving it to the advective side gives
the positive rho*R contribution used here.
[Primary Favre momentum/energy equations](https://tmbwg.github.io/turbmodels/implementrans.html).

The current reduced operator explicitly omits axial normal-stress transport.
Its numerical energy identity is valid for its declared approximation, not
proof that the omitted physical term is small. A pressure adjustment can absorb
some isotropic stress in special approximations, but no such pressure solution
has been supplied by the fixed-ambient-pressure reduced model. Do not claim
the necessary stress-work term is automatically an extra observed net flux.

## Derivation from positive semidefinite covariance

Let r=(R_sy,R_sn), s=|r|, a=R_ss and k=trace(R)/2. The Schur complement gives

    2*k >= a + s^2/a                 (a>0).

At specified k>=s the exact possible axial-variance interval is

    k-sqrt(k^2-s^2) <= a <= k+sqrt(k^2-s^2).

An explicit PSD completion with a>0 has transverse block
r*r^T/a + (2*k-a-s^2/a)*I/2. It demonstrates mathematical existence only,
not measured anisotropy or a Reynolds-stress transport closure.

Use positive advective weight w=A*rho*u*dq_area. Define
B=integral(w*s), K=integral(w*k), N=integral(w*a).
Weighted Cauchy-Schwarz then gives

    2*K >= N + B^2/N.

Consequently K<=Kcap implies N>=Kcap-sqrt(Kcap^2-B^2), when Kcap>=B.
Conversely N<=Ncap<B implies K>=(Ncap+B^2/Ncap)/2. If Ncap>=B the sharp
minimum is K>=B. These are bounds on flux VALUES, not their axial derivatives
and not the net pressure-adjusted force/energy balance.

## Fixed seven-case evaluation

Use only the completed initial_shear_tke_bounds_2026-09-06.json fine/coarse
values, retaining the frozen original geometry/rates. Reconstruct the section
only to establish u(q)>0 over0<=q<=2q0 analytically; do not repeat the large
stress-volume computation. Verify all upstream hashes.

Normalize by the existing mean axial kinetic energy flux Emean. Report:

- N/Emean minimum for K/Emean caps .12, .20, .50;
- K/Emean minimum for N/Emean caps .01, .05, .10;
- the corresponding coarse/fine bound difference, requiring<=1e-3 to call
  the diagnostic resolved. All cases retained, including incompatible caps.

These caps are transparently chosen approximation-budget illustrations, NOT
measured k, fitted constants, an acceptance threshold for observations or an
initial condition. The existing mean field may change under a full closure;
the bounds apply to the fixed shear being examined, not to every possible
future model. No default or measured score is changed.

## Limited code addition

Implement the covariance interval/completion and the associated axial
momentum and stress-work flux contributions with strict input validation.
Do not auto-insert rho*R_ss in momentum without the matching stress work,
pressure and production terms. A future coupled operator must derive these
consistently before trajectory integration.
