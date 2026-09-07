# Necessary turbulent-energy storage from the reconstructed shear stress

Registered2026-09-06 before a real-trial calculation. This is a no-fitted-constant
physical consistency diagnosis, not a k/epsilon model or accuracy improvement.
Use the completed legacy-geometry aligned witness first; keep its fields/rates/
four amplitudes fixed. A later exact-geometry result must get its own run name.

The density-weighted velocity covariance is positive semidefinite. Its trace
defines twice k; the off-diagonal Reynolds stresses are density times shear
covariances. Definitions follow the [NASA TMR Favre equations](https://tmbwg.github.io/turbmodels/implementrans.html).
The bound below is derived here, not a fitted relation taken from that source.

In normalized transverse coordinates, let S=FP-u FM, and metric lengths
Ly=sqrt(2q0)*sy and Ln=sqrt(2q0)*sn. Then R_si=Li*S_i/(A*rho). Positive
semidefiniteness implies k>=sqrt(R_sy²+R_sn²). Proof: rotate the transverse
plane toward the shear vector r. The2x2 principal minor requires
Rss*Rpp>=|r|², hence trace(R)/2>=(Rss+Rpp)/2>=|r|. It is sharp: a rank-one
covariance with diagonal (|r|, r_y²/|r|,r_n²/|r|) attains the bound.

Evaluate the actual reconstructed total stress, NOT just the new circulation
increment. Integrate Kmin=A integral(rho*u*kmin), the axial flux lower bound,
and compare with the mean axial KE flux. Report maximum pointwise lower k,
total fluctuation RMS lower bound sqrt(2kmin)/u, and the covariance metric.
This does not establish the actual tensor, actual k, isotropic intensity,
dissipation, initial turbulence level, pressure closure or extra mean transverse KE.

Use independently radial-phase-split fixed angular meshes8/48 and16/96,
matching the existing constraint-mesh construction. Require integrated Kmin
refinement<=1e-3 for a resolved diagnostic; otherwise record unresolved and
do not escalate the numerical estimate into a physical bound integrated exactly.
The pointwise algebraic inequality is independent of quadrature convergence.
No physical pass/fail cutoff on the Kmin/meanKE ratio is invented. No defaults,
observations or downstream run. If bound>0, k=0 and those nonzero stresses
cannot simultaneously represent a positive-semidefinite Reynolds covariance.

Not storing k explicitly does NOT necessarily mean a model asserts k=0:
local production-dissipation equilibrium can have nonzero k. A flux lower
bound is not a lower bound on its axial derivative, and cannot by itself
prove that the equilibrium approximation causes a specified prediction error.
