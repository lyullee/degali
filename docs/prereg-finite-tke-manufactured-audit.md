# Finite-TKE operator manufactured audit

2026-09-06. This is an operator verification, not a real-trial initialization,
trajectory, turbulence calibration or concentration/temperature score.

## Fixed manufactured inputs

Use the existing test thermodynamics: H2/air at ambient 288.65 K, 101325 Pa,
53.6666667% RH, equilibrium air condensation and temperature-dependent phase
enthalpy. The manufactured axial source is diameter .01 m, velocity 100 m/s,
density 1 kg/m3, temperature 100 K. The local state is
[rho_c,Y_c,A,theta,uc,x,z]=[1.1312,.1537,.02,.01,25,1,1.5], beta_H=1.04.
Use degree-4 center-zero square modes, zero C/H corrections, the inherited
explicit positive radial scalar mixing sampled at order 8, and thermal/species
mixing ratio 1. This is NOT a PRESLHY trial.

The exact-width/wind view uses ustar=.2, zr=.1, neutral stability, zero ambient
spread coefficients and no spread floor. Explicit manufactured Q0=2 J/m3,
all Q shape coefficients zero, ambient k=.07 m2/s2, normalized chi_k=.3 and
tau=.2 s. Fixed direction-preserving circulation amplitudes
[.0001,-.0001,.0001,0] exercise the circulation terms. These are algebraic
test inputs, NOT inferred physical values or validated coefficients.

## Independent checks and fixed gates

- All 5 global, 16 C/H weak and 9 Q weak rows retained; 32 rates with 2 fixed
  kinematic rates. Relative singular-value guard 1e-12.
- Fully phase/mixing-split streaming quadratures (radial/angular) 4/8 and 8/16.
  Scaled matrix, right-hand side and rate refinement <=1e-5, normalized by
  max(abs(fine),1), and cross-grid retained residual <=1e-5.
- At the fine rates each mean-KE, Q and thermal zeroth ledger has error
  <=1e-8 relative to max(1, individual absolute terms). Natural total-energy
  face identity and mass/source identity also <=1e-8 in their analogous scales.
- Independent energy values on coherent exact geometry: stable rectangular
  H values + paired actual EOS/mean-KE values on the union phase grid +
  smooth rectangular Q values. Do NOT use a pre-TKE direction or omit Q.
  Set h=2e-6/max(1,max(abs(fine rates))) so coordinate/shape perturbations
  remain small. Use h/2,h and kinetic grids 4/8,8/16; smooth H/Q order64.
  Each full derivative and four-way spread <=1e-5 of max(1,abs(fine E')).
- Independent 65-ray/order16 diagnostics include Q boundary-gradient mismatch,
  necessary k minus the shear-covariance bound, production, scalar/momentum
  edges, inflow and curvature. These are REPORTED even if they fail. They are
  not acceptance gates for a deliberately manufactured operator test and must
  not be hidden by its numerical-pass status. A separate physical acceptance
  flag stays false, regardless of the numerical result.

Hash all source modules, the two new tests, this specification and the audit
script before execution; verify hashes again on exit. Refuse to overwrite
an existing output or partial checkpoint. Preserve failures. The original
defaults, prior frozen results and measured scores remain unchanged.
