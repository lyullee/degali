# Rate-free boundary compatibility diagnosis

2026-09-06. After trial10's natural-ODE boundary exit and the rejected
boundary-row prototype; before this diagnostic evaluation. No new fitting.

Three boundary residuals C,H,P depend on only two local transverse fluxes
fM,fP when scalar fields, velocity and diffusivity are fixed:

    E_C = Y fM - dY
    E_H = (h - .5 Ua^2 - .5 u^2) fM + u fP - dh
    E_P = fP - Ua cos(theta) fM
    dY = A rho chi Y_a/(2q0), dh = A rho chi h_a/(2q0)  (ratio1)

Eliminating fP via E_P=0 and fM via E_C=0 gives the NECESSARY profile relation

    dh = h_eff dY/Y
    h_eff = h - .5[(u-Ua cos(theta))^2 + Ua^2 sin(theta)^2].

For the normalized residuals with positive scales S_C,S_H,S_P, a left null
vector of their3x2 flux matrix is n=(-h_eff*S_C/Y,S_H,-u*S_P). Consequently,
the sharp unconstrained minimum maximum normalized residual is

    |dh-h_eff*dY/Y| / (|h_eff*S_C/Y| + S_H + |u*S_P|).

A three-variable linear program additionally imposes fM<=0 and fP-u*fM>=0,
the latter being nonnegative inferred chi_P for the present u_q<0. Its
optimum is also only a LOWER bound for a conservative interior solution.
The locally optimal flux values are diagnostic, never inserted into the
transport equations. A bound above5% rules out all local flux adjustments
under these unchanged assumptions at that witness; a smaller bound is not
proof that a globally conservative field or downstream trajectory exists.

Evaluate1025 uniform face points including face center and corner for all7
immutable initialization states and the reconstructed8-division rejected
endpoint of trial10. Use existing normalized scales and positive original
chi_C input; use the primary amplitude update only for the moved state.
Record field values, gradients, minima, witness positions, LP residuals and
positivity. Exact physical position is sqrt(2q0)*(sy*a,sn*b).

Verify left-null and optimum residual identities, retain original failures,
hash code/input before/after, refuse overwrite. No acceptance relabeling,
diffusivity refit, downstream integration, observed score or default change.

This algebra follows from the current three reservoir conditions; it is a
diagnosis of their joint implications, not new evidence that all pressure,
Favre velocity or turbulence-energy terms of LH2 dispersion are represented.
