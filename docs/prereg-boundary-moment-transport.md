# Boundary-moment row replacement: local research screen

2026-09-06, after the trial10 downstream failure; before trial calculations.
Synthetic tests already reveal a rank-deficient zero-shape case and poor
numerical conservation for a nearly redundant case. These warnings are not
removed by loosening singular-value or conservation thresholds.

## Motivation and exact change

Weak scalar inflow conditions do not ensure small pointwise boundary errors
in the finite modal system. Test direct boundary equations in the rate solve.
Keep all5 global conservation equations, physical fields, supplied positive
chi_C, ratio1, radial conservative fM/fP, velocity shape and immediate shear
heating. Change only which finite weak equations determine the local rates.

On one symmetry face a=1,b=t, define actual physical residuals:

    E_C = Y fM - A rho chi_C Y_a/(2q0)
    E_H = h fM - A rho chi_C h_a/(2q0)
          - (.5 Ua^2+.5 u^2) fM + u fP
    E_P = fP - Ua cos(theta) fM

They are affine in the full local rates. Impose integrals of E_C/E_H against
P0,P2,P4 and E_P against P2,P4 on0<=t<=1 (eight equations). Use normalized
even Legendre tests and physical dt, NOT the pointwise diagnostic denominators.
Momentum P0 is excluded because global projected momentum/mass conservation
already determines it for the present source; verify its residual separately.

Drop4 highest ordered raw center-zero polynomial weak equations for EACH
scalar. Order raw pairs by(i+j,i,j). Whitening mixes polynomial degrees, so
first undo the test whitening; do not just remove final whitened rows.
The5 global and2(m-4) retained modal rows plus8 boundary rows give5+2m rows.
No least-squares weighting, pseudoinverse fallback or fitted physics constant.
Record ALL omitted raw weak residuals as explicit approximation defects.
Never claim the original full weak PDE is satisfied by this altered solve.

This is a boundary-bordering prototype informed by
[Dedalus tau documentation](https://dedalus-project.readthedocs.io/en/latest/pages/tau_method.html)
and [Burns et al.](https://arxiv.org/html/2211.17259), not an implementation of
their complete schemes. The former explains explicit polynomial residual
corrections/row replacement; the latter treats elliptic hypercube corner
compatibility. Neither proves this cryogenic, nonlinear conditional transport
discretization is well posed. Eight rows and their ordering are a FIXED
numerical trial, not a derived universal physical closure.

## Fixed local cases and gates

Use immutable independently passed initialization states10 and23. Trial10
was selected because of downstream edge exit;23 had the smallest initial
positive chi_P and the most marginal short-segment momentum boundary.
This selection is post-screen and not a new seven-case validation claim.

At each unchanged field evaluate fresh radial8/angular48 and16/96 meshes;
each ray splits current phase and supplied mixing knots. Compare physical
rate and matrix relative differences<=1e-5. Recompute all diagnostics using
NEW rates. Require global/retained-modal/linear/zeroth-heat/mass/P0-momentum
errors<=1e-8, all3 pointwise edges<=5%, nonnegative chi_P, inward face mass,
curvature half-width<.1; independently verify65 fixed rays/radial16.

If those local gates pass, embed exactly the same scalar fields into degree
N+2 (raw-coefficient embedding, no refit) and repeat. Compare relative C,H,
rho and u spatial rate fields on17x17 fixed points, denominator max(|fine|,1),
and base7 rates, requiring<=1e-5. Failing any stage precludes downstream use.
Passing all these is still only a candidate for separate fully phase-split
angular verification, actual-moment FD and further degree convergence.
Do not adopt, run new trajectories or field score in this screen.

Report rank/conditioning failures, pointwise defects and omitted weak
residuals for each case, including rejected runs. Freeze hashes before/after
and refuse overwrite. The existing natural-boundary ODE and its failures
remain unchanged; this prototype is in a separate module.
