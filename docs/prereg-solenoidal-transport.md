# Conservative solenoidal transverse-flux pilot

2026-09-06. Post-screen response to trial10's reproduced boundary exit and
the failed boundary-row replacement. Before real trial evaluations.

## New assumption, with all existing equations retained

The radial-only particular solution of mass/axial-momentum divergence is
not unique in two dimensions. Add explicit divergence-free vector fields,
without dropping any of the5 global or2m scalar weak transport equations.
For k=1..K let phi'_k(t)=sqrt(4k+1)P_(2k)(t), phi_k(0)=phi_k(1)=0, and

    psi_k(a,b)=a phi_k(b)-b phi_k(a)
    V_k=(a phi'_k(b)-phi_k(a), b phi'_k(a)-phi_k(b)).

Then div(V_k)=0 exactly, normal flux is zero on symmetry axes, and face
normal V_a(1,t)=sqrt(4k+1)P_(2k)(t) integrates to zero. Square symmetries
are preserved. These are normalized transverse flux shapes, not a new
velocity source or arbitrary entrainment addition.

    F_M=(a,b) fM + sum m_k V_k
    F_P=(a,b) fP + sum p_k V_k.

There are2K new UNKNOWN flow amplitudes determined together with all state
rates, not fitted physical constants. Add2K reservoir boundary moments:
C/H against even Legendre0..K/2 (K/2+1 each), P against1..K-2. For K4 this
is3+3+2; for K6 it is4+4+4. Momentum mean remains a checked global identity.
Do not refit source, fields, diffusivity or initial moments.

Use full-vector kinetic flux F_K=u F_P-.5u²F_M and production
P=-grad(u) dot(F_P-u F_M), including the new angular components. This
preserves the corresponding kinetic work identity. Include the new scalar
advection, heat production and boundary kinetic-work contributions in ALL
weak modal equations. The original global equations remain unchanged.

**This permits non-radial momentum stress.** A positive radial projection
or nonnegative shear production is NOT a scalar isotropic eddy-viscosity
closure or proof of realizability of a full Reynolds tensor. Record the
non-radial stress fraction. The old immediate-shear-heating assumption also
remains conditional; no independent TKE/pressure-work/Favre model is added.
The scope is a conditional conservative-flux pilot, never default adoption.

## Fixed cases, resolution and necessary checks

Use trial10's immutable initial field and its independently reproduced
8-division rejected endpoint at0.444045 mm. Reconstruct fields from the saved
state, not the original diagnostic's known mis-scaled physical-position
columns. Same positive chi_C32 input and primary mixing-amplitude update.

For each field, K4 uses radial8/angular48 and16/96; each ray includes current
phase and original mixing cuts. Compare rates/augmented matrix<=1e-5 relative.
Require local linear, mass, zeroth heat and mean momentum<=1e-8; all3 face
residuals<=5%; inward face mass; nonnegative sampled shear production and
radial-projected momentum diffusivity; curvature half-width<.1. Apply the
same necessary signs/edges with independent65 fixed rays/radial16 including
face centers and corners. A negative value is recorded, never clipped.

If all K4 necessary checks pass, repeat with K6. Compare base rate changes
and independent full-vector F_M/F_P values against K4, denominator
max(abs(fine),1), threshold1e-5. This convergence screen is necessary only,
not sufficient physical validation. A successful screen still requires
independent angular phase splitting, actual advective-moment FD and a
resolved momentum-stress/TKE interpretation before new field scoring.
Do not run a new downstream trajectory in this pilot.

Preserve every failure, rank and large-amplitude result; no pseudoinverse,
row omission or relaxed gate. Hash code/tests/spec/input before and after,
refuse overwrite, retain old radial model and failed row-replacement evidence.
