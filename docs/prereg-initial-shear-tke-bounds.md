# All-seven initial shear/TKE necessary-bound extension

Registered2026-09-06 before evaluating the seven initial states. Use exactly
the completed coupled_initialization_combined states, scalar coefficients and
fine rates. Reconstruct their original measured pipe/HEM sources and matched
geometry; no exact-geometry substitution, circulation, re-fit or new solve.

Apply the same coefficient-free PSD bound and physical metric conversion as
prereg-shear-tke-lower-bound.md. Set all additional circulation amplitudes to
zero. Evaluate requested8/48 and16/96 radial-phase/angle meshes. The existing
mesh includes extra diagnostic rays with zero quadrature weights; report actual
ray counts separately from the requested Gauss order. Require integrated lower
flux refinement<=1e-3 for a resolved diagnosis, never a physical accuracy pass.

Report all seven, including unresolved/failing reconstructions. A positive lower
TKE flux does not measure actual k or its axial derivative and does not alone
invalidate the existing immediate-heat/local-equilibrium approximation. The
purpose is to distinguish a trial10-specific finding from a broader omitted
storage/supply requirement before deriving a k/epsilon transport extension.
No new constants, source calibration, observed score or default change.
