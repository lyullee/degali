# Bounded three-dimensional wind candidate — 2026-09-06

This is a new opt-in extension of the frozen pressure-loss density CONTROL,
not a default replacement or fitted yaw. The first actual test is trial10,
using the independently valid matching49sample wind-FROM276.3043110304036deg
and nozzle bearing75deg. Keep the original wind magnitude/height dependence
to isolate direction, and report that magnitude/station uncertainty separately.
Trial23 far1Hz is faulty and MUST NOT be used. Coarse N/W/N brackets are not
precise contemporaneous wind and cannot be chosen by the temperature result.

Carry mass, H2, Px, Py, Pz and total energy. Density/H2 Gaussian and axial
velocity shapes and original thermochemistry stay unchanged. Velocity is
u=U dot n + uc exp(-lambda²q), and section momentum=P*n. Axial ambient
coflow only: reject U dot horizontal(n)<0 rather than clipping reverse flow.
The current implementation does not represent recirculation, ambient normal
velocity inside the finite section, Reynolds stresses or resolved TKE.

Use n=(cos(theta)cos(psi),cos(theta)sin(psi),sin(theta)). Entrained ambient
momentum is mU, energy m|U|²/2 in the existing ambient-enthalpy reference.
Normal drag is Cd*rhoa*perimeter/2 * |Uperp| Uperp, where
Uperp=U-(U dot n)n. Buoyancy is unchanged and vertical. This recovers the
original drag at psi=wind_angle=0. Local shear/density-scaled/Houf entrainment
is retained with centre velocity uc+U dot n. Generalize the original crosswind
term to alfa2*(U dot n)*|Uperp|/|U|*perimeter (zero at U=0). This is a stated
rotation-covariant empirical extension with the same coefficients, not a newly
validated universal crossflow law. Its perpendicular-coflow limitation remains.

Track horizontal travelled distance L, Cartesian X/Y/Z with dL/ds=cos(theta),
d(X,Y,Z)/ds=n. Use L for original ambient width growth; this is a modelling
assumption for yawed flow, retaining the old coplanar limit and horizontal
rotation covariance. Ground handling remains the original optional treatment.

Project unchanged source six fluxes (old Py=0), at unchanged coordinates and
lambda. Preserve existing1e-8 flux,1e-5 quadrature,2K temperature and0.05
relative halfwidth interface gates, measured relative to the ORIGINAL source
target record, not to an already approximate projected state. No new source
formation physics is claimed. Integrate only an accepted boundary. Use the
original conservative flux RK4 construction extended to6fluxes/4coordinates,
step.02m, then.01m if supported. Save10-step checkpoints and exact failures.
Require mass/species/momentum/energy ledger maximum normalized residual<=1e-5.

Receptor mapping uses the unique horizontal-normal section whose horizontal
tangent is perpendicular to receptor-minus-centre. Interpolate stored state,
then evaluate original vertical/ground-image profile using horizontal-normal
lateral offset. This preserves the old zero-yaw receptor convention exactly;
it is NOT a new exact tilted-ellipse observation operator. Require a unique
bracketed section and report missing coverage, never nearest-point extrapolate.
Angles must remain unwrapped within the supported coflow branch.

Before actual results, test zero-yaw flux/source/point limits, horizontal
rotation and reflection covariance, vector drag orthogonality, source inversion,
unmodified base object and invalid/reverse inputs. Compare all original trial10
5concentration groups/3geometry rows/19temperature sensors, minimum/p05/median
AND the newly reduced same-window means. Retain sensor-level outputs and
centre/left/right residuals. A temperature gain with concentration/geometry
damage is not overall adoption. The known near-centre cold-gap diagnosis remains
valid; yaw is not assumed to solve it. No all7-trial or full-system claim from10.
