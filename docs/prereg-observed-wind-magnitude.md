# Same-window observed wind magnitude and direction — 2026-09-06

The yaw-only trial10 run used the valid49sample vector direction but retained
the old wind magnitude/height reference deliberately for mechanism isolation.
Its actual19sensor mean/median errors decreased while minimum error and vertical
width increased. Do NOT retrospectively choose another station, direction,
averaging window or speed to optimize those outcomes.

Complete a two-factor contrast with exactly the same independently selected
49sample far sensor: vector mean norm1.9222640246473406m/s at3m and relative
yaw21.30431103040358deg. Mean vector norm, not mean scalar speed, represents
the mean momentum input. Distinct time-varying kinetic energy/entrainment effects
are NOT included by that choice. The instrument is about20m from release;
spatial representativeness and near-station disagreement remain uncertain.

Four combinations are original magnitude/direction (existing CONTROL), original
magnitude/observed direction (completed yaw pilot), observed magnitude/original
direction (new magnitude-only), and observed magnitude/observed direction (new
measured-vector). Use all four to separate input effects, not select by scores.
Update u0,reference height z0 and friction velocity together with the existing
atmospheric log-profile/stability law. Retain roughness, stability and width
coefficients, source six fluxes, lambda and thermochemistry. The adapter must
copy base/JetPlume, never mutate the original. At identical reference inputs it
recovers the original wind/flux/source exactly.

The input choice is prescribed from the raw-wind audit, with no adjustable
scalar multiplier. Handoff: unchanged original target1e-8 six-flux,1e-5 energy
quadrature,2K centre-temperature,5%halfwidth gates. Unsupported boundaries stay
unsupported. Actual trial10 .02m for both new cells and .01m measured-vector if
accepted; original19T/5C/3vertical keys, min/p05/median/mean and every per-sensor
record, same conservative1e-5 ledger. Compare .02/.01 sensor differences without
silently passing a threshold selected after seeing scores. Existing yaw-only
fine run interrupted at960 then resumed separately must retain its provenance.

This is a downstream same-source contrast, not upstream source formation with
new wind. It does not resolve initialTKE, nonstationary gust propagation,
anisotropic entrainment/second moments, or the known warm near-core residual.
Trial23 faulty far sensor remains excluded. No all-seven/default adoption from10.
