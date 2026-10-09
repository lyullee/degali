# 3-D source/wind-coupled transient dense-gas analysis

`degali.addons.transient_dense_gas_3d` adds an explicit reduced-order
finite-volume solver for a Cartesian three-dimensional domain. It is intended
for event replay and sensor timing, not for replacing a variable-density
LES/RANS solver.

The solver couples the two histories on the same event clock:

1. `SourceRateSchedule.mass_between(t, t+dt)` injects the already-atmospheric
   H2 mass during every numerical step.
2. `WindHistory` is linearly interpolated at the step midpoint and converted
   from meteorological “from” bearings to the Cartesian downwind vector.
3. The H2 field is advected by that time-varying vector and diffused with an
   explicit conservative finite-volume update.
4. By default, a prognostic velocity state relaxes toward the measured wind,
   integrates the density-driven vertical buoyancy acceleration, and feeds the
   resulting 3-D velocity back into scalar advection. Set
   `prognostic_velocity=False` only for the legacy prescribed-wind diagnostic.
   Dense source gas sinks; lighter gas rises.
5. Receptor traces, cumulative injection/outflow and mass residuals are
   emitted. Optional full field snapshots can be retained for diagnostics.

The source must already be an atmospheric H2 boundary and must end with an
explicit zero-rate record. A `SourceStateLedger` can be connected through
`SourceStateLedger.to_source_rate_schedule()` after phase routing.

## Scope boundary

This is a 3-D transient transport model with a prognostic velocity closure, but
it remains a reduced-order closure. It does not resolve pressure waves,
multiphase flashing inside the grid, turbulence-resolved eddies, or a full
variable-density Navier–Stokes pressure projection. Such a solver would
require a separate validated CFD backend. The emitted diagnostics therefore
retain the solver scope and should remain `conditional` until matched
transient field data are available.

The first validation targets are source shutoff timing, wind-direction change,
arrival-time/peak detection, dense-gas sinking, and mass conservation—not a
single tuned pointwise MAE.
