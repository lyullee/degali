# Finite gas--condensate slip: executable conservation kernel

Date: 2026-09-14

## Why this is separate from the default source

The LH2 near field can contain condensed N2/O2 and water-rich ice.  Treating
that inventory as permanently no-slip is a kinematic assumption, while
treating it as immediately stationary is the opposite limiting assumption.
Neither is a measured description of the PRESLHY releases.  In particular, a
particle size alone does not close the problem: a two-velocity calculation
also needs condensed mass, initial gas/particle slip and an applicable drag
law.

`relax_two_velocity_drag` supplies the drag-only conservation part that is
independent of those missing measurements.  For gas mass flow `m_g`, particle
mass flow `m_p`, Stokes relaxation time `tau_p`, and an explicitly chosen
constant-step duration, it advances the relative velocity exactly:

```text
d(u_g - u_p)/dt = -(1 + m_p/m_g) (u_g-u_p)/tau_p
```

The mixture momentum is preserved exactly.  The resolved relative kinetic
energy loss is returned separately as `thermalised_kinetic_energy`; it must be
inserted into a total-energy balance once, rather than silently disappearing
or being added again as latent heat.

`advance_two_velocity_phase_step` adds the immediately preceding control
volume: ambient entrainment at an explicit axial velocity, an externally
solved final condensed inventory, phase transfer, and drag. It returns **three** separate
mechanical-to-thermal terms:

1. entrainment mixing;
2. inelastic velocity mixing when condensing/re-evaporating mass joins a
   phase already moving at a different velocity; and
3. interphase drag.

Their sum equals the resolved kinetic-energy decrease.  Treating only the
last term as heat would miss the second term whenever the phases have slip;
putting any of these terms on top of an enthalpy state that already absorbed
it would double-count energy.

`schiller_naumann_relaxation_time` is also available as an opt-in local
conversion from declared sphere diameter, particle/gas density, gas viscosity
and current slip to a relaxation time.  It applies the isolated-sphere
Schiller--Naumann multiplier only for particle Reynolds number at or below
1000.  It is not an exact nonlinear trajectory: a caller must update the
local relaxation time as slip changes.  It contains no selected particle size
or rarefaction/turbulence correction.

For a constant-property interval in the same Reynolds-number range,
`relax_two_velocity_schiller_naumann_drag` instead integrates the
Schiller--Naumann relative-speed equation analytically.  It preserves mixture
momentum and records the relative kinetic-energy loss without freezing the
initial response time.  This removes a numerical approximation from the
opt-in kernel; it does not select a particle property, establish that the
isolated-sphere law is applicable to the LH2 cloud, or couple the result to
the default source.

`particle_knudsen_number`, `davies_cunningham_slip_correction` and
`cunningham_corrected_particle_relaxation_time` make the separate low-Re
rarefaction screen explicit when a caller supplies a gas mean free path. The
Davies correction increases the Stokes response time as particle-scale
continuum drag breaks down. DEGALI does not infer the mean free path of a
cold multicomponent cloud and does not silently compound this correction with
finite-Re drag; that combination needs a stated applicability model.

## What is and is not implemented

- Implemented: exact constant-drag gas/particle velocity exchange; an
  explicit-entrainment-velocity phase-transfer/drag control-volume step with
  separate mechanical-heat ledgers; momentum and kinetic-energy identity
  checks; zero-particle, re-evaporation and long-time limits; a bounded,
  declared-input finite-Re local response-time conversion and an exact
  constant-property Schiller--Naumann relative-slip update; a declared-input
  low-Re rarefaction screen.
- Already separate: particle settling and heat-limited N2/O2 sublimation
  bounds in `transported_condensed_air_source`.
- Not implemented: a default particle diameter, phase fraction, initial
  slip, variable-property/nonlinear drag trajectory through the actual plume,
  particle coalescence/nucleation, a variable-composition rarefied-gas drag
  law, or a fitted conversion of slip dissipation to a particular LH2 trial.

The final item is deliberate.  The 0.5--1.2 mm solid-air observations of
Shangguan et al. were made by feeding air into a bulk-LH2 dewar; they are not
an airborne free-jet particle distribution.  Applying them to a release jet
would be an uncontrolled transfer of scale and geometry.

The public PRESLHY E3.5 catalogue independently reinforces this restriction:
although it lists H2/O2, temperature, pressure, mass-flow, weather and video
channels and reports condensed air near the release, it lists no particle
imaging, size distribution, phase-resolved inventory, or separate particle
velocity.  The visible cloud is consequently evidence that a condensed phase
can occur, not a calibration observation for this kernel.

The independent DNV/FFI Spadeadam field campaign does not fill this gap. Its
public report documents oxygen-depletion gas sensors and thermocouples at
30/50/100 m, with 10 m/5 m mast weather measurements, but no particle image,
phase inventory, velocity profile or separate gas/particle velocity. It is an
important downstream concentration/temperature check, not a two-velocity
adoption dataset.

## Adoption gate

This kernel can be coupled to the transported source only when one matched
release supplies, at minimum, condensed N2/O2 mass fraction, a mass-weighted
particle-size distribution, gas and particle axial velocities, and their
measurement locations/time windows.  Any such coupling must conserve mass,
component mass, axial momentum and total energy through phase transfer and
drag heating; it must then be scored without refitting against both a
near-field profile and an independent downstream trajectory.

Until then it remains an opt-in, uncalibrated building block, not a claim that
DEGALI predicts particle slip for a specific LH2 release.
