# Cryo-compressed hydrogen blowdown boundary

## Purpose and scope

`degali.addons.cryogenic_blowdown` is an opt-in, zero-dimensional source-term
model for a **single-component, well-mixed, cryo-compressed hydrogen tank**.
It is not a replacement for DEGALI's atmospheric dispersion solver, and no
default DEGALI calculation invokes it.

The implementation follows the mass and tank internal-energy balance used in
Cirrone et al. (2023):

\[
\frac{d m}{d t}=-\dot m, \qquad
\frac{d(m u)}{d t}=\dot Q-h\dot m.
\]

Hydrogen properties come from CoolProp's Helmholtz-energy EOS.  The mass flow
is a one-dimensional critical-flow value multiplied by an explicit,
user-supplied `Cd`.  No regression or automatic coefficient selection is
performed.

For a single-phase, positive-flow state,
`blowdown_state_to_lh2_source()` can produce a corresponding ambient-pressure
**gas** source plane for the existing conserved LH2 near-field path.  It
expands the reconstructed critical throat while preserving mass, axial
momentum including pressure thrust, and total specific energy.  Here `Cd` is
made explicit as an effective throat area, not applied twice to the flux.
The adapter is a sequence of quasi-steady snapshots; it is **not** a transient
atmospheric transport solver.  It rejects terminal states and also rejects an
ambient-pressure two-phase expansion rather than presenting it as a gas jet.

For the latter case, the separate opt-in
`blowdown_state_to_flashing_droplet_source()` exposes the pressure-thrust
flash plane, including retained liquid and vapour flows.  It is deliberately
not a gas jet.  `blowdown_state_to_homogeneous_evaporation_source()` then uses
that same flash plane as the input to the existing common-velocity,
heat-limited complete-evaporation source.  The latter is a transparent fast
equilibrium bound; it is not a finite droplet-lifetime or gas--liquid slip
prediction.  All adapters reject a two-phase **tank** state because the
withdrawal law remains an independent source assumption.

The homogeneous-evaporation adapter defaults to 295 K specifically to match
the existing LH2 near-field solver default.  If another ambient temperature
is selected, it must be supplied consistently to both the adapter and the
near-field calculation; the source enthalpy reference is not transferable
between different ambient temperatures.

## Deliberate limitations

The published model resolves transient conduction across the tank and pipe
walls.  DEGALI now has optional one-dimensional, finite-volume thermal stores
for both: `TransientTankWall` and `TransientPipeWall`.  They use explicit
inner/outer heat-transfer coefficients, material properties, initial wall
temperature, and external boundary temperature.  The thin pipe wall is
represented by a planar equivalent with its exact cylindrical wall volume;
the flow-to-wall boundary is evaluated with the timestep-average bulk-flow
temperature.  This is a stated reduced-order approximation, not a substitute
for pipe CFD.

The older prescribed, lumped tank-wall `UA` and prescribed pipe-wall `UA`
paths remain available for simple sensitivity studies.  Each is mutually
exclusive with its corresponding transient wall, so that the same wall heat
capacity cannot be counted twice.  The pipe energy balance is solved before
the critical-nozzle calculation.  No wall parameter or temperature offset is
chosen implicitly by the package.

The E3.1 welded T4 release-line thermocouple is now read as a diagnostic when
present.  It is not automatically used as the prescribed pipe-wall
temperature: it senses the flow, whereas the model boundary is the wall.
This distinction is material: a separate published CFD study reports roughly
50 K bulk warming for a 200 bar, 4 mm cryogenic release through its own
ambient-exposed pipe, but that result depends on pipe geometry and wall
boundary.  It is evidence for an explicit sensitivity, not for inserting a
universal 50 K offset into the E3.1 source.

The decision not to infer an internal heat-transfer coefficient from this
signal or from a generic room-temperature pipe correlation is documented in
[`cryocompressed-heat-transfer-closure-audit.md`](cryocompressed-heat-transfer-closure-audit.md).

It also stops, rather than guessing, when the `(density, internal-energy)`
state reaches a two-phase tank condition.  Continuing requires an explicitly
chosen withdrawal phase, tank stratification model, and phase-equilibrium
outflow law.  A gas-only, LH2-only, or clipped continuation would not be a
conservative physical substitute.

An additional `two_phase_withdrawal="vapour"` research mode is available only
when the cylindrical tank diameter and outlet height are supplied.  It uses
equilibrium saturated-vapour critical flow and stops when the calculated
liquid level reaches the outlet.  This is not a claim that the real tank is
well mixed; it makes the alternate assumption visible and prevents an
unannounced phase switch.

`two_phase_withdrawal="homogeneous"` is a distinct research mode.  It
maximises the local-equilibrium, common-velocity two-phase mass flux along the
tank-mixture isentrope.  The throat state establishes that flux, while the
tank energy balance removes the mixture enthalpy at the tank boundary; using
the post-expansion throat enthalpy there would count nozzle kinetic-energy
conversion as tank cooling.  It does not claim phase separation, slip,
flashing kinetics, or a particular outlet phase.  It is therefore useful as
an explicit HEM envelope beside the geometry-limited vapour-withdrawal case,
never as a replacement default.  Neither two-phase path currently applies a
discharge-pipe thermal boundary.

For the HEM case only,
`hem_blowdown_state_to_flashing_droplet_source()` reconstructs the critical
two-phase throat and maps it through the pressure-thrust flash plane.  It
preserves the declared effective source area and makes the post-flash liquid
inventory available to the finite-evaporation research path.  It refuses a
vapour-withdrawal state because a separated outlet is not an HEM source.
`hem_blowdown_state_to_homogeneous_evaporation_source()` is the separate
fast complete-evaporation continuation of that same HEM flash plane; it does
not add a finite-rate evaporation or slip claim.

## First public-data check

For PRESLHY E3.1 part A, 200 bar, 4 mm, nominal 80 K, the adiabatic
single-phase calculation with `Cd=0.7`, 2.815 L and a 0.02 s step reaches the
two-phase boundary after about 1.06 s.  At that boundary it has released
about 0.0844 kg and predicts approximately 11.57 bar and 32.37 K.

This is not an accuracy claim: `Cd=0.7` is within the independently published
0.6--0.8 range, but this particular diagnostic omits the pipe-wall heat
transfer and ends before the experiment's low-pressure, mixed-phase part.
It is nevertheless a useful result: it establishes that the remaining source
discrepancy is a specific **tank phase/withdrawal** problem, not a
justification for changing the atmospheric entrainment closure.

Using the documented 160 mm tank diameter and 30 mm outlet elevation, the
explicit vapour-withdrawal continuation reaches 12.36 bar at 1.0 s and 6.63
bar at 2.0 s (with 28.84 K and equilibrium quality 0.631 at 2.0 s).  The raw
PRESLHY pressure record at the corresponding relay times is approximately
13.16 and 5.55 bar.  These two points are a transparent **non-fitted
diagnostic**, not a validation score: the calculation omits discharge-pipe
heat transfer, wall conduction, and tank stratification, and the physical
PNoz response begins about 0.065 s after the relay.

## Seven-condition public pressure screen

The complete public D=4 mm, nominal-80 K E3.1 pressure set (approximately
5--200 bar) has now been passed through the same non-fitted source envelope:
77/80/84 K, `Cd=0.6/0.7/0.8`, and explicit vapour-withdrawal/HEM end members.
Physical source start is each file's first 2-bar PNoz rise, rather than the
electrical relay. The result is recorded in
[`e31-d4-80k-blowdown-pressure-results.md`](e31-d4-80k-blowdown-pressure-results.md).

The envelope contains the 50 and 200 bar pressure records at 0.5, 1 and 2 s,
but misses some late low/intermediate-pressure stations. This is an explicit
model boundary, not a reason to select a favourable `Cd` or introduce an
unmeasured heat-transfer coefficient. The cryo-compressed source remains
opt-in and no atmospheric-dispersion default changes.

## Evidence

Cirrone, D., Makarov, D., Kashkarov, S., Friedrich, A., & Molkov, V. (2023).
*Physical model of non-adiabatic blowdown of cryo-compressed hydrogen storage
tanks*. International Journal of Hydrogen Energy, 48(90), 35387--35406.
https://doi.org/10.1016/j.ijhydene.2023.05.182

The open published version describes the non-ideal EOS, tank and pipe heat
transfer treatment, and validation against sixteen PRESLHY tests spanning
80--310 K, 0.6--20 MPa, and 0.5--4 mm nozzles.  It reports test-specific
optimum discharge coefficients in the range 0.6--0.8; this range is evidence
for a sensitivity envelope, not a package default.
