# Pre-registration: conserved near-field to crosswind handoff

Date frozen: 2026-09-04

## Question

Can the accepted conserved Gaussian LH2 near-field solution be handed to the
existing JETPLU crosswind trajectory without silently creating or destroying
mass, hydrogen, momentum, energy, or plume width at the interface?

This is an interface verification, not a fit to PRESLHY, FLADIS, or any other
downwind validation data. No empirical coefficient may be adjusted using the
outcome of this check.

## Fixed mapping

The handoff station is an explicit streamline distance `S`; if it is omitted,
the final near-field station is used. The station must be inside the computed
solution and the near-field run must pass its frozen numerical-conservation
screen.

The accepted near-field profiles are

`v = vc exp[-r^2/B^2]`

and

`rho*Y = rhoc*Yc exp[-r^2/(lambda*B)^2]`,

with `lambda = 1.16`. JETPLU uses a Gaussian standard deviation `sigma` for
contaminant concentration and a turbulent Schmidt number `Sc` for excess
velocity. Their non-fitted, shape-equivalent conversion is therefore

`Sc = lambda^2` and `sigma = lambda*B/sqrt(2)`.

The JETPLU mixing line is rebuilt from the near-field centre state at the
handoff: hydrogen mass fraction, dry-air/water ratio, temperature, and
specific enthalpy. This makes the handoff state the concentrated end of the
new dilution line rather than incorrectly treating it as pure hydrogen.

JETPLU's four free similarity variables (`cc`, `sigma_y*sigma_z`, `theta`,
and excess centre velocity) are projected against four integral targets:

1. total mass flux,
2. hydrogen mass flux,
3. horizontal momentum flux, and
4. vertical momentum flux.

Position is copied from the same near-field station. The projected state is
not sent downstream unless every acceptance condition below passes.

The source-momentum and Houf--Schefer buoyancy entrainment constants are
copied from the near-field model. This prevents an unphysical closure jump at
the interface: the crosswind equations add trajectory and atmospheric-shape
terms, but retain the same source-defined entrainment contribution.

## Frozen acceptance conditions

1. Near-field conservation screen: pass.
2. `Sc` equals `lambda^2` to a relative tolerance of `1e-10`.
3. Relative residual of each projected flux (total mass, hydrogen, horizontal
   momentum, vertical momentum): at most `1e-8`, with a `1 N` absolute scale
   for a target momentum component that is zero.
4. Hydrogen mass-density half-width mismatch: at most 5%.
5. Total-energy-flux mismatch: at most 2%. The comparison includes mixture
   enthalpy relative to ambient and absolute kinetic energy, with the same
   ambient stream subtracted from both sides.
6. Centre temperature mismatch: at most 2 K.
7. All state variables are finite and physical; the handoff concentration
   must lie inside its rebuilt mixing table.

Failure of any item makes the handoff `accepted = false`; downstream
integration must then refuse to run. A failed energy or width screen is
evidence that JETPLU needs an additional transported energy/profile state,
not permission to retune entrainment.

## Initial numerical case

The first test uses a measured-throat, ambient-pressure Hecht--Panda source
and a horizontal orientation. Crosswind is applied only in JETPLU after the
axisymmetric near-field station. This separates interface conservation from
the later experimental-validation question. PRESLHY is not used until the
interface passes unchanged.

## Interpretation

- Pass: the existing JETPLU state is sufficient at this interface and may be
  connected for independent downwind validation.
- Four balance equations pass but energy/width fails: add an explicit
  thermodynamic or profile degree of freedom before validation.
- Projection equations fail: the two similarity representations are not
  compatible under crosswind and require a larger state, not calibration.
