# Turbulence and particle-slip evidence boundary

Date frozen: 2026-09-17.

## Decision

Finite turbulent-energy transport and gas--condensed-air slip remain
**conservation-verified boundary models**, not validated LH2 field closures.
They are excluded from DEGALI's default prediction and from any quantitative
claim for a particular LH2 release until their identifying measurements are
available.  This is a model-status rule, not a claim that either physical
mechanism is absent.

The executable `TransportEvidenceBoundary` records this rule on the optional
TKE parameter conversions/transport outputs and two-velocity drag results:

```text
numerical_boundary_verified             = true
physical_closure_validated               = false
default_prediction_enabled               = false
quantitative_lh2_prediction_allowed      = false
```

No flag or coefficient in a scenario input can reverse those four fields.

## Turbulence: retained boundary

The finite-TKE branch conserves the exchange `mean kinetic energy -> Q=rho k
-> heat` and can reject an input below a positive-semidefinite Reynolds-stress
lower bound.  It accepts no default `k`, integral scale, dissipation time,
eddy viscosity, TKE transport ratio, normal stress, or pressure--strain
coefficient.  Public LES/RMS figures remain a realizability envelope, not a
gas-phase experimental TKE or epsilon boundary.

To consider a new closure, one matched release needs gas-phase velocity RMS,
Reynolds-stress information (or an explicitly stated incomplete subset), a
dissipation measurement or independently measured length/time scale, and
thermal/species observations at an independent location.  The candidate must
pass the existing total-energy and PSD screens before any field score is
calculated; it still remains non-default until independently replicated.

## Condensed-particle slip: retained boundary

The two-velocity functions exactly conserve axial momentum and return the
resolved kinetic-energy loss.  They do not select an airborne particle
diameter, phase inventory, initial slip, gas mean free path, or a drag law for
the changing cryogenic mixture.

`transported_condensed_air_source` is now explicitly labelled as one of two
kinematic limits:

- `no_slip_retained_condensate_limit`: retained condensate shares the gas
  axial velocity; optional declared diameter affects only the settling bound.
- `stationary_condensate_limit`: the opposite zero-axial-velocity deposition
  bound.

Neither result represents finite particle slip or a calibrated LH2 source.
Moving to finite-slip transport requires a co-located condensed N2/O2 phase
inventory, mass-weighted airborne size distribution, gas and particle
velocities, and gas-property/location/time-window records sufficient to test
the selected drag law.  The resulting model must close component mass,
momentum and total energy and be checked against both near-field and separate
downstream observations without refitting.

## What this freeze does not do

It does not discard the conservation kernels, set the mechanisms to zero, or
reinterpret their numerical tests as field validation.  It preserves them as
auditable bounds and prevents their present inputs from being promoted by
convenience.
