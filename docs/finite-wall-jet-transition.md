# Finite wall-jet transition research closure

`degali.addons.FiniteWallJetCrosswind` adds one bounded attachment state to
the independent-energy LH2 crosswind model. It is opt-in and does not alter
the legacy DEGADIS/JETPLU or the frozen public LH2 assessment.

The wrapper retains the five transported balances (total mass, hydrogen,
horizontal momentum, vertical momentum and energy). Ground contact adds two
effects:

1. surface stress `rho_a u_*^2` removes streamwise momentum over the contact
   chord; and
2. the positive vertical force is reduced by the dynamically attached part of
   the section.

Attachment relaxes towards a local equilibrium over a multiple of the local
above-ground plume depth. The equilibrium releases continuously as

`Ri_lift = g H max(rho_a-rho, 0) / (rho_a u_*^2)`

passes the configured critical Richardson number. The default value 30 is the
same published diagnostic already documented for `JetCoefficients`; it was
not fitted to FFI Test 6. The response-depth multiplier is exposed so that it
can be fixed from independent wall-jet or time-resolved source evidence.

Typical use after an accepted independent-energy interface is:

```python
from degali.addons import FiniteWallJetCrosswind

# Build the interface with ground_interaction="geometry" or "surface_layer".
wall = FiniteWallJetCrosswind(interface.model)
initial = wall.initial_state(interface.state)
result = wall.solve(initial, maximum_distance=100.0)
```

This is a mechanism-isolation tool, not a validated default. The single
Gaussian ellipse still does not resolve wall-normal and lateral wall-jet
profiles independently, substrate cooling is not time dependent, and the
model does not reconstruct wind-direction meander. Test 6 may therefore be
scored only as a pre-registered diagnostic until independent evidence fixes
those missing states.

The fixed Test 6 diagnostic is reported in
`prereg-ffi-test6-finite-wall-jet.md`. It rejected this Richardson-controlled
closure: the plume reached the downstream model with a lift Richardson number
already above 400, so the new state released immediately and did not repair
the 30 m concentration. The module remains available to make that negative
result reproducible; it is not a recommended Test 6 option.
