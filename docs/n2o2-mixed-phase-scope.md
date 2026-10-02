# N2/O2 condensed-phase thermodynamic scope

Date: 2026-09-17

## Decision

The existing cryogenic-air calculation can still solve its declared
component-wise mass and partial-pressure balance, but it must not be described
as a complete N2/O2 condensed-mixture EOS. `N2O2CondensedPhaseScope` now
classifies every `AirPhaseEquilibrium` state without assigning a new phase
fraction, activity coefficient, latent heat, or fit.

| Result | Interpretation | Permitted claim |
| --- | --- | --- |
| `uncondensed_bulk_air` | Neither bulk N2 nor O2 is condensed. | The component gas bookkeeping is within this particular mixture-phase scope. |
| `single_component_condensate_bound` | The separate-pure-component calculation condenses only one bulk air component. | A transparent limiting bound only; dissolved-component chemical potentials are absent. |
| `co_condensed_n2o2_mixture_unsupported` | Both N2 and O2 are condensed at or above the N2 triple temperature. | No common liquid/solid solution equilibrium or caloric claim. |
| `co_condensed_subtriple_n2o2_mixture_unsupported` | Both are condensed below the N2 triple temperature. | No mixed solid/liquid stability, composition or enthalpy claim. |

The reported condensed N2 mole fraction is an inventory ratio of the two
separately calculated condensed amounts. It is deliberately not interpreted as
the composition of any equilibrium liquid or solid phase.

## Consequence for the LH2 near-field research path

For dry air, `run_lh2_near_field_research` now scans its solved centreline
against this scope and records the classifications in its result notes. If a
bulk-air condensate occurs it emits an explicit warning: the N2/O2 treatment is
a separate-pure-component research bound, not a complete mixture closure.

For humid air, no equivalent automatic classification is reported. Water
changes the gas inventory, while a common N2/O2/H2O condensed-phase model is
also absent. The run therefore reports that the dry-air scope audit was not
performed instead of silently reusing it.

This does not change the fast gas-only `assess` route, alter any entrainment
or phase-change coefficient, or turn a research-bound result into a failed
numerical integration. It narrows the physical claim attached to a successful
conservation calculation.

## What remains required before a mixed-phase prediction

1. A common N2/O2 liquid and mixed-solid chemical-potential model over the
   required temperature/composition/pressure range.
2. Enthalpy, heat capacity and phase volumes derived from that same reference,
   including solid--liquid transitions.
3. Independent phase-inventory and particle residence/slip evidence before
   connecting the closure to a dispersion accuracy claim.

The 59--63 K local states identified from the controlled PRESLHY reduction
remain an explicit unsupported boundary. The public phase diagrams establish
that independent pure condensates are inadequate there, but do not provide the
caloric data required to create a replacement closure.
