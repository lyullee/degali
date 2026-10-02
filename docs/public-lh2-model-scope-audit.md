# Public LH2 model-scope audit: what a comparison with HyRAM can mean

Date: 2026-09-15

## Evidence examined

Sandia's *Technical Justifications for Liquid Hydrogen Exposure Distances*
(SAND2023-12548, section 3) is a public technical report describing the
HyRAM+ 4.1 calculations used for the then-current LH2 separation-distance
assessment.  It is not a like-for-like performance comparison with DEGALI;
it is used here to distinguish model scope from numerical accuracy.

The report gives three directly relevant statements:

1. its discharge model compares a homogeneous-equilibrium calculation with a
   metastable-liquid calculation and selects the former on the available mass
   flux data;
2. its unignited concentration comparison uses two horizontal DNV releases
   at 30, 50 and 100 m and a small set of PRESLHY readings, with strong wind
   variability and an H2 sensor ceiling at 4 vol-%; and
3. after the orifice, the report's dispersion and flame calculations are
   explicitly gas-only: the released fluid is assumed to retain mass,
   momentum and energy through atmospheric expansion and then to entrain and
   flash to gas, while air/humidity condensation is neglected.

The primary link and an exact provenance note are in
[`references.md`](references.md#read-directly).

## Consequence for DEGALI claims

DEGALI's default remains the validated, fast single-velocity integral path.
Its N2/O2/H2O equilibrium, finite sublimation and two-velocity helpers are
separate opt-in research branches.  They represent mechanisms excluded by the
gas-only scope above, but this **does not** establish that they are more
accurate for a particular release.  They need their own measured phase
inventory, temperature/velocity field and independent concentration test.

Conversely, a numerical agreement between DEGALI's dry default and a HyRAM
result does not validate either model's treatment of cold condensed air; it
only tests the overlapping gas-only assumptions and inputs.

## Practical comparison rule

Any future comparison must state all four items rather than reporting a single
``better than HyRAM`` label:

| Item | Required disclosure |
| --- | --- |
| Release/source | upstream state, discharge relation, area or measured mass flow |
| Near field | whether air/humidity condensation, phase equilibrium, slip and settling are active |
| Dispersion | wind treatment, receptor geometry, averaging window and sensor ceiling/censoring |
| Score | metric and whether the dataset was used to select a branch or reserved for testing |

Under this rule, the public Sandia report is an independent source-model and
far-field comparison context.  It is not a calibration target for the optional
DEGALI mixed-phase paths, because it contains neither phase-resolved airborne
particle measurements nor a matched two-velocity measurement.

## Current action

No coefficient, default setting or validation score changed as a result of
this audit.  It narrows the external claim: DEGALI can describe optional LH2
air-condensation mechanisms absent from the cited gas-only scope, but it must
not claim a quantified accuracy advantage until a phase-resolved test supports
one.
