# Pre-registration: local-shear entrainment after the PRESLHY handoff

Date frozen: 2026-09-05, before calculating this candidate's field result.

## Rationale

The seven-of-seven scalar-peak coupled run improved concentration and centre
height but underpredicted vertical width by 29%. Closing the omitted source
total-mass equation left the surviving subset's width essentially unchanged.
The remaining structural difference from the independently validated
atmospheric model is downstream entrainment.

With `momentum_entrainment_beta > 0`, JETPLU replaces its local perimeter/shear
term with the constant HyRAM source-momentum entrainment. That is appropriate
inside the source-dominated axisymmetric region, but it also disables the
Ricou--Spalding density scaling already supported by the PRESLHY field result.

## Fixed mechanism-isolation candidate

- Use the seven-of-seven scalar-peak near-field solution solely to isolate the
  downstream mechanism. Its known total-mass establishment limitation means
  this run cannot by itself promote an end-to-end model.
- Retain source-momentum entrainment through the axisymmetric `10D` near field.
- At the already fixed handoff, set JETPLU's source-momentum beta to zero and
  resume its local perimeter/shear entrainment with density scaling enabled.
- Keep the existing `alfa1=0.0875`; retain the Houf--Schefer buoyancy term and
  its existing pure-plume cap. Do not blend, fit or add both momentum terms.
- Keep every source, phase, solver, handoff, trial, arc and sensor rule in
  `prereg-preslhy-coupled-crosswind-validation.md` unchanged.

The switch is at a model-regime boundary already chosen without these field
results. It does not change the handoff state or any conserved interface flux.

## Decision

The downstream mechanism is supported only if all seven interfaces pass and,
relative to the scalar-peak/source-momentum result:

1. mean model/measured `sigma_z` is closer to one than 0.710;
2. centre-height MAE does not exceed 0.049 m;
3. concentration VG does not exceed 1.05 times 1.183;
4. FAC2 does not fall below 0.952.

MG and its interval are reported but not used to select the width mechanism.
If supported, the next test must combine local shear with a genuinely
four-flux source boundary and a compatibility-defined handoff; this isolated
scalar run is not a production default.

## Results

Supported as the downstream entrainment mechanism. All seven interfaces
passed. On the same 42 arcs, concentration changed from the scalar/source-
momentum candidate's MG/VG/FAC2 of 0.925/1.183/0.952 to
1.094/1.198/0.952. VG stayed below the frozen 1.242 ceiling and FAC2 did not
fall. On the same 17 vertical fits, centre-height MAE remained 0.049 m and
the model/measured width ratio improved from 0.710 to **0.994**.

This isolates the 29% width deficit to retaining constant source-momentum
entrainment beyond the near-field handoff. The local density-scaled shear
closure is accepted for the next fully conservative combined test, but the
scalar-peak run itself is not promoted because its establishment omits the
total-mass equation.
