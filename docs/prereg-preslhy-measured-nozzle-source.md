# PRESLHY measured-pipe source candidate: frozen protocol

Status: pre-registered before the seven-workbook source table or downstream
scores were generated.

## Physical question

The former LH2 source put the liquid on the tanker-pressure saturation curve.
PRESLHY instead reports a cold, pressure-driven liquid with little adiabatic
flash: temperature and mechanical pressure are independent observations.  The
candidate tests whether a source reconstructed at the last instrumented pipe
plane reduces the field error without a fitted source coefficient.

## Frozen population and windows

Only trials 10, 11, 12, 22, 23, 24 and 25 are admitted.  Their zero-based,
end-exclusive steady windows are those already stored in
`reference/preslhy/e35_reduced.json`: respectively 28:77, 46:273, 140:230,
34:113, 26:124, 132:156 and 12:75.  Missing or invalid source data excludes a
trial and counts as a failed source gate; it does not change the population.

## Frozen measurements

For every trial:

1. TC1, TC3, PT1, PT2 and MFM1 are read from `Flexlogger`; Coriolis mass flow,
   drive gain and calculated density are read from `Flowmeter` over the same
   frozen sample indices.
2. The converted TC columns are not used.  Raw cold-junction-corrected
   microvolts are corrected with D3.6 equation 1.  The minimum TC1 voltage in
   the window is fixed at -253 degC through the NIST ITS-90 Type-T polynomial,
   and the same gain correction is applied to TC3.
3. Source temperature is the corrected TC3 median.  Source upstream absolute
   pressure is ambient pressure plus the PT2 median gauge pressure.  PT1 and
   PT2 means, the 5th/95th TC3 percentiles and both mass-flow channels remain
   diagnostics only.

## Frozen pressure-loss mass flow

The primary mass flow is reconstructed from the liquid-filled upstream-line
pressure loss using D3.6 equation C1:

`m = 265 g/s * sqrt(((5.0 barg - median(PT1)) / 2.91 bar) * (rho / 71 kg/m3))`

Here 265 g/s, 2.91 bar and 71 kg/m3 are the report's high-pressure, 12 mm
all-liquid Trial 11 reference.  `rho` is CoolProp hydrogen density at -253
degC and the trial's median PT1 absolute pressure.  Negative pressure loss,
non-liquid state or non-positive reconstructed flow fails the source gate.
The recorded Coriolis mean is retained as the already-tested comparator, not
used to tune this candidate.

## Thermodynamic and momentum construction

The corrected TC3 temperature and median PT2 pressure define a compressed or
subcooled `T,P` state.  Its enthalpy and entropy are evaluated directly, not
from `Q=0`.  A homogeneous-equilibrium isentropic path is searched from PT2
to ambient for the maximum mass flux.  The measured/reconstructed mass flow
sets the discharge coefficient; throat advective momentum plus non-negative
pressure thrust is expanded to an atmospheric plane while conserving total
enthalpy.  That specific momentum and the same `T,P` caloric state feed the
existing condensed-air evaporation and transport ledger.

No discharge coefficient, temperature offset, pressure multiplier,
entrainment coefficient or handoff distance may be fitted to the field data.

## Gates fixed before scoring

The candidate is admissible only if all seven trials have:

- `0 < Cd <= 1`;
- finite atmospheric temperature, density and diameter;
- relative mass, momentum and energy residuals below `1e-8`;
- positive PT2 gauge pressure and corrected TC3 within hydrogen's liquid
  temperature range.

Promotion additionally requires, on the unchanged paired observations:

- concentration MG closer to 1, VG lower and FAC2 no lower than the current
  independent-energy source;
- vertical width ratio closer to 1 and centre-height MAE no larger;
- sourceFlux temperature median absolute error and RMSE both lower than the
  current sourceFlux result.

Failure of any gate leaves this route documented as a rejected physical
candidate.  It does not change the validated default.

## Source-only applicability amendment (before field scoring)

The seven-workbook source extraction showed that the three 25.4 mm open-pipe
trials (10, 22 and 25) have median PT2 from -0.030 to -0.033 barg: the small
negative value is the pressure transducer's zero offset and the outlet is at
ambient pressure.  A `T,P` pair on a two-phase boundary does not determine
quality, density or enthalpy, so applying the compressed-liquid HEM nozzle to
these trials would manufacture a state.  This was identified from source data
before any concentration, geometry or temperature prediction was run.

The frozen implementation is therefore topology-aware:

- all seven trials use the pressure-loss mass flow above;
- trials 11, 12, 23 and 24, with `median(PT2) > 0.05 barg`, use corrected TC3,
  measured PT2 and the HEM/notional-nozzle momentum construction;
- trials 10, 22 and 25 retain the existing open-pipe thermodynamic/momentum
  construction until local quality or outlet density can be independently
  constrained; their corrected TC3 remains a diagnostic;
- the seven-trial field population remains unchanged.  The `Cd` and measured
  pressure gates apply to the four nozzle trials; finite-source and balance
  gates apply to all seven.

This is an applicability split by release hardware and observed outlet
pressure, not an exclusion selected from model agreement.

## Momentum-length amendment (before field scoring)

The source-only interface audit found fully evaporated source velocity ratios
of 5.39 and 6.52 for trials 11 and 23.  Sandia's current HyRAM+ source confirms
that the measured subcritical mass flow sets the throat velocity through
`m/(Cd rho A)` and that the Yuceil--Otugen atmospheric velocity contains the
same `Cd` convention; removing it would be physically and referentially wrong.

The old `velocity / wind >= 10` screen is nevertheless incomplete for this
candidate because it discards source density and area.  Before field scoring,
the measured-source handoff is therefore frozen at

`L_m = sqrt(M_source / (rho_ambient * U_local^2))`

and `L_handoff = min(10 D_source, L_m)`, where
`M_source = rho_source * U_source^2 * A_source`.  This is the length at which
ambient dynamic pressure acting over an area `L_m^2` equals the source axial
momentum flux.  It contains no empirical coefficient and uses only quantities
already conserved at the source.  A case with `L_m < D_source` fails the
axisymmetric near-field applicability gate.  Otherwise it is retained in the
unchanged seven-trial population and transferred to the crosswind equations
at the earlier of 10D and the momentum-balance length.

## Six-constraint profile closure (before field scoring)

At the source-only handoff, trial 11 closed mass, hydrogen, two momentum
components, energy and centre temperature, but its projected H2 half-width
missed by 6.16% against the frozen 5% interface gate.  The five cross-section
unknowns are already consumed by the five flux constraints; a fixed
scalar-to-velocity spreading ratio leaves no degree of freedom for the
independently transported width.

For a measured-source handoff that fails only profile compatibility, the
velocity/scalar spreading ratio is therefore solved as a sixth internal state
over the physically bounded interval `1.0 <= lambda <= 1.5`.  Each trial
value minimizes the signed half-width residual after the five fluxes have
been closed; it is then held fixed throughout that trial's downstream ODE.
No receptor concentration, temperature or geometry observation enters this
solve.  If the resulting fluxes, width or centre temperature still fail their
existing gates, the trial and candidate remain rejected.

## Field result

The full seven-trial candidate closed every interface, but failed every field
promotion direction.  On 33 common concentration arcs, MG changed from
1.0651 to 0.6667, VG from 1.2656 to 1.3824 and FAC2 from 0.9394 to 0.8182.
The vertical width ratio moved from 1.1013 to 0.8378 and centre-height MAE
from 0.0550 to 0.0700 m.  The source reconstruction is therefore rejected as
a production model.  Its source observations and balance-safe `T,P` path are
retained for ablation: pressure-loss mass flow and measured-nozzle
thermodynamics must be separated before either can be judged.

Machine-readable result:
`reference/preslhy/measured_pipe_source_field_2026-09-05.json`.

## Frozen post-rejection ablation

Before running another field prediction, split the rejected full intervention
into exactly two coefficient-free ablations on the same seven trials and the
same gates:

- `flow_only`: pressure-loss mass flow for all seven trials, but the prior
  tank-saturation/opening thermodynamics and momentum construction;
- `nozzle_only`: the prior Coriolis window-mean mass flow for all seven trials,
  with corrected TC3/PT2 thermodynamics and measured-pressure momentum only
  on nozzle trials 11, 12, 23 and 24.

Both retain the momentum-length handoff and six-constraint internal profile
closure when the measured-source table is active.  No mixed weighting,
per-trial selection or interpolation between the two rates is allowed.

## Ablation result

`flow_only` accepted four of seven source interfaces.  Trials 12, 23 and 24
failed the frozen 2 K centre-temperature gate with residuals 2.437, 2.068 and
2.753 K.  `nozzle_only` accepted six of seven; trial 23's coefficient-free
momentum length ended at 4.79 source diameters, before the axisymmetric
self-similar development zone.  In accordance with the protocol, neither
ablation was field-scored and neither may be rescued by per-trial mixing.

Machine-readable results:
`reference/preslhy/measured_pipe_source_flow_only_interface_2026-09-05.json`
and
`reference/preslhy/measured_pipe_source_nozzle_only_interface_2026-09-05.json`.

The combined source's field overprediction together with both ablations now
points to a missing state rather than a missing source coefficient: residual
liquid hydrogen was converted to vapour at the former evaporation endpoint.
The next candidate is frozen separately in
`docs/prereg-preslhy-finite-rate-lh2-droplets.md`.
