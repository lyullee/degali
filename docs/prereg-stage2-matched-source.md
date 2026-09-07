# Stage2: matched-source thermal-profile comparison

2026-09-06. Written before the new source interfaces or field scores.

The Stage1 operational BASE and research CANDIDATE use different mass flows
and nozzle caloric states. Their temperature-error difference cannot isolate
thermal transport. The recorded pressure-loss flow follows PRESLHY D3.6 C1;
the Coriolis channel is not a new ground truth (several drive-gain readings
reach 100%). No field score will be used to select a preferred flow channel.

Use a two-by-two comparison: pressure-loss versus recorded Coriolis window
mean flow, and Gaussian bulk-density versus Gaussian volumetric-enthalpy
profiles. Both profiles at each flow use the same TC3/PT2 nozzle inputs,
open-pipe prior, HEM/phase-energy source construction, wind, humidity,
normal-hydrogen properties, consistent phase ambient, entrainment law,
source momentum-balance handoff rule, and energy transport. Their handoff
states need not be equal: projection onto a different profile is part of the
explicit model change. The velocity/scalar width parameter is held at the
same existing value (no compatibility fitting). A failed interface remains
failed; it does not trigger fitting or a relaxed acceptance threshold.

The matched Coriolis document is an analysis input, never a replacement
measurement table. Its historical pressure-loss loader slot is explicitly
labelled with the chosen flow; original values and hashes are retained.
Stored HEM diagnostic values are labelled historical, because source states
must be recalculated from the changed flow. The input file must not be
mistaken for a new pressure-loss estimate.

Preparation exposed a representation difference before any new model run:
the sealed reduced table uses three decimal places in g/s, whereas the source
table retains raw-window precision. Require identical windows and exact
agreement after rounding the raw mean to three decimals; use the sealed
reduced value to match BASE. Record the rounding difference (all seven less
than 0.0005g/s). This is not a changed physical tolerance or fitted flow.

First verify all seven interfaces for each new configuration. A field run
requires its own passed seven-case checkpoint, identical configuration and
source-input hash. Start with trials10/23 (all41 frozen thermocouples); then
extend admissible paths to the remaining five cases. Use maximum step .02m;
refine promising results to .01m before an accuracy claim.

Score the exact Stage1 concentration38/vertical17/temperature41 keys. Missing
coverage, invalid concentration or a failed source is an explicit failure,
not a reason to reduce the scored population. Per-trial pilot scores are
labelled subsets. Retain temperature minimum, p05 and median statistics
separately. These steady calculations do not predict thermocouple response
or turbulent intermittency. Report per-trial and spatial residuals as well
as aggregates. Balance residuals must remain <=1e-5.

This first calculation is mechanism isolation, not a newly validated
transport closure or an automatic promotion. It determines which errors
persist at a matched source and where the next physical correction is
warranted. Existing pressure-loss results can be reused only when the exact
executed settings, phase model, width parameter, state/flux replay, and
observation keys match. No Stage1 frozen file is changed.
