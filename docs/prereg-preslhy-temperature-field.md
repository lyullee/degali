# Pre-registration: PRESLHY near-field temperature field

Date: 2026-09-05

## Question

The trial-10 vertical-momentum audit located the remaining centre-height error
in the thermal/density history: the independent-energy branch warms earlier,
becomes positively buoyant sooner and accumulates about twice the vertical
buoyancy of the corrected JETPLU baseline.  The original reduction used only
hydrogen concentration.  The newly retrieved original PRESLHY workbooks
contain collocated type-T thermocouples, so this audit asks whether the
candidate is in fact too warm at the measured near-field points.

## Frozen data and reduction

- Use original PRESLHY workbooks for trials 10 and 23 only.  They are the two
  high-pressure horizontal trials already selected for detailed source and
  trajectory diagnosis and were retrieved before this protocol was written.
- Use sheet `Flexlogger`; select only converted near-field temperatures whose
  header ends in `TC<serial>C`.  Do not use raw thermovoltages, pipework,
  far-field stand or ambient channels.
- Map thermocouple serials and positions exactly from D3.6 Table A3.  Trial 10
  uses the 0.5 m array configuration.  For trial 23 add 1.0 m to Table A3's
  listed vertical coordinates because the array moved with the 1.5 m release.
- Freeze the release interval independently from temperature using the
  existing concentration reduction: five-point centred mass-flow smoothing,
  then the longest consecutive run above 50% of the smoothed peak.  Apply the
  same integer row interval to the 1 Hz `Flexlogger` sheet.
- The primary measured statistic is the minimum converted temperature in that
  interval, matching PRESLHY D3.6 Figure 24.  Also report the 5th percentile,
  median and number of finite samples as robustness diagnostics; none may
  replace the primary statistic after results are seen.
- Exclude a channel only if it is absent, has no finite value in the frozen
  interval or falls outside the type-T conversion range.  Do not exclude a
  cold point merely because it looks anomalous.

## Frozen model comparison

- Compare the corrected JETPLU baseline and the coefficient-free
  `source_flux` independent-energy branch, with their already frozen source,
  Li enthalpy transport, Houf velocity-width mapping and 0.02 m maximum step.
- Interpolate each steady model at the thermocouple's downwind distance and
  evaluate temperature at the **sensor height**, not at the model centre.
  This is essential because centre height is itself the disputed quantity.
- Score only downwind distances at or beyond the branch's common handoff.
  Upstream thermocouples remain source-development diagnostics and must not be
  folded into a downstream score.
- For each eligible sensor report predicted minus measured temperature in K.
  Aggregate paired errors by median absolute error and signed median bias,
  first by trial and then across the two trials.  Bootstrap uncertainty, if
  shown, resamples whole downwind sections rather than individual sensors.
- No coefficient, source state, release interval, sensor selection, trajectory
  or profile width may be fitted to these temperatures.

## Interpretation fixed before calculation

- If the source-flux branch has a positive median temperature bias at least
  10 K and larger than the baseline, its excessive buoyancy is supported as a
  thermal-closure error; the next implementation target is finite-rate
  air-condensation/heat transfer, not added drag.
- If both branches are at least 10 K colder than measurement, missing heating
  remains plausible but cannot explain why the source-flux branch rises more.
- If the source-flux branch is within 10 K while its centre trajectory still
  rises too much, temperature alone is not the cause; density composition,
  phase loading or the observation operator must be audited next.
- Trial-10 blockage of the first four concentration sampling lines does not
  remove their collocated thermocouples.  Concentration and temperature are
  not paired in time because the sampling lines have an approximately 18 s
  delay; only their separately frozen extrema may be compared descriptively.
- The report notes thermocouple uncertainty and outdoor variability.  This is
  a diagnostic mechanism test, not an independent validation campaign and not
  grounds by itself to promote either branch.

## Results

The raw workbooks supplied 24 near-field temperature channels per trial.  At
and beyond the common handoff, 42 sensor/model pairs remained.  With the
frozen mean mass-flow input and the pre-registered minimum statistic, the
corrected JETPLU baseline gave median bias +7.64 K and median absolute error
28.31 K.  The `source_flux` independent-energy branch gave median bias
+20.31 K and median absolute error 20.31 K.  Thus the candidate is too warm
against the report's cold-core statistic even though its absolute-error
ranking is mixed.

The effect is campaign-dependent and cannot be represented by one fitted
offset.  Trial 10's candidate median bias is +13.52 K; Trial 23's is +37.55 K.
On the Trial 23 centreline it is +75.34 K at 1.19 m and +85.52 K at 1.78 m,
then remains +37.55 K at 6 m.  The pre-registered robustness calculation also
shows an observation-operator limitation: against the temporal median, the
all-sensor candidate bias falls to -1.63 K and median absolute error to
14.11 K.  That apparent agreement is driven by off-centre near-ambient
sensors.  On the 12 centreline points its median absolute error remains
32.75 K, and the Trial-23 centreline median bias remains +37.26 K.  The cold
core therefore still disappears too early, but the 42-point signed bias must
not be read as a steady-mean validation statistic.  The result rejects added
drag as the primary repair and points to source unsteadiness, premature
single-temperature phase equilibration and/or premature turbulent mixing.
No coefficient was fitted and neither branch is promoted by this diagnostic.

During the audit, Table A3 parsing was also corrected to accept the seven-
character thermocouple serials and Unicode minus variants.  This restores two
ground sensors and four negative-lateral positions.  The negative-lateral
concentration channels stay below the existing signal floor, while the added
ground sensors change only the fitted vertical reductions; this is recorded
separately from the model comparison.
