# ELVHYS WP4.2 public validation boundary

The HSE ELVHYS WP4.2 record, DOI
[10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H), is a public CC0 dataset
for 48 cryogenic-hydrogen tests in a 1 m³ transfer-connection space. It is
an unusually useful independent data source: individual tests have 20 Hz
hydrogen concentration, temperature, release/nozzle pressure, ventilation
flow and weather streams. The publisher documents sensor coordinates, stated
instrument errors, and that hydrogen channels were aligned to remove their
sampling-line delay.

`degali.validation.read_elvhys_test` reads a user-held copy of a declared
test number, requiring the official `ELE402HSE[No.]...` naming and an explicit
`Time` field; nested archive folders are accepted. It preserves missing readings as `NaN`; it does not silently
subtract a baseline, align clocks, infer a release window or estimate a source
flow. The public data remain outside the package.

## What the data can test

- future confined-cryogenic-H2 observation operators: co-located temperature
  and concentration history, enclosure pressure, vent/ventilation response;
- future obstacle/wall models in the stated 1 m³ geometry, when their source
  boundary and observation operator are declared beforehand;
- repeatability and timing of a pressure/temperature/concentration chain.

## What it cannot currently score

`FLMT` records `FanFlowMeter`, i.e. **ventilation-air** flow. The data package
does not supply independently measured hydrogen mass flow. Further, the
experiment contains walls, floor impingement, finite-volume accumulation and
vent exchange. It is consequently not a calibration or quantitative score for
DEGALI's outdoor LH2-pool plume model, and the reader makes that limitation
explicit with `quantitative_outdoor_lh2_pool_validation_allowed == False`.

For this evidence to promote an enclosure branch, a source-boundary model and
a predeclared confined-space conservation/observation calculation must first
be independently tested. No wall, turbulence or source coefficient is fitted
from these data.
