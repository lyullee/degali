# Exact data request: ELVHYS Tests 10 and 11

Purpose: make the public Test-10/Test-11 records sufficient for an independent,
no-fit validation of a cryogenic-hydrogen source and near-field dispersion
model.

Dataset DOI: [10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H)

Contact listed in the dataset README: Wayne Rattigan, HSE,
`wayne.rattigan@hse.gov.uk`.

## Minimum fields needed

1. Hydrogen mass-flow rate versus time for Tests 10 and 11, with units,
   calibration date, sampling rate, uncertainty and clock alignment to the
   public files.
2. If mass flow was not measured: the effective flow area or discharge
   coefficient for the 1.0 mm nozzle, the actual bore/tolerance, upstream
   absolute pressure and temperature time series, fluid phase at the pressure
   tap, pipe internal diameter/length/roughness, valve position, and heat-leak
   assumptions needed to reproduce the flow calculation.
3. Confirmation whether `PT1ReleaseRig` and `PT2Nozzle` are gauge or absolute,
   their physical coordinates, tap geometry, response time, calibration
   uncertainty and which signal underlies D4.6 Table 5's mean nozzle pressure.
4. Confirmation of the horizontal nozzle centre coordinates. The archive
   sensor sheet gives `(20, 500, 250) mm`; D4.6 Table 1 gives
   `(20, 500, 200) mm`. Please identify the as-built value for Tests 10 and
   11 and any test-to-test repositioning.
5. Nozzle thermocouple junction type, diameter, mounting/contact method,
   exact coordinate, sampling response/time constant, radiation correction
   if any, and whether `NozzleTemp` represents wall, two-phase stream or gas
   temperature. The central-interval medians are about 158 K and 126 K.
6. Hydrogen-sampling tube internal diameter, length, pump flow, residence
   time, response/deconvolution procedure and uncertainty for `H2Bottom1` and
   `H2Bottom3`.
7. Fan-flow calibration/uncertainty and the exact supply/exhaust boundary
   arrangement for 500 L/min operation, including whether the listed value is
   standard or actual volumetric flow.
8. Confirmation of the experiment dates: the README says autumn 2024, while
   `metadata.csv` and filenames say November 2025.
9. Released-hydrogen ortho/para composition, or the storage history and
   liquefaction specification from which it can be bounded.

## Ready-to-send request

Subject: ELVHYS WP4.2 Tests 10/11 source-boundary clarification

Dear Mr Rattigan,

I am preparing an independent, no-fit validation of a cryogenic-hydrogen
dispersion model using ELVHYS WP4.2 Tests 10 and 11 from dataset
10.18710/JXJP0H. The public concentration, temperature, pressure, fan-flow
and ambient files are readable, but I cannot reconstruct the hydrogen source
mass flux from them.

Could you please provide the measured hydrogen mass-flow time series for
Tests 10 and 11, including units, uncertainty and time alignment? If flow was
not measured, could you instead provide the effective nozzle area or discharge
coefficient and the upstream thermodynamic/pipe information used to infer it?

I would also appreciate clarification of two items that materially affect the
comparison: the archive sensor sheet gives the horizontal nozzle elevation as
250 mm while D4.6 Table 1 gives 200 mm, and the README dates the campaign to
autumn 2024 while the metadata and filenames identify November 2025. Please
confirm the as-built nozzle coordinates and actual test dates. If available,
the PT1/PT2 reference basis and locations, nozzle-thermocouple construction and
response, H2 sampling-line response, and released-hydrogen ortho/para
composition would allow the uncertainty analysis to remain independent of
model fitting.

The aim is to publish all assumptions and compare the two repeats before any
model result is scored. I would be grateful for either the data or a statement
that a requested quantity was not measured.

Kind regards,

[Name / affiliation]
