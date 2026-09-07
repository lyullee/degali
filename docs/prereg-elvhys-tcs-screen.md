# Pre-registration: ELVHYS WP4.2 transfer-space screen

Date frozen: 2026-09-04, before downloading any per-test time series.

## Question and scope

Can the public ELVHYS WP4.2 data provide an independent falsification screen
for the recommended cryogenic-hydrogen near-field model?

It cannot be treated as a free-jet calibration set.  The release is inside a
1 m3 stainless-steel transfer connection space with a floor, walls, vents and
either active or passive ventilation.  Those effects are absent from the
axisymmetric model.  The screen therefore separates near-nozzle observations
from enclosure-scale observations and never uses the latter to tune an
entrainment or phase-change coefficient.

## Frozen first comparison

Use tests 10 and 11.  The dataset README describes both as nominally identical
horizontal dispersion repeats:

- 2 barg tanker pressure;
- 1.0 mm nozzle;
- 500 L/min active ventilation.

Download only `CONC`, `TEMP`, `PRES`, `FLMT` and `MISC` for these two tests.
The repository metadata lists a combined size of about 55 MB.  No other test
series is inspected until this screen is reported.

The horizontal nozzle is at `(x, y, z) = (20, 500, 250) mm`.  The bottom-axis
concentration/temperature sampling tubes are at `y = 500 mm, z = 200 mm`, so
the model must be sampled at a 50 mm radial offset, not at its centreline.  The
primary near-source stations are:

| channel | axial distance from nozzle | radial offset |
|---|---:|---:|
| `Bottom1_Front` | 0.08 m | 0.05 m |
| `Bottom3` | 0.28 m | 0.05 m |

Stations at 0.48 m and beyond are enclosure observations only because the
finite walls, floor, ventilation and recirculation can no longer plausibly be
screened out.  They may be plotted but will not score the free-jet model.

## Frozen reduction

1. Inspect headers, units, missing values and sample clocks before computing a
   result.  Record any mismatch between the README and files.
2. Identify the longest physical release interval from the nozzle pressure and
   nozzle-temperature channels.  Use a threshold defined from the file units
   and documented after header inspection; do not hand-pick a favourable time
   interval.
3. Remove the first and last 25% of that interval.  The median over the central
   50% is the frozen quasi-steady statistic.  Also retain the 10th--90th
   percentile range.
4. Compute test-10 versus test-11 repeat differences before comparing either
   test with the model.
5. Treat the published instrument errors as minimum point uncertainties:
   +/-1 vol% H2 and +/-1 K.  Sampling-line smoothing and enclosure variability
   are additional, not included in those minima.
6. Run the current recommended dry-air research configuration with the
   measured release boundary.  No coefficient is fitted to these data.  If the
   measured boundary is insufficient to define a source, report that fact and
   stop short of a quantitative prediction.

## Interpretation and adoption rule

- A discrepancy inside the repeat spread or instrument error is unresolved,
  not a model improvement opportunity.
- A consistent near-source residual in both repeats is a candidate model fault
  only if its sign cannot reasonably be explained by the 50 mm probe offset,
  sampling dynamics, active ventilation, floor or wall proximity.
- Enclosure-scale residuals diagnose a missing enclosure/ventilation model;
  they must not alter the free-jet entrainment or condensed-air closures.
- No closure change is adopted from two repeats alone.  Any candidate must be
  independently pre-registered and must not worsen the corrected 369-point
  Hecht--Panda four-slope result.
- Passing this screen is not validation of confined LH2 dispersion.  It only
  shows that the model was not falsified at the two stated near-source probes.

## Provenance warning

The downloaded dataset README says the campaign was conducted in autumn 2024,
whereas `ELE402HSEMETA.csv` gives experimental dates in November 2025 and the
per-test filenames also contain 2025 dates.  Until the publisher resolves this,
report the dataset DOI and test identifiers rather than asserting a campaign
year.

---

## Post-registered result

The ten specified files were downloaded and their official MD5 hashes
verified. All have a 0.05 s median sample interval. The longest physical
intervals selected by the frozen pressure/temperature midpoint rule are
344.20--1083.20 s for Test 10 and 112.30--902.35 s for Test 11.

Central-50% medians show that the nominal repeats do not share the same actual
source: PT2 is 1.9297 versus 1.7728 barg and `NozzleTemp` is 157.64 versus
125.58 K. Bottom-1 H2 is 35.82 versus 37.12 vol%, but Bottom-3 is 36.73
versus 44.53 vol%. Bottom-1 temperature differs by 11.64 K and Bottom-3 by
5.22 K. The fan flows, 508.21 and 507.52 L/min, are close.

The quantitative free-jet prediction was stopped as pre-registered. `FLMT`
contains only `FanFlowMeter`; no measured hydrogen flow appears in the raw
files or final-report P&ID. D4.6 also supersedes neither geometry source:
its Table 1 gives the nozzle at `z = 200 mm`, while the archive sensor sheet
used to freeze this protocol gives `z = 250 mm`. The model sampling radius is
therefore unresolved (zero or 50 mm).

Decision: no closure is adopted and these observations are not scored as
model residuals. The exact missing fields and a ready request are recorded in
[`data-request-elvhys-source.md`](data-request-elvhys-source.md); the complete
screen is in [`elvhys-tcs-audit.md`](elvhys-tcs-audit.md).
