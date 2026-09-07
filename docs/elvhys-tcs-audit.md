# ELVHYS TCS data audit

Date: 2026-09-04

## Outcome

The ELVHYS WP4.2 archive is useful for testing an eventual confined-release
model, but Tests 10 and 11 cannot yet be used as a quantitative validation of
the present free-jet model. The public files contain hydrogen concentration,
temperature, pressure, fan-flow and ambient channels, but no hydrogen mass
flow. Two official documents also disagree on whether the horizontal nozzle
is at `z = 200 mm` or `z = 250 mm`. Both quantities are first-order boundary
conditions; guessing either would convert an independent validation into an
unrecorded calibration.

No model coefficient was changed from this screen.

## Sources and integrity

- Dataset: ELVHYS WP4.2 HSE test data, DOI
  [10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H), CC0.
- Final report: Rattigan, Vizma and Welch, *D4.6 release into cold room TCS
  tests*, available from the [CORDIS ELVHYS results
  page](https://cordis.europa.eu/project/id/101101381/results).
- The cached final PDF is 35,049,436 bytes with SHA-256
  `e4efbf4d5c0e6f002588f7a2706461f555843ff1a0ea3cf3d989ea7b9bd891c0`.

Only metadata and the ten pre-registered Test-10/Test-11 streams were
downloaded. Their official MD5 values were verified:

| test | stream | MD5 |
|---:|---|---|
| 10 | CONC | `bddfc764c4793f82a9b27ca2937ca912` |
| 10 | FLMT | `10485ca53bf341a47c9afdf5ec0abcbc` |
| 10 | MISC | `af69c541a4c125abd674308354c5c2fc` |
| 10 | PRES | `1155d03dad3097badfabee97a36499af` |
| 10 | TEMP | `ebb7629d9d5c2fc88755aa7062ff3506` |
| 11 | CONC | `904756b5284fc9a3344f4ce3a2edf32e` |
| 11 | FLMT | `91e5d4b522ae75a4923583199b538c54` |
| 11 | MISC | `14d998ea9533c03f0164e6e99bce2c2a` |
| 11 | PRES | `a3db4d373885667f7d1a254129d2c758` |
| 11 | TEMP | `438034654dd079e6d440bc1501879b81` |

The archive remains under `tmp/elvhys/`; it is third-party source material
and is not vendored into the package.

## Frozen reduction and result

The protocol was fixed in
[`prereg-elvhys-tcs-screen.md`](prereg-elvhys-tcs-screen.md) before the raw
time series were read. All ten files use a 0.05 s median sample interval.
The release interval was selected by a data-driven conjunction of nozzle
pressure above its baseline/99th-percentile midpoint and nozzle temperature
below its baseline/1st-percentile midpoint. The statistic below is the median
over the central 50% of the longest detected interval.

| quantity | Test 10 | Test 11 | repeat difference |
|---|---:|---:|---:|
| detected interval | 344.20--1083.20 s | 112.30--902.35 s | - |
| PT2 nozzle pressure | 1.9297 barg | 1.7728 barg | 0.1569 barg |
| nozzle thermocouple | 157.64 K | 125.58 K | 32.06 K |
| fan flow | 508.21 L/min | 507.52 L/min | 0.69 L/min |
| H2 Bottom 1 | 35.82 vol% | 37.12 vol% | 1.30 percentage points |
| H2 Bottom 3 | 36.73 vol% | 44.53 vol% | 7.80 percentage points |
| T Bottom 1 | 259.27 K | 247.63 K | 11.64 K |
| T Bottom 3 | 274.48 K | 269.26 K | 5.22 K |

The final report's peak-enclosure summaries (18.0% and 19.2%) are reasonably
repeatable, but that statistic is not the local jet-axis concentration. The
raw Bottom-3 medians differ by 7.80 percentage points while source pressure
and measured nozzle temperature also differ. This confirms that nominally
identical test labels do not isolate model error at a local probe.

## Geometry and source ambiguities

1. The archive sensor sheet places the horizontal inlet at
   `(20, 500, 250) mm`. D4.6 Table 1 places it at `(20, 500, 200) mm`, on the
   same elevation as Bottom 1--9. The resulting radial offset is either 50 mm
   or zero. For a 1 mm jet at the 80 mm first station this is not a small
   perturbation.
2. `FLMT` contains only `FanFlowMeter`; the P&ID has no hydrogen mass-flow
   measurement. D4.6 explicitly says the 25.4 mm pressure-peaking tests used
   an estimated 180--200 g/s from previous work. That estimate does not
   define the 1.0 mm dispersion-test flow.
3. Test 10 and 11 use different actual nozzle pressures and very different
   thermocouple histories. The report itself notes that nozzle temperature
   affects mass flow.
4. The README calls the campaign autumn 2024, while `metadata.csv` gives
   5--28 November 2025 and the per-test files are dated 2025. This is a
   provenance issue, not a physics input, but should be resolved before a
   paper is submitted.

## Model implication

These tests require a source model plus wall/floor impingement, finite-volume
accumulation, active ventilation, buoyancy-driven vent exchange and sampling
line response. Treating Bottom 5 and beyond as a free jet would fold those
effects into an entrainment coefficient. The defensible next step is to
obtain the missing source boundary and coordinate clarification, then use
Bottom 1 and Bottom 3 only as a near-source falsification screen before adding
an enclosure control-volume model.
