# ELVHYS WP4.2 conditional dynamic evidence

Date: 2026-10-09

## Purpose and decision

The local evidence archive contains the public ELVHYS WP4.2 HSE dataset
(DOI [10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H)).  It is the strongest
additional local candidate for a dynamic and measured-geometry evidence lane:
the archive contains time-resolved hydrogen concentration, temperature,
pressure, ventilation and weather channels, instrument specifications and a
three-dimensional transfer-connection-space sensor/nozzle table.

This record promotes the dataset only as **conditional confined dynamic and
obstacle-context evidence**.  It is not added to the E3.5 or FFI/DNV headline
scores and it does not qualify an outdoor free-jet, pool or site-specific
obstacle-wake claim.  The reason is material: the published `FLMT` channel is
fan/ventilation flow, not a measured hydrogen mass-flow history.  A pressure
history alone is retained as a source diagnostic and is not silently converted
to a mass-flow boundary.

## Dataset boundary

The README describes 48 tests in a nominal 1 m3 stainless-steel transfer
connection space.  The campaign includes small horizontal and vertical leaks
(0.2--1.8 mm), 25.4 mm pressure-peaking tests and ignited tests, with active
or passive ventilation.  Most channels are sampled at 20 Hz; Test 30 is 100 Hz,
Test 31 is 500 Hz and ignited blast-pressure channels are 100 kHz.

The metadata record gives experimental dates 05--28 November 2025 at HSE
Science & Research Centre, Buxton, UK.  The dataset README instead describes
the campaign as autumn 2024, while the per-test filenames contain 2025 dates.
This discrepancy is retained as a provenance warning.  The manuscript cites
the DOI and test identifiers rather than asserting a campaign year.

## Measurement and geometry evidence

The archive README reports the following declared specifications:

| Channel | Declared range and specification | Treatment |
|---|---|---|
| H2 concentration | 0--100 vol% H2, Xensor XEN-5320-ALU-USB, +/-1 vol% H2 | specification metadata only; no event-level certificate promoted |
| Type-T temperature | 77--300 K, +/-1 K; five-point thermocouple calibration | specification/calibration description retained |
| Release/nozzle pressure | 0--16 barg, +/-0.08 bar | source diagnostic, not H2 mass flow |
| TCS pressure | 0--100 mbarg, +/-1.8 mbar | enclosure diagnostic |
| Ventilation flow | KROHNE flowmeter, declared fan-flow channel | not a source-rate measurement |
| Weather | wind, ambient temperature, humidity and pressure fields | contextual; weather-station indication is not used to qualify outdoor transport |

The coordinate table uses millimetres referenced to the transfer-space frame.
The horizontal nozzle is `(20, 500, 250)` and the vertical nozzle is
`(93, 500, 250)`.  Bottom-axis concentration/temperature stations are at
`z = 200 mm`, with `Bottom1_Front` at `(100, 500, 200)` and `Bottom3` at
`(300, 500, 200)`.  The 50 mm nozzle-to-probe vertical offset is therefore
explicitly retained.  The final report has a conflicting horizontal-nozzle
entry at `z = 200 mm`; the archive geometry is used for this boundary record,
and the conflict is not hidden.

## Pre-registered Test 10/11 screen

The repository protocol fixed Tests 10 and 11 before reading their raw
streams.  Both are nominally 2 barg, 1.0 mm horizontal-release tests with
500 L/min active ventilation.  The pressure/temperature rule selected the
longest physical intervals shown below; the statistic is the median over the
central 50% of each interval.

| Quantity | Test 10 | Test 11 | Difference |
|---|---:|---:|---:|
| Selected interval (s) | 344.20--1083.20 | 112.30--902.35 | -- |
| PT2 nozzle pressure (barg) | 1.9297 | 1.7728 | 0.1569 |
| Nozzle temperature (K) | 157.64 | 125.58 | 32.06 |
| Fan flow (L/min) | 508.21 | 507.52 | 0.69 |
| H2 Bottom 1 (vol%) | 35.82 | 37.12 | 1.30 percentage points |
| H2 Bottom 3 (vol%) | 36.73 | 44.53 | 7.80 percentage points |
| Bottom 1 temperature (K) | 259.27 | 247.63 | 11.64 |
| Bottom 3 temperature (K) | 274.48 | 269.26 | 5.22 |

These nominal repeats have materially different measured source conditions.
The result is therefore a useful dynamic/repeatability warning, not a clean
model-error score.  The existing ELVHYS audit retains the frozen decision:
no closure is fitted and no quantitative outdoor prediction is produced.

## Integrity record

The public archive is controlled outside this repository and is not copied to
the public snapshot.  SHA-256 values below identify the local evidence used
for this boundary record.

| Local file | SHA-256 |
|---|---|
| `00_ELVHYSWP4_2 README.txt` | `c6166236c7087ba4a457990470a69c125df436ac20e67b4b35d94e3a3df3fc42` |
| `ELE402HSEMETA.csv` | `84d71b9ee69b7814174a15a98243c7681da8cc1e9a8af7b7ebc1307e28e8d488` |
| `ELVHYS WP4.2 TCS Sensor Details.pdf` | `bedb3115c22e7f080ab10beef6ed0f3269a0997929e56a6455e80e3acd2911b6` |
| `ELE402HSE010CONC20251111.csv` | `534ddacd0b96a5cd04c46c8cce93b3d268012705e45dd28a8904563850c71480` |
| `ELE402HSE010FLMT20251111.csv` | `c09a32947c1bcdee3f5e1cfdecf41a000de4adc85a2fad201f3d0cdc263be24b` |
| `ELE402HSE010MISC20251111.csv` | `d3ed6e270af46af3059771f40509a2e2b4c2d58dc98d87ba798d7bce8301c2f5` |
| `ELE402HSE010PRES20251111.csv` | `845a555edbfd0f7e4e3c2f8e65182556d476ec739bb8633d1497e689e2a21dd0` |
| `ELE402HSE010TEMP20251111.csv` | `ea126c5ce52dc8c80a755ceaa488290040f00ea2ca4a72fed2fafd8eb64eefb2` |
| `ELE402HSE011CONC20251111.csv` | `5eaa47e1ffe1aabea53265043b976821434c5b77db110321d4023d841e86759e` |
| `ELE402HSE011FLMT20251111.csv` | `59a4c49d08ad38cc3b7ae73d62d66c96e51620c154a67dc84f541bc37315f7f5` |
| `ELE402HSE011MISC20251111.csv` | `407b276129d16cd8092e4c4818d76a8a329031fb7cf9530b8e9a1b0110cce4b5` |
| `ELE402HSE011PRES20251111.csv` | `59fd228667d0d8180edaa2f333ec5314a37bbecaf9a6dcaba2a707f5f5f243b` |
| `ELE402HSE011TEMP20251111.csv` | `e06b7621d9cda9967916641c38675e3f2a8d3fd8009fc86ff00d798ba3ff4528` |

The raw CSV digests above are local identity checks; the custodian should
recompute them from the DOI archive before any external promotion.

## IJHE use and promotion rule

The manuscript may cite this lane to demonstrate that the software can retain
time-resolved channels, calibration specifications, common-clock metadata and
measured confined-space geometry without silently pooling incompatible data.
It may not report a new MG/VG/FAC2 score from these tests, claim outdoor
obstacle-wake validation, or infer hydrogen mass flow from pressure and fan
flow.  Promotion requires a custodian-confirmed source boundary, channel-level
calibration certificate, resolved nozzle geometry and a declared enclosure
model.  Until then the lane remains `conditional` with
`promotion_allowed=false`.

Existing implementation and audit references are
[`src/degali/validation/elvhys.py`](../src/degali/validation/elvhys.py),
[`docs/elvhys-tcs-audit.md`](elvhys-tcs-audit.md) and
[`docs/prereg-elvhys-tcs-screen.md`](../docs/prereg-elvhys-tcs-screen.md).
