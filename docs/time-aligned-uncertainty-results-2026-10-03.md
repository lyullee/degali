# E3.5 common-window alignment and deterministic uncertainty results

Date: 2026-10-03

## Frozen time operator

Seven transferred PRESLHY E3.5 raw workbooks (Trials 10, 11, 12, 22, 23,
24 and 25) contain a 1 Hz Flowmeter clock and an approximately 3 Hz Xensor
clock. The analysis freezes one operator before scoring:

1. find the centre of the previously frozen sustained-flow window;
2. take an exact common 20 s source window (the shortest retained sustained
   record, Trial 24, is 24 s);
3. average Flowmeter mass flow over that recorded-clock interval;
4. shift the Xensor interval by the documented 18 s sampling-line delay and
   average each sensor over the shifted recorded-clock interval;
5. evaluate the model at the physical sensor coordinates with `averaging=20`
   and take the distance-wise sensor-arc maximum for both observation and
   prediction.

The zero-delay calculation is an observation-operator sensitivity, not an
alternative fitted alignment. Concentrations below 0.05 vol % are excluded
from logarithmic scores. Trials 20 and 21 from the earlier nine-trial score
are not included because their raw workbooks are absent from this transfer.

## Time-aligned result

The 18 s aligned, common-20-s calculation retains 46 distance-wise arcs:

| Trials | Arcs | MG (observed/predicted) | VG | FAC2 | FB | NMSE |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 46 | 0.544 | 4.601 | 0.761 | -0.318 | 0.416 |

This is materially different from the peak-concentration result. MG 0.544
means that, geometrically, the steady prediction is about 1.84 times the
20 s observed mean. The result is not evidence that the previous peak score
was incorrect; it shows that peak and common-window mean are different
observables and must not be pooled.

| Trial | Arcs | MG | VG | FAC2 |
|---:|---:|---:|---:|---:|
| 10 | 6 | 0.486 | 2.833 | 0.500 |
| 11 | 7 | 0.217 | 34.862 | 0.429 |
| 12 | 7 | 0.851 | 1.148 | 1.000 |
| 22 | 6 | 0.222 | 202.571 | 0.667 |
| 23 | 8 | 0.949 | 1.669 | 0.875 |
| 24 | 7 | 0.944 | 1.163 | 0.857 |
| 25 | 5 | 0.669 | 1.249 | 1.000 |

Trials 11 and 22 dominate the remaining scatter. They should be inspected
sensor-by-sensor before any change to the dispersion closure.

## Observation-delay sensitivity

At the primary common-window source and mean wind, removing the documented
18 s sampling-line shift changes the aggregate score to MG 0.510, VG 10.701
and FAC2 0.745 on 47 arcs. The shift therefore reduces scatter substantially
but does not remove the model-high mean-concentration tendency. Because the
report gives approximately 18 s but no probability distribution, these two
operators are reported separately rather than converted into a confidence
interval.

## Deterministic input envelope

The evidence-backed envelope uses:

- source rate: common-20-s mean to the already pre-registered sustained-window
  peak-rate bound;
- wind: campaign reported mean to reported maximum;
- observation operator: documented 18 s delay, with zero delay kept as a
  separate sensitivity.

Across all eight global combinations, MG spans 0.510–0.547, FAC2
0.745–0.766 and VG 4.587–10.876. These input/operator changes do not reverse
the aggregate model-high mean result. This is a deterministic evidence
envelope, not a 95% confidence or Bayesian credible interval.

The common-window flow coefficient of variation is 0.1%–22.3% across the
seven tests. This supports carrying source-rate variation explicitly, but it
does not supply an instrument-error distribution.

## 4 vol % distance implications

The E3.5 distance screen uses a fixed 1.5 m-high centreline receptor and the
four source-rate/wind combinations. It is limited to 40 m:

| Trial | Input-envelope LFL distance |
|---:|---:|
| 10 | 25.0–30.0 m |
| 11 | 39.1 m to >40 m (right-censored) |
| 12 | 17.9–27.0 m |
| 22 | 15.2–19.4 m |
| 23 | 15.3–18.3 m |
| 24 | 11.5–13.2 m |
| 25 | 16.2–19.2 m |

For FFI Test 6, the published low/mid/high mast wind values alone produce
34.3, 36.1 and 37.8 m at the same 1.5 m / 4 vol % definition. This is only a
wind sensitivity. The public files do not provide a defensible uncertainty
distribution for the atmospheric two-phase source state, so this 34.3–37.8 m
range is not a complete hazard-distance uncertainty interval.

## Missing uncertainty components

The current package cannot assign probabilities to sensor calibration,
source pressure/temperature, wind direction and meander, source-state model
form, or their covariance. Adding arbitrary percentages would create a
numerical interval without experimental meaning. A probabilistic propagation
should wait for calibration certificates or stated distributions and a
synchronized FFI source/wind/sensor export.

## Reproduction

```powershell
$env:PYTHONPATH = 'src'
.\.venv\Scripts\python.exe tools\audit_time_aligned_uncertainty.py `
  --output-dir tmp\time-aligned-uncertainty-reproduction
```

The generated `manifest.json` records the complete operator, omissions,
software environment and SHA-256 hashes of all seven raw workbooks and model
files. The derived CSVs and figure are in
`outputs/time-aligned-uncertainty-2026-10-03/`.
