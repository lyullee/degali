# DEGALI User Guide

This guide covers DEGALI from first installation through liquid-hydrogen
(LH2) assessments and the advanced DEGADIS-compatible workflows. DEGALI means
**Dense Gas Dispersion for Liquid Hydrogen**. It extends a verified Python
reimplementation of DEGADIS 2.1 with cryogenic-hydrogen physics.

> **Important:** DEGALI is research software for scenario comparison and
> sensitivity analysis. It is not, by itself, a certified basis for regulatory
> distances, facility layout, or safety-critical design. Always inspect the
> returned `warnings`, and independently check safety-critical results against
> experiments or an accepted consequence-analysis workflow.

## 1. Choose a workflow

| Goal | Recommended interface | Input |
|---|---|---|
| Rapid LH2 release assessment | `degali.lh2.assess()` | Physical quantities as Python arguments |
| Steady DEGADIS-compatible run | `run_steady()` or `degali steady` | `.INP` deck |
| Transient DEGADIS-compatible run | `run_transient()` or `degali transient` | `.INP` deck |
| Concentration history and dose at receptors | `TransientOutput.dose()` or `degali dose` | `.INP` deck and distances |
| JETPLU-compatible jet calculation | `run_jet()` or `degali jet` | `.INO` deck |
| Jet touchdown followed by ground dispersion | `run_jet_to_ground()` | `.INO` and `.IN` decks |

New users should begin with `assess()`, which needs no input deck. The legacy
workflows are for users who already possess lawfully obtained DEGADIS input
decks. Original FORTRAN source and third-party experimental data are not
distributed with DEGALI.

## 2. Installation

### 2.1 Recommended environment

- Python 3.10 or newer
- 64-bit Python recommended
- A fresh virtual environment recommended

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install "degali[coolprop]"
```

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "degali[coolprop]"
```

Verify the installation:

```bash
python -c "import degali; print(degali.__version__)"
degali --help
```

CoolProp is required for LH2 calculations. Installing only `degali` is enough
for the traditional DEGADIS compatibility path, but `assess()` may fail to
import CoolProp.

For development and tests:

```bash
git clone https://github.com/lyullee/degali.git
cd degali
python -m pip install -e ".[coolprop,test]"
python -m pytest -m "not slow" -q
```

## 3. First LH2 calculation

This example represents saturated liquid hydrogen stored at 6 bar absolute,
released horizontally at 0.1 kg/s through a 10 mm orifice 0.5 m above the
ground. Wind speed at release height is 2.0 m/s.

```python
from degali.lh2 import assess

result = assess(
    rate=0.10,                  # hydrogen mass release rate, kg/s
    wind=2.0,                   # wind speed at release height, m/s
    height=0.50,                # release elevation, m
    orifice=0.010,              # orifice diameter, m
    storage_pressure=6.0,       # absolute storage pressure, bar(a)
    ambient_temperature=288.15, # ambient temperature, K
    relative_humidity=65.0,     # relative humidity, percent
    ambient_pressure=101325.0,  # absolute ambient pressure, Pa
    max_distance=30.0,          # maximum integration distance, m
    at_distance=10.0,           # evaluate lowest flammable height here, m
)

print(result.report())
```

`storage_pressure` is **absolute pressure in bar(a)**, not gauge pressure. At
an atmospheric pressure of about 1 bar, 5 bar(g) is approximately 6 bar(a).

## 4. `assess()` inputs

| Argument | Unit | Meaning | Default |
|---|---:|---|---:|
| `rate` | kg/s | Hydrogen mass release rate | required |
| `wind` | m/s | Wind speed at release height | required |
| `height` | m | Release elevation; normally zero for a pool | `0.0` |
| `orifice` | m | Diameter of a pressurised-release orifice | exactly one source geometry required |
| `pool_diameter` | m | Pool or evaporation-source diameter | exactly one source geometry required |
| `storage_pressure` | bar(a) | Absolute storage pressure | `1.013` |
| `ambient_temperature` | K | Ambient temperature | `288.15` |
| `relative_humidity` | % | Relative humidity | `65.0` |
| `ambient_pressure` | Pa | Absolute ambient pressure | `101325.0` |
| `max_distance` | m | Maximum numerical integration distance | `100.0` |
| `at_distance` | m | Location for lowest-flammable-height evaluation | optional |

Specify exactly one of `orifice` and `pool_diameter`.

### 4.1 Pressurised jet

```python
from degali.lh2 import assess

jet = assess(
    rate=0.20,
    wind=1.5,
    height=0.50,
    orifice=0.0254,
    storage_pressure=3.0,
    max_distance=40.0,
)
print(jet.report())
```

The pressurised-jet path calculates a flashing equivalent source and a jet
trajectory. `height` must be positive. The present validation primarily covers
horizontal releases aligned with the mean wind.

### 4.2 Pool or low-momentum evaporation source

```python
from degali.lh2 import assess

pool = assess(
    rate=9.5,
    wind=3.0,
    pool_diameter=9.1,
    ambient_temperature=288.15,
    relative_humidity=65.0,
    max_distance=35.0,
)
print(pool.report())
```

The current pool path begins from a specified hydrogen evaporation rate and
source diameter. It does not independently predict time-dependent pool growth
and evaporation from a liquid spill rate. Obtain `rate` and `pool_diameter`
from a separate source model or a documented conservative assumption.

For a declared evaporation history, `assess_pool_history()` preserves the
quasi-steady snapshots. If a receptor and a user-supplied response time are
available, `response_time_s=...` additionally applies a causal first-order
observation kernel and returns `receptor_transient_mole_fraction`. This is a
low-order research surrogate, not transient CFD, and the response time is not
automatically inferred from validation data.

## 5. Interpreting the result

The principal `Assessment` attributes are:

| Attribute | Meaning |
|---|---|
| `distance_to_lfl` | Distance where centreline concentration falls to the hydrogen LFL of 4 mol%, m |
| `distance_to_stoichiometric` | Distance where centreline concentration falls to the stoichiometric mixture, m |
| `lowest_flammable_height` | Lowest flammable-gas height at the requested location or at the end of the flammable envelope, m |
| `regime` | One of `grounded`, `low`, or `aloft` |
| `neutral_buoyancy` | Hydrogen mole fraction where the mixture becomes lighter than air |
| `trajectory` | NumPy array with columns `[x, z, centreline mole fraction]` |
| `warnings` | Applicability or validation-range warnings |
| `notes` | Source description, model path, and post-flash hydrogen mass fraction |

`result.screening_scope` is `qualified`, `conditional`, or `out_of_scope`.
For an application that must never consume an extrapolated number, pass
`strict_scope=True`; this raises `ApplicabilityError` whenever a validation
range warning is generated. The default warning mode is useful for research
and sensitivity work but is not an approval decision.

### 5.1.1 Comparing DEGALI with another model

Use `degali.validation.compare_models` when an independently generated HyRAM,
PHAST, EFFECTS or CFD result is available. `ComparisonCase` records the source
definition, wind vector, release height, receptor operator, averaging window,
phase closure and geometry. A model with different metadata, or a prediction
fitted to the observations, is still scored for diagnosis but is excluded from
the superiority ranking. This prevents a general accuracy claim from being
made from non-equivalent runs.

For a single measured centreline point, `assess_observation_envelope()` enumerates
caller-declared source-rate and wind hypotheses and returns every one within a
specified factor of the observation. Multiple matches are explicitly marked
non-identifiable; the function does not fit or select a correction.
Sensor-height and arc-maximum observations require a separate validated
observation operator.
For declared mast coordinates `(x, y, z)`,
`project_lh2_jet_to_sensors(points_m=[...])` evaluates the jet's Gaussian
vertical/lateral profile directly and returns concentration and temperature;
points outside the integrated trajectory remain `nan`.

A `nan` result usually means that no crossing of the requested concentration
was found within the integration domain; it does not necessarily mean that the
solver failed. Check units, the trajectory, and all warnings before simply
increasing `max_distance`.

```python
import math

if math.isnan(result.distance_to_lfl):
    print("No LFL crossing was found in the calculation domain.")
else:
    print(f"Centreline LFL distance: {result.distance_to_lfl:.2f} m")

print("Cloud regime:", result.regime)
print("Model path:", result.notes["model path"])

for warning in result.warnings:
    print("WARNING:", warning)
```

### 5.1 Export the trajectory to CSV

```python
import numpy as np

np.savetxt(
    "lh2_trajectory.csv",
    result.trajectory,
    delimiter=",",
    header="distance_m,height_m,centreline_mole_fraction",
    comments="",
)
```

### 5.2 Plot trajectory and concentration

Matplotlib is not a required DEGALI dependency, so install it separately:

```bash
python -m pip install matplotlib
```

```python
import matplotlib.pyplot as plt

x = result.trajectory[:, 0]
z = result.trajectory[:, 1]
c = result.trajectory[:, 2]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(x, z)
ax1.set_ylabel("Plume centre height (m)")
ax1.grid(True)

ax2.semilogy(x, c)
ax2.axhline(0.04, color="red", linestyle="--", label="H2 LFL")
ax2.set_xlabel("Downwind distance (m)")
ax2.set_ylabel("Centreline mole fraction")
ax2.grid(True)
ax2.legend()

fig.tight_layout()
plt.show()
```

## 6. Parameter studies

The following example varies wind speed while preserving every warning in the
output. Change one input at a time when building a sensitivity study.

```python
import csv
from degali.lh2 import assess

rows = []
for wind in (0.6, 1.0, 2.0, 3.0, 4.2):
    r = assess(
        rate=0.10,
        wind=wind,
        height=0.50,
        orifice=0.010,
        storage_pressure=6.0,
        max_distance=30.0,
    )
    rows.append({
        "wind_m_s": wind,
        "lfl_distance_m": r.distance_to_lfl,
        "regime": r.regime,
        "warnings": " | ".join(r.warnings),
    })

with open("wind_sensitivity.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
```

Record model version, all inputs, units, and warnings with the results. Cache
completed scenarios rather than recalculating identical cases.

For explicit source/wind uncertainty, use the deterministic envelope helper.
The caller supplies the credible values; DEGALI does not invent a percentage
or fit a multiplier:

```python
from degali.lh2 import assess_envelope

envelope = assess_envelope(
    rates=[0.25, 0.285, 0.32],
    winds=[2.3, 2.5, 2.7],
    height=0.50, orifice=0.0254, storage_pressure=6.0,
    max_distance=6.0,
)
print(envelope.distance_to_lfl_range)
print(envelope.screening_scope, envelope.warnings)
```

This is an input-sensitivity envelope, not a confidence interval and not a
calibration to FFI Test 6.

## 7. Validation range and warnings

The public `assess()` path checks the following directly evaluated ranges:

| Path | Mass rate | Wind | Distance | Geometry |
|---|---:|---:|---:|---|
| Horizontal jet | 0.084–0.285 kg/s | 0.6–4.2 m/s | 0.79–6.0 m | Orifice 0.006–0.0254 m; height 0.5 or 1.5 m |
| Pool | 9.2–10.3 kg/s | 1.55–6.30 m/s | 0–33.8 m | One 9.1 m pool scale |

These are evidence ranges, not numerical clipping limits. The model may run
outside them, but it reports warnings and uncertainty increases. In particular,
a jet exit-velocity/wind-speed ratio below 10 can be strongly wind-steered;
steady mean-wind-aligned concentration predictions are then difficult to
defend.

The present model is not a standalone design basis for:

- buildings, obstacles, wakes, or complex terrain (the optional geometry
  screen below only detects centre-trajectory contact; it is not a wake model);
- indoor or partially confined releases;
- the default `assess()` path for arbitrary directions, strongly crosswind
  jets, or downward impinging jets. A separate `run_lh2_yawed_crosswind_research()`
  path can integrate a horizontal yawed release, but it is research-only and
  remains unvalidated for quantitative design use;
- independent pool spreading and evaporation-rate prediction;
- all possible equipment conditions or certified consequence assessment;
- ignition, flame, thermal radiation, or explosion overpressure.

### 7.1 Research-only site-geometry screen

For a separately calculated three-dimensional centre trajectory, the optional
geometry screen can identify where an **unobstructed** integral-plume result
ceases to apply. Coordinates are global `(x, y, z)` metres. A wind-frame angle
is a **towards** bearing in radians, counter-clockwise from global `+x`.

```python
from degali.addons import (
    AxisAlignedCuboid, TransverseWall, WindFrame, screen_trajectory,
)

path = [(0.0, 0.0, 1.0), (4.0, 0.0, 1.2), (8.0, 0.2, 1.8)]
obstacles = [
    AxisAlignedCuboid(3.0, 4.5, -0.8, 0.8, 0.0, 3.0, label="compressor"),
    TransverseWall(WindFrame(direction_rad=0.0), downwind_m=6.0,
                   base_height_m=0.0, height_m=2.5, half_width_m=None),
]
screen = screen_trajectory(path, obstacles)
if not screen.free_plume_prediction_applicable:
    first = screen.encounters[0]
    print(first.obstacle_label, first.entry_point_m)
    # Do not report the unobstructed concentration field beyond this point.
    assert screen.free_plume_applicable_through(first.entry_arc_length_m - 1e-6)
```

`half_width_m=None` represents an ideal laterally unbounded wall. Its reported
`required_overflight_height_m` is only the minimum geometric clearance for a
continuous centre path. It is not a calculated plume rise. A finite wall may
be passed around or over; the screen deliberately does not choose either
route. `YawedTrajectory.obstacle_screen(obstacles)` is an equivalent shortcut
after an opt-in `YawedCrosswind` calculation.

The screen never applies a wake factor, bypass fraction, wall heat-transfer
term, reflected scalar field, or obstacle-adjusted concentration. A clear
centreline does not establish clearance of the full Gaussian envelope. Use an
obstacle-resolved, separately validated analysis for a downstream building or
wall concentration prediction.

## 8. DEGADIS-compatible input decks

### 8.1 Steady release

Command line:

```bash
degali steady CASE.INP
degali steady CASE.INP --backend coolprop
degali steady CASE.INP --er1 CUSTOM.ER1 --er2 CUSTOM.ER2
```

Python:

```python
from degali import run_steady

profile, source = run_steady("CASE.INP", backend="legacy")
print("Distance to 5 mol%:", profile.distance_to(0.05))
print("Mass above LLC:", profile.mass_above_lfl)

distance = profile.column("dist")
mole_fraction = profile.column("yc")
temperature = profile.column("temp")
```

The steady `Profile` columns are `dist`, `yc`, `cc`, `rho`, `gamma`, `temp`,
`b`, `sz`, and `sy`. Prefer `column("name")` over hard-coded array indices.

### 8.2 Transient release and snapshots

Command line:

```bash
degali transient CASE.INP --snapshot 60 --snapshot 120 --snapshot 180
```

Python:

```python
import numpy as np
from degali import run_transient

run = run_transient(
    "CASE.INP",
    times=np.array([60.0, 120.0, 180.0]),
    backend="legacy",
)

for snapshot in run.snapshots:
    print(snapshot.time, snapshot.distance_to(0.05), snapshot.mass_above_llc)
```

Snapshot columns are `dist`, `dist0`, `yc`, `cc`, `ccstr`, `rho`, `gamma`,
`temp`, `sz`, `sy`, and `b`.

### 8.3 Receptor histories and dose

Command line:

```bash
degali dose CASE.INP --at 50 --at 100 --at 200
```

Python also supports crosswind and vertical receptor offsets:

```python
from degali import Receptor, run_transient

run = run_transient("CASE.INP")
receptors = [
    Receptor(x=100.0),
    Receptor(x=200.0, offsets=[(0.0, 1.5), (5.0, 1.5)]),
]

for history in run.dose(receptors):
    peak_fraction, peak_time = history.peak
    print(history.receptor.x, peak_fraction, peak_time)
    print("First-power concentration dose:", history.dose(exponent=1.0))
```

`DoseHistory.dose()` returns `(mole fraction)^n s`. Do not confuse this with
ppm-minutes or substance-specific toxic-load conventions; convert units and
select an appropriate exponent separately.

### 8.4 Jet and ground-dispersion handoff

```bash
degali jet JET.INO
degali jet JET.INO --bridge GROUND.IN
```

```python
from degali import run_jet, run_jet_to_ground

jet, deck = run_jet("JET.INO")
print(jet.touchdown, jet.distance, jet.halfwidth)

profile, jet, source = run_jet_to_ground("JET.INO", "GROUND.IN")
```

If the plume does not reach the ground at a concentration of interest, the
ground-dispersion handoff cannot run. That can be a physical model result, not
necessarily a solver failure.

### 8.5 Input-deck rules

`.INP`, `.INO`, `.IN`, `.ER1`, and `.ER2` files are free-format in appearance
but **positional**: order determines meaning. Key conventions are:

- wind in m/s, length in m, time in s, temperature in K;
- atmospheric pressure in atm in traditional decks, but Pa in `assess()`;
- concern levels as mole fractions (`0.04` means 4%);
- `CHECK4` true for steady and false for transient calculations;
- `.INO` contains jet conditions and its own property table;
- EPA example numerical coefficients are used when `.ER1` and `.ER2` are
  omitted.

The complete field sequences are documented in
[`src/degali/io/inp.py`](src/degali/io/inp.py) and
[`src/degali/io/jetdeck.py`](src/degali/io/jetdeck.py). Preserve the
original deck and change one field at a time.

## 9. Choosing `legacy` or `coolprop`

- Use `legacy` to reproduce DEGADIS 2.1 and compare historical results.
- Use `coolprop` for supported real-fluid properties and modern thermodynamic
  inversions.
- LH2 `assess()` uses CoolProp internally.

Changing backend changes the physical-property model, not merely numerical
precision. If results differ, inspect temperature, pressure, phase, and the
purpose of the comparison rather than arbitrarily selecting one as correct.

## 10. Troubleshooting

### `ModuleNotFoundError: No module named 'CoolProp'`

```bash
python -m pip install "degali[coolprop]"
```

### `give exactly one of an orifice diameter or a pool diameter`

Supply exactly one of `orifice` and `pool_diameter`.

### `a pressurised jet requires a positive release elevation`

Set `height` to a value greater than zero for a jet calculation.

### A result is `nan`

The requested concentration crossing may be outside the integration domain.
Check units, the last trajectory concentration, `max_distance`, and `warnings`.

### The calculation is slow

Real-fluid LH2 properties, flashing, and jet integration cost more than a
simple Gaussian expression. Begin a parameter study with a few cases and a
shorter `max_distance`, then expand it. Save completed results.

### An input deck cannot be parsed

The deck is positional. One omitted value can shift every subsequent field.
Check the four header lines, declared table-row counts, property-table length,
and logical flags.

## 11. Reproducible reporting

A reproducible report should preserve at least:

1. DEGALI version and DOI;
2. all inputs and units;
3. source path (jet, pool, or legacy deck) and property backend;
4. every returned warning;
5. output CSV files and the execution script;
6. the independent evidence or tool used for comparison;
7. any extrapolation beyond the validation range and applied safety margin.

Check the installed version with:

```python
import degali
print(degali.__version__)
```

The published concept DOI for DEGALI is
[`10.5281/zenodo.22646258`](https://doi.org/10.5281/zenodo.22646258). A
version-specific DOI for the current `0.3.0` source snapshot is
[`10.5281/zenodo.23256538`](https://doi.org/10.5281/zenodo.23256538).

## 12. Further documentation

- [Consolidated physics, validation, and results](docs/technical-reference.md)
- [Current model status and limitations](docs/status.md)
- [Validation process](docs/validation.md)
- [Quantitative field validation](docs/field-validation.md)
- [Evidence and claim grading](docs/claim-grading.md)
- [External-data provenance and reproduction](docs/DATA_AND_REPRODUCTION.md)
- [Public distribution scope](docs/publication-scope.md)
- [Safety and security reporting](SECURITY.md)

When opening an issue, include the operating system, Python version, DEGALI
version, minimal reproducing code, and complete error text. Do not attach
third-party raw experimental data to a public GitHub issue.
