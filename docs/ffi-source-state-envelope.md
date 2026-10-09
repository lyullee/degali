# FFI source-state envelope

`degali.validation.ffi_source_state` adds a reproducible audit boundary around
the public FFI/Spadeadam outdoor horizontal releases. It does not modify the
free-field `hydrogen_jet` equations or promote the FFI data to site acceptance.

```python
from degali.validation.spadeadam import load
from degali.validation.ffi_source_state import (
    FfiSourceState,
    ffi_source_state_envelope_report,
    run_ffi_source_state_envelope,
)

trial = next(item for item in load() if item.test == 6)
state = FfiSourceState.from_trial(trial)
readings = tuple(
    reading for reading in trial.readings
    if abs(reading.radius - 30.0) < 0.5 and not reading.over_range
)
result = run_ffi_source_state_envelope(
    trial, state, readings=readings, corrections=True,
)
report = ffi_source_state_envelope_report(result)
```

The same audit can be reproduced without embedding the extracted tables in an
output artifact:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe tools\audit_ffi_source_state.py `
  --test 6 --radius 30 `
  --rate-bounds 0.82,0.833,0.85 `
  --wind-bounds 2.2,2.3,2.4 `
  --max-cases 4 `
  --output artifacts\ffi-test6-source-state-envelope.json
```

The tool records the reference root path and the versioned report schema; raw
sensor tables remain outside the generated JSON. The CLI and audit tool also
fingerprint `conditions.csv` and `sensors.csv` with SHA-256 so a residual map
can be traced to the exact extracted tables.

For the compact Test 6 observation-operator check, use:

```powershell
.venv\Scripts\python.exe tools\audit_ffi_sensor_operator.py `
  --radius 30
```

The JSON contains two aggregate tables and no raw sensor rows. `arc_max_table`
maximises the reported and projected values over the retained heights and
bearings at each radius. `sensor_height_table` keeps the exact-coordinate
projection but groups it by reported sensor height, so the height dependence is
visible without confusing a centreline value with the finite sensor operator.
Both aggregate observed values are explicitly marked as lower bounds because a
plume can pass between the discrete sensors. The Test 6 local extract currently gives
21.0 vol% observed versus 7.4977 vol% projected at 30 m, with height-group
ratios of 2.9191, 2.8170 and 2.4808 for 0.1, 1.0 and 1.8 m respectively.

The equivalent package CLI is:

```powershell
.venv\Scripts\python.exe -m degali.cli ffi-source-state `
  --test 6 --radius 30 `
  --rate-bounds 0.82,0.833,0.85 `
  --wind-bounds 2.2,2.3,2.4 `
  --max-cases 4 --require-complete `
  --output artifacts\ffi-test6-source-state-execution.json
```

`--require-complete` returns exit code 2 if any selected sensor is withheld in
any corner. The command does not turn a partial or conditional source-state
comparison into a successful validation. Even a calculation with
`status="complete"` carries `operational_screening_allowed=false` and
`validation_qualified=false`; that status only means the requested numerical
corner calculations completed.

Every `FfiSourceState` field is a `BoundedValue` or circular direction bound
with an explicit evidence source. The runner enumerates the Cartesian lower and
upper endpoints (exact values produce one corner), evaluates the plume at each
reported sensor coordinate, and retains a sensor-level residual map. A sensor
that is upwind or outside the computed trajectory is recorded as withheld; it
is never replaced by zero. Over-range readings are excluded and counted.

The envelope report also contains `observation_operator.arc_max_table` and
`observation_operator.sensor_height_table`. These are corner-wise aggregate
views of the same residual rows: the first maximises over retained heights and
bearings at each radius, while the second maximises over bearings at each
reported height. Every row carries its exact source-state `selection`, case
status, available/withheld sensor counts, lower-bound satisfaction fraction and
an `operator_status` of `satisfied`, `violated`, `partial` or `withheld`.
The separate `lower_bound_constraint_status` retains whether the available
rows violate the one-sided bound even when a mixed available/withheld group is
reported as `partial`.
This makes a source/wind corner that underpredicts a censored arc visible
without turning that diagnostic into a symmetric validation score. A corner
with no available sensor prediction never receives a fabricated zero or ratio.

The report schema is
`degali.ffi-source-state-envelope.v1`. Its predicted lower/upper values are
deterministic sensitivity limits, not confidence bounds. FFI arc and sensor
peaks are lower bounds because the plume can pass between sensors, so the
descriptive MAE/RMSE and bias fields are not an acceptance score. No source or
model coefficient is fitted. Each sensor entry also retains the residual lower
and upper bounds, worst absolute residual and the exact corner selection that
produced it, so a large mismatch cannot be detached from its source/weather
assumptions. Because these FFI values are lower bounds, the report additionally
calculates one-sided constraint satisfaction fraction and worst lower-bound
deficit. These are censor-aware diagnostics, not symmetric accuracy scores.

The envelope currently covers source rate, orifice, release height, storage
pressure, wind speed/direction, ambient temperature/pressure, humidity and
wind reference height. It remains a steady free-field approximation applied to
a campaign whose pad includes documented structures; plume storage,
time-varying source histories, gust-induced deformation, obstacle wakes and
sensor response are outside this path. See
`docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md` and use the field
observation and measured-history contracts for those separate operators.
