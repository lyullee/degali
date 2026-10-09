# Source/weather/receptor common time history

`src/degali/addons/transient_receptor.py` now provides a causal packet-history
operator for the first comparison work package. It is intentionally narrower
than a transient CFD solver:

1. `SourceHistory` accepts an already-qualified atmospheric H2 mass-rate
   history. The final zero rate is mandatory. Optional pressure and temperature
   channels are carried as packet state, but are never converted into a new
   discharge or flash closure here.
2. `SourceHistory.to_packets()` uses deterministic midpoint quadrature. The
   packet subdivision count is explicit and the released mass is checked against
   the integrated source history.
3. `WindHistory` is sampled on the same event clock. Packet centres are advected
   by the integrated horizontal wind vector from packet release to observation;
   this retains travel-time memory instead of rotating a steady field at each
   receptor time.
4. `source_direction_to_deg` remains packet metadata and is not replaced by the
   ambient wind direction. A near-field kernel may use the fixed nozzle axis
   explicitly.
5. `replay_fixed_receptors_with_source_history()` sums a caller-supplied packet
   mole-fraction kernel, rejects negative or over-unity contributions, and then
   applies the declared sensor `t90` response. No concentration multiplier or
   source parameter is inferred from observations.

## CSV intake and provenance

`src/degali/addons/time_history_io.py` is the strict file boundary for the
source and weather channels:

- `read_source_history_csv()` requires an explicit event-clock column,
  atmospheric H2 mass-rate column, at least two rows, a zero-rate terminal row,
  and a `SourceHistoryCsvMap` containing the event/evidence IDs. Pressure,
  temperature, and fixed source direction are optional state channels and are
  retained without being converted into a new flash or discharge closure.
- `read_wind_history_csv()` requires event-clock, speed, and meteorological
  `from`-direction columns. Directions are normalized only by the existing
  `WindHistory` operator; no wind channel is selected by concentration error.
- Both readers return a SHA-256-pinned provenance record and reject duplicate
  headers, missing cells, non-finite values, non-monotone clocks, and ambiguous
  extra CSV fields. The input schema is
  `degali.source-weather-time-history-csv.v1`.

An interval table such as `start_s,end_s,speed_m_s` is not silently treated as
a node-sampled history. It must first be converted to event-clock nodes by an
explicit, documented preprocessing step that preserves the original file hash
and interval semantics. This prevents the FFI anemometer bins from being
mistaken for a synchronized receptor history.

Example:

```python
from degali.addons.transient_receptor import (
    FixedReceptor, SourceHistory, WindHistory,
    replay_fixed_receptors_with_source_history,
)

source = SourceHistory(
    time_s=(0.0, 1.0, 2.0),
    mass_rate_kg_s=(0.8, 0.4, 0.0),
    pressure_pa=(210_000.0, 205_000.0, 200_000.0),
    temperature_k=(25.0, 25.5, 26.0),
    source_direction_to_deg=(0.0, 0.0, 0.0),
)
wind = WindHistory(
    time_s=(0.0, 1.0, 2.0),
    speed_m_s=(6.7, 6.5, 6.8),
    direction_from_deg=(270.0, 275.0, 270.0),
)

def packet_kernel(packet, age_s, relative_east_m, relative_north_m, height_m):
    # Replace with an independently qualified transport kernel. This example
    # deliberately does not infer one from sensor concentrations.
    return 0.0

traces = replay_fixed_receptors_with_source_history(
    source, wind, [FixedReceptor("S01", 30.0, 0.0, 0.5, response_t90_s=2.0)],
    packet_kernel,
    packet_subdivisions_per_interval=4,
)
```

The resulting `ReceptorTrace` has a common observation clock and the measured
wind at each observation time. Its `alongwind_m` and `crosswind_m` arrays are
`NaN` by design because a packet train can contain multiple local transport
directions; a caller must not mistake one guessed instantaneous frame for the
packet cloud geometry.

This path is ready to consume a hash-pinned `FieldEvidenceManifest` package,
but it does not promote a comparison by itself. Source boundary, weather,
receptor geometry, sensor response, common-clock, and temporal-operator IDs
must still match in the existing comparison gate before a model-selection
impact is issued.

For a complete join, use `run_manifest_bound_time_history()`. It verifies the
manifest digests first, compares event/clock/receptor/operator IDs, and then
returns a schema-tagged execution record. Missing sensor/operator metadata
produces an executable but `conditional` result; digest drift, ID mismatch, or
an uncovered packet/observation clock produces `withheld` with no traces.

The current FFI Test4 inventory is recorded in
`outputs/ffi-test4-time-history-readiness-2026-10-08-v2.json`. It is deliberately
`withheld`: weather bins exist, and a one-row P04 source-state handoff provides
pressure/temperature/quality, but the source still has no `q(t)` history and the
receptor file has only 275 s aggregate means. Those files are useful
bulk/source-state diagnostics, not a synchronized source-to-sensor time
history.
