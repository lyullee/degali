# FFI/DNV Test 6 transient receptor replay

## What was added

`degali.addons.transient_receptor` is an opt-in research observation operator.
It keeps DEGALI's spatial plume equations unchanged and performs three explicit
operations after the steady calculations:

1. convert each measured meteorological wind direction from “wind from” to
   the downwind east/north vector;
2. rotate every fixed sensor into instantaneous alongwind/crosswind
   coordinates and interpolate between precomputed steady wind-speed runs;
3. optionally apply a causal first-order detector response using a declared
   90 % response time, `tau = t90 / ln(10)`.

The implementation rejects wind-speed extrapolation. Thus a 10 Hz, 120 s wind
record does not require 1,200 full plume integrations: only the wind speeds in
the small steady table are integrated.

## Evidence boundary

The public DNV report does contain time-history plots for Test 6. It states
that non-pressure channels were recorded at 10 Hz and supplies the instrument
response specification. Appendix C shows wind direction, wind speed and the
individual 30/50/100 m hydrogen traces. Therefore the former wording “time
series are unavailable” was too broad.

What is not publicly supplied is the machine-readable raw channel export. The
repository consequently does not embed a reconstructed trace. A digitized or
owner-exported CSV must identify its provenance and is hashed by the audit
tool. The expected columns are:

```text
time_s,wind_speed_ms,wind_direction_from_deg
```

Example:

```powershell
python tools/audit_ffi_test6_transient_receptors.py wind.csv `
  --table-winds 0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6 `
  --response-t90 6 --start 20 --end 140 --output artifacts/test6-replay.json
```

`6 s` is an upper-bound response case based on the reported “less than 6 s”
specification, not a fitted value. Run `--response-t90 0` as the no-lag bound.

The replay writes `degali.ffi-test6-transient-receptor-execution.v1`. It keeps
the supplied history out of the JSON, records its SHA-256 and row count, and
stores one compact row per reference receptor with true/indicated window
mean and maximum. Reported peaks remain explicit lower bounds; over-range
observations are retained as `null` rather than treated as exact values.
The input record also preserves sample count, time range, interval range,
uniform-sampling diagnostic and wind-speed range. These are diagnostics only;
the replay does not promote an irregular or sparse history automatically.
The replay is an observation operator only: `promotion_allowed` and
`validation_qualified` are always `false`.

The saved artifact can be checked without rerunning the plume calculation:

```powershell
python tools/audit_ffi_test6_transient_receptors.py --verify artifacts/test6-replay.json
```

Verification rechecks the wind-history SHA, exact CSV contract, reference
`conditions.csv`/`sensors.csv` digests, receptor identity and compact-window
invariants. It is integrity-only and reports `report_recomputed=false`.

## What this does not solve

This is not a transient dispersion solver. It does not add plume storage,
travel-time memory, gust-induced deformation, turbulence intermittency, or a
time-varying flashing source. It isolates whether wind meander plus detector
lag can explain part of the fixed-array mean/maximum difference. A model
change is justified only if the replay is reproducible from a declared input
trace and remains consistent with the other FFI arcs.

## Ground thermal bound

The independent-energy crosswind research model now also exposes a separate
fixed ground-temperature bound. Set
`ground_heat_transfer_coefficient_w_m2_k` and
`ground_surface_temperature_k` when constructing
`IndependentEnergyCrosswind`, then opt in at integration with
`include_ground_heat_transfer=True`. The added energy source is

```text
q' = h * max(T_surface - T_cloud, 0) * contact_width
```

and is exactly zero without geometric ground contact. Coefficient zero is the
adiabatic bound. This supports the EFFECTS paper's reported adiabatic/fixed
surface sensitivity without changing the default solution or inferring a
coefficient from Test 6. It is still a fixed-temperature boundary, not a
finite-conductivity slab with progressive cooling.
