# Already-atmospheric source schedule input

`degali field-source` is the strict import boundary for a generic field source
schedule that is already expressed as atmospheric H₂ vapour mass flow. It is
useful for a separately qualified post-flash, pool-evaporation, or droplet-
evaporation ledger. It is not a flash solver and does not convert a plant
historian into a leak source.

For the strict direct-flash field case, the `scenario.source` object may either
declare `mass_flow_kg_s` or explicitly declare a `pressure_driven_mass_flow`
object. The latter contains a unit-labelled bounded `ambient_pressure_pa` and
an explicit `source_id` (plus optional `allow_supercritical_gas`); the parser
then applies the existing homogeneous-equilibrium throat closure across all
declared pressure/temperature, opening-area, Cd, and ambient endpoints. The
derived `kg/s` bound and derivation metadata are retained in the request. A
case that omits both forms, or declares both, is rejected; no rate is inferred
from pressure or opening geometry implicitly.

The derived rate is a typed source boundary, not an independent interval added
after the fact. Field uncertainty envelopes propagate the ambient-pressure
corner and recompute the rate from the selected pressure/temperature,
opening-area, Cd, and ambient values, so impossible combinations of a derived
rate and its driving inputs are not manufactured.

The JSON wrapper uses schema `degali.field-atmospheric-schedule-input.v1`:

```json
{
  "schema": "degali.field-atmospheric-schedule-input.v1",
  "csv_path": "vapour_schedule.csv",
  "evidence": {
    "dataset_id": "event-2026-01",
    "path": "vapour_schedule.csv",
    "sha256": "<64 lowercase hex characters>",
    "row_count": 4,
    "source_boundary_id": "source-boundary-01",
    "common_clock_id": "clock-utc-release-01",
    "source_kind": "post_flash_atmospheric_vapour"
  }
}
```

The optional top-level `rate_operator` is `piecewise_constant` by default.
Declare `linear` only when each CSV rate is a value at its time node; the
transport then integrates each adjacent pair with a trapezoidal/linear
operator. The operator is retained in the schedule fingerprint and must match
the lower/nominal/upper corner schedules. It does not alter the default
interval-start interpretation of existing records.

For a Python `ReleaseSource` whose rate boundary is intentionally left at zero,
`pressure_driven_lh2_mass_flow()` can be applied as an explicit opt-in adapter.
It evaluates the existing homogeneous-equilibrium throat closure across the
declared pressure/temperature, opening-area, discharge-coefficient, and optional
ambient-pressure endpoint corners, then returns a provenance-tagged `kg/s`
bound. `release_with_pressure_driven_lh2_mass_flow()` attaches that bound and
records the derivation metadata, but refuses to overwrite any non-zero declared
rate. A thermodynamically unsupported corner raises an error rather than
inventing a source rate. This Python helper is not implicit in the atmospheric
schedule CSV/JSON boundary: imported schedule cases must still declare their
rate history explicitly. The strict direct-flash case can use the separate
`scenario.source.pressure_driven_mass_flow` opt-in described above.

The CSV must contain at least two rows with `time_s` and `rate_kg_s` columns.
Times start at zero and increase strictly. Rates are finite and non-negative;
the schedule must contain positive integrated mass, and the final endpoint rate
must be zero, so no source is silently continued past the declared duration.
The row count and SHA-256 must match the evidence.

The importer reads the header and every data row with exact field widths. A
truncated or extra-value row is rejected with its row number, as are empty-data
files and headers whose stripped names collide (for example `rate_kg_s` and
` rate_kg_s`). This prevents a malformed export from changing the time/rate
mapping before the source-rate physics gate.

The strict JSON wrapper also rejects duplicate or whitespace/case-normalized
object keys before resolving the CSV path and evidence. This keeps a saved
source boundary and SHA-256 pin from being changed by last-value-wins parsing.

Optional `rate_lower_kg_s` and `rate_upper_kg_s` columns must be supplied as a
pair. When present, the adapter preserves lower/nominal/upper schedules on the
same time axis and exposes deterministic `corner_schedules()`; it does not
interpret those bounds as a probability distribution.

`source_kind` is one of `declared_atmospheric_vapour`,
`post_flash_atmospheric_vapour`, `pool_vapour`, or `droplet_evaporation`. The
adapter records this declaration but does not verify the upstream phase ledger.
The in-plane distributed-source contract additionally
accepts
`in_flight_droplet_evaporation`; unknown kinds are rejected rather than routed
as an unclassified scalar source. Combining the resulting schedule with a field case
still requires the normal source, weather, sensor, obstacle, numerical
refinement, and operational applicability gates.

`declared_atmospheric_vapour` and `post_flash_atmospheric_vapour` may be passed
to `atmospheric_source_schedule=`. That argument replaces the primary direct
atmospheric-vapour rate while retaining the explicit scenario source location
and flash contract audit. A pool or droplet schedule has no location or
vertical scalar width in this file, so placing it at the nozzle would be
physically ambiguous; those kinds are rejected there and must use the explicit
distributed-source or phase-routing handoff instead.

In the Python API, pass the parsed object as
`run_field_semi_fv_envelope(..., atmospheric_source_schedule=...)`. The envelope
then combines its lower/nominal/upper schedule corners with the declared
weather/surface/sensor corners, while fixing the upstream source parameters at
their nominal values because the atmospheric schedule replaces that source-rate
boundary. Each case retains the schedule provenance in the envelope report.
For a refinement-backed operational screening, pass the same object to
`run_field_operational_uncertainty_envelope(...)`; every source corner then
receives its own numerical refinement and fail-safe decision gate. This is
still a deterministic sensitivity envelope, not a probability interval or a
design/approval conclusion.

If a nominal `direct_vapour_schedule` is attached through the Python request
API without running its matching envelope, the operational gate only treats
the schedule as replacing the upstream thermodynamic/rate boundary. Declared
source-location or source-direction bounds remain unresolved and keep the
nominal result `withheld`; the schedule does not silently relocate or
re-orient the local source plane.

For a pool or droplet schedule that already has a declared global launch
location and scalar vertical width, use `FieldDistributedVapourSource` in the
Python API instead of placing the schedule at the primary nozzle. Its optional
`position_uncertainty_m` and `vertical_sigma_uncertainty` bounds are expanded
into Cartesian source-shape corners by `run_field_semi_fv_envelope()` and the
operational envelope. The nominal run is withheld if those bounds are present
but have not been propagated. The report retains the original bounds and their
evidence sources; no probability weighting or 3-D wake correction is inferred.

The CLI composes the schedule without weakening the strict field-case parser:

```console
degali field-screen field-case.json --uncertainty-envelope \
  --atmospheric-source-case source-case.json \
  --output field-source-execution.json
```

`source-case.json` uses the schema above. The execution record retains its
resolved path, SHA-256, full schedule evidence, and the source-corner report.
This option is mutually exclusive with `measured_history` and the phase-routing
transport input, and it is refused without `--uncertainty-envelope`.

An imported `FieldAtmosphericSourceSchedule` can be materialized back into a
strict case with `write_field_source_schedule_json()`. The writer verifies that
the evidence CSV has not changed, preserves the source kind/clock and bound
column declarations, emits relative paths plus a fresh SHA-256, and re-reads
the JSON before returning. Existing case files are never overwritten; the
JSON is opened exclusively after schedule validation, so a concurrent writer
cannot replace the pinned source case after its preflight check.

After `degali field-source --output source-execution.json`,
`degali field-verify source-execution.json` rechecks the case SHA and the
fingerprinted schedule CSV. A changed source file fails verification and does
not become a new atmospheric boundary implicitly.

## Distributed source in a strict field case

When a pool or droplet ledger has an explicit global launch location and
vertical scalar width, a base `degali.field-screening-input.v1` case may carry
`distributed_vapour_sources` directly. The source object requires a label,
evidence ID, `position_m`, a unit-labelled `vertical_sigma_m` bounded object,
and a schedule object with `time_s`, `rate_kg_s`, and an explicit `source_id`:

```json
{
  "label": "pool-vapour-A",
  "source_kind": "pool_vapour",
  "position_m": [
    {"nominal": 1.0, "lower": 0.9, "upper": 1.1, "unit": "m", "source": "layout-A"},
    {"nominal": 0.0, "unit": "m", "source": "layout-A"},
    {"nominal": 0.1, "unit": "m", "source": "layout-A"}
  ],
  "vertical_sigma_m": {"nominal": 0.15, "lower": 0.10, "upper": 0.20, "unit": "m", "source": "width-A"},
  "schedule": {
    "time_s": [0.0, 0.5, 1.0],
    "rate_kg_s": [0.02, 0.01, 0.0],
    "source_id": "pool-ledger-A"
  },
  "evidence_id": "pool-evidence-A"
}
```

All three position components must use either numeric values or bounded
objects; mixed arrays are rejected. The final schedule endpoint must be zero
and the integrated mass positive. The source-position and vertical-width
corners are combined with the primary field envelope, while the report keeps
the original evidence-backed bounds. This is an atmospheric handoff, not a
flash solver or a 3-D pool/wake closure.

The distributed-source `schedule` may also carry paired `rate_lower_kg_s` and
`rate_upper_kg_s` arrays. They must use the same source ID/time axis, end at
zero, and contain the nominal rate at every row. The envelope propagates
lower/nominal/upper rate histories as deterministic source corners alongside
position and width corners; the bounds are not treated as a probability
interval. An unresolved nominal distributed schedule remains withheld by the
operational uncertainty gate until these corners are propagated.

```console
degali field-source source-case.json --output source-execution.json
```

The execution record is provenance only. It contains no hazard-distance,
design-basis, or approval conclusion.

At the lower-level Python boundary, a `SourceRateSchedule` may still be
inspected with a retained non-zero endpoint because that endpoint is not
integrated beyond the declared duration. Before finite semi-FV transport,
however, `SemiFVConfig.validate()` requires `has_zero_endpoint` for the
primary and every distributed source. This keeps direct Python use aligned
with the strict CSV and field-request termination contract; a source cannot
silently continue after its declared clock ends.

The measured-history and pool-vapour adapters apply the same rule before
returning their typed schedule wrappers. The measured-flash wrapper also
checks that its declared direct-vapour mass equals the schedule integral; the
pool wrapper checks the corresponding evaporation-ledger mass. A malformed
wrapper therefore fails before it can be attached to a field request.
