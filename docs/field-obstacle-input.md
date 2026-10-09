# Field obstacle geometry input

`degali.field-screening-input.v1` accepts a top-level `obstacles` array. All
source, detector and obstacle coordinates use the same global Cartesian metres
as `scenario.source.location_m`: `x` and `y` are horizontal, and `z` is the
declared vertical datum. A strict case must also include `coordinate_reference`
with its drawing/grid ID, origin ID, vertical datum and the site `+x` axis
bearing relative to true east. The solver rotates the meteorological wind into
that site grid, then projects each declared footprint into the wind-*to* plane
for every wind-direction sensitivity corner.

An axis-aligned cuboid is appropriate when a drawing already uses the global
`x/y` axes:

```json
{
  "x_min_m": 12.0,
  "x_max_m": 15.0,
  "y_min_m": -1.0,
  "y_max_m": 1.0,
  "z_min_m": 0.0,
  "z_max_m": 2.5,
  "label": "pump-skid-A"
}
```

For a rotated rectangular footprint, use its centre, long/short dimensions and
an explicit global-axis convention:

```json
{
  "center_m": [13.5, 0.0],
  "length_m": 4.0,
  "width_m": 1.2,
  "z_min_m": 0.0,
  "z_max_m": 2.5,
  "long_axis_bearing_deg": 35.0,
  "long_axis_bearing_convention": "math_to",
  "label": "angled-transfer-skid"
}
```

`math_to` means counter-clockwise from global `+x`; it is **not** a
meteorological wind-*from* direction. The parser rejects another convention
instead of guessing whether a drawing bearing was reversed. It also records the
pre-projection global geometry and every local wind-plane projection in the
JSON report.

### Bounded obstacle geometry

If a drawing review leaves a dimensional interval, every geometric dimension
must be an object with `nominal`, `lower`, `upper`, `unit` and an explicit
`source`, and the obstacle must add `uncertainty_evidence_id`. Axis-aligned
geometry uses `m` bounds for all six faces. Oriented geometry uses `m` bounds
for centre, length, width and vertical limits plus a `deg` circular bound for
`long_axis_bearing_deg`. Mixed numeric/bounded geometry is rejected. For each
declared bound the field envelope expands deterministic Cartesian corners,
reprojects the footprint for the wind direction and reruns the conservative
mask/transport; it never interprets the interval as a probability, lateral
bypass or wake correction.

After projection, each local rectangle must lie strictly inside the declared
semi-FV domain. An obstacle that touches or extends past the inlet or outlet is
rejected before the mesh is built; the solver never clips a declared footprint
to the finite-volume grid.

An obstacle intersecting the wind-plane centreline becomes a full-plane solid
mask. This routes blocked scalar flux above/below the mask conservatively and
reports mass residual, but remains conditional: it does not resolve lateral
bypass, building-array recirculation, drag, wall heat exchange, or wake
turbulence. A declared obstacle that is off-plane, wholly upwind, or overlaps
the source is not silently discarded; the operational decision remains
withheld until a validated 3-D analysis or an accountable geometry review
shows it lies outside the decision domain.

The operational decision exposes the same boundary as typed codes:
`obstacle_transport_conditional` identifies a represented reduced-order
obstacle, `transport_mass_residual_exceeded` identifies a failed inventory
residual gate, and `source_ledger_mismatch` identifies a per-source ledger
failure. `source_schedule_mass_mismatch` identifies a mismatch between a
declared source schedule's integrated mass and the solver's injected mass.
The numerical residuals and scope warnings remain in the report.
Typed semi-FV diagnostics also reject non-finite or negative mass/Courant
values, invalid diverted fractions, duplicate source-ledger labels and unknown
applicability states before a report can be emitted.

The local transport control `gravitational_settling_m_s` uses a signed vertical
velocity convention: positive values move scalar mass downward toward decreasing
`z`, while negative values represent an explicitly supplied upward drift. This
is a numerical transport sensitivity, not an inferred droplet or cold-cloud
closure; phase-routed liquid and droplet paths retain their own evidence and
applicability gates.
