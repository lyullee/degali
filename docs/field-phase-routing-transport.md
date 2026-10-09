# Phase-routed pool vapour in the local field transport model

`run_field_phase_routing_transport_envelope()` connects a declared LH2
phase/rainout/dynamic-pool calculation to the conservative local semi-FV
transport operator. It runs every source, wind, surface, and—when declared—
stability corner of `run_field_phase_routing_envelope()`. For each completed
corner, it retains the direct-flash vapour source and adds only the timed,
mass-conservative pool evaporation ledger as a separate internal vapour
source. If the request declares bounded ambient temperature, pressure or
air-density boundaries, the phase ledger is repeated at each temperature/
pressure corner and the matching air-density value is carried into the local
field/sensor request; the selected ambient values remain in the case record.

```python
from degali.addons import (
    FieldPhaseRoutingConfig,
    PoolVapourLaunchBoundary,
    field_phase_routing_transport_envelope_report,
    run_field_phase_routing_transport_envelope,
)

result = run_field_phase_routing_transport_envelope(
    field_request,                 # no existing atmospheric schedule/source
    FieldPhaseRoutingConfig(
        post_release_duration_s=30.0,
        puff_duration_s=5.0,
        pool_area_m2=1.0,
        pool_time_step_s=0.1,
        evaporation_coefficient_m2_s=1e-7,
    ),
    PoolVapourLaunchBoundary(
        source_height_m=0.0,
        evidence_id="pool-launch-review-01",
    ),
    pool_vertical_sigma_m=0.2,
)
report = field_phase_routing_transport_envelope_report(result)
```

The same path is available through the dedicated strict JSON schema
`degali.field-phase-routing-screening-input.v1`. Start from
[`field-phase-routing-case.example.json`](field-phase-routing-case.example.json)
and run the complete operational envelope:

```bash
degali field-screen docs/field-phase-routing-case.example.json \
  --uncertainty-envelope --allow-conditional \
  --output new-phase-routing-field-report.json
```

`--uncertainty-envelope` is mandatory for this schema. The CLI therefore
cannot present one nominal pool history as if phase, surface, weather and
detector uncertainty had been resolved. The JSON contract accepts only the
published `FieldPhaseRoutingConfig` scalar boundaries, an explicit
`ground_level_scalar` launch with evidence ID, its vertical scalar width, and
the LH2 property-table node count. The phase-routing scalar assumptions,
including pool footprint and d-squared evaporation coefficient, require a
separate evidence ID. Gas-domain/grid/time controls and the
maximum droplet-flight integration time are required in dedicated, unit-labelled
objects because the coupled solver's generic defaults are not silently assumed
to cover every release. These are whitelisted numerical-domain controls, not
wake or evaporation calibration factors. The droplet diameter classes and
liquid-mass fractions are also explicit and require their own evidence ID;
they are not inferred from the desired pool or sensor result. Fractions must
sum to one. An insufficient development, crosswind or droplet-flight domain is
reported as a withheld corner. The schema rejects a simultaneous historian
source, arbitrary imported post-flash/pool schedules, all other
subordinate-model options and undeclared wake factors.

When phase/pool scalar uncertainty is available, add a
`phase_routing.uncertainty` object. Each declared field uses the same explicit
`nominal`, `lower`, `upper`, `unit`, and `source` form as the source contract.
Supported fields are `post_release_duration_s`, `puff_duration_s`,
`pool_area_m2`, `pool_time_step_s`, and `evaporation_coefficient_m2_s`. The
pool launch `pool_vertical_sigma_m` may likewise be a bounded object (or use
the `pool_vertical_sigma_uncertainty` sibling field) with unit `m`; its
lower/upper values create separate transport source-shape corners. The
envelope evaluates the Cartesian lower/upper corners and retains the selected
values in every phase case; it does not infer bounds or assign probabilities.
A declared post-release interval makes the transport-duration requirement use
its upper bound so a late pool source cannot be clipped. Droplet-class and
mass-fraction uncertainty may be declared under
`phase_routing.uncertainty.droplet_population` as complete class-population
corners plus an evidence ID. Each corner must preserve the mass-fraction
simplex and include the nominal population; no fractions are renormalized.

`--allow-conditional` additionally requires the case's auditable
[`conditional_review` record](field-conditional-review.md). That record cannot
waive a missing pool ledger, refinement, in-plane receptor, represented
obstacle or uncertainty corner, and it never authorizes design-basis or
approval use.

For a transient case, `transport.duration_s` must equal release duration plus
`phase_routing.post_release_duration_s` (or its declared upper bound when a
phase uncertainty interval is present). A duplicate top-level
`post_release_duration_s`, when supplied, must match exactly; it cannot silently
change the pool continuation window.

The field transport duration must cover `source.duration_s +
post_release_duration_s`; a pool schedule extending beyond that duration is
withheld rather than clipped. An existing direct-vapour schedule or a
pre-existing distributed source is rejected, since this adapter must own the
unambiguous direct-flash-plus-pool source combination.

The `FieldPhaseRoutingConfig.post_release_duration_s` is authoritative for
the derived field request. It is copied into that request, so the direct-flash
source is transparently set to zero after release and the dynamic-pool schedule
uses the same continuation window. A contradictory post-release value on the
base request is not silently reused.

This is not a pool-vapour CFD closure. It leaves pool-vapour temperature,
vertical momentum, lateral footprint dilution, three-dimensional obstacle
wake, and in-flight droplet evaporation unresolved. A missing/nonconservative
or non-time-resolved pool ledger produces an explicit withheld corner instead
of a source history. The returned envelope is a deterministic physical
sensitivity report, not an operational approval, probability interval, or
design basis. Sensor calibration and numerical-refinement gates remain
separate requirements before any conditional screening use.

For the complete deterministic phase/pool/sensor disposition, use
`run_field_operational_phase_routing_transport_envelope()` after constructing
the same request, phase configuration and launch boundary. It runs every
physical phase-routing corner and every detector-calibration corner, refines
each completed field result, and returns `withheld` if any one cannot pass the
same fail-safe gate. `allow_conditional=True` remains an accountable-review
opt-in for screening only; neither API can allow design-basis or approval use.

The operational report includes the same `deterministic_sensor_envelope`
summary as the direct field path. It aggregates only completed phase/sensor
corners and keeps a phase corner with no executable trace explicitly
`withheld`. Its typed decision carries the general `transport_incomplete`
code plus `phase_routing_incomplete` when the phase ledger did not complete, or
`pool_transport_withheld` when the conservative pool-to-field handoff failed;
the aggregate therefore preserves the failure boundary without parsing the
human-readable warning.
