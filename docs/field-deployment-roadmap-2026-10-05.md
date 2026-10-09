# DEGALI field-deployment target and first implementation (2026-10-05)

## Target

The target is a broadly deployable DEGALI screening model for LH2 release
scenarios, while preserving the validated free-field DEGADIS-compatible path.
The deployable result must be auditable: source, weather, surface, sensor and
model-form assumptions are explicit, bounded uncertainty is carried through as
an interval/corner envelope, and unsupported geometry fails closed.

This is a reduced-order, semi-finite-volume extension. It is not a claim that
DEGALI has become a full finite-element or CFD solver.

## Implemented in this milestone

- `degali.addons.field_contracts` now provides `ReleaseSource`, `WeatherState`,
  `SurfaceBoundary`, `SensorModel` and `FieldScenario` contracts. Pressure,
  temperature, flow, opening area, discharge coefficient, liquid fraction,
  wind, surface heat transfer and sensor response/bias can each carry an
  explicit lower/nominal/upper interval.
- Source metadata and the declared stability-to-diffusivity map are copied into
  immutable, type-checked mappings at construction, so audit provenance and
  scalar-mixing values cannot be changed through a caller-owned dictionary.
- Direct source geometry uncertainty applies the same physical-unit contract as
  strict JSON: location bounds require `m` and direction-vector bounds require
  dimensionless `1` before deterministic corners are formed.
- Direct Python release, weather, sensor and surface boundaries reject an
  explicitly supplied non-empty unit label that disagrees with the quantity's
  expected SI unit. Empty labels remain a compatibility path for legacy direct
  `BoundedValue` callers; strict JSON still requires exact labels.
- `FieldScenario.corner_cases()` enumerates deterministic interval corners; it
  makes no probabilistic claim and has a configurable case limit.
- Python field boundaries reject boolean numeric coercion and untyped nested
  source/weather/surface/sensor/scenario/request records at construction, so a
  malformed direct API object cannot become a physical `1.0` or fail later
  through attribute access during screening.
- The lower-level atmospheric schedule, obstacle, distributed-source and
  semi-FV transport records apply the same no-boolean-coercion rule before
  solving, keeping the generic leak-source path aligned with strict case-file
  parsing.
- A `SourceRateSchedule` can be inspected as a piecewise record with its
  retained endpoint, but every finite semi-FV primary/distributed injection
  now requires that endpoint to be explicitly zero. The field request, strict
  CSV boundary and low-level solver therefore share one shutoff contract and
  cannot silently continue a source beyond its declared duration.
- `SourceRateSchedule` now keeps the legacy `piecewise_constant` operator as
  its default and accepts an explicitly declared `linear` time-node operator.
  Both modes use the same zero-endpoint and per-source mass ledger gates;
  linear histories are integrated exactly over clipped solver intervals and
  the operator is retained in strict source provenance.
- An opt-in `pressure_driven_lh2_mass_flow()` adapter now derives a bounded
  `kg/s` source rate from the declared upstream pressure/temperature, opening
  area, discharge coefficient, and optional ambient-pressure endpoints using
  the existing homogeneous-equilibrium throat closure. The release wrapper
  records the derivation/source ID and refuses to replace a non-zero measured
  rate; a failed corner remains a typed error. Strict field JSON can invoke
  the same adapter only through an explicit `pressure_driven_mass_flow` object
  carrying the bounded ambient pressure and source ID; omitting both a rate
  and that opt-in object remains invalid.
- The pressure-derived rate boundary is typed rather than left as opaque
  metadata. Field uncertainty envelopes remove the derived `mass_flow_kg_s`
  from the independent cross-product, propagate its ambient-pressure boundary,
  and recompute the throat rate at each selected pressure/orifice/ambient
  corner. This prevents physically impossible combinations of a derived rate
  with the inputs that generated it while preserving deterministic bounds.
- The serialized v1 verifier upgrades an omitted operator only in-memory to
  the unambiguous legacy `piecewise_constant` default; explicit operator
  values still participate in exact provenance matching.
- The local field request can now declare bounded ambient temperature,
  pressure and air-density boundaries in addition to the source/weather
  bounds. Those corners are transported through flash and sensor conversion,
  retained in envelope selections and reports, and an unpropagated nominal
  ambient result is withheld by the operational gate. The phase-routing/pool
  handoff now repeats its phase ledger at each ambient temperature/pressure
  corner and carries air-density corners into the corresponding field sensor
  conversion. The measured-history joint adapter now re-flashes at ambient
  temperature/pressure corners and re-applies detector conversion at ambient
  air-density corners; those selections are retained in the joint envelope
  and cleared only on the corresponding field case. Direct request bounds use
  the same `K`/`Pa`/`kg/m3` unit contract as strict JSON before any corner is
  evaluated.
- The reusable LH2 saturation table now includes the endpoint saturation
  temperatures for the full declared ambient-pressure range, so pressure
  corners are not evaluated against a nominal-only thermodynamic table span.
- A `FieldValidationEvidence` fingerprint is now explicitly treated as
  dataset identity rather than a validation result; applicability remains
  `conditional` until the matched observed-vs-model score is produced.
- `CircularBoundedValue` makes the weather direction interval circular rather
  than linear. A declared 350°–10° record crosses north as a 20° sector and
  runs those two endpoint wind planes; it is not rejected or expanded into a
  spurious 340° sector. The case/report retains this fact as
  `crosses_zero_deg`, while substantial wind movement still requires the
  separate wind-history gate or a transient model.
- `assess_field_applicability()` returns `accepted`, `conditional` or
  `blocked`, with reasons and warnings. Unresolved flash, unsupported fluid or
  unsupported geometry cannot silently produce an accepted result.
- `degali.validation.ffi_source_state` adds a separate public-FFI audit lane.
  `run_ffi_source_state_envelope()` crosses explicitly evidenced source-rate,
  orifice, release-height, storage-pressure, wind and ambient bounds through
  the existing free-field jet, evaluates the result at reported sensor
  coordinates, and preserves lower-bound-aware residual maps. It does not fit
  coefficients, assign corner probabilities or promote FFI data to site
  acceptance; upwind/unreached sensors remain withheld.
- `degali.addons.semi_fv_obstacle` provides a local x-z conservative transport
  operator. A rectangular obstacle is specified by downwind distance, width,
  height and base height. Blocked flux is routed to the nearest available cell
  above/below and the mass residual is reported. Its signed
  `gravitational_settling_m_s` boundary is explicit: positive values settle
  downward toward decreasing z; a regression fixes the centre-of-mass direction.
- `project_cuboid_to_wind_plane()` connects a global `AxisAlignedCuboid` and
  an explicit wind-*to* direction to that local obstacle. A cuboid crossing
  the centreline becomes a conditional wind-plane obstacle; an off-centre or
  source-overlapping cuboid is reported without inventing a lateral-wake or
  inlet correction.
- `run_field_semi_fv_screening()` is an end-to-end direct-vapour screening
  entry point. It maps a field LH2 flash, wind frame, optional obstacle and
  in-plane sensor into the local transport calculation; then applies declared
  sensor t90, gain, bias and averaging. A sensor outside the 2-D wind plane is
  withheld rather than assigned a fabricated off-centre concentration.
- The field flash handoff now carries its phase mass-partition, momentum and
  total-specific-energy residuals into the screening report. Non-finite or
  out-of-tolerance residuals block the field transport before semi-FV
  injection; the direct scalar branch never silently evaporates the retained
  post-flash liquid. The measured-history and cryogenic-blowdown schedule
  adapters use the same gate before producing time-varying atmospheric source
  rates.
- Semi-FV diagnostics retain a source-level injected-mass ledger for the
  primary release and every evidenced distributed source, so combined-source
  conservation can be audited per schedule as well as in aggregate.
- The source ledger now has an independent schedule-mass check: each declared
  schedule's integrated mass is compared with the solver's injected mass and
  the per-source residual plus maximum residual are serialized. A mismatch
  blocks field transport promotion instead of allowing an internally
  self-consistent but schedule-inconsistent ledger to pass.
- The field conservation gate also checks the residual between that ledger and
  total injected mass; inconsistent source bookkeeping is withheld even when
  the mesh inventory residual happens to be small.
- `SourceRateSchedule` and `FieldSemiFVRequest.direct_vapour_schedule` accept
  a declared, already-post-flash atmospheric vapour-rate history. The local
  solver integrates each schedule interval exactly across its own time step;
  it does not infer a flash from an upstream historian. A scheduled source is
  therefore marked as not time-resolved re-flash physics, and scalar corner
  uncertainty is refused until a time-aligned source envelope is supplied.
- `field_source_io.py` and `degali field-source` provide a strict external
  CSV boundary for that generic atmospheric schedule. The file is fingerprinted
  and bound to a source boundary/common clock, requires a zero final endpoint,
  and never infers phase or leak rate from plant pressure/level columns. Paired
  lower/upper rate columns are retained as deterministic same-clock corners,
  not converted into a probability band. The field envelope can combine those
  corners with weather/surface/sensor corners while fixing the replaced upstream
  source boundary at nominal values. The operational uncertainty envelope
  accepts the same adapter so each source corner receives its own refinement and
  decision gate. `field-screen --atmospheric-source-case` composes the separate
  source JSON at that envelope boundary and records its path/SHA without
  allowing it to bypass the measured-history or phase-routing gates. The direct
  composition accepts explicitly declared or post-flash atmospheric vapour;
  pool/droplet schedules remain rejected until their explicit location-bearing
  distributed handoff is used.
- `direct_vapour_schedule_from_cryogenic_blowdown()` now converts every
  positive-duration interval in an existing `CryogenicBlowdownResult` into a
  post-flash direct-vapour schedule. It supports the blowdown model's
  single-phase branch and its explicitly selected HEM withdrawal branch; an
  unsupported two-phase withdrawal is withheld. The exported ledger reports
  total discharge, direct-vapour mass and post-flash liquid excluded from the
  atmospheric scalar source.
- `MeasuredTimeSeries`, `MeasuredReleaseHistory` and
  `direct_vapour_schedule_from_measured_history()` provide the measured-data
  source route. Pressure, temperature and flow channels carry per-sample
  lower/nominal/upper values, source IDs, time offsets and response-time
  metadata. The flow clock defines release intervals; the other channels are
  linearly interpolated only within their declared coverage. Every interval is
  flashed independently into direct vapour and unrouted liquid. A missing
  channel, invalid pressure/phase state, or out-of-range interpolation blocks
  source construction. `request_with_measured_flash_schedule()` carries both
  the resulting schedule and its alignment warnings into the field report.
- `read_measured_history_csv()` is the deliberately narrow historian import
  boundary. It requires a caller-supplied SI column map, an independent
  phase/flash-state evidence ID, per-channel source, calibration,
  response-time, clock-offset and interval-error evidence;
  preserves that mapping and a SHA-256 provenance record through the field
  report; and rejects absent, non-monotonic or nonphysical P/T/mass-flow
  data. It cannot infer a release rate or a normal volume conversion from an
  ordinary pressure-only operational export.
  `direct_vapour_schedule_from_imported_history()` uses a caller-supplied
  table or, for a long selected event, constructs one across all declared
  history-temperature bounds. It retains a direct path for short records and
  records table use/automatic construction in schedule provenance.
- `direct_vapour_schedule_envelope_from_measured_history()` enumerates global
  lower/upper channel selections for deterministic time-history sensitivity;
  it is not a per-sample random draw or a confidence interval.
- `MeasuredHistoryQualityCriteria` and
  `assess_measured_history_quality()` now separate historian diagnosis from
  operational promotion. A P/T/flow history must have an evidenced limit for
  sample spacing, response time, clock offset and interval half-width before
  `request_with_measured_flash_schedule()` will attach it to transient field
  transport. The report retains every channel metric and the pass/fail basis;
  an unapproved history can be inspected or flashed diagnostically but cannot
  silently become an operational source boundary.
- `run_field_measured_history_envelope()` carries every approved global
  P/T/flow/liquid-fraction lower/upper selection through re-flash and crosses
  declared source location/direction corners before wind-plane projection.
  It reports a separate field result for each source-history/placement corner
  and serialises its quality evidence, mass partition, sensor result and
  applicability status. It deliberately does not relabel those deterministic
  source sensitivities as a probability interval or silently combine them with
  a separate weather/sensor/geometry envelope.
  Unless the request provides a table, it builds one bounded LH2 saturation
  table across every declared historian-temperature lower/upper value and
  reuses it for all source intervals and corners; the report records whether
  this accelerated property path was used.
- `PressureDrivenMeasuredHistory` is the explicit pressure-only alternative,
  not a relaxation of the historian boundary.  It requires absolute upstream
  pressure/temperature states, event and phase-evidence IDs, declared opening
  area and discharge coefficient, and a zero placeholder mass-flow boundary.  Each
  interval is recomputed with `release_with_pressure_driven_lh2_mass_flow()`
  before flash, so a pressure trend is never treated as a leak rate.  The
  coupled `direct_vapour_schedule_envelope_from_pressure_driven_history()` and
  `run_field_pressure_driven_history_envelope()` paths propagate coherent
  P/T/liquid-fraction corners into field transport; derived flow is not crossed
  independently.  Opening-area and discharge-coefficient bounds are also
  crossed coherently before throat recomputation. Declared source location and
  direction corners are fixed before wind-plane projection and receptor/
  obstacle applicability in both field-history routes. The joint
  `run_field_joint_pressure_driven_history_envelope()` path also propagates
  ambient temperature/pressure, weather, stability and detector corners before
  transport; the source-only helper rejects unresolved ambient pressure rather
  than dropping it.  A nominal pressure-history parse with unresolved
  opening-area/discharge-coefficient bounds remains diagnostic; the joint
  envelope resolves those source corners before operational gating. This programmatic route still requires independent raw historian
  SHA/event/clock/calibration evidence before any operational promotion.
  `run_field_operational_joint_pressure_driven_history_envelope()` reuses those
  corners for per-case refinement and aggregate fail-safe gating; one withheld
  corner withholds the aggregate result.
- `run_field_joint_measured_history_envelope()` now makes the valid joint
  case explicit: every approved historian P/T/flow/liquid-fraction corner is
  re-flashed at every declared ambient temperature/pressure corner and
  combined with source location/direction, ambient air-density, independent wind
  speed/direction and every in-plane detector response/gain/bias corner. Measured
  source values, opening area
  and discharge coefficient are not varied again, avoiding a false double
  count after measured flow has replaced the orifice prediction.  The case
  count is bounded before transport starts and the report keeps source versus
  field/ambient selections separate. Named supplemental-detector corners reuse the
  exact matching true transport trace instead of re-solving unchanged
  transport; an off-plane detector is refused. It deliberately does not claim a surface
  sensitivity for direct vapour with liquid unrouted: pool/phase-routing
  analysis remains required for that branch.
- `FieldSensorDeployment` adds named supplemental detectors to one shared
  conservative semi-FV transport calculation while retaining the legacy
  `FieldScenario.sensor` as `field_sensor`. Each detector receives its own
  declared response/gain/bias/averaging operator after the shared true scalar
  field is sampled. An off-plane or out-of-domain detector is recorded as
  withheld with its reason, never populated from a fabricated lateral value.
  `run_field_sensor_array_uncertainty_envelope()` evaluates every declared
  detector response/gain/bias lower/upper calibration corner against the
  exact shared true receptor histories, without rerunning transport whose
  governing physics those parameters cannot affect. Withheld detectors are
  refused from that envelope; the audit report distinguishes this deterministic
  calibration sensitivity from sensor-field calibration or a probability band.
- `field_screening_report()` and `write_field_screening_report()` produce a
  versioned JSON audit record for each requested field screen. It retains the
  input contract, applicability status, flash mass partition, wind frame,
  obstacle projection, numerical stability diagnostics and sensor extrema.
  A blocked result contains no invented transport or sensor values. The direct
  screening and sensor-result records also validate typed subrecords, duplicate
  sensor labels, obstacle projections, withheld reasons and trace/transport
  consistency at construction time; the operational gate retains a separate
  missing-sensor defense for incomplete legacy records.
  A direct post-flash `FieldSemiFVRequest` schedule also requires a zero final
  rate endpoint, matching the strict CSV/source boundary and preventing an
  ambiguous source continuation after its declared finite duration.
  Plain uncertainty envelopes, sensor-calibration envelopes and refinement
  results likewise reject empty or malformed deterministic selections before
  their reports are serialized.
- `read_field_screening_request_json()` and `degali field-screen` supply a
  strict, non-inferential case-file route for the nominal direct-flash field
  workflow. Every bounded value carries its unit and source, unknown keys and
  unit substitutions are refused, the case-file SHA-256 is emitted with the
  field/refinement report, and an existing output cannot be overwritten. This
  base format intentionally excludes uncontrolled post-flash schedules and jet
  handoff paths, which must retain their dedicated provenance gates. It now
  accepts only explicit location/width-bearing distributed atmospheric sources
  with finite zero-endpoint schedules; those source-shape corners are included
  in the field envelope. A separate `degali.field-phase-routing-screening-input.v1` contract
  exposes only the conservative rainout/dynamic-pool adapter, declared scalar
  launch boundary, property-table resolution and a strict whitelist of
  unit-labelled gas/droplet numerical-domain controls. Droplet diameter/mass
  classes are explicit, sum to one and carry a separate evidence ID; they are
  not fitted from downstream concentration or pool output. The CLI requires
  `--uncertainty-envelope` for that schema and runs the operational
phase/pool/weather/surface/sensor aggregate. For a transient case, the declared transport duration must exactly
equal release duration plus explicit post-release continuation; no ignored
duration setting is accepted.
The same strict phase object may carry evidence-backed scalar bounds under
`phase_routing.uncertainty`; pool area, evaporation coefficient, pool timestep,
and phase time-window corners are propagated explicitly, while droplet-class
fractions remain exact until separately evidenced. The pool-to-field vertical
scalar width can independently carry a bounded `m` interval, which is retained
as separate source-shape transport corners. When evidence is available,
complete droplet class populations can now be supplied as explicit simplex
corners; the implementation does not renormalize fractions or assign weights.
- `run_field_operational_uncertainty_envelope()` is the operational counterpart
  to the deterministic nominal field envelope. It refines every exact corner,
  applies the same in-plane/review gate to every result, and permits screening
  only when every case is eligible. `degali field-screen
  --uncertainty-envelope` writes that aggregate decision plus every corner;
  `--no-refinement` remains diagnostic and deliberately yields a withheld
  aggregate result.
- `FieldSemiFVEnvelope` rejects duplicate deterministic corner selections at
  construction, preventing a malformed base envelope from using set-based
  completeness checks to hide a duplicated or missing uncertainty corner.
- Measured-history flash, source-history screening and joint history screening
  envelopes now validate typed cases, bound-selection keys, non-empty warnings
  and unique source/field corner selections at construction. Invalid historian
  corners fail before field transport or report serialization.
- Phase-routing result/transport corner records, droplet-handoff refinement
  studies, imported historian artifacts and batch-export manifests now apply
  the same construction-time type, finite-value, duplicate-selection and
  complete-case gates before phase ledgers or JSON reports can consume them.
- Fingerprinted atmospheric source schedules enforce zero final endpoints for
  nominal and lower/upper rate-bound histories even when constructed directly,
  so generic source envelopes cannot bypass the CSV schedule termination gate.
- Field screening results and base envelopes reject list-valued sensor,
  projection and corner collections, preserving immutable deterministic case
  coverage through operational/report layers.
- Operational uncertainty-envelope JSON now adds a deterministic sensor
  extrema summary for peak/final/time-average true and indicated signals.
  It reports available and withheld corner counts per detector and never fills
  missing traces with zero or presents the extrema as a probability interval.
- Joint measured-history, source-plus-sensor and phase/pool operational reports
  use the same summary, including explicit withheld accounting when a physical
  phase corner has no executable screening trace.
- Standalone sensor-calibration operational envelopes use the same summary;
  direct, historian, source-sensor, phase/pool and calibration CLI branches now
  expose one detector-extrema contract.
- A finite transient `ReleaseSource` can also carry an explicit bounded release
  duration (`duration_s` as a bounded object or a separate
  `duration_uncertainty`). The nominal duration remains the source boundary,
  while lower/upper timing corners are fixed into separate envelope cases; the
  transport domain must cover the upper duration. This is a deterministic
  sensitivity bound, not an inferred timing distribution or a fabricated
  schedule.
- Source `location_m` can now be supplied as three explicit bounded metre
  coordinates. Every Cartesian lower/upper combination is inserted into the
  field envelope before wind-plane projection, so obstacle and receptor
  applicability are recomputed per location corner. Mixed numeric/bounded
  coordinates are rejected. The same explicit three-component bounded form is
  available for `direction_m`; zero-vector corners are rejected and valid
  direction corners are retained for jet/phase-routing handoff calculations.
  The direct post-flash scalar closure records those corners but does not
  invent a momentum effect it does not resolve.
- `run_field_operational_joint_measured_history_envelope()` supplies the same
  all-corner operational gate for an approved P/T/flow/(optional phase)
  historian boundary. It takes the existing joint history/weather/detector
  selections, refines every case, carries the imported CSV provenance into
  every refined report, and only then marks source-history uncertainty
  resolved. The field CLI selects this route automatically for a case with a
  controlled `measured_history` object and `--uncertainty-envelope`.
- `field_refinement_report()` and `write_field_refinement_report()` attach the
  same case's explicit grid/time-refinement diagnostics and receptor changes
  to that report. Thus a delivered screening result can show both its physical
  scope limits and whether it met the selected numerical-resolution screen.
- `evaluate_field_operational_screening()` adds the separate operational
  disposition that was previously left to each report reader. By default it
  requires completed transport, a declared in-plane receptor and a completed
  grid/time refinement study. It also withholds a nominal result with
  unpropagated source/weather/detector bounds: surface bounds are excluded
  only because the direct-vapour path does not use them while liquid is
  unrouted. A conditional physical result remains withheld unless the caller
  explicitly sets `allow_conditional=True`; even then the returned record can
  allow only screening, never a design-basis or approval decision. Missing
  refinement, an off-plane detector, unresolved relevant uncertainty or
  blocked transport is fail-safe withheld with required actions.
  Every operational envelope constructor now cross-checks its aggregate
  decision against all retained corners. A hand-built allowed aggregate cannot
  hide a withheld corner or waive per-corner refinement/uncertainty gates;
  each allowed corner must also retain completed screening and refinement
  results.
  Phase sensor/case collections are immutable and non-empty; an empty
  sensor-array aggregate remains withheld without claiming resolved evidence.
  Its decision record keeps reasons and required actions immutable, and a
  directly constructed `conditional_allowed` state must retain an explicit
  explanatory reason plus review action before it can enter a report or batch
  manifest.
  A declared cuboid that is not inserted into the local wind plane (including
  wholly upwind or off-plane geometry) is now also a non-negotiable withheld
  condition for operational screening: its three-dimensional wake cannot be
  treated as harmless by a 2-D result. A cuboid intersecting the plane remains
  a conditional full-plane obstacle representation, not a validated wake model.
- `pool_vapour_schedule_from_phase_routing()` converts a conservative,
  time-resolved dynamic-pool ledger into a separate evaporation-vapour source
  schedule, retaining its rainout impact centroid and maximum wet area. It
  rejects missing pool steps, impact centroid or conservation failure.
  `request_with_pool_vapour_schedule()` requires an explicit
  `PoolVapourLaunchBoundary` with a ground-level scalar closure and evidence
  ID. It never merges pool vapour into the nozzle source; its report records
  that launch temperature, vertical momentum and 3-D footprint are unresolved.
- `export_field_screening_batch()` runs named cases before writing any file,
  then emits one report per case plus a manifest. The output target must be
  new or empty. By default every case requests the numerical-refinement gate;
  if it cannot run, the case remains explicitly blocked in its report rather
  than being silently exported without resolution evidence. Manifest v3 also
  records the separate operational-screening disposition and its required
  actions. Skipping the report refinement does not silently waive the default
  refinement requirement for that disposition. Supplying the optional
  `atmospheric_source_schedules={case_name: parsed_schedule}` mapping routes
  those named cases through the complete operational uncertainty envelope and
  writes the envelope schema in the case report, including source evidence and
  every deterministic source/weather/sensor corner. This batch path retains
  the non-negotiable refinement, in-plane-sensor and resolved-uncertainty gates;
  it rejects a policy that attempts to waive any of them. Conditional-review
  opt-in for these source-schedule cases requires one authorization per named
  case; the record is retained in the manifest while every envelope corner
  still passes the same non-negotiable gates. The envelope also preserves declared source location/direction
  corners; a fixed schedule cannot silently erase source-geometry uncertainty
  or a non-exact release-duration bound.
- `field_batch_io.py` and `degali field-batch` provide a strict named-batch
  JSON boundary around the Python export. Every referenced case/source file is
  fingerprinted in `batch-execution.json`; historian and phase-routing cases
  are refused instead of being silently run through the direct-flash batch
  path. `--require-screening` returns exit 2 while retaining withheld reports.
  Direct `FieldBatchInput` construction also requires an immutable tuple of
  cases, and labels are checked against the safe report-name grammar before
  referenced artifacts are loaded, so malformed batches fail closed at input.
- `run_field_semi_fv_envelope()` runs every declared lower/upper corner without
  relabelling the result as a probability interval. It has an explicit maximum
  case count and reuses one bounded LH2 saturation table across the source
  temperature range when possible. `field_semi_fv_envelope_report()` now
  serialises every corner, table-use status and the explicit statement that
  the envelope is not itself an operational decision; each corner still needs
  resolution and scope review.
- The semi-FV solver now reports the resolved time step plus advection,
  vertical-diffusion and settling Courant numbers and rejects an unstable
  grid/time combination. `run_semi_fv_refinement_study()` compares receptor
  peak and dose over declared grid/time refinements, with a cell-step work
  budget so a field screening request cannot silently become a large
  CFD-like calculation. Numerical convergence remains separate from physical
  obstacle validation. `run_field_semi_fv_refinement_study()` applies that
  same test to the source, wind frame, obstacle projection and in-plane sensor
  of a completed field request, preserving all physical-scope warnings.
  Field screening now also gates the solver's worst-step mass-inventory
  residual with an explicit absolute/relative numerical tolerance. A residual
  beyond that gate keeps the transport diagnostics for audit but marks the
  field result `blocked`; obstacle routing cannot silently turn a non-
  conservative scalar field into an operational screen.
  `degali field-verify` repeats this semantic gate for serialized reports that
  cannot be replayed numerically, including nested batch and other
  option-incomplete envelope reports; the measured-history nominal and
  joint-history paths replay from their fingerprinted CSV, the phase-routing
  path replays from its serialized phase/pool options, and the standalone
  sensor-array plus joint source-sensor paths replay from their serialized
  calibration/source options. The latter also reopens an optional
  fingerprinted atmospheric schedule. The ordinary bounded `nominal_field`
  uncertainty path likewise replays its stored corner options. It
  checks the finite per-source schedule residual map and its recorded maximum,
  and rejects any over-tolerance residual without a matching blocked reason and
  operational `gate_codes` entry. It also checks compact concentration and
  sensor extrema/time/density records, retaining the no-fabricated-result
  invariant for withheld detectors. The serialized operational decision must
  likewise preserve its status/allowance/reason/action invariants and unique
  typed gate codes. Source-level injected-mass labels and totals are also
  cross-checked against the ledger and schedule residual records. The
  verifier additionally binds `primary` and `distributed:<label>` entries to
  the declared direct-vapour/distributed schedule masses retained in
  `transport_input`, so a self-consistent total cannot hide a source-level
  provenance edit. It also recomputes the closed-form final inventory identity
  `mass_injected = mass_domain + mass_outflow` within the numerical tolerance,
  so a forged maximum residual cannot hide inconsistent domain/outflow totals.
  For a strict pressure-derived source, integrity-only verification also binds
  the typed ambient-pressure boundary, source ID, uncertainty record and
  derivation metadata in saved scenario reports, including nested sensor and
  operational-envelope cases, back to the input case. This remains active even
  when the field run is not numerically replayable.
  The current field/semi-FV target is **490 passed**. The last complete
  repository baseline before the opt-in pressure-driven adapter was **1,601
  passed, 145 skipped, 0 failures**; the new adapter has its own focused
  regression coverage.
  The solver now retains the same closed-form final inventory residual in
  typed diagnostics and applies that gate before a direct field result is
  considered completed.
- `StabilityScalarMixingClosure` can map the scenario's declared stability
  class to a caller-declared effective scalar diffusivity, but only with a
  non-empty evidence ID. The selected value and provenance are included in
  the request/report audit record; missing class coverage fails at request
  construction. No undocumented atmospheric-stability correlation or
  neutral-class fallback is selected automatically.
- `FieldWindHistory` captures a provenance-tagged wind-speed/direction record
  and checks the release window against declared speed-range and
  direction-span limits before a fixed-wind-plane calculation is allowed.
  A failed gate blocks the screen; a passed gate is recorded as an explicit
  static-wind applicability check and does not silently substitute a vector
  average for the scenario meteorological boundary. Its vector-mean speed and
  direction must also agree with the nominal `WeatherState` within separately
  declared limits; a stale or mismatched weather input is blocked.
- `field_jet_scalar_handoff_from_near_field()` now converts a
  conservation-screened axisymmetric-jet station into an explicit field
  scalar source boundary. The full direct-H2 flux is checked again against
  the freshly prepared flash rate, while the jet centre height and the
  species-Gaussian vertical marginal set the local source height and scale.
  The release direction must remain within a caller-declared angle of the
  wind plane; every wind-direction corner recomputes that source position.
- `FieldSemiFVRequest.obstacles` accepts multiple labelled global cuboids.
  Each is projected independently into the current wind plane, recorded in
  the JSON audit record, and represented by the union of its local solid
  cells. The existing single-`obstacle` input remains backward compatible;
  mixing the two forms or reusing a label is rejected.
- Facility footprints can be axis-aligned or `OrientedCuboid` rectangles. An
  oriented footprint takes a centre, length, width, height and explicitly
  mathematical global long-axis bearing, so a site drawing is not inflated to
  a global-axis bounding box before each wind-plane projection. Strict JSON
  rejects an unlabelled or meteorological building bearing, and reports both
  the original global geometry and the conservative local mask.
- Strict field case files now require `FieldCoordinateReference`: drawing/grid
  ID, declared origin, vertical datum, evidence ID and the site `+x` axis
  bearing measured mathematically from true east. Meteorological wind is
  rotated into that same site grid before receptors or obstacles are projected,
  so an engineering-grid bearing is not silently treated as east/north.
- `FieldValidationEvidence` makes LH2 validation availability auditable: a
  strict case cannot set `lh2_validation_available=true` without a
  fingerprinted dataset, positive row count, common-clock ID, source/weather/
  obstacle/receptor geometry IDs and temporal-operator ID. Free-field evidence
  remains insufficient for an obstacle transport claim, and the record is
  retained in the field report.
- `degali field-validate` now consumes a separate strict observed-sensor CSV
  contract. It requires one summary row per receptor, repeats and checks the
  common clock/obstacle IDs, verifies SHA and row count, and reports only
  matched-sensor error metrics. It does not create a spatial hazard-distance
  or approval conclusion.
- `FieldSemiFVRequest.post_release_duration_s` explicitly extends a transient
  local calculation after its source ends. Declared schedules remain active
  only over their own duration; a nominal release is converted internally to
  a finite zero-after-shutoff source record. Thus conserved scalar inventory
  can reach a sensor after source shutoff without silently feeding a steady
  plume forever.
- `superpose_field_sensor_branches()` can add completed direct, jet or
  pool-vapour scalar branches at one common fixed sensor. It requires the
  exact same weather, global obstacle, resolved scalar diffusivity, sensor
  calibration and complete, identical sensor time grid. True H2 mass
  concentrations are summed first; t90, gain, bias and averaging are applied
  once to that total. Incompatible branches are rejected instead of silently
  interpolating one branch onto a different temporal operator.
- `run_field_phase_routing_envelope()` now enumerates the declared source,
  wind and surface lower/upper corners into separate rainout/pool ledgers.
  It deliberately fixes sensor calibration at nominal values because those
  settings cannot change the liquid mass ledger, and it reuses one bounded
  saturation table when possible. The returned cases remain sensitivity
  scenarios rather than a probability interval or a merged pool history.
- `run_field_phase_routing_transport_envelope()` carries each completed,
  time-resolved and conservative dynamic-pool ledger into the local semi-FV
  model as a separate distributed pool-vapour source while retaining the
  direct-flash branch. It covers declared surface and stability alternatives
  through to field transport, but it explicitly withholds corners with an
  absent/non-timed pool source and does not resolve in-flight vapour or
  pool-launch temperature/momentum/footprint physics. It is a physical
  envelope. `run_field_operational_phase_routing_transport_envelope()` adds
  every detector calibration corner, numerical refinement and the same
  fail-safe aggregate decision; it permits only conditional operational
  screening after explicit review and never a design/approval decision.
- `DistributedScalarSource` extends the semi-FV transport operator with an
  explicit internal x/z atmospheric-vapour source and its own finite schedule.
  It is conservative, rejects source locations inside a solid cell, and is
  intended as the future handoff target for a trajectory-resolved droplet
  evaporation ledger rather than a license to collapse that ledger into one
  arbitrary point source.
- `FieldDistributedVapourSource` maps an already-atmospheric global vapour
  schedule into that internal source contract. It requires an evidence ID and
  is accepted only when the source lies in the model's wind plane, downwind
  domain and height range. Off-plane, upwind and out-of-domain sources are
  blocked rather than assigned a fabricated two-dimensional concentration.
  Explicit bounded source-position coordinates and vertical scalar width are
  now propagated as separate distributed-source geometry corners through the
  field envelope; a nominal distributed-source screen is withheld by the
  operational gate until those bounds are resolved.
  Multiple distributed sources are crossed independently across their declared
  lower/nominal/upper schedules, and the operational regression retains each
  source's injected-mass ledger rather than checking only the aggregate mass.
  Its position, uncertainty, warning and request collection boundaries are
  immutable tuples, so a constructed field request cannot be mutated after
  corner coverage has been audited.
- Global axis-aligned and oriented obstacle dimensions can likewise carry
  explicit bounded geometry plus an evidence ID. The field envelope expands
  each obstacle-mask geometry corner, reprojects it for every wind corner and
  reruns conservative transport; unresolved nominal obstacle bounds remain
  withheld and never imply a probability or 3-D wake closure.
- `distributed_source_from_pool_vapour_schedule()` connects a conservative
  dynamic-pool evaporation ledger to that internal field source without
  replacing the direct flash-vapour branch. Its vertical scalar width remains
  an explicit caller input, and off-plane pool locations remain blocked.
- `transport_droplet_population(..., trajectory_segment_duration_s=...)`
  can now retain conservative, vapour-mass-weighted evaporation sections for
  each droplet class. `distributed_sources_from_droplet_evaporation()` turns
  only those recorded sections into separate, delayed internal vapour sources;
  the source rates close the in-flight-vapour ledger over the release duration.
  Endpoint-only droplet results are deliberately refused, so no in-flight
  vapour can be silently collapsed to the nozzle or to a single rainout point.
- `run_field_droplet_handoff_refinement_study()` holds a finite field release,
  direct-flash branch, weather, sensor and semi-FV mesh fixed while varying
  only the recorded droplet-trajectory section duration. It reports the
  source count, conserved in-flight vapour mass, true sensor peak and true
  sensor dose against the finest declared trajectory section. Its JSON report
  retains each case's full field-screen audit; a passing handoff-resolution
  study is not a validation of the d-squared coefficient or 3-D droplet
  turbulence closure.
- `compare_field_model_sensor_sets()` formalises DEGALI–SLABx and
  steady–transient sensor comparisons. Source boundary, weather, sensor
  geometry, temporal-operator and (when applicable) obstacle-mask/wake-closure
  identifiers are checked before a numerical difference or LFL classification
  disagreement is reported; an operator or obstacle-representation mismatch
  remains conditional rather than being called a common-input score.
- `field_model_sensor_set_from_csv()` imports a fixed-sensor external result
  only after checking its coordinate columns, concentration units and declared
  averaging window. The imported set carries the source path, SHA-256,
  prediction column and row count. The Test 4 audit now emits the resulting
  machine-readable model-form comparison alongside its paired-sensor CSV;
  native SLABx and DEGALI source mappings plus their different temporal
  operators are deliberately recorded as conditional mismatches.
- `field_model_sensor_set_from_screening()` connects a completed DEGALI field
  trace to that same fixed-sensor comparison contract. The caller must choose
  peak, final or full-trace time-average on either the true or indicated
  channel; withheld sensors and blocked executions fail closed. Execution
  applicability and the selected trace operator are retained as provenance,
  so a conditional field run cannot appear accepted merely because comparison
  basis IDs match.
- `field_model_sensor_set_from_superposition()` applies the same explicit
  trace-operator and provenance boundary to completed direct/pool/droplet
  sensor superpositions. Superposition now also refuses steady/transient
  branch mixing before summing traces.
- Superposition now inspects every branch request for unresolved source,
  sensor, ambient, obstacle, distributed-source, stability, phase-routing or
  measured-history bounds. It retains those labels and marks uncertainty
  incomplete instead of allowing a nominal source sum to masquerade as a full
  deterministic envelope.
- `field_model_decision_impact()` now translates an **accepted** matched-basis
  comparison into threshold-alert agreement/difference at the supplied sensor
  locations, including each path's farthest *monitored* exceedance. A source,
  meteorology, geometry or temporal-operator mismatch produces a withheld
  model-selection conclusion even when numerical rows are available. It never
  interpolates unmeasured points or labels the monitored value a hazard
  distance.
- `run_field_phase_routing()` connects a finite field LH2 release to the
  existing conservative droplet/rainout/dynamic-pool ledger. It requires an
  explicit pool footprint plus surface substrate/evidence ID, heat-transfer
  coefficient and surface temperature. The field surface becomes the declared
  constant flux `h(Tsurface − Tsat,H2)`; it is not mislabeled as a resolved
  soil or concrete conduction model. Pool and in-flight vapour retain their
  unresolved atmospheric launch boundaries.
- Phase-routing gas and liquid model-option mappings are copied into immutable
  snapshots at construction, preventing caller mutation from changing a
  phase-ledger control after its applicability has been assessed.
- `degali.addons.field_lh2.lh2_flash_source_from_release()` connects a nominal
  `ReleaseSource` to the existing equilibrium LH2 flash plane. It requires an
  **absolute** pressure and an **upstream** liquid fraction, applies
  `A_eff = Cd * A_opening` once, and leaves post-flash liquid explicit rather
  than treating it as an immediate atmospheric gas source.
- `prepare_field_lh2_flash()` is the corresponding failure-safe entry point:
  a blocked input returns the reasons plus no usable flash plane, preventing a
  caller from proceeding with an earlier or partially constructed source.
- `LH2SaturationTable` and `build_lh2_saturation_table_for_release()` provide
  a bounded, interpolated table for repeated saturated-H2 flash evaluations.
  The direct CoolProp route remains the reference. A saturation query outside
  the declared table span raises an error rather than extrapolating or silently
  switching property paths.
- `tools/benchmark_lh2_flash_table.py` reproduces direct-vs-tabular timing,
  accuracy and break-even-case estimates for the reference field release.
- Package-level exports and regression tests are included. The current tests
  cover source-contract validation, deterministic uncertainty corners,
  fail-closed gating, obstacle contact, mass conservation and boundary checks.

## Deliberate limits

The obstacle operator does not resolve recirculating wakes, pressure drag,
three-dimensional sheltering, or obstacle-induced turbulence. It therefore
reports `conditional` whenever an obstacle is present and carries a warning
that no site-specific LH2 obstacle validation data were supplied. Existing
`site_geometry` remains the exact trajectory/encounter screen; this new
operator is the conservative local transport closure behind that screen.

The interval corners are sensitivity bounds, not confidence intervals. A
probabilistic uncertainty statement requires a documented calibration or
measurement-error distribution.

The optional stability mixing input is a scalar effective-diffusivity closure,
not a turbulence-resolving RANS/LES treatment. It needs case-relevant evidence
or calibration before it can support an operational clearance decision.

`FieldWindHistory` is an applicability gate, not a time-varying 3-D wind
solver. A record that exceeds its declared steady limits is blocked rather
than averaged into the two-dimensional domain; such a case needs a validated
replay or transient model selection.

The jet handoff is not a 3-D jet/obstacle calculation: lateral dilution, jet
momentum after handoff and obstacle-induced jet interaction remain unresolved.
It is withheld when the conservation screen, flash-rate check, direction
alignment or source-height domain constraints fail.

Multiple projected obstacles are a conservative 2-D solid-mask union, not a
building-array wake model. Their order, recirculation, pressure losses and
lateral bypass are not resolved.

The post-release continuation is an advecting/diffusing 2-D scalar inventory,
not the separate 3-D finite Gaussian-puff formulation. It does not resolve
gravity spreading or evolving crosswind meander after shutoff.

Sensor-level source superposition is conditional even when every branch
passes its own screen. It does not represent source-source interactions,
three-dimensional lateral mixing, nonlinear cold-cloud thermodynamics or a
spatially combined concentration field.

The default screening workflow carries only direct post-flash vapour into the
local atmospheric scalar transport. Liquid droplets, rainout and pool vapour
remain explicit separate branches until they are connected through the
existing phase-routing model and a declared atmospheric launch boundary.

The opt-in pool-vapour handoff now provides a declared **ground-level scalar**
screen after a conservative dynamic-pool ledger. This is still conditional:
the pool source temperature, vertical momentum and lateral footprint are not
resolved by the local 2-D operator. In-flight droplet evaporation can enter
the local screen only through explicitly recorded trajectory segments and a
declared scalar width/evidence boundary. That path conserves its existing
d-squared-trajectory mass ledger, but it remains conditional: turbulent
droplet dispersion, cold-cloud thermodynamics, lateral dilution and any
off-plane trajectory are not resolved.

An operating P/T/flow historian is not itself a transient flash boundary. The
current time-varying path accepts an aligned measured-history flash schedule,
an established **post-flash direct-vapour** schedule, or the project's
explicit well-mixed cryogenic-blowdown model. Alignment uses declared linear
interpolation and preserves response-time metadata; it does **not** infer
unrecorded phase state or deconvolve slow instruments. A measured-history
schedule is now withheld from field transport until evidenced criteria for
its spacing, response time, timing offset and interval uncertainty pass. The
criteria are a quality gate, not a deconvolution, calibration or probability
model.

The CSV historian reader is intentionally not an Excel tag-discovery or unit
inference tool. A plant export must first be selected and normalised under
engineering control, with its calibration and conversion basis supplied in
the import map; this prevents a nominal flow or a pressure transient from
being misrepresented as atmospheric hydrogen mass release.

## Next acceptance gates

1. **외부 검증자료 필요.** The source contract is now connected to the existing
   LH2 flash/table path and the field handoff itself blocks non-finite or
   out-of-tolerance mass-partition, momentum and energy residuals while
   retaining the phase split. The remaining acceptance task is independent
   evidence that the declared source boundary and phase fractions represent
   the site release; the
   nominal and time-aligned upstream-history/blowdown bridges plus the
   trajectory-recorded in-flight-droplet launch gate, its field-sensor
   resolution study, and measured-history quality promotion gate are now
   present; the next extension is site-specific calibration evidence for the
   promoted historian criteria and a relevant LH2 obstacle validation set.
   The connected-project review is recorded in
   [field-evidence-readiness-2026-10-05.md](field-evidence-readiness-2026-10-05.md):
   current plant records can support data-readiness and operating plausibility
   review, but do not yet provide the independently measured flow, phase,
   wind and H2-receptor channels needed for release-dispersion calibration.
2. **구현 완료, 정량 site acceptance는 별도.** `LH2SaturationTable`와
   `tools/benchmark_lh2_flash_table.py`가 반복 CoolProp 조회를 bounded table로
   재사용하고, field LH2 tests가 direct-vs-table 결과를 회귀 비교한다. 개발 PC
   benchmark는 timing 진단이며 물성 정확도 검증을 대체하지 않는다.
3. **외부 matched replay 자료 필요.** Add matched DEGALI/SLABx steady and transient replay cases using the same
   sensors, averaging and wind records; import the external sensor CSV through
   the fingerprinted contract and report model-form differences rather than
   blending them into one fit.
4. **readiness audit 완료, calibration은 미충족.** Use available SLAB-LH2 and operating records for source/weather plausibility
   checks. Do not label the obstacle closure field-validated until a relevant
   LH2 obstacle dataset exists.
5. **구현 완료.** Define the operational report-review procedure and acceptance criteria for
   batch manifests (including who may resolve a conditional or blocked case).
   The strict single-case CLI now requires an auditable `conditional_review`
   record before `--allow-conditional` can be applied and preserves whether it
   was applied in the execution JSON. Batch manifest v3 likewise requires one
   authorization per named case before conditional opt-in and preserves each
   record without allowing it to waive non-negotiable gates. The strict
   `field-batch` input/CLI additionally fingerprints every referenced artifact
   and returns exit 2 when a named case remains withheld.
