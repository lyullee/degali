# DEGALI model hardening — 2026-09-21

## What was fixed

The LH₂ entry point now makes the validation boundary executable instead of
leaving it only in prose.

- `Assessment.screening_scope` classifies a result as `qualified`,
  `conditional`, or `out_of_scope`.
- `assess(..., strict_scope=True)` raises `ApplicabilityError` whenever the
  requested rate, wind, geometry, distance, or momentum regime is outside the
  evidence-backed range.
- `assess_envelope(rates=[...], winds=[...])` evaluates caller-supplied
  source/wind bounds as a Cartesian sensitivity envelope. No uncertainty
  percentage, source multiplier, or validation fit is invented.
- `assess_pool_history(..., strict_scope=True)` applies the same guard to each
  quasi-steady pool snapshot. The time-resolved interface remains explicitly
  a source-snapshot/advection label, not a transient puff solver.
- `run_lh2_yawed_crosswind_research()` now exposes the existing conserved
  six-flux yaw kernel for horizontal non-aligned releases. It is a separate
  research branch with global wind/release bearings, reverse-axial rejection,
  and `validated = False`; it is not mixed into the primary score.
- `assess_pool_history(..., response_time_s=...)` now offers an explicit
  causal first-order receptor-response surrogate for declared source history.
  The time constant is user supplied and the result remains research-only; no
  transient storage or turbulence time scale is fitted.
- `degali.validation.compare_models` now provides an apples-to-apples model
  comparison gate. It requires the same source definition, wind vector,
  release height, receptor operator, averaging window, phase closure and
  geometry. Mismatched or observation-fitted predictions are still scored for
  diagnosis, but are excluded from any superiority ranking. This makes a
  future HyRAM/PHAST comparison reproducible without bundling those programs
  or third-party raw data.
- `degali.validation.evaluate_screening` turns the scope and geometry rules
  into an operational gate. A qualified clear-path result can be used for
  conditional engineering screening; a conditional result requires explicit
  opt-in; out-of-scope results or obstacle contact are rejected. The gate
  always reports `design_basis_allowed=False` and `approval_allowed=False`.
- `assess_observation_envelope(rates=[...], winds=[...],
  observed_mole_fraction=..., distance_m=...)` now evaluates every declared
  source/wind hypothesis with the centreline observation operator. It reports the factor-bounded
  admissible set and marks multiple admissible rows as **non-identifiable**.
  No closest row is selected and no source multiplier is fitted, so the FFI
  Test 6 residual becomes a quantitative identification result rather than an
  unqualified “model failure” or an invented correction.
- `project_lh2_jet_to_sensors(points_m=[(x, y, z), ...])` exposes the same
  Gaussian vertical/lateral observation operator at exact mast coordinates.
  The centreline diagnostic is therefore no longer forced onto sensor-height
  data; points outside the integrated trajectory return `nan` instead of being
  extrapolated. Time synchronisation and meander remain separate inputs.
- `tools/audit_ffi_sensor_operator.py` replays the local public Test 6 extract
  without packaging it and records only aggregate projected values and the
  remaining applicability warnings.

The default warning mode is retained for exploratory research. Production
screening code should use strict mode and record the returned scope, warnings,
inputs, and model version.

## Evidence state after hardening

The frozen public comparison remains unchanged: PRESLHY E3.5 gives 62
downstream arc maxima with MG 1.047, VG 1.425 and FAC2 0.839; the independent
FFI/DNV six-arc screen gives MG 1.245, VG 1.373 and FAC2 0.833. FFI Test 6 at
30 m remains a 3.48-fold observed/predicted mismatch under the frozen low-mast
wind. Wind sensitivity and physical sensor geometry do not identify a unique
source or transport correction, so no fitted multiplier was added.

N₂/O₂ mixed-phase closure, particle slip, obstacle wakes, and fully transient
plume storage remain explicit evidence boundaries. They are not silently
promoted to validated physics.

The comparison gate does not create evidence that is not present. It changes
the practical status of the ``HyRAM/PHAST comparison`` item from an informal
claim to an executable protocol: when independently generated predictions are
available under one declared case, the same Hanna metrics can be compared;
until then the report must say that no superiority ranking was performed.

## Verification

- Full suite before the observation-envelope addition: **1126 passed, 139
  skipped**, no failures.
- Latest targeted scope, comparison, gate and observation-envelope tests: **18
  passed**.
- FFI/PRESLHY audit regenerated 197 concentration rows, 40 thermal rows, and
  six FFI arc rows with the frozen metrics above.

Raw third-party workbooks, sensor time series, and original DEGADIS Fortran are
not included in the distribution.
