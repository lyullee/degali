# Operational detector-calibration envelope

`run_field_sensor_array_uncertainty_envelope()` evaluates detector calibration
corners, but it is a physical/audit envelope rather than an operational gate.
Use `run_field_operational_sensor_array_envelope()` when each calibration
corner must also carry refinement and fail-safe screening evidence:

```python
from degali.addons import (
    field_operational_sensor_array_envelope_report,
    run_field_operational_sensor_array_envelope,
)

result = run_field_operational_sensor_array_envelope(
    request,
    include_refinement=True,
    refinement_factors=(1, 2),
    relative_tolerance=0.05,
)
report = field_operational_sensor_array_envelope_report(result)
```

The scalar transport is solved once and reused because response time, gain and
bias are observation operators. Each exact calibration corner receives its own
sensor traces, refinement reference and `FieldOperationalScreeningDecision`.
The aggregate is allowed only when every corner passes the ordinary completion,
in-plane receptor, refinement and resolved-uncertainty gates.
Typed envelope cases also cross-check the decision against the physical
applicability attached to that case: `screening_allowed` requires `accepted`
applicability and `conditional_allowed` requires `conditional` applicability.

This helper resolves detector calibration uncertainty only. Source rate,
weather, stability, obstacle geometry, phase routing and transport uncertainty
remain unresolved unless their dedicated envelopes are run; those conditions
therefore remain `withheld` rather than being hidden by the calibration study.
Conditional physics still requires explicit `allow_conditional=True` and never
grants a design-basis or approval decision.
