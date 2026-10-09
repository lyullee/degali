# Joint source–sensor operational envelope

When a field case contains both bounded physical inputs and supplemental
detector calibration bounds, use:

```python
from degali.addons import (
    field_operational_source_sensor_envelope_report,
    run_field_operational_source_sensor_envelope,
)

result = run_field_operational_source_sensor_envelope(
    request,
    max_cases=64,
    include_refinement=True,
    refinement_factors=(1, 2),
    relative_tolerance=0.05,
)
report = field_operational_source_sensor_envelope_report(result)
```

The physical source/weather/surface/geometry envelope is solved with
supplemental detector calibration fixed at nominal values. Each physical
corner is then crossed with every declared supplemental calibration corner;
the unchanged true receptor traces are post-processed with the exact
response/gain/bias values. Numerical refinement is attached to each physical
corner and reused across its calibration corners.

The aggregate is allowed only when every Cartesian corner passes completion,
refinement, in-plane receptor and resolved-uncertainty gates. The helper does
not infer a source history or phase route, and it does not make a probability
statement. The strict CLI equivalent is
`degali field-screen --joint-source-sensor-envelope`; it cannot be combined
with historian or phase-routing envelopes.

The report also contains `deterministic_sensor_envelope`, which gives each
detector's peak/final/time-average true and indicated extrema over completed
corners. Missing or withheld traces are counted with their reasons; they are
never replaced by zero and the extrema are not confidence intervals.
