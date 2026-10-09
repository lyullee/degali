# Field stability alternatives

The local semi-FV field operator has no embedded atmospheric-stability
correlation. A `stability` label alone is therefore only metadata unless an
explicit `stability_mixing_closure` supplies the selected scalar diffusivity.

When the observed or forecast stability class is ambiguous, a strict
`degali.field-screening-input.v1` case can propagate each declared class
through `degali field-screen --uncertainty-envelope`. Add both objects at the
top level of the case:

```json
"stability_mixing_closure": {
  "model_id": "site-reviewed-effective-diffusivity-v1",
  "evidence_id": "met-mixing-review-2026-10-05",
  "diffusivity_m2_s": {
    "neutral": 0.50,
    "stable": 0.20,
    "unstable": 0.85
  }
},
"stability_alternatives": {
  "alternatives": ["stable", "unstable"],
  "evidence_id": "met-stability-classification-2026-10-05"
}
```

The nominal class stays at `scenario.weather.stability`; it is automatically
included once, followed by each alternative. Every class must have a positive
declared diffusivity. There are no assumed Pasquill-Gifford, RANS, wake, or
probability-weight closures. Missing evidence or a missing class coefficient
is rejected before calculation.

The nominal run remains withheld for unresolved categorical stability. The
operational envelope marks that uncertainty resolved only after it has run and
refined every source, wind speed/direction, detector, and stability corner.
The same logic applies to the quality-gated historian joint envelope. This is
still a deterministic sensitivity set—not a probability interval, validated
turbulence model, or design/approval basis.
