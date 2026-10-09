# Conditional field-screening review record

`degali field-screen --allow-conditional` is an accountable opt-in for
screening use only. It is not a switch that converts a conditional model into
a validated model. A strict case must contain a `conditional_review` object
before the CLI will apply that option:

```json
{
  "conditional_review": {
    "review_id": "field-review-A",
    "reviewer_id": "process-safety-reviewer-01",
    "reviewer_role": "process_safety_engineer",
    "reviewed_at_utc": "2026-10-05T00:00:00Z",
    "evidence_id": "field-screening-review-record-A",
    "allowed_scope": "conditional_screening_only"
  }
}
```

All fields are mandatory. The time must be an ISO-8601 UTC timestamp ending in
`Z`, and `allowed_scope` has exactly one accepted value:
`conditional_screening_only`. Values such as `design_basis`, `approval` or
`unspecified` are rejected.

The execution JSON preserves the complete review record and adds
`applied_to_this_execution`. Supplying the object without
`--allow-conditional` records the review but does not activate it. Supplying
the flag without the object fails before transport starts.

## Review procedure

Before issuing the referenced review record, the accountable reviewer should
confirm that:

1. source, coordinate, weather, surface, obstacle and detector evidence IDs
   identify the records actually used by the case;
2. the matching uncertainty envelope covers every declared bound and
   categorical stability alternative relevant to the selected physical path;
3. grid/time refinement has run and meets the declared tolerance;
4. every decision receptor is represented in the supported wind plane;
5. every declared obstacle is represented, or the case remains withheld until
   an appropriate three-dimensional assessment is available; and
6. unresolved phase, pool-launch, cold-cloud, in-flight-droplet, wake and
   validation limitations in the report have been reviewed for the intended
   screening question.

The software independently enforces the non-negotiable completion,
refinement, uncertainty, receptor and obstacle gates. A review record cannot
waive them. Even an allowed conditional result always retains
`design_basis_allowed=false` and `approval_allowed=false`.

For `export_field_screening_batch()`, pass a
`conditional_review_authorizations` mapping keyed by the exact case labels.
When `allow_conditional_operational_screening=True`, every batch case must have
one authorization or the export fails before creating the output directory.
Manifest v3 preserves each review record and whether it was applied.

The example identifiers are synthetic. Replace them with the organisation's
controlled review record and role identifiers; do not treat the example as an
authorization.
