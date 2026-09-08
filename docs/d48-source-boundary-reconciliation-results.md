# D4.8 source-boundary reconciliation results

Date: 2026-09-08

## Outcome

The trial-specific pressure-loss source values are retained for the next
downstream residual decomposition. They are not replaced with the rounded
D4.8 summaries.

| trial | existing pressure-loss flow | D4.8 Appendix 1 | comparison |
|---:|---:|---:|---:|
| 10 | 284.340 g/s | approximately 298 g/s | -4.584%, pass |
| 11 | 265.275 g/s | approximately 265 g/s | +0.104%, pass |
| 12 | 105.544 g/s | approximately 90--100 g/s | +5.544% above the upper bound, near range |

Trial 10 and 11 meet the pre-registered 10% document-consistency screen.
Trial 12 lies just outside the stated range but inside its pre-registered
near-range band. This is recorded as source uncertainty and not silently
converted into a calibration.

The report explains that Test 10's mass-flow meter drive gain is saturated by
two-phase flow and reconstructs its flow from the pressure-drop ratio to Test
11. The existing pressure-loss calculation uses the trial-specific window and
measured state, so retaining 284.340 g/s preserves more information than
replacing it with the report's approximate 298 g/s.

No D4.8 numerical summary was extrapolated to Trials 22--25. No model
coefficient, observation population, or acceptance threshold changed.

## Verification

The audit is implemented in `tools/audit_d48_source_boundary.py`. Four tests
cover exact-value acceptance, range/near-range classification, refusal to
extrapolate to other trials, and rejection of nonphysical flow. All four pass.

The local result, including SHA-256 hashes of the source table, D4.8 PDF,
protocol and audit program, is stored at
`reference/preslhy/d48_source_boundary_reconciliation_2026-09-08.json` and is
excluded from distribution with the other third-party research material.

## Next model action

Use the unchanged pressure-loss source for Trials 10 and 23. Reuse the sealed
trajectories and existing energy/cold-envelope audits before running any new
field calculation. The next new calculation should add information not already
contained in those audits: a bounded transported-energy state or a matched
three-dimensional momentum handoff. It must preserve the source-flow
uncertainty as a separate sensitivity and must not fit a coefficient to the
same plume-height observations used for scoring.

