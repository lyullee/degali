# Pre-registration: D4.8 source-boundary reconciliation

Date: 2026-09-08

## Question

Do the trial-specific pressure-loss source flows already used by the research
path agree with the independent summary values in PRESLHY D4.8 closely enough
to retain them for the next residual decomposition?

## Frozen comparison

Use the existing, unchanged
`reference/preslhy/measured_pipe_source_2026-09-05.json` values. Compare only
the trials for which D4.8 Appendix 1 states a numerical high-pressure value:

- Trial 10: approximately 298 g/s;
- Trial 11: approximately 265 g/s;
- Trial 12: range 90--100 g/s.

Trial 10 and 11 pass the document-consistency screen when the pressure-loss
value differs from the stated value by no more than 10%. Trial 12 passes when
it lies within the stated range; it is marked `near_range` rather than failed
when it is no more than 10% outside the nearest bound. The tolerances represent
the report's approximate precision and are not model-fitting tolerances.

## Decision rule

- Retain the trial-specific pressure-loss value when the screen is `pass` or
  `near_range`; do not replace a time-window calculation with a rounded report
  summary.
- Stop before downstream recomputation if Trial 10 or 11 fails.
- Record Trial 12 separately because a range comparison is not equivalent to
  an exact-value comparison.
- Do not infer values for Trials 22--25 solely from geometrically similar
  Trials 10--12.

No model coefficient, observation filter, or promotion threshold may change
as a result of this audit.

