# Pre-registration: PRESLHY transient source/temperature coupling

Date: 2026-09-05

## Question

The steady validation compares one sustained-window source rate with a
temperature extreme from a release whose measured flow changes substantially.
Trial 23, for example, falls from about 185 to 75 g/s during the retained
window.  Because the Flexlogger mass-flow and converted thermocouples share a
one-second clock, this audit asks whether the remaining cold-core discrepancy
is primarily a missing time-varying source rather than a missing phase model.

## Frozen data and calculation

- Use only the original trial-10 and trial-23 workbooks and exactly the release
  windows already frozen by `read_trial`.
- Read `MFM1_Mass_Flow_Rate` and the converted centreline `TC...C` channels
  from `Flexlogger`.  Do not use the separate Flowmeter clock, pipework
  thermocouples, off-centre sensors or concentration sampling lines.
- For every available centreline channel report Spearman rank correlation
  between temperature at second `t` and source flow at `t-lag`, for every
  integer lag from 0 through 5 s.  Report all six lags; do not select only the
  best lag.
- Also split the frozen source-flow samples at their median and report the
  temperature median in the lower-flow and upper-flow halves at each lag.
  The diagnostic contrast is lower-flow minus upper-flow temperature: a
  positive value means the high-flow release carries a colder cloud.
- No correlation, lag, source multiplier or temperature offset is fitted into
  the dispersion model in this audit.

## Frozen interpretation

- Source unsteadiness is a required next model state if at least five of the
  seven trial-23 centreline channels have negative Spearman correlation at
  every lag from 0 through 3 s, and their median lower-minus-upper flow
  temperature contrast is at least +10 K for one of those common lags.
- If the sign is inconsistent or the common contrast remains below 10 K,
  source variation cannot explain the coherent cold core.  Proceed directly
  to a separate liquid-H2/vapour state and finite-rate interphase transfer.
- Passing this diagnostic does not establish that source unsteadiness is
  sufficient.  It only requires a subsequent causal, travel-time-aware model
  comparison against the complete temperature time series before any phase
  closure is changed.

## Results

Trial 23 has the expected sign: all eight centreline thermocouples have
negative flow/temperature Spearman correlation at every lag from 0 through
3 s.  The median correlation is modest (`-0.194`, `-0.139`, `-0.235`,
`-0.272`) and the corresponding median low-flow-minus-high-flow temperature
contrast is only +2.10, +1.02, +3.19 and +6.00 K.  It therefore fails the
pre-registered 10 K sufficiency screen.

Trial 10 provides an independent warning against forcing the interpretation.
Seven of eight channels have the opposite sign at every tested lag; the
median contrast is about -12 to -14 K.  Its pipework and near-field continued
to cool during the retained interval, so elapsed release history dominates a
simple flow correlation.

The measured transient source should eventually be supported for time-series
prediction, but it is not promoted as the explanation of the 30--80 K
centreline residual.  The next isolated mechanism is the hydrogen phase
state: measured or bounded nozzle vapour quality, separate liquid/vapour
enthalpy, and finite-rate interphase transfer.  Full per-channel results for
all six lags are in
`reference/preslhy/transient_source_temperature_2026-09-05.json`.
