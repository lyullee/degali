# TI021 R12 resolved-inner 60 s qualification result

Run: `outputs/thorney-ti021-inner-r12-fine-60s-dx2p5m-dt0p0625-zero-thickness-2026-10-13-v2`

Independent score: `outputs/thorney-ti021-inner-r12-fine-60s-score-2026-10-13-v2/ti021-inner-transient-score.json` and its PNG plot.

## Scope and execution

- Observation-blind resolved-inner run; the configuration records an empty list of observation values used during setup.
- 60 s simulated, `dx=2.5 m`, `dt=0.0625 s`, zero-thickness curved fence representation; 80×80×16 cells.
- Final checkpoint is `finalized` at 60 s after six 10 s segments. The progress journal records 1010.0 s cumulative wall time (about 16 min 50 s).
- Last-segment maximum global mass residual: `1.40e-11 kg/s`; momentum balance residual: `1.79e-9 kg m/s`; maximum relative species-mass residual: `2.13e-16`; zero nonconverged corrector steps.
- This is a completed qualification calculation, but it does not span the TI021 observation/dose horizon and is not grid-converged. It cannot validate hydrogen, and the 2.5 m cells alias the 0.4 m and 2.4 m receptor heights.

## Provisional comparison

The score uses 60 TI021 receptors (53 with positive observed dose) and six published radial-arc peak summaries. The run predicts threshold arrival at 43/52 receptors where arrival times are available; arrival-time MAE among detected pairs is 23.23 s. For the six arcs, the 0.6 s-averaged peak-concentration MAPE is 16.78%, while peak-time MAE is 87.10 s. Model peak times are 17.1–24.1 s versus observed 86.5–133.5 s: a clear indication that the early transient does not reproduce the observed peak timing.

The 60 s truncated dose comparison is deliberately marked provisional: sensor-dose MAE is 269.63 vol%·s and normalized mean bias is -90.05%, but the observed dose covers a much longer event, so these are not final dose errors. The 4 vol% arc comparison has six true negatives only; it provides no sensitivity evidence.

Decision remains `withheld_as_horizon_limited`; no fitting, promotion, or quantitative model-validation claim is allowed. At the time of this score, the matched frozen 1.25 m obstacle-impact/grid audit had not yet been re-evaluated.

## Scoring protocol correction

The scorer previously conflated “scoring after the solver completed” with “the simulated horizon covers the complete experiment.” Those are independent conditions. The protocol now records `scoring_opened_after_completed_run=true` for an accepted finalized qualification run and reports `observation_horizon_complete=false` separately. This changes protocol metadata only; it does not change any score or its withheld status.

## Full-horizon result already available

A separate, earlier 600 s run at the same nominal 2.5 m grid was also re-scored after completion: `outputs/thorney-ti021-coarse-600s-dx2p5m-dt0p0625-checkpoint1s-tol1e3-2026-10-11`, with score in `outputs/thorney-ti021-coarse-600s-score-2026-10-13-v2/ti021-inner-transient-score.json`. Its checkpoint is finalized at 600 s and its progress journal reports 7127.5 s wall time (about 1 h 58 min 48 s). This closes the horizon question for that older run, not for the newer zero-thickness setup.

The full-horizon six-arc peak MAPE is 52.09%; all six predicted peaks (0.52–0.72 vol%) are below the observations (1.16–1.47 vol%). Peak-time MAE is 60.65 s, with predicted peaks at 42.7–49.1 s versus observed peaks at 86.5–133.5 s. Arrival-time MAE is 22.55 s. Six-arc dose MAPE remains 79.48% and the 60-sensor dose MAE is 265.93 vol%·s, despite the full 600 s horizon. Thus the poor dose/peak result is not explained solely by the new run's 60 s truncation.

Do not interpret the 60 s and 600 s runs as a controlled before/after comparison: the full run used `corrector_tolerance=1e-3`, whereas the newer zero-thickness run used `1e-4`, and the older run record predates explicit wall-representation provenance. The full run is a completed, no-fit diagnostic at 2.5 m, but remains withheld because the near-ground receptor heights are aliased and the 1.25 m impact-window member is required for a grid audit. It also does not validate hydrogen; TI021 is an R12 dense-gas obstacle surrogate.

Subsequent matched 5/2.5/1.25 m, zero-surface, `dt=0.02 s` 20 s runs now have a controlled grid-sensitivity audit; peak/dose and 4 vol% region gates fail. Details are in `docs/thorney-ti021-zero-surface-three-grid-audit-2026-10-13.md`. This closes the missing-grid-member issue for the initial 20 s impact window but not the full 600 s event or the hydrogen-obstacle validation objective.
