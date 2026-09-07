# Stage1 limited step-sensitivity diagnostic

Registered after the .02 m candidate and .2 m BASE reproductions completed,
BEFORE running the finer steps. This is a numerical sensitivity report,
not a new physics candidate or changed promotion criterion.

On already selected representative trials10/23 only:

- BASE: reduce the output/interpolation interval from .2 m to .1 m using
  identical hydrogen_jet inputs, REACH=40 m and unchanged solver tolerances.
- CANDIDATE: rerun the existing complete source-ablation path with max step
  .01 m instead of .02 m. Keep all source, model and solver tolerances.
- Compare the same concentration arcs, vertical profiles and41 temperature
  sensors. Report maximum and mean absolute changes, concentration relative
  changes, and metric/decision sensitivity without deleting any adverse row.
- Preserve .02/.2 m historical scores as their own configuration. Fine-step
  scores must not be substituted into only one column of the old paired table.

Two levels do not prove asymptotic convergence, error bounds, or full seven-
trial step independence. If step changes materially alter an interpretation,
report that limitation and do not call the old score validated. No ad hoc
acceptance threshold will be chosen after observing the changes. This bounded
check does not restart thermal/TKE/normal-stress model development.
