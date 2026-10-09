# Matched field-validation input

`degali field-validate` scores one declared model sensor CSV against one
already-selected observed-sensor summary CSV. It is not a historian importer
and does not infer a release source, phase state, obstacle wake or sensor
response operator.

The top-level schema is `degali.field-validation-input.v1`:

```json
{
  "schema": "degali.field-validation-input.v1",
  "model": {
    "model_id": "DEGALI",
    "temporal_mode": "transient",
    "csv_path": "model.csv",
    "prediction_column": "prediction",
    "concentration_unit": "mole_fraction",
    "comparison_basis": {
      "source_boundary_id": "source-A",
      "weather_id": "weather-A",
      "sensor_geometry_id": "detectors-A",
      "temporal_operator_id": "mean-60s",
      "averaging_time_s": 60.0,
      "obstacle_representation_id": "mask-A"
    }
  },
  "dataset": {
    "csv_path": "observed.csv",
    "temporal_mode": "transient",
    "comparison_basis": {
      "source_boundary_id": "source-A",
      "weather_id": "weather-A",
      "sensor_geometry_id": "detectors-A",
      "temporal_operator_id": "mean-60s",
      "averaging_time_s": 60.0,
      "obstacle_representation_id": "mask-A"
    }
  },
  "validation_evidence": {
    "dataset_id": "lh2-obstacle-trial-A",
    "path": "observed.csv",
    "sha256": "<64 lowercase hexadecimal characters>",
    "row_count": 2,
    "source_boundary_id": "source-A",
    "weather_id": "weather-A",
    "obstacle_geometry_id": "obstacle-A",
    "receptor_geometry_id": "detectors-A",
    "temporal_operator_id": "mean-60s",
    "common_clock_id": "clock-A",
    "scope": "lh2_obstacle_transport"
  }
}
```

The observed CSV must contain exactly one summary row per sensor with
`sensor`, `x_downwind_m`, `y_crosswind_m`, `height_m`,
`observed_mole_fraction`, `observation_time_s`, `averaging_time_s`,
`common_clock_id` and `obstacle_geometry_id`. Every row must repeat the
evidence clock and obstacle IDs. The SHA-256 and row count must match the
evidence object exactly; a mismatch is rejected before scoring.

An optional `observation_kind_column` can be declared in `dataset` when the
CSV contains an `observation_kind` column. Its only accepted values are
`exact` and `lower_bound`. A lower-bound value means that the reported
concentration is a detection or censoring threshold, not an exact reading;
the model must be at least that large. Lower-bound rows are excluded from
MAE/RMSE/bias and are reported through one-sided satisfaction and deficit
metrics. A dataset containing any lower-bound rows is never `qualified`, even
when every bound is satisfied; a violated or unmatched bound is `withheld`.
If the column is omitted, all rows retain the original exact-observation
semantics.

The dataset comparison basis is also an identity boundary: its source,
weather, receptor-geometry and temporal-operator IDs must equal the matching
IDs in `validation_evidence`. A fingerprinted CSV cannot be scored under a
different source or observation operator. The physical `obstacle_geometry_id`
is kept separate from the model's `obstacle_representation_id` (for example a
site geometry record may be compared with a model mask), so that distinction
is evaluated by the model-scoring gate rather than by string equality here.

The score is restricted to matched monitored sensors. A `qualified` result is
not a spatial hazard-distance, design-basis or approval decision. Free-field
evidence cannot qualify a model basis that declares an obstacle
representation; obstacle evidence cannot qualify a model that omits that
representation. The JSON score also reports `observation_counts`, the
`lower_bound_constraint` object, and an
`exact_classification_disagreement_count` so that classification differences
from censored rows are not mistaken for exact-observation accuracy.

An in-memory `FieldValidationCase` can be materialized with
`write_field_validation_case_json()`. The writer verifies that both the model
prediction CSV and observed CSV still match their typed predictions/observations
and evidence digest, preserves custom column mappings and relative paths, then
re-reads the emitted strict JSON. Existing case files are never overwritten;
the JSON is created exclusively after both CSV fingerprints are validated, so
a competing writer cannot replace the pinned validation case.
The observed CSV reader also rejects duplicate/empty (including stripped-name)
headers, truncated or extra row values, and logical column-role collisions
before any score is computed. The model-prediction CSV uses the same exact
row-width and header checks.

The strict validation JSON reader rejects duplicate or whitespace/case-normalized
object keys before loading model and observed CSV provenance, so a saved case
cannot silently replace its evidence identity or hash.

After `degali field-validate --output validation-execution.json`,
`degali field-verify validation-execution.json` rechecks the case and observed/
model CSV provenance, then recomputes the score. It reports integrity only and
does not turn a verified score into a design-basis or approval decision.
