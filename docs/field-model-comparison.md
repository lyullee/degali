# Fixed-sensor field-model comparison

`degali field-compare` compares two externally generated, fixed-sensor CSV
files. It is for a reproducible DEGALI–SLABx, DEGALI–CFD, or steady–transient
audit; it does not execute the external model, calibrate either model, or turn
sensor points into an interpolated hazard distance.

Start from [the strict JSON template](field-model-comparison-case.example.json)
and use paths relative to that JSON file:

```bash
degali field-compare comparison.json --output new-comparison.json --require-comparable
```

Each CSV must contain, by default, `sensor`, `x_downwind_m`, `y_crosswind_m`,
`height_m`, `averaging_time_s`, and the declared prediction column. Alternative
column names can be declared individually in the JSON. Concentration units are
only `mole_fraction` or `volume_percent`; volume percent is converted to mole
fraction and retained in the report provenance. Every imported CSV receives a
SHA-256 fingerprint, absolute resolved path, row count, and declared column.
For a replay artifact that must not drift after the case is authored, each
model block may also declare `csv_sha256` and `csv_row_count`. They must be
provided together and are checked against the resolved CSV before any model
comparison is run; a changed file is rejected rather than silently acquiring
new provenance. Reports expose this distinction as
`csv_provenance.integrity_pinned`; a normal import/export fingerprint remains
observable but is not presented as an expected-hash verification.
The CSV reader also rejects empty or duplicate (including stripped-name)
header names and rows with either fewer or more values than the header, rather
than allowing a parser overwrite or implicit missing value to change which
prediction column is compared. Declared sensor-ID, coordinate,
prediction, and averaging-time columns must also be distinct single-line names;
role collisions are rejected before reading any prediction row.

The strict JSON reader applies the same fail-safe principle to object keys: an
exact duplicate or a key collision after whitespace/case normalization is
rejected before the comparison schema is mapped. A saved case therefore cannot
silently replace its model path, basis, or CSV fingerprint through JSON
last-value-wins behavior.

An already parsed `FieldModelComparisonCase` can be materialized with
`write_field_model_comparison_case_json()`. The writer preserves each model's
column mapping and comparison basis, verifies that the referenced CSV has not
changed and that its predictions still match the typed model set, writes
relative paths plus fresh hash/row-count pins, and re-reads the emitted JSON
before returning. Existing case files are never overwritten.
Exclusive creation also rejects a competing writer that wins the path after
the preflight check.

After `degali field-compare --output comparison-execution.json`, run
`degali field-verify comparison-execution.json` to recheck the case SHA,
reopen both prediction CSVs, and recompute the report. A changed CSV or edited
execution report returns a failure; successful verification is provenance only
and never promotes model selection or approval.

The two model blocks must separately declare a `comparison_basis`:

- `source_boundary_id`: release source, phase/flash closure, and boundary;
- `weather_id`: wind and stability interval used in the run;
- `sensor_geometry_id`: receptor coordinates and heights;
- `temporal_operator_id` plus optional `averaging_time_s`: snapshot, peak, or
  specified temporal/averaging operation.
- optional `obstacle_representation_id`: the declared obstacle mask/wake-closure
  representation used by the run. If an obstacle is present, both model paths
  must declare the same representation ID before a matched selection impact is
  allowed; one-sided or different IDs remain conditional.

Equal IDs are not evidence by themselves; they are an auditable assertion the
project must substantiate. Different IDs, unequal sensor IDs/coordinates, an
obstacle representation mismatch, or a steady/transient mode mismatch leave
the numerical rows in the report but withhold the model-selection impact.
`--require-comparable` then exits with code 2. When the basis is accepted, the
threshold is required to be a mole fraction in `[0, 1]`, and values outside
that physical range are rejected before rows are classified. The report still
gives only each model's farthest *monitored* threshold exceedance;
it never claims a continuous risk-distance change.

Every comparison also carries a deterministic `gate_codes` list alongside the
human-readable reasons. `matched_basis` marks an accepted operator contract;
`basis_mismatch`, `averaging_operator_unverified`, `execution_conditional`,
`execution_uncertainty_unresolved`, and `temporal_mode_mismatch` identify
conditional holdbacks; `sensor_set_mismatch`, `sensor_geometry_mismatch`, and
`execution_blocked` identify blocked inputs. The decision-facing
`model_selection_impact.gate_codes` additionally reports
`classification_disagreement` or `classification_aligned`, and adds
`model_selection_withheld` whenever the numerical rows cannot support a model
selection conclusion. These are stable automation keys; adjacent reasons
remain the audit explanation.

If material comparison evidence is outside the CSV—for example, a native
source/flash mapping or the origin of an averaging setting—add an optional
`comparison_evidence` object to the case. It has exactly `manifest_path`,
`sha256`, and `qualification`. The parser resolves the manifest relative to the
case file and refuses a changed file. This preserves the qualified evidence in
the execution report; it does not repair a mismatch or make two operators
comparable. A CSV without an averaging-time column must be explicitly declared
with `"averaging_time_column": null`, and its temporal limitation should be
recorded in this qualification. Such a model remains numerically comparable,
but the comparison is conditional and model-selection impact is withheld when
the declared basis includes an averaging time that the CSV cannot verify.

The output is not an independent validation, a ranking against a measured
release, or a design/approval decision. Pair it with aligned field observations
and the existing observation-comparison protocol before making accuracy claims.

The Python comparison records enforce the same boundary when constructed
directly: sensor predictions and comparison rows are immutable tuples, every
prediction/row has its declared typed basis, and a blocked comparison cannot
carry fabricated numerical rows or disagreement counts. These checks prevent a
hand-built report from looking like an accepted model-selection result.

## Exporting a DEGALI field run

An executed `FieldSemiFVScreeningResult` can be converted into the same
fixed-sensor comparison contract with
`field_model_sensor_set_from_screening()`. The caller must choose the trace
operator explicitly: `peak_true`, `peak_indicated`, `final_true`,
`final_indicated`, `time_average_true`, or `time_average_indicated`. The helper
does not clip an indicated value outside `[0, 1]`, does not use a withheld
sensor, and checks a declared comparison averaging time against the sensor
operator before creating predictions.

The resulting model set retains execution provenance, including the source
screening status and selected trace operator. A conditional field execution
therefore keeps a comparison conditional and withholds model-selection impact;
matching basis IDs alone cannot turn it into an accepted comparison.

To materialize a typed result for the CLI/JSON comparison path, use
`write_field_model_sensor_set_csv()`. It writes the fixed sensor geometry,
mole-fraction prediction, and (when declared by the basis) the averaging time,
then returns an equivalent model set with an absolute path, row count, and
SHA-256 fingerprint in `csv_provenance`. The comparison basis IDs are not
hidden in the CSV and must still be supplied explicitly in the comparison case.
The writer creates parent directories but refuses to overwrite an existing
file, so a regenerated artifact cannot silently replace the one already
referenced by a report. The CSV is opened exclusively after its rows are
prepared, closing the same race for concurrent writers.

The same boundary is available for a completed
`FieldSensorSuperposition` through
`field_model_sensor_set_from_superposition()`. Superposition remains
conditional because source-source interaction and three-dimensional cold-cloud
mixing are unresolved, and its provenance is carried into the comparison
report rather than being flattened into an accepted single-source result. Any
unresolved branch bounds are also retained as incomplete uncertainty; a nominal
superposition does not stand in for its deterministic source/geometry/sensor
envelope.
