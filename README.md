# DEGALI — Dense Gas Dispersion for Liquid Hydrogen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22646258.svg)](https://doi.org/10.5281/zenodo.22646258)
[![PyPI](https://img.shields.io/pypi/v/degali.svg)](https://pypi.org/project/degali/)

DEGALI stands for **Dense Gas Dispersion for Liquid Hydrogen**. It is a modern
Python model built from a verified reimplementation of the US EPA DEGADIS 2.1
dense-gas dispersion model and extended for cryogenic hydrogen releases.

> **Alpha research software.** The DEGADIS 2.1 compatibility path is strongly
> regression-tested. The liquid-hydrogen extensions are suitable for research,
> scenario comparison and sensitivity studies, but are not independently
> certified for regulatory separation distances or safety-critical design.

한국어 요약: 현재 공개본은 원본 DEGADIS 재현 경로와 액화수소 연구 확장을
함께 제공합니다. PRESLHY 수평 야외제트 범위에서 정량 비교를 마쳤지만,
미완성 TKE/압력/열폭 폐쇄와 제한된 독립 검증 때문에 설비 인허가 판단의
단독 근거로 사용하면 안 됩니다. 자세한 판정은
[1차 결과](docs/stage1-results-2026-09-06.md)를 참고하십시오.

처음 설치하거나 실제 계산을 시작하려면 상세 사용 가이드를 참고하십시오:
**[한국어](USER_GUIDE_KO.md) · [English](USER_GUIDE.md)**. LH₂ 간편
평가, 결과 해석, CSV 저장, 민감도 분석, 기존 DEGADIS 입력 덱과 CLI 사용법을
단계별로 설명합니다.

## What is included

- Python ports of all six DEGADIS 2.1 programs. The independently obtained
  Fortran oracle is retained locally and is not redistributed.
- Steady, transient, receptor-dose and jet-to-ground workflows.
- Legacy and modern thermodynamic backends.
- Opt-in LH2 source flashing, buoyant trajectory, air condensation/freezing,
  component enthalpy, ground interaction and thermal-profile research paths.
- Experimental finite-TKE, independent thermal-width and Reynolds-stress
  operators. These are deliberately marked as research-only until their
  physical closure inputs and downstream field performance are validated.

## Current evidence

The port reproduces the five EPA reference cases and was checked locally
against a source-built Fortran implementation. The first frozen LH2 assessment uses
seven horizontal outdoor PRESLHY trials, 38 concentration sections, 17
vertical profiles and 41 temperature sensors.

| Quantity | Current provided path |
|---|---:|
| Concentration MG (observed/predicted; ideal 1) | 1.102 |
| Concentration VG (ideal 1) | 1.233 |
| Concentration FAC2 | 36/38 |
| Vertical width, predicted/measured | 1.091 |
| Plume-centre MAE | 0.054 m |
| Minimum-temperature MAE | 27.25 K |

An opt-in thermal-profile candidate lowers minimum-temperature MAE to 16.47 K
but worsens other acceptance metrics, so it is not the default. These samples
are correlated observations from one campaign, not 38 or 41 independent
validation experiments.

### Current IJHE manuscript baseline

The manuscript freeze used for the IJHE submission track is reported separately
from the historical first-pass table above. On nine momentum-dominated
PRESLHY E3.5 releases (62 arc maxima), the pre-specified peak operator gives
MG = 1.047, VG = 1.425 and FAC2 = 0.839. A synchronised 20 s mean operator
on seven trials gives MG = 0.544, VG = 4.601 and FAC2 = 0.761, while the
independent six-arc FFI/DNV screen gives MG = 1.245, VG = 1.373 and FAC2 =
0.833. These are separate evidence lanes, not one pooled accuracy score; the
low-wind FFI Test 6 residual remains an explicit application limit.

The accompanying 3-D transient solver is a numerical research extension. Its
grid-convergence and mass-conservation results are documented in the
[IJHE submission-readiness record](docs/ijhe-submission-readiness.md) and do
not constitute measured field validation.

Run `python tools/audit_ijhe_submission.py` before circulating a manuscript
snapshot. It checks the claim boundaries, numerical-verification artifact and
public-snapshot safety, while retaining an explicit warning until a matched
time-resolved field event is available.

The draft upload bundle is [the manuscript](MANUSCRIPT_LH2_VALIDATION.md),
[Supplementary Information](IJHE_SUPPLEMENTARY_INFORMATION.md),
[Highlights](IJHE_HIGHLIGHTS.md) and a
[cover-letter draft](IJHE_COVER_LETTER_DRAFT.md).
The [graphical abstract](IJHE_GRAPHICAL_ABSTRACT.svg) and
[figure/table register](docs/ijhe-figure-table-register.md) keep the upload
artifacts reproducible.

An author-neutral assembled draft is generated with
`tools/build_ijhe_upload_bundle.py` in
`outputs/IJHE_UPLOAD_DRAFT_LATEST/`. Its `BUNDLE_STATUS.md` and
`bundle-manifest.json` keep the remaining author and external-evidence gates
visible; the directory must not be treated as a final upload until those gates
are closed.

The separately audited [open-channel H₂ timing boundary](docs/open-channel-h2-boundary.md)
is retained as supporting evidence only; it is not promoted to atmospheric
LH₂ validation because its weather, geometry and common-clock evidence are
incomplete.

The research-only six-flux thermal-moment marcher now passes its preregistered
numerical extension across PRESLHY Trials 11, 12, 22, 23, 24 and 25: all coarse
and refined paths reach their frozen stations, with maximum independent balance
residual `5.05e-11` and maximum terminal flux difference `1.57e-4`. This is a
conservation, solvability and step-convergence result—not a new field-accuracy
score or safety-design validation. See the
[thermal-moment extension results](docs/thermal-moment-flux-extension-results.md).

A subsequent identical-sensor audit for Trials 10 and 23 reaches 4 m with
strong conservation and step convergence under an opt-in complete-radial-core
constitutive check. It slightly improves thermal spread but worsens centre
temperature-deficit amplitude, so the field-acceptance gate fails and neither
the candidate nor the new domain option is promoted. See the
[thermal-moment observation result](docs/thermal-moment-observation-results.md).

See [model status](docs/stage1-results-2026-09-06.md),
[claim grading](docs/claim-grading.md), and
[data/reproduction notes](docs/DATA_AND_REPRODUCTION.md). The current
availability and limitations of the LH2 measurements are summarized in the
[LH2 data-readiness register](docs/lh2-data-readiness.md).

## Field-deployment extension (preliminary)

An opt-in local 2-D wind-plane, conservative semi-finite-volume screening
path now makes the release boundary, flash state, weather, surface, sensors,
obstacle geometry and deterministic uncertainty bounds explicit. A selected
historian event can enter only through an SI CSV map with calibrated P/T/mass
flow channels, a separate phase/flash-state evidence identifier and a retained
file fingerprint; an ordinary pressure trend or nominal volume flow is not
converted into a leak source. Long historian windows can reuse a bounded LH2
property table instead of calling CoolProp for every interval.

This path is suitable only for conditional operational screening and research.
It is not a CFD solver, a building-wake closure, a design-basis calculation or
an approval decision; local LH2 obstacle validation remains outstanding. See
the [deployment roadmap](docs/field-deployment-roadmap-2026-10-05.md) and
[connected-project evidence review](docs/field-evidence-readiness-2026-10-05.md).

For a declared nominal site case, start from the synthetic
[field-case template](docs/field-screening-case.example.json), replace every
evidence ID and value under engineering control, then run:

```bash
degali field-screen field-case.json --output new-field-execution.json
```

The command runs its grid/time refinement by default and records its input
SHA-256, numerical report and fail-safe disposition. `--allow-conditional` is
an explicit accountable-review opt-in for screening only after relevant
source/weather/detector uncertainty has been propagated; it never enables a
design-basis or approval result. The flag requires the strict case's
[`conditional_review` record](docs/field-conditional-review.md), which is
retained in the execution JSON. The nominal template does not include a
historian or a post-flash schedule. A selected SI historian CSV may be added
only through the [quality-gated measured-history mapping](docs/field-measured-history-input.md);
post-flash schedules still retain their dedicated API/provenance gate and must
carry an explicit non-default `source_id` when attached through the Python
field request API.

For a pressure/temperature-only historian, the programmatic
`PressureDrivenMeasuredHistory` path is the only supported source derivation:
`direct_vapour_schedule_from_pressure_driven_history()` recomputes an explicit
orifice boundary interval by interval, and
`run_field_pressure_driven_history_envelope()` propagates its coherent source
corners into the semi-FV screen. Use
`run_field_joint_pressure_driven_history_envelope()` when ambient, weather and
detector corners must be propagated with that source. A pressure trend is never interpreted as a
leak rate; declared opening area/discharge coefficient (including any bounds
that must be crossed by the joint envelope) and source location/direction
corners are fixed before wind-plane projection. Independent event, phase,
clock and calibration evidence remain mandatory.
The strict case-file equivalent is `pressure_driven_history`; its CSV map
contains P/T (and optional liquid fraction) but no mass-flow column, so the
zero source placeholder cannot be mistaken for a measured rate. `field-screen`
and `field-verify` retain the pressure-driven history kind and its CSV digest
through the joint uncertainty/refinement envelope.

Ambient temperature, absolute pressure and air-density bounds may also be
declared at the top level as `ambient_temperature_uncertainty_k`,
`ambient_pressure_uncertainty_pa` and `ambient_air_density_uncertainty_kg_m3`.
They are deterministic corners, not probability intervals; a nominal result
with unresolved ambient bounds remains withheld by the operational gate. The
phase-routing/pool transport adapter carries the same ambient corners through
the phase ledger and clears them only on the corresponding field corner;
measured-history operational runs remain withheld until their time-aligned
history and ambient boundaries are both propagated.
When an ambient-pressure bound is present, the reusable LH2 saturation table
also spans both endpoint saturation temperatures before any corner is run.
Likewise, a fingerprinted validation dataset is evidence identity only: the
field gate remains conditional until a matched observed-vs-model validation
score is generated.
Direct ambient uncertainty objects enforce `K`, `Pa`, and `kg/m3` units before
their corners enter flash, transport, or detector conversion, matching the
strict case-file boundary.
The joint measured-history envelope carries the same ambient temperature and
pressure corners through each history re-flash and applies ambient air-density
corners to detector conversion; a completed joint case clears only the ambient
corner it actually evaluated.

When the case contains relevant bounded source, weather or primary-detector
inputs, use `--uncertainty-envelope`. It runs every deterministic corner and
its refinement, then allows screening only if every corner passes the same
gate. With a quality-gated `measured_history` object, the same flag instead
runs the joint P/T/flow/phase-history, weather and detector envelope:

```bash
degali field-screen field-case.json --uncertainty-envelope --allow-conditional \
  --output new-field-envelope.json
```

The envelope is still a sensitivity set rather than a confidence interval.
For a finite transient release, `scenario.source.duration_s` may itself be a
bounded object, or a numeric duration may be paired with
`duration_uncertainty`; the nominal duration and every lower/upper timing
corner are retained in the source contract. The transport duration must cover
the upper declared release duration. No timing probability or unrecorded
source schedule is inferred.
The source `location_m` may likewise use three SI bounded coordinate objects;
their Cartesian lower/upper corners are propagated through the wind-plane
projection and transport case. Mixed numeric/bounded coordinate arrays are
rejected, and this remains a deterministic geometry sensitivity rather than a
probability distribution.
The same three-component bounded form is available for `direction_m`; every
corner must remain a non-zero vector and is preserved for jet/phase-routing
handoff calculations. The direct scalar screen records these direction
corners even where its post-flash source closure does not use momentum.
An explicit base case may also declare `distributed_vapour_sources` for
already-atmospheric pool/droplet or other ledger branches. Each source carries
its global position, a unit-labelled vertical scalar width, a finite
zero-endpoint rate schedule and evidence ID; bounded position/width values are
expanded into separate source-shape corners and are not treated as a
probability distribution. The operational gate withholds a nominal result
until those corners are propagated.
For categorical meteorological stability uncertainty, declare evidence-backed
alternatives and a diffusivity for every class; the local solver does not
invent a turbulence correlation. The [stability-alternative protocol](docs/field-stability-alternatives.md)
adds those cases to both nominal and historian joint envelopes.

`scenario.weather.direction_deg` is a circular degree interval, rather than a
linear number range. A wind record crossing north can therefore declare
`lower: 350`, `nominal: 0`, `upper: 10`; the envelope evaluates 350° and 10°
and records the short 20° sector. It does not reinterpret that input as a
340° sweep through south or average a changing wind into one wind plane.

Obstacle footprints may be global-axis-aligned or explicitly rotated; both are
reprojected for each wind-direction corner and retained in the report. Linear
dimensions (and an oriented bearing when needed) may be declared as explicit
bounded geometry values with an `uncertainty_evidence_id`. The field envelope
then recomputes every obstacle-mask corner, while the nominal operational
screen remains withheld until those corners are resolved. See the [obstacle
input contract](docs/field-obstacle-input.md) for the exact JSON forms,
required site-coordinate reference and the limits of the conservative 2-D
mask.

`lh2_validation_available: true` is not a free-standing claim. Strict case
files must also carry `validation_evidence` with a SHA-256 fingerprint,
positive row count, common-clock ID, source/weather/obstacle/receptor geometry
IDs and temporal-operator ID. A free-field evidence record cannot qualify an
obstacle case; an unbacked validation flag is rejected before transport.

For an actual fixed-sensor validation replay, `degali field-validate` accepts
the separate `degali.field-validation-input.v1` contract. It checks every
observed row's clock/obstacle IDs and file fingerprint before calculating
matched-sensor MAE, RMSE and bias. See
[field-validation-input.md](docs/field-validation-input.md); no unmeasured
hazard distance is interpolated. If an observed CSV declares an
`observation_kind` column, `lower_bound` rows are treated as censored
one-sided constraints, excluded from symmetric accuracy metrics, and prevent
`qualified` status even when satisfied.

`degali field-audit <directory>` scans an external CSV/JSON/XLSX directory read-only
for the five evidence channels needed to assemble that validation case. It is
readiness triage only: even `candidate_complete` keeps
`promotion_allowed=false` and cannot create `FieldValidationEvidence`; a scan
limit sets `scan_complete=false` and also prevents promotion. The saved record
includes candidate-file `diagnostic_counts` for near-miss categories and
recomputes them on read-back. See
[field-evidence-audit.md](docs/field-evidence-audit.md).

After an accountable reviewer reconciles the five files to one event and common
clock, `FieldEvidenceManifest.from_audit()` can record the selected paths,
SHA-256 digests and explicit source/weather/geometry IDs. New packages should
also provide `sensor_set_id`, a SHA-pinned `sensor_registry_artifact` (sensor
calibration/response metadata), and an accountable `operator_id`; missing
metadata is reported as `evidence_readiness=conditional`. The manifest
rechecks its files and remains `promotion_allowed=false`; see
[field-evidence-manifest.md](docs/field-evidence-manifest.md).
The package layout and minimum CSV/JSON fields are listed in
[field-evidence-intake-template.md](docs/field-evidence-intake-template.md).
`field-validation-input.v1` may carry the same manifest as
`validation_manifest`; the parser then rejects any path, digest or explicit-ID
mismatch between the manifest and `validation_evidence`.
The equivalent CLI boundary is
`degali field-evidence-manifest-create <audit-execution.json>` with five
repeated `--selected-path channel=path` options and explicit event/identity
metadata, plus `--sensor-set-id`, `--operator-id` and optionally
`--sensor-registry-path`; it refuses incomplete audits and existing output
paths. Verification emits machine-readable `conditional` or `withheld`
readiness instead of silently treating missing evidence as valid.

For the public FFI/Spadeadam horizontal-release lane, use
`degali ffi-source-state`. It crosses explicitly declared source, wind and
ambient endpoints through the unchanged free-field jet and evaluates every
reported sensor coordinate. The output is a deterministic residual envelope,
not a confidence interval or site-acceptance record:

```bash
degali ffi-source-state --test 6 --radius 30 \
  --rate-bounds 0.82,0.833,0.85 \
  --wind-bounds 2.2,2.3,2.4 --max-cases 4 \
  --require-complete --output ffi-test6-source-state.json
```

See [the FFI source-state envelope note](docs/ffi-source-state-envelope.md).

When a reviewer has an owner-exported Test 6 wind history, the separate
observation replay can be run with
`tools/audit_ffi_test6_transient_receptors.py`. It requires the exact
`time_s,wind_speed_ms,wind_direction_from_deg` CSV contract, a declared
steady-wind table and window, fingerprints the history, and reports compact
true/indicated sensor statistics. It does not add plume storage or infer a
time-varying source; its execution artifact always keeps promotion disabled.
Use the same tool with `--verify <replay.json>` to recheck the history and
reference-table hashes without rerunning the plume.

`degali field-source` validates a separate, already-atmospheric source-rate
CSV. Its strict evidence record fingerprints the file and binds it to a source
boundary and common clock; the final endpoint must be zero. It does not infer a
flash, liquid fraction, pool evaporation, or leak rate from pressure/level data.
See [`field-source-input.md`](docs/field-source-input.md).
The same parsed schedule can be passed to
`run_field_operational_uncertainty_envelope(...,
atmospheric_source_schedule=...)` so each deterministic source corner receives
its own refinement and fail-safe screening decision.
The resulting report also contains `deterministic_sensor_envelope`: per-detector
peak/final/time-average extrema over completed corners, with missing or withheld
traces counted explicitly rather than replaced by zero. These extrema are
deterministic sensitivity bounds, not confidence intervals.
The same field is emitted by the measured-history, source-plus-sensor and
phase/pool operational envelope reports, as well as the standalone
`--sensor-array-envelope` calibration path.
The CLI keeps that source contract in a separate file and composes it only at
the complete-envelope boundary:

```bash
degali field-screen field-case.json --uncertainty-envelope \
  --atmospheric-source-case source-case.json \
  --output field-source-execution.json
```

The source-case path and SHA-256 are retained in the execution input record.

For named batch exports, pass the parsed schedule through
`export_field_screening_batch(...,
atmospheric_source_schedules={"case-name": schedule})`. Those cases are
written as complete operational-uncertainty-envelope reports and keep the
refinement, in-plane receptor and resolved-uncertainty gates enabled; the
batch API rejects attempts to waive them. Conditional-review opt-in requires
one authorization for every named case and is retained in the manifest; the
aggregate can become `conditional_allowed` only after every deterministic
source/weather/sensor corner passes those gates.

For a file-based named batch, use the strict
[`field-batch-input.md`](docs/field-batch-input.md) contract:

```bash
degali field-batch batch.json --output-directory field-batch-output
```

The execution record fingerprints the batch file, every case JSON and every
optional source-schedule JSON. Its top-level
`degali.field-batch-summary.v1` `summary` is derived from the typed
operational decisions and reports case/status counts plus the
fail-safe `all_screening_allowed` flag. Historian and phase-routing cases are kept on
their dedicated envelope paths rather than mixed into this batch.
For a measured-history `field-screen` execution, `field-verify` also reimports
the current historian CSV and compares its SHA-256/event provenance with every
nested report, so changing the CSV without regenerating the execution fails
closed.
Complete bounded `nominal_field`, phase-routing, standalone sensor-array and
joint source-sensor executions are replayed from their serialized numerical
options; the joint source-sensor path also reopens its optional atmospheric
schedule. Only genuinely option-incomplete envelope or batch records remain
integrity/semantic checks without a fabricated replay.
`degali field-verify` also checks serialized semi-FV conservation diagnostics
inside non-replayable envelope and batch reports: per-source schedule residuals
must be finite and consistent with their recorded maximum, and any residual
outside the numerical gate must carry a matching blocked applicability reason
and operational `gate_codes` entry. Compact transport concentration and sensor
extrema/time/density records are also checked, including the invariant that a
withheld detector has warnings and no fabricated result. The serialized
operational decision itself must retain consistent status, allowance,
reason/action, and unique gate-code fields. Source-level injected-mass labels
and totals are cross-checked against the ledger and schedule residual records.
It cannot be combined with the measured-history or phase-routing source paths,
and it is rejected unless every source corner is refined and gated.
The CLI direct-envelope composition accepts only
`post_flash_atmospheric_vapour`; pool or droplet schedules require their
location-bearing distributed/phase-routing handoff and are not placed at the
nozzle by inference.

For an explicitly routed liquid/rainout branch, the [phase-routing transport
handoff](docs/field-phase-routing-transport.md) preserves direct flash and
adds only a conservative time-resolved pool-evaporation ledger as a separate
field source. It withholds missing pool source histories rather than filling
them with a fitted evaporation or wake coefficient. The dedicated
[`degali.field-phase-routing-screening-input.v1` example](docs/field-phase-routing-case.example.json)
can be executed with `degali field-screen ... --uncertainty-envelope`; the CLI
requires the complete phase/pool/weather/surface/sensor envelope for this path.
If evidence-backed bounds are available for phase/pool scalars, the strict
`phase_routing.uncertainty` object propagates their lower/upper corners (with
the upper post-release bound enforced in the transport duration); no phase
range or droplet-population distribution is inferred. The pool launch scalar
width can also carry an explicit bounded `m` interval, which is evaluated as
separate source-shape corners. Droplet population alternatives are accepted
only as complete evidence-backed class corners whose mass fractions sum to one.

To audit DEGALI against supplied SLABx/CFD predictions—or to expose a
steady-versus-transient mismatch—use the separate fixed-sensor comparison
contract. It preserves the two CSV input fingerprints, units, sensor geometry,
source/weather IDs, temporal operator and (when applicable) the obstacle-mask or
wake-closure representation rather than inferring equivalence from model names:

```bash
degali field-compare comparison.json --output new-comparison.json --require-comparable
```

See the [comparison template](docs/field-model-comparison-case.example.json)
and [comparison protocol](docs/field-model-comparison.md). A mismatch retains
numerical rows but withholds a model-selection conclusion; even an aligned
comparison reports only monitored receptors, not an interpolated risk distance.

## Installation

Install Python 3.10 or newer. For LH2 work, install DEGALI with the CoolProp
extra:

```bash
python -m pip install "degali[coolprop]"
```

Check the installed version and command-line interface:

```bash
python -c "import degali; print(degali.__version__)"
degali --help
```

CoolProp is required for LH2 calculations. The traditional DEGADIS
compatibility path can be installed with `python -m pip install degali`.

## Quick manual

### Liquid-hydrogen release without an input deck

The simplest public interface is `degali.lh2.assess`. The example below
represents a horizontal 0.1 kg/s release through a 10 mm orifice, 0.5 m above
ground, from saturated LH2 stored at 6 bar absolute:

```python
from degali.lh2 import assess

result = assess(
    rate=0.10,                  # kg/s
    wind=2.0,                   # m/s at release height
    height=0.50,                # m
    orifice=0.010,              # m; use pool_diameter instead for a pool
    storage_pressure=6.0,       # bar(a), not gauge pressure
    ambient_temperature=288.15, # K
    relative_humidity=65.0,     # percent
    ambient_pressure=101325.0,  # Pa
    max_distance=30.0,          # m
    at_distance=10.0,           # m
)

print(result.report())
for warning in result.warnings:
    print("WARNING:", warning)
```

The principal outputs are:

- `distance_to_lfl`: centreline distance to the hydrogen LFL of 4 mol%, m;
- `distance_to_stoichiometric`: centreline stoichiometric distance, m;
- `lowest_flammable_height`: lowest flammable-gas height, m;
- `regime`: `grounded`, `low`, or `aloft`;
- `trajectory`: NumPy columns `[distance, centre height, mole fraction]`;
- `warnings`: validation-range and applicability warnings. Do not discard them;
- `screening_scope`: `qualified`, `conditional`, or `out_of_scope`.

For an automated screening gate, add `strict_scope=True`; an out-of-range
request raises `ApplicabilityError` instead of returning an extrapolated value.
For explicit source/wind sensitivity, use `assess_envelope(rates=[...],
winds=[...], ...)`. The supplied values are evaluated without fitting a
correction or claiming a statistical confidence interval.

To diagnose a measured centreline point such as the FFI Test 6 residual, use
`assess_observation_envelope(...)`. It returns every declared rate/wind
hypothesis within the requested factor of the observation. Multiple matches
are reported as non-identifiable; DEGALI never chooses a fitted correction.
Sensor-height or arc-maximum values require a separate validated observation
operator and cannot be inferred from this centreline diagnostic.
For exact mast coordinates, `project_lh2_jet_to_sensors(points_m=[...])`
evaluates the Gaussian vertical/lateral profile directly and returns mole
fraction and temperature at each point.

Horizontal releases at a non-zero bearing can be explored with
`run_lh2_yawed_crosswind_research()`. This uses the conserved six-flux yaw
kernel, rejects reverse-axial flow, and is explicitly marked unvalidated; it
does not change the primary `assess()` path.

For a defensible comparison with another consequence model, use the comparison
gate rather than comparing screenshots or separately reduced tables. It
requires identical source, wind vector, receptor operator, averaging window,
phase closure and geometry; a fitted prediction is reported but cannot win a
ranking:

```python
from degali.validation import ComparisonCase, ModelPrediction, compare_models

case = ComparisonCase(
    source_mode="horizontal_jet", source_rate_kg_s=0.12,
    wind_speed_m_s=2.0, wind_direction_rad=0.0, release_height_m=1.0,
    receptor_operator="arc_max_1s_peak", averaging_time_s=1.0,
    phase_closure="measured_throat", geometry="open_horizontal",
)
report = compare_models(
    observations,
    [
        ModelPrediction("degali", degali_values, case),
        ModelPrediction("other_model", other_values, case),
    ],
    case=case,
)
print(report.ranking_allowed, report.results)
```

This protocol does not claim that DEGALI is generally more accurate than
HyRAM, PHAST or EFFECTS. It makes that claim testable when independently
generated predictions are available and otherwise reports that no ranking was
performed.

For a workflow-level safety check, `degali.validation.evaluate_screening()`
combines the result scope with an optional obstacle-geometry screen. It can
approve a qualified clear-path *screening* run, but it always returns
`design_basis_allowed=False` and `approval_allowed=False`.

Specify exactly one source geometry: `orifice` for a pressurised jet or
`pool_diameter` for a pool/evaporation source. `storage_pressure` is bar(a),
while `ambient_pressure` is Pa. All temperatures are K.

### Save the LH2 trajectory

```python
import numpy as np

np.savetxt(
    "lh2_trajectory.csv",
    result.trajectory,
    delimiter=",",
    header="distance_m,height_m,centreline_mole_fraction",
    comments="",
)
```

### Existing DEGADIS input decks

```python
from degali import run_steady, run_transient, run_jet_to_ground

profile, source = run_steady("CASE.INP", backend="legacy")
print(profile.distance_to(0.05))

transient = run_transient("TRANSIENT.INP")
profile, jet, source = run_jet_to_ground("JET.INO", "GROUND.IN")
```

The equivalent command-line workflows are:

```bash
degali steady CASE.INP
degali transient TRANSIENT.INP --snapshot 60 --snapshot 120
degali dose TRANSIENT.INP --at 50 --at 100 --at 200
degali jet JET.INO --bridge GROUND.IN
```

DEGADIS decks are positional files: a missing value changes the meaning of
every field that follows. Their pressure convention also differs from
`assess()`. Preserve the original deck, verify units, and change one field at
a time. Input decks and licensed third-party data are not bundled.

For complete installation, input, output, plotting, parameter-study, deck,
validation-range, and troubleshooting instructions, see the root-level
**[English user guide](USER_GUIDE.md)** or
**[한국어 사용자 가이드](USER_GUIDE_KO.md)**.

## Verification

Fast development checks:

```bash
python -m pytest -m "not slow" -q
```

The original Fortran oracle and raw REDIPHEM, SMEDIS and PRESLHY files are not
redistributed. Tests that need
them skip with an explicit message unless the documented environment variables
are configured. Their provenance, reduction method and aggregate validation
results are documented. See [publication scope](docs/publication-scope.md).

The consolidated [technical reference](docs/technical-reference.md) describes
the governing physics, LH2 extensions, validation process, quantitative
results and remaining model-form limitations.

For the IJHE manuscript track, use the
[submission-readiness record](docs/ijhe-submission-readiness.md) and run
`python tools/audit_ijhe_submission.py` before circulating a draft.

## Repository layout

```text
src/degali/  model and validation code
tests/       self-contained numerical, physical and regression tests
tools/       research and release verification utilities
docs/        derivations, provenance, audits, limitations and results
```

## Scope and safety

The current validated domain does not cover arbitrary equipment conditions,
obstacles, indoor releases, pool spreading, downward/strongly wind-steered
jets, or general transient source behaviour. Always compare safety-critical
results with independent experiments and an accepted consequence-analysis
workflow. See [security and safety reporting](SECURITY.md).

## Citation and license

Citation metadata, including the published concept DOI and the current
version-specific DOI, are provided in [`CITATION.cff`](CITATION.cff). The
README DOI badge uses the concept DOI so it continues to resolve to the latest
archived DEGALI release.

The exact archive for DEGALI `v0.3.0` is
[`10.5281/zenodo.23256538`](https://doi.org/10.5281/zenodo.23256538).

The complete GitHub, Zenodo DOI and PyPI release sequence is documented in the
[publication guide](docs/publication-guide.md).

The DEGALI Python implementation is available under the MIT License. Original
Fortran implementations and external experimental datasets are not bundled or
redistributed.

The chronological development log that previously occupied this front page is
preserved in [docs/README-development-log-2026-09-06.md](docs/README-development-log-2026-09-06.md)
and [CHANGELOG.md](CHANGELOG.md).
