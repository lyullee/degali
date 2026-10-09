# Observation-operator-aware validation and decision implications of a fast liquid-hydrogen dispersion model

**Article type:** Research article; **Manuscript word count (including abstract and keywords; excluding references and declarations):** 4,247; **Abstract word count:** 165

## Abstract

Liquid-hydrogen (LH2) release models must represent cryogenic source behavior
while remaining fast for engineering screening. We evaluate DEGALI with a
frozen no-fit flashing source, a receptor observation operator and free-ground
detachment. Nine PRESLHY E3.5 releases (62 arc maxima) give MG = 1.047,
VG = 1.425 and FAC2 = 0.839, where MG is geometric-mean bias, VG is geometric
variance and FAC2 is the factor-two agreement fraction. Synchronized 20 s means on seven trials give
MG = 0.544, VG = 4.601 and FAC2 = 0.761, demonstrating operator sensitivity.
A six-arc FFI/DNV screen gives MG = 1.245, VG = 1.373 and FAC2 = 0.833,
while Test 6 retains a 3.48-fold underprediction at 30 m. Two FFI cases
change 11 of 18 binary 4 vol % line classifications. A 3-D transient extension is numerically
verified for convergence and conservation, not field validation.
The result is an evidence-bounded LH2 screening framework linking source,
weather, sensors and decisions; it is not transient multiphase CFD or a
regulatory design basis.

**Keywords:** liquid hydrogen; flashing source; cryogenic dispersion; observation operator; model validation; engineering screening.

## 1. Introduction

LH2 releases span flashing, cold entrainment, possible atmospheric-condensate formation, a momentum-dominated jet and subsequent buoyant dispersion. Generic engineering models can be fast but are commonly supported by gaseous or non-cryogenic evidence. Detailed CFD can resolve site geometry and transience but is not a practical replacement for transparent rapid screening. This paper asks a narrower question: what can a fast LH2-specific path predict with current public evidence and no coefficient selected on validation data?

Recent IJHE work reinforces the relevance of digital-twin safety monitoring and
the growing LH2-release evidence base [1,2]. Physics-guided full-field
prediction has also been demonstrated for confined hydrogen experiments using
CFD-generated training data [3], while sensor-efficient source-term estimation
targets indoor industrial and community systems [4]. Those studies address learned
full-field or inverse source reconstruction; the present contribution is
complementary: it evaluates a no-fit outdoor LH2 screening path with a declared
sensor observation operator, separate evidence lanes and an explicit failure
boundary. This distinction prevents a rapid model from being presented as a
universal replacement for transient CFD or data-assimilating digital twins.

Accordingly, this study evaluates DEGALI on its own frozen basis. A cross-code
performance rank is excluded because the available calculations do not share
the same source representation, temporal operator and receptor mapping. DEGALI
instead adds a provenance and promotion layer: it records which source,
weather, geometry, sensor and temporal operator support each result, preserves
conditional mismatches, and routes evidence outside the declared envelope to
a higher-fidelity consequence study.

## 2. Frozen model path

The baseline uses a corrected expanded flashing source, reported-window mean E3.5 source rate, model evaluation at every sensor followed by the observed arc-maximum operator, free ground detachment and 0.001 m roughness. The FFI screen uses the published low-mast wind while recording a 10 m reference-height assumption. The complete baseline and the reproducible audit are specified in [`docs/lh2-paper-baseline-2026-09-20.md`](docs/lh2-paper-baseline-2026-09-20.md).

The baseline does not claim a validated N2/O2 mixed-phase closure, particle-slip calculation, obstacle wake model or transient puff solver. These mechanisms remain outside the primary path until independently validated.

The code now also exposes a separate `run_lh2_yawed_crosswind_research()`
branch based on the conserved six-flux yaw kernel. It accepts a horizontal
release bearing and a global wind bearing, rejects the reverse-axial branch,
and returns a trajectory with an explicit `validated = False` marker. This
extends the computational research envelope without silently converting an
unvalidated three-dimensional path into the primary claim.

For declared pool histories, `assess_pool_history()` can additionally apply a
user-supplied causal first-order response time at a receptor. The resulting
signal is explicitly labeled a low-order response surrogate; no transient
storage, meander or turbulence time scale is inferred from the existing
validation data.

The public API now exposes this boundary instead of leaving it to a reader's
interpretation. `assess(..., strict_scope=True)` raises
`ApplicabilityError` for an out-of-range request; the default warning mode is
retained for exploratory work. `assess_envelope(rates=[...], winds=[...])`
evaluates caller-supplied low/nominal/high input values as a Cartesian
sensitivity envelope. It is deliberately not a confidence interval and it
does not fit a source or transport multiplier.

### 2.1 Computational verification of the transient extension

The repository also contains a separate three-dimensional explicit
finite-volume transient solver. It couples a zero-terminated source-rate
schedule and an interpolated meteorological wind history on the same event
clock, advances a prognostic velocity closure with hydrogen-density buoyancy,
and reports receptor traces together with a mass ledger. This branch is a
research extension and is not used to alter the frozen E3.5 or FFI scores.

For numerical verification, a 48 x 48 x 32 calculation with a 0.0025 s time
step is treated as a declared high-resolution numerical reference. The 36 x
36 x 24 case gives an aggregate concentration RMSE of
`5.4386e-4 kg H2/m3`, a maximum absolute receptor difference of
`1.7501e-3 kg H2/m3`, and an apparent aggregate convergence order of 1.47
between the medium and fine levels. The largest absolute mass residual is
`2.22e-15 kg` across the reported coarse-to-reference runs. The fine-grid
residual itself is `2.22e-16 kg`. These are discretisation and conservation
results against a
numerical reference; they are not experimental accuracy, sensor accuracy or
evidence of obstacle-wake validity. The full table is retained in the
transient error-study record and the scope boundary is documented in the
supplementary methods.

## 3. Validation material and operators

### 3.1 PRESLHY E3.5

PRESLHY E3.5 provides concentration and temperature observations from elevated
LH2 releases [5]. The primary concentration score retains
horizontal, momentum-dominated releases: nine trials and 62 distance-wise arc
maxima. The public record contains 24 trial workbooks; the deterministic audit
scope retains Trials 10, 11, 12, 20, 21, 22, 23, 24 and 25 because they are
horizontal, have a usable source-flow record and pass the momentum-dominance
exit-ratio gate. Trials 6–9, 13–15 are non-horizontal or obstructed-source
cases, while Trials 2–5 and 16–19 are wind-steered at the source; they are not
silently pooled into the free-ground score. The receptor map has 197 sensors
but correlated sensors are not treated as independent trials. The raw
clock-aligned sensitivity lane uses the seven transferred workbooks (Trials
10–12 and 22–25), which is why its sample count is smaller than the primary
reduced-table score.

For three accepted events (Trials 10, 20 and 21), the public Flowmeter,
LocalWeather and Xensor streams were additionally reconciled on their recorded
wall clock, with the D3.6 sensor-coordinate table and the declared 18 s
sampling-line delay. The same D3.6 report gives nominal Xensor XEN-5320
accuracy of ±3% of full scale, an estimated cryogenic Coriolis flow accuracy
of ±3% of reading, source-relative receptor heights and two obstruction-test
offsets (0.16 m and 0.18 m); these are uncertainty/geometry metadata rather
than channel-level calibration or a surveyed obstacle mesh. Trial 10 is the
representative upload-facing bundle and
Trials 20 and 21 are independently replayed conditional bundles; none is
pooled into the headline score. Each resulting `FieldEvidenceManifest` is a
hash-pinned conditional intake record, not an accepted calibrated
field-validation input or a new headline score: the weather record is
five-minute, channel-level calibration certificates are not included, and the
open-pad geometry is not an obstacle survey. The extraction and gate record are
documented in
[`docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md`](docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md).
The public RADAR workbooks for Trials 20 and 21 were additionally retrieved
through their individual-file endpoints and matched byte-for-byte to the local
replay inputs; the retrieval audit is retained as a traceability artifact, not
as evidence that closes the calibration or geometry gate.

### 3.2 Thermal diagnostic

Thermal residuals are reported separately. The sealed diagnostic has two trials, two stations and five heights. It does not re-integrate the primary trajectory and is not used to choose a heat- or mass-transfer coefficient. Its purpose is to show the unresolved thermal/enthalpy structure rather than to conceal it within concentration statistics.

### 3.3 FFI/DNV far-field screen

FFI/DNV Tests 4 and 6 are horizontal 25.4 mm outdoor LH2 releases with
measurements at 30, 50 and 100 m [6]. They are physically
distinct from the near-field array and test extrapolation. Outdoor downward
releases are excluded because impingement creates a different ground-source
problem. Closed-room tests are excluded because the ventilation-mast outlet,
not the supply nozzle, is their atmospheric source.

The FFI report describes a test pad with two stacked containers, a plastic drum
and an instrument box near the release point, plus a wind mast and oxygen-based
concentration instrumentation. Those structures are not represented in the
present free-ground calculation. The six-arc lane is therefore a conditional
far-field/transport screen, not obstacle validation; the site-geometry evidence
and its SHA-256 provenance are recorded in
[`docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md`](docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md).

The FFI source-state audit evaluates each reported sensor coordinate over
declared source, wind and ambient corners. Its JSON keeps separate arc-maximum
and sensor-height tables, including available/withheld counts and one-sided
lower-bound diagnostics. These tables expose the observation-operator and
source-state interaction; they do not add weights, fit a coefficient or qualify
the model for operational validation.

### 3.4 ELVHYS WP4.2 confined dynamic and geometry boundary

The public ELVHYS WP4.2 HSE archive provides a useful conditional dynamic
evidence lane [7]. It contains 48 tests in a nominal 1 m3 transfer connection
space, with time-resolved hydrogen concentration, temperature, pressure,
ventilation and weather channels. The archive also reports concentration and
temperature sensor specifications, five-point calibration procedures and a
three-dimensional nozzle/sensor coordinate table. Most channels are sampled at
20 Hz; the archive includes higher-rate pressure-peaking and ignition records.

This evidence is not merged into the outdoor LH2 score. `FLMT` is fan-flow
ventilation rather than measured hydrogen mass flow, and the enclosure walls,
floor, vents and active/passive ventilation are part of the observed system.
The pre-registered Tests 10 and 11 screen therefore stops at a repeatability
and falsification boundary: their nominally matched conditions have different
nozzle pressures and temperatures, and Bottom-3 medians differ by 7.80 vol%.
The horizontal-nozzle height also differs between the archive geometry table
and the final-report table. These facts make the lane useful for demonstrating
time alignment, calibration metadata and measured geometry, but not for
fitting an outdoor source or claiming a validated obstacle wake. The complete
hash-pinned boundary and promotion rule are in
[`docs/ijhe-elvhys-conditional-evidence-2026-10-09.md`](docs/ijhe-elvhys-conditional-evidence-2026-10-09.md).

### 3.5 Hecht–Panda Raman benchmark

The public Hecht–Panda cryogenic free-jet benchmark supplies four aggregate
slopes (centerline mass and temperature decay, and mass and temperature
half-width) [8]. Its journal fit has an unresolved
condition-count/membership discrepancy, so it is reported as an external
boundary and not pooled with E3.5 or FFI scores.

### 3.6 Non-LH2 transfer boundary

The repository also retains a separate REDIPHEM/SMEDIS generic dense-gas
validation lane, including ammonia, propane, R-12 and passive-gas trials and
an obstacle control/fence pair. It is reported only as a transfer boundary:
the source thermodynamics, release regimes and observation operators differ
from LH2, and the available reductions do not provide a channel-calibrated,
common-clock LH2 time series. These data are therefore not pooled with the
headline scores and do not support a universal obstacle-wake closure. The
inventory and promotion decision are documented in
[`docs/ijhe-smedis-transfer-boundary-2026-10-09.md`](docs/ijhe-smedis-transfer-boundary-2026-10-09.md).

## 4. Metrics and reproducibility

MG, VG and FAC2 use conventional observed/predicted dense-gas model definitions
[9]; MG > 1 indicates underprediction.
`tools/audit_lh2_paper_baseline.py` produces a SHA-256 manifest, receptor-level
residual tables and local SVG maps. The output is ignored by Git; no
third-party workbook, raw time series or original Fortran code is redistributed.
OpenAI Codex was used as an AI-assisted software-development aid for selected
implementation, debugging and plotting tasks; all model equations, input
choices, numerical runs and interpretations were reviewed by the authors.

For each retained pair, let `r_i = O_i/P_i`, where `O_i` and `P_i` are the
observed and predicted concentrations. Non-finite and non-positive pairs are
removed before the logarithmic statistics. We report
`MG = exp(mean(ln r_i))`, `VG = exp(mean((ln r_i)^2))`, and
`FAC2 = mean(0.5 <= P_i/O_i <= 2)`. Thus MG is explicitly
observed/predicted, while FAC2 counts predictions within a factor of two of the
observation. The commonly used Hanna screening limits are descriptive
conventions, not certification criteria; trial clustering and the declared
observation operator remain the primary uncertainty controls.

**Table 1. Headline validation metrics under the declared observation
operators.** MG is observed/predicted; values above one indicate
underprediction. Arc and sensor rows are not independent trials.

| Evidence lane | Records | Operator | MG | VG | FAC2 | Role |
|---|---:|---|---:|---:|---:|---|
| PRESLHY E3.5 | 62 arcs from 9 trials | peak/arc maximum | 1.047 | 1.425 | 0.839 | primary LH2 comparison |
| PRESLHY E3.5 subset | 46 arcs from 7 trials | synchronised 20 s mean | 0.544 | 4.601 | 0.761 | temporal-operator sensitivity |
| FFI/DNV | 6 arcs from Tests 4 and 6 | peak/arc maximum | 1.245 | 1.373 | 0.833 | independent far-field screen |

The trial-cluster bootstrap and all receptor-level residuals remain in the
supplementary tables and are regenerated by the audit scripts. This compact
table is the upload-facing summary; it does not pool the evidence lanes.

As a descriptive uncertainty check, resampling complete trials (5,000 fixed-
seed draws) gives corrected E3.5 95% percentile intervals of 0.771–1.416 for
MG, 1.132–1.882 for VG and 0.677–0.983 for FAC2. The corresponding corrected
FFI intervals are 1.031–1.503, 1.075–1.754 and 0.667–1.000, but only two tests
are available. These intervals exclude source, meteorological and model-form
uncertainty; they are trial-cluster sensitivity diagnostics rather than
population-level validation claims.

The paired comparison also runs the local pre-correction
`corrections=False` reconstruction and the frozen corrected path on identical
source records and receptor positions (`tools/audit_applied_energy_comparison.py`).
This comparator is not the original DEGADIS executable or independent
software. Both model paths use a 60 s averaging setting, but field values
are release-window peaks rather than synchronized 60 s means: the
model–observation temporal operator is not fully matched. Trial-level 95%
percentile intervals resample arcs descriptively; aggregate intervals
resample complete trials (5,000 fixed-seed draws). They exclude physical
source, weather and model-form uncertainty. Timings include source,
integration and receptor extraction, after a warm-up, with three sequential
runs per trial.

## 5. Results

### 5.1 Near-field concentration

The frozen E3.5 score is MG = 1.047, VG = 1.425 and FAC2 = 0.839 (n = 62). These meet commonly reported screening thresholds, but the clustered trial population is insufficient for a universal accuracy claim. The residual map is retained because it exposes distance, height and crosswind structure that aggregate scores hide.

### 5.2 Thermal structure

The thermal map exhibits distance- and height-dependent residuals that do not collapse to a justified scalar correction. The present evidence does not support changing baseline heat or density transport. It provides a falsifiable target for a future thermal closure.

### 5.3 Far-field screen

For six FFI/DNV arc maxima, the frozen path yields MG = 1.245, VG = 1.373 and FAC2 = 0.833. Test 4 is balanced at 30 and 50 m and model-high at 100 m. Test 6 is strongly model-low at 30 m, close at 50 m and model-high at 100 m. Aggregate performance is useful but does not cancel the physical importance of the low-wind residual.

### 5.4 Test 6 decomposition

At 30 m, low/mean/high mast winds predict 6.03/6.72/7.36 vol % against a 21.0 vol % reported arc maximum, leaving a factor 2.85–3.48 mismatch. Exact fixed-wind sensor geometry changes the predicted maximum only to 5.75 vol %. The public table has no synchronized independent atmospheric two-phase source state, so source uncertainty cannot be decomposed honestly by invented ranges. Forced ground contact produces 25.50 vol % at this point but worsens the independent six-arc screen. Shape drag, vertical shear and a literature-range lift-off diagnostic do not close the residual. The full decision table is in [`docs/ffi-test6-decomposition.md`](docs/ffi-test6-decomposition.md).

### 5.5 External Raman boundary

Against the published Hecht–Panda aggregate slopes, the current branch is high by factors 2.09 (centerline mass), 2.72 (mass half-width), 2.41 (centerline temperature) and 3.46 (temperature half-width). The separately recorded conserved-energy research branch is closer, but its source-establishment and condition-membership limitations prevent promotion into the frozen crosswind path. This closes a tempting but unsupported universal free-jet claim.

### 5.6 Paired predecessor comparison, trial uncertainty and runtime

On the same 62 E3.5 arcs, the historical reconstruction yields MG 1.156,
VG 23.101 and FAC2 0.774, versus corrected 1.047, 1.425 and 0.839.
The corrected MG trial-cluster bootstrap interval is 0.771–1.416. The result
is heterogeneous: corrected Trial 20 still has MG 2.366 and Trial 21 has
MG 0.473. On six FFI arcs, the historical reconstruction gives MG 12.600 and
FAC2 0; corrected gives MG 1.245 and FAC2 0.833. Only two FFI trials are
available, precluding a stable population interval. Median warm-run time per
trial in the recorded run is 0.351 s for corrected E3.5 and 1.567 s for
corrected FFI; corresponding historical-path times are 0.373 and 1.666 s.
A prior same-workspace run was about twice as fast under different background
load, so no external-model/CFD speed advantage or stable speed ratio is
inferred. Full trial, arc, interval
and runtime tables, figures and hashes are specified in
[`docs/applied-energy-evidence-2026-10-03.md`](docs/applied-energy-evidence-2026-10-03.md).

### 5.7 Common-window time alignment and input envelope

A separate raw-workbook sensitivity replaces release-window peaks with a
pre-specified common temporal operator. For Trials 10, 11, 12, 22, 23, 24
and 25, the central 20 s of the sustained-flow record supplies the source
mean, and the physical Xensor channels are averaged over the matching clock
window after the documented 18 s sampling-line delay. The model uses the
same 20 s averaging setting and physical receptor locations. On 46 arcs this
operator gives MG 0.544, VG 4.601 and FAC2 0.761. Thus the model is high
against a common-window mean even though its peak-based aggregate is nearly
unbiased. The two scores answer different questions and are not pooled.

An evidence-bounded source/wind screen combines the 20 s source mean with
the pre-registered peak-rate bound, and reported mean with maximum wind.
Including the zero- versus 18 s delay operator sensitivity gives aggregate
MG 0.510–0.547 and FAC2 0.745–0.766. It does not reverse the model-high mean
result. At a fixed 1.5 m receptor, trial-specific 4 vol % reach ranges span
11.5–13.2 m (Trial 24) to 39.1 m/>40 m (right-censored, Trial 11). These are
deterministic evidence envelopes, not confidence intervals. Sensor
calibration, source pressure/temperature, wind direction/meander and model-
form distributions remain unavailable. Full methods and tables are in
[`docs/time-aligned-uncertainty-results-2026-10-03.md`](docs/time-aligned-uncertainty-results-2026-10-03.md).

### 5.8 Decision-oriented reference cases and evidence boundary

For this count, a line classification is true when the steady-model
concentration at a fixed 1.5 m receptor is at least 4 vol % at the declared
screening line, and false otherwise. The 18 classifications are the Cartesian
product of two FFI cases, three reported wind values and the 20, 30 and 50 m
lines; they are decision sensitivities, not safety-distance determinations.

Two transfer-like FFI cases were evaluated at a fixed 1.5 m receptor height,
three reported-wind values and declared 20, 30 and 50 m screening lines. For
Test 4, the historical reconstruction gives a 22.7–31.6 m LFL reach and the
corrected path gives 56.5–67.9 m. For low-wind Test 6, the corresponding
ranges are 12.6–14.5 and 34.3–37.8 m. Eleven of the 18 line classifications
change between model paths. This is a decision-sensitivity result, not a
setback: Test 6 retains its unresolved 30 m residual and both paths are
conditional outside their directly supported operators.

Four NASA pool-spill states provide a separate storage/spill archetype. The
model reproduces the reported grounded/low/aloft regime in all four and gives
a 2.47 m mean absolute error in lowest flammable height. Pool-axis LFL reaches
are 29.1–66.5 m, but three exceed the 33.8 m concentration-distance evidence
boundary and remain extrapolations. Pool-axis and fixed-height jet distances
are not pooled.

Published FLACS and ADREA-HF “HSL Test 6” calculations belong to a different
2010 campaign at 0.071 kg/s, not the 2019 FFI Test 6 at 0.833 kg/s, and cannot
close the present residual. No same-input PHAST, EFFECTS or CFD receptor table
is available, so this paper does not report a cross-code performance rank.
Full inputs and tables are in
[`docs/decision-oriented-scenarios-2026-10-03.md`](docs/decision-oriented-scenarios-2026-10-03.md).
The distinction is consistent with prior LH2 CFD studies that treat release
and cryogenic-dispersion conditions as scenario-specific [10,11], while
pool-spill modeling is a separate archetype [12].

### 5.9 Decision-boundary table and figure captions

**Table 2. FFI Test 6 low-wind boundary at the 30 m receptor.** The three
predictions are the declared low, mean and high mast-wind cases; the reported
arc maximum is 21.0 vol %. The source record does not support a fitted
correction.

| Case | Predicted concentration (vol %) | Ratio to reported maximum |
|---|---:|---:|
| Low mast wind | 6.03 | 3.48-fold underprediction |
| Mean mast wind | 6.72 | 3.13-fold underprediction |
| High mast wind | 7.36 | 2.85-fold underprediction |
| Reported arc maximum | 21.0 | reference |

The following captions state the evidence status so that numerical verification
and conditional boundaries cannot be mistaken for field validation.

**Figure 1. Evidence-gated LH2 screening workflow.** Source, weather,
obstacle, receptor and sensor identities are checked before a result can be
labeled qualified.

**Figure 2. PRESLHY E3.5 observed-versus-predicted residual map.** The map
shows distance, height and crosswind structure for the 62 primary arc maxima;
it is a measured-comparison residual, not an uncertainty probability map.

**Figure 3. Peak and synchronised 20 s observation operators.** The two lanes
use the same declared receptors but answer different temporal questions.

**Figure 4. FFI/DNV horizontal-release arc residual map.** The map shows the
six arc maxima from Tests 4 and 6; the low-wind Test 6 near-array mismatch is
retained as an application boundary.

**Figure 5. Decision-line sensitivity.** Historical and corrected paths are
shown as conditional model-choice sensitivities at 20, 30 and 50 m, not as
regulatory setbacks.

**Figure 6. Numerical convergence and conservation of the 3-D extension.**
The comparison is against a declared numerical reference and is not measured
transient validation.

## 6. Discussion and application envelope

DEGALI is appropriate for qualified research and engineering screening of represented momentum-dominated LH2 releases when users retain the documented source assumptions and meteorological sensitivity. It is not, on the evidence here, a replacement for transient multiphase CFD when low-wind temporal evolution, complex obstacles, walls or site-specific source geometry determine the consequence. It is therefore not a substitute for transient multiphase CFD in a safety-critical design basis. The Test 6 mismatch is reported as an application limit rather than hidden by a fitted term.

For the Test 6 wind interval, the modelled 4 vol % concentration exit at the
1.5 m-high downwind coordinate is 12.6–14.5 m for the historical
reconstruction and 34.3–37.8 m for the corrected path. Both 20 and 30 m lines
therefore receive opposite classifications at all three winds. The Test 4
case also changes classifications at the 30 or 50 m lines. These calculations
illustrate model-choice consequences, not operational setbacks: the Test 6
30 m measured arc maximum is underpredicted 3.48-fold by the corrected path,
and the peak-versus-steady time operator is unresolved.

For operational use, a result with `screening_scope == "qualified"` is inside
the frozen evidence envelope. A `conditional` result is an extrapolation that
must be carried as a sensitivity, and `out_of_scope` includes wind-steered
releases for which a steady jet concentration is not defensible. The strict
guard is an implementation safeguard, not a claim that the underlying model
has become universally accurate.

In an industrial workflow, this boundary makes DEGALI most useful as a fast
front-end for hazard identification, pre-FEED scenario ranking, sensor-array
design and emergency-response triage. The output is a provenance-linked set
of concentration and decision sensitivities that can identify cases requiring
transient multiphase CFD or a site-specific consequence study. It is not a single
automatic separation distance: each result carries its source,
weather, receptor, sensor and observation-operator identities, and a result
that fails those evidence gates remains conditional or withheld. This
positioning complements higher-fidelity transient and obstacle-resolving
studies; it does not claim that DEGALI replaces them.

## 7. Conclusions

1. A no-fit fast LH2 baseline was frozen and reproduced on 62 PRESLHY E3.5 arc maxima with MG 1.047, VG 1.425 and FAC2 0.839.
2. The independent six-arc FFI/DNV screen gives MG 1.245, VG 1.373 and FAC2 0.833, while exposing a material low-wind near-array limitation.
3. Wind choice and fixed sensor geometry do not resolve Test 6; the public source record cannot separately identify atmospheric flashing-source uncertainty.
4. Deterministic transport alternatives do not justify a new default closure. Free detachment remains the baseline and the mismatch is an explicit limit.
5. The Hecht–Panda slope boundary rules out an unqualified transfer to every cryogenic free-jet regime.
6. Peak and synchronized 20 s mean operators lead to materially different bias and scatter; both must accompany any validation claim.
7. The historical and corrected paths change 11 of 18 declared-line classifications in two transfer-like cases, showing a decision consequence without claiming a safety setback.
8. The supported claim is a reproducible LH2 fast-screening model with a validation envelope, not universal industrial certification or transient/obstacle-resolved prediction.

## Declarations

**Data availability.** The primary PRESLHY E3.5 source record is publicly
available through KITopen/RADAR (DOI [10.5445/IR/1000136281](https://doi.org/10.5445/IR/1000136281);
RADAR version DOI [10.35097/1481](https://doi.org/10.35097/1481)). The
byte-level retrieval audit for Trials 10, 20 and 21 is retained in the
submission evidence package as `preslhy-public-retrieval-index-2026-10-09.json`.
The conditional
ELVHYS WP4.2 archive is publicly available at DOI [10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H);
its raw files are not copied into this repository. Local audit products are
derived residuals; original third-party workbooks, raw time series and report
files are not redistributed. The reported summaries can be regenerated when
the corresponding evidence is available under its licence and the SHA-256
inputs are recorded.

**Code availability.** The public implementation is available at
`https://github.com/lyullee/degali`. The submission software release is
DEGALI v0.3.0 (version-specific DOI
[`10.5281/zenodo.23256538`](https://doi.org/10.5281/zenodo.23256538)); the
version-neutral concept DOI for the complete release series is
[`10.5281/zenodo.22646258`](https://doi.org/10.5281/zenodo.22646258). The
IJHE source package pins its local source and derived files in the submission
manifest; any later submission release must receive a new version-specific DOI
and archive its exact tag before that identifier replaces the one above. The audit workflows are
`tools/audit_lh2_paper_baseline.py`,
`tools/audit_applied_energy_comparison.py`,
`tools/audit_time_aligned_uncertainty.py`,
`tools/audit_decision_scenarios.py`. A publication snapshot must exclude
original DEGADIS Fortran and external raw experimental data.
The pinned `requirements-research.txt` file records the Python 3.12 runtime
used for the reported audits and document builds; it is included in the
submission manifest so that the numerical and artifact checks can be replayed
without inferring package versions from the source code alone.

**Funding.** `[AUTHOR TO CONFIRM: funder, grant number and the funder's role in
study design, data collection, analysis, writing and the decision to submit;
state "No external funding" if applicable.]`

**CRediT authorship contribution statement.** `[AUTHOR TO CONFIRM: assign
each named author only the contributions they actually made, using the CRediT
taxonomy (for example, conceptualization, methodology, software, validation,
formal analysis, investigation, data curation, visualization, writing,
supervision, project administration and funding acquisition).]`

**Declaration of competing interest.** `[AUTHOR TO CONFIRM: provide the
journal's exact competing-interest statement, including "The authors declare
that they have no known competing financial interests or personal
relationships" if that is accurate.]`

**Acknowledgements.** `[AUTHOR TO CONFIRM: acknowledgements or "None".]`

### Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work, the authors used OpenAI Codex as an AI-assisted aid for software implementation, code review, debugging, reproducible plotting scripts and language editing. The tool did not generate experimental observations or replace the declared validation calculations.
The authors reviewed and verified the resulting code, evidence provenance,
scientific interpretation and manuscript, and retain full responsibility for
the final work.

## References

[1] Naanani, H., Nachtane, M., and Faik, A. (2025). Advancing hydrogen safety and reliability through digital twins: Applications, models, and future prospects. *International Journal of Hydrogen Energy*, 115, 344–360. DOI: 10.1016/j.ijhydene.2025.02.440.

[2] Liang, Y., He, M., Qu, Y., et al. (2025). Progress and trends in liquid hydrogen release. *International Journal of Hydrogen Energy*, 144, 1147–1167. DOI: 10.1016/j.ijhydene.2025.04.026.

[3] Wang, W., Gao, W., and Li, Y. (2026). PhyGNN-LSTM: A physics-constrained spatiotemporal framework for real-time hydrogen dispersion prediction. *International Journal of Hydrogen Energy*, 235, 155179. DOI: 10.1016/j.ijhydene.2026.155179.

[4] Li, A., Lang, Z. Q., Ni, C., et al. (2026). Sensor-efficient, data-driven estimation of hydrogen leak source terms for indoor industrial and community systems. *International Journal of Hydrogen Energy*, 252, 156184. DOI: 10.1016/j.ijhydene.2026.156184.

[5] Lyons, K., Coldrick, S., and Atkinson, G. (2021). *Summary of experiment series E3.5*. PRESLHY. DOI: 10.5445/IR/1000136281.

[6] Aaneby, J., Gjesdal, T., and Voie, Ø. (2021). *Large scale leakage of liquid hydrogen (LH2) — tests related to bunkering and maritime use of liquid hydrogen*. FFI Report 20/03101. Available at: https://www.ffi.no/en/publications-archive/large-scale-leakage-of-liquid-hydrogen-lh2-tests-related-to-bunkering-and-maritime-use-of-liquid-hydrogen.

[7] Health and Safety Executive, UK (n.d.). *ELVHYS WP4.2 — Leakage into cold room/tank connection space considering barriers and obstacles*. Public dataset, DOI: 10.18710/JXJP0H.

[8] Hecht, A., and Panda, J. (2019). *Study of key parameters in modeling liquid hydrogen release and dispersion in open environment*. International Journal of Hydrogen Energy, 44. DOI: 10.1016/j.ijhydene.2018.07.058.

[9] Hanna, S. R., Chang, J. C., and Strimaitis, D. G. (1993). Hazardous gas model evaluation with field observations. *Atmospheric Environment*, 27A, 2265–2285. DOI: 10.1016/0960-1686(93)90397-H.

[10] Ichard, M., Hansen, O. R., Middha, P., and Willoughby, D. B. (2012). CFD computations of liquid hydrogen releases. *International Journal of Hydrogen Energy*, 37, 17380–17389. DOI: 10.1016/j.ijhydene.2012.05.145.

[11] Giannissi, S. G., Venetsanos, A. G., Markatos, N., and Bartzis, J. G. (2014). CFD modeling of hydrogen dispersion under cryogenic release conditions. *International Journal of Hydrogen Energy*, 39, 15851–15863. DOI: 10.1016/j.ijhydene.2014.07.042.

[12] Holborn, P., Ingram, J. M., and Benson, C. M. (2020). Modelling hazardous distances for large-scale liquid hydrogen pool releases. *International Journal of Hydrogen Energy*, 45, 23851–23871. DOI: 10.1016/j.ijhydene.2020.06.131.
