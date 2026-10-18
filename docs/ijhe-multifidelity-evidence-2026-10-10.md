# IJHE multi-fidelity evidence intake — 2026-10-10

This intake adds two public experimental component datasets and one public
multi-code CFD report without changing the frozen LH2 headline score.  The
evidence hierarchy is fixed as experiment first, independently exported model
predictions second, and published-figure context third.

## 1. USN ventilated car-park dynamics

The public USN record (DOI `10.23642/USN.14405903`) provides MATLAB v7.3 files
for four non-ignited GH2 blowdowns used here (Experiments 19, 20, 21 and 23).
The intake verifies the publisher MD5 values, reads the mass-flow,
ventilation and concentration clocks, and applies the paper's room-cloud
operator to the 24 CAN-bus sensors D139–D162.  The separate under-table sensors
D163–D168 are not silently pooled with that operator.

The reproduced first 4 vol% arrivals are 15.923, 17.832, 32.402 and 30.498 s.
The corresponding exceedance durations are 285.10, 335.70, 195.31 and 88.51 s.
Relative to the paper table, the maximum absolute arrival difference is
0.498 s and the maximum duration difference is 5.51 s.

A no-fit ideal well-mixed mass balance is retained only as a deliberately low-
fidelity comparator.  Across the four events its model/observed peak ratio is
1.542–1.869, with a mean of 1.619, and it substantially overpredicts the
flammable duration.

A separate bounded dynamic enclosure response was therefore added rather than
forcing the outdoor 3-D closure into the container. It advances H2 volume
fraction with source displacement, measured ventilation and a concentration-
dependent buoyant-opening exchange. Exp19 alone identifies effective mixed-
volume fraction (0.547503), effective opening area (0.586732 m2) and transport
delay (5.1 s). With these coefficients frozen, pooled Exp20/21/23 spatial-mean
MAE is 0.282 vol%p, RMSE is 0.332 vol%p and Pearson `r = 0.9795`, a 95.62% MAE
reduction from the same-window ideal well-mixed comparator. The maximum mass
residual is `6.44e-15 kg` and volume fraction is bounded to 0–100%.
For the three holdouts, the spatial-mean 4 vol% first-arrival MAE is 3.70 s
(maximum 5.70 s), while flammable-duration MAE is 27.1 s (maximum 45.6 s).
These are spatial-mean threshold metrics rather than any-sensor alarm scores.
The bound-span-scaled fixed-delay Jacobian condition number is 4.87. A
descriptive 1% SSE profile retains two adjacent delay values (4.8 and 5.1 s),
mixed-volume fraction 0.5475--0.5529 and effective opening area
0.5852--0.5867 m2. These are profile diagnostics, not confidence intervals.

The score is for the enclosure-average mean of non-flatlined room sensors,
not individual probes or the any-sensor alarm envelope. The model form was
selected with knowledge of the complete public record, so Exp20/21/23 are
retrospective parameter holdouts rather than prospective blind validation.
The lane is not a DEGALI outdoor free-field score and not LH2 evidence.

Audit command:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe tools\audit_usn_carpark.py `
  reference\datasets\usn-carpark-14405903\raw `
  --output outputs\usn-carpark-dynamic-audit-2026-10-10.json

.venv\Scripts\python.exe tools\validate_usn_enclosure_dynamic.py
```

## 2. HyTunnel-CS Pro-Science spatial matrices

The public Zenodo record (DOI `10.5281/zenodo.4668554`) contains four
concentration workbooks, one ventilation-velocity workbook, a metadata PDF and
an Autodesk Inventor CAD archive.  All seven published MD5 values are verified.
The adapter preserves average, maximum and minimum matrices without inventing
a transient clock or measurement height.

The workbooks contain 42 configuration sheets: 41 are labelled H2 and one is
labelled He. The apparent conflict with the published 40-case matrix has now
been reconciled condition by condition. Thirty-nine declared H2 cases are
present. The declared 1 mm, 1.5 g/s, 3.5 m/s cross-flow H2 case is represented
by a helium sheet, while two extra stagnant H2 baselines at 1.5 and 2.5 g/s are
included. These substitutions remain explicit evidence gates rather than being
silently counted as equivalent H2 validation. The parser reads 15,768
concentration values and 1,260 signed x-velocity values.
The published concentration range is -0.24–62.17 vol%; the small negative
values are retained as measurement/background-correction results rather than
clipped.  The velocity range is -8.722–6.865 m/s.  The CAD archive contains 20
Inventor part files and 10 assemblies, with no unsafe archive paths, but it has
not yet been converted into an independently verified computational geometry.

The ICHS experiment paper additionally fixes the observation protocol at nine
simultaneous concentration channels, 2 Hz sampling, and a 10 s steady averaging
window for 5 g/s releases. The official D2.3 report defines one matched
1 mm/5 g/s/3.5 m/s co-flow CFD case on a 2.5 x 0.9 x 0.9 m domain, with the
source axis at 0.45/0.45 m, 6.4 mm cells and a 7.25 mm notional nozzle. Its
published experiment/CFD 4 vol% lengths are 1.65/2.25 m. DEGALI scoring remains
withheld. The high-pressure source adapter exists, but the official CFD's
open non-reflecting condition on the five non-inlet faces does not. The
energy-conserving DEGALI
notional nozzle is `4.364 mm` rather than the official CFD's reported
`7.25 mm`, and the public report does not provide enough state variables to
reconstruct that diameter as a complete conserved plane. No concentration-
based source fitting is allowed.

The adapter entrains ambient co-flow until both Mach `<=0.30` and the receiving
active-scalar density departure `<=10%` pass. For the matched case this gives a
`0.62069 m` handoff distance, `0.18917 m` equivalent diameter, Mach `0.05588`
and `29.56` cells across the official mesh. Conservation closes to round-off.
The current boundary candidate uses a non-equilibrium prescribed-velocity
inlet and fixed ambient pressure with adjacent-plane velocity/stress on the
other five faces. Its variable-density projection, momentum boundary and
scalar/enthalpy outflow are coupled and unit verified. Fixed pressure reflects
91.72% of a linear acoustic pulse and is rejected; characteristic
incoming-wave suppression reduces the ratio to 1.64% on both streamwise and
lateral-normal faces and passes both declared 5% gates. One declared 30°
tangential mode also passes at 4.188%, but the original one-dimensional-
impedance sweep fails at large angles. The angle-aware candidate passes all
seven resolved cases with a worst reflection of 3.907% at 60°. A localised
x-max/y-max corner packet passes at 2.594%. A scaled x-y control, x-z/y-z
edges and x-y-z three-face corner pass at 4.142%, 4.154%, 4.142% and 3.644%
for the scaled 8-cell wavelength, but its 10-cell repetition fails at
6.076%, 6.073%, 6.075% and 5.999%. A coefficient-free second-order-
extrapolated angle-aware candidate then passes both matrices: its x-y/x-z/
y-z/xyz ratios are 2.418%/2.439%/2.420%/1.846% at 8 cells and
3.192%/3.179%/3.191%/3.152% at 10 cells. This closes the declared linear-
acoustic multi-wavelength 3-D matrix for the raw-extrapolation candidate only.
That candidate becomes unstable in the finite-vortex audit. A separately
named TVD-limited dilatation/curl candidate stays stable and gives 0.943%
vorticity error, but fails velocity/density contamination at 19.825%/243.230%
and fails its 10-cell acoustic matrix at 5.502--5.564%. No unified boundary is
qualified; vortical outflow remains open. A first 0.120343 s
same-grid domain audit passes at 1.092% common-field L1 and 3.865% worst
primary-output change. Explicit cached-direct pressure projection reduces its
expanded run from 477.6 s to 220.4 s without changing those metrics, passing
the 5 min development target. The official 0.707246 s/3,056-step audit is now
complete and passes at 1.2447% common-field L1 and 0.9901% worst primary
change, including enlarged-domain lateral clearance. Its expanded runtime is
1402.2 s, so the separate 5 min workflow target fails. The memory-intensive
route remains opt-in.

A no-fit component comparison has nevertheless been completed. The fixed
source/handoff plus the existing axisymmetric free-jet route was sampled at
the published average-concentration coordinates downstream of the handoff.
Across 108 pairs it gives MAE `2.688 vol%p`, RMSE `3.160 vol%p`, signed bias
`+2.687 vol%p`, Pearson `r=0.937` and FAC2 `0.519`. The high correlation does
not rescue the route: 39 of 66 model-positive 4 vol% classifications are false
positives, and the predicted 4 vol% centreline length is 3.413 m versus the
measured 1.65 m. The axisymmetric shortcut is therefore rejected.

The resolved 3-D field was also sampled without extrapolating outside its valid
overlap. Its 6% handoff and 2.52 m common domain support 54 of the 153 average-
concentration points. A deterministic endpoint replay exactly reproduces the
frozen maximum slice and its final-20%
receptor drift is at most 0.00687%, so the endpoint is a qualified steady
proxy. Endpoint-versus-average MAE/RMSE are `2.179/2.491 vol%p`, `r=0.733`,
FAC2 `0.556`, and 4 vol% sensitivity/specificity are `1.0/0.522`. The accuracy
gates still fail. It remains separate from the outdoor LH2 headline evidence.

A predeclared four-case sensitivity then separates handoff placement from
effective scalar mixing on 45 common points. At the 6% handoff,
`Sc_eff=0.5/0.7/1.0` gives MAE `1.670/2.113/2.529 vol%p`; a 5% handoff at
`Sc_eff=0.7` gives `2.185 vol%p`. This localises the larger sensitivity to the
unvalidated scalar-mixing closure, but the lowest-error case still misses FAC2
and 4 vol% specificity targets. No coefficient is selected from observations.

The resolved-source mapping has an additional fail-closed gate. A 45 mm
coarse-grid numerical preflight maps the nominal 10% mean-plane density
departure to a 13.342% local maximum; an 8% candidate reaches 11.562%.
A 6% candidate stays at 9.664%, but it is only nominated for grid/time/domain
sensitivity. The preflight reads no concentration observations and therefore
is not a calibration or validation score.

Audit command:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe tools\audit_hytunnel_ps.py `
  reference\datasets\hytunnel-ps-4668554\raw `
  --output outputs\hytunnel-ps-spatial-audit-2026-10-11-v23.json

.venv\Scripts\python.exe tools\validate_hytunnel_free_jet_handoff.py

.venv\Scripts\python.exe tools\audit_hytunnel_resolved_source_preflight.py

.venv\Scripts\python.exe tools\audit_pressure_open_boundary_reflection.py

.venv\Scripts\python.exe tools\audit_pressure_open_boundary_oblique_reflection.py

.venv\Scripts\python.exe tools\audit_pressure_open_boundary_sweep.py `
  --boundary-scheme characteristic_multidimensional_extrapolated `
  --output outputs\pressure-open-multidimensional-sweep-2026-10-10.json `
  --plot outputs\pressure-open-multidimensional-sweep-2026-10-10.png

.venv\Scripts\python.exe tools\audit_pressure_open_boundary_corner_reflection.py

.venv\Scripts\python.exe tools\audit_pressure_open_vortical_outflow.py

.venv\Scripts\python.exe tools\audit_hytunnel_plume_domain_sensitivity.py

.venv\Scripts\python.exe tools\score_hytunnel_resolved_transient.py

.venv\Scripts\python.exe tools\audit_hytunnel_source_mixing_sensitivity.py
```

The result is `component_evidence_conditional`.  It can support future
matched-coordinate obstacle/ventilation comparisons, but it is not pooled
with outdoor LH2 field metrics.

## 3. MultHyFuel D2.2 CFD context

The official MultHyFuel D2.2 report is hash-pinned at SHA-256
`30c5131e30c9158daaf14559fdcb7e35b736a338e51ed8b6f4df65e7b2a7303f`.
Its ExJet obstructed benchmark uses a 12 mm horizontal GH2 nozzle at 0.75 m,
38.0 ± 0.1 barg and 8.0 ± 1.0 °C.  The jet enters an approximately 3.5 m-long
array of 1.5 m by 0.1 m cylinders.  Concentration and velocity receptors are
reported at 1.4, 2.5 and 4.4 m with six transverse/vertical offsets.

The report compares CFX, FLACS, two OpenFOAM versions and FRED AEROPLUME.  It
states that FRED gives the closest centreline concentration for this case even
though FRED does not resolve the obstacles, while OpenFOAM 1812 gives the
closest centreline velocity decay and the CFD radial profiles expose obstacle-
induced asymmetry.  This is an important caution against treating a single
centreline match as proof of a resolved wake mechanism.

The report publishes plots but no machine-readable matched receptor/model
table, and it does not tabulate an exact measured mass-flow history.  DEGALI's
current free-field core also does not resolve the cylinder array.  Therefore
no new DEGALI-versus-CFD error score or ranking is computed.  A future ranking
requires numerical exports for the same source, geometry, receptor coordinates
and averaging operator; otherwise `compare_models` must exclude the result.

Audit command:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe tools\audit_multifidelity_cfd_evidence.py `
  reference\datasets\multhyfuel-d22\raw\MHYF_WP2_D22_Assessment_of_dispersion_for_high-pressure_H2_20221111.pdf `
  --output outputs\multifidelity-cfd-evidence-audit-2026-10-10.json
```

## 4. Manuscript use

The defensible IJHE narrative is now:

1. retain PRESLHY E3.5 and FFI/DNV as the quantitative LH2 lanes;
2. use USN for the reproduced dynamic operator and the separately scoped
   retrospective enclosure-average parameter holdout;
3. use HyTunnel as conditional spatial/ventilation/geometry evidence;
4. use MultHyFuel as a public CFD context and matched-comparison protocol;
5. do not combine GH2, LH2, confined, free-field, peak and time-average records
   into one headline metric.

This broadens the evidence base while preserving the strongest claim: DEGALI
is a fast evidence-gated LH2 screening path, and complex obstacle/transient
cases are escalated to matched CFD rather than being hidden inside a fitted
wake factor.
