# Decision-oriented LH2 scenario analysis (2026-10-03)

## Purpose and claim boundary

This audit connects a modelled concentration contour to declared 20, 30 and
50 m screening lines. The cases are reference calculations, not records from
an operating plant. They do not define statutory separation distances,
probabilistic risk, permission to transfer LH2 or a site-specific design basis.

The reproducible command is:

```powershell
$env:PYTHONPATH = 'src'
.venv\Scripts\python.exe tools\audit_decision_scenarios.py `
  --output-dir tmp\decision-scenarios-reproduction
```

The frozen copy is in `outputs/decision-scenarios-2026-10-03/`.

## Reference cases

Two FFI/DNV horizontal releases represent transfer-like cases. Test 4 uses
0.828 kg/s through a 25.4 mm opening and the reported 5.0--6.7 m/s wind
interval. Test 6 uses 0.833 kg/s through the same opening and the reported
2.3--2.7 m/s low-wind interval. The observation for the decision screen is
the steady-model concentration at the downwind centre coordinate and 1.5 m
height. The local `corrections=False` reconstruction and corrected path receive
identical inputs. The former is not an independently run original DEGADIS
executable.

Four NASA pool-spill states provide a distinct storage/spill archetype. Their
pool-axis LFL reach and the lowest flammable height at 33.8 m are reported
separately from the fixed-height jet calculation. These different observation
operators are not pooled into one accuracy score.

## Transfer-line result

| Case | Wind interval | Historical 1.5 m LFL reach | Corrected 1.5 m LFL reach |
|---|---:|---:|---:|
| FFI Test 4 | 5.0--6.7 m/s | 22.7--31.6 m | 56.5--67.9 m |
| FFI Test 6 | 2.3--2.7 m/s | 12.6--14.5 m | 34.3--37.8 m |

Across both cases, three wind values and the 20, 30 and 50 m lines, 11 of 18
binary LFL-line classifications change between the historical reconstruction
and corrected path. For Test 6, all three winds change the 20 and 30 m
classifications; neither model reaches the 50 m line at 1.5 m. For Test 4,
the corrected model reaches the 50 m line at all three winds, whereas the
historical reconstruction does not.

Here a classification is true when the steady-model concentration at the
fixed 1.5 m receptor is at least 4 vol % at the named line, and false
otherwise. The 18 values are two cases × three wind values × three lines.

This is evidence that model formulation can change a declared screening
decision. It is not evidence that any of these lines is a safe or unsafe
operating distance. In particular, the corrected Test 6 calculation still
underpredicts the reported 30 m arc maximum, so its 34.3--37.8 m interval is
conditional rather than an operational clearance.

## Pool-spill result

The model reproduces the reported grounded/low/aloft regime in all four NASA
states. The mean absolute error in lowest flammable height is 2.47 m. The
pool-axis LFL reaches are 29.1, 48.0, 66.5 and 38.7 m. Only the first remains
inside the 33.8 m concentration-distance evidence boundary; the other three
are explicitly conditional extrapolations. Their trend must not be presented
as an externally validated hazard-distance curve.

## External-model gate

The existing SLABx Test 4 calculation is the only quantitative external-model
lane presently available at the same physical sensors. SLABx gives
MAE/RMSE 0.648/1.211 vol %-point and 3/0/0 LFL true-positive/false-negative/
false-positive detections; DEGALI gives 1.539/2.679 and 3/0/3. Ranking remains
conditional because the native source mappings differ and the DEGALI steady
275 s setting is not the synchronized 275 s mean produced by SLABx.

Published FLACS and ADREA-HF studies of “HSL Test 6” do not provide a matched
comparison with the unresolved FFI Test 6. The former is a different 2010 HSL
campaign at 0.071 kg/s through approximately 26.3 mm, whereas FFI Test 6 is a
2019 campaign at 0.833 kg/s through 25.4 mm. The number coincidence is not a
case match. Those CFD studies are useful for mechanisms--two-phase source
sensitivity, pool evaporation, air condensation, humidity and wind
direction--but not for a pooled accuracy ranking.

No licensed, same-input PHAST or EFFECTS result is present. Such a result can
be added only when source state, wind reference, release geometry, receptor
coordinates and temporal operator are recorded together.

## Derived files

- `jet_scenario_summary.csv`: fixed-height LFL reach and runtime for 12 jet runs.
- `screening_line_decisions.csv`: exact concentrations and classification at each line.
- `pool_reference_scenarios.csv`: NASA pool regime and distance results.
- `external_model_metrics.csv`: derived Test 4 SLABx/DEGALI metrics.
- `external_benchmark_register.csv`: publication and comparison-eligibility register.
- `decision_lfl_distances.svg`: figure for the decision consequence section.
- `manifest.json`: frozen operators, limitations and source hashes.
