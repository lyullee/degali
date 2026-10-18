# TI021 3D modelling and validation iterations

## Objective and acceptance

Continue the existing active session goal: 600 s cumulative dose MAPE <=10%
on all 84 positive observed TI021 sensors, with all 145 sensors retained in
the output, no per-sensor gains, mass/momentum conservation, and no degradation
hidden by route/height/distance aggregation. Historical goal wording references
44.560%; the reproducible layer baseline is 29.544145611769668%.

The user also requested a 30 s end-to-end runtime. Record measured runtime
for every candidate, including the inner solve, interface and exterior solve.
Do not equate a 30 s physical prefix with 30 s wall-clock execution.

## Iteration order

1. Freeze source geometry, species identity, coordinate transforms, sensor
   heights, clock and units. Check initial inventory against the exact prism.
2. Fix implementation errors before adjusting turbulence/obstacle coefficients.
3. Use 10--30 s physical prefixes to check upward/downward velocities,
   above-wall species flux, reverse streamwise flow, positivity, conservation
   and coupling support. The workbook has no experimental concentration
   histories for these prefixes; do not score prefix integrals against 600 s
   observed doses.
4. For each physical candidate change one mechanism, execute the whole 600 s,
   score the same sensors and report height/distance/route groups. A local
   sensor subset cannot establish the full 84-sensor target.
5. Use conservative overlap of spatial cells/faces and exact time integration
   for the inner/outer interface. Export the full species, mixture, momentum
   and enthalpy inventory, including lateral/top exports and reverse fluxes.
6. Adopt only a candidate meeting the original conservation, group and runtime
   requirements. Keep the best baseline until a candidate qualifies. Target
   completion requires verified <=10%; a model limitation does not complete it.

## Findings on 2026-10-17

The previous `score_ti021_3d_obstacle.py` result of 94.1604% was an invalid
R12 vol% reference. It divided species kg/m3 by air density, overstating the
volume concentration by Ms/Mair (about 4.17). The correct factor at 291.15 K
and 100500 Pa is 19.9210360352 vol% per kg/m3. The previous solver also allowed
downward species loss through z_min. Its horizontal carrier was reset to the
wind each step, so that result did not resolve pressure-driven overtopping or
separated wake momentum.

The diffusion operator's neighbour arrays were both equal to the centre field,
so diffusion was identically zero. This has been corrected with conservative
face pairing, separately blocking each solid interface. Regression tests
verify spreading, positivity and mass preservation adjacent to an obstacle.

A sensor-envelope domain started at x=60 m, exactly through the holder centre.
An initial-cloud candidate consequently included only 1645.06 kg of the true
3290.12 kg species inventory. Source anchors now extend the domain upstream;
initial-cloud truncation fails a new inventory check. The retained diagnostic
`ti021-3d-fullcloud-angle61-2026-10-17.json` includes 3290.1158326486 kg with
relative geometry/inventory error -3.77e-14 and MAPE 91.8026976684%. It remains
a simplified scalar transport candidate, not a promoted model.

Preserve intermediate JSONs for the audit trail. The partial-inventory
`ti021-3d-repaired-initial-angle61-2026-10-17.json` is rejected. Scalar
`ti021-3d-fullcloud-vectorized-2026-10-17.json` verifies numerically identical
interpolation but was executed while the scorer report schema was being
edited; use a fresh final-code run for provenance/acceptance.

## Pressure-coupled 3D route

Reuse the existing binary finite-volume projection solver rather than
interpreting the simplified velocity closure as a wake solution. The local
domain is 140 x 140 x 40 m at 5 m isotropic spacing (6272 cells), initially
quiescent 14 m-wide / 13 m-high R12 prism, 3.4 m/s wind at 6.1 degrees,
6 m-high zero-thickness arc wall. It includes 36 measured sensor positions;
5 m vertical cells still alias the 0.4 and 2.4 m heights.

The 30 s prefix completed in 12.50697 s. It passed the execution gates,
with relative species residual 4.15e-16, continuity residual 6.54e-12 kg/m3/s,
and momentum balance residual 1.33e-9 kg m/s. The final near-wall vertical
velocity ranged from -0.873 to +1.352 m/s; two lee cells had reverse u, with
minimum -0.205 m/s. These are model predictions, not experimental validation.
41.33% of forward species mass through the downstream plane was above 6 m;
that plane fraction does not uniquely identify the wall-crossing route.

The matched 36-sensor subset of the layer baseline has MAPE 28.62156%.
Only a completed 600 s local run may be compared with that number. The new
`audit_ti021_3d_iteration.py` withholds prefix MAPE and makes subset coverage
explicit. A downstream plane alone cannot establish a full exterior mass
ledger: lateral/top losses require additional interface exports.

The 600 s local binary run completed in 98.53046 s with all internal execution
gates passing. Its 36-sensor dose MAPE is **87.75704%**, worse than the same
36-sensor layer reference (28.62156%). Height-group errors are 94.5831% (0.4 m),
88.6455% (2.4 m), 75.9030% (4.4 m), and 96.1775% (6.4 m). It is rejected for
promotion. This local subset cannot be compared directly with the scalar
91.8027% aggregate across 84 sensors.

The next paired candidate uses the already-declared stable-density vertical
diffusivity sensitivity (alpha=5, exponent=1, shear regularization=0.041/s),
with every other setting frozen. These parameters predate this scoring
iteration (see the 2026-10-12 RI-damping audit) and are a sensitivity test,
not fitted validation constants. Score the entire 600 s on the same 36 labels.
Continue with the full exterior interface only if there is verified benefit;
both this candidate and the baseline still require near-ground resolution
checks and full sensor coverage before promotion.

The frozen RI candidate completed the full 600 s in **104.63806 s** with
36-sensor MAPE **85.89580%**, an improvement of 1.86125 percentage points
relative to the paired binary baseline. All height groups improved, but
normalized bias is still -0.90733 and the layer model is better on this same
subset. Keep the damping as a research candidate, without promotion. The
next cause to examine is the strong near-ground underprediction and missing
residence, including 0.4/2.4 m aliasing in 5 m cells.

The interface adapter has an additional concrete mapping defect: its callback
returns `[mixture, species, lateral_momentum, vertical_momentum]`, whereas the
receiving layer interprets components as `[mixture, species, px, py]`. Thus
vertical momentum is treated as horizontal py and normal px is absent. Old
coupled scores cannot establish an admissible 3D handoff. Correct the donor
packet contract and receiving component names before running a coupled
acceptance candidate. Do not invent missing normal momentum from the wind.

Final-code scalar reproduction completed in **55.53873 s** including output,
with all 84-sensor dose MAPE **91.80270%**, exactly matching the prior
full-inventory scalar result. `outputs/ti021-3d-finalcode-600s-2026-10-17.json`
contains frozen input/code hashes, all sensor rows, group scores, reproduction
arguments and traces. Sixteen targeted regression tests pass. This scalar
candidate lacks a momentum equation and is not accepted. The full goal remains
active and unachieved; neither <=10% nor <=30 s has been established.

## Reproduce

From the repository with its existing virtual environment:

```powershell
& .\.venv\Scripts\python.exe -m pytest tests\test_transient_dense_gas_3d.py tests\test_score_ti021_3d_obstacle.py -q
& .\.venv\Scripts\python.exe tools\score_ti021_3d_obstacle.py --output outputs\ti021-3d-finalcode-600s-2026-10-17.json --source-mode initial_cloud --wind-angle-deg 6.1
& .\.venv\Scripts\python.exe tools\run_thorney_ti021_inner_transient.py --workbook ti21_input_ascii.xls --duration-s 600 --cell-size-m 5 --inner-domain-m 140 140 40 --handoff-site-downwind-y-m 70 --time-step-s .125 --progress-interval-s 60 --wind-speed-m-s 3.4 --wind-heading-relative-deg 6.1 --circular-arc-representation zero_thickness_surface --obstacle-height-m 6 --output-directory outputs\ti021-binary-local600s-2026-10-17
& .\.venv\Scripts\python.exe tools\audit_ti021_3d_iteration.py --directory outputs\ti021-binary-local600s-2026-10-17 --output outputs\ti021-binary-local600s-audit-2026-10-17.json
```

Inspect live process handles and progress before invoking a reproduction
command. The resolved runner rejects existing output files to prevent an
accidental restart/overwrite. Start another candidate in a new directory.
