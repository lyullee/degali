# DEGALI whole-model error audit — 2026-10-17

## Decision

The reviewed implementation is numerically reproducible and its principal
mass/momentum ledgers close, but the TI021 obstacle path is **not yet an
externally validated concentration predictor**.  The corrected physical
cavity-frame path has 27.7595566247% positive-observation dose MAPE, not the
10% target.  It is retained as a rejected internal engineering diagnostic;
the slightly lower 27.6831474976% legacy-frame result is historical replay,
not a physically preferred model.

## Verification completed

- Repository-wide Python compilation passed.
- After all repairs, including the cavity-frame, candidate-runner, legacy
  transient failure-provenance and corrected multilayer handoff, the current
  full suite completed with **2,430 passed and 145 skipped** in 1,134.88 s.
  No test failed.  The final affected-area run completed 101 focused tests
  with no failure.
- Repository-wide Python compilation and import of all 224 discovered
  `degali` modules passed.  `pip check`, the CLI entry point and
  `git diff --check` also passed.
- An isolated source/wheel build succeeded.  Installing only that wheel into
  a separate target and invoking the packaged `degali 0.3.0` CLI succeeded.
- The subsequent wall-flux provenance implementation completed 87 focused
  dense-layer/conservative-observer tests.  Replaying the promoted path with
  that code but the new option disabled produced bitwise-identical sensor
  doses.
- Five hundred randomized conservative-column reconstructions closed species
  volume to a maximum relative residual of 1.776e-15.
- The current 600 s candidate and independent replay have bitwise-identical
  time and concentration arrays.
- Regenerating the strict promotion audit directly from the three report
  artifacts reproduced every gate and score below; the saved decision is not
  a stale or hand-edited summary.
- Current-candidate conservation: maximum relative species residual within
  the 1e-8 gate, mixture residual within the 1e-6 kg audit tolerance, and the
  horizontal-momentum gate passed.
- After rejecting an experimental wall-flux trajectory, the production solver
  and scorer were restored to the exact audited SHA-256 values
  `7e0c2f21867c314b68fa882136215cf92d10f105b8a1bb4c925e9a17cb75502e`
  and `9343d6ed68d155bad1920126efc09aec212c7a62349a45a006f16bc1b9477936`.
  Ninety-seven focused solver/observer/replay/promotion tests passed after the
  restoration.  A fresh strict-audit regeneration again reports bitwise trace
  identity, all conservation/group/runtime gates passed, and external
  validation false.
- After the complete audit, a fresh baseline, corrected cavity-frame candidate
  and independent candidate reproduction were generated from the same current
  source.  The two candidate time, concentration and dose arrays are bitwise
  identical.  The strict current-code audit passes implementation identity,
  trace completeness/reintegration, conservation and the 30 s runtime gate,
  but rejects promotion because the 6.4 m group regresses by 0.31672472
  percentage point and the aggregate score remains above 10%.

## Defects found and repaired

1. **Finite-wall overtopping exponent was integrated one power too low.**
   A height-resolved openness profile `(H-z)/H` raised to exponent `p` must
   integrate to `raw_openness**(p+1)`.  The solver used `**p`, so `p=0`
   incorrectly transmitted the whole face.  Face flux and tangential bypass
   now use the integrated exponent and have exact analytical regression tests.
2. **The exterior packet path used the first wind component and a fixed R12
   molar mass.**  It now uses the declared positive carrier speed and the
   active closure's temperature, pressure, and species molar mass.
3. **Observation-component diagnostics were computed but omitted from the
   returned result.**  They are now returned, saved to the NPZ, and checked for
   pointwise closure against total concentration.
4. **Promotion evidence was incomplete.**  The audit now requires conservative
   observation-mode identity, current implementation hashes, complete 0--600 s
   traces, trace reintegration, and bitwise independent reproduction.
5. **The command-line runner silently dropped parsed physics options.**
   Top/side ambient entrainment, cavity-core entrainment reduction,
   wall-respecting horizontal diffusion, and reservoir capture mode are now
   forwarded.  Process-global option side channels were removed, and the
   unused duplicate cavity-start parameter was removed.
6. **The then-best configuration used secondary workbook inputs rather than
   the established final-report priority.**  Replacing the 6 m wall and
   3.4 m/s source wind with the report-priority 5 m and 3.9 m/s values reduced
   MAPE from 29.12388444% to 28.20152797%.  This is an input-provenance repair,
   not sensor calibration.
7. **The corrected wall law invalidated the historical exponent semantics.**
   Re-evaluating the physical ground-loading exponent after fixing the
   integral showed that `p=2` (integrated openness power 3), together with a
   0.50 m wall-impact marker width consistent with the conserved upper branch,
   reduced MAPE from 28.20152797% to 27.68314750% without degrading any
   preregistered height or distance group.
8. **True over-wall flux could not coexist with the conserved/hybrid impact
   split and lost its physical provenance.**  The arbitrary incompatibility
   guard was removed.  A conservative child marker now follows wall-crossing
   species and its age through transport, local cavity capture/return and wake
   capture without double counting the impact marker.  Its default remains
   off because the complete TI021 tests below did not pass the promotion
   gates.

## Additional defects found and repaired by the whole-model audit

1. **The residence-candidate CLI was not a baseline-preserving A/B runner.**
   `tools/run_ti021_residence_candidate.py` loaded the baseline and then
   unconditionally overwrote numerous settings with parser defaults.  Four
   early wall-flux trajectory files made through that route remain excluded
   from physical comparison because they changed more than the declared
   mechanism.  Optional overrides now default to `None` and are applied only
   when explicitly supplied.  A complete 600 s replay through the repaired
   runner produced the exact baseline time, concentration and dose arrays and
   zero configuration disagreements at 27.6831474976% MAPE.
2. **The cavity-backflow closure was tied to the global x axis.**  The mask
   could rotate with the wind while the carrier and momentum target still used
   `dx`, `dy` and `wind[0]`.  The 90-degree witness differed by 12.41--13.85%
   in mixture and 29.73--30.68% in momentum even though conservation closed.
   The closure, stagnation footprint and return geometry now share one local
   streamwise/crosswind basis, with an exact rotation-covariance regression.
   New runs default to `cavity`; `legacy_x` is explicit historical replay only.
   The corrected current-code 600 s TI021 result is 27.7595566247% in 29.586 s,
   0.076409 percentage point worse than the legacy benchmark and therefore is
   physically repaired but not promoted under the preregistered nondegradation
   rule.  An explicit `legacy_x` replay remains bitwise identical to the
   historical benchmark.
3. **Legacy transient observer failures were swallowed without provenance.**
   `src/degali/core/transient.py` caught `RuntimeError`, `ValueError` and
   `OverflowError` and returned an empty observer record.  It now records the
   failure stage, exception type and message on that observer; the aggregate
   result exposes `failures`, and the CLI warns that failed empty profiles are
   not zero-concentration predictions.  This path is not used by the TI021
   dense-layer benchmark, but the repair removes the silent safety ambiguity
   from the legacy transient route.
4. **The current TI021 scorer could not launch the multilayer observer.**
   `vertical_observation_mode='conservative_volume'` is the scorer default,
   but the input guard rejected every `vertical_layers>1` configuration even
   though the solver dispatches those configurations to the independent EOS
   volume-occupancy observer before consulting the single-layer mode.  The
   guard now rejects only the genuinely unsupported single-layer
   `handoff_gaussian` combination.  The regression test exercises the exact
   scorer-default/multilayer pairing.
5. **Multilayer ambient entrainment contradicted its documented boundary.**
   The solver states that ambient mass crosses the cloud top, but deposited
   every `dm` increment into layer zero.  New runs place it in the highest
   materially occupied layer; `legacy_ground` is retained only for explicit
   historical replay.  A two-layer witness proves that the corrected path
   dilutes the upper layer while the old path incorrectly dilutes the ground
   layer, with unchanged species and mixture ledgers.
6. **The multilayer scorer re-expanded a collapsed finite-prism handoff.**
   The integral model physically contracts the initial 13 m prism to
   6.57593994 m before the layer handoff, but the scorer passed the original
   source bottom/top back into the multilayer initializer.  That silently
   discarded the computed collapse and created a second, inconsistent source
   initialization.  Multilayer initialization now derives its centre, width
   and domain top from the actual handoff state; the report records
   `initial_vertical_gaussian_center_m`,
   `initial_vertical_gaussian_sigma_m` and
   `handoff_vertical_domain_top_m`.  Two regression tests cover both finite
   prism and Gaussian handoffs.

The explicit `NotImplementedError` branches found in boundary-only thermal
moments, unsupported phase quadrature and the unused `series()` port are
fail-fast capability boundaries, not silent numerical corruption.  They must
remain outside advertised supported configurations until implemented.

## Current result and why it is still rejected

| Quantity | Result |
|---|---:|
| Current-code baseline MAPE | 28.20152797% |
| Corrected physical candidate MAPE | 27.75955662% |
| Aggregate improvement | 0.44197134 percentage point |
| Independent replay | bitwise identical |
| Reproduction runtime including output | 28.9963 s |
| 30 s runtime gate | pass |
| Conservation gate | pass |
| 10% MAPE target | fail |
| All preregistered groups nondegraded | fail: 6.4 m +0.31672472 point |
| Promotion | rejected |

The earlier component-diagnostic run at the historical legacy-frame benchmark
shows that the unresolved column
core supplies 97.36% of predicted dose; all explicit overflow, reattachment,
marker, and slow-wake branches together supply only 2.64%.  The remaining
error cannot be repaired credibly by tuning those small branches alone.

Optimistic fixed-prediction diagnostics for that historical benchmark are
still above the target: mirrored-pair redistribution 22.78%, conservative
vertical redistribution within a tower 19.24%, one scale per downwind arc
21.46%, and one scale per horizontal tower 15.91%.  These are diagnostic
bounds, not admissible calibrated models.

The workbook also contains at least three positive-dose heights at 21
horizontal towers.  Thirteen of them (61.9%) have non-monotone arrival order
with height; their median vertical arrival range is 15 s and the maximum is
85 s.  A single instantaneous monotone column profile cannot reproduce that
structure.

## Rejected over-wall provenance candidates

The new child state identified 419--427 kg of species crossing finite wall
faces, so the mechanism is active rather than a dead option.  It was tested
without receptor-specific factors over the full 600 s event:

| Candidate | Dose MAPE | Runtime | Decision |
|---|---:|---:|---|
| Shared impact-marker crosswind and Gaussian | 27.2962% | 25.16 s | reject: 0.4/6.4 m and near/far groups regress |
| Wall-top 5.5 m narrow Gaussian | 31.6872% | 26.08 s | reject |
| Column-conservative vertical kernel with shared crosswind | 27.7536% | 28.99 s | reject: 6.4 m and middle groups regress |
| Column-conservative kernel with 5 m neutral-log upper carrier | 28.2162% | 31.28 s | reject: score, groups and runtime regress |

These results reject both a narrow fixed-height overtop Gaussian and an
independent fast upper carrier as the missing route to ten percent.  The code
path remains opt-in for other cases; it is not enabled in the promoted TI021
configuration.

An additional isolated age-dependent wall-flux trajectory was tested over the
complete 600 s event.  It produced **30.11864678% MAPE** in **30.9594 s**,
worsened the 0.4, 2.4 and 6.4 m groups and exceeded the runtime limit.  The
experimental observation mode was removed after rejection; its disabled path
does not remain in the production solver.

Broadening the already conserved impact branch was also tested without adding
a receptor factor.  A 20% split centred at 3.5 m with 1.5 m sigma reduced the
aggregate MAPE from 27.68314750% to **27.05109672%** in **28.0146 s**, but the
2.4 and 6.4 m groups regressed by 0.1337 and 0.3390 percentage point.  A 40%
split centred at 4 m with 2 m sigma scored 27.11384369% and regressed 6.4 m by
2.8572 points.  These are diagnostic evidence that vertical allocation is
involved, not promotable corrections.

## Unresolved model limitations

1. **Core material age:** the diffusive vertical-tail width uses solver clock
   time for the inherited core rather than a transported core species-age
   moment.  This is harmless only for a truly instantaneous, same-age cloud.
2. **Multilayer injection:** the current instantaneous finite-prism benchmark
   now initializes a conservative uniform inventory over the computed handoff
   height, but a later time-varying source/interface callback still does not
   carry a vertical profile.  Production multilayer source histories therefore
   need a vertically resolved injection contract; the repaired initial
   handoff does not solve that separate callback limitation.
3. **Legacy unsafe option remains representable:** wall-blind horizontal
   diffusion is retained for historical replay.  Production promotion rejects
   it, but low-level callers can still request it explicitly.
4. **Vertical dynamics are reduced:** the accepted fast path reconstructs a
   conservative column but does not solve a full vertical momentum equation or
   resolved 3-D separation/roof-overflow field.
5. **Zero-dose sensors:** MAPE is evaluated on 84 positive observations because
   percentage error is undefined at zero.  All 61 zero-observation sensors
   receive a positive predicted dose in the current run (median 0.3215,
   maximum 11.7839 and total 87.2912 vol%-s).  This is not visible in the
  historical 27.683% MAPE and still requires an explicit detection-limit/
   false-positive absolute-error gate for safety qualification.
6. **Single-case identification:** repeated TI021 use can identify mechanisms
   for this experiment but cannot establish external predictive validity.
7. **Multilayer species-bound projection:** the corrected seven-layer run
   conserves total species, but its layer-local `0 <= species <= mixture`
   repair cumulatively redistributes **318.8958 kg** between height layers.
   The same material can contribute more than once, so this is not a lost-mass
   figure; it is nevertheless too large to treat the vertical allocation as
   a passive round-off repair.  Until the mixture/species face fluxes preserve
   this invariant without projection, multilayer height predictions remain
   experimental even when the global ledgers close.

## Authoritative artifacts

- Previous strict benchmark, replayed by the current implementation:
  `outputs/ti021-report-wall5-wind3p9-e025-k002-currentcode2-600s-2026-10-17.json`
- Current candidate:
  `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-currentcode2-600s-2026-10-17.json`
- Independent reproduction:
  `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-currentcode2-independent-600s-2026-10-17.json`
- Strict audit:
  `outputs/ti021-report-wall5-wind3p9-overtop-p2-markersigma0p5-currentcode2-promotion-audit-2026-10-17.json`
- Strict audit regenerated after this review:
  `outputs/ti021-current-audit-rerun-2026-10-17.json`
- Corrected cavity-frame current-code 600 s diagnostic:
  `outputs/ti021-cavity-frame-cavity-currentcode5-600s-2026-10-17.json`
- Explicit legacy-frame bitwise replay:
  `outputs/ti021-cavity-frame-legacy-currentcode4-600s-2026-10-17.json`
- Current-code baseline for strict comparison:
  `outputs/ti021-report-wall5-wind3p9-e025-k002-legacy-currentcode5-600s-2026-10-17.json`
- Current-code internal promotion audit:
  `outputs/ti021-currentcode5-legacy-promotion-audit-2026-10-17.json`
- Fresh current-code baseline after all repairs:
  `outputs/ti021-report-wall5-wind3p9-e025-k002-legacy-currentcode6-600s-2026-10-17.json`
- Corrected physical candidate and independent reproduction after all repairs:
  `outputs/ti021-cavity-frame-cavity-currentcode6-a-600s-2026-10-17.json` and
  `outputs/ti021-cavity-frame-cavity-currentcode6-b-600s-2026-10-17.json`
- Strict current-source promotion audit:
  `outputs/ti021-cavity-frame-cavity-currentcode6-current-baseline-promotion-audit-2026-10-17.json`
- Repaired candidate-runner baseline-preservation replay:
  `outputs/ti021-residence-runner-baseline-preservation-600s-2026-10-17.json`
- Pre-handoff-repair 13-layer top-entrainment diagnostic (historical only):
  `outputs/ti021-multilayer-top-entrainment13-d028-600s-2026-10-17.json`
- Pre-handoff-repair three-layer top/legacy A/B diagnostics (historical only):
  `outputs/ti021-multilayer-top-entrainment3-d028-600s-2026-10-17.json` and
  `outputs/ti021-multilayer-legacy-ground3-d028-600s-2026-10-17.json`
- Corrected seven-layer uniform-handoff diagnostic:
  `outputs/ti021-multilayer7-handoff-top-hat-currentcode8-600s-2026-10-17.json`
- Final-code single-layer equivalence reproduction:
  `outputs/ti021-cavity-frame-cavity-currentcode8-600s-2026-10-17.json`
- Component attribution:
  `outputs/ti021-overtop-integral-fix-e025-k002-component-attribution-2026-10-17.json`
- Report-priority residual structure:
  `outputs/ti021-report-wall5-wind3p9-residual-structure-2026-10-17.json`

## Next admissible development step

Do not tune another global concentration multiplier.  Because the current
event is a single instantaneous initial cloud, a core-age tracer alone would
be identical to solver time and cannot change this score.  The next admissible
mechanism must instead carry a time-lagged vertical recirculation/core state
capable of the measured non-monotone height arrival ordering, then rerun the
same 600 s, height/distance, zero-sensor, conservation, reproduction and 30 s
gates.  If that does not materially change the 15.9--22.8% optimistic residual
structure, the 10% target is not supportable for this reduced model without a
higher-fidelity 3-D carrier or case-specific fitting.

## Multilayer audit outcome

The three repairs above make the multilayer path executable and remove the
inconsistent source re-initialization, but do not make it an accuracy
candidate.  The earlier 13-layer **79.51059416%** and three-layer
**78.19617394%** results re-expanded the already-collapsed handoff from
6.576 m back to 13 m and are retained only as pre-repair history.  They must
not be used to characterize the current model.

With the current handoff height and the finite-prism model's conservative
uniform vertical inventory, the seven-layer complete 600 s result scores
**56.95036658% MAPE in 24.7616 s**.  Species, mixture, horizontal momentum and
vertical momentum ledgers pass (maximum species relative residual
`2.194e-15`; mixture residual `-5.59e-9 kg`; vertical momentum residual
`1.86e-9`).  Its layer-bound projection nevertheless redistributes a
cumulative 318.8958 kg, so global conservation alone is not sufficient to
qualify the vertical solution.  Gaussian initialization scored 63.95--67.94%;
adding simple global settling/lift closures worsened the uniform result to
59.42--68.27%.

The important inference is structural, not a calibration result.  Source
collapse is now preserved at initialization, yet replacing the entire
single-layer carrier with fixed vertical layers still destroys the strong
conservative core and lacks the resolved vertical momentum/exchange needed
for the measured non-monotone height timing.  The current physical
single-layer reproduction remains bitwise unchanged at **27.75955662%**.
Future work must therefore confine vertical resolution to the obstacle/wake
exchange and conservatively return it to the fast outer carrier, rather than
globally replacing that carrier.  The multilayer route remains experimental
and is not promoted.

## Local wake transport audit follow-up

The local conserved cavity reservoir contained one additional control-flow
defect: its advection and crosswind convergence were nested under
`return_rate > 0`. A physically persistent wake state therefore became
stationary whenever its transfer back to the main layer was disabled. The
transport and return operators are now independent. A regression holds the
return rate at zero, transports the reservoir into the wake corridor and
checks species, mixture and momentum closure. An opt-in representative-height
mode also evaluates its transport magnitude from the neutral logarithmic wind
while retaining the instantaneous wind direction.

Complete 600 s tests separate two conclusions:

- the best all-height/all-distance nondegrading dose diagnostic remains the
  visible Gaussian-moment wake at 26.77419422% MAPE, but its output-inclusive
  independent replay runtimes were 30.84 and 30.80 s and its mean departure
  error worsened from -23.56 s to -41.45 s;
- a genuinely moving, visible conserved layer with no return repaired mean
  departure error to -1.66 s, but scored 27.56005890% with the conservative
  layer observer (32.99477944% with the Gaussian observer) and regressed the
  2.4 m/near groups;
- adding a 66.7 s return time at the maximum local-layer transport speed
  reached 26.89297295%, but left tiny 2.4 m and far-distance regressions and
  an output-inclusive runtime of 30.01 s;
- a 2.5 m neutral-log transport basis did not remove that tradeoff.

Thus the repaired transport is retained as default-off model capability, but
none of these candidates is promoted. A single well-mixed vertical wake shape
cannot simultaneously reproduce the measured height allocation and temporal
duration. The <=10% objective requires a resolved or independently validated
multi-branch vertical wake carrier; adjusting this residence state further
would be case fitting rather than a justified closure.

Focused verification after this repair is **78 passed**; `py_compile` and
`git diff --check` pass. The earlier 2,430-test result predates this small
follow-up and is not presented as a rerun of the new code. The two complete
candidate traces and doses reproduce bit-for-bit under the current hashes;
the clean strict audit rejects promotion on runtime and the unmet 10% target.

## Sequential recirculation-to-reattachment audit

The local reservoir can now return a declared conservative fraction to the
reattached carrier while the remainder rejoins the prompt core.  The same
fraction is applied to mixture, species and both momentum rows.  A 0% split
control reproduced every sensor dose and all seven preregistered group scores
exactly; focused verification is **82 passed**, with `py_compile` and
`git diff --check` clean.

Sending the entire return to a fast height-based carrier was decisively
rejected: complete-event MAPE was 44.83--47.52%, despite closed conservation,
because far-ground exposure was greatly overpredicted.  With the original
slow reattached carrier retained, a 20% split and a 1.7 m Gaussian centre
instead produced **26.68333974% MAPE**.  All four height and all three distance
groups improve relative to 26.77419422%, and an independent run reproduces the
score exactly.  Maximum relative species residual is `2.011e-15`.  The
output-inclusive independent runtime is **32.3572 s**, so this is an internal
accuracy candidate rather than a 30 s promotion.

The component audit closes dose to `2.274e-13 vol%-s` and attributes 96.06%
of positive-sensor predicted dose to the unresolved core.  An exact linear
program gives minimum MAPE values of 26.5728%, 24.9686% and 24.9078% when all
existing component amplitudes are respectively bounded to 0--1x, 0--2x and
0--5x.  The <=10% target is therefore unreachable by event-wide branch gains.
Even optimistic mirrored-pair, tower-vertical, downwind-arc and per-tower
redistribution bounds remain 21.7434%, 18.0581%, 21.3930% and 16.0133%.

Wind-history perturbations do not explain the gap.  Changing the spatial
advection speed from 3.2 to 3.9 m/s worsened MAPE to 26.8108%; +/-3 degree
direction offsets worsened it to 30.6211% and 30.2064%.  A 2.5 m/s history
cannot cover the complete 600 s layer interval and is rejected by the input
guard rather than extrapolated.  The remaining defect is structural: the
single-layer core lacks the spatially varying vertical and crosswind exchange
needed by the observed tower pattern.  No sensor-specific correction was
introduced.

New evidence:

- `outputs/ti021-sequential-strict-f020-c1p7s0p5-independent-600s-2026-10-17.json`
- `outputs/ti021-sequential-strict-f020-c1p7s0p5-components-600s-2026-10-17.npz`
- `outputs/ti021-sequential-strict-f020-c1p7s0p5-component-span-audit-2026-10-17.json`
- `outputs/ti021-sequential-strict-f020-c1p7s0p5-residual-structure-2026-10-17.json`
