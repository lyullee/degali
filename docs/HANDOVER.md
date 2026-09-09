# Handover

> **LOCAL NEXT 2026-09-08 — conservative thermal-moment downstream reach:**
> The guarded Trial 10 state-space march reaches the first profile station
> from x=.79661 to x=1.78414 m; the half-step run reaches x=1.78074 m. Terminal
> parameter/value differences .001907/.003600 pass the .005 convergence gate,
> and all local weak/positive-diffusion/inflow gates remain valid. However,
> independent six-balance residuals are 9.206e-4 and 7.465e-4, both above the
> preregistered 5e-4 reach limit. No receptor score or promotion. Read
> `thermal-moment-downstream-reach-results.md`. Next: integrate the six physical
> flux/moment values directly and invert them to section state at every RK stage;
> do not lower the balance limit or correct the residual after marching.

> **LOCAL NEXT 2026-09-08 — measured-wind displacement separated from
> intrinsic spreading:** Stored Trial 10 yaw/magnitude/vector trajectories were
> projected through the exact five-sensor operator. Direction-only yaw improves
> all four amplitude/load median log errors (thermal centre 31.4%, thermal M0
> 18.8%, H2 centre 13.3%, H2 M0 5.0%) but worsens thermal/H2 variance errors
> 36.7%/25.0%; all wind cells therefore fail the preregistered gate. Vector
> step .02/.01 converges at max 1.551e-4 versus the .005 limit, so this is not
> integration-step error. Read `wind-vector-profile-observation-results.md`.
> Next: specify and implement an opt-in conservative thermal second-moment
> state, keeping enthalpy M0 and wind displacement separate and fitting no
> coefficient. Local changes remain intentionally uncommitted/unpushed.

> **LOCAL NEXT 2026-09-08 — identical model/measurement profile operator:**
> Stored density-control and same-width enthalpy trajectories were projected at
> the exact Trial 10/23 five-point vertical crosses without reintegration. The
> enthalpy path slightly improves thermal-variance median log error
> 0.26833→0.25256 but worsens centre-deficit error 0.17924→0.18999; it fails the
> pre-registered componentwise dominance gate. All model thermal widths are
> already 1.026--1.252 times observed. Trial 10 at 4 m instead has thermal M0
> 2.30--2.40× and H2 M0 7.67--7.84× observed, so a width-only diffusivity change
> is contraindicated. Read `model-thermal-profile-observation-results.md`.
> Next: project the already stored measured-wind-vector/yaw trajectories through
> this same operator to separate lateral displacement from intrinsic thermal
> spreading; do not reintegrate or fit a yaw to temperature. Local changes are
> intentionally not committed/pushed until the user requests a batch update.

> **LATEST 2026-09-08 — paired thermal/species observation operator:**
> Official PRESLHY E3.5 Trial 10/23 workbooks were reduced read-only at the
> collocated x=1.78/4.00 m, y=0, five-point vertical crosses. The calculation
> was pre-registered, hashes every raw and executable input, fits no coefficient,
> and distributes no external data. Thermal/species widths are approximately
> equal at 1.78 m, but the thermal-deficit profile is wider at 4.00 m. Median
> squared-width growth ratios are 1.4723 (Trial 10) and 1.4501 (Trial 23).
> Do **not** use these temperature-deficit ratios as `D_h/D_C`; they reject the
> direction of the prior 0.8235 sensitivity but do not identify an enthalpy
> diffusivity. Read `observed-thermal-profile-moments-results.md` first. Next:
> project model output through this exact finite five-sensor operator and require
> simultaneous centre-amplitude and second-moment agreement before promoting a
> thermal closure. TKE remains unidentifiable from these workbooks.

> **USER RESUME 14:34 KST:** measured_vector `step001` executed to completion
> (`reference/preslhy/observed_wind_vector_step001_2026-09-06`, 1066 steps).
> Four-cell same-key re-audit complete:
> - CONTROL mean MAE = 27.0312 K.
> - yaw .02: mean MAE 19.0453 K (−29.54%), MG 0.932706, width 1.235037, centre MAE 0.177467.
> - magnitude .02: mean MAE 29.2721 K (+8.29%), MG 0.774623, width 1.093047, centre MAE 0.188774.
> - vector .02: mean MAE 23.7324 K (−12.20%), MG 0.871162, width 1.292207, centre MAE 0.200676.
> - vector .01: mean MAE 23.7332 K (−12.20%), MG 0.871186, width 1.292243, centre MAE 0.200680.
>
> vector .01 `.1e-2` behaves almost identically to vector `.02`; no new
> physics signal from step refinement alone. Next state: keep this as evidence for
> "wind-input uncertainty not in step-size," then move to the next non-repeated
> downstream/physical candidate when user re-specified.

> **USER PAUSE14:23KST:** no further work until user resumes. Original degadis
>was removed on earlier pause, recreated as degali after continue; removing
>degali now. Latest state/resume is at TOP of docs/status.md. New sealed
>yawed_crosswind.py and runner/preregs must not be edited. .02 yaw pilot and
>.01 resumed1066-step field complete; mean/median/C improved but Tmin/width
>worse, no promotion. Two observed-speed input cells .02 have terminal artifacts;
>inspect their complete/failure before any calculation. New8wind tests passed.
>Measured-vector .01 and four-cell audit NOT YET DONE; no full7 or default claim.

> **Newest13:40 Stage2:** read lateral-wind-vector-results.md FIRST.
>All runs ENDED; full regression36296 completed817PASS/128SKIP/58excluded,
>then lateral8 and ambient-vector8 tests separately PASS. Trial10 matching
>49sample mean wind-FROM276.304311deg, release75deg; trial23 far1Hz invalid.
>Do not fit yaw or pick coarse local wind by temperature. Actual fixed-field
>vector budget identifies missing Py, not field improvement. Next bounded3D
>candidate with six-flux handoff, zero-yaw consistency and same-sensor audit.
>Existing sealed artifacts immutable. Deadline15:50KST; no active calculation.

> **Newest13:01 Stage2:** read molecular-transport-scale-results.md FIRST.
>Molecular reference screen13059 ended,16tests PASS;7actual field replays,
>3247saved nodes/16235reference-q samples. Max reference width effect0.20402%,
>not physical upper bound/new field. Soret34gas evaluations/281unsupported;
>no physical alpha chosen. Two primary/reproduction source PDFs saved/checked.
>Full non-slow regression RUNNING36296: directory
>reference/preslhy/stage2_regression_2026-09-06_1300,246Python snapshot.
>Check complete.json/pytest.log FIRST. Freeze all src/tests/tools until ended;
>do not duplicate regression. All earlier solver/field sessions ended.
>Next bounded read-only task: matched left/right sensor symmetry error floor
>and wind/coordinate basis, not fitted yaw or repeated centre cold envelope.
>Then finish evidence-backed A--F delivery matrix. Deadline15:50KST unchanged.

> **Newest12:46 Stage2:** read thermochemical-ensemble-results.md FIRST.
>24970/33554 ended; original46group failures and refinement155intermediate
>objective failures retained. Final656objectives:623solver pass/581grid pass;
>finite-support diagnostic ONLY, no new dispersion/PDF parameter/adoption.
>Raw same41 TC mean audit complete, source XLSX unchanged;8new tests PASS.
>41mean-MAE18.56486K is a different statistic, not an improvement on original
>minimum17.33598K. Near-core mean residuals remain12.3/23.6K at10/23x1.19m.
>All new code/tools/preregs sealed in their result hashes: version, do not edit.
>No active calculation. Next molecular cross-transport scale screen; do not
>repeat LP/phase/sourceKE/TKE audits or choose free diffusivity/PDF variance.
>Deadline15:50KST unchanged. All running statuses below are historical.

> **Newest12:10 Stage2:** read liquid-excess-potential-results.md FIRST.
>Opt-in nonideal liquid commonG/h65--77.5K,40new/265related tests PASS,
>local16/356 complete/no field. JSON liquid_excess_potential_2026-09-06.json
>seals new src/tests/tool/prereg: do not edit them without versioned audit.
>All prior snapshot/property/tangent runs ended;134frozen/runtime unchanged.
>Next bounded D is unresolved thermochemical mixing and averaging-EOS
>noncommutation: define mass/volume/time averages and actual source constraints
>before any candidate. No fitted PDF variance, no unsupported phase extrapolation.
>Check existing work first; source energy ledger and initialTKE bounds already
>done, so do not repeat them. Deadline15:50KST unchanged.

> **Newest11:50 Stage2:** regression59000 and tangent43956 ENDED. Full non-slow
>726PASS/128SKIP/58deselected; new tangent6tests ran AFTER it, not a732full run.
>Read actual-enthalpy-tangent-results.md FIRST:315actual samples=314initial
>resolved+1supplemental near-knot resolved at unchanged1e-3 gate. No new field
>or ratio/default selected; current ratio1 has correction0. Do not repeat these
>audits or report reference-dependent composition fractions as energy shares.
>Mixed-phase transport/common excessG and actual coupling remain incomplete.
>All running statuses below are historical. Deadline15:50KST unchanged.

> **Newest11:27 Stage2:** actual energy scale audit32263 ENDED,4fields/2201nodes;
> no dispersion improvement claim. New species_carried_enthalpy operator29tests
> PASS; same-D limit unchanged, cross-flux needed for unequalD to preserve
> species-reference covariance.27new energy tests,253related PASS. Full non-slow
> regression RUNNING59000, directory reference/preslhy/stage2_regression_2026-09-06_1130.
> Check pytest.log/complete.json FIRST; do not edit hashed existing Python
> sources/tests until completion. Read mechanical-thermal-scale-results.md and
> species-carried-enthalpy-results.md. Next actual single-phase scalar-flux
> tangent/effect screen, not arbitrary k or repeated completed thermal budgets.
> All earlier running field sessions ended. Deadline15:50KST unchanged.

> **Newest11:02 Stage2:** all Python runs ENDED. Trial10 bounded-liquid533steps
> COMPLETE and audited, no sensor improvement. Do not repeat or expand this
> candidate: fine-step NOT run,23unsupported/full7incomplete, no promotion.
> Read bounded-liquid-handoff-results.md and downstream-cold-envelope-results.md
> FIRST. New D audit82rows=41sensors×2flow choices,4field replays/22sections
> exact, no new integration.21new/197related tests PASS,134frozen unchanged.
> Source evidence+9coldcomposition inventory COMPLETE in cold-mixed-phase-evidence.md.
> Next bounded D: actual meanKE/thermal flux and upstream mechanical energy
> scale versus near-core cold gap; distinguish local KE from advected upstream
> TKE. No fitted k/phase-delay, no repeated synthetic audit. Deadline15:50KST.

> **Newest10:27 Stage2:** bounded-liquid adapter22tests/176related PASS. Same
> source actual10boundary accepted,23unsupported60K. Trial10 downstream RUNNING
>19255;10-step checkpoints under reference/preslhy/bounded_liquid_trial10_2026-09-06.
> Read docs/bounded-liquid-handoff-results.md and complete.json/failure.json FIRST.
> Do not duplicate this533step run or edit code hashed by its inputs.json.
> After completion compare ONLY matched trial10 keys; whole7/38/17/41 remain
> incomplete. Probe domain errors are not proven physical crossings. Preserve
>64K guard/no subtriple phase splice; original2K gate unchanged. Deadline15:50KST.

> **Newest10:05 Stage2:** common liquid G module+audit COMPLETE,38new/154related
> tests PASS, frozen134/runtime and input hashes verified. No active process.
> Read docs/liquid-phase-potential-results.md and active log/batch FIRST.
> Next bounded step is real pressure-loss CONTROL10/23 handoff feasibility
> with the in-domain liquid candidate, same source and all5conserved fluxes.
> Stop/report unsupported lowT rather than graft pure solids at63.151K.
> Source formation remains unchanged for this partial contrast; not full EOS.
> Do not re-run completed caloric/property screens. No field improvement or
> promotion. Original deadline15:50KST unchanged; mixed-solid stability unresolved.

> **Newest09:48 Stage2:** beta-N2 common potential+corrected A-663 audit COMPLETE.
> Read docs/nitrogen-phase-potential-results.md, batch A--F and active log FIRST.
>31 new tests; related116PASS;8input hashes and frozen134/runtime verified.
> No active calculation, no field-score or default-promotion claim. Next common
> liquid enthalpy/Gibbs reference and bounded mixed-phase stability, then source
> coupling and real10/23 evaluation. Do NOT repeat completed nitrogen/oxygen
> property screens. Beta mean volume is approximate, not precise mixed-solid EOS.
> Existing heartbeat ACTIVE through15:50KST; A is partial, not wholly complete.

> **New user instruction:** all remaining work bundled for overnight execution.
> Existing degadis automation updated ACTIVE. Read
> docs/stage2-overnight-batch-2026-09-06.md A--F plus active log first.
> Continue between stages without asking for "next"; deadline15:50KST remains.
> This includes downstream heat/turbulence and final actual-field comparison,
> not just the next beta-N2 property module. Preserve partial/failed outcomes.

> **Newest09:19 Stage2:** gamma-O2 Gibbs potential+phase caloric audit COMPLETE.
> Read docs/phase-caloric-consistency-results.md and active log first.
>25 new tests; related85PASS; Stage1134 unchanged. Result
> reference/preslhy/phase_caloric_consistency_2026-09-06.json with code/source hashes.
> No active calculation. Next beta-N2 caloric/volume/fusion reference for59--63K,
> then mixed-solid/liquid consistency. Do NOT re-run completed oxygen screen
> or treat its pure-property reproduction as dispersion improvement. Deadline15:50KST.

> **Newest Stage2:** mixed-air-liquid screen COMPLETE session2268; no Python
> run active. New mixed_air_liquid.py +26tests, related60PASS; Stage1 unchanged.
> Local effect up to+3.624K at fixedh, NOT a field prediction. Read
> docs/mixed-air-liquid-screen-results.md then active overnight log FIRST.
> Next bounded task: consistent solid-liquid-vapor Gibbs/enthalpy framework,
> actual mixed-solid applicability. Do not splice at pure N2 triple point;
> NBS3921 liquid data and phase-screen work are already done. Deadline15:50KST.

> **Latest:** continuous spreading audit COMPLETE, not promoted; all prior
> sessions ENDED. Read continuous_spreading_audit_2026-09-06.json.
> Exact38/17/41 and9 candidate replays verified; step sensor delta<=.00392K.
> Next: bounded mixed-liquid N2/O2 phase screen (NOT yet a new full EOS).
> Existing closure treats species as independent pure condensates; quantify
> shared-liquid activity effects first above63.151K, with no low-T extrapolation.

> **Current:** continuous-spreading pilot65267 ENDED, full41 coverage,
> Tmin MAE17.56551K (CONTROL17.33598K), p0514.99937K, median21.84639K.
> No accuracy claim/promotion. Remaining5 .02 RUNNING63257;10/23 .01
> RUNNING60082. Finish aggregate+resolution audit, then leave this screen.
> Active log has exact paths; do not repeat completed matched-flow analysis.

> **Latest Stage2 heartbeat:** matched-flow group20918 COMPLETED; aggregate
> flow_effect_audit.json COMPLETE with14 checked replays and same38/17/41.
> New Coriolis MG.961072/VG1.169029/FAC237/38 but width.921984 and
> centerMAE.070635m; temperature minima worsen. No flow/default promotion.
> New continuous ambient-spreading option removes x<=1m derivative switch
> for the widths already evaluated there.19 independent tests pass. Original
> Fortran and Stage1 remain unchanged; sub-metre correlation validity is not
> established. Same pressure-loss density CONTROL trial10/23 field screen
> RUNNING65267; read docs/stage2-overnight-2026-09-06.md FIRST for follow-up.

> **NEW USER AUTHORIZATION: Stage2 physical error reduction + overnight work.**
> Read docs/stage2-overnight-2026-09-06.md FIRST. Deadline15:50KST2026-09-06.
> Existing degadis heartbeat updated to this scope, ACTIVE5min. No subagents.
> Stage1 frozen134 files unchanged; all hashes/runtime verified. Added matched
> source runner/protocol and15 passing tests. Coriolis density interfaces7/7;
> enthalpy5/7, failures10/25 retained (temperature2.1248/2.1287K versus2K gate).
> Pilot density10/23 COMPLETED session6814, fixed11 arcs and41 thermocouples,
> full coverage and balance5.058e-8/4.066e-8. No default promotion. Remaining
> density11/12/22/24/25 RUNNING session20918, directory
> reference/preslhy/stage2_coriolis_2026-09-06. Do not duplicate the run.
> Aggregate its complete outputs on original38/17/41 then inspect spatial
> residuals and root physics; see active log for values and follow-up details.
> Historical Stage1 stop instructions below are superseded by this request.

> **1차 COMPLETE — 2026-09-06, final delivery. ALL older blocks are HISTORY.**
> Read docs/stage1-results-2026-09-06.md and stage1-reproduction-2026-09-06.md.
>134-file stage1_delivery_manifest_2026-09-06.json verified with exact runtime.
> assessment127 upstream hashes, candidate10/23 all row/state/flux/interface
> comparisons to old .02 outputs EXACT0. Both cover fixed sensors<=6m and
> reach11.415/11.605m. BASE all7 original rows EXACT0;41 sensor correction
> disambiguates BASE vs old CONTROL. Final scope38arcs/17vertical/41temp.
>71899/92924 fine candidate runs and50842 resolution audit ENDED; candidate
> max concentration difference .0124%, temp .0295K; BASE .0634%/.210K.
> During schema review found BASE adapter double-counted stored ground image.
> Added opt-in reported_jet_trajectory.py, prereg-stage1-single-ground-image.md,
>13 manufactured tests, audit_stage1_single_image.py. No frozen core edits.
>4964 correction audit ENDED. Same38/17/41: correctedBASE MG1.102099/VG1.233174/
> FAC236/38,centerMAE.054195m,TminimumMAE27.251879K. Candidate minimumMAE16.473842K
> but MG.838826/FAC235/38/centerMAE.065565m: NOT PROMOTED. This remains true
> against both historical and corrected BASE. No same-source causal claim.
> Schema companion: original stage1_temperature output_columns named first5
> only, raw rows have12; full accurate definition in assessment and new runner.
> tools/stage1.py verify/baseline/candidate is delivery entry. Its fresh BASE
> smoke17094 ENDED all7 EXACT0 to correction audit. New delivery18 tests PASS.
> No active relevant Python process on final check. Reports/tools/tests/output
> now FROZEN by134 manifest: do not edit them in place. Status logs aren't frozen.
> Existing full448 regression and separate later7 source-coordinate tests are
> not summed with18 as a new full run. No overall plant/safety certification.
> USER allowed necessary new physics, but this bounded1st goal is finished.
> Do NOT restart TKE/normal-stress or other2nd-stage work without a new request.
> Existing degadis heartbeat PAUSED successfully by app tool after delivery.
> Final report opening was queued for this thread. Do not restart automation.

> **Stage1 06:18 KST progress (latest overrides older blocks):**
> Existing evidence audit COMPLETED: stage1_existing_evidence_2026-09-06.json,
>111 current hashes; all field leaf core/input hashes match. Only historical
> interface runner hash differs; recorded, not overwritten.38/17 aggregates
> recomputed to2.05e-16. Prior temperature CONTROL is not corrected BASE.
>23208/9488 completed candidate trial10/23 .02m end-to-end reproductions,
>534/553 states, reaching11.415/11.605m (required sensors<=6m), balances
>4.338e-8/3.779e-8. Need final row/state comparison to original, not only metrics.
>65307 COMPLETED stage1_temperature_2026-09-06.json: BASE all7 fresh runs
> reproduce concentration/vertical rows with zero error. Same41-sensor minimum
> mean absolute errors BASE27.80464, CONTROL17.33598, CANDIDATE16.47384K.
> Median-statistic MAE BASE26.26797, CONTROL21.59628, CANDIDATE21.86994K.
> Do not mix BASE with CONTROL or call median absolute error a mean MAE.
> Fresh .01m candidate resolution runs started for10/23 (sessions in current
> tool outputs); prereg-stage1-resolution-diagnostic.md limits this to2 cases.
> NEXT: verify these against .02m and BASE .1/.2m, compile full paired report,
> input/source differences and portable execution instructions/manifest.
> BASE uses Coriolis flow_mean; candidate uses pressure-loss source flow
> (trial10 189.425 versus284.340g/s). Environment/sensors identical but source
> assumptions differ; this is NOT an isolated enthalpy-only causal comparison.
> No physics/default edits in stage1 so far. Core current111 manifest frozen.
> Stage1 script files frozen by their completed audit manifests: don't edit
> audit_stage1_existing_evidence.py or audit_stage1_temperature.py in place.

> **USER SCOPE RESET06:05 KST (2026-09-06): finish FIRST MILESTONE.**
> Read docs/first-milestone-scope-2026-09-06.md FIRST. User approved all5
> delivery steps and then explicitly allowed necessary new physics. Do not
> interpret as blanket ban on physics or permission for unbounded expansion.
> Verify fixed operational baseline/latest full-field candidate, representative
> end-to-end sensor-distance rerun, paired38/17/41 errors, reproducibility and
> limits. Reuse valid completed evidence. Stop heartbeat after genuine1st delivery.
> Existing degadis heartbeat UPDATED to this bounded scope, ACTIVE5min.
>80443 ENDED448 passed,128 skipped,58 deselected991.39s. Later moving-frame7
> tests passed3.02s separately. No active28624/25377/80443 remains. Mobile
> source2/2 synthetic PASS105hashes, normal source100hashes/joint131 verified.
> No new model trajectory/rates may be inferred from those source successes.
> Current main work: inspect run_preslhy_source_ablation, completed enthalpy/
> ambient control field outputs and41-sensor temperature comparison; beware
> that its control is NOT corrected JETPLU. Complete truthful same-model pairing.
> No new first-milestone audit process has started yet. Old blocks below HISTORY.

> **IN-PROGRESS05:44 KST (2026-09-06): mobile source implemented.**
> New mobile_tke_normal_source.py and6 tests. Current28624 running; first4
> tests passed including independent actual6-direction moment Jacobian.
> Current25377 audit_mobile_tke_normal_source.py running, preserves both Q
> centers .02/2, fixed r2/3 and original6targets. Q.02 coarse error6.75e-16.
> Original fixed-center Q2 failure preserved. New all-state solve has registered
> bounds and exact geometry, moving buoyancy and Hsecond weight derivatives.
> New files frozen by25377: do NOT edit until audit ends; if changes necessary,
> preserve this output and register another version. No old rates returned.
> Finish tests/audit, verify hashes, then update this block/result docs and
> run full non-slow if relevant. Next downstream normal/pressure production
> and physically closed mixing/Q/dissipation/circulation still open. Do not
> copy synthetic inputs into real trials or score observations before gates.
> Normal-budget/source findings appended to finite-tke-transport-results.md
> and overnight log; previous statuses below are historical.

> **IN-PROGRESS05:35 KST (2026-09-06): next mobile source initialization.**
> Current heartbeat read latest logs, checked finite-Q90/source95 hashes; no
> old active processes. New primary Favre equations show normal momentum and
> stress work must not be conflated with Q advection; see new joint-budget spec.
> normal_stress_budget.py implements PSD completion and joint flux bounds.
>9 tests passed0.81s (initial exact-float equality test fixed before freezing).
>98683 joint_tke_normal_budget COMPLETED7/7 resolved,131 hashes match. At fixed
> shear: K/Emean<=.12 requires normal work N/Emean>=.06216–.09432; making N
> <=.01 Emean requires K/Emean>=.5577–.6920. These caps are illustrations,
> NOT actual k/errors; pressure compensation and derivatives NOT bounded.
> prescribed_normal_source.py adds BOTH normal momentum and work to the
> unchanged source target with explicit Rss/k and fixed-ambient-pressure choice.
>23187 tests ENDED4 passed64.90s.34432 exact manufactured source audit ENDED:
> Q0=.02 PASS; Q0=2 FAIL bounded QP, violation9.086e-4. This does not prove
> physical impossibility. New output prescribed_normal_source_2026-09-06.json.
> **NOW implement docs/prereg-mobile-tke-normal-source.md:** allow6 source
> degrees [logrho,logC,logA,loguc, Cq, Hq] to match6 original targets with
> explicit Q/normal ratio fixed. Exact geometry at every candidate. Preserve
> old failed fixed-center source, shape limits, source targets and all defaults.
> No processes currently active. Last full425, then4+9+4 separate unit tests;
> do not call this full442. No new real turbulence input or observed score.

> **CURRENT05:10 KST (2026-09-06): ALL current calculations ENDED.**
> Finite-TKE operator actual-energy audit passed (90 hashes); prescribed-Q
> source reallocation audit78086 now also completed PASS (95 hashes match).
> Source test Q flux2.606398839W is allocated WITHIN unchanged source E:
> H+meanKE changes -2.606399000W; total E change -1.612e-7W, six-moment error
>5.566e-12, refinement1.878e-14. Original projection is unchanged. This is
> MANUFACTURED source bookkeeping, not a physical initial-Q solution.
>8784 ENDED4 tests passed152.00s, AFTER full425/128skip/58deselected373.78s.
> Do not claim full429 yet. No active90354/77745/8784/78086/93133 remains.
> New frozen modules: finite_tke_transport.py, finite_tke_energy_difference.py,
> prescribed_tke_source.py; all old failed audits/defaults preserved.
> **NEXT:** determine actual initial-Q source and physical mixing/dissipation/
> circulation laws with primary-source applicability checks; implement only
> explicitly declared research candidates, then test source moments, all old
> scalar/momentum and new Q faces, PSD shear bound, positive production jointly.
> Q0=2, chi_k=.3, tau=.2 in the completed synthetic audit must NOT be copied
> to real trials. A PSD lower bound is not an actual k estimate. Retraction
> does not close Q, stress normal components or circulation amplitudes.
> Only a physically closed, gate-passing trajectory can be scored at measured
> distances. No new observed accuracy yet. Follow finite-tke-transport-results.md.
> Existing heartbeat degadis prompt updated to this next stage, ACTIVE every5min.
> No mandatory external-data/user-authority blocker established in this turn.

> **CURRENT05:00 KST (2026-09-06): finite-TKE operator implemented and verified.**
> New finite_tke_transport.py preserves ALL old C/H weak rows, adds9 Q rates
> at degree4, H source D only, Q source P-D, total H+meanKE+Q, ambient k once.
> All Q/chi_k/tau/ambient/circulation inputs remain explicit, NOT LH2 defaults.
> Session90354 ENDED: finite_tke_manufactured_2026-09-06 numerical_passed=True,
>90 SHA dependencies match. Matrix/right/rates refine<=2.38e-10, full actual
> new-direction H+meanKE+Q differences<=1.992e-7, spread1.167e-7. Ledgers pass.
> Its MANUFACTURED physical inputs FAIL Q boundary, PSD lower bound and positive
> production (see finite-tke-transport-results.md). DO NOT adopt/score these.
> Session77745 ENDED425 passed,128 skipped,58 deselected373.78s. New19 tests
> are INCLUDED in425; no active90354/77745/93133. Old defaults/failures intact.
> **Active8784:**4 new prescribed_tke_source tests; NOT included in full425.
> This boundary-only helper reallocates an EXPLICIT Q flux within the unchanged
> total source energy, retaining original six enriched moments/center/geometry.
> Targets are [M,H2,Pax,Etotal,B,Hsecond], NOT [M,H2,Px,Pz,E,B]. Never confuse.
> Next finish8784, retain results, then determine physically admissible initial
> Q/chi/tau/circulation together and recheck source, gradients, stress bound,
> all original gates. No actual measured-distance trajectory/score yet.
> Heartbeat degadis ACTIVE; follow this latest log, don't rerun old operators.

> **CURRENT HANDOFF04:30 KST (2026-09-06): ALL calculation sessions ended.**
>15967 exact_aligned_witness completed PASS: all retained rows on new exact
> geometry, edge3.5431/4.0000/1.9626%, minchiP5.49216, nonradial1.71e-15,
> cross-volume4.08e-9/8.83e-11, new-source energy8.83e-11. Actual new-direction
> E differences<=7.66e-8; four-way spread1.49e-7. This closes the numerical
> mismatch but remains a LOCAL circulation witness, NOT a physical flow law.
>69372 initial_shear_tke_bounds completed all7 resolved: minimum TKE flux
>10.514–11.722% of mean axial KE; no fitted constants, actual k or derivative
> not determined. Fresh SHA checks all match: exact geometry137, exact
> aligned148, shear bound123, initial bounds126. Old failures preserved.
> Full regression403 passed,128 skipped,58 deselected640.85s;3 later covariance
> tests passed separately. No active15967/69372/18505/12426/59111 remains.
> **NEXT concrete task:** implement and test the opt-in finite-Q=rho*k modal
> transport operator in docs/prereg-finite-tke-modal-operator.md. Add m+1 Q
> rates/weak rows, thermal source D not P, Q source P-D, total E includes Q,
> natural Q edge=k_ambient FM; independently check all ledgers and PSD bound.
> Initial Q, dissipation time, mixing and ambient k remain EXPLICIT inputs;
> do not fabricate a validated LH2 closure or score an unclosed trajectory.
> See docs/tke-transport-followup.md and docs/exact-geometry-and-tke-results.md.
> Heartbeat automation degadis stays ACTIVE every5min in this same task;
> prompt updated to this remaining work. Do actual safe work, no duplicate
> completed audits. No external blocker presently; defaults/observed scores unchanged.

> **Latest04:10 KST:** exact_geometry_energy session18505 completed PASS.
> Actual H+K errors5.10e-8–7.18e-8 on h/2,h x4/8,8/16; spread2.08e-8,
> analytic mesh1.86e-10, geometry complex-step1.33e-15. This verifies the
> COHERENT new field derivative; old failures remain untouched. Old direction
> still misses new source by-0.261643, explicitly not conservative.
> ACTIVE15967 audit_exact_aligned_witness.py re-solves ALL retained rows
> with four aligned flow amplitudes and rechecks new-direction actual energy.
> TKE audit59111 ENDED resolved=True: fixed OLD aligned witness implies a
> minimum TKE flux10.5571% of mean axial KE, from PSD shear-covariance bound,
> without fitted constants. This is not actual k/dissipation or an axial
> derivative bound; not tracking k is not automatically asserting k=0.
> TKE helper3 tests passed1.79s. Full regression12426 ENDED403 passed,
>128 skipped,58 deselected640.85s. Three TKE tests came later, not full406.
> Defaults/field scores unchanged. Geometry re-solve is the sole active audit.

> **Latest around04:00 KST:** ACTIVE18505 audit_exact_geometry_energy.py,
> final exact_geometry_energy_2026-09-06.json not yet written. Both fresh
> phase operators4/8,8/16 done; coherent expectedE'=5.121468434, refinement
>1.86e-10. Old direction is NOT re-solved and misses new energy source by
>-0.261643; this is expected and explicitly not a new conservative witness.
> Actual H+K value differences at h/2,h and two phase grids are running.
> Exact geometry primitive9 tests pass0.85s, paired exact K2 pass107.61s.
> ACTIVE12426 full non-slow regression; previous completed count385 only.
> Read-only frozen-hash audit11561 ended: all17 completed today JSON manifests
> match (aligned118, smooth131); no altered dependencies. Research TKE
> follow-up structure documented in tke-transport-followup.md using primary
> NASA TMR sources; no new k closure implemented. Defaults unchanged.

> **Latest around03:45 KST:** sessions36646/88384/35453 ended. Aligned four
> mass-mode witness necessary screen PASS: edge3.5431/4.0000/1.9626%, minimum
> chiP5.492, nonradial1.50e-15, rates refinement6.33e-10, cross-volume4.08e-9/
>2.67e-10. Still NOT circulation closure or adoption. Read new
> aligned-flux-and-energy-conditioning-results.md. Smooth high-precision H
> component with legacy wind still FAILS narrowly1.0603e-5 energy-scaled;
> exact width-root wind compared to OLD tangent differs4.865% (not field score).
> Exact root changes width5.50microns, wind7.50e-7, directional wind4.83e-6.
> NEXT: coherent exact-constraint geometry and its own derivatives in a NEW
> opt-in numerical research path; never mix new values with old tangent/default.
> mpmath1.3.0 installed in venv for diagnostics. Full regression385 passed
>226.31s; aligned2 and smooth5 passed later separately. All current sessions ended.

> **Latest around03:30 KST (2026-09-06):** all previously started audits ended.
> Paired common-grid actual-value FD still FAILS the unchanged energy gate;
> same-table inverse root polishing improves root residual1e-10 to3e-15 but
> DOES NOT resolve the FD (h/2,h errors6.06e-4/1.21e-4). Preserve both new
> JSON records and their dependencies. No active15257 or82395 remains.
> Scalar diagnosis points to cancellation in center enthalpy/wind finite
> values; this is not yet a proved complete error budget. Do not adopt witness.
> aligned_momentum_flux.py now implements a kinematic mass-circulation plus
> conservative momentum completion that retains the normalized scalar-stress
> direction. Four unit tests pass, but NO real-trial feasibility audit yet.
> Latest full suite is still376; nine later unit tests passed separately.
> Next: bounded aligned-flux feasibility screen; numerical energy discrepancy
> remains an explicit adoption blocker, not a reason to relax frozen gates.

> **Latest around03:00 KST (2026-09-06):** memory-safe streaming implemented;
> sessions51872/35153/36619 ended. Streaming coarse parity1.86e-12,
> matrix/RHS refinement1.488e-8/4.229e-9, but OLD LP witness weak error1.123e-5
> stays rejected. New phase-consistent candidate uses accurate volume matrices
> with fixed16/96 LP constraint mesh. Rates refinement6.33e-10, amplitudes
>2.23e-11; cross-volume weak errors4.08e-9/2.11e-10, heat1.05e-12, mass4.44e-15.
> Sampled edges/prod pass; stress nonradial20.71% is NOT resolved closure.
> phase_consistent_flux_witness overall False because separate-integral FD
> global ENERGY alone fails6.83e-5/1.196e-3; do not overwrite that result.
> ACTIVE15257: paired_energy_difference_2026-09-06.json audit, same actual
> fields on union±phase grid, local product differences then fsum. No analytic
> EOS derivatives, no arbitrary-precision claim. Originalh,h/2 must both pass;
>2h/4h are sensitivity only. Preserve all frozen code/spec/results while running.
> Full regression376 passed,128 skipped,58 deselected431.62s; paired2 tests
> passed later independently20.92s (not yet full378). See overnight live log.

> **Memory-safe continuation required (around02:35 KST):** session24094 was
> STOPPED deliberately after verifying its exact process identity (PID129800).
> The fully vectorized phase8/16 audit exhausted32GB physical memory and
> paged heavily. Partial JSON retained; no final result, not external blockage.
> Coarse phase4/8 already exposes witness weak residual1.123e-5>1e-8.
> NEXT: implement NEW streaming phase-angle audit in8-ray batches, preserving
> all retained weak/global matrices, RHS, advective Jacobian and fixed-witness
> diagnostics. Do not solve the rejected hard boundary moment system merely
> to audit the retained equations. Same quadrature/gates, no refit. Preserve
> the106 frozen dependencies and aborted partial evidence. All prior sessions
> currently ended. See latest overnight log. No field score/default promotion.

> **Latest around02:30 KST (2026-09-06):** read
> [conservative witness results](conservative-flux-witness-results.md) and
> overnight live log. Solenoidal hard boundary-moment pilots10 initial/failed
> both REJECTED: large nonphysical oscillatory stress despite tiny linear error.
> All original weak/global equations retained; no row omission in this variant.
> Separate admissible_flux_witness uses LP inequalities only as existence
> diagnosis, NEVER downstream closure. Initial correction0; failed field's
> fine normalized amplitude0.001502 gives independent face3.54/4.00/3.02%,
> positive production/radial chi. Rate selector change0.00718, stress nonradial
> fraction20.7%; no convergence/realizability/observed claim.
> Full regression session68446 ended371 passed,128 skipped,58 deselected203.34s.
> LP3 tests added afterward pass independently, not yet full374.
> ACTIVE session24094: audit_admissible_flux_independent.py, fully phase-split
>4/8 vs8/16 retained weak equations at FROZEN witness and actual-moment FD.
> Do not edit its hashed dependencies or duplicate the run. Other sessions
>7834/24852 ended. NEXT depends on independent audit: resolve its numerical
> failures before physical closure. A useful subsequent bounded alternative
> is divergence-free mass correction with momentum correction constrained to
> preserve the original normalized scalar-chi direction (see live log idea).
> All claims remain research-only; no full downstream/field/default change.

> **Latest completed downstream screen (2026-09-06, around02:00 KST):**
> [Results and failure diagnosis](enriched-short-segment-results.md).
> All7 initializations passed, but primary guarded segments6/7;10 fails.
> GroupB and diagnostic sessions3712/50374/49425/1250 all ENDED.
> Combined segments93hashes; independent trial10 exit93; rejected boundary-row
> prototype90; rate-free compatibility95: all match. Frozen files unchanged.
> Exact trial10 failure at0.444045mm: corner heat5.01189%, rate refinement
>6.44e-10. Boundary-row replacement (separate module) failed badly10/23;
> do NOT adopt or extend it as passed. All original weak equations were NOT
> retained by that rejected prototype. Its4 unit tests check algebra/failure
> reporting, not physical success. Rate-free compatibility6 tests pass.
> Trial10 failed-state local minimax LOWER bound only0.2977%: potential for
> flux redistribution, not proof of a global conservative solution.
> NEXT: separate divergence-free transverse mass/momentum flux prototype,
> retain ALL scalar weak/global equations, recompute production/stress direction.
> Explicitly distinguish anisotropic momentum flux from scalar viscosity.
> Original boundary-exit JSON physical position columns miss sqrt(2q0);
> see documented erratum, not used for any gate. No field score/default.

> **Latest overnight work (2026-09-06):** read
> [the live record](overnight-work-2026-09-06.md) and
> [all7 initialization results](coupled-initialization-extension-results.md).
> All7 independent local initializations pass;85 frozen hashes match.
> New enriched_segments.py implements explicit midpoint, full moving state,
> fresh current-phase mesh every RHS, explicit conditional chi_C scaling,
> and strict start/mid/end physical checks. No initial moment retraction
> downstream. Both scaling policies are hypotheses, not LH2 validation.
> Primary2/4/8 division tests and independent endpoint checks are running.
> GroupA10/11/12 completed:10 fails heat5% near0.44 mm on4/8 divisions;
>11/12 pass short primary segments. GroupB22/23/24/25 still running in3712.
>90 hashed dependencies are FROZEN; use new files for next experiments.
> Non-slow session81091 finished355 passed,128 skipped,58 deselected578.82s.
> NEXT: independently replay trial10 boundary exit, finish/combine groupB,
> distinguish boundary invariance issue from numerical resolution. No
> field score or default promotion. Prior pilot-only entries are history.

> **Latest work — passed simultaneous initialization pilots (2026-09-06):**
> read [these results](coupled-shape-initialization-results.md) first.
> `addons/coupled_shape_initialization.py` adds `FixedMeshCoupledTransport`
> and `SimultaneousShapeInitializer`. Every candidate/FD perturbation
> re-solves the full gauge-free modal-rate system; no frozen fM/fP fit.
> Fast vectorized quadrature freezes SEED phase cuts, not current EOS.
> It is an optimization mesh only. Separate moving-phase moment retraction,
> operator4/8 vs8/16 and65 fixed angles/radial16 certify final candidates.
> Preserve six physical moments, fixed centers/velocity/basis/chi_C formula,
> ratio1, radial conservative flux and immediate-shear-heating assumptions.
> Mixing interpolation32 was already checked against direct original formula.
> Objective blocks: normalized H2/H/P face residuals plus negative chi_P/chi_C,
> each scaled by sqrt(sample count); physical stress is NEVER clipped.
> Fit mesh4/24 max6 iterations, then rebuilt8/48 max3. Both pilots required
> one first-stage improvement; second-stage conservation retraction sufficed.
> The new momentum boundary5% gate strengthens local acceptance; keep prior
> failed judgments unchanged. Original all-seven0/7 refers to OLD shapes.
> `coupled_shape_initialization_pilots_2026-09-05.json` covers ONLY25/11,
> both independently passed. Start-date filename; completed09-06.
> Trial25 H2/heat/P (%):16.88/16.54/11.84→3.24/2.30/0.66;
> min checked chi_P+5.098 s^-1. Trial11:4.01/4.14/6.20→3.56/2.81/0.53;
> min chi_P changes about-6.455→+4.036. Trial23 negative sign not revisited.
> Max actual moment error1.998e-10, refinement2.341e-11; max rate/matrix
> refinement1.397e-8/6.504e-9. Linear/heat errors<=2.302e-10/1.925e-11.
> Curvature half-width<=1.993e-4. Dense face1936/2240 plus independent65
> rays and129² fields. Sampled, not global proof. Max logshape.093976/.099566;
> trial11 lies near the.1 trust bound, important for downstream continuation.
> New6 tests pass; full non-slow342 passed,128 skipped,58 deselected369.14s.
> Pilots421.49/410.66s, all78 frozen hashes match. All processes ended;
> partial JSON is a retained snapshot, not paused work. Do NOT edit hashed
> source/prereg/audit to silently rerun an altered experiment.
> NEXT: extend to10/12/22/23/24 in a NEW prereg/audit with unchanged gates;
> trial10 hasdegree6/37 states, so add operator-size coverage. Then guarded
> downstream work only with trust/phase/BC checks and explicit chi_C update
> assumptions. No field scores/default promotion. No need to add arbitrary
> redistribution/TKE terms before testing the remaining cases; pilot success
> does not establish those extensions will never be needed.

> **Latest work — rejected coupled modal transport (2026-09-05):** read
> [these results](enriched-shape-transport-results.md) first. Older entries
> are history. `addons/enriched_transport.py` adds `EnrichedModalTransport`,
> `PrescribedRadialMixing`, panel antiderivatives and actual shifted-moment
> verification. This is CONDITIONAL LOCAL RATES, not a full downstream solver.
> The q heat mode duplicates dlog(beta). beta_ref stays a coordinate origin;
> q-mode rates carry actual width evolution. Do not add independent beta_dot.
> State order:log(rho_c),log(Cc),log(A),theta,log(uc),x,z,C modes,H modes.
> Trial10 degree6 has37 states; others degree4 have23. Solve35/21 active
> rates from5 global and30/16 weak modal equations; x/z rates are kinematic.
> `old_rate_embedding` maps old gamma into heat coefficients via2*q_mode/beta².
> New angularly dependent radial conservative fM/fP integrate the actual
> BM/BP. Scalar fluxes INCLUDE angular gradients of Y and h=H/rho. Rebuild
> fK=u*fP-.5*u²*fM and P=-2*q*u_q*(fP-u*fM), never reuse old mass/stress.
> Scalar chi_C is a REQUIRED supplied phase-panel table. The audit uses old
> positive weak-baseline chi_C; there is NO downstream mixing update law.
> Energy assumption string `reduced_buoyancy_work_immediate_shear_heat` and
> explicit thermal/species ratio1; not a full pressure/Favre/TKE model.
> The radial conservative flux ansatz omits divergence-free transverse flow.
> Source/coflow reconstructs the actual pipe/HEM factory; phase force/work
> recomputed on new fields. Do not substitute synthetic source factories.
> `enriched_transport_{pilot,group_a,group_b}_2026-09-05.json` cover10,
>11/12/22,23/24/25. Original numerics4/7:10/11/24 fail chi_C interpolation8/16
> comparison, not their rate/actual-moment checks. Keep these failures.
> `enriched_mixing_input_verification_2026-09-05.json` separately passes7/7
>16/32/direct-formula checks on new samples; fine16 max error7.918e-6<1e-5.
> Combined file reports independently verified numerics7/7, adoption0/7.
> Max linear error3.331e-10, weak heat4.901e-11, rate refinement2.979e-7,
> matrix2.106e-6, independent actual-moment FD7.085e-9, step2.149e-8.
> Max matrix condition2101; no pseudoinverse selection. All faces sampled
> inward; curvature half-width<=1.993e-4. Global conservation alone is not BC.
> New edge H2/heat (%):10(15.68,15.86),11(4.01,4.14),12(8.78,8.35),
>22(16.03,16.04),23(15.09,15.82),24(6.71,6.58),25(16.88,16.54).
> Momentum edge defects6.20–20.17%. New face samples1920–2272, not the old
>4097 uniform frozen-transport sample. chi_P radial check17 samples per angle
> is necessary only. Trial11/23 chi_P minima negative; six cases fail scalar
> edge gates and11 also fails positive viscosity despite scalar<5%.
> `enriched_countergradient_witnesses_2026-09-05.json` independently checks
>9 angles/order16, same solved rates: corner-near264K points require
> chi_P=-6.455/-26.903 s^-1, reduced production=-7.70/-48.79 W/m³.
> These violate the SELECTED positive-viscosity/immediate-heating closure,
> not every possible turbulence model with backscatter. Never clip negatives.
> NEXT: simultaneous shape/transport boundary-compatible initialization;
> if required, divergence-free transverse flux redistribution with explicit
> boundary conditions. Consider velocity-shape/TKE extensions if reverse
> stress remains; current evidence does not identify one unique remedy.
> No new downstream integration or field38/17/41 scoring/default promotion.
> New12 tests;336 passed,128 skipped,58 deselected in386.52s concurrent run.
> Hashes:base66,input69,combined73,witness75; prior61 unchanged. Do not edit
> these hashed sources/preregs in-place for a new frozen experiment.
> All executions ended; partial JSONs are retained snapshots, not pending jobs.

> **Latest work — passed conservative edge refit (2026-09-05):** read
> [these results](edge-conservative-refit-results.md) first; older entries are
> history. All7 fixed-transport cross-sections now pass their independent
> moment/edge/physical gates. Do not call this downstream or field validation.
> `addons/edge_conservative_refit.py` adds `FaceSplitSquareMoments` and
> `ConservativeEdgeRefit`; prior hashed implementations remain unchanged.
> Face Ti/Y table intersections split the angular domain; each new shape and
> angle rebuilds radial phase knots. Screening1025/2049 plus8/8 vs16/16 panel
> orders verifies the original Gaussian reference in7/7 (max4.806e-10).
> Monotonic radial branches and face-event discovery are sampled, not proven
> globally. Phase thermodynamics and shape coordinates are unchanged.
> Refit uses original diagnosed degrees:10 degree6;11/12/22/23/24/25 degree4.
> Six targets/scales M,H2,axialP,totalE,Fb,physical enthalpy M2 are frozen.
> Step QP enforces their linearization and33×33 scalar bounds simultaneously;
> every nonlinear retraction uses rebuilt actual phase integrals. Final129×129
> shape checks and4097 uniform edge samples are independent. abs(log change)
> max0.0999, below original0.1; no degree/physical coefficient increase.
> Initial retraction<=4 steps; later edge improvements0..2 steps. All reach
> the4% fit target then pass5% independent gates (not optimizer stagnation).
> Independent edge(H2%,heat%):10(2.44,2.48),11(2.80,3.16),12(2.35,2.18),
> 22(1.10,0.95),23(2.32,3.63),24(2.18,2.20),25(1.66,1.16).
> Independent moment error max9.041e-11;8/16 difference max1.394e-11, gate1e-8.
> Files `edge_conservative_refit_{pilot,group_a,group_b}_2026-09-05.json`
> cover disjoint10/11,12/24,22/23/25. Combined file
> `edge_conservative_refit_combined_2026-09-05.json` has7/7 and61 matched hashes;
> each group has56. Old51 dependencies and earlier failed records are unchanged.
> `tools/audit_edge_conservative_refit.py` reconstructs actual measured pipe/HEM
> source and coflow, replaying the previous selected target; do not replace with
> synthetic source factories. `combine_edge_refit.py` rejects missing/duplicate
> trials, incomplete inputs and mismatched source versions. New7 tests;
> non-slow324 passed,128 skipped,58 deselected in288.10s (concurrent case runs).
> All executions ended. Partial JSONs are snapshots, not pending/paused jobs.
> Computational rebuilding is still research-oriented, not performance-tested
> full downstream integration. Runtime `.venv/Scripts/python.exe`, BLAS1.
> NEXT: new C/H shapes retain OLD Gaussian-derived fM,chiC,fK,u,geometry and
> gamma only conditionally. Reconstruct new 2-D mass/shear/heat transport and
> independent coefficient evolution before integrating them downstream.
> Do not insert these constant coefficients into the old8-state radial ODE or
> claim weak reservoir transport remains solved for the new shapes. No field
>38/17/41 scoring, practical safety qualification or default promotion yet.

> **Latest work — rejected fixed-transport square enrichment (2026-09-05):**
> Read [these results](edge-profile-enrichment-results.md) first. Older entries
> are history. `edge_enrichment.py` adds centre-zero, square-symmetric even
> Legendre log corrections to C/H, preserving their centre values/signs and
> the same exact phase inversion. Degree4/6 gives16/30 total coefficients.
> A,u,theta,x,z and baseline fM,chiC,fK are FROZEN. No new coefficient RHS,
> reconstructed transport, old8-state ODE insertion or receptor usage exists.
> Six intended constraints are M,H2,axialP,totalE,Fb,physical enthalpy M2.
> Flux fit65 cosine nodes, independent4097 uniform edge samples, sampled
> abs(log corrections)<=0.1. Scales are frozen; no observation fitting.
> `edge_enrichment_screen_2026-09-05.json`:14 attempts,13 fully evaluated;
> trial10 degree4 fails order96 moment retraction. Other fits stop when no
> acceptable improving step is found. Retain lower/higher-degree failures.
> Trial10 degree6,trial12 degree4,trial24 degree4 have both edge errors2–3%,
> below5%, but **all fail independent1e-8 moment gates: final0/7**.
> The96-order constraints look accurate(<=5.146e-10); independent192/384
> integrals expose conservation errors. Trial25 also fails adjacent1e-5.
> `edge_phase_quadrature.py` independently splits NEW phase rays in a1/8
> polar square wedge, rebuilding Ti/Y knots for every changed profile/angle.
> It supports sampled monotonic Ti up/Y down only, with no global proof.
> `edge_enrichment_phase_diagnosis_2026-09-05.json` selects the lower maximum
> edge-defect candidate per trial for POST-SCREEN diagnosis, NOT promotion.
> Original coefficients are unchanged.8/32,16/64,16/128 phase/angle rules:
> last changes2.282e-8..6.966e-8; actual residuals2.412e-6..5.925e-5, largest
> in net buoyancy. Mass/momentum residuals also remain; do not repair only F.
> Zero-profile reference errors9.028e-9..2.878e-8: only trials12/23 pass1e-8,
> so this independent quadrature is NOT7/7 final verified either.
> NEXT: split angular phase/face intersections or adaptive angular error
> control; then use those true moments and trust-region-constrained steps
> for conservative refitting. Do not increase degree/log limit/physical
> coefficients just to pass. Afterwards reconstruct new shape transport
> before any short segment or frozen38/17/41 observational scoring.
> Original pilot/screen hashes46; diagnostic51; old41 dependencies unchanged.
> Do not edit their hashed sources/preregs/inputs in-place for a new experiment.
> New16 tests; non-slow317 passed,128 skipped,58 deselected in58.48s.
> All executions ended; `.partial.json` files are retained snapshots, not jobs.
> Runtime `.venv/Scripts/python.exe`; use OPENBLAS_NUM_THREADS=1 if repeating.
> The diagnostic tool intentionally reads the original dated screen, not an
> arbitrary repeat path. Defaults, phase tables and prior scores are unchanged.

> **Latest work — weak reservoir thermal transport (2026-09-05):** read
> [these results](reservoir-thermal-moment-results.md) first. Earlier entries
> are history. New `addons/reservoir_thermal.py` provides a RESEARCH weak
> thermal-width RHS, not strong pointwise boundary closure or field validation.
> `ReservoirThermalMoments` requires D_h/D_C and the explicit string
> `reduced_buoyancy_work`. The audit selects D_h/D_C=1 and instantaneous shear
> heating as hypotheses, not validated constants. Core/defaults are unchanged.
> Boundary total energy is .5*u_a²*f_M; enthalpy boundary flux subtracts the
> reconstructed mean KE + stress work. Internal gradients remain unchanged;
> the natural boundary is enforced only in thermal moment weak forms.
> Source energy adds directly computed integral(u*b_s) to the existing wind
> inflow energy; never add shear production P again. This reduced convention
> is not full pressure/gravity/TKE thermodynamics. Initial added work is only
> 3.530e-6..2.313e-4 W/m, not a fitted compensation for the prior kW/m defect.
> `_PhaseForceView` evaluates the unchanged entrainment source with phase-cell
> force quadrature. Actual source/coflow reconstruction remains mandatory.
> Eight coordinates remain log(rho_c),log(Cc),log(A),theta,log(uc),x,z,log(beta).
> R2 supplies gamma (-0.633432..-0.300054 initially); beta is NOT fixed.
> The original boundary class's solve/receptor guards remain unchanged.
> `reservoir_thermal_segments_2026-09-05.json`: 7/7 initial numerical checks,
> 6/7 RK2 2/4/8-step segments. Length=min(.01 m,.05*min(sy,sn)), actually
> 1.759..4.600 mm. Trial 10 fails 8-step balance at 1.254e-5 (gate 1e-5).
> `.partial.json` is a retained mid-run snapshot, NOT an active/paused run.
> `reservoir_thermal_segment_refinement_2026-09-05.json` retains trial 10's
> 16/32/64 results: 64 balance passes but 32/64 state change 2.627e-5 fails.
> `reservoir_thermal_adaptive_verification_2026-09-05.json`: separate DOP853
> rtol 1e-8/1e-10 checks pass trial 10, balance 6.260e-8/1.558e-9 and endpoint
> change 7.388e-8. COMBINED 7/7 short segments pass; original failures remain.
> The adaptive verification script intentionally replays this frozen original
> seven-case screen. Do not treat it as a generic arbitrary-screen driver.
> Pointwise edge defects remain: heat 17.74–25.40%, species 12.36–18.19%,
> momentum 0.150–0.297%, normalized to local advective fluxes, NOT measured
> concentration/temperature errors. Total weak closure cannot certify shape.
> NEXT: enrich BOTH species and heat edge profiles / connect the ambient tail
> to reconcile specified boundary flux and internal gradients. Do not run
> field 38/17/41 scoring or promote defaults before that independent check.
> New 14 tests; non-slow 301 passed, 128 skipped, 58 deselected (59.19 s).
> Final source/input hashes: original 35, refinement 38, adaptive 41 matched.
> All main executions have ended. Runtime stays `.venv/Scripts/python.exe`.

> **Latest work — rejected shear/thermal closure (2026-09-05):** read
> [these results](shear-thermal-compatibility-results.md) first; older entries
> below are history. New `addons/shear_thermal.py` reconstructs reduced radial
> axial momentum flux, residual shear, mean KE production and stress work.
> Force density and D_h/D_C are mandatory. It is not full curved Favre RANS.
> `equilibrium_budgets` EXPLICITLY tests Q_H=P; production is not inherently
> heat. Unity diffusivity ratio is a rejected research hypothesis, NOT a
> selected model constant. No old source, EOS, boundary guard or core changed.
> `shear_thermal_compatibility_2026-09-05.json`: 7/7 numerics, 0/7 physical
> compatibility. R2 gamma=-0.306707..-0.121278 gives positive sampled mixing
> and small curvature but heat defect +17.750..27.472 kW per axial metre.
> R0 roots=46.825..136.099 lie outside the necessary positive-diffusion bands.
> Independent finite changes of actual H, KE and physical M2 pass 7/7.
> `shear_thermal_boundary_diagnostic_2026-09-05.json` is a POST-SCREEN diagnosis,
> not revised gates. Edge advection +266.8..351.5 and diffusion -328.6..-249.1
> kW/m fail to cancel. Modeled edges remain 169.8..267.4 K vs ambient near
> 289 K. The existing energy source assumes ambient enthalpy zero. This
> boundary inconsistency dominates the attempted thermal closure's defect.
> At fixed five sources, changing only Q=eta*P for any eta in [0,1] cannot
> close R0 anywhere in the necessary gamma interval. This is NOT a rejection
> of full TKE transport with a changed energy state/boundary model.
> NEXT: conservatively couple the finite Gaussian tail/ambient energy
> boundary and all heat moments. Do not just increase delta, overwrite edge
> temperature, fit D_h/D_C, or add P again to total energy. Only then attempt
> short downstream segments and the frozen 38/17/41 observation comparisons.
> New tests: 14. Non-slow suite: 287 passed, 128 skipped, 58 deselected in
> 54.53 s. No new field run/default promotion. All main runs have ended.

> **Latest work — conditional transverse mixing (2026-09-05):** read
> [these results](transverse-conservative-mixing-results.md) first.
> `addons/transverse_mixing.py` builds the affine family of eight rates:
> log(rho_c),log(Cc),log(A),theta,log(uc),x,z,log(beta_H), where Cc=rho_c*Yc.
> gamma=dlog(beta_H)/ds remains free; the origin at gamma=0 is NOT an adopted
> constant-width ODE. `ConservativeTransverseMixing` has no solve method.
> Transformed radial total-mass and H2 fluxes recover physical transverse
> velocity after adding boundary motion. Their residual defines a diagonal
> anisotropic diffusion tensor; counter-gradient conditions are not clipped.
> The seven sampled gamma intervals are nonempty NECESSARY conditions only.
> Thermal coupling requires an EXPLICIT D_h/D_C and supplies no heat source.
> `addons/phase_radial_quadrature.py` splits the monotonic phase path at table
> cells (182–200 panels). It preserves the original EOS and piecewise slopes.
> The old adaptive thermal response replays within 7.080e-7. Width changes
> invalidate a cached partition. Nonmonotonic branches/ground are unsupported.
> `tools/audit_transverse_mixing.py` reconstructs ACTUAL measured/HEM sources
> and coflow; all seven stored first source vectors and boundary moments
> replay exactly. Do not replace this with the dummy property factory RHS.
> `transverse_mixing_screen_2026-09-05.json` keeps the failed dense unsplit
> gradient reference; its all_numerics_passed=False is intentional history.
> `transverse_mixing_flux_change_verification_2026-09-05.json` independently
> differences physical flux VALUES at re-partitioned perturbed states: 7/7
> pass while all retained conservation/refinement checks also pass.
> Non-slow suite: 273 passed, 128 skipped, 58 deselected in 123.76 s.
> Next determine gamma with a thermal/mixing ratio and a consistent spatial
> mechanical-energy exchange model, then ground closure and short ODE checks.
> No gamma, D_h/D_C or thermal/mechanical source has been selected. No defaults,
> field trajectories or measured accuracy scores changed. All main runs ended.

> **Latest work — thermal moment operators (2026-09-05):** read
> [the operator results](thermal-moment-operator-results.md) first.
> New `addons/thermal_moments.py` supplies mean advective M2, its fixed-H2
> reduced derivative, the planar curved/moving rectangle weak balance, and
> the response to a SPECIFIED specific-enthalpy-gradient diffusivity.
> No actual D_h has been selected for the seven cases: reported diffusion
> values are divided by D_h, not physical heat-transfer predictions.
> The new transport/shape closure itself is NOT complete. Boundary-only
> solve/source/receptor guards in buoyancy_profile.py remain unchanged.
> `thermal_moment_operators_audit_2026-09-05.json` preserves seven GL2048
> phase-gradient failures (1.347e-5 to 1.254e-4 vs independent GK21).
> Finite edge/volume magnitudes are 0.85178–0.87422. Do not omit the edge
> or confuse its diffusion-only value with a total plume heat-loss rate.
> `addons/thermal_moment_quadrature.py` separately adapts volume and edge;
> use its status/error checks instead of treating high-order GL as certified.
> `thermal_moment_adaptive_verification_2026-09-05.json`: independent checks
> pass 7/7 (maximum component difference 2.432e-7; net 1.626e-6 vs 1e-5).
> An internal requested error is only an estimate, not a rigorous bound.
> This audit costs about 320k–367k local phase evaluations per section:
> reduce cost before putting it in every downstream stage; do not just
> revert to uncertified low-order GL. Final 17 code/input hashes match.
> 32 new tests; non-slow regression: 257 passed, 128 skipped, 58 deselected.
> Final focused/doc check: 69 passed in 12.25 s. All runs are finished.
> Next supply consistent transverse mean transport, thermal/species mixing,
> mechanical exchange moments and ground closure before any seven-case ODE.
> Dummy property/geometry replay sources MUST NOT supply entrainment RHS.
> Runtime remains `.venv/Scripts/python.exe`. Previous results below are
> historical; no defaults or observation scores changed in this step.

> **Latest work — buoyancy-constrained enthalpy width (2026-09-05):**
> [Read these results first](buoyancy-constrained-enthalpy-width-results.md).
> New `addons/buoyancy_profile.py` represents C=Cc*G and H=Hc*G^(1/beta_H^2).
> Boundary beta_H is determined from five fluxes plus near-field F_b, not
> experimental temperatures. The fixed seven-case screen passes 7/7:
> beta_H=1.03284–1.04818, maximum temperature mismatch 0.5811 K and H2
> halfwidth mismatch 2.5826%. GL force convergence is checked separately;
> independent vectorized adaptive GK21 also passes 7/7 (max change 2.760e-6).
> Artifacts: `buoyancy_constrained_enthalpy_width_interface_2026-09-05.json`
> and `buoyancy_constrained_enthalpy_width_independent_audit_2026-09-05.json`.
> The targets are the former near-field moments; upstream grid convergence
> and measured-force accuracy are NOT revalidated by this boundary fit.
> This class is deliberately boundary-only. Do not call it a new completed
> downstream model, reuse the former 5% temperature improvement claim, or
> impose dF_b/ds=0 / d beta_H/ds=0 without a new justified closure.
> Next implement the missing thermal-shape transport law following
> [this requirements/derivation note](thermal-width-transport-requirements.md),
> then short conservative/refinement tests and only then full field scoring.
> Non-slow regression: **225 passed, 128 skipped, 58 deselected**.
> The bundled runtime lost manual scientific packages during an update.
> Use `.venv/Scripts/python.exe`; scientific/test versions are pinned in
> `requirements-research.txt` and all seven old boundary values replayed.
> No experimental defaults changed. Earlier entries below are historical.

> **Latest work — Gaussian volumetric enthalpy (2026-09-05):** read
> [the full results](gaussian-enthalpy-profile-results.md) and
> [the numerical correction note](enthalpy-profile-quadrature-note.md).
> The candidate is implemented in `addons/enthalpy_profile.py` and selected
> explicitly with `downstream_thermodynamic_profile="enthalpy"`; it requires
> consistent ambient phase thermodynamics and retains velocity lambda=1.16.
> The local C,H phase inverse, non-Gaussian-density three-unknown flux inverse,
> true-density buoyancy integration and receptor profiles are implemented.
> The original energy-only tensor convergence checkpoint is superseded: use
> `gaussian_enthalpy_profile_allflux_interface_2026-09-05.json` (7/7 passed).
> The completed all-flux field result retains all seven downstream arrays.
> MG=0.838826, VG=1.190780, FAC2=0.921053, width ratio=0.953289 and centre
> MAE=0.065565 m. This worsens height against the prior consistent-ambient
> control (0.061588 m), and the joint gate fails. No default changed.
> Both models now have the same 41 temperature sensors: minimum MAE improves
> 17.3360 -> 16.4738 K, median MAE worsens 21.5963 -> 21.8699 K.
> All 63 sampled sections pass the five-flux 1x/4x quadrature gate (2.985e-6).
> A separate buoyancy diagnostic reaches 1.029e-5 N/m; do not call it a pass
> of the same 1e-5 force threshold. Full 0.01 m step refinement was not run.
> Non-slow regression: **208 passed, 128 skipped, 58 deselected**.
> Next target a thermal-shape/transport state constrained by the transferred
> buoyancy moment as well as heat/species balances. Do not infer that TKE is
> uniquely required, and do not repeat coefficient fitting or the rejected
> equal-width profile as if it were a validated law. No processes remain from
> this experiment. Earlier entries below are historical.

> **Latest work — consistent ambient phase endpoint (2026-09-05):** read
> [the ambient/profile results](phase-ambient-consistency-results.md) first.
> An opt-in research correction now uses one explicit-species ideal volume
> closure for both ambient gas and the local phase inversion. The former
> no-H2 endpoint was about 1.20 K cold. The corrected endpoint and the dilute
> limit pass independent tests; no empirical temperature offset was fitted.
> All 7 interfaces and 38 concentration arcs / 17 vertical profiles completed.
> MG=0.821070, VG=1.198787, FAC2=0.921053; width ratio=0.952346 and centre
> MAE=0.061588 m. Improvements in concentration and height are small; the
> frozen joint promotion gate still fails. Existing defaults remain unchanged.
> All seven candidate trajectories are saved. A replay audit reconstructs
> sensor temperatures and checks stored temperatures, five fluxes and arc
> concentrations before comparing observations. Old open-pipe trajectories
> were not saved, so paired before/after temperature claims must not include them.
> Trial 23's 22 paired sensors give minimum-temperature MAE 16.6861 ->
> 16.3870 K; 33--42 K warm residuals remain at 1.19--1.78 m on the fixed
> release-axis sensor height. The new two-trial temperature population is
> 41 points, not the older 42-point source experiment. Do not pool them.
> Current non-slow regression: **184 passed, 128 skipped, 58 deselected**.
> Next test the independent enthalpy/density profile closure with this
> consistent ambient option; no new thermal-shape or TKE closure exists yet.

> **Latest work — source EOS and energy ledger (2026-09-05):** read
> [the current results](preslhy-source-eos-ledger-results.md) before continuing.
> Two measured-HEM defects were corrected: missing incoming pipe kinetic
> energy and lost condensed-phase enthalpy at Gaussian establishment.
> Four independent upstream energy checks now close below 7.31e-16; 7/7
> existing interfaces pass. The ~6% former mismatch is a conservation error,
> NOT a measured concentration error reduction. The near-field mixture EOS
> and the downstream thermal-profile/TKE states have not been replaced.
> Older HEM/para field scores and regression counts below are historical.
> The source audit is quiescent; the field evaluation retains actual wind.
> New field output records per-arc predictions and (for restarted nozzle
> shards) full downstream states, enabling residual diagnosis without another
> expensive model rerun. Exact quadrature symmetry removes duplicate phase
> evaluations without changing the integration rule or tolerance.
> Current non-slow suite: **162 passed, 128 skipped, 58 deselected**.
> The complete 38-arc/17-profile field rerun is finished: MG 0.820202,
> VG 1.199058, FAC2 0.921053, width ratio 0.956049, centre MAE 0.063587 m.
> Centre MAE improves from the former HEM 0.070042 m but remains worse than
> the baseline 0.054195 m; do not promote. Next use the stored nozzle
> trajectories to localise the residual thermal/density/profile errors.
> These scores describe the earlier version without the ambient consistency option.

> **2026-09-05 coupled-field update:** the conservative near-field/JETPLU
> handoff has now been independently exercised on seven momentum-dominated
> PRESLHY horizontal releases. Adaptive energy quadrature removes a false
> numerical rejection. Retaining constant source-momentum entrainment after
> the handoff underpredicts `sigma_z` by 29%; resuming the existing local
> density-scaled shear closure gives a mean model/measured ratio of 0.994,
> centre-height MAE 0.049 m and FAC2 0.952. This downstream mechanism is
> accepted for research. The complete model is not promoted: `scalar_peak`
> does not enforce the beta-derived total-mass target at Gaussian formation,
> while the four-flux boundary passes only 5/7 JETPLU thermal/profile
> interfaces even when searched through 20D. This established the need for an
> independently transported thermal/energy profile state, not another width
> or handoff-distance adjustment. See
> [`preslhy-coupled-crosswind-results.md`](preslhy-coupled-crosswind-results.md).
> Its seven-state equations, source terms, numerical gates and unchanged field
> decision criteria are now frozen in
> [`prereg-crosswind-independent-energy-state.md`](prereg-crosswind-independent-energy-state.md)
> before implementation or scoring. The first implementation gate is now
> complete: its five-flux boundary projection passes all seven fixed 10D
> interfaces, including former failures 10 and 25, with maximum flux residual
> `6.45e-16`, width mismatch 4.761% and temperature mismatch 1.990 K. This
> confirms the missing-state diagnosis; the downstream energy ODE and field
> score have now also been completed. Direct conserved-flux integration keeps
> the worst long-range balance residual below `1.04e-7`; the 0.02/0.01 m
> results are converged. VG, FAC2 and width improve to 1.176, 0.976 and 1.066,
> but MG moves to 1.118 and centre-height MAE worsens to 0.078 m. The frozen
> all-metric decision is therefore **do not promote**. See
> [`preslhy-independent-energy-interface-results.md`](preslhy-independent-energy-interface-results.md).
> Two coefficient-free follow-ups are also complete. Re-fitting the modelled
> direct-plus-ground-image profiles at the exact PRESLHY sensor heights passes
> its fit-quality gate but makes the candidate geometry rejection stronger:
> width ratio/centre MAE are 1.420/0.137 m against baseline 1.360/0.070 m.
> Li et al. (2026)'s published enthalpy-only established-flow equation then
> passes all seven interfaces and conserved marches, but changes internal
> centre MAE by only 0.000037 m. The extra rise is therefore neither a hidden
> observation-operator artefact nor mean-kinetic-energy thermalisation. See
> [`prereg-ground-image-geometry-observation.md`](prereg-ground-image-geometry-observation.md)
> and
> [`prereg-established-flow-energy-partition.md`](prereg-established-flow-energy-partition.md).
> Direct integration of section buoyancy then agrees with the projected
> finite-Gaussian force in sign and to within 6.761% for all seven interfaces,
> rejecting an interface-force jump as well; see
> [`prereg-interface-buoyancy-moment.md`](prereg-interface-buoyancy-moment.md).
> The following source-profile audit corrects Houf `B` from the scalar to the
> velocity e-folding width by dividing by `lambda=1.16`. The repeated field
> run improves MG/VG to 1.062/1.163 and internal width ratio to 0.976, but
> centre MAE remains 0.0758 m against the 0.0542 m baseline. The correction is
> retained as the research default without production promotion; see
> [`prereg-houf-width-mapping-correction.md`](prereg-houf-width-mapping-correction.md).
> A trial-10 force budget then finds only 0.0657 N vertical-momentum difference
> at 10D but cumulative buoyancy of 5.695 versus 2.734 N by 6 m. The candidate
> is already 64.3 K/weakly buoyant at 10D while baseline is 30.7 K/dense; see
> [`prereg-trial10-vertical-momentum-budget.md`](prereg-trial10-vertical-momentum-budget.md).
> A coefficient-free metastable dry-air bound then leaves the trial-10
> near-field centre at 65.559 K and weakly buoyant. It fails the unchanged
> interface temperature gate at 2.222 K, so the frozen protocol stops before
> a seven-trial field score. Suppressing all N2/O2/Ar condensation changes the
> near-field temperature by only about 0.42 K; equilibrium condensation heat
> release is not the source of the 10D contrast. See
> [`prereg-delayed-dry-air-condensation-bound.md`](prereg-delayed-dry-air-condensation-bound.md).
> The next physical audit should isolate pre-10D mass entrainment, not retune
> width, gravity or a nucleation-delay length.
> That audit is now complete. Li's Zone-V momentum-entrainment source had been
> integrated a second time across the Zone-IV 6.2D profile transition even
> though the input was already an air-loaded Station-3 state. A four-flux
> `source_flux` lower bound removes this unsupported addition, passes 7/7 and
> improves MG/VG/FAC2 to 1.012/1.162/0.976. It reduces trial-10 10D mass from
> 1.266 to 0.748 kg/s and centre temperature from 64.3 to 47.0 K, but by 6 m
> local entrainment restores most of the warming; centre MAE is still 0.0707 m
> versus baseline 0.0542 m. Keep it research-only and do not fit `beta_A` in
> Zone IV. See
> [`prereg-source-flux-gaussian-establishment.md`](prereg-source-flux-gaussian-establishment.md).
> The final regression check passes 114 non-slow tests plus the representative
> phase-profile crosswind handoff; existing defaults and gates are unchanged.
> Ground-contact geometry and the existing DEGADIS surface-layer entrainment
> were then transferred to the independent-energy path as explicit bounds.
> Trial 10 remains 0.469 m too high at 6 m, and the 1.78 m section is unchanged
> because contact starts at 2.60 m, so the frozen all-section gate stops the
> seven-trial run. A second attempt to use Li equation 35 from Station 4 fails
> the four-flux Gaussian boundary. The remaining defensible structural target
> is an added turbulent-energy/profile or finite-relaxation phase-slip state,
> neither of which can be calibrated from the current PRESLHY point data.

> **2026-09-03 update:** the thermodynamically inconsistent flashing-source
> lookup, non-conserved initial-entrainment momentum, and the resulting
> forced-ground workaround have been corrected.
> Use [`lh2-model-improvements-2026-09-03.md`](lh2-model-improvements-2026-09-03.md)
> for current LH2 performance numbers. Older figures below are retained as
> investigation history where explicitly marked.
> A later same-day audit also found that the fitted PRESLHY e-folding width
> had been mislabeled as Gaussian sigma. Any older 0.64/0.731 or "1.4 times
> narrow" statement below is superseded by the common-definition ratio 1.033;
> see [`gaussian-width-convention.md`](gaussian-width-convention.md).
> The condensed-air follow-up is also complete as an off-default research
> path.  Stable solid N2/O2 phase balance, particle settling and a conserved
> source march were tested at pre-registered 1/10/100 um sizes; all three
> worsened concentration scatter and plume-centre height and were not adopted.
> See
> [`prereg-condensed-air-particle-transport.md`](prereg-condensed-air-particle-transport.md).
> A 2026-09-04 total-energy audit then removed kinetic-energy creation in that
> research source and tested a HyRAM+ HEM/Yuceil--Otugen pressure-thrust
> option.  Both are conserved but rejected as default improvements; three of
> nine PRESLHY sources are incompatible with nominal-area pressure thrust.
> See
> [`prereg-source-total-energy-pressure-thrust.md`](prereg-source-total-energy-pressure-thrust.md).
> Li equation 22 then identified a missing 12--78 mm evaporation-zone source
> offset.  The explicit correction is state-neutral and slightly lowers the
> plume, but fails the pre-registered campaign criteria and remains off by
> default; see
> [`prereg-evaporation-zone-origin.md`](prereg-evaporation-zone-origin.md).
> The Li stationary-condensate momentum limit has also been implemented as S0.
> It conserves separate gas/solid kinetic energy but worsens PRESLHY and
> Spadeadam; both slip limits now over-rise, so do not fit an intermediate slip
> scalar.  See
> [`prereg-stationary-condensate-bound.md`](prereg-stationary-condensate-bound.md).
> The follow-on finite-rate screen is documented in
> [`finite-rate-sublimation-audit.md`](finite-rate-sublimation-audit.md).  A
> Ranz--Marshall minimum lifetime is implemented and tested, but a coupled
> kinetic source is intentionally not selected because neither campaign
> constrains condensed-air particle number, size or new-particle nucleation.
> A later Raman validation added a separate conserved-energy axisymmetric
> free-jet model. It passes all four Hecht--Panda slopes and reproduces an
> official HyRAM 6.1 oracle within 0.75% on aggregate metrics. A later flux
> audit found 8.9--14.8% species loss in the reproduced establishment
> boundary; two exact conservative alternatives close below `5e-12` but fail
> the thermal-centreline criterion. Treat it as a benchmark/research path,
> not a production default. See
> [`prereg-hecht-panda-hyram-closure.md`](prereg-hecht-panda-hyram-closure.md).
> The Raman follow-up now corrects unequal per-case stitched-image coverage
> (369 samples), and adds a conservative ambient-water frost phase. The
> 100%-RH upper bound is numerically conservative but thermally too strong;
> an unscored 40%-RH run passes all four slopes. Since Hecht--Panda mention
> condensed moisture but publish no laboratory RH, do not present 40% as a
> calibrated validation. Request measured RH or original environmental logs.
> See
> [`prereg-humid-air-frost-upper-bound.md`](prereg-humid-air-frost-upper-bound.md)
> and
> [`prereg-raman-case-coverage-correction.md`](prereg-raman-case-coverage-correction.md).
> The exact fields, identifiers, author contact and a ready-to-send request are
> in [`data-request-hecht-panda-humidity.md`](data-request-hecht-panda-humidity.md).
> The reported 0.3 m/s honeycomb co-flow has also been implemented with
> conservative ambient momentum/energy influx. It changes every Raman slope
> by less than 0.8% and does not resolve the thermal error; see
> [`prereg-hecht-panda-coflow.md`](prereg-hecht-panda-coflow.md).
> **Superseding Raman result:** component H2/N2/O2/H2O ideal-gas enthalpies
> have now been combined with condensed-air equilibrium. The dry candidate
> conserves all boundary invariants and passes 4/4 slopes under both 549- and
> 369-point protocols. Use it as the recommended conserved axisymmetric
> research path; do not call it humidity-conditioned validation. See
> [`prereg-phase-temperature-dependent-enthalpy.md`](prereg-phase-temperature-dependent-enthalpy.md).
> The recommended flags are now wired to the public research functions
> `lh2_source_from_measured_throat` and `run_lh2_near_field_research`; legacy
> and atmospheric defaults remain unchanged. A pre-registered argon phase
> correction remains off by default because it worsens the centreline thermal
> metric despite passing 4/4 and conservation. A perfect-black radiation ODE
> bound also remains off: it improves the thermal centre by only 0.83
> percentage point while slightly degrading both mass metrics. See
> [`prereg-argon-phase-completeness.md`](prereg-argon-phase-completeness.md)
> and [`prereg-radiation-upper-bound.md`](prereg-radiation-upper-bound.md).
> A separate para-hydrogen caloric audit is complete. It improves the two
> mass slopes but makes the corrected thermal-centre error -34.79% instead of
> -21.38%; it is a conservative, rejected sensitivity, not an omitted heat
> source. See
> [`prereg-hydrogen-spin-isomer-enthalpy.md`](prereg-hydrogen-spin-isomer-enthalpy.md).
> The follow-on four-flux/two-scalar phase model is conservative and fixes the
> radial width ordering, but its corrected thermal-centre error is -29.14%
> and its temperature radial coefficient is outside the reported range. Keep
> it as a research closure only; see
> [`prereg-phase-two-scalar-four-flux.md`](prereg-phase-two-scalar-four-flux.md).
> ELVHYS Tests 10/11 have also been reduced from checksum-verified 20 Hz raw
> records. Scoring stopped before a model run because H2 mass flow is absent
> and the two official geometry documents disagree on nozzle elevation. The
> result and exact author request are in
> [`elvhys-tcs-audit.md`](elvhys-tcs-audit.md) and
> [`data-request-elvhys-source.md`](data-request-elvhys-source.md).
> **Latest Raman provenance correction:** the active observations are now the
> final-journal fits, not the 2017 conference values. The two mass slopes are
> 0.2771 and 0.07069; the recommended dry model errors become -24.86%, -15.79%,
> -21.38% and +24.56%, still 4/4 but narrowly. Figures 6 and 8 also list a
> tenth `4 bar, 45 K, 1.25 mm` series absent from Table 1 and the radial panels.
> Treat every aggregate pass count as provisional until the authors confirm
> fit membership and uncertainty. See
> [`hecht-panda-journal-benchmark-correction.md`](hecht-panda-journal-benchmark-correction.md).
> **Conservative crosswind handoff now implemented:**
> `run_lh2_crosswind_research` calculates one consistent local wind, runs the
> conserved near field, transfers its phase-aware radial thermodynamics, and
> projects mass/species/vector momentum before auditing energy, H2 width and
> centre temperature. The centre-state and flux-energy candidates are retained
> as rejected controls. The phase-profile candidate passes the frozen 0.08 m
> screen and runs JETPLU downstream; it is limited to horizontal wind-aligned
> releases. Its later PRESLHY result and the accepted local-shear mechanism
> are summarized in the 2026-09-05 update at the top of this file.
> See
> [`conservative-nearfield-crosswind-handoff-results.md`](conservative-nearfield-crosswind-handoff-results.md).
> The immediate Spadeadam pilot is complete but rejected before concentration
> scoring. The phase-manifold extension passes test 6; test 4 fails the frozen
> width/temperature interface and has source velocity/local wind = 9.091,
> below the fixed ratio 10. It requires a crosswind-aware two-phase near-field
> source or measured single-phase outlet, not relaxed tolerances. See
> [`prereg-coupled-crosswind-spadeadam.md`](prereg-coupled-crosswind-spadeadam.md).

Written at the end of a long working session. It covers what the code is, what
has been established and to what standard, every dataset and what it can and
cannot answer, the mistakes made along the way and why they matter, and what to
do next.

The process matters as much as the results here. Five times a "model defect"
turned out to be an error in how the comparison was set up, and eight
hypotheses were formed, tested and falsified. Those are recorded in full,
because the main way to waste the next session is to re-run them.

---

## 1. What this is

`degali` is two things in one package, chosen by name:

```python
from degali.presets import DEGADIS_21, LIQUID_HYDROGEN
```

**`DEGADIS_21`** is a port of the 1989 EPA model, reproducing the Fortran to
1e-12 across all six programs and all five EPA test cases. Nothing in it has
been changed. It is the reference, and every option added for hydrogen
defaults to off so that this stays true.

**`LIQUID_HYDROGEN`** is the same code with equations of state instead of the
1989 correlations, a buoyant closure so a cloud that becomes lighter than air
can leave the ground, a flashing source, a horizontal jet path, and three
corrections to that path.

**The suite passes.** Without any external data the field tests skip and the
rest pass; with REDIPHEM, PRESLHY and SMEDIS configured, all of it runs. The
count is deliberately not written down — it was, in four documents, and had
drifted in all four.

### Layout

```
src/degali/
  core/          the ported model: rkgst, steady, transient, jetplume,
                 blanket, thermo, atmosphere, numerics
  addons/        liftoff, unified closure, notional nozzle
  validation/    rediphem, smedis, preslhy, flashing, compare, statistics
  presets.py     DEGADIS_21 and LIQUID_HYDROGEN
  lh2.py         the one-call entry point, with scope warnings
reference/fortran/   the original, vendored, with a build script
tests/               the suite
docs/                21 documents; the map is in README.md
```

---

## 2. Where the model stands

### Established, deterministic (grade A)

- Reproduces DEGADIS 2.1 to 1e-12, all six programs, all five EPA cases.
- Three findings in the original, with line numbers: `ADIABAT` does not assign
  `wa` on the `ifl=1` path while assigning it on the other four, and `SZF.for`
  line 95 passes an undeclared `walay`; `PSS.for` 89 and `SSG.for` 92 write the
  layer temperature to different variables with every other argument matching;
  `GAMINC` returns the unregularised incomplete gamma **on purpose**, with a
  comment saying so — a porting trap rather than a defect.
- Liquid hydrogen is buoyant across the whole flammable range. On the adiabatic
  mixing line the density ratio is 0.98 at the lower flammable limit, 0.87 at
  stoichiometric and 0.69 at the *upper* limit of 75 mol %, with the strongest
  buoyancy at 76 mol %. This is a calculation, not a fit, and it is why the
  standard dense-gas models do not apply.

  **The crossover itself should not be quoted as a round number.** An earlier
  version of this document said 85–90 mol % — superseded; the mixing table
  puts it at 99.97. The discrepancy is not a rounding disagreement — it is that the top
  of the table is one linear segment from 0.99894 to 1.0 carrying a 50 %
  density change, in a regime where air is treated as a non-condensing ideal
  gas (`CoolPropBackend.cp_air` clamps at 100 K and says so). Quote the
  flammable limits, which are nowhere near that segment. The headline claim
  does not depend on the crossover.
- The equilibrium flammable mixture is buoyant, but the conserved cryogenic
  source can remain dense for several metres while warming. In the adopted
  Spadeadam test-6 calculation the density crossover is near 3.2 m; the older
  half-to-one-metre estimate used a different source boundary.
- EPA's published DEGADIS bias on Burro, FB −1.07, is an evaluation-height
  artefact: measurements at 1 m compared against ground-level predictions
  (their §4.3). This port gives −1.29 at ground and −0.29 at 1 m.

### Established, statistical (grade B)

**Superseded — recomputed from the reduction; the figures below were never computed by anything.** Recomputing gives MG 0.855, CI [0.672, 1.064], VG 1.30, FAC2 0.85 as shipped and MG 1.070 corrected, over 66 arcs, and the interval now includes 1. See `docs/lh2-recomputed.md`.
| claim | value | n |
|---|---|---|
| LH₂ near-field concentration, momentum-driven | MG 0.738, CI [0.591, 0.918], VG 1.41, FAC2 0.83 | 69 arcs, 9 trials |
| same bias at both release heights | 0.722 and 0.746 | 23 and 46 |
| wind-steered releases are a different population | GSD 2.94 against 1.16 | 13 and 7 trials |
| the measured LH₂ plume does not rise over 6 m | +0.00 m, 0.5 m releases | 42 vertical fits |
| the vertical spread deficit | 0.64 of measured, as shipped | same fits |
| Burro reads high at the lowest height | MG 0.811, CI [0.63, 1.02] | 61 |
| the vertical profile is too steep | MG 7.3 at 3 m, 5500 at 8 m | 59, 53 |
| the same defect on an independent dataset | MG 0.71 at 0.1 m to 0.11 at 8.5 m | 76 SMEDIS sensors |

Intervals resample **trials**, not points: readings within a trial share a
release, a wind and a source estimate. Point-level resampling gives
[0.68, 0.80] where trial-level gives [0.59, 0.92], and reporting the former
would overstate the evidence by roughly the square root of the sensor count.

### Directional only (grade C)

Lift-off height, RMS 2.9 m on four NASA spills, three of the four "measured"
values being inversions of thermocouple data through a mixing model. The
buoyancy *regime* is right 4 of 4, which is strong for what it is and still a
sample of four. Write these as "consistent with", never assert them.

---

## 3. Historical three-correction stage

> Superseded by the five-correction model. The two later source corrections
> preserve flash composition/enthalpy and conserve total momentum during
> initial air entrainment. Current numbers are in
> `lh2-model-improvements-2026-09-03.md`; this section remains as investigation
> history.

All three are published relations that the shipped model either omits or
applies at the wrong plane. Each is off by default.

**Expanded-source start.** The model was taking the density from after the
flash while keeping the concentration and the area at the orifice. That is the
error EPA's 1991 evaluation records the SLAB developer objecting to — "for jet
releases the source area should be the cross-section of the fully expanded jet
rather than the orifice". `JetPlume.expanded_source_start` builds a consistent
triple with `rho u**2 A` carried across.

**Density-scaled entrainment.** Ricou and Spalding measured entrainment
scaling as `sqrt(rho_ambient / rho_jet)`; Panda and Hecht restate it for
cryogenic hydrogen. `JETPLU` has no such scaling. **The sign is easy to
invert** — it was inverted here first, and moved the vertical spread the wrong
way, 0.70 → 0.65, before the literature was re-read.

**Plume entrainment coefficient.** `alfa1 = 0.0875` is the pure-plume value
measured by Papanicolaou and List (1988), against their pure-jet value 0.0545;
`JETPLU` ships 0.057. An LH₂ release is buoyant within a metre and spends its
measured range as a plume while being entrained as a jet. The widely cited
0.0833 is Fischer et al.'s earlier proposed value, now superseded here;
HyRAM's independent default of 0.082 remains useful corroboration.

### What they achieve

**Superseded.** The table below was never computed by anything. Recomputed
from the reduction, over the arcs at 3-6 m and the 66 arcs both
configurations reach:

| | as shipped | corrected | measured |
|---|---|---|---|
| vertical spread ratio, 0.5 m releases | 0.82 | **1.07** | 1.00 |
| vertical spread ratio, 1.5 m releases | 0.84 | **1.03** | 1.00 |
| rise, 0.5 m releases | 0.24 m | 0.23 m | −0.06 m |
| rise, 1.5 m releases | 1.08 m | 1.06 m | +0.03 m |
| concentration MG | 0.855 | **1.070** | 1.00 |

The spread correction holds and is if anything cleaner than reported. **The
trajectory result does not.** The published 1.07 → 0.19 m does not reproduce:
this configuration lifts the 0.5 m releases 0.24 m and the corrections barely
change it. The fault is real — the measurement says the plume does not rise —
and the corrections do not address it. `docs/lh2-recomputed.md`.

The superseded figures, kept because a withdrawn number should stay
visible next to the one that replaced it:

| superseded | as shipped | corrected | measured |
|---|---|---|---|
| vertical spread ratio | 0.64 | 0.94 | 1.00 |
| rise, 0.5 m releases | 1.07 m | 0.19 m | 0.00 m |
| rise, 1.5 m releases | — | 0.51 m | 0.03 m |
| concentration MG | 0.737 | 1.228 | 1.00 |

Concentration MG passing through 1.000 at the middle configuration is **two
errors meeting**, not agreement: the narrow section raises the centreline and
the excess rise lowers it. Fixing the spread unmasks the trajectory. That is
why the three quantities are judged separately.

### On "no tuning"

An earlier draft of this claimed zero tuning coefficients. That was
self-flattery. Every individual value is published, and no new free parameter
was introduced — but eight mechanisms were tried and the three that helped were
kept, and the selection was made against the data. It is model selection rather
than parameter fitting, and it should be described that way. EFFECTS uses a
1.52 correction factor on the dew-point concentration and says so; this work
should be equally plain.

---

## 4. What was tried and failed

Eight hypotheses, formed, tested, falsified. **Do not re-run these.**

| hypothesis | how it died |
|---|---|
| FLADIS channel numbering could be guessed | the guessed channels read 302 and 23.5 at 20 m; not concentrations |
| Coyote's bias is the pool radius | forcing smaller diameters makes it worse |
| LH₂ rainout feeds a ground-level source | the report: no rainout for un-impinged elevated releases |
| the notional nozzle fixes the near field | FAC2 0.69 → 0.48; against `sigma_z`, 0.64 → 0.57 |
| the section is too narrow by a factor of two | measured on mixed release heights; the real factor is 1.4 |
| the residual is an entrainment shortfall growing with distance | pre-registered; slope +0.018, CI [−0.168, +0.171] |
| the trajectory fault is the Boussinesq approximation | pre-registered; the correction *increases* the rise, 0.81 → 1.08 m |
| added mass suppresses rise for wide sources | tested elsewhere across six widths; no coefficient fits |

And four corrections that are right for their own regime and do not transfer:

| correction | source | effect here |
|---|---|---|
| additional pressure drag | Mack et al., EFFECTS | 0.02 m at their `C_d = 0.39` |
| vertical component in the shear velocity | Mack et al., EFFECTS | under 5 % |
| spread floor | Hart & Harper, UDM | floor binds in 180 of 26 646 calls |
| Richardson jet-to-plume transition | Pantokratoras | never reaches the threshold; `u_c` is too large |

The EFFECTS pair are aimed at plumes rising fast enough for a quadratic drag
and a vector shear velocity to matter — Witcofski's pool rises twenty metres.
These jets rise 0.2 to 0.5 m over six. **The papers do not say this**, so it is
worth recording.

---

## 5. The five comparison errors

Every one of these presented as a model defect and was a comparison defect.
This is the single most transferable lesson here.

| symptom | actual cause |
|---|---|
| FLADIS gave a clean MG 0.89 | the channels read were not concentrations |
| trial 2 mass flow was −0.5 g/s | the window stretched across a spike; use the longest sustained run |
| the 1.5 m releases failed completely (model 0.0 %, measured 43–98 %) | Table A3's `z` is measured from the release axis, not the ground |
| the far field was "unusable" | the plume position can be fitted out; it does not need to be known |
| Thorney Island had "no usable measurements" | it had eighteen arcs; the harness chose a pool route for a puff and then discarded the real failure reason |

**Suspect the comparison before the model.** All five were caught by looking
again rather than by accepting the result.

---

## 6. Data

### On disk

| path | size | what |
|---|---|---|
| `/home/claude/e35/10.35097-1481/` | 26 MB | PRESLHY E3.5 raw, 24 workbooks |
| `/home/claude/exp/rediphem/all/` | | REDIPHEM archive (`$REDIPHEM_ROOT`) |
| `/home/claude/exp/batch1,2,3/` | | SMEDIS spreadsheets, 30 trials |
| `/home/claude/slab/` | 328 KB | SLAB.FOR compiled and run |
| `/mnt/project/` | | reports, papers, reduced CSVs |

**None of it ships.** The readers take a path or an environment variable and
the tests skip without them.

### Readers written

| module | reads |
|---|---|
| `validation/preslhy.py` | the 24 workbooks; `fit_vertical`, `fit_arcs` |
| `validation/smedis.py` | 30 SMEDIS trials, 1250 sensors, 79 arcs |
| `validation/rediphem.py` | the REDIPHEM archive |

### What each dataset settled

| dataset | settled |
|---|---|
| DEGADIS Fortran + 5 EPA cases | the port, to 1e-12 |
| REDIPHEM Burro (8) | MG 0.811 at 1 m; the vertical profile fails above it |
| REDIPHEM Desert Tortoise (4) | MG 1.839; EPA's own score explained |
| SMEDIS FLADIS + DT (76 sensors) | the vertical defect on an independent substance |
| PRESLHY near field (66 arcs) | MG 0.855, CI [0.672, 1.064] — recomputed; superseded figures were MG 0.738, CI [0.591, 0.918] over 69 arcs |
| PRESLHY far field (18 arc fits) | the bias reverses; the trajectory is the cause |
| PRESLHY vertical (42 fits) | the plume does not rise; `sigma_z` 1.4× narrow |
| NASA Witcofski (4) | buoyancy regime, 4 of 4 |
| EPA-450/4-90-018 | an external anchor and the height artefact |

### What each refused, and why

| dataset | refused because |
|---|---|
| PRESLHY far field, point-to-point | wind direction at 22.5° against stands 10–12° apart — solved by fitting the position out instead |
| PRESLHY near field, transients | travel time under a second against 0.3 s sampling |
| PRESLHY Dräger sheet | 90 channels at 1 Hz, but values run −275 to 1970 under a `%` header; **units unresolved, and this is the one place transients could be tested** |
| REDIPHEM FLADIS | ships no `CHANDEF.DAT`; recovered from SMEDIS |
| REDIPHEM Eagle | N₂O₄ dissociates to NO₂; temperature-dependent molecular weight |
| REDIPHEM Thorney Island | the van Ulden momentum balance has no solution at H/D ≈ 1 — the original Fortran stops in the same place, and refining the grid a hundredfold does not move it |
| SMEDIS Thorney Island | peaks at 2060 under a `mean_C(%)` header |
| SMEDIS Prairie Grass | peaks at 235, same header |
| SMEDIS EMU | sensor heights read 244 553 m |
| BA-Propane (273 sensors) | tried; MG 0.168 against the REDIPHEM reduction's 3.76 — **the same trials, two reductions, a factor of twenty apart** |
| PRESLHY E3.4 pool | concentrations at 35/45/55 cm only, no downwind distance |

### Still untried

- **The Dräger time series.** Far field, travel time 3–7 s, resolvable against
  1 Hz. Blocked only on units. The one route to a transient test on hydrogen.
- **Witcofski Table 2.** Time-resolved grab bottles for Test 6 at tower 5,
  three heights. The OCR scrambles which value belongs where.
- **Chirivella and Witcofski (1986)**, AIChE Symp. Ser. **82**, 120–140. The
  only route to a lift-off sample beyond four.
- **Sandia cryogenic jet data** (1 mm, 80 K), which NCSRD used. It would span a
  range of density ratios; every campaign in hand is at one, 4.4, which is why
  the Boussinesq scaling could not be tested properly.

---

## 7. Model comparison

Four models run on one PRESLHY release: 25.4 mm, 260 g/s, 5 barg, 0.5 m
elevation, 2.5 m/s wind. Each as shipped, with its own defaults.

| | rise at 6 m | `sigma_z` at 2 m | distance to 4 mol % |
|---|---|---|---|
| SLAB 1990 | slumps to ground, run ends at 2 m | 0.043 m | — |
| DEGADIS 2.1 | +1.07 m | 0.191 m | 26.3 m |
| HyRAM 6.1 | +3.15 m | 0.278 m | 10.0 m |
| **this work** | **+0.19 m** | **0.320 m** | 24.2 m |
| measured | +0.00 m | 0.34 m | not reached by the array |

**Read the caveats before the numbers.**

- SLAB has no buoyant branch. `vertical_jet.py:145` in the Python
  reimplementation raises `BuoyantRiseNotImplemented`; the 1990 Fortran
  silently slumps the cloud instead.
- **HyRAM's `Jet` assumes a quiescent ambient.** It was given no wind, so it
  rises further than a wind-blown plume would. This is not a like-for-like
  comparison and is shown to make that visible.
- DEGADIS 2.1 and "this work" are the same code with options off and on. **Only
  this pair is a clean model comparison.**
- The array reaches 6 m. The LFL distances are extrapolations and their 2.6×
  spread is the honest statement of what an integral model settles here.

PHAST and EFFECTS could not be run. Both publish on Witcofski, qualitatively
("plume centreline z-positions agree well"), with no numbers to plot.

Figures: `lh2_five_models.png`, `lh2_centreline_decay.png`,
`lh2_model_comparison.png`.

---

## 8. What is open

**The trajectory, 0.19 m on the low releases and 0.51 m on the high ones.**
Five candidates registered and tested; what remains is slip between droplets
and the gas, which an integral model carrying one velocity per cross-section
cannot represent. `liquid_fraction` implements the mass the droplets add — a
5 % effect — and not their falling.

**Bulk-air condensation is now bounded, not solved.**  The former source
starts its gas table near 20 K with 1.07--1.55 kg dry air/kg H2.  The optional
`bulk_air_phase_safe_source` advances a conserved plug flow to the first
all-gas N2/O2/Ar state near 68 K.  It removes that invalid phase assumption,
but fails the pre-registered validation: PRESLHY VG and centre-height MAE
worsen.  It remains off by default.  The result shows that a warm single-phase
handoff removes dense condensate too early; the next implementation needs a
transported condensate mass, re-evaporation rate and slip momentum together.
See `prereg-bulk-air-phase-boundary.md`.

**Li et al. (2026) has now been obtained and reproduced.**  The exact
Zone-III calculation is in `addons.cryogenic_air`, with a phase-domain guard.
It cannot be adopted for PRESLHY: its 51--55 K validation already evaluates
liquid N2 below the 63.15 K triple point, its latent term cancels
algebraically, and it omits O2 phase change.  NBS solid N2/O2 vapour-pressure
fits have been added as the first prerequisite for a transported solid-air
state.  See `li2026-air-condensation-audit.md`.

**The vertical spread, +3.3 % mean bias.** No width correction is warranted.
The earlier apparent deficit came from comparing an e-folding width with a
standard deviation.

**Why `alfa1 = 0.0875` is used.** Papanicolaou and List (1988), §3.1,
pp. 353–354, derive 0.0545 for a jet, 0.0875 for a plume and a pure-plume
Richardson number of 0.716 from their measurements. The 0.0533, 0.0833 and
0.557 values they quote for comparison belong to Fischer et al. (1979); the
earlier attribution through two secondary sources is superseded. The primary
source therefore closes the load-bearing coefficient issue. It also supports
the direction of stronger buoyancy-generated transport, though it does not
supply the unresolved M4 coefficient.

### Targets, measured rather than argued

Any change should be judged against the quantity it was meant to affect:

- trajectory: **0.00 m** at 3–6 m for 0.5 m releases, from 42 fits
- vertical spread: `sigma_z` about **0.64 m** at 5–7 m, after converting the
  former 0.91 m e-folding width to a standard deviation
- lateral spread: the fitted standard deviation includes wind meander; the
  former 17° estimate becomes roughly 12° under the corrected convention

`fit_vertical` and `fit_arcs` reproduce all of these. **Judging a change on
concentration alone is how four of the five comparison errors happened.**

---

## 9. Next session

### 2026-09-05 current-state addendum

The older queue below is retained for history, but two items have now been
closed. A measured TC3/PT2/pressure-loss source, atmospheric flash and
micrometre LH2 droplet scale have been implemented. The homogeneous-equilibrium
fast bound passes 7/7 interfaces but is rejected on the complete 38-arc,
17-fit field score. GASFLOW equation 45 shows that retaining a mass-transfer
delay would require 41--76 um droplets, 38--53 times the correlated diameter;
do not add a free Lee coefficient or treat LH2 droplet slip as the leading
campaign-scale correction without a measured large-droplet tail.

A spin-consistent pure-para source/near-field end member is also complete. It
reduces fast-bound centre-height MAE from 0.07004 to 0.05453 m but worsens VG,
FAC2 and width, so it is an explicit uncertainty rather than the default.
The released ortho/para composition is now part of the exact data request.
The current non-slow regression result is **136 passed, 128 skipped, 58
deselected**.

The highest-value remaining structural target is the downstream
thermal/density profile or turbulent-energy state, with a smaller
coefficient-free audit still available for H2 real-gas volume in the coldest
near-field cells. Two-velocity condensed-air transport is lower priority
because both carried and stationary momentum bounds already over-rise and no
particle-size distribution exists to identify an intermediate closure.

1. **Read this document and `claim-grading.md` first.** The grades stop the
   n=4 lift-off result from being written as firmly as the n=69 concentration
   one.
2. **Run the non-slow suite** with `PYTHONPATH=src`; the current result is 162
   passed, 128 skipped and 58 deselected.
3. **Pre-register before testing.** Five registrations are in `docs/`; the
   entrainment and Boussinesq ones were both falsified, and both would have
   been written up as successes without the criteria fixed in advance.
4. Highest-value remaining work, in order:
   The exact trials, channels, units and acceptable file formats are listed in
   [`data-request-lh2-remaining-state.md`](data-request-lh2-remaining-state.md).
   - isolate why the conservative independent-energy state increases buoyant
     centre-height error without retuning the now-correct width closure
   - obtain the ELVHYS Test-10/Test-11 hydrogen mass-flow boundary and resolve
     the 200/250 mm nozzle-coordinate conflict; only then run the frozen
     Bottom-1/Bottom-3 comparison
   - obtain measured effective discharge area/feed-line pressure and a
     condensed-air particle-size or relaxation-time range for the PRESLHY
     source hardware
   - replace common axial gas/particle velocity with a finite-relaxation slip
     closure, then add persistent water/CO2/argon condensate inventories
   - the Witcofski comparison against EFFECTS on **their** scenario (note it
     runs the pool path, so the three corrections do not apply)
   - resolve the Dräger units and test transients
   - Sandia cryogenic jet data to span density ratios

### On the papers

**Paper A, evaluation methodology, is ready.** Two grade-A results with an
external anchor and reproducible code: EPA's published DEGADIS score is an
evaluation-height artefact, and conventional arc-maximum evaluation hides a
structural defect that shows up independently on SMEDIS ammonia jets.
*Process Safety and Environmental Protection* is the fit.

**Paper B, LH₂, needs the right framing.** Not "we improved the model" but
"local validation separates faults that a concentration comparison mixes, and
two of them turn out to be published relations applied wrongly". The eight
falsified hypotheses are evidence the method works, not a weakness. Mack et
al. (PSEP 176, 2023) is the natural predecessor.

**Do not combine them.** A's strong results would be dragged down by B's open
residual.

---

## 10. Standing instructions

From the person running this work, and worth keeping:

- **Do not conclude "impossible" or "cannot be done".** Twice a "structural
  limit" was declared here and twice it was wrong — the far field was
  recovered by fitting the plume position out, and the trajectory error turned
  out to be 44 % initial condition.
- **When a validation fails, suspect the comparison first.** Five for five.
- **Aggregate statistics on a model with more than one fault produce plausible
  wrong diagnoses.** Measure the sub-models.
- **Attribution needs the source and a line number**, not a plausible reading.
  A parallel SLAB reimplementation reported six defects in the original and,
  on obtaining the Fortran, found one.
