# degali — project status

> **USER RESUME 2026-09-06 14:34 KST:** measured-vector `step001` executed and
> completed (`reference/preslhy/observed_wind_vector_step001_2026-09-06`,
> 1066 steps). This is a measured-vector + observed magnitude + observed
> direction case with an explicit 0.01/0.02 cross-check.
>
> Same 19T/5C/3vertical/same-key summary:
>
> - **CONTROL:** 27.0312 K (mean MAE, baseline for this contrast).
> - **yaw .02:** 19.0453 K (−29.54%), MG 0.932706 (improved), width ratio 1.235037
>   (worse), centre MAE 0.177467.
> - **magnitude .02:** 29.2721 K (+8.29%), MG 0.774623 (∼no improvement),
>   width ratio 1.093047 (worse), centre MAE 0.188774.
> - **vector .02:** 23.7324 K (−12.20%), MG 0.871162 (improved), width ratio 1.292207
>   (worse), centre MAE 0.200676.
> - **vector .01:** 23.7332 K (−12.20%), MG 0.871186 (practically equal to .02),
>   width ratio 1.292243, centre MAE 0.200680.
>
> `measured_vector` `.01` is almost unchanged from `.02`, so step refinement does
> not recover width-trajectory trade-off. No new candidate promotion.

> **USER PAUSE 2026-09-06 14:23 KST:** stop new work and automatic continuation.
>Yaw trial10 .02 COMPLETE; interrupted .01 resumed from960 to1066 COMPLETE in
>reference/preslhy/yawed_trial10_step001_resumed_2026-09-06. Source/old134 intact.
>22yaw+3resume+8observed-wind tests PASS separately. Yaw-only mean19T MAE
>27.03115→19.04639K, but Tmin18.43483→23.36689K and vertical width ratio
>1.02950→1.23508 worse. NO adoption. .02/.01 maxsensor difference.0158211K.
>New magnitude-only and measured-vector .02 calculations have saved terminal
>artifacts under observed_wind_magnitude_step002_2026-09-06 and
>observed_wind_vector_step002_2026-09-06; read complete/failure FIRST on resume.
>Next: audit both new wind cells, measured-vector .01 only if accepted, then
>strict same-key four-cell comparison and updated delivery/regression. Do not
>restart completed runs. Latest automation id degali is being removed for pause.

> **Stage2 13:40 latest:** full non-slow regression36296 ENDED13:11:05,
>817PASS/128SKIP/58deselected,246Python snapshot unchanged. Lateral/raw-wind
>audit and conservative ambient-vector audit COMPLETE;8+8 tests separately PASS.
>No active runs. Trial10 verified wind implies missing lateral entrainment
>momentum14.67% at4m, not a new trajectory or accuracy gain. Trial23 far1Hz
>explicitly faulty and excluded. Read [results](lateral-wind-vector-results.md)
>FIRST. Next bounded opt-in3D yaw candidate, same-source conservation and
>coplanar limit before actual trial10 comparison. Deadline15:50KST unchanged.

> **Stage2 13:01 latest:** molecular reference screen COMPLETE13059,
>7fields/3247nodes/16235fixed-q path samples,16new tests PASS;134unchanged.
>Reference Fick width increase<=0.20402%, not a rigorous multiphase/T-error
>bound. Soret unit response34/315gas states;281unsupported, no alpha selected.
>Full non-slow regression RUNNING36296,246Python snapshot fixed under
>reference/preslhy/stage2_regression_2026-09-06_1300. No code changes until done.
>Read [scale/source/resume](molecular-transport-scale-results.md) FIRST.
>Next: matched left/right temperature symmetry and wind/coordinate evidence,
>without fitted yaw; then verified regression/delivery matrix.

> **Stage2 12:46 latest:** ensemble/refinement sessions24970/33554 ENDED.
>82point/support rows complete;656final objectives:623solver-certified,
>581grid-resolved,42certified but grid-unresolved,33failed;155intermediate
>failures preserved. Successful mean-T cooling<=1.04e-6K, not a new field.
>Same41raw TC sample means added read-only: CONTROL mean-MAE18.56486K,
>old minimum17.33598K unchanged; centre1.19m mean gaps10/23=12.30455/23.61503K.
>21ensemble+8mean tests PASS separately;134frozen unchanged. No active runs.
>Read [results/resume](thermochemical-ensemble-results.md) FIRST; no PDF fitting.
>Next independent screen: molecular heat/species cross-transport scale, not
>repeating failed LP refinements or completed phase/meanKE/TKE audits.

> **Stage2 12:10 latest:** common nonideal N2/O2 liquid excessG/h flash
>implemented,40new/265related tests PASS;134frozen/runtime unchanged. Same9
>NBS correlation pressure mean error6.0318%→1.76193%, not independent validation.
>16/356in65--77.5K local samples: maxfixed-h ΔT.385334K;10below/330above
>retained. No new field/adoption. [Result/resume](liquid-excess-potential-results.md).
>All calculations ended. Next D: conservative unresolved thermochemical mixing/
>mean-EOS applicability, no arbitrary PDF variance or repeated syntheticTKE.

> **Stage2 11:50 latest:** full non-slow regression ENDED59000:726PASS,
>128SKIP,58deselected; snapshot229unchanged at completion11:39:58. Not measured
>code coverage or Linux/slow completion. Actual enthalpy tangents315samples:
>314initial PASS +1near-table-knot resolved with smaller differences at original
>tolerance; original failure retained.6new tests AFTER full run PASS. All runs
>ended;134frozen/runtime unchanged. Same-D actual transport unchanged, no field
>improvement/promotion. [Tangent result and regression](actual-enthalpy-tangent-results.md).

> **Stage2 11:27 latest:** actual mechanical/thermal scale audit COMPLETE on4
> saved fields. Pressure-loss10/23 accounted posthandoff mechanical budgets
> <=.04792%/.06684% of initial thermal deficit; NOT upstreamTKE universal bounds.
> New reference-covariant species-carried enthalpy operator+29tests PASS;
> unequal diffusivities require the composition enthalpy cross-flux. Ratio1
> recovers original exactly; no field correction/adoption yet.27energy tests,
>253related tests PASS. Full non-slow regression RUNNING59000 at
> reference/preslhy/stage2_regression_2026-09-06_1130; do not edit snapshot code.
> [Energy result](mechanical-thermal-scale-results.md) · [New operator/resume](species-carried-enthalpy-results.md).

> **Stage2 11:02 latest:** bounded-liquid trial10 field COMPLETE,533steps/11.417m,
> same5C/3geometry/19T, balance1.026e-7. NO improvement/promotion: Tmin MAE
>18.434825→18.436961K; max sensor change.034344K. Trial23unsupported, fine-step
> and whole7 not done. All Python runs ended. New41sensor fixed-section cold
> envelope diagnosis COMPLETE: some near-core residuals cannot be explained
> by sensor placement alone; trial23/1.19m median gap29.942K.21new/197related
> tests pass,134frozen unchanged. Next actual mechanical/thermal energy-scale
> diagnosis, not repeated syntheticTKE or completed phase screens.
> [Field result](bounded-liquid-handoff-results.md) · [D result/next](downstream-cold-envelope-results.md)
> · [Primary mixed-phase evidence](cold-mixed-phase-evidence.md). Deadline15:50KST.

> **Stage2 10:27 latest:** bounded-liquid10boundary PASSED original5flux/width/T
> gates;23unsupported because target60.0047K cannot meet2K with64K guard.
>22new/176related tests PASS, frozen134 unchanged. Actual trial10 field RUNNING
> session19255, checkpoint directory reference/preslhy/bounded_liquid_trial10_2026-09-06.
> Read complete.json/failure.json and checkpoints before any new run. Do not
> edit hash-sealed candidate code during run. No dispersion improvement claim.
> [Results/resume](bounded-liquid-handoff-results.md).

> **Stage2 10:05 latest:** common pure/ideal-mixed liquid G module and local audit
> COMPLETE.38new/154related tests PASS;134frozen/runtime and hashes unchanged.
> Binary pressure benchmark mean/max6.0318%/10.5666%, NOT improved versus prior
> ideal screen. Caloric-only local h correction<=133.51J/kg, not field error.
> No process active. Next same-source10/23 in-domain conservative handoff
> feasibility; do not wait indefinitely for complete mixed-solid EOS or splice
> unsupported lowT. [Results](liquid-phase-potential-results.md). A remains partial.

> **Stage2 09:48 latest:** beta-N2 Gibbs potential and corrected primary-report
> reproduction COMPLETE.31 new tests, related116PASS; input8hashes and frozen134
> files/runtime unchanged. Same-source latent-table maximum discrepancy3.4045%
> →0.23333%, NOT a field accuracy result. No calculation active or promotion.
> Batch A remains PARTIAL: next common liquid h/s/G and mixed-phase stability,
> then conservative coupling and actual10/23 fields. Do not repeat N2/O2 audits.
> [Results](nitrogen-phase-potential-results.md). Deadline15:50KST unchanged.

> **Stage2 batch scope updated:** user requested all remaining work as one
> overnight run. Existing degadis heartbeat updated and ACTIVE; no duplicate.
> Read [A--F batch list](stage2-overnight-batch-2026-09-06.md) before selecting
> the next step. Includes phase/source coupling, actual fields, downstream
> transport, regression and delivery. Existing deadline15:50KST maintained.

> **Stage2 09:19:** pure-phase caloric audit and opt-in gamma-O2 Gibbs
> potential COMPLETE.25 new tests, related85PASS, frozen134files unchanged.
> Old implied O2 fusion29.4369 versus measured13.89925kJ/kg corrected in
> the separate potential. Shared-source latent table max discrepancy6.211%
> →0.1554%; NOT a new field error result. No runs active or default promotion.
> Next beta-N2 caloric/volume/fusion reference and mixed-phase consistency.
> [Results](phase-caloric-consistency-results.md).

> **Stage2 08:55:** mixed N2/O2 liquid property screen COMPLETE. Separate
> coefficient-free flash,26 new tests; related total60 pass. Same7 field replays,
>356 local samples;55 in range,10 phase changes. Fixed-h local T shift<=3.624K,
> fixed-T density change<=7.332%. NOT a new field error score or promotion.
> Next: consistent solid/liquid/vapor chemical potentials and thermal ledger.
> [Results](mixed-air-liquid-screen-results.md). No calculation remains running.

> **Stage2 08:50 latest:** continuous-spreading screen COMPLETE (7 cases,
>38/17/41, nine candidate trajectory replays, two-case .01m check).
> MG .821070→.826195, width ratio .952346→.962357; centre MAE slightly worsens
> .061588→.061837m and Tmin MAE17.33598→17.56551K. NOT promoted.
> Fine/coarse sensor-T change<=.00392K. All sessions ended. Next physical
> question: mixed N2/O2 liquid equilibrium versus independent pure condensates.

> **Stage2 current runs:** continuous spreading10/23 .02 COMPLETE; all41
> minimum-T MAE17.56551K versus CONTROL17.33598K, no improvement claim.
> Conservation<=4.331e-8 and full pilot coverage. Remaining5 .02 RUNNING63257;
> representative10/23 .01 RUNNING60082. Session65267 ENDED.

> **Stage2 latest:** matched-flow7/7 fields COMPLETE; exact38/17/41 aggregation
> and14 original/new trajectory replays verified. No source/default selection.
> Located x<=1m ambient-width derivative switch while widths still vary.
> Added opt-in continuous_ambient_spreading plus19 passing independent tests.
> Same pressure-loss CONTROL source trial10/23 field screen RUNNING65267.
> [Active log](stage2-overnight-2026-09-06.md) overrides older running statuses.

> **Stage2 STARTED by new user request — 2026-09-06.**
> User authorized physical error reduction and overnight continuation.
> See [active log](stage2-overnight-2026-09-06.md). The Stage1 134-file bundle
> remains intact. New matched-flow tool:15 tests passed. Coriolis density
> interfaces7/7 pass; enthalpy5/7, trials10/25 fail the unchanged2K gate.
> Density trials10/23 field run completed, all41 thermocouples covered;
> no accuracy improvement or default promotion claimed. Remaining5 density
> trials are running, session20918. Existing heartbeat updated and ACTIVE
> until2026-09-06 15:50KST. Older no-Stage2 instructions refer to the old scope.

> **FIRST MILESTONE COMPLETE — 2026-09-06. Latest overrides ALL older blocks.**
> BASE all7 and CANDIDATE10/23 end-to-end reproductions match frozen outputs
> exactly. Candidate reaches11.415/11.605m, including required sensors<=6m.
> Same38/17/41 comparison completed. Single-ground-image output adapter fixes
> an identified double-reflection defect without editing frozen core/defaults.
> Its portable BASE_SINGLE_IMAGE runner freshly reproduces all7 and all scores.
> Candidate NOT promoted: joint concentration/geometry criteria still fail.
>18 delivery tests passed;134-file manifest and exact runtime verified.
> Full older448 regression and later7 moving-frame tests remain separate claims.
> No active Stage1 calculation remains. Stage2 requires a new user request.
> [Results](stage1-results-2026-09-06.md) · [Run guide](stage1-reproduction-2026-09-06.md).

> **05:47 KST update:** mobile Q+normal source audit completed both synthetic
> cases PASS, six-moment residual<=1.489e-11,105 hashes match. Old fixed-center
> large-Q failure preserved. Source6 tests passed150.67s. Current full regression
>80443 running; do not count unfinished totals. Actual initial turbulence inputs,
> closed normal/pressure/dissipation transport and observation-distance scoring
> remain open. [Results](finite-tke-transport-results.md). Older blocks are history.

> **05:10 KST update:** prescribed-Q initial source energy reallocation also
> passes its independent manufactured audit;95 hashes match, six-moment error
>5.566e-12. All current calculations ended. Four extra source tests passed
> AFTER the full425 run. Next actual initial-Q/mixing/dissipation/circulation
> physics and joint boundary gates; no measured improvement claim.

> **05:00 KST (2026-09-06):** finite-TKE operator and independent full actual
> H+meanKE+Q audit completed, numerical PASS. Manufactured input PSD/boundary/
> production gates FAIL, so no physical adoption. Full425 passed,128 skipped,
>58 deselected. Explicit-Q source energy reallocation tests are running
> separately. [Current results](finite-tke-transport-results.md).

> **04:30 KST (2026-09-06):** exact aligned conservative local witness fully
> verified, including its NEW actual energy derivative; all7 initial shear/TKE
> diagnoses resolved. No running audit. Next: finite-Q modal transport operator
> with explicit dissipation/mixing/inlet inputs. Night continuation remains
> active. No new observed/default claim. [Results](exact-geometry-and-tke-results.md).

> **Latest04:10 KST (2026-09-06):** exact geometry and its own derivatives
> pass independent actual H+K directional checks (max7.18e-8). A new fully
> conservative aligned witness is being re-solved/reverified. Prior aligned
> shear requires minimum TKE flux10.56% of mean axial KE, a coefficient-free
> necessary bound, not actual k or its axial derivative. Non-slow403 passed;
> no new measured accuracy/default adoption. [Results](exact-geometry-and-tke-results.md).

> **Latest around03:45 KST (2026-09-06):** aligned circulation necessary local
> screen passes without new nonradial stress. Amplitudes are NOT physically
> closed. Energy FD remains failed; smooth high-precision component identifies
> cancellation and exact-width versus legacy-tangent inconsistency. All
> calculations ended. Next coherent exact-geometry research verification,
> not default promotion. Full non-slow385 passed; no new observed scores.
> [Latest results](aligned-flux-and-energy-conditioning-results.md).
> Older blocks below describe historical states, not currently running jobs.

> **Latest conservative-flux research (2026-09-06):** initializations7/7,
> short primary segments6/7. Trial10 failure is independently reproduced.
> Rejected hard boundary-row and solenoidal-boundary-moment prototypes stay
> rejected. A sampled admissible flux witness reduces failed-state edges to
>3.54/4.00/3.02% while retaining all weak/global equations on its mesh.
> Rate selector refinement0.00718 and nonradial stress20.7% preclude treating
> it as a validated closure. Independent phase-volume/actual-moment audit
> running. Full non-slow371 passed; no new field score/default promotion.
> [Current results](conservative-flux-witness-results.md).

> **Latest seven-case extension (2026-09-06):** all7 simultaneous
> initializations pass.85 hashes match. Scalar boundaries1.37–3.56%,
> momentum0.10–1.71%, sampled chi_P positive in every case.
> Short moving-state/mixing-update segments are being checked, without
> refit or relaxed gates. Trial10 leaves the heat boundary5% gate near
>0.44 mm; investigate independently before extending or field scoring.
> Non-slow355 passed,128 skipped,58 deselected. No default promotion.
> [Results](coupled-initialization-extension-results.md). Earlier entries are history.

> **Latest simultaneous initialization (2026-09-06):** both preselected
> pilots25/11 pass independent local numerics, six moments, H2/heat/momentum
> boundary5% and sampled positive-viscosity gates. Each trial shape re-solves
> all21 rates; no frozen-flow fit. Trial25 H2/heat16.88/16.54%→3.24/2.30%.
> Trial11 H2/heat3.56/2.81%, momentum0.53%, min chi_P+4.036 s^-1 instead
> of about-6.455. No new physical coefficient or transverse redistribution.
> Remaining10/12/22/23/24 are NOT rerun; no downstream, field or default claim.
> NEXT: extend the same method/gates to those five, then guarded short
> downstream connection and mixing-update closure. Non-slow342 passed.
> [Results](coupled-shape-initialization-results.md). Older entries are history.

> **Latest coupled shape transport (2026-09-05):** new23/37-state local rates
> remove the beta/thermal-mode gauge duplication and rebuild mass/stress fluxes.
> Scalar mixing is still explicitly prescribed, not an evolved closure.
> Original numerical4/7 plus independent fine-input verification gives7/7;
> constitutive adoption0/7. New scalar edge defects are4–17%, and trials11/23
> require negative inferred viscosity near corners. Earlier fixed-transport
>1–3% shape results do not establish this coupled transport. Next: simultaneous
> shape/transport initialization and, if needed, conservative transverse flux
> redistribution before velocity/TKE extensions or downstream integration.
> Non-slow336 passed. Defaults/scores unchanged.
> [Results](enriched-shape-transport-results.md). Older entries are history.

> **Latest conservative edge refit (2026-09-05):** all7 fixed-transport square
> sections pass independent six-moment, boundary and physical gates. New
> face/ray phase splitting fixes the prior angular reference defect; constrained
> refitting reduces internal species/heat edge residuals to0.95–3.63%.
> Maximum independent moment error9.041e-11; gate1e-8. Original failures remain.
> NEXT: reconstruct transport on these new shapes and derive independent shape
> coefficient evolution before any downstream segment/field scoring. No default
> changes. Non-slow:324 passed,128 skipped,58 deselected.
> [Results](edge-conservative-refit-results.md). Older entries are history.

> **Latest edge-shape work (2026-09-05):** square-symmetric non-Gaussian C/H
> enrichment reaches both fixed-transport edge gates in3/7 cases, but
> independent six-moment verification rejects every case: final0/7.
> A separate phase-cell diagnosis confirms under-resolved constraint
> integrals; its angular reference check itself passes only2/7 at1e-8.
> Next: finish angular error control, then conservative constrained refitting.
> No shape-coefficient transport, field score or default promotion.
> Non-slow:317 passed,128 skipped,58 deselected.
> [Results](edge-profile-enrichment-results.md). Older entries are history.

> **Latest reservoir/thermal work (2026-09-05):** an explicit ambient inflow
> energy boundary plus a reduced buoyancy-work convention closes the weak
> two-thermal-moment model. Seven initial states pass independent flux-change
> checks. Short lengths are only 1.759–4.600 mm: 6/7 pass the original RK2
> test, and trial 10 passes a separate adaptive check after retained fixed-step
> failures. Combined 7/7 is numerical/weak-model verification, not field
> accuracy. Gaussian edge species/heat defects remain 12–25%, so profile
> enrichment is next before field expansion. Non-slow: 301 passed,
> 128 skipped, 58 deselected. [Results](reservoir-thermal-moment-results.md).

> **Latest shear/thermal screen (2026-09-05):** reduced mean kinetic-energy
> production/work and thermal zeroth/second moments now couple explicitly.
> Numerics pass 7/7; unity scalar mixing plus immediate shear heating is
> physically incompatible in all seven frozen boundaries. R2 candidates
> are NOT adopted width rates. Finite-edge heat transport dominates the
> defect, so close the ambient/finite-edge energy coupling before adding
> further TKE complexity or running downstream. Non-slow: 287 passed,
> 128 skipped, 58 deselected. Defaults and observed accuracy are unchanged.
> [Detailed result](shear-thermal-compatibility-results.md).

> **Latest conditional mixing step (2026-09-05):** a radial flux ansatz in
> deforming coordinates reconstructs mass/species entrainment and inferred
> anisotropic mixing. The free log thermal-width rate remains unknown.
> Actual source/coflow reproduces all seven stored initial source terms;
> the five tangent balances and four-face mass/species identities pass.
> Independent finite changes of physical fluxes pass 7/7. Dense unsplit
> phase-gradient integration failures are retained separately. Non-slow:
> 273 passed, 128 skipped, 58 deselected. No new trajectory, score or default.
> [Results and explicit assumptions](transverse-conservative-mixing-results.md).

> **Latest thermal-transport step (2026-09-05):** mean enthalpy second-moment
> and reduced derivatives, curved/moving free-rectangle balances and the
> specific-enthalpy-gradient response are implemented. Fixed-order GL fails
> the phase-gradient 1e-5 audit on all seven frozen boundaries, despite very
> accurate Gaussian advective moments. The failed record is retained and an
> error-controlled adaptive implementation passes 7/7 independent checks
> (maximum component difference 2.432e-7; net 1.626e-6). 32 new tests; non-slow
> suite: 257 passed, 128 skipped, 58 deselected. Actual downstream closure,
> ground coupling and new field scoring are still outstanding.
> [Current detailed result](thermal-moment-operator-results.md).

> **Latest boundary result (2026-09-05):** an independent volumetric-enthalpy
> width constrained by the near-field buoyancy moment passes 7/7 interfaces.
> Thermal/species width ratio is 1.03284–1.04818, without observation fitting.
> Five fluxes, force quadrature, temperature and H2-width gates all pass;
> an independent adaptive GK21 audit confirms the fixed-section integrals.
> No downstream beta_H transport law or new field score exists, and the class
> refuses unclosed downstream use. Defaults remain unchanged.
> [Results](buoyancy-constrained-enthalpy-width-results.md) and
> [remaining transport equation](thermal-width-transport-requirements.md).

> **Latest enthalpy-profile screen (2026-09-05):** the same-width Gaussian
> H2/enthalpy-density hypothesis has been implemented and fully compared.
> MG=0.838826, VG=1.190780, FAC2=0.921053; centre-height MAE rises from
> 0.061588 to 0.065565 m relative to the ambient-consistent control.
> Paired 41-sensor minimum-temperature MAE improves 17.3360 -> 16.4738 K,
> but median-temperature MAE worsens. Research-only; no default promotion.
> Five-flux quadrature passes at 63 sampled sections; buoyancy convergence
> is separately reported. No new TKE state or mixture EOS was introduced.
> [Results and limitations](gaussian-enthalpy-profile-results.md).

> **Latest ambient/profile audit (2026-09-05):** the explicit-species phase
> model now has an optional consistent ambient density/composition endpoint.
> The hydrogen-free temperature/enthalpy defect is repaired in that option;
> 7/7 interfaces pass. On all 38 arcs and 17 profiles, concentration changes
> are small and centre-height MAE changes from 0.063587 to 0.061588 m.
> Width ratio worsens from 0.956049 to 0.952346. The joint gate still fails.
> This is a consistency repair, not completion of a mixture real-gas EOS or
> a thermal/TKE transport closure. [Results](phase-ambient-consistency-results.md).

> **Superseding HEM energy audit (2026-09-05):** explicit source phase enthalpy
> and measured-pipe kinetic energy now close the upstream-to-Gaussian ledger
> on four nozzles; all seven fixed interfaces pass. Pure-H2 real-gas screening
> is complete, but no mixture-Z multiplier or TKE model was added. Historical
> HEM/para performance figures below refer to the former boundary and cannot
> be assigned to this corrected research branch. See
> [current results](preslhy-source-eos-ledger-results.md).
> The completed 38-arc/17-profile rerun gives centre MAE 0.0636 m (former
> HEM 0.0700 m; baseline 0.0542 m), MG 0.8202 and VG 1.1991. It does not
> pass the joint promotion gate. Source-energy conservation is repaired;
> the distributed profile/TKE question remains open.

> **Current LH2 status (updated 2026-09-05):** source-state conservation and physical
> ground detachment are implemented and validated. See
> [`lh2-model-improvements-2026-09-03.md`](lh2-model-improvements-2026-09-03.md)
> for the superseding results and remaining condensed-air gap. The PRESLHY
> Gaussian-width convention has also been corrected: current
> model/measurement `sigma_z` is 1.033, not the former 0.731; see
> [`gaussian-width-convention.md`](gaussian-width-convention.md).
> A transported condensed-air source has since been implemented as a
> pre-registered 1/10/100 um sensitivity.  It closes the N2/O2 phase and
> conservation ledgers but worsens plume-height and residual-variance errors,
> so the corrected default remains unchanged; see
> [`prereg-condensed-air-particle-transport.md`](prereg-condensed-air-particle-transport.md).
> A further storage-to-source audit now closes kinetic energy and implements
> an optional HyRAM+ HEM/Yuceil--Otugen pressure-thrust source.  Prospective
> PRESLHY and Spadeadam comparisons reject both as default improvements, and
> nominal-area pressure thrust is thermodynamically incompatible with three
> eligible PRESLHY releases; see
> [`prereg-source-total-energy-pressure-thrust.md`](prereg-source-total-energy-pressure-thrust.md).
> The upstream evaporation-zone coordinate can also be preserved explicitly.
> Its 12--78 mm correction changes no source invariant and is too small to
> repair the field residuals; see
> [`prereg-evaporation-zone-origin.md`](prereg-evaporation-zone-origin.md).
> Li et al.'s zero-condensate-velocity limit is likewise available only as an
> explicit diagnostic.  It closes all source invariants but worsens both
> campaigns, ruling out field-tuning an intermediate axial-slip coefficient;
> see
> [`prereg-stationary-condensate-bound.md`](prereg-stationary-condensate-bound.md).
> A finite-rate applicability screen is now available as well.  It implements
> the Sandia Ranz--Marshall heat-transfer form and a fastest-possible
> sublimation lifetime, but does not invent the missing particle number or
> nucleation rate; see
> [`finite-rate-sublimation-audit.md`](finite-rate-sublimation-audit.md).
> A separate near-nozzle Raman follow-up now includes equilibrium freezing of
> ambient water in the conserved H2 jet. The saturated upper bound is rejected
> for excessive warming, while a 40%-RH sensitivity passes all four printed
> slopes. The source experiment does not report RH, so no humidity is fitted
> or promoted to a default; see
> [`prereg-humid-air-frost-upper-bound.md`](prereg-humid-air-frost-upper-bound.md).
> A follow-on component-enthalpy correction is accepted for the dry-air
> axisymmetric research path. It combines H2/N2/O2/H2O `h(T)` with air phase
> equilibrium, closes all conservation checks and passes all four Raman
> slopes under corrected case coverage; see
> [`prereg-phase-temperature-dependent-enthalpy.md`](prereg-phase-temperature-dependent-enthalpy.md).
> The accepted flags are now exposed through
> `lh2_source_from_measured_throat` and `run_lh2_near_field_research`, with
> conservation diagnostics and applicability warnings attached. Argon phase
> completeness and a perfect-black radiation upper bound were subsequently
> tested prospectively. Both remain explicit sensitivities: argon worsens the
> thermal-centreline error, while full-black radiation slightly worsens both
> mass metrics. See
> [`prereg-argon-phase-completeness.md`](prereg-argon-phase-completeness.md)
> and [`prereg-radiation-upper-bound.md`](prereg-radiation-upper-bound.md).
> The hydrogen spin-isomer caloric assumption has also been isolated. A
> pre-registered para-hydrogen enthalpy run is conservative but worsens the
> corrected Raman centreline-temperature error from -21.38% to -34.79%, so it
> remains an off-default sensitivity; see
> [`prereg-hydrogen-spin-isomer-enthalpy.md`](prereg-hydrogen-spin-isomer-enthalpy.md).
> A fully conservative four-flux/two-scalar phase closure was also tested.
> It reduces the temperature-width error from +24.56% to -2.27% but worsens
> centreline temperature to -29.14%, so it is not the default; see
> [`prereg-phase-two-scalar-four-flux.md`](prereg-phase-two-scalar-four-flux.md).
> Finally, verified ELVHYS Test-10/Test-11 time series have been reduced. No
> quantitative score is issued because the public package lacks H2 mass flow
> and conflicts on nozzle elevation; see
> [`elvhys-tcs-audit.md`](elvhys-tcs-audit.md).
> A source-version audit subsequently replaced the 2017 conference Raman mass
> fits with the final-journal values 0.2771 and 0.07069. The recommended dry
> model remains provisionally 4/4, but its mass-centre and temperature-width
> errors are now -24.86% and +24.56%, very near the 25% boundary. The aggregate
> plots also list an untabulated tenth 4 bar/45 K series, so condition membership
> and fit uncertainty must be obtained before claiming raw-data validation; see
> [`hecht-panda-journal-benchmark-correction.md`](hecht-panda-journal-benchmark-correction.md).
> The formerly missing conservative near-field/crosswind interface is now
> implemented for horizontal wind-aligned releases. Two pre-registered
> reduced-thermodynamic mappings were rejected; transferring the full
> near-field N2/O2/H2O phase and component-enthalpy radial profile passes at
> 0.08 m with 0.0585% energy, 1.57% H2-width and 0.94 K centre-temperature
> mismatch while all four native balance residuals remain below `6e-16`.
> A strict 81-point/0.25-mm repeat retains the pass; see
> [`prereg-handoff-grid-convergence.md`](prereg-handoff-grid-convergence.md).
> Failed interfaces cannot run downstream. This is interface verification,
> not yet PRESLHY/Spadeadam validation; see
> [`conservative-nearfield-crosswind-handoff-results.md`](conservative-nearfield-crosswind-handoff-results.md).
> PRESLHY testing of that coupling now selects seven momentum-dominated
> horizontal releases. Adaptive energy quadrature removes a numerical
> 32/64-point false rejection and all seven scalar-peak interfaces pass.
> Concentration improves (MG 1.074 to 0.925; VG 1.213 to 1.183; FAC2 0.952
> unchanged), but vertical width worsens from 1.091 to 0.710 of measured, so
> the coupled path is not promoted. A four-flux source boundary closes total
> mass as well but does not repair width and passes only five of seven
> handoffs. See
> [`prereg-preslhy-coupled-crosswind-validation.md`](prereg-preslhy-coupled-crosswind-validation.md)
> and
> [`prereg-preslhy-four-flux-establishment.md`](prereg-preslhy-four-flux-establishment.md).
> A subsequent mechanism-isolation test identifies the downstream width loss:
> switching from constant source-momentum entrainment to JETPLU's local
> density-scaled shear term at the handoff moves model/measured `sigma_z` from
> 0.710 to 0.994 while retaining 7/7 interfaces, 0.049 m centre-height MAE and
> FAC2 0.952. This mechanism is supported but the scalar source boundary is
> not a fully four-flux closure. Combining it with four-flux establishment
> passes only 5/7 interfaces; a 10D--20D compatibility search does not rescue
> trials 10 and 25. The remaining structural requirement is an independent
> transported thermal/energy profile state. See
> [`prereg-preslhy-posthandoff-local-shear.md`](prereg-preslhy-posthandoff-local-shear.md)
> and
> [`prereg-preslhy-compatible-handoff-search.md`](prereg-preslhy-compatible-handoff-search.md).
> The resulting seven-state model structure is now frozen before coding or
> field scoring in
> [`prereg-crosswind-independent-energy-state.md`](prereg-crosswind-independent-energy-state.md).
> Its first implementation gate now passes all seven fixed 10D boundaries.
> Former single-scalar failures 10 and 25 are recovered without changing a
> threshold; maximum five-flux residual is `6.45e-16`, H2-width mismatch 4.761%
> and centre-temperature mismatch 1.990 K. This boundary result validated the
> missing-state diagnosis; see
> [`preslhy-independent-energy-interface-results.md`](preslhy-independent-energy-interface-results.md).
> The downstream five-balance ODE is now complete and step-refined. It retains
> 7/7 interfaces and improves VG 1.213 to 1.176, FAC2 0.952 to 0.976 and width
> ratio 1.091 to 1.066, with maximum balance residual `1.04e-7`. It is not
> promoted because MG worsens 1.074 to 1.118 and centre-height MAE 0.054 to
> 0.078 m. The remaining error is buoyant rise, not width or numerical energy
> closure.
> A common sensor-height ground-image observation audit now passes 17/17 fits
> but confirms the rejection (candidate width ratio/centre MAE 1.420/0.137 m,
> baseline 1.360/0.070 m). The Li (2026) enthalpy-only established-flow balance
> also passes 7/7 and conservation but changes centre MAE by less than 0.04 mm.
> A direct section-buoyancy audit then finds the near-field and projected
> forces have the same sign and differ by at most 6.761% at all seven 10D
> interfaces. A HyRAM source-profile audit then finds that Houf `B` must be the
> velocity e-folding width, smaller than the stored scalar width by 1/1.16.
> Correcting it improves MG/VG to 1.062/1.163 and internal width ratio to
> 0.976, but internal centre MAE remains 0.0758 m versus 0.0542 m baseline.
> The corrected equation is retained research-only; the remaining target is
> downstream buoyancy accumulation, especially trial 10. A direct trial-10
> budget now shows the candidate accumulates 5.695 N buoyancy from 10D to 6 m
> versus 2.734 N for baseline, while their 10D vertical-momentum difference is
> only 0.0657 N. The next separation is early mass entrainment versus
> instantaneous equilibrium air-condensation heat release.
> The coefficient-free metastable dry-air bound now resolves that separation:
> trial 10 stays at 65.559 K and weakly buoyant, while its projected temperature
> differs by 2.222 K and fails the unchanged interface gate. No seven-trial
> field score was run. Equilibrium condensation heat release is not the source
> of the 10D thermal contrast; pre-10D total-mass entrainment is the next target.
> See
> [`prereg-delayed-dry-air-condensation-bound.md`](prereg-delayed-dry-air-condensation-bound.md).
> A source-zone audit then finds that `entrained_mass` had also integrated the
> Zone-V `beta_A` source across the Zone-IV 6.2D profile transition, although
> the PRESLHY input is already an air-loaded Station-3 state. The conservative
> `source_flux` lower bound removes this second application, passes 7/7 and
> improves MG/VG/FAC2 to 1.012/1.162/0.976. Centre MAE improves to 0.0707 m
> but remains above baseline 0.0542 m, and width error is slightly worse than
> baseline, so it is retained research-only. See
> [`prereg-source-flux-gaussian-establishment.md`](prereg-source-flux-gaussian-establishment.md).
> The implementation is regression-clean: the complete non-slow suite passes,
> and the 13 selected long-running phase/energy conservation tests also pass.
> Two final coefficient-free bounds do not close the remaining trajectory
> residual. Ground contact begins only after 2.60 m in trial 10 and lowers its
> 6 m centre from 0.737 to 0.695 m versus 0.227 m measured. Applying Li equation
> 35 from Station 4 rather than 10D fails the conservative Gaussian boundary
> on the representative source. Neither proceeds to a seven-trial score; see
> [`prereg-independent-energy-ground-contact.md`](prereg-independent-energy-ground-contact.md)
> and [`prereg-zone-v-enthalpy-nearfield.md`](prereg-zone-v-enthalpy-nearfield.md).
> Direct thermocouple validation is now available from the original trial-10
> and trial-23 workbooks. Against the report's minimum-temperature statistic,
> the 42-point `source_flux` branch has +20.31 K median bias and 20.31 K median
> absolute error; trial 23 remains +37.55 K warm in the median. Raising the input from sustained mean
> to the measured-window peak lowers the aggregate error to 16.29 K and moves
> MG to 1.000, but worsens VG from 1.162 to 1.179. It therefore fails its
> pre-registered promotion gate and remains an uncertainty bound. The default
> is unchanged. A temporal-median robustness check removes the aggregate
> signed bias but not the centreline error: centreline median absolute error
> is still 32.75 K and Trial 23 remains +37.26 K warm. The next model target
> is synchronized source unsteadiness followed by a separate liquid-H2
> thermal/mass state with finite-rate interphase transfer, not an empirical
> source-rate multiplier. See [`prereg-preslhy-temperature-field.md`](prereg-preslhy-temperature-field.md)
> and [`prereg-preslhy-source-rate-consistency.md`](prereg-preslhy-source-rate-consistency.md).
> A synchronized 0--5 s lag audit then finds the expected source-flow sign in
> all eight Trial-23 centreline channels, but the common temperature contrast
> reaches only 6.00 K and fails the frozen 10 K gate. Trial 10 has the opposite
> sign because the apparatus continues cooling through the window. Source
> unsteadiness is therefore a future time-series feature, not the primary
> cold-core repair; see [`prereg-preslhy-transient-source-temperature.md`](prereg-preslhy-transient-source-temperature.md).
> Separating liquid temperature from tanker pressure confirms a second source
> defect but cannot be applied as a temperature-only patch.  The 20.3689 K
> boiling-liquid bound makes all seven pressure-thrust-free source states fail
> the fixed velocity/wind threshold (0.929--5.892 versus 10), because their
> pressure work is discarded.  The candidate is rejected before field
> scoring; the next implementation reconstructs a measured two-phase nozzle
> state from pressure, temperature/quality, density and effective area. See
> [`prereg-preslhy-subcooled-liquid-source.md`](prereg-preslhy-subcooled-liquid-source.md).
> The subsequent seven-workbook reconstruction now separates measured TC3
> temperature and PT2 pressure, retains post-flash liquid H2, and applies the
> Yellow Book atomisation range without fitting. The four identifiable nozzle
> sources remain 98.08--100% liquid immediately after flashing. A collective
> phase-correct fast bound passes all seven interfaces and improves campaign
> VG and width, including a large trial-23 concentration improvement, but
> worsens overall MG, FAC2 and centre-height MAE. It is therefore retained as
> a research bound rather than the default; see
> [`prereg-preslhy-measured-nozzle-source.md`](prereg-preslhy-measured-nozzle-source.md)
> and
> [`prereg-preslhy-finite-rate-lh2-droplets.md`](prereg-preslhy-finite-rate-lh2-droplets.md).
> A subsequent spin-consistent source audit applies `ParaHydrogen` from the
> compressed liquid through flash, collective evaporation and conserved
> near-/far-field calorics. It passes 7/7 interfaces and reduces the
> fast-bound centre-height MAE from 0.07004 to 0.05453 m, but worsens VG,
> FAC2 and width ratio. It is therefore an explicit composition uncertainty,
> not a promoted correction; see
> [`prereg-preslhy-liquid-spin-source.md`](prereg-preslhy-liquid-spin-source.md).
> The exact experimental channels needed to resolve these states, including
> trial numbers, units and minimum useful subsets, are in
> [`data-request-lh2-remaining-state.md`](data-request-lh2-remaining-state.md).
> A pre-registered independent Spadeadam pilot was then attempted with the
> existing 1 um conserved source reconstruction. Test 6 passes the extended
> phase-manifold interface, but test 4 fails the width/temperature screens and
> starts below the velocity/wind ratio limit (9.091 versus 10). Per protocol,
> no partial sensor score was issued. See
> [`prereg-coupled-crosswind-spadeadam.md`](prereg-coupled-crosswind-spadeadam.md).

Written at the point where the DEGADIS work is finished and the liquid
hydrogen work begins.

## What exists

A Python reimplementation of **DEGADIS 2.1**, the dense gas dispersion model
of Spicer and Havens released by EPA in 1989 and still named in 49 CFR
193.2059 for LNG exclusion zones. There was no other open Python
implementation: a GitHub-wide search returns one repository, a C# front end
that shells out to the original DOS executables without porting the physics.

All six of DEGADIS's programs are ported and all five EPA test cases run end
to end in Python. 136 tests, 35 modules, about 8 000 lines.

| stage | original | status |
|---|---|---|
| source blanket | `DEG1`, `SRC1`, `SRC1O` | parity |
| source thinning | `CRFG`, `SZF`, `NOBL` | parity |
| steady downwind | `DEG2S`, `PSS`, `SSG` | parity |
| transient observers | `DEG2`, `SSSUP`, `OB`, `UIT`, `TUPF` | parity |
| cloud snapshots | `DEG3`, `GETTS`, `SORTS` | parity |
| receptor histories | `DEG4`, `GETTD`, `DOSOUT` | done |
| jet/plume | `JETPLU`, `SETJET`, `ELLIPS` | parity |
| jet-to-ground bridge | `DEGBRIDG` | parity |
| input/transfer files | six formats | parity |

```python
from degali import run_steady, run_transient, run_jet_to_ground
profile, source = run_steady("B9.INP")
profile.distance_to(0.05)      # metres to the LNG lower flammable limit
```

## How it was validated against the original

The Fortran is built from source **inside the repository** and run as part of
the test suite. Six purpose-written probes link against its object files and
dump internal state at full double precision, because printed output carries
five significant figures and cannot distinguish a correct port from a close
one.

Two layers: the patched Fortran must still reproduce EPA's 2012 golden
listings (it does, byte for byte apart from timestamps), and the Python must
reproduce the Fortran at 1e-12.

Eight portability patches were needed for gfortran. One was load-bearing:
`COMMON /ERROR/` interleaves eighteen `REAL*8` with an `INTEGER*4`, and
`ESTRT1`'s `EQUIVALENCE` overlay made gfortran pad that block differently
there than in every other unit, so the air-entrainment coefficient was
silently read as zero. Intel's `/align:dcommons` had been masking it.

### Three things that trip up a port

Reproduced rather than fixed, because the point is to reproduce DEGADIS.

- **`ADIABAT` computes `wa` for `ifl` of 0, −1, −2 and 2, but not for 1**
  (`TPROP.for` lines 300–320), so on that path it uses whatever the caller
  passed. `SZF.for` line 95 passes `walay`, which appears nowhere else in that
  file — no declaration, no assignment. Under `/noauto` that is static,
  zero and never written, so `SZF` runs every lookup as though the mixture
  contained no dry air.
- **`SURFAC` receives a mix of heated and adiabatic layer properties**, and
  `PSS` and `SSG` differ by one letter (`temlay` against `temlam`) in which
  temperature they pass.
- **`GAMINC` returns the unregularised incomplete gamma.** Reading it as the
  regularised form leaves the flammable-mass derivative low by exactly
  Γ(1/(1+α)) while everything else in the routine stays correct.

### One reachable model limit

Thorney Island's cylinder has a height-to-diameter ratio near one, which turns
on the van Ulden momentum balance — a branch no EPA test case reaches. It does
not converge, and `SRC1` stops. The original Fortran stops in the same place
with `STOP SRC1 velocity loop`, confirmed by building the deck and running the
1989 executable on it.

## Two backends

`legacy` is bit-faithful to 1989 including its approximations: the adiabatic
mixing table is built on a grid computed in single precision, `ALPH`'s
quadrature error moves the fitted wind exponent more than its own root-finder
tolerance, `GSERIES` clamps its argument, and `ZBRENT` cannot converge below
`EPS = 3e-8`.

`coolprop` replaces the 1989 correlations with equations of state and solves
the thermodynamic inversions accurately. On Burro 9 it moves the contaminant
enthalpy by 6 % and the final distances by half a per cent — the density
change enters gravity slumping as a square root and entrainment suppression
through Φ, and the two largely cancel. **Against measurement the two backends
are indistinguishable**, which is itself a result: the 1989 correlations are
not what limits DEGADIS's accuracy on LNG.

## What the field evaluation found

Against the Risø REDIPHEM database. The headline is not that DEGADIS is
accurate or inaccurate, it is that **the conventional way of evaluating it
hides a structural error**.

At the lowest instrumented height — what the published evaluations report —
Burro gives MG = 0.81 with a 95 % interval of [0.63, 1.02]. That straddles
Hanna's bias bound rather than satisfying it.

Using every instrumented height instead:

| Burro, by elevation | n | MG | FAC2 |
|---|---|---|---|
| 1 m | 61 | 0.81 | 0.56 |
| 3 m | 59 | 7.3 | 0.19 |
| 8 m | 53 | 5500 | 0.09 |

The measured cloud is close to well mixed over the bottom eight metres from
50 m downwind; the model confines it under a metre and a half. The lateral
spread runs the other way, one to three times too wide against Desert
Tortoise. The single-height statistic is carried by the one elevation where
those two errors cancel.

Two results were withdrawn during this work, both worth recording:

- **FLADIS.** The series ships no `CHANDEF.DAT`; a fallback list of channel
  numbers selected channels reading 302 and 23.5 at 20 m, which are not
  concentrations. The result looked clean and meant nothing. The reader now
  refuses to guess.
- **The Coyote explanation.** The bias was attributed to a pool radius taken
  from a spill pond. Forcing smaller diameters makes it worse; the hypothesis
  is falsified and no replacement is offered.

## Repository

```
src/degali/          core/ io/ validation/, run.py, cli.py
reference/fortran/     patched DEGADIS 2.1 + build.sh + six probes
reference/testcases/   the five EPA cases with their golden output
tests/                 136 tests
docs/validation.md            against the Fortran
docs/field-validation.md      against measurement
docs/data-inventory.md        what data exists and what state it is in
```

Field-trial tests need `$REDIPHEM_ROOT` and skip cleanly without it, since the
database is not redistributable. `CITATION.cff`, a CI workflow and a
changelog are in place; wheel and sdist both build, with the reference
implementation shipped in the sdist because the validation claim is not
reproducible without it.

## Where this goes next

The LNG and ammonia work is done. The open question is **liquid hydrogen**,
for which the project folder holds the PRESLHY and HSL datasets and the
supporting literature. LH2 is a genuinely different problem from LNG —
hydrogen vapour is buoyant at ambient temperature but the cold cloud is
initially dense, so a release passes through a dense phase into a buoyant one,
and DEGADIS has no buoyant-rise physics at all. See
`docs/data-inventory.md` for what is available and
`docs/lh2-plan.md` for what can and cannot be attempted.
