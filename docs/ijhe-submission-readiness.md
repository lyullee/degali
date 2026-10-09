# IJHE submission-readiness record

**Target:** *International Journal of Hydrogen Energy*
**Scope:** the observation-operator-aware LH2 validation manuscript and its
reproducible DEGALI research software
**Last reviewed:** 2026-10-09

This record separates what is already defensible in a manuscript from what
still requires a new experiment or an accountable evidence package. It is a
submission-control document, not a claim that the model is certified for
safety-critical design.

## Journal fit

Elsevier describes IJHE as publishing original analytical and experimental
research across hydrogen energy, including storage, transmission, utilization
and enabling technologies. The current journal listing also includes a
Safety/Sensors topic. The present manuscript is therefore framed around LH2
release safety, sensor observation operators and evidence-bounded engineering
screening rather than around a generic software release. See the
[official journal description](https://shop.elsevier.com/journals/international-journal-of-hydrogen-energy/0360-3199)
and the [current journal topic listing](https://www.sciencedirect.com/journal/international-journal-of-hydrogen-energy/vol/210/suppl/C).

## Publisher submission controls checked

The package also follows the current general Elsevier author resources: the
[Highlights guidance](https://www.elsevier.com/researcher/author/tools-and-resources/highlights)
specifies three to five bullets with no more than 85 characters each; the
[graphical-abstract guidance](https://www.elsevier.com/researcher/author/tools-and-resources/graphical-abstract)
asks for an original, separate visual file and defers journal-specific artwork
rules to the journal guide; and the [research-data guidance](https://www.elsevier.com/researcher/author/tools-and-resources/research-data)
recommends an explicit data-availability statement where the journal requests
one. The repository audit checks the first and third items and records the
graphical abstract as a separate SVG source; final file-format and pixel-size
acceptance must still be confirmed in the submission portal.
The current graphical-abstract SVG is 1400 x 760 px, which exceeds the
Elsevier general minimum of 1328 x 531 px. Its 1.84:1 layout is intentionally
kept unclipped for the three-panel workflow; the preferred 500:200 aspect
ratio is recorded by the figure audit and should be confirmed in the portal
before final upload rather than silently distorting the artwork.
The manuscript also uses a separate AI declaration immediately before the
references, consistent with Elsevier's [generative-AI journal policy](https://www.elsevier.com/en-gb/about/policies-and-standards/generative-ai-policies-for-journals);
the code-development use is additionally described in the reproducibility
section. The manuscript body uses sequential numeric citations in square
brackets, and the submission audit verifies that the numbered bibliography is
sequential and fully cited.

## Portal preflight matrix

| Upload item | Current package state | Final action |
|---|---|---|
| Article type and title page | Research-article framing; author-neutral title page with main-text and abstract word counts | Confirm article type, author order, affiliations, postal addresses and corresponding-author details in the portal |
| Abstract and keywords | 165-word abstract; six keywords; machine-checked | Confirm the portal's live word limit and paste the final author-approved text |
| Highlights | Five bullets; maximum 76 characters; standalone Word file plus source text | Upload `Highlights.docx` as the portal requests; retain the source text for audit |
| Graphical abstract | Separate SVG and PDF; 1400 x 760 px; minimum-dimension check passes | Confirm the portal accepts the aspect ratio and final file type |
| Main figures and tables | Six inline manuscript figures, two tables, seven audited PDF derivatives | Confirm figure permissions and portal artwork format; no re-export is currently required |
| Supplementary information | Separate five-page DOCX and PDF derivatives | Upload as supplementary files, not as additional main-text pages |
| Data and code statements | Conditional raw-data access, public code DOI and the pinned `requirements-research.txt` runtime list are stated | Confirm licence/custodian wording and the exact release DOI after author approval |
| Declarations | Funding, CRediT, competing interests, acknowledgements and AI fields are present but author-neutral | Replace all author-confirmation placeholders before final export |
| Submission history and permissions | Exclusive-submission, preprint/prior-publication and third-party permissions fields remain open | Complete the cover-letter and portal declarations |

This matrix is a submission control, not a substitute for the journal's live
portal validation. The current Elsevier resources require a separate graphical
abstract where requested and recommend an explicit data statement; the journal
page also lists Safety/Sensors among its current hydrogen-energy topics.

## Recommended paper position

The paper should be presented as:

> A provenance-controlled, observation-operator-aware fast screening framework
> for momentum-dominated outdoor LH2 releases, with an explicit evidence
> boundary and a decision-oriented comparison of peak and time-synchronised
> concentration operators.

The central contribution is an independently framed DEGALI validation and
application envelope. The source state, wind record, sensor operator and
applicability status are kept explicit so that a numerical result cannot
silently become a safety-distance claim. No cross-code performance rank is
included in the manuscript.
In practice, this supports fast HAZID/pre-FEED scenario ranking, sensor-array
planning and emergency-response triage, while routing evidence-gated or
site-specific cases to transient multiphase CFD or a qualified consequence
study. The industrial role is therefore a provenance-controlled front-end and
decision aid, not an automatic separation-distance generator.

The revised Introduction now positions the work against recent IJHE activity:
the digital-twin safety review by [Naanani et al.](https://doi.org/10.1016/j.ijhydene.2025.02.440),
the LH2-release review by [Liang et al.](https://doi.org/10.1016/j.ijhydene.2025.04.026),
the physics-guided confined-dispersion model by [Wang et al.](https://doi.org/10.1016/j.ijhydene.2026.155179),
and sensor-efficient indoor source-term estimation by [Li et al.](https://doi.org/10.1016/j.ijhydene.2026.156184).
The manuscript explicitly distinguishes those data-rich learned or inverse
systems from its no-fit outdoor LH2 screening and evidence-boundary focus.

The 3-D transient dense-gas solver is a separate research extension. Its
numerical verification may be reported as an implementation subsection or
supplement, but it must not be pooled into the field-validation score. A
conditional common-clock PRESLHY replay lane now covers accepted Trials 10, 20
and 21.  The associated D3.6 report gives nominal sensor accuracies,
source-relative receptor coordinates/heights, common-logging details and
obstruction-test offsets, which supports an explicit uncertainty and geometry
boundary; it does not remove the event-level calibration-certificate and
surveyed-site requirements for an operational claim. See
[`docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md`](ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md).
The public RADAR downloads for Trials 20 and 21 are byte-for-byte matched to
the local replay inputs in
`outputs/preslhy-public-retrieval-index-2026-10-09.json`; this improves input
traceability but does not close the event-level calibration, absolute-clock or
surveyed-geometry gates.
That index also pins the non-redistributed D3.6 report used for the replay
(6,518,995 bytes; SHA-256
`865b2b9f966e05023fbb581a4ed68f95f2353d6d6db2ad137aae7a9cf17c4f67`; pages
14–15, 43–46 and 50), so the nominal instrument and geometry provenance is
auditable without overstating validation completeness.
The local evidence archive also contains the public ELVHYS WP4.2 dataset,
which supplies time-resolved confined H2 channels, declared calibration
specifications and measured TCS sensor/nozzle geometry. It is added as a
conditional dynamic/geometry evidence lane, not as a new outdoor accuracy
score, because the archive has no measured H2 mass-flow history and includes
walls, vents and active/passive ventilation. See
[`docs/ijhe-elvhys-conditional-evidence-2026-10-09.md`](ijhe-elvhys-conditional-evidence-2026-10-09.md).
The separate inspection-centre workspace was also inventoried. Its 17 XLSX
and 3 CSV files contain process PT/TT/DPT histories and operating sequences,
but no atmospheric wind, fixed-receptor H2 concentration, surveyed obstacles,
or approved source-rate history. It is therefore retained as operational
source context only; see
[`docs/ijhe-center-operational-data-audit-2026-10-09.md`](ijhe-center-operational-data-audit-2026-10-09.md).

## Claim-to-evidence matrix

| Claim | Status | Evidence that must be cited | Remaining action |
|---|---|---|---|
| The frozen LH2 path is reproducible without a fitted validation coefficient | **Ready** | `docs/lh2-paper-baseline-2026-09-20.md`; `tools/audit_lh2_paper_baseline.py` | Freeze the exact input/provenance manifest used for submission |
| E3.5 peak-operator score is MG 1.047, VG 1.425, FAC2 0.839 on 62 arc maxima | **Ready with scope** | `MANUSCRIPT_LH2_VALIDATION.md`; baseline audit output | State trial count and correlation explicitly; do not call arcs independent trials |
| Peak and synchronised 20 s mean operators give different conclusions | **Ready with scope** | `docs/time-aligned-uncertainty-results-2026-10-03.md` | Keep both operators in the main table and explain the 18 s sampling-line delay |
| FFI/DNV is an independent far-field screen | **Conditional** | `docs/ffi-test6-decomposition.md`; FFI audit outputs | Retain the low-wind Test 6 mismatch as a limitation; do not fit a correction |
| Cross-code performance ranking | **Excluded from submission** | manuscript, supplementary information and graphical abstract | Keep the paper focused on DEGALI; disclose the separate related manuscript to the editor in the cover letter |
| The 3-D transient solver is numerically conservative and convergent | **Ready as numerical verification** | `docs/transient-dense-gas-3d.md`; `outputs/transient-3d-error-study-2026-10-08/error_study.md` | State that the reference is numerical, not measured validation |
| The model predicts site-specific obstacle wakes or complex facilities | **Not supported** | `docs/field-deployment-roadmap-2026-10-05.md`; `docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md`; obstacle gate records | The FFI pad has known structures, but the six-arc lane does not resolve them; do not use this as an IJHE headline claim without matched measured obstacle data |
| The field workflow is operationally deployable | **Conditional** | `docs/field-evidence-manifest.md`; `docs/field-conditional-review.md` | Require source/weather/receptor/sensor/operator identities and SHA-256 fingerprints |
| A common-clock LH2 event can be reconciled without hidden joins | **Conditional** | `docs/ijhe-preslhy-trial10-conditional-evidence-2026-10-09.md`; derived `FieldEvidenceManifest` | Obtain channel calibration and measured site geometry before operational promotion |
| ELVHYS WP4.2 supplies time-resolved confined H2 and measured TCS geometry | **Conditional boundary** | `docs/ijhe-elvhys-conditional-evidence-2026-10-09.md`; `docs/elvhys-tcs-audit.md` | Resolve source mass-flow boundary and enclosure geometry discrepancy before quantitative promotion |
| The model is suitable as a regulatory separation-distance basis | **Withheld** | `docs/publication-scope.md`; README warning | Keep the prohibition in the abstract, conclusions and software README |

The upload bundle also contains separate [Highlights](../IJHE_HIGHLIGHTS.md)
source text and the standalone `Highlights.docx`,
[Supplementary Information](../IJHE_SUPPLEMENTARY_INFORMATION.md),
[graphical abstract](../IJHE_GRAPHICAL_ABSTRACT.svg) and
[cover-letter draft](../IJHE_COVER_LETTER_DRAFT.md) files. The figure/table
generation register is in [docs/ijhe-figure-table-register.md](ijhe-figure-table-register.md).
An author-neutral assembled hand-off is generated by
`tools/build_ijhe_upload_bundle.py` into
`outputs/IJHE_UPLOAD_DRAFT_LATEST_2026-10-09_DEGALI_ONLY_V2/`; its bundle manifest records hashes and
keeps `final_upload_allowed=false` until author and evidence gates are closed.
The supplied legacy related-submission Word files are explicitly separated in
[docs/ijhe-word-artifact-boundary.md](ijhe-word-artifact-boundary.md) and are
not the DEGALI upload source.
The required author, funding, declaration and permission fields are collected
in [docs/ijhe-author-submission-intake.md](ijhe-author-submission-intake.md).
The combined hand-back sheet for those fields and the matched-event evidence
gate is [docs/ijhe-closeout-worksheet.md](ijhe-closeout-worksheet.md).
The final source and derived-figure hash record is generated by
`tools/build_ijhe_submission_manifest.py` into the ignored
`outputs/ijhe-submission-manifest-2026-10-09.json` file.
The requirement-by-requirement completion record is in
[docs/ijhe-final-gap-audit-2026-10-09.md](ijhe-final-gap-audit-2026-10-09.md).
The exact export request for the remaining calibration- and geometry-complete
matched event is in
[docs/ijhe-field-evidence-request-2026-10-09.md](ijhe-field-evidence-request-2026-10-09.md).
The non-LH2 SMEDIS/REDIPHEM candidate package and its promotion decision are
recorded separately in
[docs/ijhe-smedis-transfer-boundary-2026-10-09.md](ijhe-smedis-transfer-boundary-2026-10-09.md);
those generic-gas records are not pooled with the LH2 scores.
The Highlights and Supplementary Information are content-ready. The cover
letter contains explicit author placeholders and must be reviewed before
submission; it is not treated as a completed declaration by the automated
audit.

## Minimum IJHE revision package

Before submission, the paper package should contain the following items.

1. **Frozen manuscript inputs.** Archive the exact model configuration,
   source-state identifiers, weather assumptions, receptor geometry and
   observation-operator definitions used in every table.
2. **Independent validation separation.** Keep E3.5, FFI/DNV and Hecht–Panda
   in separate evidence lanes. Do not pool them into one accuracy number or
   introduce a cross-code performance rank.
3. **Matched transient validation.** Promote at least one event with a
   calibration-complete common source/weather/receptor/sensor clock. The
   preferred unit of holdout is the
   release event or trial, not an individual sensor row. ELVHYS can document
   dynamic and measured-geometry handling, but its ventilation-flow channel is
   not a measured H2 source rate and must remain conditional.
4. **Uncertainty accounting.** Report source, wind, clock/operator and sensor
   calibration uncertainty separately from discretisation error. A grid
   convergence result is not an experimental accuracy interval.
5. **Failure boundary.** Retain the Test 6 underprediction, the
   Hecht–Panda transfer boundary and the obstacle/2-phase limitations in the
   main manuscript, not only in supplementary material.
6. **Reproducibility.** Provide the code snapshot, environment, SHA-256
   manifest and commands below. Do not redistribute third-party raw workbooks,
   PDFs or proprietary Fortran sources unless their licences permit it.

## Reproduction commands

The commands below are intentionally independent. The first two require only
the public repository; the LH2 evidence audits additionally require the
locally controlled data paths described in their documentation.

```powershell
.venv\\Scripts\\python.exe tools\\audit_ijhe_submission.py
.venv\\Scripts\\python.exe tools\\build_ijhe_docx.py `
  --output outputs\\ijhe-manuscript-draft-YYYY-MM-DD.docx
.venv\\Scripts\\python.exe tools\\check_publication.py
.venv\\Scripts\\python.exe -m pytest -q -m "not slow" --ignore=tests/test_reference_parity.py
.venv\\Scripts\\python.exe tools\\benchmark_transient_dense_gas_3d.py `
  --output-dir outputs\\transient-3d-error-study-2026-10-08
.venv\\Scripts\\python.exe tools\\plot_transient_3d_error_study.py `
  outputs\\transient-3d-error-study-2026-10-08\\error_study.json `
  --output outputs\\transient-3d-error-study-2026-10-08\\transient-3d-convergence.svg
```

The 3-D error study reports a fine-grid aggregate RMSE of
`5.4386e-4 kg H2/m3`, a maximum absolute receptor error of
`1.7501e-3 kg H2/m3`, an apparent aggregate order of approximately 1.47 and a
maximum mass residual near machine precision. These values quantify numerical
resolution and conservation only; they are not field accuracy or sensor
accuracy.

## Reviewer-risk controls

- **Novelty risk:** describe the observation operator, provenance gate and
  decision consequences as the contribution; do not claim a new universal
  dense-gas closure.
- **Validation risk:** show trial-level and operator-specific results, and
  keep the low-wind FFI residual visible.
- **Overclaim risk:** use “screening”, “conditional” and “evidence-bounded”
  consistently; reserve “validated” for the declared operator and release
  envelope only.
- **Reproducibility risk:** make every reported table traceable to a script,
  input fingerprint and software version; reject local outputs containing raw
  third-party data from the public snapshot.
- **Industrial-use risk:** state that the workflow supports scenario
  screening and emergency-response analysis, but is not a standalone design
  basis or regulatory approval tool.

## Decision

The package is ready for an **IJHE manuscript draft and internal co-author
review**. It is not yet ready for a claim of universal transient or obstacle
validation. The blocking item for a stronger operational claim is external,
matched, time-resolved field evidence—not another fitted coefficient or a
finer numerical grid.
