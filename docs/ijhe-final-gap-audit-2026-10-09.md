# IJHE final gap audit — 2026-10-09

This is a completion audit for the current DEGALI IJHE package. It records
what is proved by local evidence and what still requires an author or an
external evidence change. A green software test does not close a missing
field-data or submission-metadata requirement.

| Requirement | Current evidence | Status | Remaining action |
|---|---|---|---|
| IJHE hydrogen-energy fit | Official journal scope; LH2 safety, dispersion and sensor/operator framing in the manuscript | Ready | Keep the title/abstract focused on hydrogen safety and evidence-bounded screening |
| Scientific novelty | No-fit outdoor LH2 path; explicit peak/20 s observation operators; source/weather/sensor provenance gate; recent IJHE literature positioning | Ready with scope | Do not recast the work as a universal CFD or AI-surrogate replacement |
| Primary validation | 62 E3.5 arc maxima from 9 trials; MG/VG/FAC2 = 1.047/1.425/0.839 | Ready with scope | Preserve trial clustering and operator definition |
| Independent far-field check | Six FFI/DNV arcs; MG/VG/FAC2 = 1.245/1.373/0.833; Test 6 3.48-fold residual retained | Conditional | Keep the low-wind limitation visible; no fitted correction |
| Sensor uncertainty | Supplementary S9 separates calibration, delay, display precision, geometry and observation operator; D3.6 records nominal Xensor/Dräger/flow specifications separately from event calibration | Ready as evidence boundary | Add channel-level calibration only if supplied by the data custodian |
| Statistical uncertainty | Trial-cluster percentile intervals now accompany MG, VG and FAC2 for the paired-path run; the corrected E3.5 intervals are MG 0.771–1.416, VG 1.132–1.882 and FAC2 0.677–0.983 | Ready with scope | Keep intervals descriptive; do not present them as source or model-form confidence bounds |
| 3-D transient extension | Numerical convergence figure, RMSE `5.4386e-4 kg/m3`, order 1.47, mass residual near machine precision | Ready as numerical verification | Do not call this experimental validation |
| Obstacle/wake claim | FFI report geometry is now explicitly documented (containers, drum and instrument box), but DEGALI does not resolve it and no matched LH2 obstacle event is promoted | Withheld | Obtain machine-readable measured geometry with a matched source/weather/receptor/sensor package |
| FieldEvidenceManifest gate | Five channels and SHA-256s are pinned; sensor-set/operator identities are explicit; non-certified or nominal-only sensor metadata now remains `conditional` | Ready as evidence gate | Supply a custodian-confirmed certificate and declare `sensor_calibration_status=certified` for an operational package |
| Reproducibility | Audit scripts, figure generators, public snapshot check, SHA-256 submission manifest, the stable DEGALI concept DOI (`10.5281/zenodo.22646258`) and the v0.3.0 release DOI (`10.5281/zenodo.23256538`) | Ready with release boundary | Regenerate the manifest after final author edits and retain the exact v0.3.0 archive reference |
| Figures and tables | Figure/table register; seven local SVGs, seven PDF derivatives and six correctly mapped inline review images regenerated and visually inspected | Ready for author export | Confirm third-party data permissions and final image format in the portal; no further export is currently required |
| Manuscript text | 165-word abstract with metric definitions, six keywords, 12 cited numeric references, separate AI declaration immediately before references and methods disclosure for code assistance | Ready | Final language proof after Word conversion |
| Highlights upload | Five bullet-only highlights, maximum 76 characters, with a standalone Word artifact verified by the submission audit | Ready | Upload `Highlights.docx` and retain the source text for the audit trail |
| Author declarations | Manuscript now has explicit Funding, CRediT, competing-interest and Acknowledgements fields; a prior PSEP package provides candidate metadata but is a different manuscript and still requires confirmation | Pending author input | Confirm the candidate against `docs/ijhe-author-metadata-source-2026-10-09.md` and complete `docs/ijhe-author-submission-intake.md`; do not infer omitted co-authors |
| Cover letter | Quantitative claims, scope boundary and preprint/prior-publication disclosure field are drafted; author/declaration placeholders remain | Pending author input | Replace placeholders and confirm exclusive submission and publication history |
| Word/PDF upload | Author-neutral DEGALI DOCX rebuilt with six inline figures, Word-exported to a 12-page PDF and visually checked page-by-page; strict export now also refuses unresolved declaration placeholders | Draft verified | Fill author metadata and declarations, then rerun strict export and final portal checks |
| Draft upload assembly | 35 hash-pinned upload and hand-off files assembled without third-party raw data in `outputs/IJHE_UPLOAD_DRAFT_LATEST_2026-10-09_DEGALI_ONLY_V2/`, including the portal-compatible graphical-abstract, standalone Highlights Word file and supplementary-information PDF/DOCX derivatives, the pinned research-runtime list, the local operational-data audit, the external-evidence search and byte-level selective-retrieval records, the three-event conditional PRESLHY manifest index and a single closeout worksheet; bundle status remains author-neutral | Ready as draft hand-off | Complete the closeout worksheet, replace placeholders, rerun strict export and rebuild the bundle after author edits |
| Matched transient field event | [Conditional PRESLHY E3.5 Trial 10, 20 and 21 bundles](ijhe-external-evidence-search-2026-10-09.md); D3.6 additionally supplies nominal sensor accuracy, receptor coordinates/heights and 0.16/0.18 m obstruction offsets; public H2SAFE/HyWAM records were also checked but are surrogate/context-only; the inspection-centre workspace has process PT/TT/DPT histories but no atmospheric/receptor package; field audits `candidate_complete`, manifests `conditional`, promotion `false` | Conditional | Obtain custodian-confirmed calibration, source-rate history, absolute clock and surveyed site geometry before strengthening the operational or obstacle claim |
| ELVHYS WP4.2 dynamic/geometry evidence | Public DOI archive with 48 confined tests, time-resolved H2/temperature/pressure channels, declared calibration specifications and measured TCS geometry; no measured H2 mass-flow history | Conditional boundary | Use for dynamic/geometry traceability only; resolve source boundary and enclosure/report geometry discrepancy before quantitative promotion |

## Verification record

The current local audit reports **75 pass, 0 fail and 3 warnings** after the
ELVHYS evidence-lane, cover-letter boundary and Word/PDF audit updates. The
warnings correspond exactly to the two author-finalisation rows and the
conditional matched-event gate. The latest
public snapshot contains 808 files and is 7.7 MiB. The figure-artifact audit
also confirms that all seven locally available SVG/XML artifacts have valid
dimensions and SHA-256 records.

Reference integrity was rechecked on 2026-10-09: all eleven DOI-bearing
references in the manuscript returned HTTP 200 from `doi.org`, and the FFI
report's DOI-free official publication page also returned HTTP 200. No DOI was
invented for the FFI report.

The software DOI chain was checked against Zenodo's public record API. The
previously reserved v0.2.0 DOI `10.5281/zenodo.23105451` is now published with
the exact `v0.2.0` artifacts. The current `v0.3.0` release is pinned to its own
reserved DOI `10.5281/zenodo.23256538`, while the resolving concept DOI remains
`10.5281/zenodo.22646258`. The submission audit now requires both the exact
release DOI and the stable concept DOI in public-facing metadata.

A repository-wide Markdown link scan then found 69 broken local targets in the
archived development log, all caused by an extra `docs/` prefix from that
file's own directory. The relative targets were corrected, and a repeat scan
of 342 Markdown files returned zero broken local links. The check is now also
covered by `tests/test_publication.py` and the publication checker itself.

The graphical-abstract audit now records the Elsevier general minimum of
1328 x 531 px explicitly. The current 1400 x 760 px SVG passes that minimum;
its preferred-aspect-ratio difference is retained as a portal confirmation
item rather than corrected by a potentially distorting resize.

The supplementary-information source was exported to a five-page DOCX and PDF
after adding a specification-only 2/4 vol% sensor-error table; all five PDF
pages were visually inspected, and the submission audit checks both the DOCX
container and the non-empty PDF derivative.

After clarifying the abstract's transfer-case sentence, the main DOCX and its
Word-exported 12-page PDF were regenerated. All 12 pages were inspected from
the latest PNG render; the bundled LibreOffice renderer was unavailable on the
host, so the controlled Word export plus Poppler rasterization was used as the
rendering fallback.

The manuscript-facing source was then cleaned of internal `Artifact:` file
paths and the upload-oriented subsection title was replaced with a
reviewer-facing decision-boundary/figure-caption heading. The DOCX/PDF were
regenerated again and all 12 pages were re-inspected; figure provenance remains
in the separate figure/table register and hash manifest. The abstract was also
updated to define MG, VG and FAC2 at first use; it remains within the checked
submission limit at 165 words.
The author-neutral title page now also records the research-article type and
the 4,247-word main-text count (abstract and keywords included; references and
declarations excluded).

The abstract and Section 5.9 now define the reported 11/18 decision changes as
binary 4 vol% classifications at a fixed 1.5 m receptor across two FFI cases,
three wind values and the 20/30/50 m lines. The main DOCX/PDF were regenerated
after this clarification and all 12 rendered pages were re-inspected; the
decision definition is visible without relying on the supplementary record.

The DOCX audit now directly fails if journal-facing text contains local
`docs/`, `tools/`, `outputs/` or `tmp/` path markers, including table cells.
The current export passes this check, so path hygiene no longer depends only
on the builder regression test.

The submission audit also now fails on accidental `TODO`, `FIXME` or `TBD`
production markers in the manuscript, supplementary information, highlights
or cover-letter source. The current journal-facing files contain none.

The latest language pass standardised the manuscript-facing prose to American
spelling (`centerline`, `labeled`, `modeling`) without changing cited titles or
headline results. The main DOCX/PDF were regenerated once more and all 12
rendered pages were re-inspected; the document audit remains warning-only for
the unresolved author fields.

The DOCX builder now maps repository-only `docs/` and `tools/` path labels to
reviewer-facing descriptions in the journal-facing DOCX/PDF. Reproducibility
links remain in the Markdown source and hash records, while the rendered
manuscript no longer exposes local repository paths. The supplementary DOCX/PDF
were rebuilt and all five rendered pages were re-inspected after the same
cleanup.

Figure 4 was also reconciled with its generated artifact: the SVG is a six-arc
FFI/DNV map covering Tests 4 and 6, while the separate Test 6 source/wind
decomposition remains a boundary table. The manuscript caption, figure
register and audit check now use the same description.

The latest non-slow regression command completed with **1,589 passed, 6
skipped and 19 deselected** in 1,232.33 s. This
confirms the current software and submission helpers in the checkout; neither
run closes the author or external-field evidence gates.

The paired-path uncertainty audit was regenerated after adding trial-cluster
VG intervals; the frozen seed, 5,000 draws and input hashes are retained in
`outputs/applied-energy-evidence-2026-10-09/manifest.json` and the derived
aggregate table. The update does not change any headline score or promotion
status.

The PRESLHY Trial 10 evidence manifest was also regenerated under the current
manifest schema. Its selected files verify successfully, its explicit
`sensor_calibration_status=specification_only` yields
`evidence_readiness=conditional`, and `promotion_allowed` remains `false`.
The submission audit now checks that gate directly, and the upload-bundle
builder refuses to copy a stale manifest that lacks the explicit calibration
status or has contradictory promotion flags.

The same public-data extraction and manifest verification was then replayed for
accepted Trials 20 and 21. Both manifests verify all selected files and retain
`evidence_readiness=conditional` with `promotion_allowed=false`; they are
workflow-reproducibility evidence, not additional headline validation scores.

The post-cleanup upload hand-off was rebuilt with 35 copied files and no
missing or hash-mismatched entries. Its machine-readable gate remains
`final_upload_allowed=false`, `preslhy_evidence_gate.status=conditional` and
`sensor_calibration_status=specification_only` until the author and custodian
requirements are closed.

The completion state is therefore **draft-ready with external validation and
author-finalisation pending**, not final-upload-ready. The Word/PDF artifact
row is now verified; the state should only be promoted after the remaining
author and evidence rows have authoritative evidence and the final Word/PDF
pages have been visually checked again after any author edits.

The public KITopen/RADAR PRESLHY E3.5 research-data record was rechecked on
2026-10-09 as a stronger custodian lead. Its metadata identifies 25 elevated
LH2 releases and names mass-flow, pressure, H2 concentration, temperature,
near/far-field weather, humidity, O2 and video channels; the archived package
is approximately 11.3 GB and is published under CC BY-SA 4.0
([KITopen record](https://publikationen.bibliothek.kit.edu/1000136281),
DOI `10.5445/IR/1000136281`; RADAR version DOI `10.35097/1481`). The landing
record does not itself provide channel-level calibration certificates, a
machine-readable obstacle survey or a compact matched-event manifest. The
selectively retrieved D3.6 report does provide nominal accuracy bounds,
source-relative sensor coordinates/heights and the two obstruction offsets,
but not the missing event-level calibration and surveyed 3-D geometry. It
therefore does not change the conditional gate or headline score. Selective
retrieval and custodian-confirmed provenance remain the highest-value external
action.

The individual public Trial 10 workbook was then selectively retrieved and
verified byte-for-byte against the local controlled copy (575,705 bytes;
SHA-256 `17987bc4e2f0a8d7bec40900063e3a3c2ce7a9ad6f9ab4014a5715ab04924aa`).
Its five sensor/time-series sheets confirm the source, H2, temperature and
weather channels but contain no separate calibration certificate, absolute UTC
clock, surveyed obstacle sheet or event-level geometry manifest. This is a
reproducibility confirmation of the conditional lane, not a promotion to
certified field validation. The same check now covers Trials 20 and 21 in
`outputs/preslhy-public-retrieval-index-2026-10-09.json`.

The manuscript and supplementary information now state the complete E3.5
selection rule and the 24-workbook population explicitly. The main DOCX/PDF
and five-page supplementary DOCX/PDF were regenerated after that edit; all 12
main pages and all five supplementary pages were raster-rendered and visually
checked. The focused IJHE regression suite remains 10 passed with `git diff
--check` clean.
