# IJHE closeout worksheet

This is the single hand-back sheet for the two remaining IJHE gates. It is
not an authority for inventing metadata or filling missing measurements. The
corresponding author or accountable data custodian should complete the fields,
attach the supporting evidence, and return the signed version to the package
owner. Blank or uncertain fields remain open in the automated audit.

## A. Author and submission metadata

| Field | Confirmed value | Evidence or approver |
|---|---|---|
| Manuscript title |  |  |
| Article type |  |  |
| Complete author list and order |  |  |
| Affiliations for each author |  |  |
| ORCID for each author |  |  |
| Corresponding author, email, address |  |  |
| All-author approval obtained |  |  |
| Exclusive submission confirmed |  |  |
| Preprint/thesis/conference overlap or “none” |  |  |
| Suggested/excluded reviewers (optional) |  |  |

## B. Declarations to copy into the manuscript

| Declaration | Final text or decision | Approver/date |
|---|---|---|
| Funding and grant number |  |  |
| Funder role in design, analysis and submission |  |  |
| CRediT contribution for every named author |  |  |
| Competing-interest statement |  |  |
| Acknowledgements or “None” |  |  |
| Generative-AI disclosure approved for this paper |  |  |
| Third-party data and figure permissions |  |  |

Candidate values in `docs/ijhe-author-metadata-source-2026-10-09.md` are only
provenance for confirmation. They must not be copied automatically from the
prior PSEP submission.

## C. Matched-event evidence gate

Complete one row per event. A single event is sufficient for a qualified
case-study lane only when every required channel is supplied and reviewed.

| Channel | File/record supplied | SHA-256 | Custodian / status |
|---|---|---|---|
| Stable `event_id` and common-clock record |  |  |  |
| Source location, phase, geometry and mass-flow history |  |  |  |
| Weather speed, direction, height, stability and station |  |  |  |
| Measured obstacle/site geometry and coordinate reference |  |  |  |
| Receptor coordinates and H₂ concentration time series |  |  |  |
| Sensor registry, calibration certificate and response metadata |  |  |  |
| Operator/data-custodian identity and processing version |  |  |  |

The minimum machine-readable layout is described in
`docs/ijhe-field-evidence-request-2026-10-09.md`. A nominal instrument
specification, a numerical convergence study or a nearby process historian
does not substitute for an event-linked calibration and source/weather/
receptor record. Until this section is complete, the event remains
`conditional` or `withheld` and the manuscript must retain its current
operational-claim boundary.

## D. Final sign-off

- [ ] I confirm the author list, order, affiliations and ORCIDs for this exact
  DEGALI manuscript.
- [ ] I confirm the declarations, permissions, exclusive-submission status and
  prior-publication disclosure.
- [ ] I confirm that every promoted event has an accountable custodian,
  common-clock evidence and SHA-256-pinned source files.
- [ ] I understand that the current paper does not claim regulatory separation
  distances, obstacle-resolved CFD replacement or certified alarm thresholds.

**Name / role:**
**Signature or accountable approval record:**
**Date (UTC):**
