# IJHE author-submission intake

This form is the final human-input gate for the DEGALI IJHE package. It is
deliberately separate from the scientific manuscript so that author metadata
cannot be inferred from the software repository or copied from a separate
submission package.

For a single hand-back sheet covering both author declarations and the matched
event evidence gate, use [`docs/ijhe-closeout-worksheet.md`](ijhe-closeout-worksheet.md).

## Required before upload

## Candidate identity found in the repository

The public `CITATION.cff` and `pyproject.toml` name **Ugwiyeon Lee** with
**Korea Gas Safety Corporation** as an affiliation. A prior, separate PSEP
upload package also contains a candidate four-author list, ORCIDs, funding,
CRediT and competing-interest wording. These are only candidate inputs for
author confirmation; they are not copied into the manuscript or cover letter
automatically. The source and exact confirmation gate are recorded in
[`docs/ijhe-author-metadata-source-2026-10-09.md`](ijhe-author-metadata-source-2026-10-09.md).

- **Manuscript title:** Observation-operator-aware validation and decision
  implications of a fast liquid-hydrogen dispersion model
- **Article type:** Research article (confirm in the submission system)
- **Corresponding author:** full name, institutional email, postal address,
  telephone if required by the system, and ORCID
- **Complete author list:** spelling, order, affiliations, ORCID, and approval
  from every author; compare the candidate source note but do not accept it
  without explicit confirmation
- **CRediT contributions:** conceptualization, methodology, software,
  validation, formal analysis, investigation, data curation, visualization,
  writing, supervision, project administration, and funding acquisition as
  applicable to each named author
- **Funding:** funder name, grant number, and whether the funder had a role in
  study design, data collection, analysis, writing, or the decision to submit
- **Competing interests:** the exact journal-form wording, including a
  statement if none exist
- **Data permissions:** confirmation that every third-party report, workbook,
  sensor record, and derived table is used under its permitted terms
- **Exclusive submission:** confirmation that the manuscript is not under
  consideration elsewhere
- **Preprint/prior-publication disclosure:** identify any preprint, thesis,
  conference version or overlapping publication, or explicitly confirm that
  none exists and that the submitted work is not previously published
- **Suggested/excluded reviewers:** optional; provide only after conflict and
  institutional-affiliation checks

## Final file gate

After the fields above are confirmed, replace the placeholders in
`IJHE_COVER_LETTER_DRAFT.md` and the declarations at the end of
`MANUSCRIPT_LH2_VALIDATION.md`. Then regenerate the Word/PDF submission files
from the canonical Markdown sources and perform page-by-page visual QA. The
legacy Word files supplied from the separate workspace are not valid substitutes;
see `docs/ijhe-word-artifact-boundary.md`.

## Scientific claims to confirm at upload

The corresponding author should explicitly approve that the submitted version:

1. keeps E3.5, FFI/DNV and Hecht–Panda evidence in separate lanes and excludes
   cross-code performance ranking from the submission;
2. reports the Test 6 low-wind residual rather than fitting it away;
3. labels the 3-D transient result as numerical verification, not measured
   validation; and
4. does not present the model as a regulatory design basis, obstacle-resolved
   CFD replacement, or certified safety-distance tool.
