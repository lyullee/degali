# Distribution and commit plan — 2026-09-19

This is a release-preparation record, not a release.  No commit, tag, push,
upload or DOI metadata change is made by this plan.

## Current local gate

The complete clean-cache suite completed on 2026-09-19 with **1108 passed,
139 skipped**.  The skipped cases are deliberately gated on optional
dependencies or non-distributed reference material.  `git diff --check` was
clean.

`python -m build --wheel --sdist` built:

- `degali-0.1.0-py3-none-any.whl`
- `degali-0.1.0.tar.gz`

Both artifacts were inspected by filename.  The combined 730 archive entries
contain no `reference/` path and no `.for`, `.f90`, `.f`, `.xlsx`, `.xls`,
`.zip` or `.pdf` payload.  Build products stay in ignored `dist/` and are not
part of a commit.

## Commit boundary

Stage only source, tests, tools and documentation that explain the public
model and its evidence boundaries.  Before staging, review every currently
modified or untracked path against these exclusions:

1. Original DEGADIS Fortran and any reference executable/source stay under
   ignored `reference/` and never enter Git.
2. PRESLHY, SMEDIS, REDIPHEM, DNV, Sandia, FFI/AIJ and any other third-party
   raw measurement, image, spreadsheet, archive or digitised derivative stay
   outside the repository.
3. `dist/`, virtual environments, caches, local manifests and generated audit
   JSON stay untracked.
4. A provenance link and the observation method may be committed; a copied
   measurement value or calibration series may not.

## Proposed release sequence after review

1. Inspect the staged diff and confirm the above exclusions.
2. Commit the reviewed public source/tests/docs as one evidence-boundary
   change; do not mix it with unrelated local edits.
3. Push only after the complete suite is still green on the intended commit.
4. Create a GitHub release tag for `0.1.0`, then allow the configured PyPI and
   Zenodo workflows to run.  Check workflow logs before treating either
   publication as complete.
5. After Zenodo assigns a DOI, add the version-specific DOI only in a follow-up
   metadata/documentation commit.  Never borrow another repository's DOI.

## Scientific release statement

The public package provides a verified DEGADIS reimplementation plus
research-grade LH2 source/dispersion extensions.  E3.4 trial identity,
time-resolved source-to-plume timing, low-temperature N2/O2 condensed-phase
scope, and obstacle/wall contact are explicitly bounded.  They are not
presented as fitted field-accuracy, atmospheric transient-puff, mixture-EOS,
or quantitative wake closures.
