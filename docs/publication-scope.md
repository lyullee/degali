# GitHub publication scope

## Public repository

The GitHub repository contains the DEGALI Python source, self-contained tests,
research tools, physical derivations, data provenance, validation methods and
aggregate validation results needed by software users.

The public snapshot is intentionally an alpha research release. The stable
claim is source-level DEGADIS 2.1 reproduction. Liquid-hydrogen extensions are
opt-in and carry their own validation and applicability statements.

## Kept local

The following remain on the workstation and are excluded by `.gitignore` and
the source-distribution manifest:

- the original DEGADIS Fortran source and executable oracle;
- all external experimental data, including reduced or selected copies;
- third-party journal and technical-report PDFs;
- solver checkpoints, partial JSON files and large numerical field dumps;
- virtual environments, caches, compiled Fortran objects and scratch clones;
- machine-specific regression XML and temporary probes.

Nothing is deleted by this policy. Large local results remain traceable through
the dated documentation and can be archived separately in an access-controlled
or DOI-backed data deposit if their licences permit it.

## Before the first GitHub release

1. Choose the final GitHub owner and verify the URLs in `pyproject.toml` and
   `CITATION.cff`.
2. Replace the author placeholder in `CITATION.cff` with the desired public
   name, ORCID and affiliation.
3. Decide whether the first tag remains a development snapshot or becomes an
   alpha release; keep `pyproject.toml` and `CITATION.cff` versions identical.
4. Run `python tools/check_publication.py`, build wheel/sdist, and inspect both.
5. Run the self-contained Python suite in GitHub Actions and retain the full
   Fortran-oracle/experimental-data audit as a controlled local release record.
6. Reserve a Zenodo DOI before inserting a DOI and release date, then archive
   the exact Git tag rather than a later working tree.

## Claims that must remain visible

- Conservation and source-oracle agreement are numerical verification, not
  experimental validation.
- The main LH2 quantitative assessment is PRESLHY-heavy and does not cover all
  release orientations, source conditions or facility geometries.
- Finite-TKE, pressure/stress and independent thermal-width work is not yet a
  validated downstream default.
- Safety-critical use requires independent review and comparison with accepted
  engineering practice.
