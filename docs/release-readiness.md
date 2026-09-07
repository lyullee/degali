# Release readiness

Checks applied before publishing, taken from a parallel reimplementation of
SLAB that had already been through the process and recorded what went wrong.
Each is now a test, so it is checked on every run rather than once.

## Patterns that hide errors

**A bare `raise` outside an `except` block.** The SLAB audit found one that had
drifted below the handler that gave it meaning: a genuine property-lookup
failure surfaced as `RuntimeError: No active exception to reraise`, with the
original error and its message gone. This package has none, and a test walks
the AST to keep it that way.

**`except Exception: pass`.** One site did this — the vapour enthalpy lookup,
which needs a fallback above the critical point where there is no saturation
line. It also swallowed an unknown fluid name, returning a plausible number
for a substance that does not exist. It now distinguishes the two: a missing
critical temperature raises, a missing saturation line falls back.

## What ships

| check | state | enforced by |
|---|---|---|
| third-party data in the package | none | `test_the_package_carries_no_third_party_data` |
| suite without any external data | passes; the field tests skip | CI, with no data configured |
| absolute paths in the package *and the tests* | none | `test_no_absolute_paths_in_the_package` |
| superseded figures in the prose | none | `test_no_retired_figure_survives_in_the_prose` |
| printed evidence matches the suite's | yes | `test_the_reported_evidence_is_the_measured_evidence` |
| reference Fortran or external experimental data in the sdist | none | `pyproject.toml`, publication check |

**Counts are not recorded here.** They were, in four documents, and they had
drifted in all four — 103, 167, 168/19 and 192 against an actual 202. A number
in prose has nothing checking it. What is recorded instead is the *property*,
and a test that holds it.

REDIPHEM, the SMEDIS spreadsheets and the PRESLHY workbooks are all someone
else's to distribute. None ships. Every one of them is located from an
environment variable — `REDIPHEM_ROOT`, `DEGALI_E35_ROOT`,
`DEGALI_E35_REPORT`, `DEGALI_SMEDIS_ROOT` — and the tests skip without it,
so a clone with no data still passes.

That was not true when this document first claimed it. The REDIPHEM tests read
a variable; the PRESLHY and SMEDIS tests held paths from the machine this was
written on, and the guard that was supposed to catch that only walked `src/`.
The suite passed on a clone by skipping for the wrong reason, which is the
failure mode this whole document is about.

The original Fortran and every external experimental-data copy remain local.
The public documentation records provenance, methods and aggregate results;
authorized reviewers can rerun the controlled local validation separately.

## Regression tiers

Use `pytest -m "not slow" -q` while developing. Whole-model LH2 assessments,
campaign reductions, strict multiphase integrations and full comparison runs
carry `@pytest.mark.slow` and run with `pytest -m slow -q`. This split is by
physics cost, not importance: both tiers are required before release. The
complete source-built Fortran oracle tier remains a controlled local release
check and is not part of public CI.

## Versioning, from someone else's mistake

The SLAB release history contains a version that must not be cited: its
archived README and `CITATION.cff` point at the *previous* version's DOI rather
than their own. DEGALI uses Zenodo's GitHub integration, which assigns a DOI
only after the first GitHub release. The first archive therefore carries the
correct title, author, version, release date and repository URL without a DOI;
the newly assigned DOI is added in the first post-release metadata update.
Never insert a DOI belonging to another project or version.

That is worth knowing before the first release rather than after it.
