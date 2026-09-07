# Contributing

DEGALI combines a compatibility port with opt-in research models. Changes
must keep those claims separate.

1. Create a focused branch and add tests for every numerical or physical
   change.
2. Run `python -m pytest -m "not slow" -q` during development.
3. Run the full suite, including the Linux `gfortran` oracle, before a release.
4. Do not tune a physical coefficient against a validation result without a
   pre-registered hypothesis, an untouched comparison set and documented
   applicability limits.
5. Do not commit raw third-party datasets, literature PDFs, credentials,
   machine-specific paths, checkpoints or generated bulk fields.
6. Mark research-only candidates explicitly; passing a conservation test is
   not evidence of predictive accuracy.

Bug reports should include the input deck or a minimal synthetic case, Python
version, operating system, selected backend/preset, full traceback and the
expected physical behaviour. Remove confidential plant details before posting.
