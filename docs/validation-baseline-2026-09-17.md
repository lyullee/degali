# Validation baseline — 2026-09-17

This document fixes the reproducible starting point for the next downstream
thermal/density and pool-to-plume investigations.  It is a regression and
evidence-scope record, not a new accuracy claim.

## Regression result

Command:

```powershell
$env:PYTHONPATH='src'
.\.venv\Scripts\python.exe -m pytest -q -m "not slow"
```

Result: **1049 passed, 127 skipped, 58 deselected** in 686.97 s.
The deselected tests are explicitly marked `slow`; skipped tests are
dependency- or non-distributed-data-gated.  `git diff --check` was clean.

Two independent invariants that initially failed in this working tree were
also rerun after correction:

- the pipe → flash → phase-plane → Gaussian source energy ledger for normal
  and para hydrogen, with both ambient-closure variants;
- exact replay of the frozen PRESLHY trial-23 downstream states, fluxes and
  sensor arcs.

The former correction preserves the actual post-flash specific enthalpy for a
subcooled liquid plane; reconstructing it from a saturated quality had
discarded sensible enthalpy.  The latter does not change the present water-ice
property model.  It selects the historical water-ice table only when replaying
an unversioned, already frozen 2025-09 field, so that archived states are not
silently reinterpreted by later physics changes.

## Validation-portfolio result

`tools/audit_validation_portfolio.py` reported:

- automatic default promotion: **not allowed**;
- composite accuracy score: **not calculated**;
- third-party observations distributed by this repository: **no**.

The qualified field-screen branches remain `fast_single_velocity_lh2`
(PRESLHY E3.5 and the independent reported Spadeadam arcs) and
`cryocompressed_source` (PRESLHY E3.1 pressure history).  They are qualified
screens, not field-accuracy or automatic-default claims.  The legacy
reconstruction, cryogenic-gas axisymmetric, finite-TKE, finite-rate droplet
and equilibrium-air-condensation branches remain reconstruction-only or
research-only according to their documented evidence limits.

## Consequence for the next work

The baseline supports diagnostic comparison of downstream thermal and density
profiles.  It does **not** authorize fitting a new closure, promoting an
unvalidated turbulence/slip branch, or claiming obstacle-wake prediction.
Any proposed physics change must preserve this regression set and be accepted
only against independent observations using a pre-declared observation
operator.

## Current-source confirmation — 2026-09-19

After the E3.4 substrate-to-pool heat-balance source path, the E3.5
measurement-clock alignment path, and the time-resolved pool-source interface
were added, the entire *current* test collection was rerun with a clean
pytest cache:

```powershell
.\.venv\Scripts\python.exe -m pytest --cache-clear -q
```

Result: **1105 passed, 139 skipped** in **1884.65 s**.  The skipped cases are
explicit external-reference or optional-dependency gates; no test was
deselected.  `git diff --check` was clean.

This confirms code consistency, not a new field-accuracy score.  In
particular, the downward-impingement Spadeadam screen is now explicitly a
bounded diagnostic: using the physically relevant Gaussian axis concentration
rather than a conserved section mean reveals an approximately 20-percent
dependence on the unresolved impingement footprint.  It must not be promoted
as a tuned quantitative validation branch.

## Scope-closure regression — 2026-09-19

After the final E3.4 trial-identity gate and the public hydrogen-channel
timing intake were added, the complete current collection was rerun again
with a clean cache:

```powershell
.\.venv\Scripts\python.exe -m pytest --cache-clear -q
```

Result: **1108 passed, 139 skipped** in **2453.13 s**.  No test was
deselected and `git diff --check` was clean.

This baseline closes four claim boundaries without adding a fitted
coefficient: (1) a D3.5 E3.4 report example cannot be numerically compared
against a workbook with a different trial identifier; (2) the public FFI
hydrogen-channel record may audit separately clocked source/sensor timing but
cannot validate an atmospheric LH2 transient plume; (3) co-condensed N2/O2
states remain a flagged mixture-EOS boundary; and (4) cuboid/wall contact
remains an invalidation screen, not a quantitative LH2 wake prediction.
