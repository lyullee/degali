# Direct-H₂ cylinder-impingement benchmark audit (2026-10-13)

## Decision

The 2025 unignited-H₂ cylinder-impingement experiment is a useful direct-H₂
obstacle physics target, but it is **not yet an executable quantitative DEGALI
validation case**. The public article record supports the test envelope and
qualitative concentration trends. It does not expose a machine-readable
concentration table, a time history, or enough boundary data to reconstruct a
unique transient source. Keep it in the conditional evidence queue; do not
digitize an abstract claim into a validation score or call it a CFD comparison.

## What the primary article record establishes

The experiment concerns an unignited hydrogen release impinging on a cylinder.
It reports a 1.5 mm nozzle, release pressures from 0.15 to 0.30 MPa, and
nozzle-to-cylinder distances from 0.20 to 0.40 m. The bottom of the cylinder
was approximately 0.16 m above the floor. The article describes concentration
decay/profile measurements after impingement and reports trends with pressure
and spacing; a longitudinal profile is described as Gaussian in a
dimensionless hydrogen mass fraction, with a separate stagnation-point vertical
correlation. This is direct-H₂ obstacle evidence and is physically closer to
the sought impingement mechanism than neutral-gas obstacle tests.

Primary source: Yao et al., *International Journal of Hydrogen Energy* 100
(2025), 388–397, DOI [10.1016/j.ijhydene.2024.12.307](https://doi.org/10.1016/j.ijhydene.2024.12.307),
[ScienceDirect article record](https://www.sciencedirect.com/science/article/abs/pii/S0360319924055204).
The public indexed record was readable during this audit; opening the article
page itself returned HTTP 403. No article PDF, supplementary file, or dataset
was located in the scoped local evidence package or the targeted searches for
the PII/DOI and data-availability terms. This is a search result, not proof
that no author-held or subscription data exist.

## Missing evidence before a score can be computed

| Requirement | Current state | Consequence |
|---|---|---|
| Exact cylinder diameter, height, orientation, and full installation geometry | Not confirmed from the accessible primary record | Do not infer geometry from a related ignited-jet paper; it is a different experiment |
| Sensor count, coordinates, response/calibration, and sampling cadence | Not registered from accessible primary record | No defensible receptor operator or uncertainty weighting |
| Raw concentration measurements and synchronized time histories | Not found | No sensor-history score, arrival-time/dose score, or measured contour comparison |
| Complete source boundary (mass-flow history, reservoir state/temperature, discharge coefficient) | Pressure range/nozzle size reported; source history not found | Pressure cannot be silently substituted for measured mass flow; transient score withheld |
| Ambient conditions, wind/confinement, and test-by-test conditions | Incomplete | Cannot set reproducible outer boundaries or attribute mismatch |
| Machine-readable observation tables / uncertainty | Not found | Plot digitization, if later performed, must be labeled figure-derived and uncertainty-tagged |

## Appropriate use if fuller data become available

This case should be a separate unignited, near-field H₂ impinging-cylinder
benchmark, not pooled with Thorney Island or FFI. Start with one nominally
complete run only after source and geometry registration. Preserve measured
vol% as vol%; if the paper provides mass fraction instead, convert explicitly
using the H₂/air mixture molecular weights and document the paper's exact
definition. Compare wall/stagnation/side profiles at the reported stations,
then conduct distance and pressure sensitivity without tuning to the scored
profiles. A steady or time-averaged profile comparison must remain distinct
from transient validation. Report peak and profile error only with pointwise
observations and their uncertainty; report the 4 vol% exceedance region only
if spatial support and a defensible interpolation rule are available.

The experiment does **not** by itself validate cryogenic thermodynamics: its
accessible abstract does not identify a cryogenic release state. Do not merge
its results with the 2023 cryogenic vertical-wall experiment without a separate
case manifest and verified source conditions.

## Next acquisition action

Seek the full text and any supplementary files via the authors' institutional
repository or a legitimately available library copy, then inspect the actual
figures/tables and data-availability statement. If no raw data are released,
requesting them from the corresponding author would require explicit user
authorization. Until then, the existing FFI bundle remains the best
machine-readable direct-H₂ obstacle time-series candidate, although its
cross-beam geometry and source-clock coverage gates remain unresolved.
