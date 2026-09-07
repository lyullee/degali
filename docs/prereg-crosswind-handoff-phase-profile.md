# Pre-registration addendum: phase-profile thermodynamic handoff

Date frozen: 2026-09-04

## Reason for the third candidate

The centre-state candidate failed the frozen energy, width, and temperature
screens. The subsequent flux-energy candidate closed energy and all four
native JETPLU balances, and reduced width mismatch below 1%, but its smallest
centre-temperature mismatch over 0.08--0.13 m was 2.06 K, just outside the
unchanged 2 K screen.

Inspection identified a model-form discontinuity: the accepted near-field
model includes equilibrium N2/O2/H2O phase change and component
temperature-dependent enthalpy, while JETPLU's hydrogen table includes only
the historical air/water mixture rules. The missing phase/caloric physics,
not an entrainment parameter, is the third candidate.

## Candidate

At the selected near-field station, sample the complete conserved radial
profile already calculated by the near-field model:

- hydrogen mass concentration,
- mixture density,
- hydrogen/dry-air/water composition,
- equilibrium temperature, and
- specific enthalpy relative to ambient.

Reverse and deduplicate those samples into the monotonically increasing
concentration table required by JETPLU, with an exact ambient row prepended.
This transfers the N2/O2/H2O phase equilibrium and component `h(T)` model
without fitting a coefficient or inventing a virtual centre temperature.

Project total mass, hydrogen, horizontal momentum, and vertical momentum as
before. Because enthalpy is no longer a linear function of hydrogen mass
fraction, evaluate its cross-sectional flux by deterministic Gaussian
quadrature. The kinetic-energy part uses the analytic Gaussian factors.

## Frozen acceptance conditions

The original, unchanged interface screens apply:

- each native balance residual at most `1e-8`,
- total-energy mismatch at most 2%,
- hydrogen mass-density half-width mismatch at most 5%,
- projected centre-temperature mismatch at most 2 K,
- energy quadrature change on doubling the order at most `1e-5` relative,
- finite positive state inside the transferred table.

The phase-profile candidate is evaluated over the already computed
0.08--0.13 m interval. The earliest passing station is reported, but the
public default remains the caller's explicit station or the final computed
near-field station; no downwind experimental result selects the boundary.

## Interpretation

Passing this screen permits crosswind integration but is not validation of
that integration. PRESLHY/FLADIS concentration and ground-contact comparisons
remain a separate independent stage.
