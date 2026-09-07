# Conserved LH2 near-field to crosswind handoff: results

Date: 2026-09-04

## Outcome

The accepted conserved axisymmetric LH2 near field can now be connected to
JETPLU without silently copying incompatible centreline variables. The
interface projects total mass, hydrogen mass, horizontal momentum and vertical
momentum as integral constraints, audits total energy and profile width, and
refuses downstream integration on any failed screen.

The accepted interface is the pre-registered **phase-profile** candidate. It
transfers the near-field equilibrium N2/O2/H2O radial thermodynamics and
component temperature-dependent enthalpy into JETPLU's concentration table.
No entrainment or dispersion coefficient was fitted.

## Representative calculation

The numerical screen used the Hecht--Panda 5 bar-class measured throat:

- throat diameter 1.0 mm,
- absolute throat pressure 2.422 bar,
- throat temperature 37.4 K,
- throat density 1.65 kg/m3,
- throat velocity 498.2 m/s,
- horizontal release at 0.5 m,
- 2.5 m/s wind at 10 m, neutral stability and 1 mm roughness,
- resulting local wind 1.687380257 m/s, supplied consistently as near-field
  co-flow,
- dry air at 295 K and 101325 Pa.

At the first tested phase-profile station, `S = 0.080 m`, the frozen audit was:

| quantity | relative mismatch | limit | result |
|---|---:|---:|---|
| total mass | 4.18e-16 | 1e-8 | pass |
| hydrogen mass | 5.04e-16 | 1e-8 | pass |
| horizontal momentum | 1.67e-16 | 1e-8 | pass |
| vertical momentum | 1.02e-20 | 1e-8 on a 1 N scale | pass |
| total energy | 5.85e-4 (0.0585%) | 2% | pass |
| H2 mass-density half-width | 1.571e-2 (1.57%) | 5% | pass |
| centre temperature | 1.045 K | 2 K | pass |
| doubled-order energy quadrature | 1.25e-6 | 1e-5 | pass |

The accepted state then integrated through the crosswind equations past
2.08 m in the smoke test without a non-physical state or solver failure.

A separately pre-registered strict-grid repeat (81 radial points, 0.25 mm
maximum step, `5e-8` tolerance) retained the pass: 0.05853% energy, 1.571%
H2 half-width, 0.943 K centre temperature and `1.916e-6` doubled-order energy
quadrature residual. See `prereg-handoff-grid-convergence.md`.

Handoff-position uncertainty was then isolated. The instantaneous
phase-profile table passes at 0.08 and 0.10 m but reaches its table ceiling at
0.13 m. The separately pre-registered phase-manifold extension passes all
three and changes 0.5--2 m centreline concentration by at most 0.305%, sigmas
by at most 0.347%, and centre height by 0.000307 m. The public measured-source
default remains the strictly confirmed 0.08 m phase-profile boundary; see
`prereg-handoff-location-invariance.md`.

## Falsified alternatives retained in code

The original centre-state dilution line fails at 0.025 m: energy 13.46%, H2
half-width 5.05%, and centre temperature 3.974 K. Its four native balance
residuals are nevertheless below `1e-15`, showing why checking only mass and
momentum would have admitted a thermodynamically discontinuous handoff.

The flux-averaged energy line makes energy an exact fifth constraint. At
0.080 m it closes energy to 1.26e-11 and width to 0.958%, but centre
temperature still differs by 2.823 K and it is rejected without relaxing the
2 K criterion.

Both candidates remain selectable by name for reproducibility, but a rejected
handoff cannot call the downstream integrator.

## Public research path

`run_lh2_crosswind_research` now assembles one consistent path from a measured
or ambient-pressure source. It calculates the local log-law wind first,
supplies that velocity to the conserved near field, performs the audited
phase-profile handoff, and returns an object whose `run()` method is enabled
only after acceptance.

The current coupled path deliberately accepts only a horizontal release
aligned with the mean wind. An angled or cross-axis release needs a
three-dimensional near-field formulation; pretending that transverse wind is
axisymmetric would reintroduce the momentum defect the handoff removes.

## What this does not establish

This is an interface-verification result, not independent dispersion
validation. The next stage is to rerun eligible PRESLHY/Spadeadam cases with
the coupled path and score concentrations, vertical spread, centre height and
ground contact without retuning. Ground heat transfer and atmospheric
unsteadiness also remain separate uncertainty tasks.

The two public Spadeadam horizontal sources cannot yet be used directly for
that test. Their published flow, nominal orifice and P04 pressure produce
conserved HEM atmospheric states with H2 vapour qualities only 0.1007 (test 4)
and 0.1142 (test 6), at 20.369 K and densities 11.33/10.18 kg/m3. They are
two-phase jets, whereas this coupled path starts from a single-phase gas
plane. The public runner now rejects that misuse explicitly.

To connect these trials without an assumption, the required boundary is the
first atmospheric-pressure single-phase station after flash/evaporation:
effective diameter, bulk velocity, temperature, density, hydrogen mass
fraction, and preferably separate residual liquid-H2 and condensed-air mass
flows. A validated two-phase source-zone calculation providing the same
quantities would also close the gap. The existing P04 pressure and total flow
alone do not identify them.

A pre-registered attempt with the existing 1 micrometre transported-air,
total-energy source reconstruction was stopped before sensor scoring. A
phase-manifold table made test 6 interface-compatible, but test 4 still failed
the fixed width (6.94%) and centre-temperature (5.257 K) screens and had a
source-velocity/local-wind ratio of only 9.091, below the independently fixed
momentum-dominated limit of 10. See
`prereg-coupled-crosswind-spadeadam.md`. This narrows the blocker: test 4 needs
a crosswind-aware three-dimensional/two-phase near field or measured
single-phase outlet data, not a looser handoff tolerance.
