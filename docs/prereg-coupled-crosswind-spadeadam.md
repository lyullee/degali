# Pre-registration: coupled phase-profile path on Spadeadam 4 and 6

Date frozen: 2026-09-04

## Purpose

Test the newly conserved near-field/crosswind handoff on independent outdoor
field concentrations. No Spadeadam value is used to form the handoff or tune a
coefficient.

The published Spadeadam flow/orifice/P04 data identify a two-phase atmospheric
H2 state, not the single-phase gas plane required by the axisymmetric model.
The candidate therefore uses the already implemented, previously
pre-registered 1 micrometre transported-condensed-air source with total source
energy conservation and no pressure-thrust option. This source choice predates
the coupled handoff and was already found not to improve the old JETPLU path;
it is not selected from the new result.

The 1 micrometre particle size remains a model-form sensitivity, not a measured
Spadeadam input and not a proposed facility default.

## Fixed cases and conditions

- Outdoor horizontal tests 4 and 6 only.
- Published H2 mass flow, 25.4 mm nominal orifice and P04 pressure.
- Low-mast wind: 5.0 m/s for test 4 and 2.3 m/s for test 6.
- Release height 0.50 m.
- Ambient 277.15 K and 90% RH, following the DNV modelling assumption already
  used by the existing comparison.
- Neutral stability, 1 mm roughness, 10 m wind reference.
- Source position begins at the transported source-zone distance.
- Conserved axisymmetric integration ends at `S = 3.0 m`, using 41 radial
  points, 0.01 m maximum step and `2e-6` relative tolerance.
- Default phase-profile handoff; every interface screen must pass before a
  trajectory is scored.
- JETPLU integrates to at least 250 m without coefficient changes.

## Observations and scoring

Use the existing extracted sensor table. At 30, 50 and 100 m, compare the
maximum predicted concentration over the same three sensor heights with the
maximum reported peak over bearing at each height. This retains the current
six-arc lower-censored comparison.

Report Hanna geometric mean bias (MG), geometric variance (VG), FAC2, the six
observed/predicted maxima, centre heights and vertical sigmas. The reference
corrected-default result is MG/VG/FAC2 = `0.869/1.476/0.667`.

## Frozen decision rule

The candidate passes this pilot only if:

1. both handoffs pass their independently frozen conservation, energy, width,
   temperature and quadrature screens;
2. all six arcs are reached;
3. MG is closer to one than 0.869;
4. VG is no greater than 1.476;
5. FAC2 is at least 0.667; and
6. no non-physical state, solver failure or forced-ground workaround occurs.

If it passes, repeat at 81 radial points and 0.0025 m maximum step before
quoting the field result. If it fails, retain the coupled interface for
measured single-phase sources but do not promote this uncertain two-phase
source reconstruction.

## Results after pre-registration

The 1 micrometre total-energy source endpoints were:

| test | T (K) | velocity (m/s) | diameter (m) | density (kg/m3) | H2 mass fraction | source-zone x (m) |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 68.745 | 30.680 | 0.3537 | 1.4700 | 0.18689 | 0.7970 |
| 6 | 68.735 | 35.180 | 0.3310 | 1.4638 | 0.18799 | 0.7446 |

The original instantaneous phase-profile table hit its concentration ceiling
and failed the native balances for both sources. The separately
pre-registered phase-manifold domain extension then gave:

| test | source velocity/local wind | max native residual | energy | H2 width | centre T | interface |
|---|---:|---:|---:|---:|---:|---|
| 4 | 9.091 | 1.68e-16 | 0.606% | 6.94% | 5.257 K | fail |
| 6 | 22.662 | 1.64e-16 | 0.180% | 1.42% | 0.885 K | pass |

Test 4 fails the fixed width and temperature screens and is already below the
independent momentum-dominated applicability ratio of 10. Because both
handoffs were required to pass before any sensor prediction, no arc
concentration or MG/VG/FAC2 was calculated.

## Decision

Reject this reconstructed-source pilot. The phase-manifold extension is useful
for diagnosing the representation domain and passes test 6, but it cannot turn
a wind-steered two-phase test-4 source into an axisymmetric near field. Do not
weaken the interface or velocity-ratio thresholds and do not quote a partial
three-arc field score.
