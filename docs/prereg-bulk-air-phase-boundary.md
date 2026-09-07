# Pre-registration: bulk-air phase-safe jet boundary

Written before running the PRESLHY or Spadeadam comparisons with the new
source option.

## Fault under test

The current flashing source stops when the last liquid hydrogen evaporates.
For the PRESLHY pressure range that endpoint is about 20.0 K and contains
1.07--1.55 kg of entrained dry air per kg of hydrogen.  The downstream table
then treats that dry-air mass as an ideal gas.  Nitrogen and oxygen cannot
both remain gaseous at their resulting partial pressures, so the source state
is outside the thermodynamic model's stated domain.

Sandia's HyRAM+ technical reference places an initial plug-flow entrainment
and heating zone before its Gaussian jet.  It conserves mass, momentum, and
energy until a minimum usable temperature is reached.  The current Sandia
implementation leaves the minimum user-selected and disabled by default.

This test removes that arbitrary temperature.  It advances the adiabatic
plug-flow source until the three bulk constituents of dry air (N2, O2, and
Ar) can all be gaseous at their own mixture partial pressures.  At each
candidate temperature the entrained-air ratio follows from the component
enthalpy balance.  The handoff temperature is the first point at which

```text
max_i(p_i / p_sat,i) = 1,  i in {N2, O2, Ar}.
```

Ambient water and trace CO2 are deliberately excluded from that criterion.
DEGADIS already carries water condensate separately, while CO2 is below
0.1 % of dry-air mass.  No condensation or slip coefficient is fitted.  The
boundary is the no-loss/no-slip limit: condensate may form inside the skipped
zone, but all bulk dry-air condensate has re-evaporated by its exit.

## Predictions fixed in advance

- The 1--8 barg hydrogen endpoint will move from about 20 K to roughly
  65--75 K.
- The required dry-air ratio will increase and the hydrogen mass fraction
  will fall.
- N2, O2, and Ar saturation ratios will be no greater than one at the new
  source; at least one will equal one to root-solver tolerance.
- With entrained air initially at rest, `m_total u_out = m_H2 u_in` will
  continue to close to round-off.
- The source will be wider and slower.  The directions of the concentration,
  centre-height, and vertical-width changes are not assumed.

## Fixed comparisons and decision rule

- PRESLHY: peak-flow, momentum-filtered sensor-arc MG/VG/FAC2 on the arcs
  reached by both the adopted and candidate models.
- PRESLHY geometry: the same well-constrained, momentum-filtered vertical
  fits, using the common standard-deviation convention.
- Spadeadam horizontal tests 4 and 6 at the reported low mast winds.

The option is promoted only if it removes the invalid phase state without a
material aggregate regression: FAC2 must not decrease, VG must not increase
by more than 5 %, the mean vertical-width ratio must remain within 0.8--1.2,
and centre-height MAE must not rise by more than 20 %.  Independent Spadeadam
behaviour is reported as a guard against selection on PRESLHY alone.

The trial flow, entrainment coefficients, wind treatment, and Gaussian-width
definition remain fixed.  If the candidate fails, the phase-safe calculation
stays available as a model-form sensitivity bound rather than becoming the
corrected default.

## Results recorded after the fixed comparisons

At 1, 5, and 8 barg (288.6 K ambient), the former and candidate dry-air
source states are:

| barg | former T (K) | candidate T (K) | former air/H2 | candidate air/H2 | candidate H2 mass fraction |
|---:|---:|---:|---:|---:|---:|
| 1 | 20.027 | 68.612 | 1.550 | 4.240 | 0.191 |
| 5 | 20.086 | 68.235 | 1.265 | 3.865 | 0.206 |
| 8 | 20.127 | 67.956 | 1.073 | 3.611 | 0.217 |

Oxygen sets the boundary in all three cases and its saturation ratio closes
to one within `3.3e-13`.  Thus the thermodynamic invariant and the first
three predictions pass.

The repository's headline validation uses `flow_mean_gs`, despite an older
sentence in `nearfield.py` that incorrectly said peak.  Both choices were
therefore evaluated.  Results below are on the arcs reached by both models:

| PRESLHY source flow | model | n | MG | VG | FAC2 |
|---|---|---:|---:|---:|---:|
| window mean | adopted source | 45 | 1.068 | 1.611 | 0.778 |
| window mean | phase-safe candidate | 45 | **1.032** | 2.027 | **0.844** |
| window peak | adopted source | 44 | 0.988 | 1.532 | 0.773 |
| window peak | phase-safe candidate | 44 | **0.985** | 1.925 | **0.795** |

The mean bias and FAC2 improve, but VG worsens by about 26 % under either
flow convention.  The enlarged established-source plane also reduces the
headline mean-flow comparison from 62 reachable arcs to 45; the omitted
near-source measurements are not extrapolated backwards.

On the 20 common, well-constrained vertical fits using the window-mean flow,
the mean modelled/measured standard-deviation ratio improves slightly from
1.058 to 1.041.  However, centre-height MAE increases from 0.160 to 0.274 m
(71 %), and the mean signed error moves from +0.038 to +0.259 m.  The
candidate makes the 1.5 m releases rise too early.

The independent Spadeadam tests 4 and 6 also give a mixed-to-negative result:
MG/VG/FAC2 move from 0.869/1.476/0.667 to 1.367/1.633/0.833.  Centre heights
increase at every arc; for test 6 they move from 3.35/6.78/13.65 m to
4.72/8.69/16.14 m at 30/50/100 m.

## Decision

The candidate fails the pre-registered VG and centre-height criteria and is
**not promoted**.  It remains an off-by-default, coefficient-free upper-bound
path selected with `bulk_air_phase_safe=True` in `equivalent_source`, or
`bulk_air_phase_safe_source=True` in `hydrogen_jet` and the PRESLHY helpers.

This rejection is informative: merely skipping forward until all bulk dry
air has re-evaporated removes the invalid phase state but also removes the
dense condensed phase too early.  A production correction needs a transported
condensed-air mass, finite-rate re-evaporation, and phase-slip momentum term.
Those coupled equations cannot be replaced by a warmer single-phase source
without increasing trajectory and residual variance errors.
