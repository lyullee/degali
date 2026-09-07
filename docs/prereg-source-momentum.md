# Pre-registration: momentum of the flash-entrained air

> **Later width-normalisation audit:** the 0.933 and 0.731 spread ratios below
> used a measured e-folding width under a sigma label. Their
> common-standard-deviation values are 1.319 and 1.033. Momentum and
> concentration results are unaffected; see `gaussian-width-convention.md`.

Written before running the PRESLHY or Spadeadam comparisons with the new
option.

## Fault under test

The mixed equivalent source contains `1/Y` kg of hydrogen plus entrained air
per kg of hydrogen.  Its area is calculated from that total mass, but the
current source velocity is increased according to the density change,

```text
u_old = u_2 sqrt(rho_2/rho_3).
```

Consequently its total momentum is `m_H2 u_old/Y`, between 3.3 and 6.2 times
the incoming value over the 1 and 5 barg PRESLHY source states.  Station 2 to
3 is an entrainment and heating zone, not a pressure-expansion zone.  Both the
Sandia initial-entrainment model and Zhang et al.'s public condensation-model
precursor impose

```text
m_2 u_2 = m_3 u_3 + m_condensed u_condensed.
```

The testable no-slip limit is therefore

```text
u_3 = Y u_2
A_3 = m_H2 / (Y rho_3 u_3).
```

This is not yet a statement about the preceding under-expanded zone.  If
pressure thrust is material, its Station-2 momentum must be calculated
separately rather than manufactured during air entrainment.

## Prediction fixed in advance

- The total momentum invariant will close to round-off.
- The source speed will fall and its diameter will increase.
- The wider source should increase initial vertical spread, but the slower
  source may strengthen wind steering and dilution.  The direction of MG and
  trajectory change is not fixed in advance.
- A campaign-wide improvement is not assumed.  If the concentration variance
  or trajectory deteriorates materially, the option will remain experimental
  and the missing Station-2 pressure thrust will be recorded explicitly.

## Fixed comparisons

- PRESLHY: the same 69 momentum-driven sensor arcs, MG/VG/FAC2.
- PRESLHY geometry: the same 23 filtered vertical fits, mean
  `modelled_sigma_z/measured_sigma_z` and centre-height errors.
- Spadeadam tests 4 and 6: centre height at 30 and 100 m, without forced ground
  contact.

No coefficient will be fitted.  The only new switch selects the conservation
law used between the already-defined source planes.

## Results recorded after the fixed comparisons

The invariant closes to round-off in the regression test.  Relative to the
table-consistent but momentum-inconsistent source, the new boundary reduces
total source momentum by factors of about 3.28 at 5 barg and 6.20 at 1 barg.

Seven of the 69 PRESLHY arcs at 0.35/0.53 m lie upstream of the enlarged
established plane.  On the 62 common downstream arcs:

| source boundary | MG | VG | FAC2 |
|---|---:|---:|---:|
| mass/enthalpy consistency only | 1.442 | 1.563 | 0.69 |
| plus momentum consistency | **1.047** | **1.425** | **0.84** |

For the 23 fixed vertical fits, the mean spread ratio changes from 0.933 to
0.731.  That is a regression relative to the intermediate source but remains
above the historical 0.64.  Mean/median signed centre-height errors improve to
0.036/0.002 m, and mean absolute error is essentially unchanged at 0.145 m.

At the low reported Spadeadam winds, the six-arc free-detachment comparison
changes from MG/VG/FAC2 0.795/1.711/0.67 to 1.245/1.373/0.83.  Five of six arcs
improve.  At the mean mast winds, test 4 remains low at 1.43 m at 30 m, while
test 6 rises to 4.00 m at 30 m and 16.55 m at 100 m, preserving the observed
grounded/lifted regime split.

## Decision

The option is adopted as the corrected LH2 default.  It is a conservation law,
has no fitted parameter, and improves concentration statistics in both the
PRESLHY and independent Spadeadam comparisons.  The narrower-than-measured
vertical spread is reported as the remaining cost.  Station-2 pressure thrust
and condensed-phase slip are not hidden inside this balance and remain explicit
source-model uncertainties.

The phrase "narrower-than-measured" above is withdrawn after the Gaussian
convention audit. The adopted source has a +3.3% mean vertical-width bias.
