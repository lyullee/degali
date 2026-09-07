# Pre-registration: flashing-source thermodynamic state consistency

> **Later width-normalisation audit:** width ratios in this historical test
> used a measured e-folding width under a sigma label. Multiply them by
> `sqrt(2)` for comparison with JETPLU's standard deviation; see
> `gaussian-width-convention.md`.

> Historical result. The later momentum-conservation correction supersedes
> the final performance numbers here; see `prereg-source-momentum.md` and
> `lh2-model-improvements-2026-09-03.md`. This record is retained unchanged
> below because it was written before that test.

## Defect found before the comparison

For Spadeadam test 6 the equivalent source is

```text
H2 mass fraction       0.401087
density                2.76899 kg/m3
H2 mole fraction       0.905840
```

but the first `JetPlume.derivatives` lookup maps its contaminant
concentration to

```text
H2 mass fraction       0.984672
density                1.12789 kg/m3
H2 mole fraction       0.998915
```

The source plane is therefore not thermodynamically self-consistent even
though its density, fraction and area are mutually consistent in the initial
condition.  The lookup table is still the pure saturated-hydrogen-to-air
mixing line.  `from_concentration` also assumes concentration is monotone on
that line; for cryogenic hydrogen it is not, so it selects the high-H2 branch
for the already air-entrained source.

A second defect is exposed at the same state.  The CoolProp backend queries
pure hydrogen at the *total* ambient pressure.  At 20.04 K and 1 atm that is
liquid (`71.2 kg/m3`), although hydrogen in the mixture is at its partial
pressure and is vapour (`about 1.36 kg/m3` on the metastable gas branch).
The plume thermodynamics are gas-phase thermodynamics; silently taking the
liquid branch is not a valid equation-of-state improvement.

## Frozen correction

1. Force the CoolProp contaminant property calls used by the **gas-plume
   paths** onto the gas branch (`T|gas`).  The explicitly modelled hydrogen
   droplet mass in `JetPlume` remains separate; this change does not turn a
   liquid source into gas by assumption.  The switch is scoped to the jet
   builders so the independently validated historical pool/dense-gas path is
   not silently changed before air-condensation physics is available there.
2. For `corrections=True`, treat the equivalent source as a secondary source,
   just as DEGADIS `DEG2S` treats a source that has already entrained air.
3. Convert the flash calculation's dry-air ratio `r` to a humid source:

```text
M = 1 + r (1 + humidity)
w_H2 = 1/M;  w_dry_air = r/M;  w_water = r*humidity/M
```

4. Evaluate the source enthalpy at the flash temperature with that composition
   and rebuild the adiabatic table from `(w_H2, w_dry_air, h_source)` to
   ambient.  Use the rebuilt endpoint density and H2 fraction in the expanded
   source geometry.  This makes the state entering the ODE and the state
   returned by its first thermodynamic lookup identical.
5. Keep the old pure-source table when `corrections=False`, preserving the
   Fortran oracle and the published baseline.

No coefficient or experimental result enters this correction.

## Predictions and decision rule

1. The rebuilt table's endpoint and the ODE initial state must agree in H2
   mass fraction, density, mole fraction and temperature to interpolation
   tolerance.  Failure is a coding error, not a model outcome.
2. Retaining the already entrained air should reduce initial H2 concentration
   and the excessive buoyancy impulse.  The modelled rise must not increase.
3. On the PRESLHY reduced campaign, MG must move closer to one without a lower
   FAC2; geometric variance and the measured/modelled vertical-width ratio
   are reported separately and may veto adoption.
4. On Spadeadam tests 4 and 6, the low-wind centre height at 30 m must fall and
   its 1 m concentration must move toward the measurement.  Test 4 may not be
   degraded outside a factor of two at 30 or 50 m.
5. The Houf--Schefer option remains off during this test.  Its failed hybrid
   cap is not allowed to hide or amplify the source-table effect.

## Results and adoption

The invariant passes.  For Spadeadam test 6 the state entering the ODE and
the first table lookup now both return H2 mass fraction `0.400001`, density
`3.050197 kg/m3`, mole fraction `0.905219` and temperature `20.03847 K`.
The former lookup returned mass fraction `0.984672` and density
`1.12789 kg/m3`, numerically deleting almost all flash-entrained air.

On 66 common PRESLHY arcs (gas-branch correction present in both sides), the
source-table correction changed:

| state table | MG | VG | FAC2 | model/measured sigma-z |
|---|---:|---:|---:|---:|
| pure-H2 endpoint | 1.409 | 1.834 | 0.788 | 0.971 |
| mixed flash endpoint | 1.414 | **1.522** | 0.712 | 0.933 |

The geometric-mean criterion and FAC2 criterion do not both pass, so the
change is not described as a uniform statistical win.  It is nevertheless
adopted as a correctness fix: conservation at the source is non-negotiable,
VG crosses the commonly used 1.6 threshold, mean centre-height bias falls
from `0.245 m` to `0.113 m`, and the alternative state violates its own
specified mass fraction and density before any experimental comparison.

The correction also invalidates the earlier decision to force permanent
ground contact.  With the consistent source, allowing detachment gives the
observed Spadeadam split between grounded high-wind test 4 and lifting
low-wind test 6 and improves the six-arc peak comparison from
`MG/VG/FAC2 = 0.310/12.81/0.50` to `0.795/1.71/0.67` at the lower measured
wind.  The public liquid-hydrogen path has been changed accordingly.
