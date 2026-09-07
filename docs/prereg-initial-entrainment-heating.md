# Pre-registration: initial entrainment and heating to the dry-air dew point

**Frozen before the first model result.** Do not edit above the `RESULTS`
line after a candidate result is known.

## Hypothesis

The source-conservative Gaussian boundaries become colder than both the
ambient-pressure source and the validity range of their gas-only air model.
HyRAM+ technical manual section 3.2.3, equations 58--64, publishes an omitted
plug-flow zone that entrains ambient air while conserving fuel mass, advective
momentum and total energy until a specified minimum temperature. The GUI
default is 0 K, so the zone is normally inactive; that default is not physical
for a 27--46 K cryogenic source mixing with air.

## Candidate fixed before running

Implement equations 58--64 before Gaussian establishment. Do not fit
`T_min`. Set it to the ideal dry-air dew temperature at the experiment's
ambient pressure, obtained from

```
y_N2 P / Psat_N2(T_dew) + y_O2 P / Psat_O2(T_dew) = 1
```

with `y_N2=0.78848834`, `y_O2=0.21151166`, `P=101325 Pa` and CoolProp
liquid-vapour saturation pressures. The fixed result is
`T_min=82.1508333771 K`. At and above this temperature, an infinitesimal first
drop of dry air is no longer stable at equilibrium.

For an ambient-pressure input plug, solve the published mass, momentum and
total-energy equations for fuel mass fraction `Y`; calculate density from the
published additive-volume equation and diameter from the dimensionally
correct mass identity `d=sqrt(4 mdot/(pi rho v))`. Calculate zone length from
equation 64 using the published source-momentum entrainment with
`beta_A=0.28`. The output plug then enters the four-flux `entrained_mass`
Gaussian boundary and the existing density-profile conserved ODE.

Component enthalpy is the same constant-`cp*T` model on both sides of this
isolated test. No Raman value changes a coefficient, source state or length.
The original source fuel mass, advective momentum and energy relative to
ambient air must be preserved. All defaults and the production JETPLU path
remain unchanged.

## Decision rule

1. Initial-zone fuel mass, momentum and total-energy residuals below `1e-10`;
   Gaussian boundary residuals below `1e-8` for all nine cases.
2. Downstream species and energy drift below `2e-4`.
3. All four printed Raman slopes within 25%, with median radial coefficients
   inside 33--64 (mass) and 21--49 (temperature).
4. The initial zone has `0<Y<1`, positive diameter/length and exactly
   `T_min`; all downstream states remain physical.
5. A 161-point radial rerun changes every slope by less than 0.2%.

Failure of any criterion rejects the candidate. `T_min`, latent heat,
entrainment, source temperature and virtual origin will not be tuned.

---

## RESULTS

Rejected. The 121-point result over all 549 frozen samples is:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21625 | 0.06302 | 0.01945 | 0.07817 | 53.49 | 33.93 | 2/4 |

Relative errors are -17.7%, -3.1%, **-31.2%**, and **+26.3%**. The initial
zone repairs the gas-only temperature domain and both mass metrics pass, but
centreline warming and temperature width fail criterion 3. The option is not
adopted alone.

A 161-point rerun gives 0.21627, 0.06302, 0.01945 and 0.07817; every slope
changes by less than 0.01%. This demonstrates that physically heating the plug
does not remove the downstream one-scalar profile limitation. Because the
initial zone and independent `Y`--`T` profiles address different equations and
were each fixed from external physics, their combination is tested in a new
pre-registration rather than altering this rejected candidate.
