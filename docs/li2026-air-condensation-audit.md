# Li et al. (2026) air-condensation model: direct audit

## Outcome

The final paper was obtained and equations 13--24 were reproduced in
`degali.addons.cryogenic_air.li2026_zone3`.  For the paper's Test 2
Station-2 values, the implementation gives a Station-3 hydrogen mass fraction
of 0.790481 and velocity 473.878 m/s, consistent with Fig. 4 and satisfying
the printed mass, momentum and area balances to numerical precision.

It is **not a physically admissible replacement for the current PRESLHY LH2
source**.  Three independent reasons prevent adoption:

1. Test temperatures 51--55 K are below nitrogen's 63.15 K triple point, yet
   equation 17 uses a liquid-N2 saturation pressure and the text calls the
   collected phase liquid.  The physical phase there is solid.  PRESLHY's
   saturated-LH2 source, roughly 23--28 K, is further outside this domain.
2. Equations 15 and 16 add both condensed-N2 enthalpy and a latent term, with
   `gamma = h_N2 - h_LN2`.  They therefore reduce algebraically to assigning
   the condensed mass the **gas** enthalpy.  The advertised latent-heat effect
   cancels instead of appearing as the gas-to-condensate enthalpy difference.
3. The model drops nitrogen at zero axial velocity, retains all oxygen in the
   gas, and fixes Station 3 at the pure-N2 normal boiling point, 77.35 K.
   At LH2 temperature both N2 and O2 freeze/condense, and the appropriate
   boundary depends on component partial pressure rather than pure-component
   pressure at one atmosphere.

The code therefore rejects sub-triple-point use by default.  An explicit
`allow_subtriple_liquid_extrapolation=True` switch exists only to reproduce
the published calculation; it is not forwarded to the production or
validation paths.

## Exact paper closure

With entrained-air flow `m_air`, nitrogen condensation fraction `sigma`, and
zero condensed-phase velocity, equations 13--20 reduce to

```text
m3 = mH2 + m_air * (1 - YN2,air * sigma)
v3 = mH2 * v2 / m3
YH2,3 = mH2 / m3
sigma = 1 - Psat,N2(T2) / (0.78 * Pamb)
```

The implementation keeps the paper's 78/22 N2/O2 mole split and converts it
to mass fractions.  Pure-component gas enthalpies are used at both inlet and
outlet.  This is important: combining CoolProp's pseudo-pure `Air` enthalpy
with pure N2/O2 enthalpies mixes arbitrary reference states and changes Fig. 4
from about 0.79 to about 0.85.

## Stable solid-vapour pressure

`air_saturation_pressure` adds a separate phase-domain correction.  Below the
triple point it interpolates the tabulated pressure in Clausius-Clapeyron
coordinates,

```text
ln(P / Pa) = a / T + b
```

to the NBS Circular 564 solid-vapour tables:

| species | table range | points retained | interpolation residual |
|---|---:|---:|---:|
| N2 | 52--63.156 K | 6 | zero at table points |
| O2 | 36--54.36 K | 11 | zero at table points |

At 23 K the N2 solid-vapour extrapolation gives `1.89e-6 Pa`.  CoolProp's metastable
liquid ancillary gives about `1.94e4 Pa` and, more decisively, decreases as
temperature rises from 23 to 31 K.  That extrapolation cannot be used for an
LH2 air-freezing calculation.

The low-temperature tail is a constant-latent-heat extrapolation outside its
tabulated range.  It is sufficient to enforce the phase domain and show that
essentially all newly entrained N2/O2 is initially supersaturated, but not to
close condensate energy, particle residence time, settling or re-evaporation.

## What can reduce error next

The next defensible model is a distributed, three-component near-field state:

- transport gaseous H2, N2 and O2 separately;
- carry condensed/solid N2 and O2 mass as additional states;
- use stable solid-vapour pressure below each triple point and liquid-vapour
  pressure above it;
- put phase-change enthalpy into one energy balance exactly once;
- give condensate a finite slip/settling velocity and allow re-evaporation;
- hand off to the existing DEGADIS jet only after condensed mass is negligible.

The missing load-bearing input is not another paper equation but particle
kinetics: a measured or independently justified droplet/solid-particle size
distribution (or relaxation/settling time) for the PRESLHY release hardware.
Until that is available, zero slip and instantaneous full dropout should be
implemented only as uncertainty bounds and judged against both PRESLHY and
Spadeadam, never fitted to one campaign.

## Primary property source

NBS Circular 564, *Tables of Thermal Properties of Gases*, Table 7-11,
reports the solid-nitrogen vapour pressures and both species' phase points:
<https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular564.pdf>.

The solid-oxygen points are the Aoyama and Kanda measurements transcribed in
Georgia Tech Project A-593, Table 35:
<https://repository.gatech.edu/bitstreams/3b16b36e-0022-4035-b21f-a0bba98886b4/download>.
