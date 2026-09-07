# Finite-rate condensed-air sublimation audit

Date: 2026-09-04.  This is an applicability audit, not a field-data-fitted
candidate and not a pre-registration.  No PRESLHY or Spadeadam concentration
residual is used to select a particle size or transfer coefficient.

## Question

The transported condensed-air source currently imposes N2/O2 phase
equilibrium at every axial station.  The next physical question after the
carried and stationary momentum bounds is whether a particle can actually
heat and sublime over that station's residence time.

A complete answer needs two additional transported states per condensed
species (particle mass and particle enthalpy/temperature) and a particle
number flux.  The number flux is load bearing: it sets the interfacial area
available for both condensation and sublimation.  Neither PRESLHY nor
Spadeadam reports condensed-air particle number or size, and the initial
particles do not close the problem because newly entrained air can nucleate
while the plume is below the N2/O2 phase boundaries.

## Source-backed lower bound

Sandia's thin-skin spray formulation writes the droplet mass, momentum and
energy balances separately and uses

```
Nu_0 = 2 [1 + Re_p^(1/2) Pr^(1/3) / 3]
```

before its Stefan-flow correction.  See SAND2002-3419, sections 2.1--2.2:
<https://www.osti.gov/servlets/purl/807053>.

The 2023 single-moving-LN2 experiment found that all six tested classical
gas-phase models underpredicted the high evaporation rate and that omitting
the Stefan-blowing suppression gave the better direction; its proposed model
was within 1.1% for relative velocity above 0.5 m/s:
<https://doi.org/10.1016/j.ijheatmasstransfer.2022.123584>.

For a deliberately fastest possible screen, hold the particle at its initial
temperature, assign every watt of convective heat to latent heat, and omit
Stefan resistance.  A spherical balance then gives

```
d(d_p^2)/dt = -4 Nu_0 k_g (T_g - T_p) / (rho_p L)
t_min = rho_p L d_p^2 / [4 Nu_0 k_g (T_g - T_p)]
```

`minimum_heat_limited_sublimation_time` implements this equation with the
stable solid density and sublimation enthalpy already used by the phase
ledger.  Because every assumption maximises evaporation, the result is a
**lower bound on survival time**, not a kinetic prediction.

## Campaign-condition screen

The lower-bound `d_p^2` erosion was integrated along the already-computed
equilibrium-source temperature, velocity and distance histories.  The
particle was held at 20.4 K, hydrogen properties were used for the
gas-dominated film, terminal settling supplied relative Reynolds number, and
the fixed sizes remained the previously declared 1/10/100 um cases.

| source | equilibrium handoff (m) | 1 um N2 disappears (m) | 10 um | 100 um |
|---|---:|---:|---:|---:|
| PRESLHY 10 | 0.532 | 0.023 | 0.151 | after handoff |
| PRESLHY 11 | 0.252 | 0.030 | 0.241 | after handoff |
| PRESLHY 12 | 0.125 | 0.021 | after handoff | after handoff |
| PRESLHY 20 | 0.477 | 0.008 | 0.061 | 0.416 |
| PRESLHY 21 | 0.239 | 0.007 | 0.049 | after handoff |
| PRESLHY 22 | 0.533 | 0.020 | 0.141 | after handoff |
| PRESLHY 23 | 0.252 | 0.024 | 0.166 | after handoff |
| PRESLHY 24 | 0.126 | 0.019 | after handoff | after handoff |
| PRESLHY 25 | 0.533 | 0.020 | 0.141 | after handoff |
| Spadeadam 4 | 0.797 | 0.040 | 0.255 | after handoff |
| Spadeadam 6 | 0.745 | 0.038 | 0.261 | after handoff |

At equal diameter and thermal history, the solid-density-times-latent product
for O2 is about 1.41 times the N2 value, so its lower-bound lifetime is longer.
The qualitative size classification is unchanged.

The screen does **not** say that the whole condensed phase disappears at the
listed location.  It follows an initial particle only.  The equilibrium path
continues creating condensed material from newly entrained air, and a
finite-rate model must predict the nucleation/growth surface for that material.

## Decision

Do not add a fitted sublimation-delay multiplier.  For 1 um particles, even
the deliberately fastest bound lasts only 7--40 mm; a more complete model may
lengthen that, but this scale cannot be promoted as the explanation for the
metre-scale trajectory error.  At 10--100 um, finite-rate survival is plainly
important, but the answer is dominated by an unmeasured particle size and by
nucleation surface area.  Selecting either from PRESLHY residuals would be
calibration without an independent observable.

The code therefore gains the auditable transfer-number and minimum-lifetime
functions, while the validated default and the equilibrium research source
remain unchanged.  A two-temperature/two-phase plume should only be promoted
after obtaining at least one of:

1. N2/O2 particle-size or number-density measurements within the first metre;
2. time-resolved optical extinction/Mie scattering that constrains particle
   surface area and disappearance;
3. gas and condensed-phase temperature or velocity histories through the
   20--80 K region.

Sandia's current cryogenic-release facility explicitly lists Raman imaging
and particle-image velocimetry as relevant diagnostics, but its public page
does not provide the needed condensed-air size distribution:
<https://crf.sandia.gov/research/experimental-capabilities/laser-diagnostics-for-cryogenic-gas-releases/>.

The closest public temperature/concentration field is Hecht and Panda's
48--63 K Raman study, OSTI record **1529288** / SAND2018-7244J:
<https://www.osti.gov/biblio/1529288>.  The accepted manuscript and plotted
median profiles are public, but the OSTI record links no raw image or numeric
field data product.  Those measurements can independently test jet warming
and scalar spreading at 20--60 mm, but they neither reach the saturated-LH2
20 K source state nor identify condensed N2/O2 particle size.
