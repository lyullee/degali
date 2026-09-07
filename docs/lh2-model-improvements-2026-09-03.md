# LH2 model physics correction — 2026-09-03

This note supersedes the liquid-hydrogen performance numbers in older working
notes where they conflict. Historical calculations remain reproducible with
both `source_table_consistency=False` and
`source_momentum_consistency=False`.

## 2026-09-04 near-nozzle thermal follow-up

The separate conserved-energy Gaussian jet now supports explicit ambient
humidity with local H2/N2/O2/H2O equilibrium. Water vapor sensible enthalpy,
freezing/vaporisation latent heat and condensed-water volume all remain
inside the audited total-energy flux. The 100%-RH Hecht--Panda upper bound
has maximum nine-case boundary/species/energy residuals of `1.08e-15`,
`4.80e-6` and `9.71e-5`, respectively. It overpredicts the printed
temperature centreline slope by 44.8% and is not adopted.

This falsifies the earlier implicit assumption that laboratory moisture is a
minor correction: the dry equilibrium model is about 31% low, whereas a
40%-RH sensitivity is about 9% low and places all four slopes inside the
pre-registered 25% band. Hecht--Panda describe condensed moisture in the
images but do not state RH, so the result is parameter sensitivity rather
than predictive validation. The experimental environmental log or raw-data
metadata is the next required item. See
`prereg-humid-air-frost-upper-bound.md` and
`prereg-raman-case-coverage-correction.md`.

The percentages in the preceding humidity paragraph are retained from the
then-active 2017 conference benchmark. A later source audit made the final
2019 journal fits authoritative and exposed an unresolved tenth legend entry;
see `hecht-panda-journal-benchmark-correction.md` before quoting any 4/4 count.

The next isolated correction replaces constant ambient-temperature heat
capacities with low-density Helmholtz ideal-gas component enthalpies inside
the dry N2/O2 equilibrium calculation. It is accepted: the corrected
369-point slopes have final-journal relative errors -24.86%, -15.79%, -21.38%
and +24.56%,
with both radial coefficients inside their source ranges. Maximum boundary,
species and energy residuals are `1.50e-14`, `1.53e-5` and `3.77e-5`. This is
the first Raman candidate to combine end-to-end conservation with a
provisional 4/4 aggregate-fit result. It is the recommended dry-air
axisymmetric research configuration,
not a humidity-conditioned production validation; see
`prereg-phase-temperature-dependent-enthalpy.md`.

## Adopted correction 1 — preserve the flashing-source state

The equivalent-source calculation already contains air entrained while the
flashed liquid evaporates. The jet ODE nevertheless looked the initial
contaminant concentration up on a *pure saturated-H2-to-air* table. Because
that table is non-monotone in H2 concentration, the lookup selected its
high-H2 branch and numerically removed most of the entrained air.

Spadeadam test 6 exposed the discontinuity:

| quantity | equivalent source | former first ODE lookup | corrected lookup |
|---|---:|---:|---:|
| H2 mass fraction | about 0.400 | 0.984672 | 0.400001 |
| density, kg/m3 | about 3.05 on the humid table | 1.12789 | 3.050197 |
| H2 mole fraction | about 0.905 | 0.998915 | 0.905219 |
| temperature, K | 20.038 | — | 20.03847 |

The corrected path rebuilds the adiabatic table from the mixed flash state.
For flash dry-air ratio `r` and ambient absolute humidity `h`,

```text
M = 1 + r(1+h)
w_H2 = 1/M
w_dry_air = r/M
w_water = rh/M
```

The source enthalpy is evaluated at the flash temperature and composition.
The table endpoint density and fraction then define the expanded source plane.
No experimental coefficient enters this correction.

CoolProp contaminant properties are forced onto its metastable gas branch only
on the corrected gas-jet paths. Remaining H2 droplets are already carried as
separate mass. This switch is deliberately not global: changing the
independently validated NASA pool path before adding air-condensation physics
spoils two of its four regime classifications.

## Adopted correction 2 — allow a grounded plume to detach

The former public assessment forced `JetPlume.ground_effect=True`. That option
removes the buried entrainment perimeter and scales upward buoyancy by the
fraction of the ellipse that has cleared the ground. With the old source table
this masked excessive rise. With the corrected source it permanently holds
both Spadeadam horizontal plumes down.

That is contrary to the experiment. Mack et al. report that high-wind test 4
remained grounded and low-wind test 6 lifted; their EFFECTS calculation and an
independent CFD study reproduced the split. DEGALI now does too without a
forced ground condition:

| mean of reported mast winds | test 4 | test 6 |
|---|---:|---:|
| centre height at 30 m | 1.43 m | 4.00 m |
| centre height at 100 m | 4.67 m | 16.55 m |
| observed interpretation | grounded | lifted |

At the lower reported wind, the six peak-concentration arcs change as follows:

| treatment | MG | VG | FAC2 |
|---|---:|---:|---:|
| forced ground contact | 0.359 | 3.93 | 0.33 |
| free to detach | **1.245** | **1.37** | **0.83** |

Five of six arcs improve. Only the test-6 30 m peak reading favours the
attached solution. These are sparse,
meandering-plume measurements, and the 100 m arc maxima are lower bounds when
the narrow plume passes between sensors, so they are not used for coefficient
fitting.

## Adopted correction 3 — conserve momentum across initial air entrainment

The expanded-source area used the total mixed mass flow, but its velocity was
increased as `u2 sqrt(rho2/rho3)`. The resulting total mixture momentum was
not the incoming hydrogen momentum: it was about 3.28 times larger at 5 barg
and 6.20 times larger at 1 barg. The model therefore created momentum while
stationary air was entrained.

The adopted no-slip balance is

```text
m_H2 u2 = (m_H2/Y) u3
u3 = Y u2
A3 = m_H2 / (Y rho3 u3)
```

where `Y` is the hydrogen mass fraction at the established source plane. This
is the same initial-entrainment balance used in the current Sandia HyRAM jet
implementation. It contains no fitted coefficient. Any pressure thrust from
the preceding under-expanded zone must be calculated explicitly at Station 2;
it cannot be inserted by accelerating already entrained air. Condensed-air
slip would add a separate momentum term and remains unresolved.

## Campaign result after correction

PRESLHY E3.5, momentum-driven horizontal trials, model evaluated at each sensor
and maximised over the same arc as the measurement:

| model | n | MG | trial-bootstrap 95% CI | VG | FAC2 |
|---|---:|---:|---:|---:|---:|
| historical reconstruction, full array | 69 | 1.121 | 0.736–2.264 | 16.85 | 0.80 |
| historical reconstruction, common range | 62 | 1.156 | 0.720–2.521 | 23.10 | 0.77 |
| mass/enthalpy consistency only, common range | 62 | 1.442 | 1.072–1.849 | 1.563 | 0.69 |
| mass/enthalpy/momentum consistency, common range | 62 | **1.047** | **0.759–1.402** | **1.425** | **0.84** |

The seven 0.35/0.53 m arcs are upstream of the enlarged established-source
plane. They are not extrapolated backwards; the corrected validation range
therefore begins at the next instrumented arc, 0.79 m. On this common basis
the complete correction passes all three Hanna screening criteria and moves
MG close to unity. The trial-level residual has no significant monotone
association with pressure, flow, wind, ambient temperature, humidity, release
height or orifice size over the nine trials, so a universal scalar coefficient
is not justified.

For the 23 well-constrained, momentum-driven vertical fits, mean
`modelled_sigma_z/measured_sigma_z` is **1.033** (historical model 0.901).
The earlier values 0.731 and 0.64 compared JETPLU's Gaussian standard
deviation with a measured e-folding width that had been mislabeled `sigma_z`.
Those width parameters differ by exactly `sqrt(2)`. The fitting code and the
reduced data now use `exp(-0.5*((z-zc)/sigma_z)^2)` consistently; fitted
curves, centres, concentrations, R2 values and the 23-point selection are
unchanged. See `gaussian-width-convention.md`.

Mean and median signed centre-height errors remain 0.036 and 0.002 m; mean
absolute error remains 0.145 m. In the 3–7 m band, MG/VG/FAC2 remain
0.904/2.115/0.67; excluding trial 20 they remain 0.736/1.775/0.69. The former
27% vertical-width deficit was a validation-normalisation error, not missing
plume entrainment. The corrected mean width bias is +3.3%.

## Investigated but not adopted

### Circular cross-section and extra entrainment

A no-coefficient upper-bound diagnostic preserved `sigma_y sigma_z` but
forced `sigma_y=sigma_z`. On the corrected Gaussian convention it moves the
vertical-width ratio from 1.033 to 1.097, so it is unnecessary and worsens the
small mean bias. It was not added to production code. See
`prereg-cross-section-split.md`.

Increasing total entrainment is not supported as a generic cryogenic-jet fix
either. Hecht and Panda's controlled cryogenic-hydrogen jets had a measured
mass-fraction half-width growth near 0.065, below the roughly 0.10--0.11
room-temperature values they cite. The former apparent need for about 40%
more vertical spread was instead the exact Gaussian-width normalisation error.

Source: Sandia publication record and public manuscript,
<https://www.sandia.gov/research/publications/details/mixing-and-warming-of-cryogenic-hydrogen-releases-2019-04-02/>,
<https://h2tools.org/sites/default/files/2019-09/123.pdf>.

### Houf–Schefer buoyancy entrainment

The official HyRAM equations and current Sandia implementation were added as
an off-by-default option. Its pure-plume cap conflicts with the existing
DEGADIS shear closure: it improves some concentration statistics, but moves
the corrected standard-deviation width ratio from 1.033 to about 0.91. No
hybrid coefficient was fitted, so it remains experimental. See
`prereg-houf-entrainment.md`.

Sources:

- SAND2015-10216, *HyRAM 1.0 Technical Reference Manual*:
  https://www.osti.gov/servlets/purl/1814840
- current Sandia HyRAM source:
  https://github.com/sandialabs/hyram/blob/master/src/hyram/phys/_jet.py

### Ground-layer top entrainment

`ground_effect` omitted DEGADIS's friction-velocity/Richardson-number
top-surface entrainment. The existing equation was added over the
ground-contact chord as an off-by-default option. In its pre-registered test
on the earlier table-consistent source, it moved the six-arc Spadeadam
comparison only from MG/VG/FAC2 0.323/13.46/0.50 to 0.333/12.75/0.50 and
changed PRESLHY negligibly. It was not promoted after the momentum correction:
the public path does not force ground contact, and the closure cannot act
without it. See `prereg-ground-layer-entrainment.md`.

## Largest remaining physics gap — condensed air and phase slip

The corrected equivalent source conserves its declared gas-mixture mass,
enthalpy and momentum, but still assumes the entrained dry air remains
gaseous.  Its 1--8 barg endpoint is only about 20.0--20.1 K and contains
1.07--1.55 kg dry air per kg H2, so that phase assumption is not physically
available. Two primary studies show why this is the next source-model item:

- Li et al., *Modeling of cryogenic hydrogen jets with air condensation*,
  International Journal of Hydrogen Energy 250 (2026), article 156128,
  DOI 10.1016/j.ijhydene.2026.156128
  (<https://www.sciencedirect.com/science/article/pii/S0360319926027667>).
  Its improved integral model adds
  nitrogen-condensation latent heat to the energy balance and momentum loss
  from liquid/gas slip in the initial entrainment/heating zone.
- Sun et al., *Phase change modeling of air at the liquid hydrogen release*,
  International Journal of Hydrogen Energy 50 (2024), 717–731,
  DOI 10.1016/j.ijhydene.2023.06.201. It identifies two opposing effects:
  condensed particles raise mixture density, while released latent heat warms
  the cloud and increases buoyancy. During a turbulent continuous pool
  release the condensed-air volume is small and sensitivity to its Lee-model
  coefficient is weak; after release termination it grows.

Zhang et al.'s public conference precursor, *Analytical Model of Cryogenic
Hydrogen Releases*
(<https://hysafe.info/uploads/papers/2023/189.pdf>), was also read. It prints
the Zone-II mass, momentum and
energy balances, but its printed equations contain dimensional and sign
problems.  The final Li et al. (2026) article was supplied on 2026-09-03 and
read directly.  It corrects several conference-print errors and specifies
`v_LN2 = 0`, `T3 = 77.35 K`, equations 13--24 and `beta_A = 0.281`.

Direct reproduction also exposes limits that prevent transplantation to
PRESLHY.  Equation 17 evaluates liquid-N2 saturation below the nitrogen triple
point (all four cold validation cases are 51--55 K), equations 15--16 make
`h_LN2 + gamma = h_N2` so the latent term algebraically cancels, and O2 is
forced to remain gaseous.  A phase-domain guard, exact paper reproduction and
NBS solid-vapour correlations are now in `addons.cryogenic_air`; none is
enabled in the production path.  See `li2026-air-condensation-audit.md`.

Sandia's April 2025 HyRAM+ 6.0 technical manual (SAND2025-04942, equations
58--64) and the corresponding open-source `_dev_plug` implementation were
also checked.  They conserve mass, momentum and energy through an initial
plug-flow heating zone until a user-specified minimum temperature.  The
manual explicitly says the default is 0 K, so this zone is disabled in the
GUI; the current 6.1 source retains the user-selected boundary.

DEGALI now contains a stricter, coefficient-free sensitivity path.  Instead
of choosing 65 K, it solves the component enthalpy balance and advances to the
first temperature at which the N2, O2 and Ar partial pressures do not exceed
their saturation pressures.  For 1--8 barg LH2 this gives 67.96--68.61 K,
3.61--4.24 kg dry air/kg H2, with oxygen setting the boundary.  Select it with
`bulk_air_phase_safe=True` in `equivalent_source` or
`bulk_air_phase_safe_source=True` in the jet helpers.

It is **not adopted as the corrected default**.  On the 45 common PRESLHY
mean-flow arcs it moves MG/FAC2 from 1.068/0.778 to 1.032/0.844, but VG worsens
from 1.611 to 2.027.  On 20 common vertical fits the width ratio improves from
1.058 to 1.041 while centre-height MAE worsens from 0.160 to 0.274 m.  The
independent Spadeadam MG/VG also moves from 0.869/1.476 to 1.367/1.633.  A
single warm all-gas handoff removes the invalid state by removing the dense
condensed phase too early.  It therefore remains a model-form sensitivity
bound; see `prereg-bulk-air-phase-boundary.md`.

The formal paper is no longer missing.  A transported N2/O2 candidate has now
also been implemented and tested prospectively at fixed 1, 10 and 100 um
particle diameters.  It starts from a pressure- and enthalpy-closed solid-air
H2 endpoint, conserves condensed component mass, permits re-evaporation and
settling, removes dropped-particle momentum and energy, and preserves the
source-zone axial distance at the JetPlume handoff.  All 11 PRESLHY/Spadeadam
source marches close relative mass, momentum and energy residuals below
`1e-10`.

It is **not adopted**.  On common PRESLHY arcs the three candidate
MG/VG/FAC2 values are respectively 0.968/2.058/0.846,
0.982/2.150/0.846 and 1.362/8.946/0.775.  Centre-height MAE rises from
0.187 m on the common fits to 0.285, 0.295 and 0.482 m.  Independent
Spadeadam MG/VG/FAC2 moves from 0.869/1.476/0.667 to
1.353/1.598/0.833, 1.357/1.600/0.833 and 1.460/1.694/0.833.  Fine retained
particles reproduce the over-buoyant warm-source failure; coarse dropout
makes the remaining source more H2-rich and worse.  Full equations,
pre-registered criteria and ranges are in
`prereg-condensed-air-particle-transport.md`.

Production use still needs measured particle-size/residence data and a
two-axial-velocity slip closure.  The present research path uses a common
axial velocity and does not yet transport argon, ambient water ice or CO2
frost, so no particle diameter may be selected by whichever field residual
happens to score best.

### Storage-to-source total energy and pressure thrust (2026-09-04)

A subsequent conservation audit found that the solid-air endpoint spent all
storage enthalpy on phase enthalpy and the transport march then added kinetic
energy.  The research source can now solve phase and kinetic energy together
with `source_total_energy_consistency=True`.  A separate option evaluates the
HyRAM+ homogeneous-equilibrium critical throat and energy-conserving
Yuceil--Otugen pressure-thrust plane.  Both require the explicit transported
condensed-air source and remain off by default.

The total-energy correction is real but is **not adopted as a field-model
improvement**.  On 38 common PRESLHY arcs E changes MG/VG/FAC2 from
1.056/1.630/0.789 to 0.965/2.094/0.842, while centre-height MAE worsens from
0.187 to 0.284 m.  Independent Spadeadam MG/VG/FAC2 worsens from
0.869/1.476/0.667 to 1.336/1.583/0.833.

The pressure-thrust candidate reduces centre-height MAE from 0.186 to
0.145 m on 15 common PRESLHY fits, but fails every campaign-wide guard:
common-arc MG moves 1.182 to 1.201, width ratio remains 1.083 and Spadeadam
MG/VG becomes 1.466/1.688.  More fundamentally, PRESLHY trials 10, 22 and 25
have no ambient-pressure energy state under the geometric-orifice thrust
model.  This diagnoses line/orifice incompatibility rather than a parameter
to clip.  Full prospective criteria, source ranges and results are in
`prereg-source-total-energy-pressure-thrust.md`.

Li et al. equation 22 also exposed a smaller coordinate omission: the
evaporation endpoint already contains entrained air but had been placed at
zero source distance.  The optional `evaporation_zone_distance=True` now adds
the paper's 0.012--0.050 m PRESLHY / 0.072--0.078 m Spadeadam Zone-III
length.  It changes no source state.  Prospective testing rejects it as an
error-reduction mechanism: E's centre MAE improves only 0.284 to 0.280 m,
PRESLHY MG worsens 0.965 to 0.954, and Spadeadam MG remains 1.331.  See
`prereg-evaporation-zone-origin.md`.

The opposite condensed-phase momentum limit has also been closed.  Li et al.
assume negligible condensed-phase axial velocity; the S0 research bound now
applies that assumption to stable solid N2/O2 while conserving separate gas
kinetic energy and the enthalpy carried into deposition.  All eleven sources
close below `6.2e-12`, but S0 is decisively **not adopted**.  Against the
corrected default, PRESLHY common-arc MG/VG changes from 1.047/1.425 to
1.285/2.308, centre MAE / width ratio from 0.145 m / 1.033 to
0.321 m / 1.269, and Spadeadam MG/VG from 0.869/1.476 to 1.563/1.827.
Because both the fully carried and stationary limits over-rise, no
intermediate velocity-relaxation scalar will be selected from these residuals.
See `prereg-stationary-condensate-bound.md`.

The next finite-rate audit now adds a source-backed Ranz--Marshall transfer
number and a deliberately fastest heat-limited sublimation lifetime.  On the
eleven eligible PRESLHY/Spadeadam sources, an initial 1 um N2 particle has a
7--40 mm minimum survival distance, while 10 um often persists 0.05--0.26 m
and 100 um usually survives beyond the equilibrium handoff.  This proves that
kinetics can matter for coarse particles, but also that the unmeasured
particle number/size and new-particle nucleation determine the answer.  No
delay multiplier or diameter is fitted; see `finite-rate-sublimation-audit.md`.

### Raman-validated conserved-energy free jet (2026-09-04)

Hecht and Panda's nine simultaneous Raman concentration/temperature jets add
a direct near-nozzle test that the outdoor campaigns cannot provide. A
traceable transcription uses only Table 1 and the numerical fits printed in
Figures 5--8; no plotted points were digitized. The first JETPLU comparison
showed that changing a single entrainment coefficient could not reproduce
both centreline dilution and radial width.

DEGALI now contains an independent implementation of the published
axisymmetric Gaussian conservation model. It integrates total mass,
horizontal and vertical momentum, hydrogen species and total enthalpy plus
kinetic energy, uses the published `lambda=1.16`, `beta_A=0.28` and
`alpha<=0.082`, and converts the plug source to its established Gaussian
profile before integration. The physical fluxes are integrated directly and
their Jacobian is evaluated numerically, rather than transcribing another
program's expanded equations.

Across all 549 frozen points, relative errors in the four printed slopes are
**-3.1%, -11.2%, -19.9% and +14.7%**; all pass the prospective 25% band.
The same values from a separately installed official HyRAM 6.1 oracle differ
from the independent implementation by at most 0.75%. Species and total-energy
conservation *within the established-flow ODE* is pinned by tests. A later
end-to-end audit found that the published algebraic plug-to-Gaussian boundary
loses 8.9--14.8% species and 20.4--29.0% relative energy before that ODE
begins. The original Fortran-compatible path remains unchanged.

This model is retained as an **official-model reproduction and Raman
benchmark**, not yet accepted as an end-to-end source-conservative model. Two
coefficient-free conservative boundary replacements close below `5e-12` but
each fails the centreline-temperature slope at about -31%, implicating
missing low-temperature mixing-layer physics rather than a tunable virtual
origin. The formerly missing cross-wind handoff has now been implemented
separately. It projects mass, hydrogen and vector momentum exactly and
transfers the accepted near-field phase/enthalpy radial profile before
screening total energy, H2 width and centre temperature. The phase-profile
candidate passes the frozen representative interface test; two simpler
mappings are rejected and retained. Independent PRESLHY and Spadeadam pilots
are now complete. On seven PRESLHY releases, resuming JETPLU's local density-
scaled shear entrainment after the handoff corrects the mean vertical-width
ratio from 0.710 to 0.994 while retaining FAC2 0.952. The complete path is not
promoted: its fully four-flux source passes only five of seven single-scalar
JETPLU interfaces, and searching 10D--20D worsens the temperature/profile
mismatch. The pre-registered independent-energy crosswind state has since
recovered all seven interfaces and completed a step-refined field run. It
improves VG, FAC2 and width but worsens MG and centre-height MAE, so it remains
research-only. See `preslhy-independent-energy-interface-results.md` and
`prereg-crosswind-independent-energy-state.md`.

## Other residual limitations

- Ground heat transfer is steady and cannot represent progressive substrate
  cooling during a long release. Mack et al. show this can switch test 6
  between grounded and lifted regimes.
- The atmospheric boundary layer is steady. Spadeadam wind speed varied
  20–30% with larger peaks and wind direction by roughly 10–15 degrees with
  larger excursions; a single deterministic centreline cannot reproduce plume
  meander or peak/mean ratios.
- The corrected default has no transported N2/O2 phase.  The off-default
  research source assumes instantaneous phase equilibrium and includes
  settling, but not finite-rate condensation/sublimation or a separate axial
  particle velocity.  Those require additional transported state variables,
  not another scalar correction to entrainment.
- The measured flow and nominal opening do not define a unique pressure-thrust
  source when upstream line loss and effective discharge area are unknown.
  Three PRESLHY sources fail total-energy admissibility if the full nominal
  area is assigned the HEM throat pressure.

These limits should be propagated as scenario/sensitivity ranges in a safety
assessment rather than hidden inside one best-fit coefficient.
