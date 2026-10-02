# LH2 data readiness for the next model stage

Last reviewed: 2026-09-15

This register separates data that can constrain the model from material that
is useful only for context or bounds. Third-party raw files are kept under the
ignored `reference/` tree and are not distributed with DEGALI.

## Decision summary

| question | usable evidence now | readiness |
|---|---|---|
| PRESLHY E3.5 source mass flow | raw meter histories plus the D4.8 pressure-drop reconstruction | usable with trial-specific quality flags |
| PRESLHY E3.5 nozzle state | pipe/nozzle pressure and nozzle-flow/wall thermocouples | usable, subject to two-phase sensor interpretation |
| PRESLHY E3.5 atmospheric boundary | ambient temperature, RH, wind speed and direction | usable at the reported reference height; no local 3-component wind |
| PRESLHY E3.5 plume temperature and H2 | thermocouple arrays, Xensor and Drager channels | usable with alignment, saturation and coordinate caveats; paired 1.78/4 m vertical thermal/species moments now reduced for Trials 10/23 |
| transported turbulent energy | no turbulence channel in E3.5; a public 5 bar/50 K LES axial-RMS contour now supplies a non-experimental realizability envelope | conservation/realizability boundary only; epsilon remains unidentified |
| condensed-air relaxation | the E3.5 public catalogue lists thermal, pressure, H2/O2, weather and video channels but no particle-size, phase-fraction or separate particle-velocity instrument | no-slip/stationary kinematic limits only; finite slip unvalidated |
| ELVHYS Tests 10 and 11 geometry | archive dictionary separates the nozzle at 250 mm from bottom samplers at 200 mm | usable for archive reduction |
| ELVHYS Tests 10 and 11 hydrogen flow | `FLMT` is ventilation-air flow only | missing |
| Sandia cryogenic free jets | Raman concentration and temperature profiles are available in the paper/local reduction | usable for mean thermal/species profiles, not a complete TKE boundary |
| Sandia 2025 steady-crosswind LH2 pools | 16 real-LH2 spill tests; public condition tables and downwind H2 contour plots | usable as a separately declared pool-vapour validation screen after plot reduction; raw histories remain request-only |

## 1. PRESLHY E3.5 outdoor releases

Primary archive: HSL/PRESLHY E3.5, DOI
[10.35097/1481](https://doi.org/10.35097/1481).

The seven trials currently selected for the downstream thermal/buoyancy work
are 10, 11, 12, 22, 23, 24 and 25. Their local workbooks contain the same five
data groups:

| sheet | relevant measurements | intended use |
|---|---|---|
| `Flexlogger` | mass-flow indication, pipe and nozzle pressure, nozzle-flow and nozzle-wall temperature, ambient temperature, RH, wind, thermocouple arrays | source history and thermal trajectory |
| `Flowmeter` | mass flow, indicated density, meter temperature, drive gain and totalizer | source-flow quality and phase-instability screening |
| `Draeger` | H2 ppm, vol% and %LEL at the stand network | far-field concentration; values at the 4 vol% ceiling are censored |
| `Xensor` | fast hydrogen-sensor outputs | near-source timing/profile reduction after calibration mapping |
| `LocalWeather` | temperature, humidity, pressure, wind and dew point | atmospheric boundary and condensation bounds |

### Source-flow hierarchy

Do not select one flow definition silently. Preserve three values where they
exist: raw meter history, stable/peak meter estimate, and the report-derived
flow. The D4.8 appendix gives the most important reconstruction:

- Test 11, 12 mm nozzle at 5 bar tanker pressure: approximately 265 g/s.
- Test 12, 6 mm nozzle: approximately 90--100 g/s.
- Test 10, 25.4 mm open pipe: approximately 298 g/s, inferred from
  `sqrt(3.8/3.0)` times the Test-11 flow because the two-phase meter output is
  unreliable.

The report-derived value is the preferred boundary when drive gain, density,
or void fraction shows that the Coriolis meter is outside reliable operation.
The meter time series remains useful for release timing and stability.

### What E3.5 cannot identify

None of the seven workbooks contains velocity profiles, three-component
velocity fluctuations, Reynolds stresses, turbulent kinetic energy, or a
dissipation rate. Concentration and temperature alone cannot independently
identify a transported-TKE closure. A new turbulence coefficient derived only
from plume height would therefore be a calibration, not an independent
physical validation.

The public archive catalogue was rechecked on 2026-09-15.  It lists
near/far-field thermocouples, pipe/tank pressure, mass flow, H2 concentration
and ppm channels, O2 depletion, near/far-field weather, humidity and video;
it does **not** list particle imaging, size distribution, condensed N2/O2
mass, or a separate particle-velocity measurement.  Its description records
visible condensed air near the release/impingements, but that observation
cannot determine a two-velocity or particle-size closure.  The full archive
is about 11.3 GB and includes video; it was not copied into this project.

## 2. ELVHYS WP4.2 confined releases

Primary archive: DataverseNO, DOI
[10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H), CC0.

Only the README, metadata, and ten Test-10/Test-11 streams were selected from
the 198-file, approximately 1.43 GB archive. They provide 20 Hz concentration,
temperature, release/nozzle pressure, weather and fan-flow histories.

### Coordinate interpretation

The archive data dictionary gives separate coordinates:

- horizontal LH2 inlet: `(x, y, z) = (20, 500, 250) mm`;
- Bottom 1--9 temperature/concentration samplers: `z = 200 mm`.

Accordingly, the archive should be reduced with a 250 mm nozzle elevation and
200 mm bottom-sampler elevation. The 200 mm nozzle entry in D4.6 is retained
as a report discrepancy, but it no longer makes the archive coordinate origin
ambiguous: the machine-readable data dictionary explicitly distinguishes the
two physical objects.

### Remaining source limitation

`ELE402HSE010FLMT...csv` and `ELE402HSE011FLMT...csv` contain only
`FanFlowMeter`, around 500 L/min. This is ventilation-air flow, not hydrogen
mass flow. The pressure and nozzle-temperature histories can constrain a
source calculation, but they do not constitute an independently measured H2
mass-flow boundary. These tests remain a confined-release falsification case,
not a clean calibration case for the outdoor free-jet model.

## 2b. DNV/FFI Spadeadam large outdoor releases

The primary FFI report, **FFI-RAPPORT 20/03101**, was visually reviewed again
on 2026-09-15. Its seven outdoor tests use 10 m and 5 m mast wind sensors,
pad/field thermocouples, and oxygen sensors converted to H2 concentration by
oxygen depletion. The public report presents the outdoor field array at 30,
50 and 100 m, at 0.1, 1.0 and 1.8 m height. It does not describe particle
imaging, a condensed-air mass measurement, velocity profiles, or separate gas
and particle velocities.

This makes the series a valuable independent downstream concentration/
temperature/trajectory test but not an adoption dataset for a two-velocity
condensed-air closure. Reported wind is also insufficient to reconstruct a
resolved release-height velocity field: the document states that only the
high-mast average is included.

### Directly reproducible public screening anchors

The report narrative itself supplies two particularly useful horizontal,
outdoor anchors.  They are retained here as *reported maxima*, not as a
replacement for the unpublished time series or a fitted target:

| test | release | reported wind, high / low mast (m/s) | reported peak H2 at 30 m | reported peak H2 at 50 / 100 m | use |
|---:|---|---:|---:|---:|---|
| 4 | horizontal, 25.4 mm, 0.828 kg/s | 6.7 / 5.0 | 17.2 vol % at 0.1 m | 6.5 / 1.1 vol % | high-wind grounded-cloud lower bound |
| 6 | horizontal, 25.4 mm, 0.833 kg/s; values before ignition | 2.7 / 2.3 | 21.0 vol % at 0.1 m | <= 1.8 / <= 1.8 vol % | low-wind grounded-cloud lower bound |

Those values were checked directly against the primary report on 2026-09-15.
They show the important physical discrimination without inventing a particle
closure: comparable horizontal releases have materially different downstream
behaviour under their distinct wind conditions.  The local validation reader
uses the report's static sensor-table transcription only for a reproducible
comparison and keeps it in the ignored `reference/` tree; neither that table
nor any external raw time series is part of the distributed package.  A
reported arc maximum is lower-censored by the finite sensor array, so it must
not be treated as an observed centreline value or averaged against an
unmatched model timescale.

The same primary-source audit corrected two local, ignored transcription
values used only by the separate downward-release regression: Test 3 is
43.8 kg/min (0.730 kg/s), not 0.630 kg/s; Test 5 is 42.9 kg/min (0.715 kg/s),
not 0.739 kg/s.  The correction changes the deliberately non-promoted
downward `LiftoffPlume` diagnostic from MG/VG `4.79/12.32` to `4.74/11.94`;
its conclusion is unchanged (strong dilution/over-rise bias and no
quantitative validation claim).  No coefficient was changed to absorb this
source-data correction.

## 3. Sandia Hecht--Panda cryogenic jets

The local paper and its reduced Raman data cover near-liquid-temperature
hydrogen jets through 1.0 and 1.25 mm nozzles at pressures up to 5 bar. They
are suitable for checking mean centreline concentration, temperature decay,
radial profile width, concentration RMS and intermittency. Sandia also
describes simultaneous Raman/PIV capability, but a complete public numerical
set of `u`, `v`, `w`, RMS, Reynolds-stress and `epsilon` profiles has not been
located. A 2024 public LES paper supplies mean and axial-RMS velocity contours
for the 5 bar/50 K/1 mm case. Those contours now constrain a TKE realizability
envelope, not experimental model accuracy; see
[`public-turbulence-evidence-audit.md`](public-turbulence-evidence-audit.md).

## 3b. Sandia steady-crosswind LH2 pooling experiments

Hecht's public Sandia report, [*Liquid Hydrogen Pooling and Vaporization
Experiments*](https://doi.org/10.2172/2585576), reports 16 controlled,
unignited LH2 spills on concrete and steel in a large tunnel. It is
particularly valuable because the cross-wind is steady and unidirectional,
unlike the outdoor field series. The public report gives approximately 8--44
L/min LH2 onto approximately 1 m square substrates at 1.7--4.3 mph
cross-wind, plus per-test condition tables, pool-radius summaries, and
downwind concentration contour figures.

The report's Figures 3-27--3-31 have labelled downwind distance, height and
cross-wind coordinates and 0.04-mole-fraction contour increments. They can
therefore be digitised into a *local*, traceable validation reduction without
receiving a private file. The contour is a spatial interpolation of a finite
sensor array, so a reduction must record the figure, panel, axis calibration,
contour level, and at least one-grid/line-width spatial uncertainty. It must
not claim an interpolated contour is an observed centreline or instantaneous
maximum.

This is a pool-vapour source experiment, not an orifice-jet or cryo-compressed
blowdown experiment. It can falsify or validate a future pool/steady-crosswind
source handoff and downstream thermal/species trajectory, but cannot select a
nozzle two-phase discharge law or a finite-TKE dissipation coefficient. The authors
report that co-located thermocouple-derived mole fractions and extractive HyWAM
measurements had excellent temporal agreement, while minor scaling/offset
discrepancies remained. This supports retaining the existing
temperature--composition check as a measurement diagnostic, not treating it
as exact phase-equilibrium calibration. The report says raw data can be
provided on request; it is not redistributed in this repository.

## 4. Evidence classes

### Ready for direct model use

- PRESLHY E3.5 trial-specific pressure, temperature, mass-flow indicators,
  weather, H2 and thermocouple histories.
- D4.8 reconstructed flow bounds for unreliable two-phase meter cases.
- ELVHYS Test-10/Test-11 pressure, temperature, concentration, ventilation
  and clarified sensor/nozzle coordinates.
- Spadeadam 30/50/100 m H2-equivalent concentration and field-temperature
  observations, with their oxygen-depletion and high-mast-wind limitations.
- Sandia mean concentration and temperature profiles.
- Sandia concentration RMS, intermittency and radial PDF-shape evidence.
- Sandia 2025 steady-crosswind LH2 pool condition tables and contour-derived,
  uncertainty-labelled concentration geometry, reduced locally from public
  figures only under the contained-pool criterion.

### Ready only for physical bounds

- ambient RH/dew point for water/air-condensation limits;
- condensed-air equilibrium calculations without particle size or slip;
- the E3.5 catalogue's qualitative observation of condensed air near the
  release/impingements, without particle-resolved channels;
- the public PRESLHY E3.1 80 K cryo-compressed-H2 Excel/BOS package as a
  potential gas-jet thermal/optical comparison: the source pressure/trigger
  and 0.25/0.75/1.75 m axis-temperature channels are now identified, while
  the reported +7 K thermocouple cold-bath bias and H2-tube delay keep it out
  of direct transient or particle-closure scoring;
- literature-only plume-width and centreline-decay trends;
- published LES mean/axial-RMS velocity contours used only as a realizability
  envelope;
- pressure-based ELVHYS source estimates without measured hydrogen flow.

### Still missing

1. Radial mean velocity and three-component RMS or section-integrated
   `rho*u*k` at a matched cryogenic-jet station.
2. Condensed N2/O2 phase fraction, mass-weighted particle-size distribution,
   and gas-particle slip velocity.
3. Independently measured hydrogen mass flow for ELVHYS WP4.2 Tests 10 and 11.
4. Local three-component wind at PRESLHY release height.
5. Trial-specific hydrogen ortho/para composition or liquefier conversion
   history.
6. Particle-resolved mass/size/slip data for the DNV/FFI large releases;
   their public field array does not provide it.

## 5. Implementation consequence

The source-boundary reconciliation and the coefficient-free thermal-moment
observation test are now complete. The latter is numerically converged through
4 m but improves thermal variance while worsening centre amplitude; see
[`thermal-moment-observation-results.md`](thermal-moment-observation-results.md).
Further thermal-width tuning from the same four profiles would be calibration,
not independent validation.

The next DEGALI step no longer has to wait for a private turbulence file. The
public LES axial-RMS contour gives a defensible lower TKE envelope, and the
public Raman fluctuation paper gives an independent scalar-intermittency
check. A matched experimental velocity profile and epsilon are still needed
before a finite-TKE closure can be validated or promoted as default. Until
then the run remains an explicit bounded sensitivity and cannot support a
field-accuracy claim. Likewise, finite condensed-air relaxation requires an
independent particle-size or slip datum rather than a coefficient inferred
only from plume temperature.

The D4.8 source comparison has now been completed; see
[`d48-source-boundary-reconciliation-results.md`](d48-source-boundary-reconciliation-results.md).
The trial-specific pressure-loss values were retained.

The sealed downstream evidence has also been combined without reintegration;
see [`downstream-residual-map-results.md`](downstream-residual-map-results.md).
It rejects resolved mean-kinetic-energy thermalisation as the primary
cold-core explanation and fixes independent thermal-profile transport as the
next bounded implementation target.

The cross-branch claim boundary is now executable in
[`validation-integration.md`](validation-integration.md). It deliberately
keeps reconstruction, qualified field screens, numerical bounds and unmeasured
mixed-phase/TKE states separate rather than producing a composite accuracy
score.

The public Sandia steady-crosswind pool reduction is registered in
[`prereg-sandia-lh2-pool-contour-screen.md`](prereg-sandia-lh2-pool-contour-screen.md).
It uses a one-sided contour-reach falsification gate. Its initial short-reach
signal exposed a section-mean/point-contour comparison error; after the
existing Gaussian axis profile was used for reported centreline concentration,
the Test-8 lower bound is not falsified. A second, contained Test-5 lower
bound does falsify the present pure-vapour pool path. The earlier Test-3 trace
was excluded because the report says its pool ran off the substrate. This does
not embed the external figure reduction or promote a replacement pool model.
