# LH2 data readiness for the next model stage

Last reviewed: 2026-09-08

This register separates data that can constrain the model from material that
is useful only for context or bounds. Third-party raw files are kept under the
ignored `reference/` tree and are not distributed with DEGALI.

## Decision summary

| question | usable evidence now | readiness |
|---|---|---|
| PRESLHY E3.5 source mass flow | raw meter histories plus the D4.8 pressure-drop reconstruction | usable with trial-specific quality flags |
| PRESLHY E3.5 nozzle state | pipe/nozzle pressure and nozzle-flow/wall thermocouples | usable, subject to two-phase sensor interpretation |
| PRESLHY E3.5 atmospheric boundary | ambient temperature, RH, wind speed and direction | usable at the reported reference height; no local 3-component wind |
| PRESLHY E3.5 plume temperature and H2 | thermocouple arrays, Xensor and Drager channels | usable with alignment, saturation and coordinate caveats |
| transported turbulent energy | no `u/v/w`, velocity RMS, Reynolds stress, `k`, or `epsilon` in the E3.5 workbooks | not identifiable from current data |
| condensed-air relaxation | literature establishes condensation; no matched particle size, phase fraction, or slip history | bounded only |
| ELVHYS Tests 10 and 11 geometry | archive dictionary separates the nozzle at 250 mm from bottom samplers at 200 mm | usable for archive reduction |
| ELVHYS Tests 10 and 11 hydrogen flow | `FLMT` is ventilation-air flow only | missing |
| Sandia cryogenic free jets | Raman concentration and temperature profiles are available in the paper/local reduction | usable for mean thermal/species profiles, not a complete TKE boundary |

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

## 3. Sandia Hecht--Panda cryogenic jets

The local paper and its reduced Raman data cover near-liquid-temperature
hydrogen jets through 1.0 and 1.25 mm nozzles at pressures up to 5 bar. They
are suitable for checking mean centreline concentration, temperature decay,
and radial profile width. Sandia also describes simultaneous Raman/PIV
capability, but a complete public numerical set of `u`, `v`, `w`, RMS,
Reynolds-stress and `epsilon` profiles has not been located.

## 4. Evidence classes

### Ready for direct model use

- PRESLHY E3.5 trial-specific pressure, temperature, mass-flow indicators,
  weather, H2 and thermocouple histories.
- D4.8 reconstructed flow bounds for unreliable two-phase meter cases.
- ELVHYS Test-10/Test-11 pressure, temperature, concentration, ventilation
  and clarified sensor/nozzle coordinates.
- Sandia mean concentration and temperature profiles.

### Ready only for physical bounds

- ambient RH/dew point for water/air-condensation limits;
- condensed-air equilibrium calculations without particle size or slip;
- literature-only plume-width and centreline-decay trends;
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

## 5. Implementation consequence

The next DEGALI step can proceed without another bulk download. The source
boundary has been reconciled with the D4.8 flow hierarchy, so the next new
work is the downstream temperature, density, buoyancy and energy residual
analysis for PRESLHY Trials 10 and 23. A transported-TKE or finite-relaxation parameter
must not be fitted until an independent velocity/TKE or particle-slip datum is
obtained. Until then, such physics should be implemented as a bounded
sensitivity branch and reported as structurally unvalidated.

The D4.8 source comparison has now been completed; see
[`d48-source-boundary-reconciliation-results.md`](d48-source-boundary-reconciliation-results.md).
The trial-specific pressure-loss values were retained.

The sealed downstream evidence has also been combined without reintegration;
see [`downstream-residual-map-results.md`](downstream-residual-map-results.md).
It rejects resolved mean-kinetic-energy thermalisation as the primary
cold-core explanation and fixes independent thermal-profile transport as the
next bounded implementation target.
