# Pre-registration: PRESLHY subcooled-liquid source state

Date: 2026-09-05

## Question

The current release builder infers liquid temperature from tanker gauge
pressure as though the tanker contained saturated LH2 at that pressure.  At
5 barg this gives 28.27 K.  D3.6 instead describes pressure-driven liquid
near the LH2 boiling temperature in the measurement pipe, with low intrinsic
vapour pressure and little adiabatic flash.  Driving pressure and liquid
thermodynamic temperature are therefore independent measured source
variables.  This audit tests the coefficient-free lower-enthalpy bound before
adding any finite-rate droplet parameter.

## Frozen implementation and cases

- Add an optional explicit liquid storage temperature to `hydrogen_jet`.
  Omission must reproduce the existing tanker-pressure saturation assumption
  exactly.
- Add a named `ambient_boiling` source mode to the reduced PRESLHY entry
  points.  It uses CoolProp saturated-liquid temperature at the recorded
  ambient pressure (about 20.37 K), matching D3.6's corrected boiling-liquid
  interpretation.  Tanker pressure remains the mechanical driving-pressure
  input; it is not replaced by ambient pressure.
- Compare `tank_saturation` and `ambient_boiling` on the frozen seven trials
  10, 11, 12, 22, 23, 24 and 25, using sustained-window mean flow.
- Keep source-flux Gaussian establishment, Li enthalpy transport, equilibrium
  dry-air phase closure, Houf velocity-width mapping, local-shear crosswind
  entrainment, free ground geometry and 0.02 m maximum step unchanged.
- Recompute the 42 raw thermocouple points using the pre-registered minimum,
  5th-percentile and temporal-median diagnostics.  Centreline results remain
  a required separate report so off-centre ambient sensors cannot carry the
  thermal conclusion.
- No source temperature, flow multiplier, phase coefficient, sensor subset or
  profile parameter may be fitted to the field measurements.

## Frozen decision rule

The subcooled source can replace the present source assumption only if:

1. all seven five-flux interfaces and downstream balance residuals pass the
   existing limits;
2. concentration absolute log-MG decreases, with VG no larger and FAC2 no
   smaller;
3. centre-height MAE does not increase and vertical-width ratio moves no
   farther from unity;
4. both all-sensor and centreline minimum-temperature median absolute error
   decrease, without reversing either trial's signed centreline bias beyond
   10 K; and
5. the maximum downstream balance residual remains below `1e-6`.

Even if all numerical gates pass, promotion is provisional until the
per-trial corrected nozzle-liquid temperature or vapour quality can be
reconstructed from the pipework channels.  Failure remains useful: it bounds
how much the present pressure/temperature conflation can contribute before a
finite-rate two-phase state is introduced.

## Result

The candidate stopped at the source-applicability gate; no field integration
or temperature scoring was run.  At 101325 Pa the explicit boiling-liquid
temperature is 20.3689 K, the equilibrium flash fraction is zero and the
liquid density is 70.8483 kg/m3.  With the existing pressure-thrust-free
source reconstruction this produces the following established axisymmetric
states:

| trial | diameter (m) | velocity (m/s) | local wind (m/s) | velocity/wind |
|---:|---:|---:|---:|---:|
| 10 | 0.31731 | 1.9809 | 2.0966 | 0.9448 |
| 11 | 0.14986 | 11.9548 | 2.2949 | 5.2092 |
| 12 | 0.07486 | 13.5225 | 2.2949 | 5.8924 |
| 22 | 0.31707 | 1.7067 | 1.7000 | 1.0039 |
| 23 | 0.14988 | 5.6761 | 1.7000 | 3.3389 |
| 24 | 0.07491 | 9.7485 | 1.9000 | 5.1308 |
| 25 | 0.31699 | 1.7187 | 1.8500 | 0.9291 |

All seven are below the frozen momentum-dominated threshold of 10.  The
candidate is therefore rejected, but the low ratios are not evidence that the
experiments themselves were weak jets.  They show that lowering the liquid
temperature without reconstructing the nozzle converts the measured flow
through a high liquid density while discarding the pressure work that must
become velocity and/or pressure thrust.

For trial 23, the old 5-barg saturation assumption gives 28.2674 K,
23.4760% equilibrium flash and 5.3469 kg/m3 orifice density; the boiling-
liquid bound gives 20.3689 K, zero flash and 70.8483 kg/m3.  That order-of-
magnitude density change makes a temperature-only substitution mechanically
inadmissible.  The default remains unchanged.

The next source model must close mass, momentum and energy from independent
pipe/nozzle measurements: static pressure, corrected temperature or vapour
quality, two-phase density (or its calibration method), and effective area or
discharge coefficient.  Trials 10 and 23 already contain nozzle-pressure and
calculated-density channels in their original workbooks.  The remaining five
workbooks are to be retrieved individually from the public TAR by HTTP byte
range, rather than downloading the 11.3 GB archive.  The complete numerical
record is `reference/preslhy/subcooled_liquid_source_screen_2026-09-05.json`.
