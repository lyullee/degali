# Data request: remaining LH2 trajectory states

Date: 2026-09-05

## Purpose

The remaining PRESLHY error is excessive buoyancy accumulation between 10D
and 6 m. Source-flux bookkeeping, equilibrium versus delayed dry-air
condensation, interface buoyancy, ground contact and total-energy versus
enthalpy switches have now been separated. None closes the trajectory error.
The next defensible model needs either a transported turbulent-energy/profile
state or finite-relaxation condensed-air momentum. The current concentration
and temperature peaks cannot identify either one.

## Priority 1 — turbulent-energy/profile boundary

For PRESLHY E3.5 trials **10, 11, 12, 22, 23, 24 and 25**, request any
available measurements or CFD reductions at Station 3, Station 4, 10D, and
the 0.79, 1.78, 4 and 6 m arrays:

1. Mean axial, radial and vertical gas velocity, `u`, `v`, `w` in m/s, versus
   radius/height; include the coordinate origin and whether velocity is gas,
   droplet, mixture or tracer velocity.
2. RMS fluctuations `u_rms`, `v_rms`, `w_rms` in m/s and Reynolds stresses
   `u'v'`, `u'w'`, `v'w'` in m2/s2. If only turbulent kinetic energy is
   available, request `k` in m2/s2 and its averaging definition.
3. Dissipation rate `epsilon` in m2/s3, or the integral length/time scale used
   to infer it.
4. Synchronized H2 mole fraction, temperature and velocity profiles with raw
   sample rate, averaging window, calibration uncertainty and channel time
   offset. A profile averaged at a different time from the release flow is not
   sufficient.
5. If only CFD is available, request the area integrals at each section:
   `integral rho*u dA`, `integral rho*Y_H2*u dA`,
   `integral rho*u^2 dA`, `integral rho*u*k dA`, and
   `integral rho*u*(h-h_amb) dA`, in SI units, plus mesh/time-averaging
   uncertainty.

Minimum useful subset: radial mean velocity plus one of three-component RMS
velocity or section-integrated `rho*u*k` at 10D for trials 10 and 23 (different
release height/wind). Without this, a turbulent-energy relaxation coefficient
would be fitted to plume height and is not identifiable.

## Priority 2 — condensed-air finite relaxation

For the same source hardware, or a geometrically and thermodynamically matched
LH2 jet, request:

1. Condensed N2 and O2 mass flow or mass fraction separately, kg/s or kg/kg,
   versus axial station and time.
2. Particle/droplet size distribution as number- and mass-weighted
   `D10/D50/D90` in micrometres, with instrument lower detection limit.
3. Gas velocity and particle velocity measured independently in m/s; at
   minimum axial slip `u_particle-u_gas` and vertical settling velocity.
4. Particle and gas temperature separately, K, or a documented assumption of
   thermal equilibrium with its response-time estimate.
5. Deposition/rainout mass collected on the pad versus distance, split by
   water frost and dry-air condensate if chemically analysed.
6. Ambient RH/dew point and CO2 mole fraction during each release, since water
   and CO2 solids change both particle mass and latent heat.

Minimum useful subset: mass-weighted `D50` with a credible 10--90% range,
condensed-air mass fraction at one station, and axial slip velocity. A total
visible-cloud image alone cannot distinguish water frost from N2/O2 condensate
and is not a usable calibration target.

## Priority 3 — source and coordinate closure

For every requested trial, include:

- hydrogen mass-flow time series in kg/s at the highest available sample rate;
- absolute upstream/nozzle pressure in Pa and fluid temperature in K;
- effective flow area in m2 or internal diameter in m, discharge coefficient,
  nozzle length and edge geometry;
- vapour quality/void fraction at the flow meter and nozzle, or the method used
  to classify liquid/two-phase/gas flow;
- release-axis elevation above the pad and the exact axial origin of every
  sensor array;
- local three-component wind at the release height, not only a 10 m mast mean.
- released-hydrogen ortho/para mole fractions. If they were not measured,
  provide the liquefier conversion specification/certificate, conversion
  catalyst and temperature, delivery date, tank fill and release timestamps,
  boil-off history, and whether hydrogen contacted a conversion catalyst
  after liquefaction. A label such as "LH2" alone does not identify the spin
  composition used by the caloric equation of state.

The spin item is quantitatively material: the pure-para end member changes
the seven-trial fast-bound centre-height MAE by 15.5 mm, despite changing
measured-nozzle density and velocity by less than 0.21%. Until the composition
is supplied, normal and pure-para hydrogen remain reported uncertainty bounds.

For residual liquid-H2 droplets, the minimum useful distribution is
mass-weighted D10/D50/D90 at the nozzle and one downstream station, in
micrometres, with the instrument cutoff and sampling-bias correction. The
current GASFLOW screen requires about 41--76 um droplets to retain phase lag
over the collective source transit, versus 0.78--2.01 um from the fixed TNO
correlation; a number-weighted mean without the large-droplet mass tail cannot
resolve this question.

For ELVHYS Tests 10 and 11 specifically, the minimum boundary remains hydrogen
mass flow versus time plus the resolution of whether the published 200/250 mm
coordinate is nozzle elevation, sensor elevation or a separate origin.

## Delivery format

Preferred: CSV or HDF5 with one channel per named column and a separate data
dictionary containing SI units, clock origin, sample rate, missing-value code,
sensor coordinates and uncertainty. Native TDMS/MAT files are acceptable if
the channel map is included. PDF plots without underlying values are useful
for bounds but not for parameter identification.
