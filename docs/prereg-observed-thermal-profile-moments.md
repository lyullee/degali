# Pre-registration: observed thermal/species profile moments

Date frozen: 2026-09-08, before computing the paired profile moments.

## Question

The model already transports an independent enthalpy-profile moment, but its
downstream thermal/species diffusivity ratio remains explicit.  Can the public
PRESLHY E3.5 measurements constrain the *direction and scale* of differential
thermal/species spreading without fitting a Gaussian profile or a model
coefficient?

This is an observation-operator audit.  It does not identify a turbulent
Prandtl or Schmidt number and cannot promote a downstream closure by itself.

## Frozen source and scope

- Source dataset: Lyons, Coldrick and Atkinson (2023), PRESLHY E3.5,
  [DOI 10.35097/1481](https://doi.org/10.35097/1481), CC BY-SA 4.0.
- Read the original Trial 10 and Trial 23 workbooks in place.  Record their
  SHA-256 hashes and never copy them into the distributed package.
- Use `Flexlogger` converted near-field thermocouple columns ending in `C`,
  its `Ambient_Temperature` column, and `Xensor` hydrogen-output columns.
- Keep the previously frozen inclusive Flexlogger rows 30--78 for Trial 10
  and 28--125 for Trial 23.  These are the mass-flow-selected release windows.
- Use only the two complete cross-shaped stations, x=1.78 and 4.00 m, and only
  the five vertical points at y=0 and z-axis=-0.50,-0.25,0,+0.25,+0.50 m.
  Temperature and hydrogen serials differ, so pair them by the D3.6 Table A3
  physical coordinates, not by channel names.

## Frozen reduction

1. Convert both sheet clocks to seconds after midnight.  The Flexlogger points
   are the target 1 Hz clock.  Linearly interpolate Xensor output to those
   target times only when bracketed by raw Xsensor points; never extrapolate.
2. At each section and time define nonnegative measured scalars

       theta_i = max(T_ambient - T_i, 0)
       c_i     = clip(X_H2,i, 0, 100)

   The clipping is an observation convention for small electrical noise, not
   a model state correction.  Count every clipped value.
3. On the measured vertical interval use trapezoidal quadrature to calculate

       M0 = integral phi dz
       zbar = integral z phi dz / M0
       sigma_trunc^2 = integral (z-zbar)^2 phi dz / M0

   separately for theta and c.  No Gaussian extrapolation outside the one-metre
   instrument span is allowed.
4. A scalar profile is usable only when M0>0, its centre-point value is at
   least 5 K for theta or 1 vol-% for c, and all five coordinates are present.
   These are signal gates, not adjustable accuracy criteria.
5. Estimate the integer 0--5 s downstream delay independently from the
   centre-point theta and c series by the lag having the largest Pearson
   correlation.  A common parcel-pair screen passes only when the two selected
   lags differ by at most 1 s and both correlations are at least 0.5.
6. Pair upstream moments at time t with downstream moments at t+the rounded
   mean of the two lags.  Report every usable pair, the thermal/species
   truncated-width ratio at each section, and positive squared-width growth.

## Fixed interpretation

- At least ten common usable time pairs are required for a trial-level result.
- If the common-delay screen fails, report the station moments but issue no
  differential-spreading inference.
- If squared width decreases or changes sign, retain it.  Do not take an
  absolute value or discard it when summarising growth signs.
- A positive median growth ratio

      R = (sigma_theta,4m^2 - sigma_theta,1.78m^2)
          / (sigma_c,4m^2 - sigma_c,1.78m^2)

  is only a finite-window temperature-deficit/species proxy.  It is not
  `D_h/D_C`, because density, enthalpy nonlinearity, phase change, lateral
  truncation and unsteady advection remain unresolved.
- Agreement with the independently published `Sc_t/Pr_t=0.70/0.85` range may
  support carrying that value as a research sensitivity.  Disagreement or
  poor common-delay evidence leaves the ratio explicit.  No value from this
  audit becomes a default and no field trajectory is re-integrated.
