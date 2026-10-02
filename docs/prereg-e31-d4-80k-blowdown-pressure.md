# Pre-registration: E3.1 4-mm, nominal-80-K blowdown pressure audit

Date frozen: 2026-09-17, before calculating the seven-run comparison.

## Question

Can the explicit adiabatic, well-mixed cryo-compressed hydrogen vessel model
reproduce the pressure-history *envelope* of the seven public 4-mm,
nominal-80-K PRESLHY E3.1 releases without selecting a discharge coefficient,
wall-transfer coefficient, or phase-withdrawal mode from their results?

This is a source-boundary audit. It does not validate atmospheric dispersion,
particle slip, a finite pipe heat-transfer closure, or an LH2 spill source.

## Frozen public population and synchronization

- Use every public `d4_80k_data/*.xlsx` workbook: the nominal 5, 10, 23, 50,
  100, 150 and 200 bar runs. Read each in place; distribute neither workbooks
  nor derived per-sample histories.
- Use its recorded 4-mm diameter and initial vessel-pressure metadata.
- Define physical source start as the first PNoz value at least 2 bar above
  its pre-relay median. The 2-bar criterion is declared before model results.
- Initialize the tank pressure to the *measured vessel pressure at that source
  start*. This is an observed initial condition, not a fitted pressure offset.
- Compare pressure at 0.5, 1.0 and 2.0 s after that same source start.

## Frozen model envelope

- Fixed vessel volume: 2.815 L; tank inside diameter: 160 mm; outlet height:
  30 mm; nozzle diameter: recorded 4 mm; ambient pressure: 101325 Pa.
- Use no tank-wall or discharge-pipe heat transfer. Their geometries and
  coefficients are not identified by the public pressure traces.
- Temperature sensitivity is the uncorrected nominal 80 K plus the documented
  closed-thermocouple cold-bath offset bracket: 77, 80 and 84 K. No member is
  chosen as the actual storage temperature.
- Discharge-coefficient sensitivity is the independently reported
  `Cd = 0.6, 0.7, 0.8` range. No per-run selection is permitted.
- Continue the two-phase tank state separately with (a) equilibrium vapour
  withdrawal limited by outlet liquid level and (b) homogeneous-equilibrium
  common-velocity withdrawal. These are labeled end members, not fitted modes.
- Integrate to 2.0 s with 0.01-s state steps. Record termination, discrete
  energy residual, pressure and temperature at each station.

## Decision rules

- Report every 3 x 3 x 2 end-member result. Do not collapse them to an RMS or
  choose the closest case.
- At each station report whether the measured pressure falls in the full
  declared envelope, only as a consistency screen. Coverage is not an
  accuracy score because the members omit pipe/wall heat transfer and a real
  stratified withdrawal process.
- If a model stops before a station, retain the stop and do not extrapolate or
  substitute a different phase mode.
- No result changes a DEGALI default. A failure identifies the missing source
  boundary measurement or model physics; it does not authorize a fitted UA or
  `Cd`.
