# Pre-registration: established-flow energy-partition audit

Date: 2026-09-05

## Competing published equations

The independent-energy model currently follows HyRAM+ 6.1 equation 82 and
transports the ambient-relative total enthalpy-plus-mean-kinetic-energy flux,

```text
Q_total = integral rho v (h - h_amb + v^2/2) dA.
```

Li et al. (2026), equation 35, use a different established-flow closure for
the same seven variables and Gaussian density/species profiles:

```text
Q_h = integral rho v (h - h_amb) dA.
```

The distinction is a physical partition, not a numerical method. The total-
energy closure converts loss of resolved mean kinetic-energy flux into the
represented equilibrium thermal state. The Li closure transports mixture
enthalpy independently, leaving that mean-flow loss outside the resolved
thermal state. The latter can therefore remain colder and denser downstream,
which is a directional candidate for the excessive rise of the completed
total-energy calculation.

## Frozen implementation

- Add an explicit `energy_transport` choice, `"total"` or `"enthalpy"`, to
  the independent-energy research model. The default remains `"total"` so
  the already recorded result is reproducible.
- For `"enthalpy"`, remove only `v^2/2` from the fifth section flux and set
  its entrained-ambient source to zero because the transported enthalpy is
  relative to ambient. Mass, H2, vector momentum, trajectory, Gaussian
  profiles, phase equilibrium, all entrainment terms and coefficients remain
  unchanged.
- Construct the 10D target with the matching enthalpy-only integral of the
  existing conserved near-field profile. Do not approximate it from centre
  temperature.
- Use direct conserved-flux RK4 at the already converged 0.02 m maximum step.
- Use the same seven trials, 42 concentration arcs and 17 vertical fits.
- Report both internal-state geometry and the pre-registered sensor-height
  direct-plus-image Gaussian observation operator. The latter is the primary
  comparison with the fitted PRESLHY geometry.

## Gates and decision

1. All seven 10D interfaces must pass the existing `1e-8` five-flux,
   `1e-5` quadrature, 5% H2-width and 2 K temperature gates.
2. Every downstream run must remain positive and the maximum direct selected-
   flux balance residual must remain below `1e-6`.
3. All 42 concentration arcs and 17 vertical sections must remain present;
   the minimum sensor-profile R-squared must remain at least 0.85.
4. The Li closure replaces total energy as the independent-energy research
   base only if no concentration metric worsens (`abs(log MG)`, VG, FAC2),
   and both the sensor-fit width-ratio error and centre MAE improve relative
   to the total-energy candidate.
5. Production promotion remains harder: the same all-metric rule must also
   beat the corrected JETPLU baseline. No tolerance or coefficient will be
   changed after inspecting the result.

The expected direction is lower centre temperature or slower warming, higher
centre density and less rise. The magnitude and the concentration/width
trade-off are deliberately not predicted.

## Result

Both published closures pass all seven 10D interfaces. For the Li enthalpy
closure the maximum interface flux residual is `1.13e-15`, quadrature change
`8.42e-6`, H2-width mismatch 4.762% and centre-temperature mismatch 1.941 K.
All seven downstream integrations remain positive and the maximum direct
enthalpy-flux balance residual is `5.18e-8`.

The expected direction is confirmed but the magnitude is negligible:

| independent-energy candidate | HyRAM+ total energy | Li enthalpy |
|---|---:|---:|
| MG / VG / FAC2 | 1.118 / 1.176 / 0.976 | 1.118 / 1.176 / 0.976 |
| internal mean `sigma_z` ratio | 1.065947 | 1.065476 |
| internal centre MAE, m | 0.077939 | 0.077901 |
| sensor-fit mean `sigma_z` ratio | 1.420033 | 1.419304 |
| sensor-fit centre MAE, m | 0.137151 | 0.136869 |
| minimum sensor-profile R-squared | 0.971273 | 0.971219 |

The Li branch satisfies the frozen research-continuation direction at reported
precision and is retained as the published alternative for the next audit.
It is not a production promotion: its concentration bias and both geometry
errors remain worse than the corrected JETPLU baseline. The sub-millimetre
change in aggregate centre MAE also falsifies mean-kinetic-energy thermalisation
as the cause of the excessive rise. The HyRAM+ total-energy branch remains the
reproducibility control and the API default.
