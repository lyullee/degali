# Pre-registration: Houf velocity-width mapping correction

Date: 2026-09-05

## Defect found before calculation

The crosswind independent-energy state stores the scalar/density Gaussian
standard deviations. Its profiles are

```text
rho-rho_a proportional to exp[-r^2/(2 sigma_s^2)]
v-v_a       proportional to exp[-lambda^2 r^2/(2 sigma_s^2)] .
```

Sandia HyRAM+ 6.1 commit
`b45abf9a6d995951311be6aad836f1874e4d420b`, `_jet.py` lines 279--280 and
515--520, defines `B` in the Houf local-Froude and entrainment equations as
the **velocity** e-folding width in `v=V_cl exp(-r^2/B^2)`. Therefore the
elliptical area-preserving mapping is

```text
B_velocity = sqrt(2 sigma_y sigma_z) / lambda,
```

not the previously implemented scalar e-folding width
`sqrt(2 sigma_y sigma_z)`. With `lambda=1.16`, the former value is 16% too
large; before capping, the Houf buoyancy contribution is consequently too
large by `lambda^2=1.3456`, and the pure-plume cap is too large by `lambda`.
This supersedes that one conversion sentence in `prereg-houf-entrainment.md`.

## Frozen implementation

- Add an explicit `houf_width_mapping` choice to the independent-energy
  research model: `"velocity"` uses the corrected equation and is the new
  default; `"scalar"` reproduces the 2026-09-05 historical audit.
- Change no source, profile, thermodynamics, force, coefficient, cap, handoff,
  integration setting or observation operator.
- The five-flux boundary must be identical because the mapping occurs only in
  the downstream entrainment source.
- First run a trial-10 directional check, then the same seven trials, 42 arcs
  and 17 internal and sensor-fit vertical sections at a 0.02 m maximum step,
  using Li equation 35 enthalpy transport.

## Gates and interpretation fixed before calculation

1. A unit test must recover `B_velocity=sqrt(2*sysz)/lambda` and the scalar
   compatibility option must recover `sqrt(2*sysz)`.
2. All seven interfaces, downstream positivity, 42/17 coverage, balance below
   `1e-6`, and sensor-fit R-squared above 0.85 remain mandatory.
3. Physical correctness is independent of validation direction: the velocity
   mapping becomes the research default even if the metrics worsen. The old
   scalar mapping remains only for result reproduction.
4. Production promotion still requires the pre-existing all-metric rule
   against the corrected JETPLU baseline. No coefficient will be retuned after
   seeing the result.

The sign of the field response is not pre-claimed. Reduced entrainment delays
both dilution and ambient heating; for an LH2 cloud crossing from dense to
buoyant, those effects act in competing directions.

## Result

The implementation passes the analytic mapping test and both retained
mapping paths pass the representative conservative interface/integration
test. The corrected 0.02 m field run retains all seven interfaces, 42 arcs,
17 vertical sections and positive states. Its maximum direct balance residual
is `5.07e-8`, and the minimum sensor-profile R-squared is 0.9756.

| model | MG | VG | FAC2 | internal sigma ratio | internal centre MAE, m | sensor-fit sigma ratio | sensor-fit centre MAE, m |
|---|---:|---:|---:|---:|---:|---:|---:|
| corrected JETPLU baseline | 1.074 | 1.213 | 0.952 | 1.091496 | 0.054195 | 1.359686 | 0.069557 |
| Li enthalpy, old scalar-width mapping | 1.118 | 1.176 | 0.976 | 1.065476 | 0.077901 | 1.419304 | 0.136869 |
| Li enthalpy, corrected velocity-width mapping | **1.062** | **1.163** | **0.976** | **0.976421** | 0.075770 | **1.291707** | 0.141674 |

Correcting the equation materially improves all concentration metrics, moves
the internal width to within 2.36% of measurement, and reduces internal centre
MAE by 2.13 mm relative to the old mapping. The sensor-fit width also improves
substantially. The physical correction is therefore retained as the research
default. Production promotion still fails only the centre-height condition:
both the internal and sensor-fit centre MAEs remain above the baseline, and
the sensor-fit MAE is 4.81 mm worse than the old mapping.

Trial 10 at 6 m remains the dominant warning. The corrected internal centre
falls slightly from 0.7742 to 0.7649 m against a measured fit at 0.2266 m, but
the sparse-height direct-plus-image refit returns 1.7621 m despite R-squared
0.9756. This is reported rather than used to retune the physics. The next
audit must separate cumulative vertical momentum/force error from the
ill-conditioned observation fit.

Machine-readable aggregates are in
`reference/preslhy/houf_width_mapping_correction_2026-09-05.json`.
