# TI021 600 s dose root-cause diagnostics — 2026-10-15

## Purpose

The acceptance metric remains the 600 s cumulative dose MAPE on the fixed
`ti21_input_ascii.xls` workbook. The historical baseline is 84.901%. Short
prefixes and altered domains are diagnostics only; they are not promoted as
the final model unless the full dose score improves and the conservation and
route/height checks remain acceptable.

## Evidence from full 600 s runs

The historical 2.5 m local patch with analytic exterior packet continuation
gave 79.885% dose MAPE. Its 84 positive-observed sensors split into 33 sensors
inside the layer patch and 51 sensors on the exterior packet route. The
predicted/observed dose sums were:

| route | positive sensors | observed dose sum (vol% s) | predicted dose sum (vol% s) | ratio |
|---|---:|---:|---:|---:|
| layer patch | 33 | 13,796.749 | 1,886.879 | 0.1368 |
| exterior packets | 51 | 3,520.380 | 993.308 | 0.2822 |

The species residual was about `1.25e-15` relative and the unmapped species
fraction was 0.786%. The dominant problem is therefore spatial transport and
receptor exposure, not disappearance of conserved species.

### Sensor-envelope domain diagnostic

The solver was extended with `--domain-mode sensor_envelope`. This puts the
full 145-sensor footprint and the source/fence neighbourhood inside the
finite-volume patch, removing the analytic exterior route without changing
the source or layer coefficients.

| case | cells | dose MAPE | runtime |
|---|---:|---:|---:|
| local patch, 5 m | 44 x 56 | 79.885% | 4.5 s |
| sensor envelope, 5 m | 110 x 120 | 75.367% | 11.7 s |
| sensor envelope, 2.5 m | 220 x 240 | 77.657% | 63.9 s |

The 5 m envelope improvement of 4.52 percentage points confirms that the
external continuation is a material error source. The 2.5 m envelope is worse
than the 5 m envelope, so uniform grid refinement is not a sufficient fix and
the single-layer closure is not numerically converged.

The 5 m envelope still underpredicts near sensors: S001 is 184.74 versus
1,832 vol% s observed and S002 is 46.71 versus 364.936 vol% s. Thus external
coverage is only one part of the error.

### Vertical observation diagnostic

The scorer was then connected to the existing diffusive vertical-tail
operator. This operator is deliberately diagnostic: it changes the receptor
observation of unresolved vertical exchange without adding conserved species
to the layer state.

| domain | vertical operator | dose MAPE | low-height MAPE | mid-height MAPE | upper-height MAPE |
|---|---|---:|---:|---:|---:|
| local, 5 m | top-hat | 79.885% | 64.45% | 81.95% | 90.59% |
| envelope, 5 m | top-hat | 75.367% | 56.80% | 68.66% | 91.90% |
| local, 5 m | diffusive tail, 10 m, D=1 m²/s | 66.581% | 62.02% | 68.85% | 63.98% |
| envelope, 5 m | diffusive tail, 10 m, D=1 m²/s | 46.848% | 53.97% | 47.11% | 38.25% |

When the same tail is normalized so that the resolved top-hat plus tail has
the original vertical integral, the envelope score returns to **81.169%**.
This rejects the unnormalized 46.848% result as a production correction: its
improvement came partly from exposing extra concentration above the layer
without transferring a corresponding conserved vertical state. The
mass-conserving version passes the species residual gate but is not yet a
useful accuracy improvement. A real multilayer state is still required.

The large response confirms that missing height-dependent exposure is a major
error source. It is not yet a production correction: the tail is an
observation operator rather than a conserved multilayer transport solution.
The next implementation must transfer the same vertical species profile into
multiple conserved layers and repeat the 600 s score. The 46.848% value is
therefore a diagnostic lower bound/target, not a validated accuracy claim.

### Conserved multilayer transport

The layer solver now has an opt-in conservative split transport with vertical
species/mixture/momentum layers. Horizontal FV transport is applied to the
column, the column update is distributed over the existing vertical mass
fractions, ambient entrainment enters the ground layer, and a zero-flux
vertical diffusion step redistributes all four conserved components. The
receptor reads the layer state rather than applying an extra tail.

| case | vertical layers | dose MAPE | species residual | runtime |
|---|---:|---:|---:|---:|
| envelope, 5 m, top-hat | 1 | 75.367% | 1.2e-15 relative | 11.7 s |
| envelope, 5 m, conservative multilayer | 5 | **53.283%** | 1.1e-15 relative | 22.1 s |
| envelope, 5 m, conservative multilayer | 10 | 53.607% | 3.5e-15 relative | 39.5 s |

Adding a bounded directional separated-flow closure (50 m wake length, 10 m
base width, 0.2 spread slope, 0.5 recirculation strength) to the same
conservative multilayer case gives **52.441%** at 5 layers.  The 10-layer
rerun gives **52.751%**, a 0.310 percentage-point change, and is therefore a
stable practical candidate rather than a one-layer numerical artifact.

The no-wake 5-versus-10-layer change is 0.325 percentage points, which is a useful
initial convergence check. For the accepted recirculation candidate, the
corresponding change is 0.310 percentage points. Height-group MAPE for the
accepted 5-layer result is
54.33% (0--1 m), 46.32% (1--3 m), and 54.84% (3--6 m), all lower than the
79.885% local baseline groups. This is the first structural candidate that
improves the full dose score without arbitrary source scaling and while
passing the mass gates. It remains experimental: the vertical split has no
vertical momentum equation or resolved wake, and the outer domain still uses
the envelope rather than a recoupled density-current exterior.

Using the historical local-route partition (S001--S033 in the resolved patch,
S034 onward in the exterior packet route), the same candidate also improves
both route groups: local-patch MAPE changes from 85.03% to 63.61%, and the
exterior-route group changes from 76.56% to 45.22%. Thus the aggregate gain is
not caused by sacrificing one route for the other.

The latest reruns also pass the explicit momentum ledger.  For the 5-layer
case the maximum momentum-budget residual is `9.31e-10` kg m/s and
the maximum vertical-mixing residual is `4.66e-10` kg m/s; the
10-layer values are `9.31e-10` and `4.66e-10` kg m/s.
Both runs report `momentum_conservation_gate=true` alongside the species and
mixture mass gates.

An isotropic conserved horizontal diffusivity of 20 m²/s was also tested as a
cheap wake-mixing surrogate. It worsened the multilayer score from 53.283% to
68.526%, despite passing the mass residual gate. This rejects uniform scalar
diffusion as a substitute for an obstacle wake; a directional separated-flow
or recirculation closure is still required.

A directional lee-wake velocity-deficit diagnostic (strength 0.5, length 50 m,
base width 10 m, spread slope 0.2) gave **53.964%**, slightly worse than the
53.283% no-wake multilayer result. It is therefore not promoted either; a
prescribed speed deficit is not enough to reproduce the missing 3-D wake and
recirculation structure.

When the same closure includes a bounded reverse-flow recirculation component,
the sensor-envelope 5-layer run gives **52.441%** (strength 0.5); every route
and height group improves relative to the no-wake 53.283% candidate.  Its
momentum and mass gates pass.  This is accepted as the session's best
experimental closure, while remaining explicitly unqualified as full 3-D CFD.

### Exterior density-current continuation diagnostic

The exported boundary packets now carry the resolved hydrostatic gravity-current
speed.  A separate exterior continuation mode projects that speed toward each
receptor and adds a conservative shallow-current spreading term; the historical
ambient-wind route remains unchanged for the reference comparison.  This is a
diagnostic closure only: it is not a recoupled outer grid and does not invent or
renormalize packet mass.

| local 5 m, 5-layer case | exterior route | density-current strength | dose MAPE |
|---|---|---:|---:|
| multilayer reference | ambient wind | 0.00 | **71.177%** |
| multilayer diagnostic | density current | 0.25 | 73.408% |
| multilayer diagnostic | density current | 0.50 | 73.281% |
| multilayer diagnostic | density current | 0.75 | 73.154% |
| multilayer diagnostic | density current | 1.00 | 73.027% |

All cases retained the same species residual (about 2.4e-15 relative) and
mixture residual (about 1e-9 kg), but every density-current strength worsened
the full 600 s dose score.  The route is therefore rejected as a production
correction.  The result is useful evidence: the remaining exterior error is not
fixed by changing packet travel speed alone; recoupled 3-D bypass, separation,
and wake re-entry are still missing.

### Source inventory sensitivity diagnostic

`--source-mass-scale` was added as an explicitly diagnostic-only perturbation.
It scales the source prefix and handoff inventory consistently; it does not
change the production default and must not be interpreted as calibration.

| source mass scale | dose MAPE |
|---:|---:|
| 1.0 (historical local baseline) | 79.885% |
| 1.5 | 69.173% |
| 2.0 | 63.907% |
| 3.0 | 65.391% |

This non-monotone response shows that source inventory is a major sensitivity,
but it does not identify a defensible physical mass. The workbook's 4,962 kg/s
for 1 s is a rate-times-duration proxy, while the current initial cloud is
about 3,290 kg. The measured release mass-flow history must decide this input;
the factor 2.0 result is not promoted.

Using the workbook proxy ratio 1.5091 together with the sensor envelope and a
mass-conserving vertical tail gave **71.661%** MAPE. It is also not promoted:
the proxy is not a verified atmospheric source history, and the conservative
tail still lacks a transported multilayer momentum state.

### Coefficient ablations rejected

At the 2.5 m local domain, Parker bulk suppression produced 81.923% MAPE and
uniform wind produced 82.062%, both worse than the 79.885% baseline. Further
coefficient tuning is therefore not the current priority.

## Current root-cause ranking

1. Source inventory/time history is underidentified and strongly affects dose.
2. Analytic exterior continuation underrepresents long-range receptor exposure;
   a speed-only density-current continuation did not cure it.
3. The one-layer vertical profile and missing 3-D obstacle wake/re-entry bias
   near- and multi-height sensor dose.
4. Grid refinement alone does not cure the structural closure error.

## Temporal residual decomposition

For the two full-event NPZ traces, the model arrival operator was defined as
the first time the trace reached 5% of its own local peak. This is comparable
to the workbook's documented arrival convention, but remains diagnostic because
the workbook does not contain the underlying concentration traces.

| case | arrival MAE | arrival bias | departure MAE | departure bias |
|---|---:|---:|---:|---:|
| 52.441% retained candidate | 31.77 s | -27.50 s | 121.24 s | -121.24 s |
| 44.560% earlier hand-off candidate | **28.57 s** | **-23.76 s** | **106.87 s** | **-103.32 s** |

The new candidate improves arrival timing at every reported height group
(0.4/2.4/4.4/6.4 m: 31.62/29.07/23.28/29.32 s versus
38.16/30.65/24.53/30.94 s). Departure remains systematically early, which
identifies missing residence/recirculation duration rather than a simple source
arrival error. No temporal values were used to alter the dose score.

## Next implementation gate

The next physical change must be a conservative exterior density-current
transport/re-entry closure or a domain-coupled outer patch, with explicit
height-resolved receptor sampling. The current 44.560% candidate beats the
52.441% predecessor and passes the active species, momentum, route, and
height-group gates, but it remains experimental until the early departure bias
is addressed with a recoupled wake residence model.

