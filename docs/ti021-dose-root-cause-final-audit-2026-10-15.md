# TI021 dose root-cause and repair audit

## Final criterion

The authoritative metric is the cumulative 600 s dose MAPE over the 84
positive-observed receptors. Arrival time and 30 s dose are diagnostic
indicators only. A candidate is promotable only when its full-event MAPE is
below the reproducible baseline and its conservative mass/momentum gates and
route/height checks do not regress.

## Previously retained candidate

Artifact:
`outputs/integral-local-patch-2026-10-15/dose-envelope5m-multilayer5-wake-recirculation05-v2.json`

| item | result |
|---|---:|
| historical reference MAPE | 84.901% |
| reproducible local baseline MAPE | 79.885% |
| accepted candidate MAPE | **52.441%** |
| grid | 5 m, 110 x 120 sensor envelope |
| vertical state | 5 conservative layers, (D_z=1) m²/s |
| runtime | 22.1 s |
| species relative residual | 2.49e-15 |
| mixture mass residual | -6.52e-9 kg |
| maximum momentum budget residual | 9.31e-10 kg m/s |
| maximum vertical momentum residual | 4.66e-10 kg m/s |
| momentum gate | pass |

The 10-layer convergence check with the same wake closure gives 52.751% MAPE
and changes the score by 0.310 percentage points, so the 5-layer result is
retained as the practical candidate. Historical route partitions both improve:
local-patch MAPE 85.03% -> 63.61%, exterior-route group 76.56% -> 45.22%.
Height groups (0--1 m, 1--3 m, 3--6 m) are 54.33%, 46.32%, and 54.84%, each
below the 79.885% baseline group score.

## New lower-error candidate under the same acceptance gates

The remaining-physics sweep added an earlier source-to-layer overlap and a
declared vertical observation tail while preserving the workbook wind speed.
The full 600 s rerun is:

`outputs/integral-local-patch-2026-10-15/dose-envelope5m-physically-declared-h40.json`

| item | result |
|---|---:|
| full-event MAPE | **44.560%** |
| runtime including output | 14.22 s |
| hand-off radius | 40 m |
| workbook wind speed | 3.4 m/s |
| vertical observation | diffusive tail, 10 m, 2 m²/s |
| horizontal diffusion | 1 m²/s |
| drag / entrainment | 0.05 / 0.10 |
| species relative residual | 6.92e-16 |
| mixture mass residual | 8.38e-9 kg |
| maximum momentum residual | 1.40e-9 kg m/s |
| momentum gate | pass |

The candidate is lower than 52.441% without using observed dose values during
the run. All 84 positive receptors remain on the same `layer_patch` route.
The height-group MAPEs (0--1 m, 1--3 m, 3--6 m) are 52.44%, 45.95%, and
37.18%, respectively, so none regresses against the retained candidate's
54.33%, 46.32%, and 54.84%. This candidate is therefore the current accepted
diagnostic candidate for the active root-cause goal. It remains experimental:
the diffusive tail is an observation closure, not a resolved 3-D vertical
momentum field.

## Root-cause candidates and decisions

| candidate | 600 s result | decision |
|---|---:|---|
| source scale 1.5 | 69.173% | reject: unverified mass history |
| source scale 2.0 | 63.907% | reject: diagnostic only, non-physical calibration |
| source scale 3.0 | 65.391% | reject: non-monotone sensitivity |
| isotropic horizontal diffusion 20 m²/s | 68.526% | reject: worsens score |
| prescribed directional wake | 53.964% | reject: worse than 53.283% |
| multilayer + bounded wake recirculation 0.5 | **52.441%** | **retain** |
| exterior density-current strength 0.25 | 73.408% | reject |
| exterior density-current strength 0.50 | 73.281% | reject |
| exterior density-current strength 0.75 | 73.154% | reject |
| exterior density-current strength 1.00 | 73.027% | reject |
| conservative multilayer + sensor envelope | 53.283% | intermediate |
| earlier 40 m hand-off + vertical tail + horizontal diffusion | **44.560%** | **retain as diagnostic candidate** |

The evidence indicates that the dominant correctable error was the one-layer
vertical observation/transport closure combined with the too-small local
domain and analytic exterior route. Source scaling and speed-only density
current continuation are sensitivities, not defensible production corrections;
the simple wake deficit is insufficient to represent separated 3-D flow.

## Scope boundary

The accepted result remains experimental and is not a field-qualified CFD
replacement. A resolved vertical momentum equation, separated wake/re-entry,
and a recoupled outer density-current grid remain future physics work. They are
not silently treated as validated by this audit.

## Verification

The implementation and diagnostics pass the relevant focused suite:

```text
70 passed
```

The code also passes Python compilation and `git diff --check`.

## Same-day strict recheck

The workbook-wall-6 configuration was rerun after the physics sensitivity
work, with the same 600 s horizon and 84 positive receptors:

`outputs/integral-local-patch-2026-10-15/recheck-workbook-wall6.json`

It gives dose MAPE **40.884245%**, prediction runtime 20.95 s, maximum species
residual `8.19e-12 kg`, and maximum momentum-budget residual
`6.98e-10 kg m/s` (gate true). A mass-conserving ground-concentration plus
diffusive-tail observation candidate was also tested at exponents 0--2; its
MAPE was 82.20--79.30% and was discarded. Thus no tested physically declared
candidate reaches the 10% criterion. The remaining gap is not a numerical
conservation failure; it is unresolved 3-D release/head, wall-impact and
lee-side residence physics, compounded by the lack of synchronized raw sensor
time series.
