# PRESLHY independent-energy interface result

Date: 2026-09-05

## Outcome

Separating centre density and centre H2 mass fraction supplies the missing
thermodynamic degree of freedom. The pre-registered five-flux projection passes
all seven fixed PRESLHY 10D interfaces without changing a physical coefficient,
handoff location or acceptance threshold. The former single-scalar failures,
trials 10 and 25, both pass.

This was a successful structural diagnosis. The subsequent fifth energy ODE
has now also been implemented and field-tested; its result is recorded below.

## Frozen-gate results

| trial | five-flux max residual | energy quadrature change | points | H2-width mismatch | temperature mismatch, K | former JETPLU interface |
|---:|---:|---:|---:|---:|---:|---|
| 10 | 3.51e-16 | 6.71e-6 | 64 | 4.637% | 1.668 | fail |
| 11 | 3.87e-16 | 9.98e-7 | 64 | 3.073% | 1.832 | pass |
| 12 | 4.45e-16 | 5.90e-7 | 64 | 3.374% | 1.990 | pass |
| 22 | 6.44e-16 | 4.02e-6 | 64 | 4.267% | 1.503 | pass |
| 23 | 2.82e-16 | 4.62e-6 | 64 | 1.728% | 1.142 | pass |
| 24 | 4.01e-16 | 4.78e-6 | 64 | 3.014% | 1.789 | pass |
| 25 | 3.15e-16 | 8.37e-6 | 128 | 4.761% | 1.721 | fail |

The unchanged gates are `1e-8` for every flux, `1e-5` for energy-quadrature
change, 5% for H2 half-width and 2 K for centre temperature. Aggregate maxima
are `6.45e-16`, `8.37e-6`, 4.761% and 1.990 K, respectively.

## Interpretation

The earlier 5/7 outcome was not evidence that the four-flux near field or its
10D station was physically inadmissible. It was a rank mismatch: four JETPLU
physical states were already consumed by mass, species and vector momentum,
leaving no state with which to satisfy total energy independently. Adding the
HyRAM-style density/composition separation resolves precisely those two
failures while preserving the independently checked width and temperature
profiles.

## Downstream field result

The seven-state model integrates the five conserved fluxes directly and
recovers a positive physical state at every Runge--Kutta stage. All seven
trials completed. At the primary 0.02 m step the worst direct flux-minus-source
residual is `5.18e-8`; at the 0.01 m repeat it is `1.04e-7`. The latter is
slightly larger only because twice as many nonlinear state inversions are
accumulated, and remains negligible.

On the same 42 arcs and 17 vertical fits:

| model | MG | VG | FAC2 | mean sigma ratio | centre MAE, m |
|---|---:|---:|---:|---:|---:|
| corrected JETPLU baseline | 1.074 | 1.213 | 0.952 | 1.091 | 0.054 |
| independent energy, 0.02 m | 1.118 | **1.176** | **0.976** | **1.066** | 0.078 |
| independent energy, 0.01 m | 1.118 | **1.176** | **0.976** | **1.066** | 0.078 |

The step-refined width ratio changes by only `6.9e-5` and centre MAE by
`2.7e-6 m`; all reported concentration metrics are unchanged at three decimal
places. The result is numerically converged.

The model is **not promoted** under the frozen all-metric rule. It lowers
variance, increases FAC2 and moves width closer to one, but `abs(log(MG))`
worsens and centre-height MAE rises from 0.054 to 0.078 m. The independent
energy state is therefore retained as the correct structural research base,
not as the production default. The next physics audit must isolate the extra
buoyant rise without retuning the already supported local-shear width closure.

Machine-readable results are in
`reference/preslhy/independent_energy_interfaces_2026-09-05.json` and
`reference/preslhy/independent_energy_field_2026-09-05.json`.

## Follow-up audits

Applying the same single-Gaussian fit to direct-plus-ground-image model values
at the actual sensor heights does not explain the rejection. All 17 fits pass
the frozen R-squared gate, but candidate sensor-fit width ratio/centre MAE are
1.420/0.137 m against baseline 1.360/0.070 m. See
`prereg-ground-image-geometry-observation.md`.

Li et al. (2026) equation 35 was then tested as a coefficient-free published
enthalpy-only alternative to the HyRAM+ total-energy equation. It passes all
interfaces and balances, and changes the internal centre MAE only from
0.0779385 to 0.0779011 m. This direction is retained for research continuity,
but the negligible change rules out resolved kinetic-energy thermalisation as
the cause of the rise defect. See
`prereg-established-flow-energy-partition.md`.

Finally, direct radial integration of `g (rho_amb-rho)` was compared with the
exact finite-Gaussian force after projection. All seven signs agree and the
largest magnitude difference is 6.761%, below the frozen 10% gate. Thus the
five-flux handoff does not inject a material vertical-force discontinuity;
interface density shape is also rejected as the rise-error mechanism. See
`prereg-interface-buoyancy-moment.md`.

The Houf entrainment-width conversion was then corrected against the actual
HyRAM velocity profile. The required `/lambda` factor reduces the uncapped
buoyancy entrainment by 25.7%. On the repeated seven-trial field run, MG/VG
improve to 1.062/1.163, internal width ratio to 0.976 and internal centre MAE
to 0.0758 m. It remains research-only because centre MAE is still worse than
the 0.0542 m baseline. This physically corrected run supersedes the earlier
scalar-width field numbers; see `prereg-houf-width-mapping-correction.md`.

A trial-10 conservative budget then locates the remaining centre-height
defect. The candidate and baseline differ by only 0.0657 N vertical momentum
at 10D, but accumulate 5.695 versus 2.734 N buoyancy by 6 m. The candidate is
already 64.3 K and weakly buoyant at 10D while the baseline is 30.7 K and
dense. The next target is therefore establishment-to-downstream thermal and
density evolution, not an initial vertical-momentum correction; see
`prereg-trial10-vertical-momentum-budget.md`.

The first frozen separation, complete suppression of dry-air condensation,
does not account for that thermal contrast. Trial 10 remains 65.559 K at the
transferred centre and weakly buoyant; the interface temperature mismatch
increases to 2.222 K and fails closed. Per protocol no seven-trial score was
run. Finite nucleation delay is therefore rejected as the next correction,
and pre-10D mass entrainment is the remaining isolated target; see
`prereg-delayed-dry-air-condensation-bound.md`.

The next source audit identifies a genuine zone-bookkeeping issue. The
PRESLHY source is already the air-loaded Station-3 endpoint, while the
`entrained_mass` Gaussian boundary also integrated the Zone-V `E_mom` over
the Zone-IV 6.2D development length. A coefficient-free `source_flux` lower
bound removes that unsupported second application and conserves the four
Station-3 fluxes. It passes 7/7 and improves MG/VG/FAC2 to
1.012/1.162/0.976; centre MAE falls from 0.0758 to 0.0707 m. It is still not
promoted because baseline centre MAE is 0.0542 m and its width error is
slightly larger. See `prereg-source-flux-gaussian-establishment.md`.

Transferring JETPLU's ground-contact geometry and the existing DEGADIS
surface-layer entrainment to this research path does not close the residual.
In trial 10, contact starts at 2.60 m and the strongest bound lowers the 6 m
centre only from 0.737 to 0.695 m, still 0.469 m above the fitted measurement.
The 1.78 m section is unchanged and the frozen all-section gate stops the
seven-trial run. See `prereg-independent-energy-ground-contact.md`.

Applying Li equation 35 from Station 4 rather than only after 10D is also
rejected before field scoring. With a `source_flux` boundary the representative
plug-to-Gaussian solve cannot conserve mass, species, momentum and enthalpy
simultaneously; residuals remain 27--58% on three balances. This exposes the
need for an additional turbulent-energy/profile state across the development
zone rather than a total-energy/enthalpy switch. See
`prereg-zone-v-enthalpy-nearfield.md`.
