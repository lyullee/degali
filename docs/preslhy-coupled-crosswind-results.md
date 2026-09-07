# PRESLHY coupled near-field/crosswind results

Date: 2026-09-05

## Outcome

The new conserved thermodynamic near field can be coupled to JETPLU and gives
good concentration and centre-height performance on independent PRESLHY data.
The field test also found and repaired a specific downstream width defect:
constant source-momentum entrainment must end at the near-field boundary and
JETPLU's local density-scaled shear entrainment must resume.

The complete coupled model is not yet promoted. The scalar-peak Gaussian
boundary does not enforce the total-mass target implied by its own beta-based
formation-zone air influx. The four-flux boundary removes that inconsistency,
but two of seven large-source cases cannot be represented by JETPLU's single
thermodynamic scalar while simultaneously meeting the unchanged temperature
and H2-width screens.

## Frozen population

The condition-only filters select PRESLHY trials 10, 11, 12, 22, 23, 24 and
25. All are unobstructed horizontal releases with both orifice and atmospheric
source velocity/wind ratios at least 10. Model predictions use measured-window
mean H2 flow, the same reduced JSON, the same sensors and only arcs downstream
of the handoff. No coefficient was fitted to these observations.

## Results

Because rejected interfaces change the common population, each row states its
own sample size. Baselines are recomputed on exactly the same arc/fitting keys.

| candidate | interfaces | concentration n | MG | VG | FAC2 | vertical n | mean sigma ratio | centre MAE, m | decision |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| scalar boundary, source-momentum downstream | 7/7 | 42 | 0.925 | 1.183 | 0.952 | 17 | 0.710 | 0.049 | reject: width |
| scalar boundary, local shear downstream | 7/7 | 42 | 1.094 | 1.198 | 0.952 | 17 | **0.994** | 0.049 | entrainment mechanism supported; boundary not four-flux |
| four-flux boundary, source-momentum downstream | 5/7 | 33 | 0.970 | 1.197 | 0.970 | 12 | 0.713 | 0.040 | reject: 2 interfaces |
| four-flux boundary, local shear, 10D--20D search | 5/7 | 33 | 1.165 | 1.231 | 0.970 | 12 | **1.019** | 0.040 | reject: 2 interfaces |

For the seven-trial/42-arc population, the corrected atmospheric baseline is
MG 1.074, VG 1.213, FAC2 0.952, mean sigma ratio 1.091 and centre-height MAE
0.054 m. Thus local shear repairs width and slightly improves VG and centre
height, but its MG is marginally farther from one than the baseline; the
original all-metric promotion rule remains failed.

For the five-trial/33-arc four-flux subset, the corrected baseline is MG
1.146, VG 1.240, FAC2 0.970, sigma ratio 1.149 and centre MAE 0.040 m. The
combined candidate gives good geometry and slightly lower variance, but a
subset result cannot override the pre-registered 7/7 interface gate.

## Numerical and conservation audit

The original fixed 32/64-point energy quadrature falsely rejected five scalar
handoffs. Adaptive doubling to at most 512 points retains the same `1e-5`
criterion; every scalar case converged with 64 or 128 points. Native handoff
mass, H2 and vector-momentum residuals are around machine precision. Physical
energy mismatch remains below 0.46%, H2-width mismatch below 4.49% and centre
temperature mismatch below 1.25 K for the accepted scalar interfaces.

The four-flux source boundary closes total mass, H2 mass, momentum and total
energy below `1.4e-15` in all seven trials. At fixed 10D, trials 10 and 25
miss the 2 K centre-temperature screen by only 0.069 and 0.106 K; trial 25
also misses the 5% width screen by 0.036 percentage point. Moving the handoff
out to 20D worsens rather than repairs both mismatches, reaching roughly
6.5 K and 8.3%.

The scalar boundary's 20.58% representative discrepancy is not H2 loss. Its
Gaussian total mass is below the source-plus-air target calculated from
`beta_A=0.28`; H2, momentum and energy close. It therefore implies a different
formation-zone air influx from the one subsequently used by the ODE, which is
why it is not called a fully four-flux boundary.

## Model decision

- Keep the legacy atmospheric `assess()` result as the validated default.
- Keep `crosswind_entrainment="local_shear"` as the supported coupled research
  option; do not continue constant source-momentum entrainment downstream.
- Keep `nearfield_establishment="entrained_mass"` as the physically complete
  source-boundary option, with fail-closed interface behavior.
- Do not tune spreading ratio, entrainment coefficient, handoff distance or
  the 2 K/5% screens to admit trials 10 and 25.
- The next model extension requires a separately transported thermal/energy
  scalar or profile width in the crosswind state. JETPLU's present single
  concentration table cannot represent the multi-valued phase/temperature
  history exposed by the four-flux boundary.

That extension has since been implemented and prospectively tested. A HyRAM-
style state with independent centre density and H2 fraction represents all
five near-field fluxes and passes all seven unchanged 10D width/temperature
screens, including trials 10 and 25. Its downstream energy ODE improves VG,
FAC2 and width but worsens MG and centre-height MAE, so it is retained as a
research base and not promoted; see
`preslhy-independent-energy-interface-results.md`.

Reproduction is provided by
`degali.validation.nearfield.coupled_from_reduced`; its arguments explicitly
select Gaussian establishment, crosswind entrainment and compatible-handoff
search.
