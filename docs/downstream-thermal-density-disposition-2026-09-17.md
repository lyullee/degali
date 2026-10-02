# Downstream thermal/density disposition — 2026-09-17

## Evidence rerun

The public PRESLHY E3.5 thermal-envelope audit was rerun locally without
trajectory reintegration or coefficient selection.  The owner-held raw
workbooks were read in place and remain excluded from the package. Their
source is Lyons, Coldrick and Atkinson (2023), PRESLHY E3.5, DOI
[10.35097/1481](https://doi.org/10.35097/1481).

For the density-profile control, the 4 m centre temperature deficit exceeds
the recorded mean by 38.33 K in Trial 10 and 13.96 K in Trial 23. Both values
lie near the upper end of their respective finite-record envelopes. The
same-width enthalpy-profile candidate is colder still at that station.

The paired 4 m / 1.78 m deficit ratio is above the observed 95th percentile
in Trial 10 and near its upper tail in Trial 23. The persistent direction after
dividing by the upstream deficit means that a source-amplitude change alone
cannot explain the downstream error.

## Decomposition and disposition

| Candidate explanation | Evidence disposition |
|---|---|
| Loss of post-flash sensible enthalpy | Corrected. The flash source now carries its physical post-flash enthalpy instead of reconstructing a subcooled liquid from saturated quality. |
| Mean resolved kinetic energy turning into heat | Rejected as primary cause: the required enthalpy gaps are orders of magnitude larger than the local mean kinetic-energy scale. |
| Different scalar/thermal width with the current closure | Rejected for promotion: it improves width error but worsens centre amplitude under the frozen joint gate. |
| Humid-air condensation absent from the calculation | Rejected: the active thermodynamic path already includes the stated equilibrium phase terms; the record does not isolate an omitted latent source. |
| Unresolved turbulent energy and heat transport | Plausible but unidentifiable from these records: no matched velocity RMS, Reynolds stress, length scale or dissipation observation is available. |
| Minimum-versus-steady sampling | Material uncertainty, but not a full explanation: the normalised downstream decay retains the same error direction. |

## Resulting model decision

No heat-transfer multiplier, turbulent Prandtl number, buoyancy multiplier, or
source-rate fit is enabled. The physically correct source-enthalpy preservation
is retained and regression-tested. The remaining downstream thermal/density
residual is recorded as a boundary on the present fast model until an
independent velocity/thermal profile can identify a transport closure.
