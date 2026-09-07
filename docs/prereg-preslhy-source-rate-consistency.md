# Pre-registration: PRESLHY source-rate consistency

Date: 2026-09-05

## Question

The reproducible reduction contains both a sustained-window mean and its
peak, while the independent-energy entry points previously hard-coded only
the mean.  The legacy `conditions.csv` also contains `flow_used_gs`, but that
column mixes report set-points and a different reduction and is explicitly
non-authoritative in `DATA_AND_REPRODUCTION.md`.  It is therefore not used to
select the model input.  Instead, this audit treats the measured peak as a
deliberately high source-rate uncertainty bound.  In trials 10 and 23 the
bound raises the hydrogen source from 0.189 to 0.285 kg/s and from 0.121 to
0.186 kg/s, respectively.  The question is whether this bound alone can
explain the early dilution, heating and excessive buoyancy.

## Frozen change

- Add an explicit `source` selector to the independent-energy reduced-data
  entry points, matching the selector already exposed by the legacy field
  path.  Keep `flow_mean_gs` as the backward-compatible default until the
  comparison is complete.
- Compare `flow_mean_gs` against `flow_peak_gs` from the same sustained raw
  mass-flow window.  Do not use or reinterpret `flow_used_gs`.
- Keep every other setting fixed: source-flux Gaussian establishment, Li
  enthalpy transport in both regions, equilibrium dry-air phase closure,
  Houf velocity-width mapping, local-shear crosswind entrainment, free ground
  geometry and 0.02 m maximum step.
- Run all seven pre-selected trials 10, 11, 12, 22, 23, 24 and 25.  Do not
  select trials by their agreement.
- Use the same concentration arcs and vertical fits available in the bundled
  reduction for this mechanism isolation.  The separately discovered Table
  A3 parser omissions will be corrected and rescored as a distinct data audit
  so their effect cannot be credited to the source-rate change.
- On the raw-temperature subset, use the already frozen 42 post-handoff
  thermocouple points and their minimum-temperature statistic from
  `prereg-preslhy-temperature-field.md`.

## Frozen decision rule

The selected-rate branch is a better source specification only if all of the
following hold without tuning:

1. all seven five-flux handoffs pass the existing conservation, quadrature,
   H2-width and centre-temperature gates;
2. concentration moves closer to unity in absolute log-MG, with VG no larger
   and FAC2 no smaller;
3. centre-height MAE does not increase and vertical-width ratio moves no
   farther from unity;
4. the 42-point temperature median absolute error falls, and its absolute
   median bias falls; and
5. the maximum downstream balance residual remains below `1e-6`.

If only the temperature path improves, retain the rate as a source-uncertainty
bound rather than promoting it.  D3.6 Table 4 is reported as an external
cross-check but is not substituted after results are seen.

## Results

All seven peak-rate runs passed the five-flux interface screens and the
largest downstream balance residual was `4.23e-8`.  Relative to the frozen
mean-rate `source_flux` result, the peak bound improved concentration MG from
1.0117 to 1.0003, vertical-width ratio from 0.8986 to 0.9461, centre-height
MAE from 0.0707 to 0.0605 m, and the 42-point temperature median absolute
error from 20.31 to 16.29 K.  However, VG worsened from 1.1616 to 1.1790;
therefore decision rule 2 fails and the peak rate is not promoted.

The thermal mechanism remains unresolved rather than repaired.  Trial 23's
peak-rate `source_flux` median bias is still +35.02 K, with centreline errors
of +76.72 K at 1.19 m and +86.16 K at 1.78 m.  The peak rate is retained only
as an input uncertainty bound.  The full numerical record is
`reference/preslhy/temperature_source_rate_audit_2026-09-05.json`.
