# TI021 named-interface repair and paired 600 s verification

The original goal remains active and unachieved: all 84 positive-observation
TI021 sensors must reach cumulative dose MAPE <=10% with the required
conservation and group checks. Keep the 29.5441456118% layer baseline.

## Repairs

- The donor now exports an optional named marked-branch packet with mixture,
  species and x/y/z advective momentum. Momentum uses the same branch mass
  and the solved donor velocities. Background anomalies remain distinct.
- HDF5 v2 histories retain compatibility through a separately versioned
  `degali.marked-layer-donor-packet.v1` extension. The extension survives
  checkpoint resume/truncation. Old files lacking normal momentum cannot be
  silently fed to the layer.
- The adapter returns `[mixture, species, px, py]`; pz is audited separately.
  It integrates the piecewise-linear time history exactly over each receiving
  step and remaps lateral face overlaps conservatively. Wider target planes
  get zero exterior flux, without renormalization. Incomplete time histories
  raise an error instead of supplying zero for the missing simulation.
- Tests include donor velocity/component identities, exact coarse/fine time
  integration, fractional face overlap, discarded vertical-impulse accounting,
  checkpoint persistence, and a receiving-layer mass/xy-momentum ledger.
  The selected handoff/layer/preflight regression suite passes 46 tests.

## Frozen comparison

The inner binary donor has the same geometry, wind, source, 600 s horizon,
5 m grid and RI parameters as the previous candidate. It was regenerated to
export the missing named packet fields. Its full runtime was 122.6082 s.
The layer configuration is copied from the existing baseline, with only the
interior flux-replacement interface at x=130 m added. All 145 sensor rows are
retained and the same 84 positive observations are scored.

| Candidate | All-84 MAPE | Normalized bias | Receiving calculation incl. output |
|---|---:|---:|---:|
| Existing layer baseline | 29.54415% | -0.26513 | historical 17.74 s |
| Named packet | 54.64259% | -0.38451 | 13.28 s |
| Same-obstacle control subtracted | 50.81658% | -0.22117 | 14.97 s |

These receiving times exclude the new donor and air-only control solves.
They do not satisfy the 30 s end-to-end requirement for the coupled model.

The named-packet receiver has species residual 1.23e-11 kg and xy momentum
residual 1.40e-9 kg m/s. The control-subtracted receiver has species residual
9.55e-12 kg and xy momentum residual 5.82e-10 kg m/s. This establishes the
numerical receiver ledger, not full physical conservation across a vertically
collapsed interface.

## Additional cause found and tested

Subtracting the nominal inlet carrier from a resolved obstacle plane includes
obstacle-induced background velocity differences in the mixture anomaly.
The positive marked branch therefore kept injecting ambient mixture long
after most species had passed. The complete 0--600 s plane contained
1,432,385 kg marked mixture for only 3295.33 kg species.

An air-only control with the identical holder momentum policy, geometry,
wind and time grid was executed. Subtracting this plane from the donor before
forming the positive marked branch reduced its mixture integral to
26,590.54 kg. The receiving interval starts at 2.575 s, so its delivered mass
is slightly lower (26,582.20 kg). All species packets remain unchanged.

The control simulation completed and flushed all 4801 frames through 600 s,
but its first report attempt used an incorrect result attribute after saving
the histories. The report was recovered by reopening and verifying the three
complete HDF5 histories, without restarting the simulation. The recovery audit
explicitly records that the in-memory control solver diagnostics were not
persisted. Future executions use the correct `domain_mass_kg` attribute and
save control diagnostics before postprocessing. No promotion relies on the
recovered control result.

## Decision and next requirements

Both candidates are rejected: neither improves the layer baseline or its
group requirements. The interface-cycle audit independently reintegrates all
saved concentration traces and verifies the recorded dose, common horizon,
sensor order and observations.

The next coupling cannot use a positive-only anomaly as if it were a full
mixture flux. The control pair contains 23,220.59 kg of negative mixture
anomaly that this closure discards. The control-subtracted receiver also drops
3900.83 kg m/s absolute vertical impulse, and the native coarse upstream layer
still supplies the near-field sensor predictions. These must be resolved in
a signed full-mixture/ambient-state interface with retained vertical state.
Lateral and roof donor exports remain uncoupled. Improving these diagnostics
does not establish the 10% goal or justify calling this result full 3D validation.

## Evidence and reproduction

- `outputs/ti021-interface-cycle-audit-2026-10-17.json`
- `outputs/ti021-named-packet-coupled600s-2026-10-17.json` and `.npz`
- `outputs/ti021-control-subtracted-coupled600s-2026-10-17.json` and `.npz`
- `outputs/ti021-same-obstacle-control-2026-10-17/audit.json`
- `outputs/ti021-binary-local600s-ri5-packets-2026-10-17/`

```powershell
& .\.venv\Scripts\python.exe tools\run_ti021_named_packet_candidate.py --donor outputs\ti021-binary-local600s-ri5-packets-2026-10-17\ti021-inner-handoff.h5 --interface-x-m 130 --output outputs\ti021-named-packet-coupled600s-2026-10-17.json
& .\.venv\Scripts\python.exe tools\run_ti021_named_packet_candidate.py --donor outputs\ti021-same-obstacle-control-2026-10-17\control-subtracted.h5 --interface-x-m 130 --output outputs\ti021-control-subtracted-coupled600s-2026-10-17.json
& .\.venv\Scripts\python.exe tools\audit_ti021_interface_cycle.py
```
