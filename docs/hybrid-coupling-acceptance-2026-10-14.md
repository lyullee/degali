# DEGALI hybrid coupling acceptance

## Implemented

- Cartesian solver-face coordinate generation.
- Strict patch-to-face alignment checks; no silent interpolation.
- `face_flux_callback` in the fast 3-D transport step and solver.
- Patch boundary face integration and `PatchFluxAccumulator` history.
- Conservative overlap mass remap and shared-clock exchange executor.
- Face-resolved packet extraction for Cartesian x/y/z interfaces, including
  positive outward total/species mass, sensible enthalpy and momentum fluxes;
  incoming packets can be distributed to a boundary cell strip.
- Bidirectional state-callback contract with donor-mass limiting.
- Bidirectional conservative field packet contract for total mass, H2 species
  mass, sensible enthalpy, and momentum; packets above donor inventory are
  rejected rather than clipped.
- Solver-result time alignment and total-species conservation audits.

## Quantitative gate

The obstacle result remains `conditional` until one coupled run demonstrates,
on the same time grid, that the coarse and patch states are advanced with
equal-and-opposite interface face fluxes and that the total hydrogen mass
residual is below the declared tolerance. Independent runs, geometric weights,
or a residual inferred after the fact do not satisfy this gate.

The 4 vol% criterion must use the declared H2/air pressure-temperature closure;
an ambient, zero-enthalpy conversion is diagnostic only and is not a
cryogenic-facility validation.

## Implementation checkpoint (2026-10-14)

- Fast obstacle transport now exposes scalar face-flux, interface-state, and
  conservative x-min species-inlet-flux callbacks in the timestep loop.  The
  latter is included in the solver mass ledger, so a handoff is not a field
  overwrite.
- The shared bridge also exposes explicit state-packet mass/enthalpy ledger
  callbacks.  State application and ledger accounting are separate hooks, so
  a nonzero packet cannot silently appear as an unaccounted concentration
  overwrite.
- Resolved obstacle transport exposes callbacks for the hydrogen scalar and sensible-enthalpy loops, before each local transport step.
- Both 3-D solvers now expose an optional six-array momentum interface callback;
  the resolved solver additionally exposes a conserved-state callback that is
  invoked before the thermal corrector, so incoming patch packets are not
  skipped by the variable-density path.
- Transactional bidirectional smoke test: 4 steps, initial/final field mass both `0.125 kg`, with no unaccounted transfer.
- Conservative field exchange smoke: shared coarse/patch clock preserves all
  four ledgers to `1e-10` absolute tolerance across four patch steps.
- Regression subset: 98 tests passed; syntax compilation and whitespace checks passed.
- This is an interface/ledger checkpoint only; no facility validation accuracy claim is made.
- The TI021 local-patch runner records this as
  `coupling_mode=one_way_handoff_frame_export` until the real outer and patch
  callbacks are wired to this packet contract.
- The callback hooks are now exercised by a shared-clock driver, but no
  facility-quality run is promoted until face-resolved total-mass, species,
  enthalpy and momentum packets use machine-registered geometry and the
  resulting field ledger closes under the external validation gates.
- A concurrent two-solver smoke now exercises that driver path for two shared
  steps: fast outer + resolved patch, two exchanged packets, zero reported
  solver mass residual, wall time about `0.024 s`.  It remains conditional
  because the smoke uses overlapping toy domains, nominal outer density and
  normal momentum flux only.
- The TI021 fast outer screening runner now reads the durable HDF5 handoff
  lazily and injects its remapped species flux into a 5 m obstacle grid.  The
  600 s screening replay completed in `15.47 s` with maximum mass residual
  `2.46×10⁻11 kg`; it is still conditional because the flow closure is
  potential-flow and the current handoff species is R12 rather than a
  calibrated H2 event.

## FFI O1 production-shaped witness (2026-10-14)

`tools/run_ffi_o1_hybrid_patch_prefix.py` now wires the declared callback
contract to three non-overlapping O1 blocks (fast/resolved/fast), rather than
the overlapping two-step toy smoke.  The 60 s run completed with 3,000 shared
clock exchanges on each bridge, runtime `56.26 s`, and maximum partitioned mass
residual `1.56×10⁻17 kg`.  The resolved patch uses active thermal buoyancy and
the new characteristic momentum option: an incoming packet is mapped through
`u_in=p_in/m_in` to the boundary cells and normal face, while species mass and
sensible enthalpy remain explicit conservative packet ledgers.  Its maximum
patch Mach was `0.125`; the prior
additive cell-impulse mode exceeded Mach `0.30` over the same window.
The resolved patch reported zero energy residual, continuity residual
`1.97×10⁻15 kg m⁻³ s⁻¹`, and maximum pressure correction `0.147 Pa`.

This is the first executable outer/patch/downstream witness, not promotion:
FFI beam CAD, source/nozzle closure, clock registration and sensor calibration
remain conditional, and the patch has ambient-enthalpy/no-local-source input in
this run.

The runner also has an explicit isothermal time-convergence mode.  A matched
3 s `dt=0.02/0.01 s` pair changed the common sensor means by only
`2.40×10⁻4 vol%-points` mean absolute (`0.00377` maximum), with field RMSE
`0.00357 vol%` and mass residuals below `2.2×10⁻19 kg`.  This numerical audit
does not replace the active thermal run because the latter is tied to the
D3Q19 physical timestep.
