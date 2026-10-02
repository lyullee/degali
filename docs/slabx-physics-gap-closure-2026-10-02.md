# SLABx / SLABx-LH2 physics-gap closure (2026-10-02)

## Decision

The comparison target is DEGALI's opt-in LH2 research path, not the verified
DEGADIS-compatible core.  Existing baseline results are therefore unchanged.

| Physics | DEGALI state after this change | Boundary |
|---|---|---|
| pressure source / flash / jet | already implemented by the notional nozzle, post-flash droplet source, conserved Gaussian jet and independent-energy crosswind | no automatic selection among incompatible source descriptions |
| source input handoff | post-flash liquid rate, breakup diameter, density and vector jet velocity now map directly to droplet transport | alternative diameter populations must be declared |
| droplet crosswind and gravity | class-resolved 3-D mean-wind transport with Schiller--Naumann drag and a declared d-squared evaporation coefficient | no stochastic turbulent dispersion or secondary breakup |
| rainout and pool formation | impact rate/centroid/RMS radius plus an exact H2 component ledger; deposited mass feeds a declared disk/annulus and is advanced by a conservative axisymmetric shallow-layer spreading solver with local-contact substrate evaporation | the deposition footprint remains an input; slope, curbs, drains, obstacles, wind shear on the liquid and non-axisymmetric flow are unresolved |
| steady crosswind / yaw | already implemented by the independent-energy and yawed crosswind solvers | reverse axial branch and resolved gust dynamics remain outside scope |
| time-varying wind | a circular-direction/speed applicability gate is automatically applied to a measured pre-transition window; measured wind vectors drive puff momentum after transition with exact sample-boundary splitting | the pre-transition plume still uses one nominal vector and is rejected/marked conditional when the measured window is not steady |
| finite plume-to-puff switch | the material clock is integrated from conserved fluxes; at source duration a full-section mass/momentum/energy puff handoff is produced and continued by a native 3-D puff | field validation of the puff concentration branch is still pending |
| puff continuation | H2 inventory conservation, ambient entrainment, steady or measured-vector wind momentum relaxation, buoyancy, form drag, ground contact, anisotropic horizontal gravity spreading, thermodynamic reconstruction and fixed-receptor histories | no obstacles, terrain, combustion or stochastic meander |

## Conservation rules

Droplet transport separates each liquid class into airborne vapour, airborne
liquid at the trajectory time limit, and ground liquid.  These three streams
must sum to the injected liquid rate.  Post-release pool coupling adds direct
flash vapour, pool vapour, remaining pool inventory and any radial-domain
escape to the same hydrogen ledger. Positivity correction is reported as a
separate numerical adjustment in every dynamic-pool ledger row.

Finite-release transition uses the H2 species clock.  The declared H2 source
rate must equal the transported H2 flux within a declared tolerance.  Total
carrier-mixture flux is used only for total puff inventory and bulk velocity;
using it as the H2 clock is rejected.  DEGALI integrates the full Gaussian
cross-section, so SLAB's quarter-cloud factors are not copied.

## Intentional non-claims

This closes the previously silent interfaces and adds a native conservative
puff rather than extending the steady plume beyond source cessation.  The
puff is a DEGALI Gaussian full-section formulation, not a line-for-line SLAB
quarter-cloud port.  It supports downstream fixed-receptor histories, but its
concentration branch remains a research result until independent finite-LH2
field validation is completed.

The default pool path supports concurrent constant rainout, radial spreading
and evaporation. It uses a declared deposition disk/annulus because a small
deterministic droplet sample cannot identify a turbulent impact footprint.
The existing fixed-area calculation remains available only through the
explicit ``pool_model='fixed'`` option. A moving, non-axisymmetric deposition
footprint still requires a different model.

## Code and verification

- `degali.addons.droplet_rainout`: source adapter, class transport, impact
  footprint statistics and post-release pool ledger.
- `degali.addons.dynamic_pool`: conservative radial shallow-layer spreading,
  local substrate-contact evaporation and domain/numerical mass ledgers.
- `degali.addons.finite_release`: finite-source puff handoff and steady-wind
  applicability gate.
- `degali.addons.finite_puff`: native conservative 3-D puff and receptor
  observation operator.
- `degali.lh2.run_lh2_finite_release_research`: one-call jet, yawed plume,
  exact source-clock transition and puff continuation.
- `degali.lh2.run_lh2_rainout_pool_research`: one-call post-flash droplet,
  crosswind/gravity transport, rainout and concurrent dynamic spreading pool.
- `tests/test_dynamic_pool.py`, `tests/test_droplet_rainout.py` and
  `tests/test_finite_release.py`: mass, source-vector, pool, H2-clock,
  transition-distance, domain escape and circular wind tests.

The complete repository regression after these additions was **1158 passed,
145 skipped, 0 failed** (2026-10-02); after aligning the validated evaporation-
momentum default and adding the constant-flux water boundary, the focused
jet/rainout/pool/puff regression was **36 passed, 0 failed**. A representative
5 s flash/rainout run
closed the dynamic-pool mass ledger to 7.11e-15 kg and the end-to-end hydrogen
ledger to -1.04e-14 kg. These are implementation/conservation checks, not an
independent field validation of spreading radius or puff concentration.

For a like-for-like 0--5 s component case using the same CoolProp properties,
inflow, radial grid and validated zero-radial-momentum-vapour closure, the
DEGALI solver and the DynamicLH2PoolX route used by SLABx-LH2 agreed at every
reported time: maximum liquid-inventory difference 0 kg, maximum cumulative-
evaporation difference 0 kg and final reported-radius difference 0 m. This is
cross-implementation parity, not field validation.

The DEGALI implementation was also run directly against the local open JUEL-
3155 figure transcription and the frozen DynamicLH2PoolX Stage-C protocol.
It reproduced the holdout radius metrics: water Trial 4 RMSE 0.13744 m and
aluminium Trial 6 RMSE 0.10159 m, with 60 s radii 0.40 and 0.46 m. Maximum
mass-ledger residuals were 4.81e-13 and 3.49e-13 kg. This supports only the
declared-ground-inflow, horizontal smooth-axisymmetric component scope; it
does not validate rainout footprint prediction or non-axisymmetric spills.
