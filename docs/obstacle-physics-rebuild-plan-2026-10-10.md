# Obstacle-physics rebuild plan

Date frozen: 2026-10-10

## Objective

Rebuild DEGALI as a validated multi-fidelity system for transient hydrogen
dispersion in facilities containing finite obstacles.  The legacy non-local
scalar-routing operator is not an admissible quantitative obstacle model.  A
fast screen and an obstacle-resolved tier must instead share explicit
conservation laws, solid boundary conditions, evidence gates and comparable
outputs.

This is a model-development objective.  A journal submission is neither the
definition of completion nor evidence that the physics is correct.

## Concrete product target and boundary

DEGALI is to become a **dynamic local-facility hydrogen consequence model**,
not a second free-field slab model and not a general-purpose CFD replacement.
Its intended decision window is the region where finite equipment, walls,
roofs and sensor locations materially change a time-varying hydrogen cloud.
SLABx-type free-field/far-field dispersion remains complementary rather than a
competitor.

The declared input contract is one `FieldEvidenceManifest` containing a common
event clock, time-varying release boundary, wind/weather history, finite 3-D
obstacles, receptor/sensor geometry, calibration evidence and immutable file
hashes.  The declared output contract is:

1. transient 3-D hydrogen concentration in vol%;
2. arrival, peak, dose and duration above declared thresholds at sensors;
3. 4 vol% instantaneous and maximum envelopes;
4. windward, roof, side-bypass and lee/wake diagnostics;
5. mass, energy, continuity, grid/time, applicability and runtime evidence;
6. an uncertainty envelope and an automatic `accepted`, `conditional` or
   `withheld` decision.

Ignition, combustion, flame acceleration, blast, two-phase source formation
inside process equipment and approval/design-basis decisions remain outside
the current product target.  Those boundaries must not be hidden by a good
dispersion score.

## Definition of done

The rebuild is complete only when all of the following gates pass on frozen,
reproducible cases.  A manuscript or release tag is not a substitute.

### G1 -- numerical core

- relative hydrogen-mass and sensible-energy residuals are at most `1e-6`;
- impermeable-solid normal mass/species flux is zero to solver tolerance;
- projected local continuity residual is at most `1e-7 kg/m3/s` on the
  declared reference cases;
- no negative or above-unity mole fraction is produced;
- lattice-Mach, density-departure and nonlinear-corrector gates all pass;
- two adjacent grid and time refinements change every primary safety output
  (sensor peak/dose, 4 vol% volume and threshold duration) by at most 10%, and
  full-field differences decrease monotonically.  Non-monotonic refinement is
  `withheld`, even when absolute differences are small.

### G2 -- obstacle-flow component evidence

- analytic channel/Poiseuille, wall behaviour and finite-cube wake tests pass;
- AIJ Case H is evaluated on a pre-frozen common observation set with no
  closure tuning on the evaluation points;
- a candidate SGS model replaces the simpler baseline only if it reduces both
  `U` and `W` RMSE by at least 5%, does not reduce either correlation by more
  than 0.01, preserves recirculation-sign accuracy and passes the same
  grid/time and flow-residual gates;
- inlet/outlet and wall treatment are named and independently tested.

### G3 -- dynamic hydrogen evidence

- at least two independent hydrogen datasets are evaluated: one finite-
  obstacle/open-or-semi-confined case and one enclosure/channel case;
- source, weather, geometry, sensors, calibration and clock evidence are
  complete in `FieldEvidenceManifest`; otherwise the case remains
  `conditional` or `withheld`;
- at sensors whose observed peak is at least 1 vol%, the frozen targets are
  median absolute peak error <= 1 vol%-point, median absolute dose error
  <= 30%, and median normalized arrival-time error <= 20%;
- 4 vol% event detection has sensitivity and specificity of at least 0.90,
  reported together with false negatives rather than hidden in an average;
- scores are reported by region (windward/roof/side/wake) and dataset, not only
  as one pooled number.

These are engineering development targets, not claims about current accuracy.
They must be frozen before evaluation and may be revised only with a versioned
rationale, never after seeing a held-out result.

### G4 -- usable multi-fidelity workflow

- Tier 1 returns a conservative screen and evidence decision in <= 10 s on a
  committed reference case and named reference machine;
- Tier 3 returns the converged reference transient in <= 5 min on the same
  machine, or explicitly reports that the runtime target was missed;
- both tiers consume the same manifest and emit compatible sensor/threshold
  outputs, so escalation is automatic and auditable;
- a single documented command reproduces every acceptance record and figure.

### Promotion states

- **Research-ready:** G1 and G2 pass; external hydrogen results may still be
  conditional.
- **Field-pilot-ready:** G1--G4 pass on the two frozen hydrogen datasets and a
  prospective, non-calibration field event.
- **Design-basis/approval-ready:** explicitly outside this rebuild; it would
  require a separate governance, uncertainty and authority review.

## Ordered execution path

1. Stabilise the resolved flow horizon and qualify inlet/outlet boundaries;
   do not choose an SGS closure until it materially beats the simpler model.
2. Pass the coupled species/energy grid-and-time gate at roof and side
   receptors, where the current solver still changes by 25--35%.
3. Freeze and run the two hydrogen validation packages without calibration on
   their evaluation observations.
4. Connect the validated tiers to `FieldEvidenceManifest`, uncertainty
   propagation and automatic fidelity escalation.
5. Only then freeze an industrial pilot release and prepare any publication
   around the evidence actually earned.

## Finding that triggered the rebuild

The legacy semi-FV operator removes scalar flux from the last windward fluid
cell and deposits it in the nearest fluid cell found by a non-local search.
For a ground-mounted block this is normally the first cell over the roof.  In
the 2026-10-10 diagnostic case the windward maximum was 4.65 vol%, whereas the
global maximum was 37.41 vol% exactly at the first roof cell.  Integral mass
was conserved, but local momentum, pressure, residence time and scalar-flux
continuity were not represented.  Those concentration and hazard-distance
numbers are withdrawn from quantitative use.

## Evidence map

| Evidence | Physical use | Prohibited inference |
|---|---|---|
| AIJ Case H, DOI 10.5281/zenodo.15429872 | Neutral single-building velocity, turbulence and tracer-wake observations | Cryogenic or universal LH2 wake coefficient |
| Tominaga et al. (2008), DOI 10.1016/j.jweia.2008.02.058 | CFD-quality and building-flow benchmark practice | Validation of DEGALI without matched outputs |
| Nielsen (1991), DOI 10.1016/0950-4230(91)80004-E | Continuous dense-gas wall response, including spatially signed concentration changes | Scalar wall attenuation or enhancement factor |
| Davies and Singh (1985), DOI 10.1016/0304-3894(85)85043-3 | Full-scale dense-gas obstruction regimes | Hydrogen-specific buoyancy law |
| Akhter and Kaye (2025), DOI 10.3390/fluids10080194 and 10.3390/fluids10020046 | Upstream dense release, windward impact, wake trapping and Richardson-number regimes | Direct atmospheric-scale calibration without similarity analysis |
| Olvera and Choudhuri (2006), DOI 10.1016/j.ijhydene.2006.02.022 | Hydrogen buoyancy, roof and lee cavities around a cuboid | Experimental validation of DEGALI |
| Henriksen et al. (2024), DOI 10.1016/j.ijhydene.2024.10.038 | Hydrogen sensor histories with and without channel obstruction; LES comparison | Open-atmosphere facility validation |
| Angot, Bruneau and Fabrie (1999), DOI 10.1007/s002110050427 | Navier--Stokes/Brinkman obstacle penalisation | Turbulence closure or validation data |
| Kolomenskiy and Schneider (2012), DOI 10.1016/j.jcp.2012.02.009 | No-slip flow and no-flux scalar volume penalisation | Automatic engineering accuracy on coarse grids |
| NIST Chemistry WebBook SRD 69, molecular hydrogen | Molecular mass and gas-phase thermochemical reference | Constant-cp accuracy below the tabulated Shomate range or two-phase LH2 closure |
| Lu and Pantano (2020), DOI 10.1016/j.jcp.2019.109132 | Discrete variable-density zero-Mach mass constraint and pressure solvability | Validation of DEGALI momentum/species/energy coupling by continuity alone |
| Zou and He (1997), DOI 10.1063/1.869307 | Pressure/velocity LBM flow boundaries and boundary canonical tests | Automatic suitability for the present D3Q19 obstacle outlet |
| Guo, Zheng and Shi (2002), DOI 10.1088/1009-1963/11/4/310 | Non-equilibrium extrapolation velocity/pressure boundary candidate | Promotion without the present-domain comparison and canonical verification |

## Required governing physics

The resolved tier must advance, in a documented low-Mach or incompressible
approximation, mixture mass, momentum, hydrogen species and sensible energy:

1. a pressure--velocity coupling that constrains the resolved volume flux;
2. no penetration and no slip at solid surfaces;
3. zero normal hydrogen flux at impermeable, non-absorbing walls;
4. temperature- and composition-dependent density or a declared Boussinesq
   limit;
5. gravity acting on the density anomaly;
6. molecular plus modelled turbulent momentum, heat and species transport;
7. finite-width three-dimensional over-roof and lateral bypass;
8. separated roof and lee shear layers, recirculation and reattachment at the
   resolved tier;
9. time-varying source and wind histories on one event clock.

Combustion, ignition, flame acceleration and blast are outside this objective.

## Fidelity ladder

### Tier 0 -- geometry gate

Retain exact obstacle-contact and coordinate-frame checks.  Contact invalidates
an unobstructed prediction; it does not create a concentration correction.

### Tier 1 -- fast conservative screen

Solve a geometry-dependent, divergence-free potential-flow field around a
ground-mounted wind-plane obstacle, then advect and diffuse hydrogen locally
with zero normal solid flux.  This tier represents stagnation and continuous
over-obstacle diversion without a non-local scalar jump.  It deliberately
does not claim no-slip separation or a recirculating wake.

### Tier 2 -- resolved two-dimensional audit

Solve viscous Navier--Stokes with pressure correction or a validated
Navier--Stokes/Brinkman formulation.  Couple species and energy transport and
resolve roof separation, the lee cavity and reattachment.  Use this tier to
falsify Tier-1 assumptions and to exercise grid/time convergence cheaply.

### Tier 3 -- resolved three-dimensional facility model

Extend the pressure-projected variable-density solver to finite cuboids and
multiple obstacles.  Resolve both lateral and vertical bypass and retain the
same source, receptor, conservation and evidence contracts as the lower tiers.

## Work packages and gates

### WP0 -- quarantine legacy routing

- Mark semi-FV routing as legacy and non-quantitative.
- Keep it only for regression, conservation and failure-demonstration tests.
- Remove it from any automatic industrial/model-selection path.
- Preserve the 2026-10-10 roof-cell spike as a regression that the replacement
  must not reproduce by non-local transfer.

### WP1 -- Tier-1 potential-flow transport

- Sparse Laplace solve for a wind-plane streamfunction.
- Velocity derived from the streamfunction, with exact zero scalar flux across
  every solid face.
- Conservative local finite-volume species transport.
- Windward, roof and lee receptor histories, mass ledger and runtime.
- Open-field identity, solid impermeability, divergence, positivity and
  refinement tests.

### WP2 -- Tier-2 pressure/Brinkman solver

- Momentum advection, viscosity and pressure projection.
- No-slip solid and ground, declared top/inlet/outlet conditions.
- Thermal and species transport with independent turbulent Prandtl/Schmidt
  numbers.
- Boussinesq first, then low-Mach variable-density audit.
- Canonical channel, lid-driven cavity and bluff-body tests before dispersion.

### WP3 -- Tier-3 finite-width three-dimensional solver

- Axis-aligned and oriented cuboid masks.
- Finite lateral bypass and arbitrary horizontal wind bearing.
- Multiple obstacles and facility-domain boundary conditions.
- Conservative receptor interpolation outside solid cells.

### WP4 -- external validation

Use immutable train/evaluation splits.  Geometry, boundary conditions and
observation operators are frozen before reading evaluation concentrations.
Report vol% absolute error, MAE, RMSE, fractional bias, NMSE, geometric mean,
geometric variance, FAC2, arrival time, peak, dose and threshold duration.
Neutral, dense-gas and hydrogen evidence remain separate; success in one does
not silently validate the others.

## Numerical acceptance criteria

- cumulative hydrogen mass residual <= 1e-6 of injected mass;
- no scalar inventory in declared solid cells above floating-point tolerance;
- no advective or diffusive scalar flux through an impermeable solid face;
- pressure-corrected normalised divergence below a predeclared tolerance;
- non-negative hydrogen concentration without post-hoc loss of material mass;
- adjacent grid/time refinements change predeclared receptor peak and dose by
  no more than 10% before quantitative promotion;
- every reported field includes grid, time step, conservation, applicability,
  execution time and evidence identifiers;
- failure of a physics or evidence gate yields `conditional` or `withheld`,
  never an unqualified hazard distance.

## Completion evidence

Completion requires all fidelity tiers, tests and external comparisons above.
A visually plausible contour, a passing conservation test, or agreement at one
centreline sensor is not completion evidence.

## Implementation status at 2026-10-10

- WP0 is enforced in code: a legacy semi-FV run that contacts an obstacle sets
  `quantitative_prediction_allowed=False`, and the operational decision layer
  returns `withheld` even if a reviewer enables conditional screening.
- WP1 is implemented in `degali.addons.fast_obstacle_transport`, including a
  divergence-free potential-flow field, local conservative species fluxes,
  zero solid-face flux, vol%, receptors, mass ledger and runtimes.
- WP2 has an explicitly conditional prototype in
  `degali.addons.resolved_obstacle_flow_2d`. It uses a low-Mach D2Q9 viscous
  momentum solver with no-slip link bounce-back, an isothermal lattice
  pressure, and a Helmholtz projection that maps the velocity field onto the
  conservative scalar-control-volume faces. It resolves an emergent negative
  lee velocity and exposes pressure, velocity, vorticity, concentration,
  conservation and runtime diagnostics.
- The WP2 prototype is not a completed Tier 2: its effective viscosity is
  prescribed, thermal energy and hydrogen buoyancy feedback are absent,
  canonical channel/cavity and external building-flow validation are still
  required, and `quantitative_prediction_allowed` remains false.
- The matched diagnostic figure and machine-readable record are generated by
  `tools/compare_obstacle_fidelity_tiers.py`. They are a model-development
  comparison, not experimental validation.
- The D2Q9 momentum and no-slip wall implementation now passes an analytic
  Poiseuille verification: relative L2 profile error 0.220%, relative maximum
  error 0.161%, and mass-flux error 0.240%.
- A four-grid obstacle audit (`60x20`, `90x30`, `120x40`, `150x50`) preserves
  one 4.5 s physical flow horizon and uses bilinear fluid-cell receptor
  interpolation. Peak and dose changes on the finest pair are below 10% at
  windward, roof, lee-low and lee-shear receptors; 4 vol% area and maximum
  downwind extent also change by less than 10%.
- Quantitative promotion is nevertheless withheld because roof duration above
  4 vol% changes by 44.721% on the finest adjacent grids. Threshold duration
  is therefore not numerically converged even though peak, dose and hazard
  envelope geometry meet their current gates.
- The official AIJ Case H package is now ingested with immutable SHA-256
  fingerprints and its measured approach-flow profile is imposed directly at
  the Tier-2 inlet (with an explicit no-slip completion from the lowest
  measurement at 0.01 m to the ground). All 220 available centre-plane mean
  `U/W` points are evaluated without calibration. The comparison gives
  `U/U_H` MAE 0.356, RMSE 0.486 and correlation 0.828; `W/U_H` MAE 0.256,
  RMSE 0.434 and correlation 0.721. Streamwise recirculation-sign accuracy is
  79.5%, and one run takes about 29 s on the current machine.
- This AIJ result fails the quantitative gate. The experiment has
  `Re_H=43,000`, whereas the stability-qualified laminar/effective-viscosity
  run has `Re_H=32.3`; the flow residual remains 0.102; and the two-dimensional
  blockage forces all bypass over the roof instead of allowing the finite
  building's lateral route. AIJ tracer concentration is therefore withheld,
  as are all hydrogen-accuracy claims. The next physics priority is Tier 3
  finite-width momentum with a turbulence treatment, not coefficient fitting
  to this failed two-dimensional comparison.
- The reproducible AIJ run and plot are generated by
  `tools/validate_aij_case_h_flow.py`. Its machine record explicitly separates
  accepted conditional mean-flow evidence from prohibited hydrogen and
  concentration inferences.
- A finite-width three-dimensional fast baseline is now implemented in
  `degali.addons.fast_obstacle_transport_3d`. It projects the prescribed inlet
  flow onto a divergence-free Cartesian face field around one or more finite
  cuboids, so roof and both lateral bypass routes coexist. Hydrogen transport
  is local, conservative, and has exactly zero flux across every solid face.
  It remains a Tier-1 geometry model: slip potential flow has no boundary
  layer, separation, recirculation, turbulence, heat or buoyancy feedback.
- On AIJ Case H, adding finite width without fitting reduced `U/U_H` RMSE from
  0.486 to 0.357 (26.46%) and `W/U_H` RMSE from 0.434 to 0.137 (68.52%). The
  3-D projection used 405 CG iterations, achieved core divergence
  `7.43e-9 s-1`, and ran in about 0.81 s, versus 26.49 s for the current 2-D
  viscous run. This directly supports the diagnosis that the 2-D full-span
  blockage is a major magnitude-error source.
- The 3-D potential result is not a promotion: `W` correlation fell from 0.721
  to 0.629 and the model cannot produce the measured viscous wake. The next
  resolved implementation must combine finite width with no-slip momentum and
  an independently verified turbulence treatment. Reproducible comparison:
  `tools/compare_aij_case_h_2d_3d.py`.
- A first finite-width viscous prototype is now implemented in
  `degali.addons.resolved_obstacle_flow_3d`. It advances the isothermal
  low-Mach D3Q19 Navier--Stokes limit, applies link bounce-back on finite
  cuboids, optionally adds a declared Smagorinsky strain viscosity, and uses a
  separate Helmholtz-projected face field for conservative local hydrogen
  transport. It generates both lateral bypass and negative lee velocity.
- The D3Q19 core independently passes periodic 3-D Poiseuille verification:
  relative L2 velocity error 0.712%, relative maximum error 0.525%, mass-flux
  error 0.776%, horizontal nonuniformity near machine precision, and relative
  density variation below `4e-16`.
- On the common 195 AIJ points available to every tier, the 3-D viscous
  prototype gives `U/U_H` RMSE 0.169 and correlation 0.934, and `W/U_H` RMSE
  0.086 and correlation 0.808. Relative to the matched 2-D result, RMSE falls
  67.11% for `U` and 80.15% for `W`; recirculation-sign accuracy rises from
  80.2% to 92.2%. Runtime is about 11.6 s. These are evaluation results, not a
  promotion or a calibrated fit.
- Quantitative status remains `withheld`: the AIJ run is coarse (`0.02 m`
  cells), has `Re_H=32.3` rather than 43,000, ends at residual 0.00227, and
  uses unverified equilibrium-inlet/copy-outlet and Smagorinsky wall treatment.
  Energy, composition-dependent density and buoyancy feedback are still
  absent. The scalar face projection reduces the raw face divergence from
  `41.6 s-1` to `1.09e-9 s-1`, but this numerical conservation result does not
  cure the remaining physical and convergence gates.
- The first three-grid AIJ audit holds the physical flow horizon at 0.2786 s
  and evaluates only the 176 observations available on all three grids. From
  `0.020 m` to `0.0125 m`, normalized `U` RMSE changes 8.85%, but `U` MAE
  changes 11.28%, `W` RMSE changes 12.08%, and the minimum lee velocity
  changes 23.13%. The finest grid also misses the residual tolerance. The
  spatial-convergence gate is therefore `withheld`; the improved finest-grid
  RMSE is not treated as validation.
- The independent flow-horizon audit holds the `0.020 m` grid fixed and runs
  0.1393, 0.2786 and 0.5573 s. Between the last two horizons, every declared
  AIJ error/correlation/sign metric changes by less than 1.84%, while the
  minimum lee velocity changes 4.58%. This clears the 10% sensitivity part of
  the time gate, but the longest run still has residual 0.00160 and therefore
  remains `withheld` rather than being labelled converged.
- Reproducible records and figures are generated by
  `tools/audit_aij_resolved_3d_grid.py` and
  `tools/audit_aij_resolved_3d_time.py`. The next implementation gate is not
  more presentation work: it is a Reynolds-consistent resolved treatment,
  qualified inlet/outlet and near-wall SGS behaviour, followed by coupled
  energy, species, composition-dependent density and buoyancy.

## Active thermal-scalar checkpoint

- `degali.addons.hydrogen_air_thermodynamics` now reconstructs temperature,
  hydrogen mole fraction, ideal-gas mixture density and vertical buoyancy from
  two independently conserved fields: hydrogen mass concentration and signed
  sensible-enthalpy density. The force sign has an independent canonical
  check: density above ambient produces downward acceleration and density
  below ambient produces upward acceleration.
- The ideal-gas neutral-density map gives 282.240 K at 4 vol% H2 for the
  declared 293.15 K, 101325 Pa, 1.204 kg/m3 ambient. Thus 4 vol% is not
  intrinsically buoyant: below that transition temperature it is denser than
  ambient; above it, it is lighter. The map and machine record are generated
  by `tools/audit_hydrogen_air_buoyancy_transition.py`.
- `resolved_obstacle_flow_3d` now has an opt-in coupled transient path. It
  advances D3Q19 momentum with Guo buoyancy, predicts hydrogen mass and
  sensible enthalpy, reconstructs EOS density, projects **face mass flux** to
  the variable-density continuity constraint, and repeats scalar transport to
  a fixed point. The corrected velocity is reconciled with the D3Q19 momentum
  while retaining its non-equilibrium stress. Receptor output includes H2
  vol% and temperature histories.
- In the reproducible 3.0 s finite-block development witness (0.5 m/s wind,
  80 K source declaration, 0.20 m cells), the current observed local runtimes
  are about 1.1 s for the neutral path and 6.3 s for the corrected active path
  after reusing the immutable pressure-projection matrix and factorization
  (18.17 s before those optimizations).
  The active run reached 274.49 K, 5.21% maximum relative density departure
  and remained below lattice Mach 0.093. All 250 corrector steps converged in
  at most five iterations (mean 3.47). The worst local variable-density
  continuity residual was `7.68e-15 kg/m3/s`, the worst global residual was
  `2.93e-15 kg/s`, and hydrogen-mass and sensible-energy absolute residuals
  were below `7e-18 kg` and `2e-11 J`. The neutral/active 4 vol% maximum-field
  volumes are both 0.816 m3; the largest local maximum-field change is 2.340
  vol%-point. These are development differences, not accuracy evidence.
  Reproduction:
  `tools/compare_obstacle_thermal_buoyancy.py`.
- This checkpoint completes the discrete mass-flux continuity link, but does
  **not** complete the target thermodynamic model. Air carrier inventory is
  still represented by a bounded active-scalar closure, heat capacities are
  constant, and phase change/real-gas source enthalpy are not solved. A
  per-cell relative-density gate, a corrector-convergence gate, a lattice-Mach
  gate and an unconditional quantitative prohibition remain active.

## Variable-density mass-flux projection checkpoint

- `degali.addons.variable_density_projection_3d` now enforces the cellwise
  finite-volume constraint
  `(rho_new-rho_old)/dt + div(rho*u) - S_rho = 0`. Inlet and impermeable faces
  retain their prescribed mass flux, while the positive-x pressure outlet is
  corrected so cellwise continuity and the whole-domain mass budget close
  together. The correction potential is also returned as a pressure increment
  `delta_p = phi/dt` for the next momentum-coupling step.
- Four independent manufactured cases were frozen: random provisional flux at
  uniform density, uniform density decrease, a localized mass source around
  an impermeable block, and a nonuniform density front. The largest post-
  projection local residual is `1.03e-13 kg/m3/s`; the largest global balance
  residual is `1.69e-14 kg/s`; all pass predeclared `1e-8` and `1e-9`
  limits. Total local runtime is about 0.017 s. Cached-direct and CG solutions
  are also compared directly in the unit suite.
- Solid-face mass flux remains exactly zero. A disconnected fluid pocket with
  no pressure-outlet path is rejected instead of sending a singular system to
  the iterative solver.
- The reproducible record and figure are generated by
  `tools/verify_variable_density_mass_flux_projection.py`. Acceptance applies
  only to this continuity operator.

## Coupled low-Mach thermal-corrector checkpoint

- `degali.addons.low_mach_thermal_corrector_3d` now couples hydrogen species,
  sensible enthalpy, ideal-gas mixture density and the mass-flux pressure
  projection in one fixed point. The integrated Tier-3 path feeds the returned
  velocity and pressure correction back into D3Q19 momentum.
- Independent manufactured checks cover ambient uniform flow and a local cold
  hydrogen source upstream of an impermeable finite block. Both converge; the
  cold-source case takes three iterations, closes local continuity to
  `1.04e-14 kg/m3/s`, closes the global balance to `8.18e-16 kg/s`, conserves
  species and sensible enthalpy to the printed precision, and has exactly zero
  solid-face penetration. Reproduction:
  `tools/verify_low_mach_thermal_corrector.py`.
- Status is `accepted_numerical_operator_only`. Physical promotion still
  requires grid studies for this newly coupled path, temperature-dependent/
  real-gas property checks, and external hydrogen obstacle concentration data.
  The AIJ velocity comparison alone cannot validate cryogenic density/species/
  energy coupling.
- Tier 3 now defaults to a bounded MUSCL-Hancock reconstruction with the
  van-Leer limiter; first-order upwind remains an explicit comparison option.
  On a smooth 1-D-in-3-D Gaussian translation, RMS error drops from 0.02794 to
  0.00501 (82.1%) while total mass, non-negativity and the upper bound are
  retained. Slopes are forced to zero where their stencil touches a solid, so
  the no-flux obstacle boundary is not crossed.
- The first coupled time-step audit holds the 0.20 m grid, 7.2 s physical flow
  spin-up, 3.0 s transient and source history fixed while halving `dt` from
  0.012 to 0.006 to 0.003 s. All 1,750 corrector steps converge and retain the
  continuity gates. With bounded MUSCL transport, domain maximum H2 changes
  21.141 -> 21.216 -> 21.239 vol%, and the 4 vol% volume changes
  0.816 -> 0.800 -> 0.800 m3. However, successive
  full-field RMS differences do **not** decrease: maximum-H2 differences are
  0.0300 then 0.0335 vol%-point, final-H2 differences are 0.0342 then 0.0356
  vol%-point, and final-temperature differences are 0.0723 then 0.0752 K.
  Formal time convergence is therefore withheld despite small absolute
  sensitivity. Optimized coarse/medium/fine runtimes are 6.16, 12.35 and
  24.43 s. Reproduction: `tools/audit_active_thermal_corrector_time.py`.
- A separate fixed-Mach joint grid/time audit uses exactly grid-aligned
  obstacles and a constant 1.5 refinement ratio (`dx=0.300, 0.200, 0.1333 m`).
  All correctors and conservation gates pass, but the grid gate fails. The
  bounded-MUSCL finest-pair domain maximum and 4 vol% volume change by 5.86%
  and 0.48%, and the windward peak/dose by 5.84%/5.21%; however roof peak/dose
  change by 28.12%/34.53% and side peak/dose by 25.10%/31.48%. Full-field RMS
  differences also increase instead of decreasing. The status therefore
  remains `withheld_grid_refinement_gate_failed`. Relative to the upwind
  baseline, the safety-threshold volume and windward metrics improve, while
  low-concentration roof/side bypass remains unresolved. The optimized finest
  runtime is 131.88 s (352.26 s before pressure-solver optimization).
  Reproduction: `tools/audit_active_thermal_corrector_grid.py`; the upwind
  baseline can be regenerated with `--scheme upwind`.
- The source kernel and obstacle mask are now discretised consistently across
  grids. Cuboids use half-open cell membership, and the separable Gaussian
  release is integrated analytically over each Cartesian control volume rather
  than sampled at its centre. A nested-grid unit check confirms that eight fine
  cells aggregate to the matching coarse-cell source weight to round-off.
  This correction does **not** cure the coupled concentration-grid failure. In
  the repeated MUSCL audit the finest-pair changes are 7.47% for domain maximum,
  1.56% for 4 vol% volume and 6.50%/5.73% for windward peak/dose, but
  roof peak/dose still change by 29.58%/36.42% and side peak/dose by
  26.29%/33.18%; full-field RMS differences remain non-monotonic. The result is
  therefore still `withheld_grid_refinement_gate_failed`. Reproduction artifact:
  `outputs/obstacle-physics-rebuild-2026-10-10/active_thermal_corrector_grid_audit_cell_integrated_source.json`.
- The next closure gate is WALE, following Nicoud--Ducros and the LBM
  channel/wall-mounted-cube verification sequence of Weickert et al. It must be
  implemented as an opt-in candidate and pass near-wall, channel and cube
  checks before replacing Smagorinsky. A turbulent Schmidt/local scalar-flux
  model is a separate external-data gate and will not be fitted to this grid
  audit.
- WALE is now implemented as an opt-in D3Q19 closure and passes its local
  invariant/unit checks, including zero response to uniform flow and pure
  shear, non-zero response to rigid rotation, solid masking and relaxation-
  time limiting.  It has not replaced the default Smagorinsky option.
- A fixed `0.020 m`, 195-point AIJ Case H comparison gives normalized `U/W`
  RMSE of `0.1680/0.0860` without SGS, `0.1685/0.0863` with Smagorinsky
  `Cs=0.12`, and `0.1684/0.0863` with WALE `Cw=0.50`.  Relative to
  Smagorinsky, WALE reduces the four MAE/RMSE measures by only
  `0.02--0.15%`, slightly reduces `W` correlation, leaves sign accuracy
  unchanged, and increases local runtime from `16.9 s` to `32.4 s`.
  No-SGS is marginally best on this coarse/effectively viscous case, and all
  three flow residuals remain unconverged.  Therefore no closure is promoted;
  the record remains `withheld_candidate_comparison`.  Reproduction:
  `tools/compare_aij_case_h_subgrid_models.py`.
- A Guo--Zheng--Shi non-equilibrium-extrapolation velocity inlet/fixed-density
  outlet was implemented as an opt-in boundary candidate with exact interior-
  plane moment tests. It is not promoted: on the same no-SGS AIJ case all four
  MAE/RMSE measures increase by `2.78--5.97%`, the final instantaneous residual
  increases from `0.00229` to `0.04288`, and the wake is substantially weakened.
  Spatial localisation places the largest residual at `x=0.67--0.69 m`, the
  last two cell centres. Reproduction:
  `tools/compare_aij_case_h_flow_boundaries.py` and
  `tools/diagnose_aij_flow_residual_regions.py`.
- Moving the legacy copy outlet downstream by extending the domain from
  `0.70` to `0.90` and `1.10 m`, while preserving `0.020 m` cells and 1.286
  nominal flow-throughs, changes every matched AIJ metric by less than 1% on
  the longest pair and the minimum lee velocity by 2.55%. Thus the evaluated
  observation region is insensitive to downstream fetch at the 5% gate, but
  the instantaneous residual remains non-monotonic and unconverged. The
  overall audit remains `withheld`. Reproduction:
  `tools/audit_aij_case_h_outlet_fetch.py`.
- The resolved solver now supports an explicit spin-up plus finite-window
  time-mean observation operator. It compares the first and second halves of
  the averaging window rather than requiring an unsteady instantaneous field
  to become stationary. With no SGS, 3000 total steps and the final 1000 steps
  averaged, the mean-field residual is `0.008015`, passing the predeclared 1%
  component-stationarity gate. Relative to the 1500-step instantaneous field,
  `W` MAE/RMSE improve by about 1.26% while `U` RMSE worsens by 0.38%; this
  mixed response indicates observation-operator alignment rather than a tuned
  accuracy claim. The mean `U/W` RMSE is `0.16861/0.08493`, correlations are
  `0.93384/0.80969`, and runtime is 35.5 s. Status is
  `accepted_component_stationarity`; grid and hydrogen concentration gates
  remain open. Reproduction: `tools/audit_aij_case_h_time_averaged_flow.py`.
- The first time-averaged three-grid rerun exposed a geometry-discretisation
  defect rather than a turbulence-closure defect. The historical middle
  `0.020 m` grid placed the declared `x=0.15/0.25 m` block faces through cell
  centres. Inclusive endpoint masking therefore voxelized the nominal
  `0.10 m` streamwise block as `0.12 m`, while the coarse and fine blocks were
  exactly `0.10 m`. On that inconsistent family, `W` model-field change was
  non-monotonic (`15.64% -> 16.21%`) and the grid gate remained withheld.
- Cuboid voxelisation now uses half-open cell-centre membership and reports,
  for every obstacle, declared/effective bounds, signed face and size errors,
  cell counts, maximum error in metres/cells and an exact-alignment flag. This
  prevents a face lying on a cell centre from claiming both endpoint layers;
  it does not pretend that a staircase grid is a cut-cell method.
- A new all-face-aligned grid family (`dx=0.025, 0.0166667, 0.0125 m`) holds
  the total flow horizon at `0.5573 s` and the averaging duration at `0.1858 s`.
  All three mean fields pass the 1% stationarity gate. Common-point `U`
  relative-L2 changes decrease `2.88% -> 0.95%`; `W` changes decrease
  `12.54% -> 5.11%`. On the finest pair, all AIJ error/correlation/sign
  changes are below 10%, minimum lee velocity changes by only 0.60%, and
  projected divergence passes. The audit therefore records
  `accepted_component_grid`. Fine-grid mean `U/W` RMSE is
  `0.14536/0.06614`, correlation is `0.94802/0.89522`, and runtime is
  `186.3 s`. This acceptance is limited to the AIJ carrier-flow component;
  the deliberately effective-viscosity Reynolds number, hydrogen scalar,
  active thermodynamics and concentration grid gates are not promoted.
  Reproduction:
  `python tools/audit_aij_case_h_time_averaged_grid.py --grid-family aligned`.
- The active concentration failure has now been decomposed rather than hidden
  behind a looser gate.  On one fixed physical carrier field, the conservative
  species/enthalpy, EOS and variable-density mass-flux corrector converges at
  essentially first order: maximum-H2 RMS differences decrease
  `0.01801 -> 0.00894 vol%-point`, final H2 `0.001213 -> 0.000606 vol%-point`,
  and final temperature `0.002832 -> 0.001416 K`.  The independent D3Q19
  carrier-flow time audit is also monotonic.  The failure occurs only when the
  projected velocity is fed back into the isothermal BGK populations.
  Reproduction: `tools/diagnose_active_time_fixed_carrier.py` and
  `tools/diagnose_active_time_carrier_flow.py`.
- Four pressure-coupling variants were tested at `Sc_eff=0.7`, `Pr_eff=0.85`.
  Removing the preliminary incompressible face projection does not cure the
  failure, zero gravity does not cure it, and replacing instantaneous
  equilibrium reconstruction with a deferred Guo pressure impulse does not
  cure it.  Applying the variable-density projection to conservative scalar
  transport but not feeding it back to D3Q19 *does* pass the component time
  gate: maximum-H2 RMS differences decrease `0.01956 -> 0.01042`, final H2
  `0.003972 -> 0.002774 vol%-point`, and temperature
  `0.01093 -> 0.00770 K`.  That `transport_only` result is diagnostic evidence,
  not a promotable model, because its pressure-corrected velocity is not the
  next momentum state.  It identifies a structural incompatibility between
  the present isothermal BGK pressure dynamics and the required variable-
  density dilatation.  The active Tier-3 path remains withheld; the next
  implementation gate is a single, consistent variable-density low-Mach
  momentum/pressure integrator rather than another coupling coefficient.

- That integrator is now implemented as the opt-in
  `active_momentum_solver="finite_volume_projection"` candidate. D3Q19 is used
  only to initialise the carrier field; the active release then advances
  conservative cell momentum with the same projected face mass flux used by
  hydrogen and sensible enthalpy, followed by one variable-density pressure
  projection. The release momentum velocity is an explicit, audited input and
  cannot be silently used by the legacy BGK path.
- The first implementation exposed a time-step-independent face-filtering
  defect: reconstructing absolute face velocities from averaged cells every
  step was not the identity as `dt -> 0`. The corrected map preserves the
  previous projected faces and adds only the new cell-momentum increment. A
  zero increment is now an exact identity under obstacle and outer-boundary
  masks.
- A four-level fixed-carrier momentum audit gives approximately first-order
  time convergence. Successive relative-L2 differences decrease
  `2.88e-6 -> 1.40e-6 -> 6.89e-7` for U and similarly for V, W and speed;
  observed orders are `1.01--1.04`. Reproduction:
  `tools/diagnose_fv_momentum_time.py`.
- With species, enthalpy, EOS, buoyancy and the pressure corrector active,
  maximum-H2 RMS differences decrease `0.01794 -> 0.00891 vol%-point`, final
  H2 `0.001205 -> 0.000602 vol%-point`, and final temperature
  `0.002804 -> 0.001403 K`; all are first order. The aligned three-grid audit
  also passes its declared component gate. On the finest pair, maximum H2,
  4 vol% volume and windward/roof/side/lee peak changes are respectively
  `5.83%`, `4.08%`, `3.20%/0.90%/0.38%/3.59%`. Momentum-balance residuals are
  below `9e-15 kg m/s`; the fine case runs in `37.34 s`. Reproduction:
  `tools/audit_active_thermal_corrector_time.py` and
  `tools/audit_active_thermal_corrector_grid.py` with
  `--active-momentum-solver finite_volume_projection`.
- This is accepted **numerical component evidence**, not measured hydrogen
  validation. The candidate and `Sc_eff=0.7`, `Pr_eff=0.85` remain opt-in until
  an independent obstacle hydrogen experiment passes the external-validation
  gate. No default or quantitative facility claim is promoted yet.
- The HyTunnel public record is now reconciled against its published 40-case
  test matrix rather than treated as an unexplained 42-sheet mismatch. It
  contains 39 matching H2 cases, one declared cross-flow H2 case supplied as
  He, and two extra stagnant H2 baselines. A single matched reference case is
  frozen at 1 mm, 5 g/s and 3.5 m/s co-flow, including the official
  2.5 x 0.9 x 0.9 m CFD domain, 6.4 mm mesh, 7.25 mm notional nozzle and
  published 4 vol% experiment/CFD lengths of 1.65/2.25 m.
- The official boundary has been rechecked rather than inferred from the
  project name. The experiment was a fan-driven jet in a 6 m diameter vessel.
  The CFD prescribed 3.5 m/s air at `x=0` and open non-reflecting conditions
  on the remaining five faces; the 0.9 m square is not a no-slip tunnel.
- The high-pressure source route is now implemented without injecting the
  approximately `2000 m/s` atmospheric notional-nozzle state into the low-Mach
  solver. The existing HyRAM+/Yuceil--Otugen mass/momentum/energy nozzle is
  mixed with explicitly accounted ambient air until both physical Mach
  (`<=0.30`) and the active-scalar density-departure bound (`<=10%`) pass.
  For the frozen HyTunnel case the controlling density gate places the handoff
  `0.62069 m` downstream, with Mach `0.05588`, `-10%` density departure and a
  `0.18917 m` equivalent diameter. On the official `6.4 mm` grid this is
  `29.56` cells across. Hydrogen mass, full plane momentum and total energy
  close to round-off; the resolved source receives only the H2 mass plus the
  source increments of momentum and sensible enthalpy because entrained air is
  already part of the Eulerian carrier.
- This does not make the source definition identical to the official CFD.
  DEGALI's energy-conserving nozzle is `4.364 mm`, versus the report's
  `7.25 mm` Molkov notional nozzle. Diameter alone cannot reconstruct the
  latter's complete mass/momentum/energy plane, so both definitions are kept
  and no concentration-based diameter tuning is permitted.
- The resolved source now supports a cell-integrated anisotropic Gaussian:
  radial width matches the second moment of the equivalent uniform disk and
  axial width matches one uniform grid-cell slab. A four-cell diameter gate
  fails closed on under-resolved sources. The no-slip tunnel option remains
  available for genuinely confined facilities but is withdrawn from this
  matched case. The former all-face symmetry/free-slip surrogate remains a
  diagnostic. A new five-face pressure-open candidate now couples fixed
  ambient pressure, adjacent-plane velocity/non-equilibrium stress,
  multi-face mass projection and scalar outflow. Reflection and domain-size
  tests were then added: characteristic incoming-wave suppression reduces the
  1-D linear acoustic reflection ratio from 91.72% to 1.64% on both the
  streamwise and lateral-normal faces and passes the 5% gates. Oblique
  incidence at one declared 30° Fourier mode also passes at 4.188%. The
  expanded angle/wavelength sweep fails the 5% gate: the resolved long-wave
  mode gives 8.268% at 45° and 12.774% at 60°. A separate angle-aware
  characteristic candidate estimates incidence cosine from local normal and
  tangential acoustic velocity. On the unchanged matrix all seven resolved
  cases pass, with a worst reflection of 3.907% at 60°. This is component
  evidence only. A localised 45° packet aimed at the intersecting x-max/y-max
  faces passes at 2.594% using a conservative window beginning immediately
  after the direct pulse. A scaled x-y control and x-z/y-z edges pass at
  4.142%, 4.154% and 4.142%; the x-y-z three-face corner passes at 3.644%.
  A second 10-cell wavelength gives 6.076%, 6.073%, 6.075% and 5.999% for
  those four orientations and exposes a failure in that first candidate. A
  coefficient-free second-order-extrapolated angle-aware candidate gives
  2.418%/2.439%/2.420%/1.846% at 8 cells and
  3.192%/3.179%/3.191%/3.152% at 10 cells. Both declared 3-D matrices pass
  the 5% gate for that raw-extrapolation candidate, but it becomes unstable in
  a finite-vortex audit. A separately named TVD-limited local dilatation/curl
  candidate stays stable and gives 0.943% vorticity error, but velocity and
  density contamination are 19.825% and 243.230%; its 10-cell acoustic matrix
  is also 5.502--5.564%. No single pressure-open candidate passes both gate
  families. A first same-grid 0.120343 s plume-domain audit
  does pass: common maximum-field normalised L1 is 1.092%, worst primary-output
  change is 3.865%, and 4 vol% volume changes by 0.235%. The expanded run takes
  220.4 s with the explicit cached-direct pressure projection, down from
  477.6 s, and passes the 5 min development target without changing the
  concentration metrics. The official `0.707246 s`/3,056-step audit is now
  complete. Worst primary-output change is `0.9901%`, common maximum-field L1
  is `1.2447%`, and the enlarged-domain lateral clearance gate passes. The
  enlarged case takes `1402.2 s`, so the 5 min G4 target fails. This
  memory-intensive route remains opt-in. Unified vortical boundary treatment,
  source sensitivity and independent obstacle-H2 scoring remain open.
- A no-fit falsification run has been completed before spending the full 3-D
  cost. It joins the fixed source/handoff to the existing axisymmetric
  source-flux Gaussian jet and compares only the 108 published average
  concentration points downstream of the handoff. It obtains MAE/RMSE
  `2.688/3.160 vol%p`, bias `+2.687 vol%p`, Pearson `r=0.937`, MG `0.478`,
  VG `2.042` and FAC2 `0.519` at the declared `0.10 vol%` logarithmic floor.
  At 4 vol%, sensitivity is 1.0 but specificity is only 0.519 because there
  are 39 false positives. The model 4 vol% centreline length is 3.413 m versus
  the experiment's 1.65 m (`+106.9%`) and the official CFD's 2.25 m.
- This result rejects the **axisymmetric free-jet shortcut**, despite its high
  correlation. It neither scores nor rejects the resolved 3-D open-box model.
  No observation was used to fit the source, entrainment or spreading terms.
- The first resolved overlap score uses only the 54 measurements downstream of
  the 6% handoff and inside the 2.52 m common domain; 99 unsupported points are
  not extrapolated. A deterministic replay exactly reproduces the maximum
  slice and has only 0.00687% maximum receptor drift over its final 20%, so its
  endpoint is accepted as a steady proxy. Endpoint-versus-average MAE/RMSE are
  `2.179/2.491 vol%p`, `r=0.733`, FAC2 `0.556`, and 4 vol% sensitivity/
  specificity are `1.0/0.522`. This still fails the accuracy gates.
- A four-case, predeclared source/mixing sensitivity uses 45 points common to
  every case. At the 6% handoff, `Sc_eff=0.5/0.7/1.0` produces MAE
  `1.670/2.113/2.529 vol%p`; the 5% handoff at `Sc_eff=0.7` produces
  `2.185 vol%p`. The mixing closure dominates this handoff perturbation, but
  the lowest-error case still has FAC2 `0.689` and 4 vol% specificity `0.675`.
  It is diagnostic only: measured concentration may not select `Sc_eff` or a
  handoff threshold. Independent mixing evidence is required.
- A subsequent closure audit found that the sensitivity above used the base
  numerical momentum viscosity in a uniform `nu_base/Sc` surrogate rather
  than the local SGS eddy viscosity. It remains reproducible but is not
  physical evidence for selecting `Sc_t`. Schema v15 adds the separately named
  local closure `D_molecular + nu_sgs/Sc_t` (and the analogous thermal
  `alpha_molecular + nu_sgs/Pr_t`), conservative cellwise-diffusivity face
  fluxes, and active finite-volume SGS refresh from the evolving velocity.
  The first 100-step applicability check rejects the 6% handoff at 10.46%
  density departure and accepts the already declared 5% handoff at 9.25%; the
  decision uses no observed concentration. The full 3,056-step local-SGS case
  reaches a stationary, admissible endpoint but fails the external accuracy
  gates: on 45 common points MAE/RMSE are `3.425/3.976 vol%p`, bias is
  `+3.380 vol%p`, `r=0.396`, FAC2 `0.378`, and 4 vol% specificity `0.25`.
  Relative to the same-handoff legacy surrogate, MAE worsens by 1.239 vol%p.
  The outcome rejects local WALE alone on this coarse steady carrier and also
  prevents the legacy numerical compensation from being called validated
  turbulence. The next closure gate requires independently specified inlet/
  source turbulence or resolved unsteadiness plus scalar-flux evidence; it is
  not another fit of `Sc_t` to this sheet.
- Schema v16 can impose a complete `z-y` inlet plane on the active
  finite-volume momentum, pressure-correction and scalar path. The public
  HyTunnel fan workbook's first one-minute mean plane was transferred from
  measured `x=0.25 m` to the computational inlet without concentration-based
  fitting. The resulting full-duration score is effectively unchanged from
  uniform-in-y local SGS: MAE changes `3.425 -> 3.432 vol%p`, RMSE
  `3.976 -> 3.986 vol%p`, correlation `0.396 -> 0.392`, and 4 vol%
  specificity remains `0.25`. This rejects mean-profile nonuniformity as the
  missing correction. The workbook has no time series, RMS, integral scale,
  spectrum or dissipation, so the turbulence-inlet gate remains withheld.
- The public PRESLHY E3.4 package has now been audited as a local dynamic
  component dataset. Its four wind workbooks contain corrected pool mass and
  synchronized 10 Hz three-component sonic/H2 histories; the workbook flags
  all CH35 channels and the Sand02 wind channels defective. The remaining
  three workbooks can support conditional dynamic replay, but not obstacle-
  transport qualification: the vapour source flux is not measured directly,
  site obstacles are only approximately visible in a report figure, and the
  sensor record has nominal specifications rather than event calibration.
- A concentration-blind source-window policy now opens only the corrected
  mass sheet and retains every prominent post-fill decay. It freezes 13 source
  windows across the four records before any H2 observation or prediction is
  scored. The reviewed v2 artifact records the policy, raw-file hashes and all
  boundaries. Liquid-mass loss remains a conditional proxy rather than a
  qualified atmospheric source; porous substrates retain an artefact warning.
- Schema v17 adds a three-component time-dependent inlet contract to the
  finite-volume active momentum path. The vector history is interpolated on
  the event clock and carried through momentum, pressure correction and scalar
  transport. It deliberately cannot be used on the older D3Q19-active path,
  cannot reverse the current x-min inlet, and does not synthesize a spatial
  turbulence field. Before E3.4 replay, the sonic-to-site coordinate transform
  and treatment of instantaneous normal-flow reversal must be frozen without
  looking at H2 error.
- PRESLHY reproduction:
  `tools/audit_preslhy_e34_dynamic_evidence.py`,
  `tools/freeze_preslhy_e34_source_windows.py`,
  `outputs/preslhy-e34-dynamic-evidence-2026-10-11-v2/`, and
  `outputs/preslhy-e34-source-window-freeze-2026-10-11-v2/`.
- The first Sand03 dynamic component replay is now configuration- and
  scoring-policy-frozen before H2 access. Its 268.5 s, 26,850-step execution
  passes the numerical gates in 241.42 s, but fails every predeclared external
  accuracy target: median peak error is 21.53 vol%p, median dose error 124.4%,
  normalized arrival error 150%, and pooled 4 vol% sensitivity/specificity
  0.594/0.787. CH45/CH55 MAE are 4.280/2.305 vol%p. This is a genuine failed
  primary rather than a tunable sensitivity.
- Post-score diagnosis shows two structural evidence mismatches. The observed
  window starts at 15.86/3.52 vol% at CH45/CH55 because the post-fill interval
  omits plume history generated during filling, whereas the model is
  initialized at zero. Both model peaks occur at 109.96 s during the lowest
  inlet-U quartile; inverse-U correlation is 0.889/0.858 in the model but only
  0.085/0.057 in the observations. The next admissible work is a source-only
  prehistory sensitivity and independently supported spatial inflow, while
  preserving the failed primary and prohibiting observation-based time shift,
  source scaling or closure selection. E3.4 still cannot qualify obstacle
  transport.
- Dynamic replay reproduction:
  `tools/build_preslhy_e34_dynamic_replay_package.py`,
  `tools/preflight_preslhy_e34_dynamic_component.py`,
  `tools/score_preslhy_e34_dynamic_component.py`,
  `tools/diagnose_preslhy_e34_dynamic_component_score.py`, and
  `outputs/preslhy-e34-sand03-window01-preflight-2026-10-11-v1/`.
- A post-primary, non-validating fill-prehistory diagnostic now tests the
  initial-condition hypothesis without H2-based onset or rate fitting. The
  mass-only 10% fill crossing occurs at 308.12 s, 129.5 s before the primary
  window. Extending the frozen post-fill rate through that interval reduces
  median peak error 21.53 -> 14.63 vol%p and arrival error 150% -> 55%, but
  median dose error remains 116.3% and pooled 4 vol% sensitivity/specificity
  remain 0.590/0.792. This falsifies missing prehistory as the dominant error
  source. The next physics/evidence gate is spatial inflow plus source/mixing
  evidence, not another initial-condition or scalar coefficient fit.
- Fill-prehistory reproduction:
  `tools/build_preslhy_e34_fill_prehistory_diagnostic.py`,
  `tools/run_preslhy_e34_fill_prehistory_diagnostic.py`, and
  `outputs/preslhy-e34-sand03-window01-fill-prehistory-2026-10-11-v1/`.
- Thorney Island Trial 021 has now been recovered as a second obstacle-shape
  evidence package rather than forced through the generic SMEDIS reader. The
  dedicated reader retains 145 dose sensors (84 positive), six arc summaries,
  the workbook SHA-256, and the distinction between dose, peak concentration,
  arrival, departure and peak time. The final report takes precedence where
  it conflicts with the workbook: the impermeable semicircular fence has
  radius 50 m and height 5 m, event wind is 3.9 m/s and its heading is -6.1
  degrees relative to site +y. Workbook alternatives (6 m, 3.4 m/s and
  -15.8 degrees) are retained as declared sensitivity values rather than
  silently averaged. Physical fence thickness remains unreported.
- `CircularArcWall` is now a native solid in both the fast 3-D tier and the
  resolved 3-D tier. A concentration-blind three-grid geometry audit represents
  the reported wall as a one-cell radial band whose physical thickness shrinks
  under refinement. It remains one connected wall, leaves the source cell
  fluid, preserves the exact 5 m height, and reduces arc-length error from
  5.68% at 1 m to 1.22% at 0.5 m and 0.27% at 0.25 m. This passes the frozen
  geometry-component gate only; it is not a dispersion-validation result.
- Schema v18 now contains the required source-free full binary-mixture path.
  It transports conserved species mass with total projected mass flux,
  `-rho D grad(Y)` binary diffusion and conservative sensible enthalpy, recovers
  constant-pressure ideal-gas mixture density, and feeds that density and
  buoyancy into the finite-volume momentum/continuity path. It is separate
  from the bounded H2 added-scalar approximation.
- The twelve-sided, 14 m across, 13 m filled holder is represented as a regular
  dodecagonal prism using exact polygon/cell and vertical-overlap integrals.
  Interpreting the ambiguous word "across" as flat-to-flat gives 2048.204 m3,
  only -0.0876% from the reported 2050 m3. At 32 mol% R12 the reconstructed
  effective molar mass is 58.389 g/mol and relative density is 2.0158, matching
  the workbook 58.5 g/mol and report 2.02 without concentration fitting.
  Geometry volume and initial R12 mass are conserved to round-off on 1.0,
  0.5 and 0.25 m grids.
- A reduced coupled canonical step passes species-mass and total-continuity
  gates while exceeding 100% density departure, proving that the new path is
  not silently using the former 10% active-scalar limit. TI021 is therefore
  promoted only to `conditional` full-scale preflight.
- The molecular R12-air input is now independently frozen from Matsunaga,
  Hori and Nagashima (1994, DOI `10.1299/kikaib.60.504`):
  `D=4.77e-10 T^1.75 m2/s`, measured from 273.15--453.15 K at atmospheric
  pressure with reported +/-2% uncertainty. At 291.15 K this is
  `9.78865e-6 m2/s`. CoolProp 8.0.0 revision
  `ae81610e7d23efc57f9d051c8e70a4d66e87537f` gives R12 cp
  `606.0749 J/kg/K` at 291.15 K and 100500 Pa; the exact version and NIST
  REFPROP DOI are recorded. Molecular evidence does not select a turbulent
  scalar closure. Quantitative scoring remains withheld pending an
  independently frozen turbulent-mixing policy, full-scale grid/time/domain
  convergence and the no-fit dose/arc comparison. Treating the one-second
  workbook field as a continuous jet remains prohibited.
- The first report-scale execution exposed a real numerical incompatibility:
  physical air viscosity on a 1--2 m atmospheric grid puts the optional BGK
  relaxation effectively at 0.5. Schema v19 therefore separates the
  finite-volume carrier initial condition from D3Q19. Its `projected_inlet`
  path maps the declared vector wind around impermeable faces and projects it
  to discrete continuity while preserving `1.5e-5 m2/s` for subsequent
  finite-volume momentum. It is explicitly an unconverged initial condition,
  not a wake solution or turbulence closure.
- A concentration-blind 1 m full-site single-step preflight now covers
  291,200 cells, the exact 5 m wall height and the 2048.204 m3 initial cloud.
  It completed in 32.4 s with AMG-CG; R12 mass residual was at round-off, maximum
  local continuity residual was `3.700e-11 kg/m3/s`, all correctors converged
  in at most 11 iterations, and every declared preflight gate passed. This
  promotes only the execution evidence to
  `accepted_full_scale_single_step_preflight`. The TI021 package remains
  `conditional_pending_turbulent_closure_and_no_fit_run`; no dose, arc,
  long-horizon or grid/time/domain accuracy claim follows from 0.001 s.
- A cached symmetric smoothed-aggregation hierarchy now preconditions the
  variable-density CG projection. A frozen same-case 2 m comparison reduces
  summed projection iterations from 4,313 to 107 (40.31x) and total runtime
  from 4.40 to 1.53 s (2.87x), while maximum mole fraction and density
  departure match within `1e-10` and both continuity residuals remain below
  `1e-7 kg/m3/s`. This is an accepted solver-acceleration component, not a
  relaxation of a conservation gate or a physical-validation result.
- A non-promoting 2 m time-step diagnostic advances the complete coupled
  report-scale domain for 1 s. `dt=0.05 s` reaches momentum CFL 0.878 and
  runs in 34.93 s; `dt=0.025 s` reaches CFL 0.441 and runs in 60.06 s. Their
  final maximum R12 mole fractions differ by 0.0001783 (0.0178 vol%-point),
  and both retain conservation and corrector convergence. This is not the
  required three-level, full-field/receptor monotonic time audit, so the
  smaller step is only the next candidate rather than an accepted time grid.
- The turbulent-mixing policy is now frozen before any TI021 concentration,
  dose or arc score is read. Schema v20 adds a neutral Monin--Obukhov
  background `nu_t=kappa*u_star*z` using the event evidence
  `u_star=0.205 m/s` and `z0=0.005 m`. The primary prior is
  `Sc_t=0.7`, `Pr_t=0.74`, but the complete Cartesian envelope
  `Sc_t={0.2,0.7,1.3}` by `Pr_t={0.74,1.0}` is mandatory and all six outputs
  must be retained. The reported D/E class does not supply a measured
  Monin--Obukhov length, so no stable correction is invented. This closure
  provides atmospheric background mixing but no obstacle-generated turbulence
  production. A 1 m, 291,200-cell primary-member coupled preflight passes all
  geometry and numerical gates in 34.0 s. The evidence status is therefore
  `conditional_pending_full_transient_convergence_and_no_fit_run`; the
  full-transient three-level grid/time/domain audit and no-fit dose/arc score
  remain mandatory.
- A sensor-coordinate-only coverage audit then invalidated use of the small
  preflight box as a full validation domain: it contains 33/145 sensors, while
  the network reaches 500 m downwind and +/-300 m crosswind. The frozen
  200 x 200 x 40 m resolved inner domain contains the complete fence and all
  six 75/100 m arc centreline observations; the 600 x 640 x 40 m outer domain
  contains all 145 sensors. At the measured solver cell-step rates, a 600 s,
  1.25 m monolithic outer-domain run requires at least 1,730 h per prior
  member. A monolithic six-member/three-grid campaign is therefore withheld.
  The next architecture gate is an overlapping conservative handoff of species
  mass, sensible enthalpy and face mass flux from the resolved inner tier to
  the outer transport tier, with <=1e-6 integral flux mismatch, <=1e-7
  kg/m3/s projected continuity residual and <=10% overlap peak/dose change.
  The package status is now
  `conditional_pending_conservative_nested_handoff_and_full_transient`.
- The first nested operator now uses exact rectangular face overlaps to remap
  time-dependent total-mass anomaly, species-mass, sensible-enthalpy and
  lateral/vertical momentum flux densities from the inner y-z plane to the
  wider outer plane. It retains signed flux and local backflow without clipping
  or post-hoc normalisation. The TI021-shaped analytic 1.25 m to 5 m audit
  preserves all five integrated flux histories to `1e-12` relative tolerance.
  This is accepted numerical remap evidence only. Time interpolation into the
  outer finite-volume boundary, addition of the independently declared ambient
  flux and an overlap-domain peak/dose comparison are still required. The
  package therefore advances only to
  `conditional_pending_outer_boundary_coupling_and_full_transient`.
- Binary-mixture corrector schema v2 now accepts a complete prescribed x-min
  total/species/enthalpy flux contract. The remapped history is linearly
  interpolated, added to a separately declared ambient total mass flux and
  inverted from flux ratios to boundary composition, specific enthalpy,
  temperature, density and three velocity components. The x-min total mass
  flux remains fixed while projection corrects the downstream pressure-open
  face. A 5 m, 8,192-cell one-step audit closes species to `1e-10 kg`, sensible
  enthalpy to `1e-6 J` and local total continuity to `1e-7 kg/m3/s`. This is an
  accepted outer binary boundary component, not the full handoff: transverse/
  vertical momentum application and an overlapping inner/outer peak-dose
  comparison remain. The package status is
  `conditional_pending_outer_momentum_overlap_and_full_transient`.
- Binary-mixture corrector schema v3 and the nested outer execution component
  now apply the reconstructed `u/v/w` to the finite-volume momentum predictor,
  apply the same conserved total/species/enthalpy inlet to the pressure/
  thermodynamic corrector, and feed the projected faces back to the next time
  step. The 48-step same-5 m narrow-reference/wide-outer comparison passes all
  conservation and corrector gates; maximum peak/dose change over eight common
  receptors is `7.35e-9`/`2.15e-8`. This isolates nesting from grid
  convergence. The separately retained 2.5-to-5 m diagnostic fails at 18.81%
  peak and 45.12% dose, so the 5 m grid remains unqualified. The evidence state
  is `conditional_pending_full_transient_convergence_and_no_fit_scoring`; the
  architecture permits the 145-sensor run, but does not claim it has run.
- Full-transient plan schema v2 corrects the physical coordinate join before
  any long run. The original 600 x 640 x 40 m box remains a monolithic
  sensor-coverage/runtime reference, not the nested outer grid. The actual
  handoff is site-local downwind y=100 m, represented at inner solver x=160 m
  and outer x=0. The resulting outer domain is 440 x 640 x 40 m and spans
  local y=100--540 m. Inner and outer cover 60/94 sensors, share 9 overlap
  sensors, and cover 145/145 in union without reading dose or concentration.
- The actual-inner streaming path is now exercised rather than inferred from
  synthetic planes. It records every total/species/enthalpy/lateral-momentum/
  vertical-momentum frame at inner x=160 m to append-only HDF5, reopens it and
  area-remaps it to the outer `q=-site_x` plane with `1.94e-16` maximum
  integrated residual. The HSE report's p.150 stationary-cloud description is
  now the explicit TI021 initial-momentum contract. The report-based dx=5 m
  run completes 600 s and 60 inner sensors in 501.1 s with mass, energy,
  continuity and nonlinear-corrector gates passed. Its no-fit score is a failed
  accuracy diagnostic, not a promotion: six-arc peak MAPE is 47.8%, dose MAPE
  is 81.5%, and detected arrival MAE is 25.2 s. The one-cell fence and aliased
  near-ground sensor heights make 2.5/1.25 m grid checks mandatory before any
  quantitative claim. Evidence v10/v13 therefore reports
  `conditional_pending_outer_transient_and_grid_convergence`.
- Thorney reproduction:
  `tools/audit_thorney_ti021_evidence.py`,
  `tools/audit_thorney_ti021_curved_wall_geometry.py`,
  `tools/audit_thorney_ti021_initial_cloud.py`,
  `tools/audit_thorney_ti021_r12_properties.py`,
  `tools/preflight_thorney_ti021_full_scale.py`,
  `tools/audit_thorney_ti021_projection_acceleration.py`,
  `tools/audit_thorney_ti021_turbulent_mixing_policy.py`,
  `tools/audit_thorney_ti021_full_transient_plan.py`,
  `tools/audit_thorney_ti021_nested_handoff.py`,
  `tools/audit_thorney_ti021_outer_inlet_coupling.py`,
  `tools/audit_thorney_ti021_outer_momentum_overlap.py`,
  `tools/audit_thorney_ti021_handoff_streaming_preflight.py`,
  `tools/run_thorney_ti021_inner_transient.py`,
  `tools/score_thorney_ti021_inner_transient.py`,
  `outputs/thorney-ti021-evidence-audit-2026-10-11-v13/`,
  `outputs/thorney-ti021-handoff-streaming-preflight-2026-10-11-v1/`,
  `outputs/thorney-ti021-inner-transient-dx5m-600s-quiescent-2026-10-11-v2/`,
  `outputs/thorney-ti021-inner-transient-score-dx5m-2026-10-11-v1/`,
  `outputs/thorney-ti021-full-transient-plan-2026-10-11-v2/`,
  `outputs/thorney-ti021-nested-handoff-audit-2026-10-11-v1/`,
  `outputs/thorney-ti021-outer-inlet-coupling-audit-2026-10-11-v1/`,
  `outputs/thorney-ti021-outer-momentum-overlap-audit-2026-10-11-v1/`,
  `outputs/thorney-ti021-turbulent-mixing-policy-2026-10-11-v1/`,
  `outputs/thorney-ti021-full-scale-preflight-dx1m-mixing-2026-10-11-v1/`,
  `outputs/thorney-ti021-full-scale-preflight-dx1m-2026-10-11-v1/`,
  `outputs/thorney-ti021-projection-acceleration-audit-2026-10-11-v1/`,
  `outputs/thorney-ti021-initial-cloud-audit-2026-10-11-v1/`, and
  `outputs/thorney-ti021-curved-wall-geometry-audit-2026-10-11-v1/`, and
  `outputs/thorney-ti021-r12-property-audit-2026-10-11-v1.json`.
- Reproduction: `tests/test_high_pressure_source_handoff_3d.py`,
  `tests/test_validate_hytunnel_free_jet_handoff.py`,
  `outputs/hytunnel-ps-spatial-audit-2026-10-11-v23.json`,
  `outputs/hytunnel-resolved-transient-score-2026-10-11-v2.json/.png`, and
  `outputs/hytunnel-free-jet-handoff-comparison-2026-10-10.json`.
- A separate mapped-source preflight prevents the 3-D run from inheriting an
  invalid mean-plane assumption. On a declared 45 mm coarse grid, the nominal
  10% top-hat density-departure handoff produces a 13.342% maximum departure
  after cell-integrated Gaussian mapping. Tightening the mean-plane gate to 8%
  still gives 11.562%; 6% gives 9.664%, zero out-of-range cell-steps and
  round-off mass balance. No measured concentration or CFD output is read.
  Thus 10% and 8% are rejected for this grid, while 6% is only a numerical
  candidate for grid/time/domain sensitivity—not a promoted physical value.
  These values use the characteristic pressure-open candidate; oblique
  planar non-reflection and the raw-extrapolated 3-D multi-wavelength linear-
  acoustic matrix are component-qualified, but no candidate also passes the
  nonlinear vortical gate;
  plume-domain independence is demonstrated only for the 0.120343 s horizon.
- Preflight reproduction:
  `tools/audit_hytunnel_resolved_source_preflight.py` and
  `outputs/hytunnel-resolved-source-preflight-2026-10-10-v6.json`.
- Canonical reflection reproduction:
  `tools/audit_pressure_open_boundary_reflection.py` and
  `outputs/pressure-open-boundary-reflection-audit-2026-10-10-v2.json/.png`.
- Single-angle oblique reflection reproduction:
  `tools/audit_pressure_open_boundary_oblique_reflection.py` and
  `outputs/pressure-open-oblique-reflection-audit-2026-10-10.json/.png`.
- Angle/wavelength sweep reproduction:
  `tools/audit_pressure_open_boundary_sweep.py` and
  `outputs/pressure-open-angle-wavelength-sweep-2026-10-10.json/.png`.
  Add `--boundary-scheme characteristic_multidimensional_extrapolated` to reproduce
  `outputs/pressure-open-multidimensional-sweep-2026-10-10.json/.png`.
- Intersecting-face corner reproduction:
  `tools/audit_pressure_open_boundary_corner_reflection.py` and
  `outputs/pressure-open-corner-reflection-audit-2026-10-10.json/.png`.
- Finite-vortex outflow reproduction:
  `tools/audit_pressure_open_vortical_outflow.py` and
  `outputs/pressure-open-vortical-outflow-audit-2026-10-10.json/.png`.
- Z-edge and three-face corner reproduction:
  `tools/audit_pressure_open_boundary_3d_corner_reflection.py` and
  `outputs/pressure-open-3d-corner-reflection-audit-2026-10-10.json/.png`.
  Add `--pulse-wavelength-cells 10` and distinct output paths to reproduce
  `outputs/pressure-open-3d-corner-reflection-audit-2026-10-10-wavelength10.json/.png`.
  Add `--boundary-scheme characteristic_multidimensional_extrapolated` and
  the matching wavelength to reproduce the two passing
  `-extrapolated-wavelength8` and `-extrapolated-wavelength10` artifacts.
- Short-horizon plume-domain reproduction:
  `tools/audit_hytunnel_plume_domain_sensitivity.py` and
  `outputs/hytunnel-plume-domain-sensitivity-2026-10-10-v3.json/.png`.
