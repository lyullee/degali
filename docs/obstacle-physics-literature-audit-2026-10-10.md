# Obstacle-flow and hydrogen-dispersion literature audit

Date: 2026-10-10

## Audit question

Which physical mechanisms and numerical models are required before DEGALI can
make a quantitative concentration claim for hydrogen striking and passing a
finite obstacle?

## Direct comparison

| Study or dataset | Equations/evidence | Observed or resolved behaviour | DEGALI legacy semi-FV | New Tier 1 | Requirement carried forward |
|---|---|---|---|---|---|
| Olvera and Choudhuri, *International Journal of Hydrogen Energy* 31 (2006), DOI 10.1016/j.ijhydene.2006.02.022 | 3-D finite-volume mass, momentum and energy; standard k-epsilon; building-flow validation | Upper-windward-edge separation, roof cavity, lee cavity and reattachment; hydrogen buoyancy changes impact | No pressure, energy, turbulence or recirculation | Continuous over-roof flow only | Tier 2/3 must solve pressure-coupled momentum, energy and separated wake |
| Henriksen et al., *International Journal of Hydrogen Energy* (2024), DOI 10.1016/j.ijhydene.2024.10.038 | Experiments with 29 H2 sensors plus LES in a 5.8 m channel; cases with/without a 36.8% blockage | Time histories, ceiling accumulation, buoyant vertical gradient and obstruction response | No channel flow, ceiling buoyancy or sensor-matched velocity field | No thermal/buoyant feedback | Use as the first hydrogen obstacle time-series benchmark after Tier 2 thermal coupling |
| Akhter and Kaye, *Fluids* 10 (2025), DOI 10.3390/fluids10020046 and 10.3390/fluids10080194 | Water-flume LIF experiments patterned on Thorney Island; continuous and finite upstream dense releases | Windward impact, delayed lee response, Richardson-dependent over/around flow and wake trapping | Teleports windward flux to roof and cannot trap a wake | Restores local windward impact but has no trapping | Tier 2 must reproduce distinct low-/high-Ri regimes and lee residence time |
| Davies and Singh, *Journal of Hazardous Materials* 11 (1985), DOI 10.1016/0304-3894(85)85043-3 | Full-scale Thorney Island dense-gas obstruction trials | Finite cloud interaction with solid/porous barriers and isolated building | No resolved velocity or finite-width bypass | Two-dimensional over-roof path only | Tier 3 required before field-scale isolated-building comparison |
| Nielsen, *Journal of Loss Prevention in the Process Industries* 4 (1991), DOI 10.1016/0950-4230(91)80004-E | Continuous liquefied-propane field experiments with a transverse wall | General ground-level dilution but elevated increases caused by cloud lifting | A single positive roof relocation and no signed spatial response | Can lift locally but cannot separate/mix turbulently | Validation must score every sensor; no scalar wall factor is admissible |
| AIJ Case H, DOI 10.5281/zenodo.15429872 | Wind-tunnel U,V,W, RMS velocity, TKE, mean tracer and concentration RMS around a 1:1:2 block | Neutral single-building flow and tracer wake with a near-ground source on the leeward side | No predicted wake velocity/TKE | No recirculating wake/TKE | Tier 2 neutral-flow validation; not direct LH2 validation |
| Tominaga et al., *Journal of Wind Engineering and Industrial Aerodynamics* 96 (2008), DOI 10.1016/j.jweia.2008.02.058 | Cross-code/experiment CFD benchmark guidance | Grid, inflow turbulence and numerical choices materially affect building-flow accuracy | No flow-field validation | Only potential-flow consistency checks | Freeze grid, inlet, observation operator and evaluation data before calibration |
| Angot, Bruneau and Fabrie, *Numerische Mathematik* 81 (1999), DOI 10.1007/s002110050427 | Navier--Stokes/Brinkman penalisation with convergence analysis | Cartesian representation of no-slip obstacles through vanishing permeability | Not used | Not used | Selected basis for Tier 2/3 solid momentum treatment |
| Kolomenskiy and Schneider, *Journal of Computational Physics* 231 (2012), DOI 10.1016/j.jcp.2012.02.009 | Volume penalisation for incompressible flow and scalar advection--diffusion | No-slip flow and no-flux scalar with immersed obstacles | Non-local target-cell routing | Explicit zero solid-face scalar flux | Use no-flux penalisation/conservative cut-face audit in Tier 2/3 |
| NIST Chemistry WebBook SRD 69, molecular hydrogen | Molecular weight and temperature-dependent ideal-gas thermochemistry | Reference enthalpy/heat-capacity basis and stated temperature ranges | No thermal state | No thermal state | Do not extrapolate the current constant-cp development closure into a validated cryogenic property model |
| Lu and Pantano, *Journal of Computational Physics* 404 (2020), DOI 10.1016/j.jcp.2019.109132 | Variable-density zero-Mach mass constraint, pressure solvability and manufactured solutions | Temperature/density and velocity cannot be selected independently; discrete compatibility controls accuracy | No mixture continuity equation | Divergence-free constant-density volume flux only | Enforce density evolution, source mass and boundary mass flux in the pressure equation before strong-density-ratio claims |
| Nicoud and Ducros, *Flow, Turbulence and Combustion* 62 (1999), DOI 10.1023/A:1009995426001 | WALE eddy viscosity from the square of the velocity-gradient tensor | Responds to strain and rotation and recovers the proper near-wall cubic scaling without an explicit damping distance | No wall-aware SGS model | No viscous SGS model | Replace undamped Smagorinsky as the next wall/wake candidate; verify before promotion |
| Weickert et al., *Computers & Mathematics with Applications* 59 (2010), DOI 10.1016/j.camwa.2009.08.060 | WALE implemented in an LBM framework; channel and wall-mounted-cube tests through Re=40,000 | Demonstrates an algebraic wall-adapting LBM route and transitional-to-turbulent obstacle tests | No LBM flow | No LBM flow | Use the paper's channel/cube sequence as the implementation-verification order |
| Longo et al., *Building and Environment* 185 (2020), DOI 10.1016/j.buildenv.2020.107066 | Multi-fidelity local turbulent-Schmidt modelling across single-building, array and canyon cases | Concentration accuracy depends on spatially varying scalar/momentum mixing rather than one universal scalar diffusivity | One prescribed scalar diffusivity | Local molecular-plus-SGS candidate and explicit legacy surrogate | Do not tune a single Schmidt number to the current roof/side grid audit; qualify a local scalar-flux closure against external data |
| Bae, Stanford Report TF-163 (2016), filtered scalar transport formulation | SGS scalar flux is closed by gradient diffusion with eddy diffusivity `D_e = nu_e/Sc_e`, while molecular and SGS diffusion remain separate | Direct basis for `D_eff = D_molecular + nu_sgs/Sc_t` | Uniform `nu_base/Sc` surrogate | Spatial molecular-plus-SGS field, with legacy route retained only for reproduction | Re-run concentration evidence without selecting `Sc_t` from the evaluation data |

## Governing-model decision

### Rejected

1. **Non-local target-cell redistribution.** It conserves the inventory but
   violates the local transport path and creates a first-roof-cell spike.
2. **A scalar obstacle multiplier.** SMEDIS and field literature show both
   increases and decreases depending on position.
3. **A prescribed wake factor without a velocity solution.** It cannot meet
   the windward, roof and lee observations simultaneously.
4. **Calling a potential-flow contour CFD validation.** Potential flow has no
   no-slip separation, turbulent production or recirculating wake.

### Selected

1. Tier 1: sparse streamfunction/potential-flow solve plus local conservative
   scalar finite volumes and zero solid-face flux.
2. Tier 2: transient incompressible Navier--Stokes with pressure projection
   and Brinkman penalisation, followed by species and thermal transport.
3. Tier 3: three-dimensional low-Mach/variable-density extension with finite
   cuboids, lateral bypass and multiple obstacles.

## Dimensionless controls that must be reported

- obstacle Reynolds number, `Re_H = U_H H / nu`;
- source and building Froude numbers;
- source Richardson number or reduced-gravity Richardson number;
- Schmidt and turbulent Schmidt numbers;
- Prandtl and turbulent Prandtl numbers;
- blockage ratio and obstacle aspect ratios;
- source distance, height and diameter normalised by obstacle height;
- atmospheric stability and roughness/inlet boundary-layer description.

Agreement obtained at one unreported combination of these quantities is not
transferable to a different hydrogen or cryogenic regime.

## Current implementation result

`degali.addons.fast_obstacle_transport` now supplies Tier 1.  It solves a
sparse Laplace problem for the obstacle-dependent streamfunction.  Face
velocities are a discrete curl, so their cell divergence is zero to roundoff.
Scalar advection and diffusion are local and every fluid/solid face is forced
to zero normal flux.  The result exposes pressure coefficient, receptor time
histories, concentration in kg/m3 and vol%, mass ledger, divergence and wall
clock time.

In the first 12 m by 6 m regression, the windward receptor peak exceeded the
roof peak, which exceeded the lee-ground peak.  Maximum divergence was below
`3e-15 s-1`; cumulative mass residual was below `8e-17 kg/m` for a `0.006
kg/m` finite source.  This corrects the non-local roof spike but remains
conditional because potential flow cannot satisfy the resolved-wake evidence.

## Next implementation gate

Tier 2 may start only with canonical verification cases and independent
diagnostics for:

1. pressure-projection divergence reduction;
2. no-slip decay inside a Brinkman solid;
3. Poiseuille-channel convergence;
4. lid-driven cavity vortex location;
5. bluff-body recirculation and reattachment under grid refinement;
6. no-flux passive-scalar mass conservation;
7. Boussinesq sign and energy/species conservation.

Only after those pass should AIJ velocity/TKE/concentration observations be
used.  Hydrogen obstacle concentrations remain withheld until the thermal and
buoyancy terms pass the channel experiment comparison.

## Tier-2 prototype checkpoint

The first viscous prototype now solves the isothermal low-Mach D2Q9
lattice-Boltzmann hydrodynamic limit around a no-slip ground-mounted block.
Unlike the legacy routing rule, the lee recirculation is generated by the
momentum solution. A separate Helmholtz projection makes the face field used
by the finite-volume hydrogen solver divergence-free in every non-outlet fluid
cell. The scalar operator retains exactly zero obstacle-face flux.

For the matched 16 m by 6 m diagnostic case, the prototype produced 648
negative-streamwise-velocity cells in the defined lee region, a minimum lee
velocity of -0.333 m/s, core divergence below `2e-13 s-1`, hydrogen mass
residual below `3e-15 kg/m`, and a 4.90 s runtime. The legacy roof-routing
maximum was 42.54 vol%; the Tier-1 and Tier-2 maxima were 7.40 and 6.76 vol%
respectively and remained near the source rather than appearing in the first
roof cell.

These are numerical-development observations, not accuracy claims. The
prototype still lacks the energy equation, composition-density feedback,
buoyancy, a turbulence closure, finite-width lateral bypass, grid/time
convergence and comparison against AIJ, SMEDIS/Thorney Island or hydrogen
obstacle measurements. Its quantitative facility flag therefore remains
false.

## Canonical and grid-convergence checkpoint

The common D2Q9 equilibrium, BGK collision, Guo forcing and no-slip
bounce-back were independently exercised in a periodic Poiseuille channel.
The numerical velocity profile produced 0.220% relative L2 error, 0.161%
relative maximum error and 0.240% mass-flux error against the analytic
parabola. This verifies the laminar channel building block; it does not verify
the obstacle inlet/outlet boundaries or hydrogen dispersion.

The obstacle calculation was then repeated on four matched grids from 0.20 m
to 0.08 m cells at one physical flow duration. Receptor observation changed
from nearest-cell sampling to bilinear interpolation over fluid cells only.
For the finest pair, peak changes were 4.28%, 2.85%, 4.47% and 5.76% at the
windward, roof, lee-low and lee-shear receptors. Corresponding dose changes
were 4.81%, 4.37%, 9.27% and 4.66%. The 4 vol% envelope area changed 5.10%
and its maximum downwind position changed 9.88%.

The overall numerical gate remains withheld: roof duration above 4 vol%
changed from 0.222 s to 0.402 s, a 44.721% difference. The code records this
failure rather than promoting the otherwise converged peak and dose results.
The reproducible audit is generated by
`tools/verify_resolved_obstacle_numerics.py`.

## AIJ Case H external-flow checkpoint

The official AIJ Case H Zenodo files are now hash-pinned and read by a strict
loader. The 17-point measured approach profile is applied at the inlet; a
declared `(z,U)=(0,0)` point completes the unmeasured 0--0.01 m ground layer.
The evaluation uses all 220 centre-plane points having both mean streamwise
and vertical velocity. No Case H output is used to fit viscosity or another
coefficient.

The conditional run gives correlations of 0.828 for `U/U_H` and 0.721 for
`W/U_H`, but normalised RMSE values remain 0.486 and 0.434. Recirculation-sign
accuracy is 79.5%. These correlations show that the model recognizes the
gross building-flow pattern; the magnitude errors are far too large for a
quantitative validation claim.

The failure is physically informative. The finite AIJ building permits
sideways bypass, while the present two-dimensional block occupies the full
span and drives excessive over-roof acceleration and vertical motion. The
stability-qualified effective-viscosity run also has `Re_H=32.3` rather than
the experimental 43,000, contains no inlet turbulence or turbulence closure,
and did not reach its steady residual tolerance. Fitting the effective
viscosity to the same evaluation data would not repair these missing physics
and is prohibited. AIJ concentration and TKE remain withheld. This checkpoint
therefore promotes the finite-width Tier-3 momentum solver and independent
turbulence verification to the next implementation gate.

## Finite-width geometry isolation

To separate dimensional blockage from viscosity/turbulence errors, a
three-dimensional potential-flow baseline was added without calibrating to
Case H. It preserves the measured inlet profile and the published
0.10 x 0.10 x 0.20 m finite building, while allowing both lateral and roof
bypass. A conservative local 3-D scalar operator is attached for later fast
screening, but concentration is disabled in this flow-only evaluation.

Relative to the 2-D viscous prototype, the finite-width baseline reduces the
normalised streamwise RMSE by 26.46% and vertical RMSE by 68.52%, while running
in about 0.81 s. This is strong evidence that treating the building as an
infinite transverse wall materially corrupts the centre-plane magnitude.
It is not evidence that potential flow is an adequate wake model: it has slip
walls, no separation or negative-velocity cavity, and its `W` correlation is
only 0.629. The comparison therefore fixes the model architecture: retain the
3-D potential solver as a rapid geometry screen, but require finite-width
no-slip momentum plus turbulence for the resolved tier.

## D3Q19 finite-width viscous checkpoint

The next resolved prototype now uses the standard D3Q19 equilibrium and
weights, BGK collision, link bounce-back no-slip solids, lattice pressure and
an optional explicit Smagorinsky--Lilly strain viscosity. The coefficient is a
declared input (`0.12` in the AIJ development run), not fitted to the Case H
outputs. The scalar operator receives a separately projected divergence-free
face field; the reported validation velocity remains the unprojected D3Q19
momentum field so that scalar bookkeeping cannot silently alter the quantity
being compared.

Before external comparison, the same D3Q19 equilibrium, Guo body force and
bounce-back links were checked against analytic 3-D Poiseuille flow. Relative
L2, maximum and mass-flux errors were respectively 0.712%, 0.525% and 0.776%.
The finite-cuboid regression then showed simultaneous lateral bypass, roof
flow, 143 negative-streamwise lee cells, zero solid-face scalar flux and a
hydrogen mass residual below `1e-18 kg`.

The no-fit AIJ comparison materially improves: on 195 common observation
points, `U/U_H` RMSE/correlation become 0.169/0.934 and `W/U_H` become
0.086/0.808. The matched 2-D errors are 0.512 and 0.435. This is the first
evidence in the rebuild that finite width and viscous wake physics improve the
same independent flow field together.

Promotion is still prohibited. The run is at effective `Re_H=32.3`, not the
experimental 43,000; the grid and flow duration have not converged; the final
residual is 0.00227; 25 near-solid points are excluded by the coarse trilinear
stencil; the inlet/outlet boundary treatment and SGS wall behaviour have not
been independently qualified; and no energy or density/buoyancy feedback is
present. AIJ concentration, TKE and every cryogenic-hydrogen inference remain
withheld pending those gates.

## Spatial and flow-horizon audit

The three-grid audit uses spacings of 0.025, 0.020 and 0.0125 m at the same
0.2786 s physical flow horizon. It compares the same 176 AIJ points on every
grid. The finest pair changes normalized `U` MAE by 11.28%, `W` RMSE by
12.08%, and minimum lee velocity by 23.13%, all above the predeclared 10%
limit; the finest residual is 0.00167. Spatial convergence is therefore
`withheld` even though the finest-grid `U` and `W` RMSE values improve to
0.145 and 0.0679.

The separate time audit uses the 0.020 m grid for 0.1393, 0.2786 and 0.5573 s
physical horizons and retains 195 common points. Between the two longest
runs, normalized `U` RMSE changes 0.64%, `W` RMSE 1.83%, both correlations
less than 0.07%, recirculation-sign accuracy 0%, and minimum lee velocity
4.58%. These metric-sensitivity checks pass the 10% limit, but the 0.5573 s
residual is still 0.00160, so the flow is not declared converged. The audit
status remains `withheld` and no concentration evidence is promoted.

Machine-readable records and plots are written to
`outputs/obstacle-physics-rebuild-2026-10-10/` by
`tools/audit_aij_resolved_3d_grid.py` and
`tools/audit_aij_resolved_3d_time.py`.

## Thermal and buoyancy coupling checkpoint

The first coupled path now transports hydrogen mass concentration and signed
sensible-enthalpy density as separate conservative finite-volume quantities.
The recovered ideal-gas H2/air density supplies a vertical Guo force to every
subsequent D3Q19 momentum step. This is genuine scalar-to-momentum feedback,
not a post-processed plume lift, and mass and sensible-energy ledgers are
reported independently.

The canonical neutral-density map passes both sign checks. At the declared
ambient state, the 4 vol% neutral-buoyancy temperature is 282.240 K: colder
mixture states are dense and receive downward acceleration, while warmer
states are light and receive upward acceleration. Yellow applicability lines
in the generated map mark the 10% relative-density bound of the active-scalar
approximation.

The finite-block development case now uses the same 3-D obstacle and source for
a neutral run and a mass-flux-corrected active run. The current active result
retains hydrogen-mass residual below `7e-18 kg` and sensible-energy residual
below `2e-11 J`; all 250 fixed points converge in at most five iterations. The
worst local variable-density continuity residual is `7.68e-15 kg/m3/s` and the
worst global residual is `2.93e-15 kg/s`. After reusing the immutable
projection matrix and factorization, observed local runtime is about 6.3 s
(18.17 s before the optimizations) versus about 1.1 s for the one-way path.
The neutral and active 4 vol% maximum-field volumes are both 0.816 m3, while
the largest local maximum-field change is 2.340 vol%-point. This demonstrates
that the coupling is live; it does not establish physical accuracy.

Promotion remains prohibited because the thermodynamic carrier representation
is still a bounded constant-cp active-scalar approximation. Continuity,
species, energy, EOS density and pressure correction are now discretely linked,
but temperature-dependent properties and real-gas/phase-aware source-state
handoff remain missing. External hydrogen obstacle concentration data must be
scored before the coupled path is used quantitatively.

## Variable-density continuity-operator checkpoint

The first replacement for the constant-volume constraint now projects face
mass flux rather than velocity. Its target divergence contains both the
discrete density time derivative and the local total-mass source. A zero-
potential positive-x outlet adjusts its boundary mass flux, which avoids the
previous practice of excluding the outlet cell from the continuity norm.

Random-flux, uniform-expansion, localized-source/obstacle and nonuniform-front
manufactured cases all pass. The worst cell residual after projection is
`1.03e-13 kg/m3/s` and the worst global residual is `1.69e-14 kg/s`. The
localized-source case maintains zero mass flux across every obstacle face;
sealed fluid components without access to the pressure outlet fail explicitly.

This is accepted numerical evidence for the projection operator only.

## Coupled low-Mach corrector checkpoint

The time step now iterates hydrogen transport, signed sensible-enthalpy
transport, ideal-gas mixture reconstruction and the variable-density mass-flux
projection. The returned face velocity and pressure correction are reconciled
with the D3Q19 momentum populations without discarding the post-collision
non-equilibrium stress.

Two independent manufactured cases are frozen by
`tools/verify_low_mach_thermal_corrector.py`. Ambient uniform flow remains
unchanged and converges in the required two iterations. A local 80 K hydrogen
source upstream of a finite block converges in three iterations, closes local
continuity to `1.04e-14 kg/m3/s`, closes global balance to `8.18e-16 kg/s`,
conserves species and sensible enthalpy to printed precision, and retains
exactly zero normal velocity on solid faces. The machine status is
`accepted_numerical_operator_only`.

The resolved scalar flux now defaults to MUSCL-Hancock reconstruction with a
van-Leer limiter; first-order upwind remains available as a declared comparison
path. In a smooth Gaussian translation canonical, RMS error decreases from
0.02794 to 0.00501 (82.1%) without violating mass conservation, non-negativity
or the initial upper bound. Reconstruction slopes are zeroed when their stencil
touches a solid, preserving the obstacle no-flux condition.

The resolved tier remains `withheld`. The coupled step still needs a temporal
and grid sensitivity study with active thermodynamics, stronger property/source
canonical tests, and external cryogenic-hydrogen obstacle validation. The AIJ
velocity benchmark is valuable for the carrier flow but cannot qualify
hydrogen thermodynamics or concentration error.

The first three-level time audit is now available. It fixes the 0.20 m grid,
7.2 s physical flow spin-up, 3.0 s transient and source history, then halves
the coupled step from 0.012 to 0.006 to 0.003 s. All 1,750 fixed points converge
and pass local/global continuity gates. With bounded MUSCL transport, the
domain maximum approaches 21.239 vol% H2 and the 4 vol% volume is identical
(0.800 m3) on the medium and fine
runs. Nevertheless, the full-field RMS differences plateau rather than
decrease: 0.0300 -> 0.0335 vol%-point for maximum H2, 0.0342 -> 0.0356
vol%-point for final H2, and 0.0723 -> 0.0752 K for final temperature. The
machine status is therefore `withheld_time_refinement_not_monotonic`. The small
absolute sensitivity is useful for engineering triage, but it is not a formal
order-of-accuracy result; spatial, boundary and LBM relaxation errors must be
separated in the next coupled grid/time audit. Optimized coarse/medium/fine
runtimes are 6.16, 12.35 and 24.43 s.

The first active joint grid/time audit now uses a 1.5 refinement ratio with the
finite-block boundaries exactly aligned on all three grids. Cell size decreases
from 0.300 to 0.200 to 0.1333 m while lattice Mach, physical spin-up, source
history and transient duration remain fixed. Conservation and corrector gates
pass on every grid. With bounded MUSCL transport, the finest-pair maximum H2
and 4 vol% volume changes are 5.86% and 0.48%, and the windward sensor
peak/dose changes are 5.84%/5.21%.
Those global and windward metrics meet the 10% development tolerance, but the
roof peak/dose changes by 28.12%/34.53% and the side peak/dose by
25.10%/31.48%. Successive full-field RMS differences increase rather than
decrease. The grid status is `withheld_grid_refinement_gate_failed`.

This result localizes the remaining numerical weakness: coarse cells provide
too much cross-stream/top diffusion, so low-concentration bypass predictions
are not grid independent even where the windward hazard metric is relatively
stable. The higher-order flux materially improves the 4 vol% volume and
windward finest-pair sensitivity relative to upwind, but it does not cure the
low-concentration bypass dependence. Reusing a cached sparse factorization
reduces the finest observed runtime from 352.26 to 131.88 s. Further promotion
requires a finer asymptotic study and independently qualified wall/wake and
turbulent scalar-flux closures; tolerances must not be relaxed to hide this
behavior.

The next turbulence implementation gate is therefore WALE, not an arbitrary
retuning of the present Smagorinsky constant. WALE must first reproduce its
near-wall limiting behaviour and pass channel plus wall-mounted-cube tests in
the D3Q19 formulation. Only after that may the active obstacle grid audit be
repeated. Turbulent scalar diffusivity remains a separate gate: a fixed or
local Schmidt closure cannot be selected from the current no-observation grid
study, because that would turn numerical sensitivity into calibration.

WALE has now been implemented as an opt-in candidate and passes the local
invariant checks, but the first external carrier-flow comparison does not
justify promotion. On the same `0.020 m` AIJ Case H grid and 195 common mean
velocity observations, no SGS, Smagorinsky (`Cs=0.12`) and WALE (`Cw=0.50`)
give normalized `U` RMSE `0.1680`, `0.1685` and `0.1684`; normalized `W` RMSE
is `0.0860`, `0.0863` and `0.0863`. WALE's four MAE/RMSE reductions relative
to Smagorinsky are only `0.02--0.15%`, while runtime rises from `16.9 s` to
`32.4 s`; no SGS is marginally best. All three flow residuals are still above
the convergence tolerance. This is consistent with a coarse, deliberately
high-effective-viscosity development case in which an SGS term has little
room to improve resolved dynamics. It is not evidence against WALE at a
Reynolds-consistent resolution, but it is decisive evidence not to make WALE
the default now. Reproduction:
`tools/compare_aij_case_h_subgrid_models.py`.

## Open-boundary and mean-observation audit

Zou and He (1997, DOI `10.1063/1.869307`) show that equilibrium-only flow
boundaries can introduce significant error and establish pressure/velocity
boundary canonical tests. Guo, Zheng and Shi (2002, DOI
`10.1088/1009-1963/11/4/310`) decompose a boundary population into prescribed
equilibrium and an adjacent-fluid non-equilibrium extrapolation. That NEE
operator was implemented without coefficient fitting and its imposed velocity
and density moments are exact on interior boundary nodes.

Izquierdo and Fueyo (2008, DOI `10.1103/PhysRevE.78.046707`) show why exact
fixed-pressure moments are insufficient: a subsonic outlet must suppress the
incoming Euler characteristic while allowing the outgoing one to leave. Their
one-dimensional pressure-wave test is mirrored here with a short domain and a
no-return long-domain reference.

The present obstacle case falsifies a generic promotion, however. On the
`0.020 m`, no-SGS AIJ run, NEE increases `U` MAE/RMSE by 5.97%/2.78% and `W`
MAE/RMSE by 3.67%/4.16%. Its final 100-step velocity-change residual is
`0.04288`, versus `0.00229` for the equilibrium/copy baseline, and the peak is
at the last cell centre. Reduced raw flux imbalance and density variation do
not compensate for worse wake and observation error. NEE therefore remains a
rejected opt-in candidate rather than a default selected from general theory.

A downstream-fetch audit then separates outlet proximity from the observation
region. Extending the domain from `0.70` to `0.90` to `1.10 m`, with fixed
`0.020 m` cells and equal nominal flow-through count, changes all longest-pair
AIJ metrics by less than 1% and the lee minimum by 2.55%. The observation
region therefore passes a 5% outlet-distance sensitivity gate, even though the
instantaneous residual does not monotonically vanish.

Finally, AIJ reports mean velocities, whereas the previous DEGALI comparison
used one final instantaneous field. The new observation operator performs a
declared spin-up and averages a finite trailing window. For a 3000-step no-SGS
run, the last 1000-step mean differs by `0.8015%` between its two halves, below
the predeclared 1% stationarity tolerance. Mean `U/W` RMSE is
`0.16861/0.08493`; the mixed changes relative to the shorter instantaneous
case are small and not uniformly favourable. This is accepted only as
component stationarity/observation-operator alignment. It does not cure the
failed spatial refinement gate or validate hydrogen concentration.

The time-averaged refinement audit then showed that the old failure was partly
geometric. At `dx=0.020 m`, both streamwise block faces crossed cell centres;
the inclusive mask occupied six x layers and represented `0.12 m` rather than
the declared `0.10 m`. Coarse and fine grids represented the correct width, so
the three solutions did not share one geometry. Half-open voxel membership and
explicit effective-bound diagnostics now make this discrepancy visible.

Repeating the audit on `dx=0.025, 0.0166667, 0.0125 m`, for which every block
face lies on a cell face, changes the conclusion. Common-point model-field
relative-L2 differences decrease `2.88% -> 0.95%` for U and
`12.54% -> 5.11%` for W. Every finest-pair error/correlation/sign metric is
within 10%, lee minimum velocity changes 0.60%, and all statistical mean-field
residuals are below 1%. This supports an `accepted_component_grid` result for
the low-Re/effective-viscosity AIJ carrier-flow calculation. It still cannot
be used as evidence for Reynolds-consistent turbulence, scalar mixing or
hydrogen concentration.

## Active pressure--momentum coupling diagnosis

The individual numerical components now separate cleanly.  A fixed-carrier
test gives approximately first-order time convergence for species, enthalpy,
the ideal-mixture EOS and variable-density mass-flux projection.  The D3Q19
carrier alone is also time-monotonic.  In contrast, every tested route that
feeds the variable-density pressure correction back into the isothermal BGK
momentum state remains non-monotonic in final H2 and temperature: the original
equilibrium replacement, variable-density projection without a preceding
incompressible projection, zero-gravity isolation, and a deferred Guo pressure
impulse.

A diagnostic `transport_only` mode retains cellwise variable-density
continuity for species/enthalpy fluxes but does not alter the next D3Q19
momentum state.  It passes the three-level component time gate.  This is useful
causal evidence but is intentionally not promoted: it breaks the desired
pressure-to-momentum feedback.  The evidence points to a formulation mismatch,
not a tunable Schmidt number or source-splitting defect.  The isothermal BGK
step tends toward its own near-incompressible lattice-pressure state, while the
thermodynamic continuity equation requires non-zero mixture expansion and
contraction.  A production active solver therefore needs one consistent
variable-density low-Mach momentum and pressure discretisation; stacking a
second pressure projection on the existing BGK pressure dynamics is no longer
an accepted development route.

## Unified finite-volume low-Mach candidate

The replacement follows the conservative projection structure of Almgren et
al. (1998, DOI `10.1006/jcph.1998.5890`), the time-accurate variable-density
predictor/corrector formulation of Lessani and Papalexandris (2006, DOI
`10.1016/j.jcp.2005.07.001`), and McDermott's mass/energy/EOS-consistent
low-Mach divergence constraint (2014, DOI `10.1016/j.jcp.2014.06.019`). These
references motivate the equation structure; they do not validate DEGALI's
obstacle or hydrogen predictions.

The opt-in candidate advances conservative momentum, species and sensible
enthalpy with one shared projected face mass flux. Viscosity and buoyancy are
included in the provisional momentum step and one variable-density projection
supplies pressure. D3Q19 supplies only the initial carrier field, avoiding two
competing active pressure dynamics. A declared source-momentum velocity is
included in the conservative balance; a non-zero value is rejected on the
legacy BGK active path rather than ignored.

Manufactured balances and four-level temporal refinement pass. The fully
coupled time audit is first order, and an aligned three-grid audit passes its
predeclared component thresholds with sub-6% finest-pair changes in every
reported H2 peak and the 4 vol% volume. This resolves the numerical coupling
failure identified above. It does **not** resolve turbulence-closure or
experimental-validation uncertainty: external obstacle H2 sensor/contour data
remain the promotion gate.

## High-pressure source to low-Mach obstacle handoff

The atmospheric notional-nozzle construction already follows the declared
HyRAM+/Yuceil--Otugen sequence: measured mass fixes the discharge coefficient,
pressure thrust contributes to momentum, and the atmospheric state closes
total energy. That plane is not itself a permissible low-Mach inlet. In the
matched HyTunnel case it is approximately `2000 m/s`, so inserting it into the
resolved obstacle grid would violate the model equations even if a small time
step prevented numerical failure.

The new adapter performs an explicit entraining control-volume reduction. It
adds ambient air with its declared co-flow momentum and energy, conserves H2
mass, plane momentum and total energy, and selects the first plane satisfying
both Mach and receiving-EOS density-departure limits. The distance follows the
same source-momentum entrainment closure used by the existing axisymmetric jet;
it is not fitted to HyTunnel concentration. The finite-width source mapping is
also deterministic: its radial Gaussian matches the second moment of the
equivalent uniform disk and its axial Gaussian matches a one-cell uniform slab.

The matched result is controlled by the `10%` active-scalar density bound, not
the Mach limit: distance `0.62069 m`, diameter `0.18917 m`, Mach `0.05588` and
relative density departure `-0.10`. This makes the limitation explicit. DEGALI
does not resolve the earlier compressible region or concentrations above the
receiving closure's density range; its reported downstream coordinate must
retain the handoff offset.

Confinement is a separate boundary choice. `all_no_slip_tunnel` uses full
population reversal at sidewalls and roof and the same half-cell viscous
no-slip gradient in the finite-volume momentum predictor. It is appropriate
to test actual confined facilities, but it is not the HyTunnel boundary. The
official HyTunnel CFD prescribed an inlet velocity at `x=0` and open non-
reflecting conditions on all five remaining faces. The implemented
`all_symmetry_free_slip` mode removes the ground drag but remains impermeable
and is retained only as a diagnostic surrogate. The new `all_pressure_open`
candidate applies adjacent-plane velocity/non-equilibrium stress and fixed
ambient pressure on the five non-inlet faces, opens the variable-density
projection on those faces, removes wall viscous drag, and accounts for scalar
and enthalpy outflow. Uniform-state and conservation tests pass. The first
fixed-pressure implementation nevertheless reflects 91.72% of the incident
linear density-pulse amplitude and therefore fails the predeclared 5% gate.
Replacing its macroscopic target with the isothermal incoming/outgoing Riemann
invariants reduces the ratio to 1.64% at both the streamwise and independently
isolated lateral-normal faces, passing both 5% tests. One declared 30°
tangential Fourier mode also passes at 4.188%, but the broader one-dimensional-
impedance sweep fails at 45° and 60° long waves. A separate angle-aware
candidate uses the local normal/tangential acoustic velocity ratio and passes
all seven resolved sweep cases, with a worst value of 3.907% at 60°. One
localised x-max/y-max corner packet passes at 2.594%. The first angle-aware
candidate passes the scaled 8-cell x-y/x-z/y-z/xyz matrix at
4.142%/4.154%/4.142%/3.644%, but fails the 10-cell matrix at
6.076%/6.073%/6.075%/5.999%. A coefficient-free second-order-extrapolated
candidate passes both matrices at 2.418%/2.439%/2.420%/1.846% (8-cell) and
3.192%/3.179%/3.191%/3.152% (10-cell). Axis rotation and the declared linear-
acoustic matrix are supported for that raw extrapolation only. It becomes
unstable for a finite Gaussian vortex. A TVD-limited local dilatation/curl
candidate is stable and keeps vorticity error to 0.943%, but fails its
velocity/density contamination limits at 19.825%/243.230% and its 10-cell
acoustic matrix at 5.502--5.564%. No unified pressure-open candidate is
qualified. A same-grid
0.120343 s plume-domain audit passes with 1.092% common-field L1 and 3.865%
worst primary-output change. The official 0.707246 s audit also passes its
same-grid domain gate at 1.2447% common-field L1 and 0.9901% worst primary
change; its expanded runtime is 1402.2 s, so the separate 5 min workflow gate
fails. This domain result does not qualify the failed nonlinear vortical
boundary component.

### No-fit free-jet route falsification

Before a full resolved open-box run, the fixed source and handoff were joined to
the existing conserved axisymmetric free-jet route. The model coordinate was
shifted so that its source-flux establishment plane begins at the physical
handoff coordinate; this changes neither state nor flux and uses no measured
concentration. Only 108 published steady-average points downstream of that
plane were scored. Upstream measurements, maximum/minimum matrices, finite-box
boundary and inlet-profile effects, and unavailable raw transient signals were
excluded explicitly.

The outcome demonstrates why correlation is not an accuracy claim. Pearson
`r=0.937`, yet MAE/RMSE are `2.688/3.160 vol%p`, the signed model-minus-data
bias is `+2.687 vol%p`, and FAC2 is only `0.519`. All 27 observed 4 vol% points
are detected, but 39 additional points are false positives, giving precision
0.409. The predicted centreline 4 vol% length is 3.413 m, compared with 1.65 m
in the experiment and 2.25 m in the official CFD. The 10 vol% length is
1.280 m versus 0.60 m measured. The route therefore fails the predeclared FAC2
and threshold-length gates and is marked `diagnostic_component_route_rejected`.

This rejection is scoped to the axisymmetric shortcut. It does not prejudge
the resolved 3-D equations, which must be run with the official open boundary
or a qualified domain-insensitive surrogate, identical observation operator,
grid/time audits and both declared source models. The result is recorded in
`outputs/hytunnel-free-jet-handoff-comparison-2026-10-10.json` and its derived
contour figure; model execution took 2.51 s on the current workstation.

### Resolved-source profile admissibility

The first coarse resolved-source preflight exposed a distinction that the
mean-plane handoff alone cannot test. The nominal plane satisfies a 10%
density-departure bound, but finite-width Gaussian mapping concentrates the
local scalar source. With the characteristic pressure-open candidate, a 45 mm
grid and 50 source steps produce a 13.342% maximum density departure for the
10% mean-plane case. An 8% candidate reaches 11.562%; a 6% candidate reaches
9.664% with no out-of-range
cell-steps. Corrector convergence and mass balance pass in all three cases.

This is a numerical admissibility screen, not coefficient calibration: it
reads no measurement or CFD result. Its `0.04 m2/s` coarse-grid viscosity is
explicit numerical regularisation, the carrier flow resolution is not
independently qualified, and the 6% result may change with grid, time step,
duration and domain. The nominal
10% mapped source therefore fails closed, while the 6% candidate may proceed
only to those sensitivity studies. The audit is
`outputs/hytunnel-resolved-source-preflight-2026-10-10-v6.json`. The declared
planar angle/wavelength sweep exposed a 10-cell failure in the first angle-
aware candidate. The raw second-order-extrapolated candidate passes both
declared 8-cell and 10-cell 3-D matrices but is unstable under finite-vortex
passage. Its stable TVD-limited counterpart still fails velocity/density and
10-cell acoustic gates. A unified open boundary remains a separate
prerequisite. The official-duration domain artifact is
`outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json`.
It passes numerical domain sensitivity but its 1402.2 s enlarged run misses
the 5 min development target.

### Resolved transient overlap score

The 6% density-departure handoff begins at `x=1.08625 m`; the common resolved
box ends at 2.52 m. Therefore only the 54 public average-concentration points
at x=1.25--2.50 m are supported. The 81 upstream and 18 beyond-domain points
are excluded rather than extrapolated. A deterministic replay reproduces the
frozen maximum slice at normalised L1=0 and gives only 0.00687% maximum
receptor drift over the final 20%, qualifying its endpoint as a steady proxy.
Against the experimental average, MAE/RMSE are `2.179/2.491 vol%p`, bias
`+2.170 vol%p`, `r=0.733` and FAC2 `0.556`. The 4 vol% confusion matrix is
TP/TN/FP/FN = 8/24/22/0, hence sensitivity 1.0 and specificity 0.522. It still
fails the concentration gate. Reproduction:
`tools/score_hytunnel_resolved_transient.py` and
`outputs/hytunnel-resolved-transient-score-2026-10-11-v2.json/.png`.

The next no-fit sensitivity separates effective scalar mixing from the
receiving-model handoff. All four cases are scored on the same 45 points.
At the 6% handoff, `Sc_eff=0.5/0.7/1.0` gives MAE
`1.670/2.113/2.529 vol%p`; moving to the more conservative 5% handoff at
`Sc_eff=0.7` gives `2.185 vol%p`. Therefore the current discrepancy is more
sensitive to the unvalidated gradient-diffusion closure than to that handoff
change. The apparent best perturbation is not adopted: FAC2 is only `0.689`
and 4 vol% specificity `0.675`, and observation-based coefficient selection
would be calibration leakage. Artifact:
`outputs/hytunnel-source-mixing-sensitivity-2026-10-11-v1.json/.png`.
