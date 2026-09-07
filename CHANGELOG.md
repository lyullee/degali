# Changelog

## Unreleased - GitHub alpha snapshot

- Renamed the project, Python package and command from `degadisx` to `degali`.
  DEGALI means **Dense Gas Dispersion for Liquid Hydrogen**. This is a naming
  migration only; it does not change the frozen physical or numerical results.
- Prepared a compact public repository while retaining the full local research
  ledger. Third-party PDFs, raw workbooks, checkpoints and large generated
  fields are excluded without deleting them.
- Replaced the development-log front page with a concise status, installation,
  evidence and safety overview. The former front page is preserved under
  `docs/README-development-log-2026-09-06.md`.
- Added publication, contribution and security policies plus an automated
  repository size, restricted-file and credential check.
- The current package remains `0.1.0.dev0`. No research-only LH2 candidate was
  promoted and no physical or numerical result was changed by packaging.

## Research history

- Added joint PSD bounds showing that TKE and axial normal-stress work cannot
  both be assumed small for the frozen7 initial shear fields; all7 resolved.
  Added explicit normal momentum/work source allocation. Preserved fixed-center
  Q2 failure; new exact-geometry six-DOF mobile initialization passes both .02/2
  synthetic cases without changing source targets or fitting Q/r (max1.489e-11).
  Mobile source6 tests passed,105 hashes match. Added straight-frame stress/
  production transformation primitives; no full normal/pressure closure or
  observed/default adoption (2026-09-06,05:47 KST).

- Completed the explicit-Q source energy reallocation audit: six retained
  moments close to5.566e-12 while Q is paid from the unchanged total source
  energy. Original projection preserved;95 hashes match. Four new source
  tests passed separately AFTER full425. Synthetic-only, no actual Q selected
  or physical initial-boundary acceptance (2026-09-06,05:10 KST).

- Implemented opt-in finite Q=rho*k modal transport with thermal dissipation
  D, turbulent production P-D, and total H+meanKE+Q energy. All original scalar
  rows retained. Manufactured phase-refinement and actual energy gates pass
  (max1.992e-7), but its physical input gates FAIL and it is not adopted.
  Full non-slow425 passed. Added explicit-Q boundary source energy reallocation
  helper for separate testing. [Results](docs/finite-tke-transport-results.md)
  (2026-09-06,05:00 KST). Defaults/field scores unchanged.

- Completed exact-geometry conservative aligned witness re-solve and its new
  actual-energy differences (max7.66e-8); all retained equations and independent
  boundary gates pass. Extended necessary shear/TKE bounds to all7 initial
  cases (10.51–11.72% of mean axial KE flux), all numerically resolved.
  Next finite-TKE operator specified, not yet implemented. No physical closure,
  downstream observation or default promotion (2026-09-06,04:30 KST).

- Added opt-in exact-constraint transverse geometry with consistent tangents;
  actual full-energy directional differences pass all four registered checks,
  max7.18e-8. The old direction is not reused as a new conservative solution.
  Added sharp PSD shear/TKE lower-bound diagnostic (no fitted coefficient):
  minimum TKE flux10.56% of mean axial KE for the fixed prior aligned witness.
  Non-slow403 passed;3 later bound tests passed separately. No field/default
  adoption. [Results](docs/exact-geometry-and-tke-results.md) (2026-09-06).

- Added bounded-memory phase operators and a direction-preserving four-mode
  circulation feasibility pilot: failed trial10 edges3.543/4.000/1.963%,
  nonradial stress1.5e-15, all original weak/global rows retained. Not a closure.
  Preserved failed raw/paired/root-polished energy checks; new high-precision
  smooth enthalpy diagnosis separates cancellation and width/wind inconsistency.
  No field/default change. Non-slow385 passed;7 later unit tests separately pass.
  [Results](docs/aligned-flux-and-energy-conditioning-results.md) (2026-09-06).

- Completed seven guarded new-shape segment tests:6/7 pass short intervals,
  trial10's corner heat5.01189% failure independently reproduced. Added
  rate-free compatibility bounds, a rejected boundary-row prototype and
  conservative divergence-free flux modes retaining all original weak rows.
- Hard boundary-moment solenoidal pilots fail; separate small-amplitude
  inequality FEASIBILITY witness satisfies sampled gates but needs independent
  phase-volume verification and a physical flow/stress closure. No downstream
  adoption or observed score. Full non-slow371 passed;3 later LP tests pass.
  [Results](docs/conservative-flux-witness-results.md) (2026-09-06).

- Extended simultaneous initialization to all7 unchanged cases; all pass
  independent local gates. Added degree6/37-state operator coverage.
  New guarded research segment driver evolves fields without endpoint
  refitting and requires explicit mixing-amplitude scaling.12 guard/state
  tests pass; full non-slow355 passed,128 skipped,58 deselected578.82s.
  Downstream tests ongoing; trial10 heat boundary failure is preserved.
  No observed score/default promotion.
  [Results](docs/coupled-initialization-extension-results.md) (2026-09-06).

- Added research-only simultaneous shape/rate initialization: vectorized
  fixed-mesh iterations re-solve all modal rates for every candidate and
  difference sample, with hard six-moment constraints. Final candidates
  are retracted/verified by separate moving-phase quadrature and65 rays.
- Both preselected pilots25/11 pass local gates, including the newly added
  momentum boundary5% requirement. Trial25 scalar defects fall to3.24/2.30%;
  trial11 no longer requires negative inferred viscosity on checked samples.
  No fitted diffusivity, new flow redistribution, default or observed score.
  Added6 tests;342 passed,128 skipped,58 deselected. Other five pending.
  [Results](docs/coupled-shape-initialization-results.md) (2026-09-06).

- Removed redundant thermal-width/shape coordinates in a new research path;
  added23/37-state conditional modal transport with newly reconstructed
  mass/axial-momentum fluxes, stress work, angular scalar diffusion and weak
  reservoir boundaries. The supplied scalar mixing profile remains explicit.
- Preserved the original4/7 numerical screen and independently verified the
  fine mixing input for a combined7/7 numerical result. Constitutive adoption
  fails all seven: renewed edge mismatches and two negative-viscosity cases.
  Higher-order witnesses locate those signs near finite-section corners.
  Added12 tests; non-slow336 passed,128 skipped,58 deselected. No default,
  new downstream segment or field score. [Results](docs/enriched-shape-transport-results.md).

- Added face/phase-intersection angular splitting, rebuilt phase-ray
  quadrature, six-moment derivatives and hard-constraint scalar shape
  refitting. All seven frozen cases now pass independent section gates.
  Internal edge defects are0.95–3.63%, not observed prediction errors.
- Kept prior failed evidence, original shape degrees, log correction limit,
  transport, physical parameters and core/defaults unchanged. Added7 tests;
  non-slow:324 passed,128 skipped,58 deselected. No new shape-coefficient
  transport law, field score or default promotion.
  [Results](docs/edge-conservative-refit-results.md).

- Added opt-in square-symmetric C/H shape enrichment with fixed transport,
  six nonlinear moment constraints, exact phase inversion and independent
  edge/moment verification. Three of seven cases reach both edge gates;
  zero pass all conservation gates, so no candidate is adopted.
- Retained both basis degrees and failures; added a separate phase-cell
  polar integration diagnosis without changing original coefficients.
  Angular reference reproduction still needs refinement in five cases.
  Added16 tests; non-slow:317 passed,128 skipped,58 deselected.
  No new field score, downstream shape transport or default changes.
  [Results](docs/edge-profile-enrichment-results.md).

- Added opt-in ambient-reservoir total-energy flux boundaries and a reduced
  buoyancy-work ledger for weak thermal zeroth/second moments. Added a free
  thermal-width RHS and short research segment driver without changing the
  core, original boundary guards, source coefficients or default options.
- All seven initial weak balances and independent physical flux differences
  pass. Original 2/4/8-step short-segment checks pass 6/7; retained trial 10's
  fixed-step failures and separately verified it with adaptive DOP853.
  Combined short-segment checks pass 7/7, not pointwise/field validation.
  Edge species/heat gradient defects remain 12–25%. Added 14 tests; non-slow:
  301 passed, 128 skipped, 58 deselected. [Results](docs/reservoir-thermal-moment-results.md).

- Added reduced axial shear-stress/work reconstruction and an explicit
  thermal equilibrium compatibility screen. All seven numerical audits pass;
  unity thermal/species diffusion plus immediate shear heating fails the
  independent zeroth heat balance in all seven frozen boundary cases.
- Separated finite-edge advection/diffusion and ruled out merely reducing
  the heat-conversion fraction while keeping the current five source terms.
  Preserved failed physical candidates. Added 14 tests; non-slow regression:
  287 passed, 128 skipped, 58 deselected. No core/default or field-score change.
  See [the results](docs/shear-thermal-compatibility-results.md).

- Added conditional conservative transverse mass/species reconstruction,
  actual-source/coflow replay, a free thermal-width tangent family and an
  inferred anisotropic mixing tensor. Negative diffusivity is rejected;
  thermal/species ratio and width rate are never silently supplied.
- Added monotonic phase-cell-split quadrature and independent finite-change
  verification of all seven boundary tangent families. Retained failures of
  the dense unsplit gradient reference. Added 16 tests for diffusion,
  moving/curved sections, phase splitting and invalid/underdetermined inputs.
  Non-slow: 273 passed, 128 skipped, 58 deselected. No default promotion.

- Added research thermal-shape operators: mean advective enthalpy second
  moment and reduced derivative, planar curved/moving rectangular balance,
  and variable-density specific-enthalpy diffusion response. No diffusivity,
  mechanical heating distribution or downstream beta_H closure is invented.
- Preserved a failed seven-case fixed-order phase-gradient quadrature audit;
  added separate adaptive volume/edge integration with explicit error and
  nonconvergence reporting. Gaussian diffusion and curved manufactured
  balances plus negative input/closure tests are covered by 32 new tests.
  Adaptive verification passes 7/7 fixed sections, including the net response;
  no observed accuracy score or downstream transport closure is claimed.
  Non-slow suite: 257 passed, 128 skipped, 58 deselected. No default change.

- Added a boundary-only independent enthalpy width constrained by five fluxes
  and the near-field buoyancy moment. All seven fixed boundaries pass both
  conservation/force and existing temperature/width gates, with independent
  adaptive GK21 confirmation. No new field accuracy score or promotion.
- Added implicit phase/moment derivatives, fixed-flux warm starts and a
  force-cancellation rejection regression. Unclosed downstream transport and
  inappropriate five-constraint use fail explicitly. Pinned the recovered
  scientific environment in a project-local virtual environment.

- Added an opt-in Gaussian volumetric-enthalpy profile with strict local C,H
  phase inversion, a dedicated conservative downstream inverse and actual
  non-Gaussian density buoyancy integration. The seven-case / 38-arc /
  17-profile rerun is complete. Same-41-sensor minimum-temperature MAE drops
  about 5%, but height and median-temperature errors worsen; not promoted.
- Added exact polar integration over the original finite square, checking
  all five fluxes instead of energy alone for the new profile. Retained the
  original failed numerical checkpoint and rejected it for subsequent runs.
  Added high-order downstream audits, thermodynamic-profile/quadrature merge
  guards and hashed reuse of frozen temperature reductions.

- Added opt-in explicit-species ambient/phase EOS consistency, including RH
  conversion and an exact ambient table node. No-H2 ambient gas no longer
  acquires a spurious approximately -1.20 K temperature offset in this mode.
  The 7-interface / 38-arc / 17-profile rerun retains the existing gates and
  defaults. Centre-height MAE changes from 0.06359 to 0.06159 m; concentration
  changes are small and vertical width worsens slightly, so no promotion.
- Added frozen-trajectory replay with temperature, five-flux and sensor-arc
  identity checks before comparing raw thermocouple observations. Workbook
  hashes and fixed sampling windows are retained without modifying sources.
  Field merging now rejects incompatible ambient closures as well as source
  physics/isomers, and newly generated field shards include input/code hashes.

- Corrected the measured-LH2 HEM research source's pipe-to-Gaussian energy
  ledger: carry incoming pipe kinetic energy and explicit condensed-phase
  enthalpy with a checked component reference. Independent upstream checks
  remove a 5.99–6.22% energy mismatch on four measured nozzles; this is not an
  experimental accuracy claim. Historical HEM field results remain archived.
  The completed 38-arc/17-profile rerun reduces HEM centre MAE from 0.0700
  to 0.0636 m, but concentration and the joint promotion gate do not improve.
- Added stable pure-H2 HEOS volume/caloric departure screening, deliberately
  separate from a mixture EOS. The 20 K source requires nonideal treatment;
  multiplying the full H2/air mixture volume by pure-H2 Z is not adopted.
- Preserved the crosswind tensor quadrature while merging exactly duplicate
  radial exponents by symmetry. Added source-physics/isomer guards to field
  shard merging, pointwise field output and downstream trajectory storage.

- Added a spin-consistent LH2 source/near-field selector.  The pure-para
  PRESLHY end member preserves normal hydrogen as the default and passes all
  seven conservative interfaces. It reduces the fast-bound centre-height MAE
  from 0.0700 to 0.0545 m, but worsens VG, FAC2 and vertical width, so it is
  retained as an explicit composition uncertainty rather than promoted.
- Added a measured-pipe LH2 flash/droplet source that conserves mass,
  pressure-thrust momentum and total energy while retaining separate liquid
  and vapour rates.  It reproduces the published ESREL/TNO micrometre
  atomisation scale and fixes `C_ds=10/15/20` as a non-fitted sensitivity.
- Implemented PRESLHY D3.1 equation 45 for the physical GASFLOW-MPI phase
  relaxation coefficient and inverse critical-diffusivity/diameter screens. The
  measured-source droplets would require effective diffusivity below 0.071%
  of the 20.4 K NBS H2 reference value, or droplets 37.6--52.9 times larger,
  to survive the collective source transit, so no empirical phase-delay
  coefficient is introduced.
- Added a phase-correct homogeneous-equilibrium evaporation-source bound and
  connected it to the seven-trial PRESLHY conservative interface/field path.
  It passes all interfaces and improves VG and vertical width, but worsens MG,
  FAC2 and centre-height MAE, so the validated default is unchanged.
- Made the downstream five-flux state inversion retry strict failed solves
  from neighbouring physical width/velocity states without relaxing its
  tolerance.  This recovers the formerly fragile trial-23 path while leaving
  established solutions unchanged; the complete non-slow suite passes.

- Added a conservative axisymmetric-to-JETPLU LH2 handoff with transferred
  phase thermodynamics, adaptive energy quadrature, explicit crosswind
  entrainment selection, and reproducible PRESLHY coupled validation. The
  local density-scaled shear transition repairs the diagnosed vertical-width
  loss, while the fully four-flux path now fails closed where JETPLU needs an
  independent thermal-profile state.
- Pre-registered the coefficient-free seven-state crosswind extension that
  separates centre density from centre H2 fraction and transports total energy
  alongside mass, species and vector momentum.
- Implemented its first fail-closed gate: independent density/composition
  Gaussian profiles and an adaptive-quadrature five-flux boundary solve. All
  seven frozen PRESLHY 10D interfaces pass, including the two former
  single-scalar failures.
- Added the downstream five-balance ODE using conservative-flux RK4 with a
  positive state inversion at every stage. The 0.02/0.01 m PRESLHY runs are
  converged and close long-range balances below `1.04e-7`; VG, FAC2 and width
  improve, but MG and centre-height MAE fail the frozen promotion rule.
- Added a common PRESLHY geometry observation operator that fits modelled
  direct-plus-ground-image values at the actual sensor heights. It confirms
  that the independent-energy rise/width rejection is not a hidden-state
  comparison artefact.
- Added Li et al. (2026)'s enthalpy-only established-flow balance as an
  explicit alternative to HyRAM+ total energy. It passes all seven interfaces
  and conserved marches but changes centre MAE by less than 0.04 mm, ruling
  out resolved kinetic-energy thermalisation as the rise-error mechanism.
- Added a direct near-field versus projected crosswind buoyancy-moment audit.
  All seven interfaces preserve the force sign and differ by at most 6.761%,
  rejecting an interface density-shape discontinuity as the rise-error cause.
- Corrected the Houf local-Froude width conversion from scalar to velocity
  e-folding width (`B=sqrt(2 sigma_y sigma_z)/lambda`). The repeated 7-trial
  run improves MG/VG to 1.062/1.163 and internal width ratio to 0.976, but
  centre-height MAE remains above baseline, so the path stays research-only.
- Added a trial-10 vertical-momentum budget. It finds only 0.0657 N handoff
  momentum difference but 5.695 versus 2.734 N cumulative buoyancy by 6 m,
  locating the remaining rise error in the early thermal/density evolution.
- Added a coefficient-free metastable dry-air condensation bound. Trial 10
  remains warm and weakly buoyant and fails the unchanged 2 K interface gate,
  so delayed nucleation is rejected before any seven-trial field scoring and
  pre-10D mass entrainment becomes the next isolated mechanism.
- Added a four-flux `source_flux` Gaussian-establishment bound after tracing
  the source-zone definitions in Li et al. (2026). It removes a second
  application of Zone-V momentum entrainment inside Zone IV, passes 7/7 and
  improves concentration to MG/VG/FAC2 1.012/1.162/0.976, but does not beat
  baseline centre-height and width errors and remains research-only.
- Added coefficient-free `geometry` and `surface_layer` ground-contact bounds
  to the independent-energy crosswind model. Trial 10 remains 0.469 m too high
  at 6 m and fails the frozen all-section gate, so both remain off-default and
  no seven-trial score is calculated.
- Added an explicit Li-equation-35 enthalpy transport option to the conserved
  near field. Its `source_flux` plug-to-Gaussian boundary is incompatible on
  the representative source and now fails closed before downstream scoring.
- Added direct PRESLHY thermocouple validation from the original trial-10 and
  trial-23 workbooks, including local sensor-height temperature observation.
  The 42-point audit finds the independent-energy candidate systematically
  too warm, especially in trial 23, and locates the next gap in phase/thermal
  development rather than drag.
- Corrected the D3.6 Table A3 reader for Unicode minus signs and seven-character
  thermocouple serials.  Added an explicit sustained-mean/peak source selector;
  the peak bound improves temperature and geometry but fails the frozen VG
  gate, so the mean-rate baseline remains unchanged.
- Added a synchronized raw-workbook source/temperature audit. Trial 23 shows
  the expected colder-at-higher-flow sign at all eight centreline channels,
  but the common contrast is too small to explain the cold-core residual;
  trial 10 is dominated by continued apparatus cooling. This prevents a
  transient source correction from being mistaken for a phase-model repair.
- Separated LH2 caloric temperature from mechanical tanker pressure and
  pre-registered a boiling-liquid source bound.  The temperature-only form
  fails the source applicability gate in all seven PRESLHY trials because a
  pressure-thrust-free reconstruction loses the available pressure work.  It
  remains off-default while measured nozzle pressure, quality/density and
  effective area are reconstructed.  Added a range-selective RADAR TAR
  utility so individual workbooks can be obtained without the 11.3 GB archive.

### Finite-rate condensed-air applicability is now quantified

The cryogenic-air add-on now exposes Sandia's uncorrected Ranz--Marshall
transfer number and a heat-limited minimum N2/O2 particle sublimation time.
The latter assigns all convective heat to latent heat and omits Stefan
resistance, so it is an explicit fastest-transfer bound rather than a fitted
kinetic model.  Across the eligible LH2 sources, 1 um N2 particles disappear
within 7--40 mm on this bound while 100 um particles usually survive past the
equilibrium handoff.  No finite-rate default is adopted because the field
campaigns do not constrain particle size, number density or nucleation.

## Unreleased

First working version.

### The conserved LH2 near field now has an audited crosswind handoff

A new research path carries the accepted axisymmetric near field into JETPLU
by projecting total mass, hydrogen and horizontal/vertical momentum as exact
integral constraints. Total energy, hydrogen half-width and centre temperature
are independent fail-closed screens. JETPLU now exposes its exact integral
fluxes, including analytic Gaussian kinetic-energy factors.

Two pre-registered reduced thermodynamic mappings fail unchanged. The accepted
candidate transfers the full near-field N2/O2/H2O equilibrium and component
enthalpy radial profile. On the frozen 0.08 m representative boundary it
closes native balances below `6e-16`, energy to 0.0585%, H2 width to 1.57% and
centre temperature to 0.94 K on a strict-grid confirmation.
`run_lh2_crosswind_research` assembles the
consistent local-wind/near-field/handoff path for horizontal wind-aligned
releases; rejected handoffs cannot integrate downstream. Independent
PRESLHY/Spadeadam rescoring remains outstanding.

A pre-registered Spadeadam 4/6 pilot using the existing 1 micrometre conserved
source reconstruction stops before sensor scoring: test 6 passes a
phase-manifold table-domain extension, while test 4 fails fixed width and
temperature screens and lies below the momentum-dominated velocity ratio.
The coupled runner now makes that applicability condition fail-closed.

### The Raman benchmark now follows the final journal source

The active Hecht--Panda observations now use the final 2019 article's printed
mass-decay and mass-width fits, 0.2771 and 0.07069, while preserving the 2017
conference values as versioned provenance. Model predictions are unchanged.
The recommended dry conserved model remains inside all four 25% bands, but its
mass-centre error is now -24.86% and its temperature-width error +24.56%.

The validation record also no longer silently resolves a source contradiction:
Table 1 and radial panels contain nine conditions, whereas the aggregate-fit
legends list an additional 4 bar, 45 K, 1.25 mm series. Pass counts are marked
provisional pending the authors' fit membership, reduced data and uncertainty.

### Para-hydrogen calorics and ELVHYS source identifiability are audited

The conserved axisymmetric model now accepts an explicit normal-, para- or
ortho-hydrogen component enthalpy table without changing the species mass or
density relation. The public LH2 research entry point exposes the same
off-default spin-isomer sensitivity. A pre-registered nine-case para-hydrogen
run closes every conservation threshold, but worsens the corrected Raman
thermal-centreline error from -21.38% to -34.79%; normal hydrogen remains the
recommended dry baseline.

Checksum-verified ELVHYS Test-10/Test-11 raw records were independently
reduced. A quantitative model score is deliberately withheld: the public
archive contains no hydrogen mass-flow channel, the nominal repeats have
different pressure/temperature source histories, and official documents give
conflicting 200 and 250 mm nozzle elevations. The exact missing fields and a
ready-to-send HSE request are documented.

A second pre-registered candidate combines independent thermal/species radial
profiles with phase equilibrium through a fully conservative four-flux
boundary. It reduces the Raman temperature-width error from +24.56% to -2.27%
but worsens centreline temperature to -29.14% and is retained only as a
research closure. The former unsafe phase/two-scalar combination remains
blocked unless four-flux establishment is selected.

### The accepted Raman model now has a public research entry point

Measured throat conditions can now be expanded through the same conserved
HyRAM+ boundary used in validation and passed directly to the recommended
dry-air axisymmetric model. The returned object carries source state,
trajectory, boundary residual, species drift, external-energy-corrected drift
and applicability warnings. Existing atmospheric and DEGADIS defaults are
unchanged.

Atmospheric argon phase change was added as a pre-registered explicit
sensitivity. It remains 4/4 and grid-converged but worsens the thermal-centre
error, so it is not enabled in the recommended configuration. A grey
radiative-absorption energy source and separate external-heat ledger were also
implemented. Even the perfect-black `5B` upper bound slightly worsens both
mass metrics, so its absorptivity stays zero by default.

### Humid-air frost is now included in the conserved Raman model

The axisymmetric research model can now carry explicit ambient absolute
humidity and solve local H2/N2/O2/H2O phase equilibrium, including water
freezing, vaporisation-plus-fusion latent heat and condensed-water volume.
It remains off by default. A pre-registered saturated-air upper bound closes
all boundary invariants and keeps nine-case energy drift below `1e-4`, but
overpredicts centreline warming and is rejected. A 40%-RH sensitivity puts
all four Hecht--Panda slopes inside 25%; the experiment does not report RH,
so this is a diagnosed input uncertainty rather than a fitted improvement.

The Raman comparison now follows the unequal stitched-image coverage in the
source paper (369 rather than 549 samples). Fixed turbulent Prandtl/Schmidt
spreading, temperature-dependent ideal enthalpy, dry-air-dew-point initial
heating and their two-scalar combination were also tested prospectively and
rejected without altering production defaults.

The reported 0.3 m/s honeycomb co-flow can also be represented explicitly.
Entrained co-flow mass, vector momentum and kinetic energy are conserved and
the entrainment speed is relative to the moving ambient. Its effect on all
four Raman slopes is below 0.8%, so it is retained as a boundary option but
rejected as an explanation of the thermal-centreline error.

The dry equilibrium phase model now optionally replaces constant ambient
heat capacities with low-density Helmholtz ideal-gas enthalpy tables for H2,
N2, O2 and H2O. This candidate closes boundary/species/energy residuals at
`1.50e-14`, `1.53e-5` and `3.77e-5`, and passes all four Raman slopes under
both the 549-point audit and corrected 369-point protocol. It is the first
end-to-end conservative Raman candidate to do so and is the recommended
dry-air axisymmetric research configuration; atmospheric production defaults
remain unchanged pending humidity-conditioned data.

### A conserved-energy cryogenic free-jet model is Raman-validated

An independent implementation of the published axisymmetric Gaussian
mass/momentum/species/energy balances now covers quiescent cold-gas near
fields. It includes plug-to-Gaussian flow establishment, source-momentum and
buoyancy entrainment, kinetic energy, and separate density/species profiles.
On all nine Hecht--Panda Raman releases the official-establishment variant
passes all four pre-registered centreline and half-width slopes; its aggregate
metrics agree with a separate official HyRAM 6.1 oracle to within 0.75%.
An end-to-end audit then exposed 8.9--14.8% species loss in that published
boundary. Two exactly conservative alternatives are retained as rejected
research variants because both fail the centreline-temperature criterion.

### Source total energy and pressure thrust are explicit research options

The transported condensed-air source can now include kinetic energy in the
storage-to-evaporation balance, removing an energy-creation defect at its
first station.  A separate HyRAM+ 6.0 path finds the homogeneous-equilibrium
critical throat, infers `Cd` from measured flow, conserves geometric-orifice
pressure thrust with the Yuceil--Otugen convention and solves the atmospheric
state from total enthalpy.  Unit tests close throat energy and notional-nozzle
mass, momentum and energy.

Neither option changes the corrected default.  The energy-only candidate
worsens PRESLHY variance/centre height and independent Spadeadam bias.  The
pressure-thrust candidate improves one PRESLHY centre-height subset but
worsens campaign-level bias and width, fails Spadeadam, and has no physical
ambient-pressure state for three of nine eligible PRESLHY releases.  Those
incompatibilities are reported rather than clipped.

The source can also preserve the Li et al. equation-22 evaporation-zone
distance instead of placing an already air-loaded state at the orifice.  The
12--78 mm coordinate correction changes no handoff state and slightly lowers
the predicted plume, but does not improve the two validation campaigns enough
for adoption.

Li et al.'s stationary-condensate limit is now available as a second
off-default research bound with separate gas/particle endpoint kinetic energy.
Deposited solid N2/O2 removes its enthalpy and zero axial momentum, with mass,
momentum and energy closed explicitly.  It is rejected: PRESLHY variance and
centre/width geometry and independent Spadeadam bias are all worse than the
corrected default.  Together with the fully carried bound, this prevents an
unsupported intermediate slip coefficient from being fitted to field errors.

### Condensed-air transport is available as a rejected research bound

The LH2 evaporation endpoint now permits stable solid N2 and O2 below their
triple points, closes component partial pressures and enthalpy, and feeds an
off-default source march that conserves retained/dropped mass, axial momentum
and energy.  The march retains H2 non-ideality at its component partial
pressure, transports re-evaporating condensed air, applies
Schiller--Naumann settling, and preserves its physical distance when handing
off to JetPlume.  Unit and integration tests cover pressure, energy,
continuity and the explicit source option.

Fixed 1, 10 and 100 um cases were registered before comparison.  None is
adopted: PRESLHY common-arc VG becomes 2.058, 2.150 and 8.946, centre-height
MAE becomes 0.285, 0.295 and 0.482 m, and independent Spadeadam MG/VG becomes
1.353/1.598, 1.357/1.600 and 1.460/1.694.  Fine retained condensate reproduces
the too-buoyant warm-source result; coarse dropout makes the remaining source
more H2-rich and worse.  The public corrected default is unchanged.

### The flashing source is now conserved through the first ODE station

The equivalent-source flash calculation already included the air needed to
evaporate the remaining hydrogen, but the first jet lookup placed that state
on a pure saturated-hydrogen mixing curve. For Spadeadam test 6 this changed
the hydrogen mass fraction from about 0.400 to 0.985 and discarded most of the
entrained air. The corrected path rebuilds the thermodynamic table from the
mixed flash composition and enthalpy, and a regression test now requires mass
fraction, temperature and density to survive the first lookup.

Gas-branch CoolProp evaluation is scoped to corrected jet sources. Applying it
globally damaged two independently checked NASA pool regimes. The public LH2
jet preset no longer forces permanent ground contact: with the conserved
source, test 4 remains low while low-wind test 6 detaches, matching the
reported Spadeadam regime split.

The expanded source also used to manufacture total mixture momentum while
stationary air was entrained: about 3.28 times the incoming value at 5 barg
and 6.20 times at 1 barg.  The corrected source now enforces
`u_out = Y_H2 u_in`, matching the no-slip initial-entrainment balance in the
current Sandia HyRAM jet model.  Pressure thrust from an under-expanded zone
and condensed-phase slip remain separate, explicit uncertainties.

Seven 0.35/0.53 m PRESLHY arcs lie upstream of the enlarged established source
plane. On the 62 common downstream arcs, the complete corrected model has MG
1.047, VG 1.425 and FAC2 0.84. On 23 filtered vertical fits it raises the mean
spread ratio from 0.64 to 0.731 and leaves mean / median signed centre-height
errors of 0.036 / 0.002 m. The earlier table-only source reached a spread ratio
of 0.933 but did so with non-conserved momentum, so that number is superseded.
Houf–Schefer buoyancy entrainment and DEGADIS surface-layer top entrainment are
available as pre-registered, off-by-default experiments because neither gave a
campaign-wide improvement without a competing regression.

The one-call report now cites the corrected jet evidence (MG 1.047, n=62)
instead of the historical baseline. Pool LFL distance no longer borrows the
PRESLHY jet statistic; it states that no direct pool-concentration validation
is available.

### The LH2 statistics are now computed

Reduced tables of the PRESLHY E3.5 campaign and the REDIPHEM archive make the
field results reproducible without either archive, and the suite recomputes
them rather than quoting them. `validation/nearfield.py` holds the liquid
hydrogen jet configuration, which had never been committed: every LH2 result
in this project had been produced by assembling the thermodynamics, the
boundary layer, the flash and the coefficients by hand in a session.

The configuration is checked against a full parameter dump of a reference run
and matches it to four significant figures. Getting there needed four fixes,
one of which is worth stating on its own: **`distmx` is the integration step,
not the limit.** Passing the 40 m limit there integrated the whole plume in a
single stride and produced a trajectory that was monotone, smooth, plausible
and wrong by a factor of two. The averaging time (60 s) and roughness
(0.001 m) were nowhere in the records and had to be recovered by matching
`deltay` and `ustar` against the dump; the source flow is the window mean, not
the peak.

**The published results reproduce.** Provenance exactly — 24 workbooks to 9
trials to 69 arcs. Concentration to two per cent: MG 0.722, CI [0.574, 0.903],
VG 1.44, FAC2 0.82 against a published 0.738, [0.591, 0.918], 1.41, 0.83. The
vertical spread ratio at the band it was quoted over, 0.64, exactly.

The vertical spread reproduces exactly, 0.64 -> 0.97 against a published
0.64 -> 0.94, once three unrecorded choices are pinned: the basis is 23
momentum-filtered fits rather than the 42 used for the trajectory, the
statistic is the mean of per-fit ratios rather than the ratio of medians, and
there is no distance filter.

**One published result is withdrawn.** The trajectory figures `1.07 -> 0.19 m`
come from configurations that cannot be reconstructed, and were taken over
different subsets besides. From the two dumped configurations the pair is
**0.877 -> 0.658 m** on trial 10 at 6 m, and 1.00 -> 0.85 m over the 5-7 m
band: the corrections take about a quarter off the rise, not four fifths. The
fault itself is confirmed and is the largest open defect in the LH2 path --
the measurement puts the plume below the nozzle where the model puts it a
metre above.

Every discrepancy in this exercise traced to an aggregation reported without
its subset, so every statistic now carries its distance band, release heights,
filter and n. `docs/lh2-recomputed.md` has all of it, and superseded figures
are marked in place rather than deleted.

### Review pass: the reported numbers and the computed ones

A review compared what the package prints and what its documents claim against
what the suite actually produces. The physics was unaffected; the reporting was
not.

**Two defects in `degali.lh2`.** The scope warning checked `max_distance` —
the integration limit, which defaults to 100 m — instead of the furthest
distance any reported answer relies on, so every call taking the default warned
that 100 m was out of range while every number it returned was inside it. And a
single `VALIDATED` range conflated two campaigns: the wind bounds came from the
NASA pool spills and the distance bound from their 33.8 m tower row, but both
were applied to jets whose corrected evidence spans 0.79 to 6 m. A jet asked about 30 m
was reported as inside the checked range when nothing had checked it. The range
is now per release type, and a "range" that is really one point says so.

**`degali.evidence`, new.** Every validated number, with its sample size,
grade and source, in one module. `lh2.py` printed a superseded statistic —
`MG 0.74, VG 1.22, n=53` — from a hard-coded f-string, while every document
carried the current `MG 0.738, VG 1.41, n=69`;
the neutral-buoyancy concentration appeared as four different values across the
README, the handover, a test assertion and the program's own output. Both are
now read from one place, and the report carries each claim's grade so that an
n=4 result cannot be read as firmly as an n=69 one.

**Guards, so this does not recur.** A retired figure reappearing in prose fails
the suite; the printed evidence is checked against the module; the scope-warning
regression has a test naming it; and the absolute-path guard now walks the
tests as well as the package — it had only ever walked `src/`, and the tests had
accumulated eleven machine paths.

**Test data is located from the environment.** `DEGALI_E35_ROOT`,
`DEGALI_E35_REPORT` and `DEGALI_SMEDIS_ROOT` join `REDIPHEM_ROOT`. The
reduced PRESLHY tables are derived products of this work rather than the
dataset, so they are vendored under `reference/preslhy/` with their provenance;
two far-field results now reproduce on a clone with no third-party data at all.

**Citation metadata.** `CITATION.cff` had no `authors` block, which CFF 1.2.0
requires and Zenodo rejects. Added, with `version`, and a test that blocks the
placeholder once the version is no longer a development one. The REDIPHEM entry
cited a title that is not the title of anything: REDIPHEM is the CEC DG XII
ENVIRONMENT project, and the database report is Nielsen and Ott, *A collection
of data from dense gas experiments*, Risø-R-845(EN), 1995.

**Documentation corrections.** The NASA lift-off table in `docs/liftoff.md` was
stale by up to 4 m against what the code produces, and quoted a mean error of
−0.1 m where the four values give −1.0. The suite regenerates it. The
handover's neutral-buoyancy figure of 85–90 mol % is superseded; the mixing table
gives 99.97, and the reason the two are so far apart — one linear segment at the
top of the table, in a regime where air is treated as a non-condensing ideal gas
— is now recorded alongside it. Test counts have been removed from prose
entirely: they were quoted in four documents and had drifted in all four.

### The model

All six DEGADIS 2.1 programs are ported: the source blanket (`DEG1`), the
steady and transient downwind models (`DEG2S`, `DEG2`), the cloud snapshots
and receptor histories (`DEG3`, `DEG4`), and the jet/plume model with its
touchdown bridge (`JETPLU`, `DEGBRIDG`). All five EPA test cases run end to
end in Python.

### Validation against the original

The Fortran is built from source in this repository and run as part of the
test suite. Purpose-written probes extract its internal state at full double
precision, because printed output carries five significant figures and cannot
distinguish a correct port from a close one. Every module below the
integration drivers matches to round-off.

Eight portability patches were needed to build the 2012 sources with gfortran.
One was load-bearing: `COMMON /ERROR/` mixes `REAL*8` with an `INTEGER*4`, and
`ESTRT1`'s `EQUIVALENCE` overlay made gfortran pad it differently there than
elsewhere, so the air entrainment coefficient was silently read as zero.

### Three things that trip up a port

- `ADIABAT(ifl=1)` reads its `wa` argument rather than computing it; `SZLOCAL`
  passes an uninitialised local, so `SZF` runs its lookups as though the
  mixture contained no dry air.
- `SURFAC` receives a mix of heated and adiabatic layer properties, and `PSS`
  and `SSG` differ by one letter in which temperature they pass.
- `GAMINC` returns the *unregularised* incomplete gamma, and says so in the
  comment above the assignment. That is deliberate, not a defect — but
  mapping it onto `scipy.special.gammainc`, which is regularised, leaves the
  flammable-mass derivative low by Gamma(1/(1+alpha)).

### Field evaluation

Against the Risoe REDIPHEM database, with bootstrapped confidence intervals.

At the lowest instrumented height Burro gives MG = 0.81 with a 95 % interval
of [0.63, 1.02], which straddles Hanna's bias bound rather than satisfying it.
Over the whole cloud nothing passes: DEGADIS gets the ground-level centreline
roughly right and the vertical distribution badly wrong, so the conventional
single-height statistic is carried by the one elevation where the vertical and
lateral errors cancel. The lateral spread runs one to three times too wide.

The reader refuses to guess channel type numbers, which are per series. An
earlier version fell back on a plausible list for FLADIS, which ships no
definitions, and produced a clean-looking result computed on channels that
were not concentrations. That result is withdrawn.
