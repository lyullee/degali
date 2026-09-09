# Pre-registration: six-trial extension of the flux-space thermal march

Date frozen: 2026-09-09, before calculating any of the six extension trials.

## Question

Does the Trial 10 direct-flux result extend to the other six already selected
PRESLHY boundary cases without loss of inverse solvability, local physical
validity, conservation or step convergence?

This is an applicability and numerical-conservation extension. It does not
fit observations, change a source, select a new diffusivity, or promote the
research closure.

## Frozen population and targets

Retain trials 11, 12, 22, 23, 24 and 25. For each trial, use the first
well-constrained vertical concentration profile strictly downstream of its
stored 10D boundary:

| Trial | Boundary x (m) | Target x (m) |
|---:|---:|---:|
| 11 | 0.7259876095 | 0.79 |
| 12 | 0.3585012082 | 1.78 |
| 22 | 0.7960846902 | 1.78 |
| 23 | 0.6106914874 | 0.79 |
| 24 | 0.4066452935 | 1.78 |
| 25 | 0.7959056170 | 1.78 |

The targets must be re-derived from the hashed reduced PRESLHY record at run
time and must equal this table. Trial subsets may run in separate immutable
files, but no six-trial conclusion is allowed until all rows are combined.

## Frozen physics and numerics

- Replay each accepted independent-enthalpy-width boundary.
- Reconstruct the same measured-pipe/HEM source and ambient state.
- Retain `thermal_species_ratio=1`, not fitted to these trials.
- Retain `mechanical_work=reduced_buoyancy_work`.
- Retain the free planar geometry; do not add yaw or ground interaction.
- Use the already verified six-flux inverse and classical RK4 implementation.
- Coarse arc step: 0.005 m; refined arc step: 0.0025 m.
- Maximum arc length: 2.0 m; accepted-step ceiling: 4000.
- Six-flux inverse tolerance: `1e-8`; target-x tolerance: `1e-8` m.
- Stop a trial at its first inverse, phase, weak-budget, inward-flow,
  positive-diffusion or curvature failure. Do not clip or extrapolate.

## Acceptance and decisions

For every selected trial:

- boundary moment replay error must be at most `1e-9`;
- both marches must reach the frozen target;
- independently reconstructed order-16 terminal balance must be at most
  `1e-5` for both steps;
- maximum coarse/refined scaled difference must be at most `0.005` for both
  physical section parameters and the six order-16 fluxes;
- hashes must be unchanged before and after each calculation.

Retain all failures. No sensor/profile score, coefficient fit, default change,
or package release is permitted in this stage. If all six pass, a later,
separately pre-registered observation comparison may use their calculated
states.
