# Pre-registration: near-field/crosswind buoyancy-moment audit

Date: 2026-09-05

## Question

The five-flux interface conserves mass, H2, vector momentum and the selected
energy flux. The first downstream vertical-momentum derivative is instead
driven by a non-conserved section moment,

```text
F_b = integral g (rho_amb - rho) dA.
```

The near-field density profile follows local N2/O2/H2O phase equilibrium and
need not be Gaussian. The crosswind state replaces it with one Gaussian
density-excess profile. Existing interface screens constrain H2 half-width and
centre temperature but do not constrain `F_b`. An exact five-flux projection
can therefore begin with the wrong vertical force.

## Frozen calculation

- Use the seven fixed momentum-dominated PRESLHY releases and their fixed 10D
  handoff stations.
- Use the Li enthalpy transport selected for research continuation in
  `prereg-established-flow-energy-partition.md`; repeat total energy only if
  the result is unexpectedly sensitive to that choice.
- Integrate the near-field force directly on its existing radial grid and use
  the exact finite-Gaussian area factor for the crosswind force.
- Do not change the state projection, a coefficient, handoff distance or any
  acceptance threshold.
- Report signed forces, signs, absolute differences and
  `projected/near-field` ratios when the denominator is not near zero.

## Interpretation fixed before calculation

- A sign reversal or a greater than 10% magnitude difference in a trial is a
  material interface defect and motivates a buoyancy-preserving profile
  representation before any more downstream coefficient work.
- If all seven signed forces agree within 10%, the interface density shape is
  rejected as the cause of the excessive rise.
- This is a diagnostic gate, not permission to add a sixth algebraic
  constraint to five states. If it fails, the next model must supply a real
  density/thermal profile degree of freedom or transport closure; it must not
  rescale gravity or fit a force multiplier to PRESLHY.

## Result

All seven interfaces preserve the sign of the section-integrated buoyancy and
remain inside the frozen 10% magnitude gate:

| trial | near-field `F_b`, N/m | projected `F_b`, N/m | projected/near-field | difference |
|---:|---:|---:|---:|---:|
| 10 | 0.0184834 | 0.0172337 | 0.932387 | 6.761% |
| 11 | 0.00975757 | 0.0101377 | 1.038962 | 3.896% |
| 12 | 0.00277225 | 0.00285855 | 1.031129 | 3.113% |
| 22 | 0.0123868 | 0.0123574 | 0.997627 | 0.237% |
| 23 | 0.00597760 | 0.00634850 | 1.062048 | 6.205% |
| 24 | 0.00203882 | 0.00215952 | 1.059201 | 5.920% |
| 25 | 0.0159395 | 0.0151558 | 0.950834 | 4.917% |

The worst discrepancy is trial 10 at 6.761%, and no sign reversal occurs.
The interface-density-shape hypothesis is therefore rejected under the frozen
rule. The five-flux projection does not begin its downstream march with a
materially wrong vertical force; adding a fitted gravity multiplier or an
extra algebraic boundary constraint is not justified. The next diagnosis must
follow the downstream accumulation of buoyancy, particularly trial 10, where
the longest-range centre-height error dominates the aggregate result.

Machine-readable values are recorded in
`reference/preslhy/interface_buoyancy_moment_2026-09-05.json`.
