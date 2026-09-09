# Preregistered analytic enthalpy-section width split

Date frozen: 2026-09-10 (Asia/Seoul), before implementing the analytic split.

## Observation

At the retained C1 Trial 24 failure, the reduced M2 Jacobian predicts an area
column of `23.8892219176`, while central differences below `3e-6` give
`-0.7715352408`.  The discrepancy is not thermodynamic.  The inherited JETPLU
width splitter solves the two algebraic geometry constraints with a Brent root
absolute tolerance of `1e-4`; the section widths therefore remain on a
numerical plateau during the much smaller changes required by the six-flux
inverse.

The existing constraints are

`sigma_y * sigma_z = A`

and

`sigma_y**2 - sigma_z**2 = sigma_ya**2 - sigma_za**2 = D`.

Writing `t = sigma_y**2` gives the positive quadratic root

`t = (D + sqrt(D**2 + 4*A**2))/2`.

This is an exact solution of the established geometry, not a new fitted model.

## Frozen implementation rule

1. Add a width-split hook to `IndependentEnergyCrosswind`; its default remains
   the legacy `jetplume._split` behavior.
2. Override only `GaussianEnthalpyCrosswind` with the positive analytic root.
   Evaluate the root in a cancellation-safe form when `D < 0`.
3. Preserve the existing `spread_floor` rule exactly: when enabled,
   `sigma_y = max(sigma_y, sigma_ya)` and `sigma_z = A/sigma_y`.
4. Use the hook consistently in both wind averaging and section geometry.
5. Do not change the core/original JETPLU splitter, phase data, transport
   equations, source terms, coefficients, inverse tolerances, or defaults
   outside the independent-enthalpy add-on.

## Frozen checks

1. Without the spread floor, the analytic widths satisfy both algebraic
   constraints to relative `2e-14` across positive randomized inputs spanning
   six decades.
2. The legacy hook remains bit-for-bit the existing splitter for ordinary
   `IndependentEnergyCrosswind` instances.
3. At the retained Trial 24 failure, the reduced M2 Jacobian agrees with central
   differences at steps `1e-5`, `3e-6`, and `1e-6` to relative `2e-5` or
   absolute `2e-6` in every column.
4. Existing thermal and reference-parity tests pass without tolerance changes.
5. The already preregistered C1 Trial 24 pilot is rerun first with unchanged
   physics and acceptance settings, retaining the prior failed evidence and
   writing a new immutable result.

Failure of the derivative check rejects the analytic split as the Trial 24
resolution even if the march happens to proceed.
