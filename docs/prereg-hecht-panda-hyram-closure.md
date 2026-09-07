# Pre-registration: HyRAM momentum and energy closures in JETPLU

**Frozen before the first run of either candidate.**  Do not edit above the
`RESULTS` line after a candidate result is known.

## Why this is not coefficient tuning

The first Raman comparison rejected all existing DEGALI variants.  The
closest, unmodified JETPLU, overpredicted the H2 half-width slope by 77% while
overpredicting the centreline mass-decay slope by only 31%.  Hecht and Panda
identify this exact combination as evidence that the velocity/scalar width
ratio and entrainment closure cannot simply retain their warm-jet values.

Sandia's independently published HyRAM+ model supplies both missing pieces:

- `lambda = 1.16`, the density/scalar width relative to velocity width.  In
  JETPLU's exponent convention this is `Sc = lambda^2 = 1.3456`.
- `beta_A = 0.28`, giving source-momentum entrainment
  `E_mom = beta_A sqrt(A_exp rho_exp v_exp^2 / rho_ambient)`.

These are the current official defaults in Sandia's public model, not values
estimated from the Raman curves.  HyRAM+ also conserves the cross-sectional
flux of `rho v (h + v^2/2 - h_ambient)` (technical manual equation 82), while
JETPLU has no energy balance and obtains temperature from concentration alone.

## Frozen candidates

1. `hyram_momentum`: replace only JETPLU's local perimeter/shear entrainment
   with the source-momentum expression above and use `Sc=1.3456`.  No
   Ricou--Spalding density multiplier is applied because source density is
   already present in `E_mom`.
2. `hyram_momentum_energy`: add a conserved total-energy flux to candidate 1.
   The energy term must include kinetic energy and the enthalpy of entrained
   ambient air.  It may use numerical radial quadrature, as HyRAM+ does, but
   may not introduce a fitted heat-transfer coefficient.

Both switches remain off by default, preserving the Fortran oracle exactly.

## Predictions and decision rule

1. `hyram_momentum` should reduce both mass-decay and mass-width slope errors
   from the legacy values, without changing source mass, momentum or total
   energy.
2. Because it slows excessive mixing, it should also reduce both temperature
   slope errors, but may make the centreline too cold without kinetic-energy
   dissipation.
3. `hyram_momentum_energy` should leave the two mass metrics nearly unchanged
   (within 5%) and reduce the centreline-temperature error relative to
   candidate 1.  A larger mass change indicates unintended thermodynamic-
   buoyancy feedback in this momentum-dominated range and must be reported.
4. Eligibility requires all four printed slopes within 25%, both equivalent
   radial coefficients inside the printed case ranges, exact mass/momentum/
   energy source residuals below `1e-8`, and no regression in the existing
   PRESLHY comparison or Fortran parity tests.
5. If candidate 1 fails the mass metrics, candidate 2 is not eligible: energy
   conservation cannot repair the wrong mixing field.  No coefficient will be
   adjusted to force a pass.

---

## RESULTS

The regression interpretation correction recorded in
`prereg-hecht-panda-raman.md` was applied before these candidate results were
calculated: centreline slopes have fitted intercepts; half-width slopes are
forced through zero.  All nine Table-1 releases and all 61 axial samples per
release are included.

| model | mass decay | mass width | T decay | T width | `A_Y` | `A_T` | slopes within 25% |
|---|---:|---:|---:|---:|---:|---:|---:|
| observed | 0.26260 | 0.06503 | 0.02826 | 0.06188 | 49 (33--64) | 42 (21--49) | - |
| JETPLU + HyRAM momentum | 0.20827 | 0.06533 | 0.02748 | 0.08208 | 52.0 | 31.8 | 3/4 |
| conserved-energy Gaussian jet | 0.25456 | 0.05774 | 0.02264 | 0.07098 | 64.2 | 41.2 | **4/4** |
| official HyRAM 6.1 oracle | 0.25266 | 0.05782 | 0.02255 | 0.07103 | 64.0 | 41.1 | **4/4** |

Relative errors for the conserved-energy implementation are **-3.1%**,
**-11.2%**, **-19.9%**, and **+14.7%**, respectively.  Its cross-case relative
RMSE is 0.0385 for centreline mass decay and 0.0109 for centreline temperature
decay.  The median `A_Y=64.2` is 0.3% above the largest individually printed
case coefficient (64), which is reported as a strict miss rather than rounded
down; `A_T=41.2` lies inside the printed 21--49 range.

The separately installed official HyRAM 6.1 package was used only as an
oracle.  The DEGALI implementation evaluates the five published integral
fluxes and their numerical Jacobian independently; no GPL implementation code
is included.  For the first release at 40--100 mm, centreline mass fraction,
temperature and width agree with the oracle to about 0.5%.  Across the frozen
metrics, the independently implemented slopes differ from the oracle by
0.75%, 0.14%, 0.41%, and 0.07%.

### Decision

Candidate 1 is retained as a useful diagnostic but rejected as the final
free-jet closure because temperature width misses the frozen band by 32.6%.
Candidate 2 passes all four primary slopes and independently reproduces the
official model.  It is therefore accepted as the new **quiescent
axisymmetric near-field** model.  It is not silently substituted into the
cross-wind JETPLU path: a conservative handoff between these different profile
definitions is required first.  The original path and all Fortran parity
defaults remain unchanged.

### Post-result establishment audit — acceptance narrowed

The statement immediately above is retained as the decision made from the
registered Raman metrics, but a subsequent end-to-end flux audit found a
boundary defect before the ODE: the published algebraic plug-to-Gaussian
conversion loses 8.9--14.8% of hydrogen species flux and 20.4--29.0% of
relative total-energy flux across these nine cases. The ODE then conserves
the reduced flux accurately, which is why an ODE-only residual test did not
expose it.

Two coefficient-free conservative replacements were frozen and tested. Both
close their boundary fluxes below `5e-12`, but each misses the centreline
temperature slope by about 31% and therefore scores only 3/4. Consequently,
the published-establishment model is retained as an **official-model
reproduction and Raman benchmark**, not accepted as an end-to-end
source-conservative production model. The two failed repairs are recorded in
`prereg-entraining-gaussian-establishment.md` and
`prereg-scalar-constrained-establishment.md`.

### Final-journal source correction

The final-journal mass fits change the conserved-energy errors to -8.13% and
-18.33%; its temperature errors remain -19.88% and +14.71%. The 4/4 result
therefore remains, but only as a provisional aggregate-fit check because the
published fit legends contain an untabulated tenth series. See
`hecht-panda-journal-benchmark-correction.md`.
