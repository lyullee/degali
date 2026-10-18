# TI021 early-arrival / low-concentration bias: stratification audit plan

Status: predeclared physical diagnosis; no coefficient has been selected from
TI021 concentration, dose, arrival or departure observations.

## Evidence that motivates the audit

The completed 2.5 m, 600 s baseline passes conservation and nonlinear-solver
checks but predicts the six arc peaks and doses low by 51.9% and 79.5% MAPE,
respectively, and predicts their peak times 60.7 s early in absolute mean. Its
60-sensor dose ordering is nevertheless strong (Pearson 0.847; Spearman
0.914). A strict `1e-4` versus `1e-6` peak-window replay passes every numerical
tolerance gate, while the 5-to-2.5 m maximum-field relative-L1 change remains
28.9%. Numerical grid convergence is therefore still unresolved, but loose
nonlinear tolerance is not a viable explanation for the large physical bias.

The current TI021 configuration applies a neutral surface-layer mixing field
through local eddy-viscosity scalar and thermal diffusion. That field depends
on height and the declared neutral `u*`/roughness inputs, but it is not reduced
inside a locally stable, dense R12 layer. This can spread the low-concentration
front too quickly while diluting the near-ground peak and dose.

This is a physically plausible hypothesis, not a conclusion:

- Briggs et al.'s dense-plume experiment reports that local gradient
  Richardson number best describes changes in mean flow and turbulence, with
  marked dense-gas effects beyond a value near 0.25
  ([DOI 10.1016/S0304-3894(98)00162-9](https://doi.org/10.1016/S0304-3894%2898%2900162-9)).
- Kaltenbach, Gerz and Schumann find that stable stratification suppresses
  vertical overturning and that strong-stratification vertical diffusivity is
  much smaller than cross-stream diffusivity
  ([DOI 10.1017/S0022112094002831](https://doi.org/10.1017/S0022112094002831)).
- The EPA DEGADIS 2.1 guide separates negative-buoyancy, stably stratified and
  passive regimes and states that dense-gas stratification damps turbulence
  and reduces vertical mixing
  ([EPA report](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=2000J5GU.TXT)).

## Observation-blind implementation order

1. Finish the revised 5/2.5/1.25 m impact-window baseline. Do not modify the
   running fine member or use TI021 observations to choose a closure.
2. Add a diagnostic-only local gradient Richardson field,
   `Ri_g = max(0, -(g/rho) drho/dz / (|dU_h/dz|^2 + epsilon_s^2))`, with the
   density and velocity gradients evaluated by the same finite-volume geometry
   operators used by transport. Record the regularisation scale explicitly.
3. Emit volume-weighted `Ri_g` distributions and neutral eddy-diffusivity
   exposure separately for the fixed front, roof, side and wake regions and
   for the 60 receptor sampling cells. This stage must not alter fluxes.
4. Only after the diagnostic is reproducible, implement an opt-in anisotropic
   candidate that damps the vertical turbulent diffusivity while retaining a
   separately reported cross-stream component. Molecular diffusion is never
   damped. A hard `Ri_g=0.25` cutoff is not an acceptable hidden default.
5. Freeze at least a neutral baseline and a literature-bounded smooth-damping
   candidate before opening any evaluation observations. Keep all density,
   species and sensible-enthalpy fluxes conservative and include the candidate
   choice in restart/configuration hashes.
6. Require numerical time/grid checks first, then report TI021 peak, dose and
   arrival changes. A candidate can become a dense-gas component model only if
   it improves independent field structure without breaking conservation,
   boundedness, obstacle-region or 4 vol% diagnostics. TI021 remains R12
   evidence, not hydrogen validation.

## Falsification conditions

The hypothesis is rejected or retained as non-default if any of the following
occurs:

- the revised fine-grid baseline removes the early/low bias without the new
  operator;
- diagnostic `Ri_g` is small in the dense layer and receptor-relevant path;
- damping improves peak amplitude only by delaying or trapping mass in a way
  that worsens dose, arrival, wake transport or global mass closure;
- the effect is not robust to the predeclared time/grid refinement path; or
- independent hydrogen/helium obstacle evidence contradicts the resulting
  front/roof/side/wake redistribution.

Until those gates are evaluated, the production 3-D default remains the
current neutral surface-layer closure and quantitative promotion remains
withheld.

## Implementation status

Step 2 and the region-summary primitive are now implemented in
`dense_gas_stratification_3d.py`. Vertical gradients use arithmetic
fluid-fluid face reconstruction and the transport operator's zero-normal-
gradient convention at solid/domain faces. The shear regularisation is a
required explicit input; solid cells are zeroed and excluded from summaries.
The companion summary reports `Ri_g` quantiles, the fraction at or above 0.25,
and the *undamped* neutral scalar eddy-diffusivity exposure separately. It has
no damping function, changes no flux, and carries no observation-derived
coefficient. It will be applied to the checkpointed baseline only after the
running 1.25 m member is complete, preserving item 1 above.

The post-run driver is
`tools/audit_thorney_ti021_stratification_diagnostic.py`. It refuses an
incomplete checkpoint journal, verifies every segment SHA-256, requires the
previously frozen concentration-blind mixing-policy artifact, and samples the
same fixed obstacle regions plus all coordinate-only inner receptor cells.
The shear regularisation is a required command-line value rather than a hidden
default. The output remains non-promoting even when every diagnostic gate
passes.

Before any diagnostic output was opened, the regularisation envelope was
frozen as `epsilon_s = {0.5, 1, 2} u*/H_wall`, using only the accepted mixing
policy's `u*=0.205 m/s` and the report-priority `H_wall=5 m`. The required
values are therefore `0.0205`, `0.041`, and `0.082 s^-1`; the middle member is
reported first but all three must be retained and selection after execution is
prohibited. The hash-pinned plan is
`outputs/thorney-ti021-stratification-diagnostic-plan-2026-10-11-v1/`.
The audit driver requires that plan as an input and rejects either a mixing-
policy hash mismatch or a regularisation value outside the three frozen
members.
`tools/audit_thorney_ti021_stratification_envelope.py` executes all three
members, hash-pins every member record and figure, and reports time-resolved
minimum--maximum envelopes for each fixed region and receptor cell. It does
not expose a best-member selection.
