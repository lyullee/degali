# Pre-registration: smooth phase lookup trial 24 continuation pilot

Date frozen: 2026-09-10, after the independent C1 lookup confirmation and
implementation tests, before projecting or marching the trial 24 boundary.

## Question

Does the nodally exact C1 phase lookup remove the arbitrary lookup-cell corner
stagnation in trial 24 while retaining the same six boundary moments, local
physical gates, flux conservation and step convergence?

## Frozen calculation

- Use trial 24 only and the same measured-pipe source, ambient state and stored
  10D boundary moments as the retained failed extension.
- Select `phase_interpolation="c1_hermite"`; retain the existing bilinear
  path as the package default and negative-control history.
- Reproject the stored six boundary moments on the smooth phase section from
  the stored state. Require maximum moment residual `1e-8`, quadrature
  difference `1e-5`, and thermal width strictly inside 0.5--2.0.
- March only to x=0.43 m, beyond the retained x=0.41398 m cell-corner failure.
- Retain six primary fluxes, reduced buoyancy work, thermal/species ratio 1,
  internal inverse tolerance `2e-5`, endpoint tolerance `1e-8`, coarse/refined
  maximum steps 0.005/0.0025 m, minimum step 0.000078125 m, eight-success
  regrowth, and every existing local physical gate.
- Independently reconstruct terminal order-16 fluxes and require scaled
  balance `1e-5` for both marches.
- Require both terminal physical-parameter and independent-flux
  coarse/refined differences below `0.005`.
- Record rejections, growths and endpoint fallback use. Hash inputs and write
  a new immutable JSON file. Do not overwrite any failed evidence.

Passing permits a separately recorded full trial 24 extension. It does not
change defaults, fit observations, score sensors or promote the downstream
model.
