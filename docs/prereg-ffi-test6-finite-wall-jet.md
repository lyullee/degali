# Pre-registration: FFI Test 6 finite wall-jet diagnostic

Date frozen: 2026-09-21, before running the finite wall-jet field trajectory.

## Purpose and claim boundary

Test whether retaining a finite ground-attachment state can bridge the known
Test 6 30 m bracket between free detachment and permanent ground geometry.
This is a mechanism diagnostic only. Test 4 failed the earlier joint source
gate, so a Test 6-only result cannot be promoted as validation and no partial
MG/VG/FAC2 score will be quoted as a model-validation result.

## Frozen input

- Test 6 published flow 0.833 kg/s, orifice 25.4 mm and P04 2.53 barg.
- Low mast wind 2.3 m/s, release height 0.50 m.
- Ambient 277.15 K, 90% RH, neutral stability, 1 mm roughness and 10 m wind
  reference height.
- Existing pre-registered 1 micrometre transported-condensed-air source with
  total source energy, no pressure thrust and the transported source-zone
  distance retained.
- Conserved axisymmetric near field to S=3.0 m, 41 radial points, 0.01 m
  maximum step, 2e-6 relative tolerance and phase-manifold handoff.
- Independent-energy projection must pass the existing five-flux, energy,
  width, temperature and quadrature gates before field integration.

## Frozen comparison

Three downstream calculations use the same accepted boundary:

1. free interaction;
2. permanent ground geometry; and
3. `FiniteWallJetCrosswind` with critical Richardson number 30,
   response length 1 local plume depth and wall stress `rho_a u_*^2`.

No value may be changed after inspecting concentrations. Integrate past the
100 m arc with a 0.1 m maximum step and repeat the finite wall-jet calculation
at 0.05 m. Report the maximum relative five-flux balance residual and the
coarse/refined change.

## Observations and reported diagnostics

At 30, 50 and 100 m, report:

- maximum model concentration over the reported 0.1, 1.0 and 1.8 m sensor
  heights, compared with the maximum reported peak over bearing;
- centre height and vertical Gaussian sigma;
- finite wall-jet attachment fraction, geometric contact fraction and lift
  Richardson number.

The finite wall-jet mechanism is directionally useful only if it raises the
30 m sensor-height maximum above the free result while remaining below the
permanent-ground result, and subsequently releases rather than staying
permanently attached. It remains research-only regardless of outcome because
the response-depth closure and substrate thermal history lack independent
Test 6 measurements.

## Results after pre-registration

The reconstructed source reproduced the earlier endpoint: 68.7353 K,
35.1798 m/s, 0.330998 m diameter, 1.46376 kg/m3, H2 mass fraction 0.187993
and source-zone x=0.74456 m. The velocity/local-wind ratio was 22.6617.
The phase-manifold interface passed with maximum native flux residual
`3.13e-16`, energy-quadrature residual `1.37e-6`, H2 half-width residual
2.85% and centre-temperature residual 1.706 K.

The pre-registered 0.1 m flux-space step failed its strict flux-to-state
inversion (`3.42e-5` maximum residual). The existing state-space solver could
produce numbers but accumulated 0.23--0.31 relative balance errors and was
therefore not used for interpretation. The conservative recovery used 0.05 m
and 0.025 m steps without changing a physical coefficient. Its maximum
balance residual was `3.50e-7`; maximum coarse/refined concentration change
was `2.18e-6` relative and centre-height change was `1.63e-6 m`.

| path | 30 m max (%) | centre (m) | 50 m max (%) | centre (m) |
|---|---:|---:|---:|---:|
| observed peak | 21.0 | -- | 1.8 | -- |
| free | 4.079 | 5.176 | 0.762 | 10.718 |
| permanent ground geometry | 11.847 | 2.261 | 3.401 | 5.583 |
| finite wall jet | 3.737 | 5.374 | 0.606 | 11.195 |

At 30 m the finite model's geometric contact fraction was 0.141, but its lift
Richardson number was 432.3 against the fixed threshold 30. The equilibrium
and transported attachment fractions therefore fell to approximately `1e-9`;
the vertical force was effectively fully released. At 50 m there was no
geometric contact. The candidate moved in the wrong direction and fails the
declared directional criterion.

At the 3 m handoff itself the section already overlapped the ground by 25.2%,
but `Ri*=137.3`, so the registered equilibrium initialization was detached.
As a post-result bounding diagnostic, forcing the handoff attachment state to
one (without changing its one-depth response length) still relaxed to about
`1e-9` by 30 m and predicted 3.769% at 30 m and 0.609% at 50 m. Thus merely
remembering contact at the handoff does not change the decision; an unsupported
long attachment-memory length would be required.

The conservative five-state Gaussian inversion stopped before the 100 m arc;
no 100 m value is quoted. This is an additional representation-domain failure,
not permission to substitute the non-conservative state integration.

## Decision

Reject Richardson-controlled downstream attachment as the Test 6 correction.
Retain the code as an opt-in mechanism-isolation tool, but do not enable it in
the public model. The diagnostic narrows the missing state to contact history
formed before the 3 m handoff, time-dependent substrate/source thermodynamics,
or an observation-time intermittency operator. A downstream lift-off threshold
cannot recreate a wall jet after the interface already classifies it as free.

## Post-result observation-operator check

The reported 21.0% is a temporal peak. At the 90-degree sensors on the 30 m
arc, the three mean/peak/standard-deviation values are 13.5/21.0/4.3,
15.4/20.6/3.7 and 13.3/18.6/3.8%. Permanent ground geometry predicts
11.85/11.63/11.10% at the same heights: 0.38, 1.02 and 0.58 observed standard
deviations below the time means. Thus the phase-manifold grounded branch is
not a factor-of-two error against the 30 m time mean even though it remains
below the peaks.

At 50 m the same sensors report means 0.2/0.7/0.9% and peaks 1.0/1.8/1.8%.
The permanent-ground prediction 3.21/3.27/3.40% is five to ten standard
deviations above the means and exceeds every peak. The free prediction
0.659/0.691/0.762% is close to the means. This independently confirms the
required qualitative sequence: ground-influenced near 30 m, free by 50 m.
It does not identify a transition distance or response time without fitting
the same observations, so no new default parameter is inferred from it.
