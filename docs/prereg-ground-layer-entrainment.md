# Pre-registration: ground-layer top entrainment

> Historical result. These figures use the table-consistent source before the
> later momentum-conservation correction. See `prereg-source-momentum.md` and
> `lh2-model-improvements-2026-09-03.md` for the current default and metrics.

## Defect found before the comparison

When `JetPlume.ground_effect` detects that the lower edge of a horizontal
plume intersects the ground, it removes the buried part of the curved
entrainment perimeter and scales upward buoyancy by the fraction of the
plume that has cleared the ground.  It does **not** add the top-surface
entrainment used by DEGADIS for a cloud in the atmospheric surface layer.

This makes the present option one-sided: ground contact reduces rise, but it
also reduces dilution.  After correcting the flashing-source thermodynamic
state, this is visible directly in Spadeadam tests 4 and 6.  The 30 m plume
height is brought close to the ground, while concentrations at 50 and 100 m
remain much too high.

DEGADIS already contains the missing closure.  In neutral conditions its
surface-layer entrainment velocity is proportional to the friction velocity,
and stable density stratification suppresses it through `PHIF`:

```text
w_e = kappa u_* (1 + alpha_wind) / Phi(Ri*)
```

The same neutral form is already used by `addons.UnifiedClosure`; the
Richardson correction is the existing `core.entrainment.phi` implementation.
No new empirical coefficient is needed.

## Frozen implementation

The correction is active only while `ground_effect` is active and the
elliptical cross-section intersects the ground.

For ellipse semi-axes `a = delta sigma_y`, `b = delta sigma_z` and the cut at
`z_cut = z_j / cos(theta)`, the horizontal contact chord is

```text
L_ground = 2 a sqrt(1 - (z_cut/b)^2)
```

The local power-law exponent `alpha_wind` is the exponent already calculated
by `JetPlume._wind`; it is not re-fitted to the trials.  The layer Richardson
number uses the above-ground depth and the Gaussian-to-layer conversion that
is already represented by `delta`:

```text
H = z_j + delta sigma_z cos(theta)
Ri* = g max(rho_bulk-rho_air, 0) H
      / (rho_air u_*^2 delta)
w_e = VKC u_* (1 + alpha_wind) / PHIF(Ri*, 0, 3)
E_ground = w_e L_ground
```

Only positive density excess is passed to `PHIF`.  That matches the existing
DEGADIS observer path, which applies stable suppression to a dense layer but
does not reinterpret a buoyant plume as a negative-stability surface layer.
`E_ground` is added to the total-mass entrainment balance.  It vanishes
continuously as the plume clears the ground.  The existing free curved
perimeter, cross-flow and passive-dispersion terms are unchanged.

The switch is off by default in `JetCoefficients` and is enabled explicitly
by the corrected liquid-hydrogen entry point.  The Fortran oracle therefore
remains unchanged.

## Predictions and decision rule fixed before the run

1. The term is exactly zero for an airborne plume and positive for a grounded
   plume; the free-plume oracle must remain bit-for-bit unchanged.
2. It must increase total entrainment and lower the excessive 50 and 100 m
   concentrations in Spadeadam tests 4 and 6 without increasing their centre
   heights materially (more than 0.25 m at 30 m).
3. Across the six horizontal Spadeadam arcs, geometric mean bias and geometric
   variance must both move toward one, or one may remain within 2% if FAC2
   increases.  A change that helps only one selected arc is rejected.
4. The corrected flashing-source endpoint invariant must remain exact.
5. PRESLHY common-arc MG, VG and FAC2 will be reported.  Since those near-field
   plumes have little ground contact, an absolute FAC2 loss greater than 0.03
   or an MG change greater than 5% vetoes adoption.
6. Houf--Schefer entrainment remains off, and no coefficient is fitted after
   seeing the results.

## Results

The missing term was real but not the dominant error.  With the corrected
flashing source and both horizontal Spadeadam trials run at the mean of the
two reported mast winds, adding the term changed the six-arc comparison as
follows:

| ground-layer top term | MG | VG | FAC2 |
|---|---:|---:|---:|
| off | 0.323 | 13.46 | 0.50 |
| on | 0.333 | 12.75 | 0.50 |

The direction is favourable but small: 100 m concentration falls by 5.1% in
test 4 and 1.4% in test 6.  PRESLHY is unchanged to the shown precision
(`MG 1.3849 -> 1.3850`, `VG 1.4905 -> 1.4908`, `FAC2 0.739 -> 0.739`).
The term is therefore retained as an explicit, off-by-default experimental
closure; it is not used to claim a material accuracy improvement.

The comparison exposed a larger issue.  Once the flashing-source state is
made self-consistent, the old `ground_effect` no longer improves the model:
it prevents the low-wind plume from detaching.  At the lower measured wind,
with no ground-layer add-on:

| trajectory treatment | MG | VG | FAC2 |
|---|---:|---:|---:|
| forced ground contact | 0.310 | 12.81 | 0.50 |
| free to detach | **0.795** | **1.71** | **0.67** |

The qualitative result also becomes correct: the high-wind test 4 stays
close to the ground, while the low-wind test 6 rises from 2.46 m at 30 m to
12.47 m at 100 m.  Mack et al.'s independent EFFECTS analysis and the DNV
interpretation report the same grounded/lifted split.  Consequently the
public liquid-hydrogen assessment no longer forces `ground_effect=True`.
