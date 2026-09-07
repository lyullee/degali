# Pre-registration: Houf--Schefer buoyancy entrainment

## Why this mechanism

The corrected liquid-hydrogen path still over-predicts plume rise, especially
at low wind, and the separately modelled Spadeadam ventilation mast passes
above every detected external sensor.  The existing optional
`buoyant_entrainment` term has the right qualitative effect but uses a chosen
velocity scale and an unfounded multiplier.  It cannot be adopted on that
basis.

DEGADIS `JETPLU` has shear, cross-flow and passive-atmospheric entrainment, but
no explicit buoyancy-driven entrainment.  Houf and Schefer's hydrogen plume
model does.  HyRAM 1.0 documents it in SAND2015-10216, equations 36--42, and
the current Sandia HyRAM+ 6.1 implementation retains the closure.  This is a
specific omitted physical term with published coefficients, not another
free fit.

## Frozen implementation

The implementation will follow Sandia HyRAM+ commit
`b45abf9a6d995951311be6aad836f1874e4d420b` (version 6.1),
`src/hyram/phys/_jet.py`, lines 210--212, 278--285 and 390--399 as retrieved on
2026-09-03.

At the source plane:

```text
Fr_den = v0 / sqrt(g D0 |rho_a-rho_0| / rho_0)
a = 17.313 - 0.11665 Fr_den + 2.0771e-4 Fr_den^2,  Fr_den < 268
a = 0.97,                                             Fr_den >= 268
```

At each integration point:

```text
v_cl = max(u_c + u_a cos(theta), 0)
B = sqrt(2 sigma_y sigma_z)
Fr_L = v_cl^2 rho_cl / (g B |rho_a-rho_cl|)
E_b = (a/Fr_L) (2 pi v_cl B) max(sin(theta), 0)
```

`B = sqrt(2 sigma_y sigma_z)` is the only profile conversion: HyRAM writes
velocity as `exp(-r^2/B^2)`, while DEGALI's reported widths are Gaussian
standard deviations, `exp(-r^2/(2 sigma^2))`.  The geometric mean preserves
area for the elliptical section.

The DEGADIS shear contribution and `E_b` will be summed, then limited by the
Houf--Schefer pure-plume cap:

```text
E_shear+b <= 0.082 (2 pi B v_cl)
```

Cross-flow entrainment and passive atmospheric spread are outside the
axisymmetric Houf--Schefer closure and will be added after this cap.  The
buoyancy contribution is zero for a descending trajectory: the current
HyRAM source uses signed `sin(theta)`, but importing a negative value would
remove DEGADIS entrainment from a downward release and is not a defensible
interpretation of *buoyancy-driven entrainment* in the cross-wind model.

The option is off by default so the Fortran oracle remains exact.  No
coefficient will be fitted to Spadeadam, PRESLHY or NASA data.

## Predictions and decision rule fixed before the run

1. The term will be weak while a horizontal jet is momentum-dominated, then
   grow as the plume turns upward and its local Froude number falls.
2. It will reduce the excessive low-wind rise more than the high-wind rise.
3. Near-field concentration variance should fall without moving the
   geometric mean farther from unity or reducing FAC2.
4. The 30 m centre-height ratio for Spadeadam test 6 versus test 4 should
   move down from the corrected baseline, and test 6 concentration at 1 m
   should move toward the observed 15.4 vol %.  A change in the opposite
   direction falsifies the mechanism for the outdoor horizontal path.
5. In the mast comparison, the median vertical separation between the model
   centre and detected sensors should fall.  The existing zero false-LFL
   result for non-detect sensors must be retained.
6. A mechanism that improves only one selected case while materially
   degrading the campaign statistics will remain experimental.  It is
   adopted into `corrections=True` only if directions 3--5 all hold and the
   existing parity and regression suite remains green.
