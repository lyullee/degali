# FFI Test 6 near-field underprediction: physical decomposition

## Comparator

At 30 m in FFI/DNV Test 6, the reported arc maximum is 21.0 vol % H2. The frozen fast LH2 path predicts 6.03 vol % at the same measurement heights: observed/predicted = **3.48**. Test 6 is a horizontal 25.4 mm release at 0.5 m elevation, with reported mass flow 0.833 kg/s and P04 = 2.53 barg.

The FFI report describes two stacked containers, a plastic drum and an
instrument box on the test pad. The present calculation does not resolve those
structures; this decomposition therefore remains a free-ground transport
diagnostic, not an obstacle-wake validation. The site-geometry evidence and
hash record are kept in
`docs/ijhe-ffi-site-geometry-boundary-2026-10-09.md`.

This is a decomposition of what public evidence can quantify, not a calibration exercise.

## Results

| Category | Tested/publicly determined quantity | 30 m consequence | Decision |
|---|---|---|---|
| Source | The mass flow, orifice and P04 are reported, but the local public table has no independent atmospheric two-phase source state or machine-readable synchronized source history | Source contribution is not identifiable without inventing a range | Keep uncertainty explicit; do not fit flow, flash fraction or thrust |
| Wind | Low/mean/high mast = 2.3 / 2.5 / 2.7 m/s | predicted maximum = 6.03 / 6.72 / 7.36 vol %; residual = 3.48 / 3.13 / 2.85 | Material but insufficient |
| Observation operator | 15 physical 30 m sensors projected from the reported wind direction | largest steady prediction at exact sensor positions = 5.75 vol % | Fixed array geometry does not explain the residual; transient sampling remains unresolved |
| Baseline transport | free detachment | centre = 4.41 m; predicted maximum = 6.03 vol % | Retain |
| Forced ground contact | existing optional geometry switch | centre = 1.35 m; predicted maximum = 25.50 vol % | Reject: improves this point by overpredicting it and degrades the six-arc screen |
| Shape drag = 2 | deliberately strong diagnostic | centre = 4.10 m; predicted maximum = 6.79 vol % | Reject: only 13 % concentration change |
| Vertical shear | existing diagnostic | centre = 4.38 m; predicted maximum = 6.07 vol % | Reject: negligible change |
| Lift-off Richardson = 30 | literature-range diagnostic | centre = 4.41 m; predicted maximum = 6.03 vol % | Reject: negligible change |

## Conclusion

No demonstrated correction is adopted. The low-wind mismatch is real, but the public evidence does not uniquely separate an unmeasured flashing source state, wind/time variability, and unresolved transient/two-phase transport. Tested deterministic transport changes do not solve it without damaging independent balance.

The correct claim is therefore a fast quasi-steady LH2 dispersion model with a documented low-wind Test 6 limitation, not a universally accurate transient plume model. A future default change requires an independently documented boundary such as a time-resolved atmospheric source state, machine-readable synchronized wind/sensor record, or a separate low-wind flashing-release data set. The public report does contain plotted 10 Hz histories; these are now usable through the separate replay path documented in `ffi-test6-transient-receptor.md`, provided the digitization/export and its provenance are supplied explicitly.

## No-fit source/wind admissible-set diagnostic

The public values were also passed through the **centreline** operator of
`assess_observation_envelope()`. This is a source/wind diagnostic, not a claim
that a centreline value is identical to a sensor-height arc maximum.
For the reported 0.833 kg/s source and the measured 2.3--2.7 m/s mast-wind
range, the 30 m prediction remains approximately 0.0608--0.0787 mole
fraction against the 0.21 observation; none is within a factor of two. This
is not a fitted result. It is an explicit enumeration of the measured source
and wind hypotheses.

As a diagnostic only, extending the rate grid to 1.5--3.0 kg/s produces
factor-two matches, but every one of those rates is outside the frozen jet
evidence range (0.084--0.285 kg/s). The admissible set is therefore
**non-identifiable and out of the validated source envelope**, not a basis for
changing the default mass flow. The remaining residual must be treated as a
source-state/observation-operator/transport uncertainty until a synchronized
source record or low-wind release dataset is available.

The newly exposed sensor operator was then applied to the fifteen local
30-metre Test 6 sensor coordinates (without copying that table into the
package). Using the reported low mast wind, winter ambient assumptions and
the measured source, the exact-coordinate projection reaches **7.50 vol %** at
the highest sensor, versus the 21.0 vol % reported arc maximum. Thus the
sensor-height/lateral operator changes the point estimate but does not explain
the residual; it is now quantified rather than silently folded into a
centreline comparison.
