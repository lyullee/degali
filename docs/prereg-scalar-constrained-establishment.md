# Pre-registration: scalar-constrained conservative establishment

**Frozen before the first result from this candidate.** Do not edit above the
`RESULTS` line after a candidate result is known.

## Candidate

The fixed-velocity boundary has no conservative solution, while integrating
constant source entrainment over the entire empirical development length adds
too much air. Retain instead the published plug-to-Gaussian scalar-peak
relationship,

`Y_cl / Y_plug = (lambda^2 + 1) / (2 lambda^2)`,

and solve centreline velocity, Gaussian width and centreline density so that
the established profile exactly preserves plug hydrogen-species, axial-
momentum and total-energy flux. Total mass is then a prediction, not a fitted
constraint. Development length, location, angle, radial profiles and every
entrainment constant remain frozen.

## Decision rule

1. The three normalized boundary residuals must be below `1e-8` for all nine
   sources and every solved state must be positive.
2. All four Raman slopes must be within 25%.
3. Median radial coefficients must be within 33--64 for mass and 21--49 for
   temperature.
4. The downstream species and energy tests and all historical defaults must
   remain unchanged.
5. No profile factor, virtual origin or regression definition may be changed
   after the run.

---

## RESULTS

Rejected. All nine sources have positive states and boundary species,
momentum and energy residuals below `5e-12`, but the four-slope result is only
3/4:

| mass decay | mass width | T decay | T width | `A_Y` | `A_T` | passes |
|---:|---:|---:|---:|---:|---:|---:|
| 0.21622 | 0.05962 | 0.01942 | 0.07456 | 60.3 | 37.5 | 3/4 |

Relative errors are -17.7%, -8.3%, **-31.3%**, and +20.5%. This independently
repeats the previous candidate's thermal-centreline failure while avoiding
its assumed total-air integral. Therefore the failure is not repaired by a
different conservative Gaussian boundary constraint. No virtual origin,
profile factor or acceptance band is changed.

Both conservative boundaries leave the same omitted physical process in the
40--80 K mixing layer: N2/O2 condensation and its latent heat. That mechanism
is investigated separately rather than using a heat-transfer multiplier.
