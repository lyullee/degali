# Pre-registration: trial-10 vertical-momentum budget

Date: 2026-09-05

## Question

After correcting the Houf velocity-width mapping, concentration and width pass
their production comparison but centre-height MAE does not. Trial 10 at 6 m
is the largest internal centre residual. The remaining distinction must be
located before another physical term is changed:

1. excess vertical momentum already present at the fixed 10D handoff;
2. excess downstream integrated buoyancy relative to vertical drag; or
3. comparable momentum histories but a different velocity/residence-time
   conversion from momentum to trajectory.

## Frozen calculation

- Use PRESLHY trial 10 only, with the same corrected source, 10D near-field,
  Li enthalpy transport, corrected Houf velocity-width mapping, local-shear
  entrainment and 0.02 m maximum step as the accepted research branch.
- Compare it with the corrected JETPLU baseline built from the same reduced
  trial conditions.
- At the 10D handoff and at 1.78, 4 and 6 m, report centre height, trajectory
  angle, vertical momentum flux and section-integrated buoyancy.
- Integrate buoyancy and the signed vertical drag separately over streamline
  distance from the common 10D x-coordinate. Verify their sum against the
  computed vertical-momentum change. Do not infer force from centre height.
- JETPLU output stores ground-imaged concentration. Recover its un-imaged
  Gaussian state exactly from the stored centre height and vertical width
  before evaluating flux and force. Record the image correction magnitude.
- Change no model coefficient, state, source, handoff or comparison threshold.

## Interpretation fixed before calculation

- If most of the 6 m vertical-momentum difference is already present at 10D,
  the next target is the source-to-10D near-field trajectory/force history.
- If the handoff momenta are comparable but the force integrals diverge, the
  next target is downstream density/thermal-profile evolution.
- If vertical momenta remain comparable but heights diverge, the next target
  is velocity profile, added mass or the momentum-to-trajectory kinematics.
- The sparse sensor-fit centre at 6 m is not used to choose among these
  mechanisms; the primary diagnostic is the internal conservative budget.
- This audit may locate an error but cannot promote a new term or fitted
  coefficient.

## Result

The 10D station is 0.7966 m downstream because the conserved expanded source
diameter is 79.66 mm, not the 25.4 mm physical orifice. The two branches are
already thermodynamically different there:

| at 10D | corrected JETPLU baseline | independent-energy candidate |
|---|---:|---:|
| total mass flux, kg/s | 1.0270 | 1.2658 |
| centre temperature, K | 30.71 | 64.31 |
| centre density, kg/m3 | 2.7123 | 1.1826 |
| vertical momentum, N | -0.0631 | 0.00265 |
| section buoyancy, N/m | -0.3392 | 0.01723 |

The candidate has entrained 0.239 kg/s more total mass and has already crossed
from negatively to positively buoyant, while the baseline remains cold and
dense. Yet the absolute handoff vertical-momentum difference is only 0.0657 N.
The later momentum separation is generated primarily after the handoff:

| at 6 m | corrected JETPLU baseline | independent-energy candidate |
|---|---:|---:|
| internal centre height, m | 0.5407 | 0.7649 |
| vertical momentum, N | 2.8105 | 5.6169 |
| cumulative buoyancy from 10D, N | 2.7336 | 5.6949 |
| cumulative vertical drag, N | -0.0141 | -0.0806 |

The candidate's net integrated vertical force is 2.895 N above the baseline;
the observed 6 m vertical-momentum difference is 2.806 N. Thus initial
vertical momentum is not the material cause. The sign and magnitude point to
the thermal/density path that begins at establishment and persists downstream.

The candidate's direct conservative budget closes to `4.47e-6 N`. JETPLU's
diagnostic flux change differs from the integrated explicit buoyancy-plus-drag
source by 0.154 N (5.4% of its momentum change), even after reducing output
spacing to 0.02 m. This is reported as a legacy wind/profile-derivative
closure residual and does not explain the factor-of-two buoyancy difference.
The ground image is negligible at 10D (`1.1e-80`) and therefore cannot explain
the establishment-state contrast.

The next physical audit is now narrowed to two linked terms: the larger
axisymmetric mass entrainment before 10D and instantaneous equilibrium
air-condensation/latent-heat release. These must be separated prospectively;
the previous transported-particle sensitivities cannot be assumed to test
delayed nucleation without checking their actual phase source terms.

Machine-readable values are in
`reference/preslhy/trial10_vertical_momentum_budget_2026-09-05.json`.
