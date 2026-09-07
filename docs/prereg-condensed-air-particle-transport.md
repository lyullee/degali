# Pre-registration: transported condensed-air source zone

Date fixed: 2026-09-03, before any PRESLHY or Spadeadam run with this model.

## Question

Can a near-field source zone that conserves N2/O2 phase mass and permits
finite particle slip reduce LH2 concentration and trajectory error without
fitting an entrainment coefficient?

Li et al. (2026) demonstrates that condensed-air momentum changes the
Station-3 source, but its zero-velocity N2 dropout, fixed 77.35 K boundary and
sub-triple liquid saturation cannot be used for saturated-LH2 releases.  The
new candidate will instead use component partial pressures and stable solid
vapour pressure, transport both N2 and O2 condensed mass, and count phase
enthalpy once.

## Fixed particle-size cases

The volume-equivalent particle diameters are fixed before model comparison:

| case | diameter | role |
|---|---:|---|
| fine | 1 um | rapidly velocity-coupled aerosol |
| intermediate | 10 um | transition case |
| coarse stress bound | 100 um | appreciable settling/dropout |

A liquid-nitrogen surface experiment reported abundant cryogenic aerosol in
the 1--5 um range, but it is not an LH2 jet measurement.  Therefore 1 um is a
literature-supported order of magnitude, not a fitted PRESLHY value; 10 and
100 um are declared sensitivity bounds.  No diameter may be tuned after
viewing validation residuals.

## Fixed physics

1. At each plug-flow station, gaseous H2 is conserved and the N2/O2 gas-solid
   split satisfies `p_i <= p_sat,i(T)` using component partial pressure.
2. N2/O2 saturation below their triple points comes from the tabulated
   solid-vapour correlations in `addons.cryogenic_air`; liquid saturation is
   used only above the triple point.
3. Particle velocity relaxation uses drag with the Schiller--Naumann
   Reynolds correction.  Settling removes only the portion whose cumulative
   fall distance reaches the local lower plume boundary.
4. Retained condensate may re-evaporate as the gas warms.  Dropped condensate
   cannot re-enter downstream.
5. Gas plus retained-particle momentum is conserved across lateral
   entrainment; ambient air enters with zero streamwise momentum.
6. Gas, retained condensate and dropped condensate share one energy ledger.
   No separate latent term is added on top of condensed-phase enthalpy.
7. Handoff to JetPlume occurs at the first station with less than 1% retained
   condensed mass relative to the largest condensed mass encountered, or at
   2 m if that condition is not reached.  A 2 m failure remains a multiphase
   result and must not be passed to the single-phase dispersion ODE.

## Evaluation population

- PRESLHY E3.5 reduced table, horizontal momentum-driven trials, window-mean
  release flow, sensor-position arc maximum; comparisons use common arcs.
- PRESLHY well-constrained vertical Gaussian fits with the same momentum
  filter; comparisons use common fits.
- Independent Spadeadam horizontal outdoor tests 4 and 6, low-mast wind and
  the existing censored-arc convention.

The as-built corrected model remains the comparator.  Current references are
PRESLHY MG/VG/FAC2 `1.047/1.425/0.840` on 62 arcs, width ratio `1.033`, common
fit centre MAE about `0.160 m`, and Spadeadam `0.869/1.476/0.667`.

## Predictions fixed in advance

- The 1 um case will remain almost fully carried, approach the already-tested
  phase-safe warm-source bound, and is expected to worsen PRESLHY residual
  variance.
- The 100 um case can lose condensed mass before re-evaporation and should
  reduce the excessive buoyant rise of a fully retained warm source, but may
  over-concentrate the plume.
- The 10 um case should lie between those limits.  A non-monotone result must
  be explained by phase-boundary or handoff changes, not selected as a hidden
  best fit.

## Adoption rule

A candidate is adopted only if all of the following hold on common samples:

- PRESLHY MG remains 0.7--1.3 and moves no farther from 1.0;
- PRESLHY VG does not exceed 1.455 and FAC2 does not fall below 0.82;
- vertical centre MAE does not exceed 0.160 m and width-ratio error does not
  exceed 0.033;
- Spadeadam MG moves closer to 1.0, VG does not exceed 1.506, and FAC2 does
  not fall below 0.65;
- mass, component mass, axial momentum and energy residuals are each below
  `1e-6` relative at every source-zone step;
- the same fixed diameter passes both campaigns.  A diameter that helps one
  and harms the other remains a sensitivity case, not a correction.

## Evidence still missing

No direct size distribution for condensed N2/O2 in the PRESLHY hardware has
been found.  The decisive experimental input would be Mie/PDA measurements
within 0--0.5 m of the outlet reporting D32 or d10/d50/d90 together with
distance, nozzle diameter, reservoir/exit temperature, pressure and flow.  If
such data become available, they replace the three-case bracket prospectively;
they are not used to choose whichever existing case scored best.

## Results after pre-registration

The implementation first solves a saturated-liquid-H2 evaporation endpoint
with solid N2 and O2, then marches the retained condensed mass and entrained
air to the declared handoff.  For a 27 K, 288 K ambient example this replaces
the impossible gaseous-air endpoint by `T = 20.3689 K`, `m_air/m_H2 = 0.70008`
and essentially all N2/O2 in the solid phase.  The pressure sum closes to one
atmosphere and the relative energy residual is `1.3e-16`.

Across the nine eligible PRESLHY sources and Spadeadam tests 4 and 6, all 11
source marches reached the handoff at 0.117--0.793 m.  The largest relative
mass, axial-momentum and energy residuals over every step and size were
`7.5e-15`, `1.5e-15` and `9.5e-11`.  The final source ranges were:

| diameter | temperature (K) | H2 mass fraction | cumulative dropout / initial condensed air |
|---:|---:|---:|---:|
| 1 um | 67.999--68.725 | 0.188--0.217 | 0.00001--0.00068 |
| 10 um | 68.001--68.659 | 0.189--0.217 | 0.0012--0.0645 |
| 100 um | 64.130--68.350 | 0.199--0.360 | 0.046--1.222 |

The last ratio can exceed one because newly entrained air can condense and
drop after the initial evaporation plane.

PRESLHY concentration results below use only arcs reached by both the adopted
source and each candidate.  The changing `n` is a physical range change, not
an outcome-based exclusion.

| diameter | n | adopted MG/VG/FAC2 | candidate MG/VG/FAC2 |
|---:|---:|---:|---:|
| 1 um | 39 | 1.061 / 1.612 / 0.795 | **0.968 / 2.058 / 0.846** |
| 10 um | 39 | 1.061 / 1.612 / 0.795 | **0.982 / 2.150 / 0.846** |
| 100 um | 40 | 1.098 / 1.676 / 0.775 | **1.362 / 8.946 / 0.775** |

On the 17 common vertical fits, the adopted-source centre-height MAE is
0.187 m and its mean modelled/measured width ratio is 1.113.  Candidate
centre MAE / width ratio are 0.285/1.022, 0.295/1.025 and 0.482/1.113 for
1, 10 and 100 um.  Fine particles improve width while moving the whole plume
too high; coarse dropout worsens the centre without improving width.

The independent six-arc Spadeadam comparison is:

| diameter | MG | VG | FAC2 |
|---:|---:|---:|---:|
| adopted source | 0.869 | 1.476 | 0.667 |
| 1 um | 1.353 | 1.598 | 0.833 |
| 10 um | 1.357 | 1.600 | 0.833 |
| 100 um | 1.460 | 1.694 | 0.833 |

All candidates increase centre height at all six Spadeadam arcs.  Test 6 at
30/50/100 m moves from 3.35/6.78/13.65 m to 4.64/8.66/16.19 m for 1 um and
4.87/8.96/16.59 m for 100 um.

## Decision and implementation boundary

No fixed diameter passes the pre-registered concentration and centre-height
criteria in either campaign.  The candidate is therefore **not adopted** and
the corrected default is unchanged.  It is available only through the
explicit `condensed_air_particle_diameter` research option; its source state
and physical axial distance are preserved in the JetPlume handoff.

This is not yet a complete finite-slip aerosol model.  The current controlled
implementation uses a common axial gas/particle velocity and a well-mixed
settling loss with Schiller--Naumann terminal speed.  The 100 um stress case
has a relaxation time large enough that a two-axial-velocity closure is needed
before quantitative use.  Argon, ambient water ice and CO2 frost are also not
transported, and the sub-1%-of-peak residual condensate is folded into the
single-phase handoff.  These limitations independently bar promotion even if
one concentration statistic were favourable.

The rejection diagnoses the structural error: warming nearly all retained
condensed air to an approximately 68 K handoff makes the source too buoyant
and raises the plume too early.  Letting coarse material fall out makes the
remaining source still more H2-rich and performs worse.  A particle diameter
must therefore not be tuned to the field residuals; the next defensible step
is a measured size distribution plus two-velocity transport and persistent
water/CO2/argon solids.
