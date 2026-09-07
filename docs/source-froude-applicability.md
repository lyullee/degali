# Source Froude number: what can and cannot be claimed

Date: 2026-09-03

## Decision

Do **not** replace the current descriptive exit-speed/wind screen with the
proposed `Ri < 1e-4` filter. Keep the existing reported sample unchanged until
a cryogenic source plane and an applicability threshold are pre-specified.

This is not because a source Froude number is unhelpful. It is because the
candidate Richardson number and the primary-source Froude number are different
quantities, and the source experiment does not match a flashing LH2 jet in
crosswind.

## Primary evidence

Schefer, Houf and Williams, *Investigation of small-scale unintended releases
of hydrogen: Buoyancy effects*, International Journal of Hydrogen Energy 33
(2008) 4702–4712, DOI
[`10.1016/j.ijhydene.2008.05.091`](https://doi.org/10.1016/j.ijhydene.2008.05.091),
tested round gaseous-hydrogen jets at `Fr = 58, 99, 152, 268`. The `Fr = 268`
case followed momentum-dominated centreline decay; buoyancy effects increased
as Fr fell. Sandia identifies the work as SAND2007-7287J in its
[official publication record](https://www.sandia.gov/research/publications/details/investigation-of-small-scale-unintended-releases-of-hydrogen-buoyancy-effec-2007-11-01/).

The corresponding densimetric source number is

    Fr = u / sqrt(g d |rho_a - rho_j| / rho_j)

where `u`, `d` and `rho_j` must describe the same jet plane. The candidate in
`prereg-applicability-filter.md` was

    Ri_candidate = g d (rho_a - rho_j) / (rho_a u^2)

It uses `rho_a` rather than `rho_j` in the density ratio. Consequently its
threshold is not the reciprocal of the Schefer number and `Ri < 1e-4` is not a
primary-source criterion.

HyRAM 1.0 uses the Schefer/Houf formulation to blend momentum and buoyancy
entrainment, with a coefficient branch at source `Fr = 268`; see Sandia report
[SAND2015-10216](https://h2tools.org/sites/default/files/SAND2015-10216-HyRAM-1.0-Technical-Reference-Manual.pdf).
That branch changes an entrainment correlation. It is not a rule for deleting
all lower-Fr measurements from validation.

## Diagnostic calculation on the PRESLHY sample

Using the published trial conditions, CoolProp hydrogen properties, the
absolute density contrast and internally consistent values at each source
plane gives:

| trial | Fr at orifice | Fr at equivalent plane |
|---:|---:|---:|
| 10 | 160 | 162 |
| 11 | 1,400 | 1,418 |
| 12 | 2,236 | 2,265 |
| 20 | 29.8 | 24.4 |
| 21 | 70.3 | 57.7 |
| 22 | 137 | 139 |
| 23 | 665 | 673 |
| 24 | 1,614 | 1,634 |
| 25 | 138 | 140 |

This is a useful regime label: trials 20 and 21 are clearly the least
momentum-dominated sources. It is not yet a defensible exclusion rule. Trial
20 changes from 29.8 to 24.4 solely by moving from the orifice to the
equivalent plane, while the experiment being modelled is a flashing cryogenic
release whose mixture warms from initially dense to strongly buoyant.

## Why `Fr = 268` is advisory here

1. Schefer's source was a 1.91 mm gaseous jet; PRESLHY uses a 25.4 mm flashing
   liquid-hydrogen source.
2. Schefer's jet was vertical in a quiescent ambient; PRESLHY is horizontal in
   crosswind.
3. `Fr = 268` identified momentum-like source behaviour. Buoyancy becomes more
   important downstream as local velocity falls, so a source label cannot
   guarantee far-field momentum dominance.
4. The LH2 density changes sign relative to air as the cryogenic mixture warms.
   A single source density does not describe that evolution.

## Required next test

Pre-register all of the following before recomputing performance:

- source plane: orifice or fully expanded equivalent plane;
- density and concentration at that same plane;
- signed versus absolute density contrast;
- whether Fr is a stratification label, a model-branch switch, or an exclusion
  criterion;
- a threshold supported by flashing-LH2/crosswind evidence rather than chosen
  after seeing trial 20's error.

Until then, report the current population and show trial 20 sensitivity as a
secondary analysis. Do not present its removal as an accuracy improvement.
