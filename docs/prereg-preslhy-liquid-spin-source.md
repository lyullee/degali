# Pre-registration: PRESLHY liquid-hydrogen spin-consistent source

Date frozen: 2026-09-05, before the first para-hydrogen source/interface or
field result.

Do not edit above the `RESULTS` line after a candidate result is known.

## Physical question

The measured-pipe LH2 reconstruction currently uses CoolProp `Hydrogen`
(normal hydrogen) from the compressed liquid state through atmospheric flash,
then also uses normal-hydrogen calorics in the conserved near field.  Stored
liquid hydrogen near 20 K is normally converted toward para-hydrogen during
liquefaction, but the released spin composition of PRESLHY E3.5 is not
reported.  NIST supplies separate reference equations for both fluids.

A property-only screen on the four measured nozzle states found that replacing
normal by para hydrogen changes upstream density and atmospheric velocity by
less than 0.21%, but increases the source-to-ambient enthalpy deficit by
12.90--12.97%.  This is large enough to test because the rejected fast-bound
candidate becomes positively buoyant too early.  No in-flight ortho--para
conversion heat is added; the spin composition is frozen over the release
residence time.

## Candidate fixed before running

1. Add an explicit `Hydrogen`/`ParaHydrogen`/`OrthoHydrogen` property selector,
   leaving `Hydrogen` as the public default.
2. For the para candidate, use `ParaHydrogen` consistently for liquid density,
   viscosity and surface tension, pressure-thrust flash, saturation quality,
   H2 evaporation endpoint, gas density, and near-/far-field component
   enthalpy.  Do not mix a para liquid source with normal gas calorics.
3. Keep all measured rates, TC3/PT2 inputs, geometry, `C_ds=15`, entrainment,
   spreading, handoff locations, numerical tolerances and field observations
   unchanged.  The three open-pipe trials use the same para property selector
   on their existing source topology; no missing nozzle pressure is invented.
4. Fit no spin fraction or field coefficient.  The candidate is an explicit
   pure-para end member because the experiment did not report composition.

## Staged decision rule

1. Unit tests must recover the normal-hydrogen default exactly, close all
   source balances below `1e-8`, and demonstrate reference-state invariance.
2. All seven conservative interfaces must pass their existing flux, H2-width
   and centre-temperature limits before field scoring.
3. On the same 38 common concentration arcs and 17 vertical fits as the
   normal-hydrogen fast bound, para hydrogen may be promoted only if MG moves
   closer to one, VG does not increase, FAC2 does not decrease, width ratio
   moves closer to one, and centre-height MAE does not increase.
4. Regardless of score, the default cannot change until the released
   ortho/para composition or an adequate storage/liquefaction history is
   documented.  A failed candidate remains an off-default uncertainty bound.

Sources: Leachman et al., *Fundamental Equations of State for Parahydrogen,
Normal Hydrogen, and Orthohydrogen*, DOI `10.1063/1.3160306`; NBS Monograph
168, section 4; CoolProp `ParaHydrogen` fluid documentation.

---

## RESULTS

Implemented as an explicit selector with normal hydrogen unchanged as the
default.  The four measured nozzle sources close all balances below
`2.1e-16`.  Relative to normal hydrogen, para hydrogen changes upstream
density by -0.181% to -0.206%, atmospheric velocity by +0.039% to +0.083%,
and the correlated droplet diameter by +0.776% to +0.865%.  Trials 12 and 24
gain only 0.00215 and 0.00234 absolute flash quality.  The collective
all-vapour air/H2 ratio decreases 0.806--0.847% and its formation distance
decreases 1.255--1.347%.

The para candidate passes all seven fixed interfaces.  Maximum five-flux
residual is `7.97e-14`, maximum H2-width mismatch is 2.211%, and maximum
centre-temperature mismatch is 0.710 K.

On the same 38 arcs, the normal fast bound has MG/VG/FAC2
`0.82075 / 1.19741 / 0.92105`; para hydrogen gives
`0.82149 / 1.20316 / 0.89474`.  Across the same 17 vertical fits, para
hydrogen materially improves centre-height MAE from 0.07004 to 0.05453 m,
nearly the 0.05420 m legacy baseline, but moves width ratio away from one,
from 0.98014 to 0.94841.  Every downstream balance remains below `4.51e-8`.

The candidate therefore fails the frozen VG, FAC2 and width conditions and is
not promoted.  The result is still physically informative: spin-dependent
calorics can account for about 15.5 mm of the fast-bound centre-height error,
whereas spin-dependent liquid density, flash and atomisation are only percent-
level effects.  Released ortho/para composition is now an explicit source
uncertainty, not an implicit normal-hydrogen assumption or a field-fitted heat
term.

Machine-readable results:
`reference/preslhy/measured_pipe_droplet_source_para_2026-09-05.json`,
`reference/preslhy/measured_pipe_droplet_equilibrium_para_interface_2026-09-05.json`
and
`reference/preslhy/measured_pipe_droplet_equilibrium_para_field_2026-09-05.json`.
