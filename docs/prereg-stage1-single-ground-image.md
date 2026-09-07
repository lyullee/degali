# Stage1 reporting-interface correction: apply the ground image once

Registered after the fixed reproduction/step diagnostic and BEFORE evaluating
corrected field scores. This is a demonstrated output-adapter defect, not a
new fitted physical closure. Existing frozen core, adapter and evidence stay
unchanged. The new opt-in adapter is named ReportedJetTrajectory.

## Derivation from executed code

`JetPlume.run.outp` stores column2 as C_report=C_bare*(1+exp(-2*z_c^2/sigma_z^2)).
It also stores the corresponding imaged density, temperature and mole fraction.
`Trajectory` currently treats column2 as C_bare and multiplies by the direct
Gaussian plus its ground image again. At a stored node and y=0,z=z_c this
reapplies the same factor. It is not a new entrainment or thermodynamic law.

Recover C_bare at EACH stored node by division by the known factor, then
interpolate the five primitive columns and apply the direct+image expression
once at the receptor. Do not divide by a factor evaluated only after
interpolating an imaged concentration. Source, ODE and geometry are unchanged.

## Scope and checks

- Only positive-elevation, non-touchdown JetPlume reported rows. Explicitly
  reject touchdown/invalid rows: the legacy touchdown interpolation has a
  separate convention requiring a separate audit. Do not apply to the
  independent-energy candidate, whose stored state has a different meaning.
- Manufacture reported rows from a known bare Gaussian. Verify centreline
  recovery, lateral/vertical receptors, symmetry, zero normal gradient at the
  ground, interpolation order, source-row immutability and rejection of bad
  inputs. No field-derived coefficient or tolerance is introduced.
- Reuse freshly reproduced BASE all7 output rows with their verified hashes;
  rebuild exactly the same thermodynamic tables. Check uncorrected scoring
  first, then compare fixed38 arcs/17 profiles/41 temperature receptors.
- Preserve BASE/CONTROL/CANDIDATE original scores. Label the additional result
  BASE_SINGLE_IMAGE; no concentration sensor, temperature statistic or trial
  may be removed to improve it. Geometry must remain unchanged.
- Report coarse/fine .2/.1m changes on10/23 with the corrected adapter too.
- Correctness is judged by the exact reporting identity, not by whether error
  scores improve. If worse, retain the worse scores. Do not claim this repairs
  the independent-energy closure or makes the model generally validated.

This bounded correction is included because the user explicitly authorized
necessary model corrections/new equations. It does not reopen Stage2 TKE work.
