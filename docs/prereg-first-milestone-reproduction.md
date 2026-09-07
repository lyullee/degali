# Stage1 reproduction and paired performance decision

2026-09-06. This fixes the new reproduction/aggregation checks BEFORE the new
end-to-end reruns, not a claim that historical performance was unseen.

## Fixed variants and populations

- BASE: corrected measured-mean-source JETPLU as instantiated by
  validation.nearfield.from_reduced(corrections=True). Its actual code uses
  STEP=.2 m and REACH=40 m; do not substitute the .02 m candidate step or
  the historical preregistration's wording for the executed baseline setting.
- CONTROL (context only): source-energy/ambient-consistent density profile.
- CANDIDATE: last completed full-field source-energy/ambient-consistent
  volumetric-enthalpy Gaussian profile, normal H2, full measured pipe source,
  collective equilibrium source bound, fixed velocity exponent1.16^2,
  polar_square_all_flux_v1, candidate max step .02 m and tolerance1e-5.
- EXCLUDED from field adoption: all newly enriched thermal/TKE/normal-stress
  modules that have only initial or short-segment/synthetic verification.

Keep the frozen seven trials10,11,12,22,23,24,25; concentration38 common arcs
and vertical17 common profiles. Use their exact recorded sensor populations.
The41-temperature set is trial10's19 and trial23's22 sensors from
gaussian_enthalpy_profile_temperature_audit_2026-09-05.json. Do not merge it
with the different42-point historical comparison. Match by trial/channel/
serial/coordinates, not list index; preserve minimum,p05,median separately.

## Provenance and reproduction

Trace complete-field aggregation hashes to all leaf shards. Check core code,
inputs and raw-workbook hashes. Record every mismatch, including the earlier
interface checkpoint's runner-script hash (core/input hashes still match).
Never rewrite that earlier evidence to match current files. Store a fresh
manifest of the current source, required tools, relevant input/outputs and
runtime versions. Do not claim historical hashes covered files not listed.

Reintegrate trials10 and23 from their unchanged sources through the latest
candidate using the existing source-ablation runner and .02 m max step.
It also recomputes the baseline. Compare the new matching concentration and
vertical rows with the frozen field outputs at relative/scale tolerance1e-8.
Compare interface and saved final-state quantities where the arrays coincide;
report coverage and every mismatch. Both must reach their required sensors,
and retain the original interface/energy balance acceptance criteria.

For BASE temperature, run the exact same hydrogen_jet configuration and
STEP/REACH as from_reduced, then evaluate at all41 fixed sensors, asserting
domain coverage instead of accepting the temperature fallback outside range.
Verify that this same trajectory reproduces the stored BASE concentration
and vertical rows on these trials before using its temperature. Replay the
candidate from its stored states, with the existing temperature/flux/arc
checks. Use actual corrected BASE, CONTROL and CANDIDATE labels throughout.

Independent aggregate checks recompute MG=exp(mean(log(O/P))),
VG=exp(mean(log(O/P)^2)), FAC2 and vertical width ratio/centre bias/MAE
from the paired rows. Max scaled difference from stored aggregates<=1e-10.
Report all eligible paired rows and rejected/missing counts. No after-the-fact
coefficient changes or sensor screening. An improvement in one metric does
not override a failure of the original joint promotion criteria.

## Scope of decision

Deliver working settings, reproduction instructions, paired error tables and
known limits. Reusing fully hashed evidence plus two representative reruns is
not seven fresh integrations or an independent holdout campaign. No Linux
oracle, slow/full release certification, general plant safety approval or
SCI readiness is inferred. Necessary new physics may be developed if a
specific demonstrated defect requires it; it must be isolated, tested and
recompared. A large unresolved new theory is not silently made a condition
of stage1 completion. Failure/no improvement is a valid stage1 result.
