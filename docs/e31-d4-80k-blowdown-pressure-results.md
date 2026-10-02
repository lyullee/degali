# E3.1 4-mm, nominal-80-K cryo-compressed source-pressure audit

Date: 2026-09-17

## Scope and result

The public PRESLHY E3.1 part-A 4-mm, nominal-80-K set contains seven
cryo-compressed hydrogen releases: approximately 5, 10, 23, 50, 100, 150 and
200 bar. All seven have now been evaluated using the same source-start rule
and the same declared, non-fitted blowdown envelope.

This completes **source-boundary validation**, not atmospheric-dispersion
validation. It establishes which part of the pressure history the present
well-mixed, adiabatic source can represent and where an unmeasured thermal or
two-phase-withdrawal closure would be required.

## Fixed comparison

- Dataset: Lyons, Coldrick and Atkinson, PRESLHY E3.1 part A,
  [DOI 10.35097/1187](https://doi.org/10.35097/1187), CC BY 4.0.
- The raw workbooks stay local and excluded from the package. Only their
  cryptographic hashes are retained in ignored local audit records.
- Source start is the first recorded PNoz value 2 bar above its pre-relay
  median. This avoids treating the electrical relay as instantaneous flow.
- At each physical source start, measured vessel pressure initializes the
  source model. It is an observed initial condition, not a fitted offset.
- Every case uses the full predeclared `T0 = 77, 80, 84 K`, `Cd = 0.6, 0.7,
  0.8`, and two-phase vapour-withdrawal/HEM end-member set. No best member is
  selected. Tank/pipe wall heat transfer remains zero.

The detailed pre-registration is
[prereg-e31-d4-80k-blowdown-pressure.md](prereg-e31-d4-80k-blowdown-pressure.md).

## Timing evidence

The pressure-derived PNoz response delay is not constant across the set:

| nominal initial pressure (bar) | 5 | 10 | 23 | 50 | 100 | 150 | 200 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| PNoz 2-bar response after relay (s) | 0.1010 | 0.0905 | 0.0810 | 0.0655 | 0.0780 | 0.0675 | 0.0660 |

The 0.0355-s range is large enough that relay time must not be silently used
as a common mass-flow start. DEGALI’s E3.1 audit reader now uses a direct
numeric XML path for these unusually large public pressure sheets while
requiring exact row alignment across time, vessel pressure, PNoz and relay.
It reduces one 662,000-sample pressure read from an impractical object-heavy
scan to about 13 s in the present environment; small or structurally unusual
workbooks retain the ordinary safe reader.

## Declared pressure-envelope screen

Each entry reports whether the observed vessel pressure is inside the complete
18-member envelope at a fixed time after measured PNoz response. It is not an
accuracy score and does not identify an individual `Cd`, temperature or
withdrawal mechanism.

| pressure (bar) | 0.5 s | 1.0 s | 2.0 s |
|---:|---|---|---|
| 5 | outside | outside | no member remained above ambient |
| 10 | inside | outside | no member remained above ambient |
| 23 | inside | inside | outside |
| 50 | inside | inside | inside |
| 100 | inside | inside | outside |
| 150 | inside | inside | outside |
| 200 | inside | inside | inside |

At 50 and 200 bar, the pressure record is consistent with the declared
adiabatic/end-member source envelope through 2 s. At low pressure the recorded
vessel empties or depressurizes faster than every declared model member; at
23, 100 and 150 bar this divergence appears by 2 s. The result is informative
because it is **not** repaired by changing `Cd` within the independently
published 0.6--0.8 range.

For orientation only, the fixed central sensitivity (`80 K`, `Cd=0.7`) gives
at 2 s 2.21/4.46/5.86/6.69 bar for the 50/100/150/200-bar vapour-withdrawal
cases, whereas the synchronized records are 1.61/3.21/4.38/5.86 bar. The
central values are not fitted and are not promoted over their envelopes.

## Physical decision

1. The source model now has a public, seven-condition pressure-history screen
   rather than a single 200-bar diagnostic.
2. The initial single-phase tank path, explicit vapour withdrawal and HEM are
   retained as research source bounds. They are not default dispersion sources.
3. The late-time misses cannot justify a hidden heat-transfer coefficient:
   the pressure record alone does not identify tank-wall material, internal
   convection, pipe geometry, pipe-wall temperature, phase stratification or
   actual outlet phase.
4. The default DEGALI atmospheric dispersion solver remains unchanged. The
   pressure boundary is available only when the user supplies an explicit
   cryo-compressed source model and its assumptions.

The next source improvement would require a geometrically specified tank and
pipe thermal boundary, or independent time-resolved mass-flow/outlet-phase
data. Neither is inferred from the pressure trace here.
