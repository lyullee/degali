# What the data has been asked, and what is left in it

An interim record. Every dataset in hand, what it was used for, what it
refused, and what remains untried.

> **Gaussian-width correction:** the `sigma_z is 1.4x narrow` finding below
> compared a fitted e-folding width with JETPLU's standard deviation. On a
> common definition the corrected model/measurement ratio is 1.033. See
> `gaussian-width-convention.md`.

## Used, and settled

**Superseded — recomputed from the reduction; the figures below were never computed by anything.** Recomputing gives MG 0.855, CI [0.672, 1.064], VG 1.30, FAC2 0.85 as shipped and MG 1.070 corrected, over 66 arcs, and the interval now includes 1. See `docs/lh2-recomputed.md`.
| dataset | what it settled |
|---|---|
| DEGADIS 2.1 Fortran | the port reproduces it to 1e-12 across all six programs |
| EPA test cases (5) | end-to-end reproduction, byte for byte |
| REDIPHEM Burro (8) | MG 0.811 at 1 m, CI [0.63, 1.02]; the vertical profile fails above it |
| REDIPHEM Desert Tortoise (4) | MG 1.839; and EPA's own DEGADIS score explained |
| SMEDIS FLADIS + DT (76 sensors) | the vertical defect on an independent substance and release type |
| PRESLHY E3.5 near field (69 arcs) | MG 0.738, CI [0.591, 0.918], VG 1.41, FAC2 0.83 |
| PRESLHY E3.5 far field (18 arc fits) | the bias reverses; the trajectory is the cause |
| PRESLHY vertical fits (42) | the plume does not rise; `sigma_z` is 1.4× narrow |
| NASA Witcofski (4 spills) | the buoyancy regime, 4 of 4 |
| EPA-450/4-90-018 | an external anchor, and the height artefact behind it |

## Asked and refused, with the reason

| dataset | asked for | refused because |
|---|---|---|
| PRESLHY far field, point-to-point | concentrations at fixed stands | wind direction recorded to 22.5° against stands 10–12° apart — solved instead by fitting the plume position out |
| PRESLHY near field | transient behaviour | travel time under a second against 0.3 s sampling |
| PRESLHY Dräger sheet | far-field time series | 90 channels at 1 Hz, but values run −275 to 1970 under a `%` header; units unresolved |
| REDIPHEM FLADIS | anything | ships no `CHANDEF.DAT`; recovered from SMEDIS instead |
| REDIPHEM Eagle | N₂O₄ dispersion | the substance dissociates to NO₂, so its molecular weight is temperature-dependent |
| REDIPHEM Thorney Island | puff dispersion | the van Ulden momentum balance has no solution at H/D ≈ 1; the original Fortran stops in the same place, and refining the grid a hundredfold does not move it |
| SMEDIS Thorney Island | 263 sensors | peaks at 2060 under a `mean_C(%)` header |
| SMEDIS Prairie Grass | 200 sensors | peaks at 235 under the same header |
| SMEDIS EMU | 224 sensors | sensor heights read 244553 m |
| PRESLHY E3.4 pool trials | LH₂ pool dispersion | concentrations at 35, 45 and 55 cm only, no downwind distance |
| PRESLHY E3.1a | dispersion | near-orifice discharge, 1.2 GB in nested archives |

## Tried this session and recorded as unusable

**BA-Propane, 273 sensors.** Eight of the ten SMEDIS trials have downwind
sensors and a REDIPHEM counterpart. Paired at sensor level they give
**MG 0.168, VG 1581, FAC2 0.11** — the model reading six times high.

The REDIPHEM reduction of the same trials gives MG 3.76: the model reading
nearly four times *low*. **Two reductions of the same experiments, differing
by a factor of twenty and in opposite directions.**

Lathen was already recorded as a series not to use — small low-elevation
propane jets where release momentum dominates the near field, dispersed from
an equivalent source that discards exactly that. This confirms it and adds
something: when a comparison is this badly posed, the *reduction* decides the
answer more than the model does. That is EPA's 1991 conclusion arriving from a
third direction.

**EEC17 is genuinely empty.** REDIPHEM places every one of its 36
concentration channels at or upwind of the release, and SMEDIS independently
places all 34 of its sensors between −83.8 and −23.8 m. Two sources agree, so
the "no usable measurements" message on that trial is true.

## Still untried

**The Dräger time series.** 90 channels at 1 Hz across the far field, where
travel time is 3 to 7 seconds — long enough to resolve against 1 Hz sampling,
unlike the near field. Blocked only on the units. If the header can be
resolved this is the one place the transient path could be tested on hydrogen.

**Witcofski Table 2.** Time-resolved grab-bottle concentrations for Test 6 at
tower 5, three heights. The OCR scrambles which value belongs to which height
and time; the original table would give a single trial's time development,
an axis this work has none of.

**Chirivella and Witcofski (1986)**, AIChE Symposium Series 82, 120–140. The
follow-up data paper to the 1984 one, not in hand. It is the only route to
enlarging the lift-off sample beyond four.

## What no dataset here can settle

The trajectory fault. Five candidates were registered and tested; what remains
is slip between droplets and gas, which needs a second velocity field. No
amount of the data above distinguishes a model that lacks that from one that
has it wrong, because every campaign in hand is at one density ratio.
