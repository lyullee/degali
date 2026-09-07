# The liquid hydrogen statistics, recomputed

> **Superseded for the current corrected model (2026-09-03).** This document
> preserves the earlier reconstruction and its historical configurations.
> Current source-state invariants and campaign statistics are in
> [`lh2-model-improvements-2026-09-03.md`](lh2-model-improvements-2026-09-03.md).
> Width ratios in this historical reconstruction used unlike Gaussian width
> definitions. See [`gaussian-width-convention.md`](gaussian-width-convention.md)
> for the correction from 0.731 to 1.033.

Every figure here is produced by the test suite from
`reference/preslhy/e35_reduced.json`. Nothing in this document is typed in.

That distinction is the reason the document exists. An audit found that no
statistical claim in this package was computed by anything: the parity results
were enforced to twelve decimals and the grade-B results, which are what the
papers rest on, were numbers a working session had produced and a person had
transcribed. Two of them had already drifted.

Run it:

```
pytest -k "near_field_statistic_recomputed or vertical_sub_models or exclusion_chain"
```

---

## The provenance reproduces exactly

| | count |
|---|---|
| workbooks in the dataset (DOI 10.35097/1481) | 24 |
| excluded — not a horizontal release | −7 |
| carried forward | 17 |
| excluded — wind-steered, exit speed below ten times the wind | −8 |
| **trials in the reported statistic** | **9** |
| arcs, as shipped | 69 |
| arcs, both configurations reach | 66 |

Every number matches what `docs/lh2-results.md` published. The exclusion
criteria are computed from the flash and the measured flow, never from
agreement, and `test_the_exclusion_chain_reproduces_the_published_provenance`
holds them.

The drop from 69 to 66 is new and worth stating: the expanded-source
correction starts the plume further downstream, so the 0.35 m arc falls
outside the trajectory for trials 10, 22 and 25. Comparing 69 against 66 makes
part of the difference a change of sample, so both are restricted to the 66
they share.

---

## The configuration is now validated, not guessed

The jet setup had never been committed. It was rebuilt, and the session
holding the archives supplied a full parameter dump for one run — trial 10,
corrections on — to check the rebuild against.

**It matches to four significant figures.** `ustar`, `rhoa`, `rhoe`, `deltay`,
`deltaz`, `betaz`, `gammaz`, the starting vector, and the trajectory at every
tabulated station: `z` 1.156 against 1.158 at 6.04 m, `sigma_z` 0.6917 against
0.6909.

Four things had to be corrected, and each was worth a factor:

| | was | is |
|---|---|---|
| `distmx` | 40.0 — passed as if it were the limit | **0.2, the integration step**; `smax=40` is the limit |
| `rhoe` with the expanded source | expanded density | **saturated-vapour density, unchanged** |
| averaging time | 18.4 s | **60 s** |
| roughness | 0.01 m | **0.001 m** |
| source flow | `flow_peak_gs` | **`flow_mean_gs`**, the window mean |

The first is the one worth remembering. `distmx` is the **step**, not the
limit, and passing 40 there integrated the whole plume in one stride. The
result was monotone, smooth, and wrong by a factor of two — nothing in the
output said it was under-resolved. Neither the averaging time nor the
roughness had been written down anywhere; both were recovered by matching
`deltay` and `ustar` against the dump.

---

## Concentration

Arc maximum against the model centreline, on the 66 arcs both configurations
reach, from nine trials.

Italic rows are the superseded published figures, kept beside the ones that
replaced them.

| config (italics superseded) | MG | 95 % CI | VG | FAC2 |
|---|---|---|---|---|
| as shipped | **0.722** | [0.574, 0.903] | 1.44 | 0.82 |
| *published* | *0.738* | *[0.591, 0.918]* | *1.41* | *0.83* |
| corrections on | **1.113** | [0.873, 1.395] | 1.30 | 0.89 |
| *published* | *1.228* | — | — | — |

**The baseline reproduces to two per cent on all four statistics.** The
earlier conclusion that the interval had moved to include 1 was an artefact
of the wrong source flow and is withdrawn; the model does read high, and the
interval does exclude unity, exactly as published.

The corrected figure is 1.113 against a superseded 1.228. That pair
was computed on 69 arcs against 66 — the expanded source starts the plume
downstream of the 0.35 m arc for three trials — so part of that gap is the
sample rather than the model. On the arcs both reach, the correction takes MG
from 0.722 to 1.113 and both pass all three of Hanna's criteria.

### The convention still moves the answer

| as shipped | MG | VG |
|---|---|---|
| against the model centreline | 0.722 | 1.44 |
| against the model at each sensor | 1.132 | 19.1 |

The session with the archives confirmed the published statistic used the
centreline, **and agreed the at-sensor convention is the defensible one**: by
6 m the model centreline sits above the topmost sensor, so it is a
concentration the array could not have measured.

That change is not free. VG goes from 1.44 to 19.1, and the reason is
mechanical rather than modelling: once the plume centre climbs past the array
the arc maximum is a tail value, and tails are where small differences in
`sigma_z` become large differences in concentration. The centreline
convention does not flatter the model, it smooths the comparison.

**This is a decision for the paper, not a calculation.** Whichever is chosen,
it has to be stated once and applied to the baseline as well.

---

## The vertical sub-models

The published basis is 42 fits: every well-constrained fit on a horizontal
trial, no momentum filter, where well-constrained is `r² > 0.85`, at least
four points, `sigma_z < 3`. The momentum filter applies to the concentration
statistics only, because the fits are measurements and the filter is a
statement about model applicability. That was confirmed, and it had never been
written down as a criterion.

### The spread correction holds

| | ratio to measured |
|---|---|
| as shipped | 0.49 |
| corrections on | **0.77** |
| *published* | *0.64 → 0.94* |

Neither endpoint is matched but the shape is: the section is too narrow as
shipped, and the corrections remove **55 per cent of the deficit** against a
published 83 per cent. The corrected run is the one with a reference dump and
it reproduces exactly, so the discrepancy sits on the baseline start, which
has none.

### Both configurations are validated against dumps

| | as shipped | corrected |
|---|---|---|
| `rho_exit` | 5.34689 (orifice, two-phase after flash) | 2.54663 (expanded) |
| `concentration` | 1.0 | 0.44141 |
| `diajet` | 0.02540 m (the orifice) | 0.04602 m (equivalent bore) |
| exit velocity | 69.9 m/s | — |
| `alfa1` | 0.057 | 0.0833 (historical dump; superseded default) |
| density-scaled entrainment | off | on |

The baseline source plane is worth stating on its own, because it is what the
expanded-source correction exists to fix: **the density is the two-phase value
after the flash while the concentration and the area are the orifice's.** One
plane for one quantity and a different plane for the other two. The correction
makes the three consistent, and that is a coherence argument rather than a
fitted improvement.

Both trajectories reproduce to four significant figures, with identical row
counts and integration ranges.

The reference dumps predate the primary-source correction. The production
configuration now defaults to Papanicolaou and List's measured plume value,
`alfa1 = 0.0875`; only the dump-matching calls retain 0.0833 explicitly.

### The trajectory: a quarter, not four fifths

`1.07 → 0.19 m` is **withdrawn at both ends**. Neither figure comes from a
configuration that can now be reconstructed, and they were taken over different
subsets besides — 5–7 m over all release heights against 3 m and beyond for the
0.5 m releases alone, four fits against seven.

| | as shipped | corrected | improvement |
|---|---|---|---|
| rise at 6 m, trial 10 | **0.877 m** | **0.658 m** | 25 % |
| rise over the 5–7 m band | 1.00 m | 0.85 m | 15 % |
| measured, 5–7 m | — | — | **−0.12 m** |

**The corrections take about a quarter off the rise, not four fifths.**

The fault is untouched and is the largest open defect in the liquid hydrogen
path. At 5–7 m the measurement puts the plume *below* the nozzle and the model
puts it a metre above; a quarter off a metre is not a fix.

### The vertical spread reproduces exactly

| | ratio |
|---|---|
| as shipped | **0.64** — published 0.64 |
| corrections on | **0.97** — published 0.94 |

Three things had to be identified to get this, none of them written down:

**The basis is 23 fits, not 42.** The momentum filter *is* applied to the
spread ratio and is *not* applied to the trajectory fits. Two different
populations in the same table, distinguished nowhere.

**The statistic is the mean of per-fit ratios**, not the ratio of medians. On
this data those differ by 0.11 — 0.64 against 0.53 — which is larger than the
correction being measured on some subsets.

**The distance band is all of it**, no filter.

The corrections remove 92 per cent of the spread deficit. That result is
unchanged from the publication and is now computed rather than transcribed.

### The pattern behind all of it

Every discrepancy in this exercise came from an aggregation reported without
its subset, and there is more to it than reproducibility.

The trajectory figures were three medians over three different bands quoted as
one quantity. The spread ratio was over a filtered 23 while the trajectory was
over an unfiltered 42, in adjacent rows of the same table. The statistic form —
mean of ratios against ratio of medians — was worth more than the effect being
measured. The source flow was a window mean in one place and a peak in another.

None of it was wrong. All of it was unlabelled, and the labels are what make a
number comparable to the next one.

**This is not the first time in this project.** "The section is twice too
narrow" was two release heights averaged together. "The entrainment shortfall
grows with distance" was a median table over a seven-trial subset. Both fell
in pre-registered regression tests. With those, three cases make an argument
rather than a caveat: **a statistic reported without its aggregation is not
merely unreproducible — it misleads its own author**, because the next
comparison is made against it as though it were the same quantity.

**Every statistic here states its distance band, release heights, filter,
statistic form and n.**

## Burro cannot be recomputed yet

`MG 0.811` at the lowest instrumented height is now the only grade-B claim in
the package with nothing behind it, and `reference/rediphem_reduced.json` does
not close it.

The reduction is sound as far as it goes — 56 trials, 840 arc maxima, the
series counts and the three Burro heights all check out, and
`test_the_rediphem_reduction_cannot_yet_recompute_burro` verifies that much.
But arc maxima are *measurements*, and a statistic needs a prediction beside
them. Building one calls `trialcase.to_case`, which reads thirteen fields from
the trial specification:

```
pool diameter          release rate            site average windspeed
surface roughness      Monin-Obukov length     ambient temperature
ambient pressure       relative humidity       exit temperature
nozzle diameter        initial concentration   release duration
reference height for wind
```

The reduction carries five keys per trial and none is one of these.

**A second reduction with those thirteen fields per trial would close the last
one.** It would be a few kilobytes.

---

## What this changes for the papers

**Paper A (evaluation methodology) is unaffected.** Its two results are
grade A: EPA's published DEGADIS score is an evaluation-height artefact, and
arc-maximum evaluation hides a structural defect that appears independently on
SMEDIS ammonia. Neither touches these numbers.

**Paper B (LH₂) needs two edits, and gains a result.**

1. **The near-field concentration stands.** MG 0.722 [0.574, 0.903] against a
   published 0.738 [0.591, 0.918] — reproduced to two per cent, now computed
   by the suite rather than transcribed. The corrected figure becomes 1.113
   on the arcs both configurations reach.
2. **The trajectory correction is withdrawn.** The three corrections move the
   rise by three centimetres, not by 0.88 m. The published pair compares
   against a baseline that cannot be reconstructed.
3. **The vertical spread correction stands** in shape if not in magnitude:
   55 per cent of the deficit removed, against a published 83 per cent.

Losing the trajectory fix costs the paper its tidiest result and gives it a
better one. The corrections demonstrably fix the section and demonstrably do
not fix the trajectory, so the two faults are **independent** — which is the
argument for local validation the paper is making anyway, now carried by a
negative result rather than confused by a claimed fix.

The framing holds unchanged: *local validation separates faults that a
concentration comparison mixes.* One of them turned out not to be fixed. That
is the method working.

**A methodological point worth making in the paper.** Every discrepancy in
this exercise traced to an aggregation that was reported without its subset —
a rise quoted as one number that was three medians over three different
distance bands and release-height groups, a spread ratio over 42 fits
described in the same breath as a concentration statistic over 9 trials. None
of it was wrong; all of it was unlabelled. Every statistic in this document
now carries its distance band, release heights, filter and n, and that is a
cheaper discipline than the alternative.

---

## Two corrections from an independent check

Both were checked back against the reduction and both stand.

### The far-arc variance is one trial, not a distance effect

Reported here as "beyond 3 m the scatter goes out of control, VG 50.7". That
number is real and the reading of it was wrong.

| 3-7 m | n | MG | VG | FAC2 |
|---|---|---|---|---|
| all trials | 15 | 2.171 | 50.72 | 0.67 |
| without trial 20 | 13 | 1.107 | 1.35 | 0.77 |

Every other trial in the campaign sits between VG 1.02 and 1.59. Trial 20 sits
at 9724. It is the lowest wind, 0.57 m/s, and the model lifts its plume to
4.9 m at 6 m downwind against a top sensor at 1.75 m.

**The applicability filter should have excluded it and cannot.** It keeps a
trial when the exit velocity is at least ten times the wind, and trial 20 has
the lowest exit velocity in the campaign -- 9.8 m/s, a factor of forty below
trial 11 -- yet passes at 17.3 because the wind is 0.57 m/s. Dividing by a
small wind makes a weak jet look momentum-dominated.

That is a defect in the criterion visible without knowing the outcome, but the
outcome is what prompted the look, so the replacement is pre-registered in
`docs/prereg-applicability-filter.md` rather than applied. Both numbers are
pinned in the suite.

### "Everything except the trajectory reproduces" was too strong

`sigma_y` and far-field dilution are not separately established.

- **`sigma_y`.** The far-field arc fits give 0.8 to 2.2 m once meander is
  removed, against a modelled 1.14 m. That is inside a range, not an
  agreement, and the meander angle used to remove it comes from the same fits.
- **Dilution.** The band where the model is exact, 0 to 1 m, is where the
  centreline falls only from 100 to about 88 vol %: almost no dilution has
  happened. Dilution is exercised from 1 to 6 m, where MG runs 1.1 to 1.5.

The defensible statement is: **the source state and the near-field `sigma_z`
are established, the trajectory is the principal fault, and `sigma_y` and
far-field dilution have not been separated from each other.**

### And the trajectory does not explain everything

The overlap between failing arcs and arcs where the modelled centre is above
the array is real but partial: 0 of 27 in the first metre, 2 of 26 from 1 to
3 m, 5 of 13 beyond 3 m. The two extreme arcs are both centre-above-array, so
it explains the tail; eight of thirteen arcs beyond 3 m are inside the array
and it does not explain those.
