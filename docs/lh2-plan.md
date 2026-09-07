# Liquid hydrogen: what can be attempted

Written before any LH₂ modelling, so that the limits are on record in advance
rather than discovered as excuses afterwards. The LNG and ammonia work is in
`docs/field-validation.md`; the data is described in `docs/data-inventory.md`.

## Why hydrogen is not just another dense gas

DEGADIS is a *dense* gas model. Its whole structure — the gravity-slumping
source blanket, the Richardson-number entrainment suppression, the transition
from a dense phase to a passive Gaussian one — assumes the cloud starts
heavier than air and gets lighter as it dilutes, monotonically.

Liquid hydrogen does not do that. At 20 K the vapour is about 1.3 kg/m³
against ambient air at 1.2, so a cold hydrogen cloud *is* initially dense. But
hydrogen's molecular weight is 2, so as it warms and dilutes it becomes
violently buoyant — a 4 % hydrogen-air mixture at ambient temperature is
lighter than air. The cloud passes through neutral buoyancy and then lifts
off.

**DEGADIS has no lift-off physics.** Its Gaussian phase disperses a passive
tracer at ground level; nothing in it can make a plume rise. The presence of
`aeatliftoffmodelling.pdf`, `aeatliftoffreview.pdf`, the Mack buoyant-plume-
rise preprint and the HGSYSTEM documentation in the project folder says the
people who worked on this reached the same conclusion.

So the honest framing is not "does DEGADIS predict LH₂ dispersion" — it
cannot, past the point of neutral buoyancy. It is **"how far does the dense
phase carry, and where does DEGADIS stop being applicable"**, which is a
question worth answering because it bounds where a dense-gas model may
legitimately be used for hydrogen.

## What the data can carry

From `docs/data-inventory.md`, the HSL E3.5 trials are the only LH₂ dataset
here with dispersion measurements. Three properties of it govern everything
below.

**The sensors saturate at 4 % H₂, which is the lower flammable limit.** 57 of
686 readings are at the ceiling and the maximum recorded value is exactly
4.00. Those are censored observations. Treating them as measurements biases
every statistic toward the model, and it is the direction that flatters a
model that reads high.

**The far field ends at 14 m.** Burro's arcs run to 800 m. For a 285 g/s
release DEGADIS's own secondary source would be a substantial fraction of
that fetch, so there may be little or no region where the downwind model is
even the thing being tested.

**Stand coordinates are in a figure, not a table.** Table A4 gives stand
number and height but not `x` and `y`, and two layouts were used depending on
wind direction. Without those, 534 detections have no downwind distance
attached.

Per-trial signal, from the reduced CSVs:

| | trials |
|---|---|
| with a usable flow rate | 21 of 24 |
| with signal at 8 or more stands | 11 |
| with at least one saturated sensor | 11 |
| with essentially no signal | 1 (trial 14) |

## Plan

**Step 1 — recover the stand geometry.** Nothing downwind can be done without
it. In order of preference: find the coordinates in one of the other reports
(`rr985/986/987`, the E3.5 report's figures, `xxiiipaper65`); or reconstruct
them from the raw `.xlsx` if the sheets carry positions; or, failing both,
digitise Figure A4 and record it as an estimate with an explicit uncertainty.
If none works, the far-field comparison is not possible and that should be
stated rather than worked around.

**Step 2 — establish where hydrogen stops being dense.** A pure calculation,
no measurements needed: along the adiabatic mixing line for LH₂ vapour into
ambient air, find the mole fraction and temperature at which the mixture
density crosses ambient. `Thermo.build_adiabatic_table` already does this for
any substance, and the answer sets the concentration below which DEGADIS is
structurally inapplicable. If that crossover sits above 4 %, the model is
inapplicable everywhere the flammability question is asked, and that is the
result.

**Step 3 — source term.** These are two-phase flashing releases from 1 or
5 barg, which is the situation `validation/flashing.py` handles, though it was
built and checked on ammonia. Hydrogen differs in that the flash fraction is
large and rainout is the phenomenon the trials were designed to study, so the
"all liquid evaporates" assumption needs checking against the report rather
than being carried over.

**Step 4 — comparison, with censoring handled.** Saturated readings enter as
"at least 4 %", not as 4 %. That means the ordinary geometric statistics do
not apply and a censored-data treatment is needed — at minimum, reporting how
many pairs are censored and giving bounds; better, a survival-analysis
estimator. Reporting MG over a set where a sixth of the observations are
capped would be misleading in the model's favour.

**Step 5 — say where the model fails, and why.** Given steps 2 and 4, the
expected outcome is that DEGADIS over-predicts near the source and cannot
represent the far field at all. Establishing that quantitatively is a more
useful contribution than a set of statistics that pretend the comparison was
valid.

## What will not be attempted

- **Lift-off and buoyant rise.** Outside DEGADIS. Implementing it would be
  writing a different model and calling it a validation.
- **Rainout and pool spreading on the E3.4 substrates.** The concentration
  measurements there span 20 cm of height and no downwind distance; the
  dataset is about the source term, and DEGADIS represents substrate through a
  single heat-transfer coefficient.
- **The E3.1a discharge trials.** Near-orifice physics, 1.2 GB in nested
  archives, not extracted. Relevant to a source term, not to dispersion.

## The result this is likely to produce

Not "degali predicts LH₂ dispersion". More probably: a quantitative
statement of the concentration and distance beyond which a dense-gas model
cannot be applied to hydrogen, obtained with a code that has been shown to
reproduce the reference model exactly — so that the failure is attributable to
the model rather than to the implementation.

That is worth having. Regulatory guidance for hydrogen still borrows
dense-gas tooling built for LNG, and a defensible boundary on where that
borrowing stops being legitimate is more useful than another set of
performance measures.
