# LH2 pool heat balance and PRESLHY E3.4 observation boundary

## What changed

The ordinary `degali.lh2.assess()` pool route needs a **vapour** mass rate.
It does not turn a rainout rate into vapour automatically.  The new
`degali.addons.pool_evaporation` module supplies a transparent post-rainout
source boundary: a liquid pool held at hydrogen saturation temperature draws
heat through a homogeneous solid substrate, and that heat is divided by the
declared latent heat to obtain its vapour rate.

\[
 \rho_s c_{p,s}\frac{\partial T_s}{\partial t}
 = k_s\frac{\partial^2T_s}{\partial z^2},\qquad
 q''=-k_s\left.\frac{\partial T_s}{\partial z}\right|_{z=0},\qquad
 \dot m_v=\frac{Aq''}{h_{fg}}.
\]

The conduction equation uses a backward-Euler, one-dimensional cell-centred
grid.  Its surface is fixed at the declared LH2 saturation temperature while
liquid remains, and its deep boundary stays at the declared initial substrate
temperature.  The result reports the depth-to-thermal-penetration ratio, so a
finite numerical depth cannot silently be described as semi-infinite.

This is not a fitted evaporation correlation.  No E3.4 rate, material
constant, heat-transfer coefficient, wind correction, or visual pool estimate
is embedded in the package.

The conduction boundary is also independently consistent with the
non-spreading LH2 concrete-pool experiment of Zhou *et al.*: that study
separates container and concrete heat inputs and reports the container term as
small relative to the concrete term for its apparatus.  It supports the
*structure* of this boundary, not transfer of its apparatus-specific material
values or an interfacial boiling coefficient into DEGALI.  See
[Processes 11, 1415 (2023)](https://doi.org/10.3390/pr11051415).

## Use it as a source boundary

```python
from degali.addons import SolidSubstrate, substrate_conduction_evaporation
from degali.lh2 import assess

concrete = SolidSubstrate(
    conductivity_w_m_k=1.5,
    density_kg_m3=2200.0,
    heat_capacity_j_kg_k=900.0,
    initial_temperature_k=293.15,
    depth_m=1.0,
    cells=96,
)
pool = substrate_conduction_evaporation(
    concrete,
    area_m2=0.25,
    duration_s=300.0,
    time_step_s=1.0,
    initial_liquid_mass_kg=1.0,
)

# `assess` is a steady plume assessment.  Choose a declared time window or
# snapshot; do not average a rapidly changing source and call it transient.
rate = pool.steps[-1].vapour_rate_kg_s
answer = assess(rate=rate, wind=2.0, pool_diameter=0.50)
```

For a changing pool source, `assess_pool_history(pool, ...)` evaluates each
positive evaporation step as a separate **quasi-steady** plume.  It can label
a receptor arrival interval only when the caller declares a propagation speed.
This preserves the source-step duration instead of misrepresenting it as an
instantaneous release.  It is a
source-to-plume sequence, not a claim of validated transient cloud storage or
puff transport.

Use independently supported substrate properties and report them.  The
default latent heat is `4.46e5 J kg-1`; it too is an explicit argument.  For
early-time mesh checks, compare the numerical result with
`semi_infinite_heat_flux()` at a positive elapsed time.  Refine `cells` until
the selected observation time is mesh-insensitive; the supplied model does
not claim that `cells=48` is universally sufficient.

## What this model does not claim

It is limited to a substantially fixed-footprint, post-rainout pool above a
homogeneous solid.  The following require a separately declared model and are
not inferred by this module:

- source flashing, rainout fraction, incoming-liquid enthalpy and pool spread;
- water/ice phase change, porous-media imbibition, sand or gravel displacement;
- wind-driven surface transfer, film boiling/nucleate boiling, and radiation;
- a moving liquid level, finite rim overflow, or local impingement heat input;
- particle size, gas--particle slip, and turbulent transport.

If a finite initial liquid mass is supplied, the module stops generating
vapour when it is exhausted.  It deliberately has no hidden replenishment.

## E3.4 raw observation operator

The optional `degali.validation.preslhy_e34` reader consumes a user-held copy
of the public PRESLHY E3.4 workbooks (Lyons, Coldrick and Atkinson, DOI
[10.35097/1319](https://doi.org/10.35097/1319), CC BY-SA 4.0).  Raw files are
not distributed by DEGALI.

```python
from degali.validation.preslhy_e34 import (
    read_e34_pool_history, evaporation_window,
)

history = read_e34_pool_history("/local/PRESLHY/20200324-Concrete02-Final.xlsx")
observed = evaporation_window(
    history,
    start_s=950.0,
    end_s=1050.0,
    thermocouples=("TA000-00", "TA000-05"),
)
print(observed.evaporation_rate_kg_s)
```

The reader uses the campaign-provided corrected `m(LH2) [g]` channel, rather
than differentiating total facility weight.  It preserves the workbook time
axis and keeps all thermocouples alongside that mass history.  The caller
must predeclare the evaporation interval.  E3.4 contains filling, sloshing,
multiple release and sometimes substrate-loss intervals; choosing the
steepest or best-agreeing mass slope would be fitting the source to the model.

The returned slope `R²` and selected-temperature span are diagnostics, not
automatic pass/fail metrics.  Sand and gravel are marked for manual review
because the report and corrected channels themselves identify material-loss
artefacts.  Water is retained as an observation case but is outside this
solid-substrate model until a separately validated water/ice boundary exists.

`tools/audit_e34_pool_source.py` makes this comparison reproducible without
placing either a workbook or its derived observations in the repository.  Its
external JSON declaration must name the workbook, its exact `trial_id`,
pool-contact time, fixed comparison window, pool area and independently
selected substrate properties.  A `report_trial_id` can be supplied only when
it is identical to `trial_id`; the command rejects a report example applied to
a merely similar workbook name.
It writes a new JSON report only when the caller provides an unused output
path, and records that no coefficient was fitted.

### Provenance gate for the local public copy

The detailed D3.5 heat-balance narrative identifies its worked example as
`20200423-Concrete02`.  The local public workbook collection presently holds
`20200320-Concrete01`, `20200324-Concrete02`, and
`20200417-Concrete03Wind`, but not a workbook with the report example's date.
Therefore the example's plotted time intervals must **not** be applied to the
similarly named `20200324-Concrete02` workbook by assumption.  The audit tool
remains usable for a separately declared, internally self-consistent window,
but no report-to-workbook numerical comparison is promoted until the source
identity is demonstrated from the records themselves.

The public D3.3 release-and-mixing report independently confirms that the
worked case was filled repeatedly and documents the mass/temperature plot for
`20200423-Concrete02`.  It is useful corroboration of the experiment identity
and of the pool fill/dry-out sequence, but its figure is not a substitute for
the missing time-series workbook: digitising it would not establish the scale
correction, sample clock or the report-to-workbook mapping.  It is retained as
provenance support only, not converted into a fitted numerical data series.

## Validation sequence

1. Freeze each E3.4 interval from the mass and thermal record before looking
   at model output.
2. Run a material-property sensitivity stated from independent sources;
   do not optimise conductivity, heat capacity or initial temperature to E3.4.
3. Compare mass loss and its time dependence first.  Only then pass the
   resulting vapour rate to the atmospheric pool plume.
4. Keep the E3.4 0.35--0.55 m local hydrogen probes separate from a far-field
   validation score.  They do not observe the complete downstream plume.
5. Use E3.5 as the independent outdoor release/transport data set.  Its
sampling and gas-tube timing must remain in its own observation operator.

This separation prevents downstream entrainment or lift-off parameters from
compensating for an unspecified pool energy balance.
