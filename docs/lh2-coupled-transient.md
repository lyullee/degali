# LH2 coupled transient research path

`LH2_COUPLED_TRANSIENT` is an opt-in orchestration layer. It does not replace
or modify the frozen `DEGADIS_21` compatibility path.

The model starts from one pressure-thrust flash plane and routes its phases at
the same source clock:

1. direct flash vapour -> conserved near field -> yawed crosswind plume ->
   finite Gaussian puff;
2. residual liquid -> finite-size droplet flight -> in-flight evaporation,
   airborne liquid or ground rainout;
3. ground rainout -> concurrent fixed or axisymmetric spreading pool ->
   time-binned pool evaporation.

The direct-vapour adapter preserves its flash-plane mass flow, common phase
velocity and equilibrium vapour enthalpy. The liquid inventory is not
evaporated or discarded to create that gas source.

## Deliberate boundary

In-flight droplet vapour and pool vapour have conserved masses and time bins,
but the current evidence does not uniquely provide their atmospheric launch
temperature, vertical momentum and concentration profile. They are therefore
returned as `LH2TransientAtmosphericSource` entries with
`dispersion_model=None`. They are not silently merged into the nozzle jet.

This distinction appears directly in the result:

- `result.conservative` tests the end-to-end hydrogen ledger and the native
  jet/puff conservation diagnostics;
- `result.atmospherically_complete` is false while any vapour source lacks a
  dispersion boundary;
- `result.accepted` additionally applies the existing gas-branch
  applicability screen;
- `result.unresolved_atmospheric_mass_kg` quantifies the missing atmospheric
  continuation instead of hiding it in a warning.

Passing the conservation screen is not field validation.

## Public API

```python
from degali.addons import SolidSubstrate, flashing_hydrogen_droplet_source
from degali.lh2 import run_lh2_coupled_transient_research

flash = flashing_hydrogen_droplet_source(
    mass_flow=0.285,
    orifice_diameter=0.006,
    upstream_temperature=28.2550342766,
    upstream_pressure=600000.0,
    upstream_quality=0.0,
)
ground = SolidSubstrate(
    conductivity_w_m_k=1.4,
    density_kg_m3=2200.0,
    heat_capacity_j_kg_k=850.0,
    initial_temperature_k=293.15,
    depth_m=0.5,
)

result = run_lh2_coupled_transient_research(
    flash,
    release_duration_s=5.0,
    post_release_duration_s=60.0,
    puff_duration_s=60.0,
    release_position_m=(0.0, 0.0, 1.5),
    release_azimuth_rad=0.0,
    release_elevation_rad=0.0,
    wind_speed_m_s=2.5,
    wind_to_angle_rad=0.0,
    evaporation_coefficient_m2_s=1.0e-9,
    pool_area_m2=1.0,
    pool_time_step_s=0.1,
    substrate=ground,
)
print(result.report())
```

The gas branch currently accepts a horizontal nozzle. Arbitrary horizontal
bearing and global x/y source offsets are retained. Measured vector wind can
be supplied with `puff_wind_history`. Existing subordinate controls remain
available through `gas_model_options` and `phase_model_options`; shared
scenario values cannot be overridden through those dictionaries.

## Current scope

The calculation remains research-only. It does not resolve non-axisymmetric
pool spread, terrain slope, curbs, drains, wind shear on the liquid surface,
turbulent droplet dispersion, combustion, or obstacle interaction. The
unresolved source entries are the explicit interface for a future measured or
independently validated pool-vapour and distributed-droplet dispersion
closure.
