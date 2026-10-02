# Obstacle-wake evidence boundary

## What is now available

The geometry screen remains the package's operative obstacle/wall result: it
uses exact intersections with axis-aligned cuboids and transverse walls and
marks an unobstructed plume prediction invalid after contact. It does not
invent a wake concentration change.

On 2026-09-17, a public quantitative starting benchmark was identified:
AIJ UWE Case H on Zenodo, DOI
[10.5281/zenodo.15429872](https://doi.org/10.5281/zenodo.15429872). It is a
neutral wind-tunnel experiment around a 0.10 m × 0.10 m × 0.20 m cuboid, with
a near-floor ethylene point source. Its published `RS_caseH.csv` contains
three-component mean velocity, velocity RMS, TKE, mean concentration and
concentration RMS; its approach-flow file supplies the inflow profile. The
record is CC BY 4.0.

`degali.addons.read_aij_case_h` is the first quantitative-wake intake. It
reads only a user-supplied local CSV, validates the published schema and
preserves missing values. The benchmark is not packaged or redistributed.

## Strict scope

Case H permits a future **neutral-gas, single-cuboid** wake observation test.
It does not identify a cryogenic density correction, an LH2 source boundary,
liquid/condensed-phase interaction, or a universal wake coefficient. Its
benchmark object therefore explicitly reports that an LH2 quantitative-wake
prediction is not allowed.

For dense-gas obstacles, the public [SMEDIS database](https://admlc.com/smedis-dataset/)
provides downloadable Excel batches containing source, meteorology and
concentration measurements, including Thorney Island and other complex cases.
The openly supplied EEC360/EEC361 propane pair is a checked control/single
2 m-fence pair: release conditions, ambient metadata and all 26 sensor
channels correspond. Its nonzero control observations include both increases
and decreases behind the fence. This falsifies a scalar "wall attenuation"
factor before one can be added to the model.

`degali.addons.read_smedis_fence_pair` is a no-fit intake for a user-held copy
of such a declared pair. It rejects a changed source, meteorology, fence count
or sensor layout and retains repeated colocated channels. It only exposes the
published paired observations; it does not yield an LH2 obstacle prediction.
The spreadsheets are not bundled or redistributed.

The public obstacle field literature also rules out a scalar ``wall factor``:
a wall normal to the wind can dilute ground-level dense gas while increasing
some elevated concentrations through plume lifting.  Thus even a correct
ground-level attenuation would be wrong for a three-dimensional concentration
claim.  This reinforces the separation between the current contact screen and
a future, observation-constrained wake model.

## Required progression

1. Verify the AIJ local CSV against its published checksum and reproduce only
   its stated neutral wake observations.
2. Define a no-fit neutral wake observation operator for the one cuboid.
3. Use the paired SMEDIS observations to reject any proposed reduced wake
   closure that cannot reproduce its spatially signed effect without fitting.
4. Only if the neutral and dense regimes pass their distinct predeclared
   checks may a clearly labelled, non-default reduced wake closure be exposed.

Until then, `site_geometry` is the only operative site feature: it is a
validity screen, not a CFD surrogate or a facility-design concentration model.
