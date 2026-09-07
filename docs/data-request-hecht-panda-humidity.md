# Exact data request: Hecht--Panda Raman provenance, humidity and reduced data

## Critical fit-provenance request

The final journal manuscript has a source inconsistency that must be resolved
before treating the four aggregate slopes as a strict validation target.
Table 1 and the radial-profile figures contain nine conditions, but the
legends of Figures 6 and 8 also list `4 bar, 45 K, 1.25 mm`. The conclusion
states a 50--61 K range, and the prose continues to refer to nine conditions.

Please request:

1. whether the `4 bar, 45 K, 1.25 mm` run was included in each of the four
   aggregate fits, or whether the legend entry is erroneous;
2. its Table-1-equivalent record: number of imaging heights, nozzle pressure
   and temperature, throat temperature/pressure/density/velocity, run ID and
   acquisition range;
3. the run IDs, point counts and weights used in each centreline and half-width
   regression, plus fitted intercepts and parameter uncertainty/covariance;
4. confirmation of why the final journal mass slopes (`0.2771`, `0.07069`)
   differ from the conference values (`0.2626`, `0.06503`) while the two
   temperature slopes are unchanged; and
5. which atmospheric temperature was used for each centreline-temperature
   normalization, given the 292--300 K per-condition radial fits.

## Why this item is now blocking

The paper's dry-air equilibrium calculation underpredicts the printed
centreline temperature-decay slope by about 31%. The new conservative
100%-RH frost bound overpredicts it by about 45%, while an unscored 40%-RH
sensitivity is about 9% low and puts all four printed slopes inside the
pre-registered 25% band. Relative humidity is therefore a first-order model
input, not a minor nuisance parameter.

The Sandia publication page, OSTI record 1529288, accepted manuscript and
ICHS 2017 paper provide no RH, dew point, absolute humidity, test timestamp or
separate environmental-data attachment. Do not infer or fit an RH from the
four aggregate published regressions.

## Minimum useful request

Ask Sandia for the following, in priority order:

1. Relative humidity **or dew point at the 19 cm co-flow honeycomb inlet** for
   each of the nine Table-1 release conditions. A single laboratory value is
   still useful if all runs occurred in one stable session.
2. Ambient/co-flow air temperature paired to that humidity measurement and
   the sensor accuracy, location and sampling interval.
3. Whether the 0.3 m/s co-flow was untreated room air, facility compressed
   air, bottled air, or air passed through a dryer/humidifier. If dried, ask
   for the specified outlet dew point.
4. The timestamp or run identifier mapping each Table-1 condition to the
   environmental log.

Any one of RH, dew point, water-vapour partial pressure, humidity ratio
`kg H2O/kg dry air`, or a specified dryer dew point is sufficient. The model
already accepts humidity ratio directly.

## Strongly preferred validation data

To replace regression-on-four-printed-slopes with a defensible uncertainty
analysis, also request:

- reduced median H2 mole-fraction, N2 mole-fraction and temperature arrays for
  all 2--6 stitched heights in each condition, including pixel/mm coordinates;
- per-height 25th/75th percentile or the 400 processed single-shot arrays;
- nozzle pressure and temperature time series aligned to image acquisition;
- calibration constants, masks and background/flat-field arrays used for both
  cameras;
- the applied lateral slice shifts and any flags for suspected out-of-plane
  motion or nozzle ice;
- a statement on whether visible condensed moisture originated mainly in the
  co-flow, room air entrainment, or external nozzle frosting.

Preferred machine-readable formats are HDF5, NetCDF, NumPy NPZ, MATLAB MAT,
or CSV plus a field dictionary. Original image TIFFs are not required for the
first humidity-conditioned model test if reduced arrays are available.

## Record identifiers and contact

- Article: *Mixing and warming of cryogenic hydrogen releases*
- DOI: `10.1016/j.ijhydene.2018.07.058`
- Sandia report number: `SAND-2018-7244J`
- OSTI ID: `1529288`
- ICHS paper: `ID123`, 2017
- Corresponding author printed in the paper: `ehecht@sandia.gov`
- OSTI record: <https://www.osti.gov/biblio/1529288>
- Sandia record:
  <https://www.sandia.gov/research/publications/details/mixing-and-warming-of-cryogenic-hydrogen-releases-2019-04-02/>

## Ready-to-send request

> Subject: Environmental humidity and reduced Raman data for SAND-2018-7244J
>
> I am validating a mass/momentum/species/total-energy-conserving cryogenic
> hydrogen jet model against the tabulated conditions in Hecht and Panda,
> DOI 10.1016/j.ijhydene.2018.07.058. The published paper notes condensed
> moisture and nozzle ice, but I could not locate the co-flow/laboratory
> relative humidity or dew point. A conservative phase-equilibrium sensitivity
> shows this input materially changes the predicted centreline warming.
>
> I also found that Table 1 and the radial-profile panels contain nine runs,
> while the aggregate-fit legends list an additional 4 bar, 45 K, 1.25 mm
> series. Could you confirm whether that run was included in each of the four
> regressions and provide its Table-1-equivalent throat/source record? The
> run IDs, point counts, weighting, intercepts and fit uncertainty would let
> me reproduce the published 2019 slopes; confirmation of the change from the
> 2017 mass slopes would also be valuable.
>
> Could you provide, if retained, the RH or dew point and paired temperature at
> the 19 cm co-flow inlet for the nine Table-1 runs, the co-flow air source or
> dryer specification, and the run/timestamp mapping? Reduced H2/N2/temperature
> arrays for each stitched imaging height and their uncertainty/percentile
> information would also enable a much stronger validation than the four
> aggregate published fits. Any native machine-readable format is welcome.
