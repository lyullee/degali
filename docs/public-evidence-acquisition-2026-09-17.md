# Public LH2 evidence acquisition — 2026-09-17

## Purpose and handling boundary

This is an acquisition record for the next physical-model branch, not a
validation result.  Downloaded source archives remain outside version control
under `reference/`; no external raw measurements, report figures, or source
code from other projects are distributed with DEGALI.

The immediate physics priority is a pool source boundary: liquid inventory,
pool footprint, substrate cooling, and vapour mass flux must be separated
before a pool-fed atmospheric plume can be assessed.  A downstream transport
coefficient must not be fitted to compensate for an unspecified pool source.

## Public sources found

| Source | What is available | Direct use in DEGALI | Access state |
|---|---|---|---|
| PRESLHY E3.4, DOI [10.35097/1319](https://doi.org/10.35097/1319) | Ten LH2 pool tests on concrete, sand, gravel, and water; pool mass, 79 thermocouples, local H2 readings, and three-component wind at about 10 Hz.  401.8 MB, CC BY-SA 4.0. | Primary independent test of substrate-conduction/evaporation balance and crosswind sensitivity.  Its concentration probes cover only 0.35--0.55 m height, so it is **not** a far-field dispersion validation. | Already present, complete, in the local SLABX validation archive: ten unique XLSX test records plus the technical report.  Do not duplicate the raw measurements in DEGALI. |
| PRESLHY E3.5, DOI [10.35097/1481](https://doi.org/10.35097/1481) | 24 HSL outdoor rainout tests: release flow/pressure/temperature, thermocouples, H2/O2, weather, and video.  Archive is 11.34 GB because it includes video. | Independent near-/far-field concentration and temperature checks after source formation.  Existing local reductions and seven selected workbooks already cover the current thermal audit. | Already present, complete, in the local SLABX validation archive: 24 XLSX records and 24 videos.  Do not duplicate the raw measurements in DEGALI. |
| PRESLHY E3.1a (DisCha), DOI [10.35097/1187](https://doi.org/10.35097/1187) | More than 200 hydrogen blow-down experiments, including roughly 80 K cryogenic cases; selected mature-set Excel records for 0.5, 1, 2, and 4 mm nozzles and all tested pressures.  Published archive size is 1.3 GB, CC BY 4.0. | Independent source-boundary screen: reservoir pressure/temperature evolution, two-phase blow-down, and source temperature before any pool calculation. | Public download; accept the repository's ordinary download terms and keep the archive outside version control.  It is not an external-data request. |
| PRESLHY E3.1b (Cryostat), DOI [10.35097/1317](https://doi.org/10.35097/1317) | The accompanying report describes twelve reproducibility-tested blow-downs, including actual LH2 cases at about 20 K.  The published archive is 1.06 GB, CC BY-SA 4.0, but its exposed numeric workbooks are only the five 290 K reference cases (2/4 mm, 3.2--5 bar). | Measurement-operator and pressure-synchronisation reference; the report figures may be digitised as bounded evidence. | **Do not label this archive as LH2 raw validation.**  No 20--30 K numeric workbook was found in the published archive, so the low-temperature branch remains unvalidated by E3.1b raw data. |
| Sandia SAND2019-9998, [open report](https://h2fcp.org/sites/default/files/LH2_ColdPLUME_SAND_REPORT.pdf) | Calibrated planar-Raman concentration/temperature fields for 50--64 K, up-to-5-bar cryogenic hydrogen jets; the report also describes PIV velocity measurements. | Independent, controlled near-jet shape and thermal-profile screen.  It is not an LH2-pool test and cannot identify a pool evaporation law. | Already present locally as a reference report. |
| Sandia SAND2025-08847, pooling/vaporisation report | Pool containment, spread and 4% contour figures for the crosswind tests. | Falsification bounds for the pool-to-plume handoff; Test 5 already rules out the current all-inflow-as-immediate-vapour source. | Report is public; numeric raw series are stated to be available on request. |
| FFI/DNV Spadeadam large releases, [FFI report 21/03101](https://www.ffi.no/en/publications-archive/large-scale-leakage-of-liquid-hydrogen-lh2-tests-related-to-bunkering-and-maritime-use-of-liquid-hydrogen) | Public report of large outdoor releases: pool radius, field H2/temperature observations and directional limits. | Large-scale qualitative and contour/lower-bound checks, separate from the small controlled pool tests. | Report is public; raw DNV channels require an access enquiry/selection, not an automatic bulk download. |
| ELVHYS WP4.2, DOI [10.18710/JXJP0H](https://doi.org/10.18710/JXJP0H) | CC0 data from 48 cryogenic-H2 tests in a 1 m³ transfer-connection space, including small leaks and a 25.4-mm LH2 pressure-peaking class. | Future wall/obstacle and confined-source validation; not a substitute for an outdoor pool source boundary. | Already present, complete, in the local SLABX validation archive (196 CSV records).  Defer analysis until the obstacle branch has a pre-registered observation operator. |
| PRESLHY D4.8, [open report](https://hysafe.info/wp-content/uploads/sites/3/2020/12/PRESLHY_D4.8-report-4th-Dec.pdf) | HSE analysis of condensed air and the E3.5 source measurement.  It gives traceable source-flow ranges for 5-bar/1-bar and 6/12/25.4-mm configurations, and explains when the mass-flow meter is invalid in two-phase service. | A physics constraint on source-flow reconstruction and the N2/O2 condensed-phase branch; it prevents treating an invalid Coriolis trace as an observation. | Public report; use its reported ranges and uncertainty structure, not a fitted condensed-particle coefficient. |
| DNV/FFI large-release presentation, [open slides](https://hysafe.info/wp-content/uploads/sites/3/2021/05/3_1_Allason_NPRA-Large-Volume-Liquid-Hydrogen-Releases.pdf) | Test-by-test mean outflow, source vapour-quality ranges, release duration, observed pool evidence, concrete-temperature traces and LFL limits for the seven Spadeadam outdoor releases. | Independent large-release inequality/bound validation.  Plot values can be digitised with an explicitly recorded pixel-to-axis uncertainty. | Public secondary evidence only; it does **not** replace the unavailable original time-series channels. |
| Verfondern & Dienhart, [ICHS 2005 open paper](https://h2tools.org/sites/default/files/2019-09/110078.pdf) and [2007 open technical report](https://juser.fz-juelich.de/record/1311/files/Energie%26Umwelt_10.pdf) | LAuV shallow-layer pool equations, solid/water-ground treatment and legacy LH2-pool evidence. | Equation and range cross-check only.  Do not transplant LAuV coefficients into DEGALI. | Public substitutes for the paywalled 1997 article. |

## What is genuinely still unavailable

The following original measurements have not been publicly located.  They are
not a reason to request access: the associated branches remain bounded,
non-default physics until public evidence exists.

1. Cryogenic-H2 radial mean velocity together with three-component RMS (or a
   section-integrated `rho*u*k`) and an integral scale/dissipation rate.
2. Condensed N2/O2 mass fraction, mass-weighted particle-size distribution,
   and independently measured gas--particle slip at a matched LH2 release.
3. Sandia 2025 per-test numeric pool mass/area/evaporation and concentration
   records for the contained low-flow pool cases.
4. FFI/DNV Spadeadam machine-readable source, weather, pool and concentration
   channels for the unignited outdoor LH2 releases.  The local FFI report is
   complete, but no accompanying instrument archive is public.  Use the
   report and the public DNV slides only as digitised lower/upper bounds.
5. The nominal-20--30 K E3.1b Cryostat workbooks described by the technical
   report.  The public E3.1b package currently exposes its 290 K reference
   workbooks only; this is a publication-coverage fact, not a licence to
   infer a cryogenic result from the ambient records.

No public field report found here supplies those three quantities.  A visible
cloud, an infrared image, or a computed CFD velocity field cannot replace them
as calibration evidence.

## Public-only route

**No paid purchase is needed before the pool-boundary implementation starts.**
The only directly relevant paywalled paper found was Verfondern & Dienhart,
*Experimental and theoretical investigation of liquid hydrogen pool spreading
and vaporization*, *International Journal of Hydrogen Energy* 22 (1997),
649--660, DOI
[10.1016/S0360-3199(96)00204-2](https://doi.org/10.1016/S0360-3199(96)00204-2).
It would add original water/aluminium graphs, but the public 2005/2007
materials above cover the governing shallow-layer formulation and the new
PRESLHY E3.4 raw data are a better independent target.  Purchase it only if
we decide to validate water/metal substrate extrapolation quantitatively.

No external data-access request is part of this work.  The execution order is:

1. extract the local E3.4 pool mass and substrate histories.  This is now
   implemented as `degali.validation.preslhy_e34`; the new conduction-only
   pool boundary is `degali.addons.pool_evaporation`.  Both require declared
   intervals/properties and contain no raw record or fitted coefficient.
2. use local E3.5 raw records with the public HyWAM measurement description
   to preserve sensor transit-time and position constraints.  The E3.5 reader
   now maps Flowmeter and Xensor windows by their recorded clocks before it
   falls back to record fractions.  D3.6 documents approximately 18 s for the
   30 m pumped near-field lines, so 18 s is the documented default; a caller
   may replace it only with a separately declared measurement-line value,
   never one inferred from a concentration peak.
3. ingest the public E3.1a source archive as a separate source-boundary
   validation set;
4. digitise only the public DNV/FFI presentation plots needed for pre-declared
   large-release bounds; and
5. leave Sandia numeric pool records, direct turbulence statistics, and
   gas--particle slip outside calibration and report them as evidence gaps.
