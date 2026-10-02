# LH2 paper baseline — 2026-09-20

## Frozen configuration

This fixes the fast LH2 model path for the validation paper. It is not a new public API and contains no validation fit.

| Item | Frozen choice |
|---|---|
| Source path | corrected expanded flashing source |
| E3.5 source rate | reported window mean (`flow_mean_gs`) |
| Concentration operator | evaluate every sensor, then use the same arc maximum as the measurement |
| Ground contact | free detachment (`ground_effect=False`) |
| Surface roughness | 0.001 m |
| FFI wind | published low-mast value; 10 m reference height is recorded assumption |
| Fitted coefficients | none |

This scope is limited to momentum-dominated horizontal LH2 releases. It excludes impinging downward jets, ventilation-mast outlets, obstacle wakes and fully transient clouds.

The one-call API now makes the boundary executable:

```python
from degali.lh2 import assess, assess_envelope

screen = assess(
    rate=0.285, orifice=0.0254, storage_pressure=6.0,
    wind=2.5, height=0.5, max_distance=6.0,
    strict_scope=True,
)

sens = assess_envelope(
    rates=[0.25, 0.285, 0.32],
    winds=[2.3, 2.5, 2.7],
    orifice=0.0254, storage_pressure=6.0, height=0.5,
    max_distance=6.0,
)
print(sens.distance_to_lfl_range)
```

The envelope is an explicit input-sensitivity calculation, not a fitted
uncertainty band. `strict_scope=True` refuses any evidence-range warning;
warning mode remains available for exploratory extrapolation.

## Reproducible audit

The audit writes only local derived residuals and maps. It does not copy a workbook, report PDF or time series into the distribution.

```powershell
$env:PYTHONPATH = 'src'
.\.venv\Scripts\python.exe tools\audit_lh2_paper_baseline.py `
  --output-dir tmp\lh2-paper-validation-2026-09-20
```

It creates SHA-256 fingerprints, concentration and thermal residual maps, FFI arc and sensor outputs, an FFI Test 6 decomposition, and a separate Hecht–Panda slope-boundary map.

## Frozen results

The primary E3.5 score uses 62 arc maxima from nine momentum-dominated horizontal releases:

| n | MG (observed/predicted) | VG | FAC2 | FB | NMSE |
|---:|---:|---:|---:|---:|---:|
| 62 | 1.047 | 1.425 | 0.839 | +0.068 | 0.122 |

The local residual map retains 197 sensor records; 190 have positive model concentrations. These receptors show distance, height and crosswind structure, but are not independent validation trials.

The independent FFI/DNV horizontal-release screen contains six arc maxima (Tests 4 and 6 at 30, 50 and 100 m): MG 1.245, VG 1.373 and FAC2 0.833. It is a small far-field screen, not a general field-accuracy claim.

The thermal map re-expresses a sealed E3.5 diagnostic with two trials, two stations and five heights. It does not re-integrate the primary trajectory and cannot promote a thermal-transport correction.

The public Hecht–Panda Raman benchmark is kept as an external boundary rather than merged into the primary score because the journal fit contains an unresolved condition-count/membership discrepancy. The current LH2 branch is above the reported aggregate slopes by factors 2.09 (centreline mass), 2.72 (mass width), 2.41 (centreline temperature) and 3.46 (temperature width). This is evidence against a universal warm/free-jet transfer claim, not a licence to fit those slopes.

## Provenance

- Lyons, Coldrick and Atkinson, *Summary of experiment series E3.5*, DOI [10.5445/IR/1000136281](https://doi.org/10.5445/IR/1000136281).
- Aaneby, Gjesdal and Voie, *Large scale leakage of liquid hydrogen (LH2)*, FFI Report 20/03101, [public report](https://www.ffi.no/publikasjoner/arkiv/large-scale-leakage-of-liquid-hydrogen-lh2-tests-related-to-bunkering-and-maritime-use-of-liquid-hydrogen/21-03101.pdf).
- Hecht and Panda, *Study of key parameters in modeling liquid hydrogen release and dispersion in open environment*, IJHE 44 (2019), DOI [10.1016/j.ijhydene.2018.07.058](https://doi.org/10.1016/j.ijhydene.2018.07.058).

The repository contains methods, provenance and derived summaries. It does not redistribute third-party raw workbooks, raw sensor time series or original DEGADIS Fortran source.
