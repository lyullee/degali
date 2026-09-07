# DEGALI — Dense Gas Dispersion for Liquid Hydrogen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22646258.svg)](https://doi.org/10.5281/zenodo.22646258)
[![PyPI](https://img.shields.io/pypi/v/degali.svg)](https://pypi.org/project/degali/)

DEGALI stands for **Dense Gas Dispersion for Liquid Hydrogen**. It is a modern
Python model built from a verified reimplementation of the US EPA DEGADIS 2.1
dense-gas dispersion model and extended for cryogenic hydrogen releases.

> **Alpha research software.** The DEGADIS 2.1 compatibility path is strongly
> regression-tested. The liquid-hydrogen extensions are suitable for research,
> scenario comparison and sensitivity studies, but are not independently
> certified for regulatory separation distances or safety-critical design.

한국어 요약: 현재 공개본은 원본 DEGADIS 재현 경로와 액화수소 연구 확장을
함께 제공합니다. PRESLHY 수평 야외제트 범위에서 정량 비교를 마쳤지만,
미완성 TKE/압력/열폭 폐쇄와 제한된 독립 검증 때문에 설비 인허가 판단의
단독 근거로 사용하면 안 됩니다. 자세한 판정은
[1차 결과](docs/stage1-results-2026-09-06.md)를 참고하십시오.

처음 설치하거나 실제 계산을 시작하려면 상세 사용 가이드를 참고하십시오:
**[한국어](USER_GUIDE_KO.md) · [English](USER_GUIDE.md)**. LH₂ 간편
평가, 결과 해석, CSV 저장, 민감도 분석, 기존 DEGADIS 입력 덱과 CLI 사용법을
단계별로 설명합니다.

## What is included

- Python ports of all six DEGADIS 2.1 programs. The independently obtained
  Fortran oracle is retained locally and is not redistributed.
- Steady, transient, receptor-dose and jet-to-ground workflows.
- Legacy and modern thermodynamic backends.
- Opt-in LH2 source flashing, buoyant trajectory, air condensation/freezing,
  component enthalpy, ground interaction and thermal-profile research paths.
- Experimental finite-TKE, independent thermal-width and Reynolds-stress
  operators. These are deliberately marked as research-only until their
  physical closure inputs and downstream field performance are validated.

## Current evidence

The port reproduces the five EPA reference cases and was checked locally
against a source-built Fortran implementation. The first frozen LH2 assessment uses
seven horizontal outdoor PRESLHY trials, 38 concentration sections, 17
vertical profiles and 41 temperature sensors.

| Quantity | Current provided path |
|---|---:|
| Concentration MG (observed/predicted; ideal 1) | 1.102 |
| Concentration VG (ideal 1) | 1.233 |
| Concentration FAC2 | 36/38 |
| Vertical width, predicted/measured | 1.091 |
| Plume-centre MAE | 0.054 m |
| Minimum-temperature MAE | 27.25 K |

An opt-in thermal-profile candidate lowers minimum-temperature MAE to 16.47 K
but worsens other acceptance metrics, so it is not the default. These samples
are correlated observations from one campaign, not 38 or 41 independent
validation experiments.

See [model status](docs/stage1-results-2026-09-06.md),
[claim grading](docs/claim-grading.md), and
[data/reproduction notes](docs/DATA_AND_REPRODUCTION.md).

## Installation

Install Python 3.10 or newer. For LH2 work, install DEGALI with the CoolProp
extra:

```bash
python -m pip install "degali[coolprop]"
```

Check the installed version and command-line interface:

```bash
python -c "import degali; print(degali.__version__)"
degali --help
```

CoolProp is required for LH2 calculations. The traditional DEGADIS
compatibility path can be installed with `python -m pip install degali`.

## Quick manual

### Liquid-hydrogen release without an input deck

The simplest public interface is `degali.lh2.assess`. The example below
represents a horizontal 0.1 kg/s release through a 10 mm orifice, 0.5 m above
ground, from saturated LH2 stored at 6 bar absolute:

```python
from degali.lh2 import assess

result = assess(
    rate=0.10,                  # kg/s
    wind=2.0,                   # m/s at release height
    height=0.50,                # m
    orifice=0.010,              # m; use pool_diameter instead for a pool
    storage_pressure=6.0,       # bar(a), not gauge pressure
    ambient_temperature=288.15, # K
    relative_humidity=65.0,     # percent
    ambient_pressure=101325.0,  # Pa
    max_distance=30.0,          # m
    at_distance=10.0,           # m
)

print(result.report())
for warning in result.warnings:
    print("WARNING:", warning)
```

The principal outputs are:

- `distance_to_lfl`: centreline distance to the hydrogen LFL of 4 mol%, m;
- `distance_to_stoichiometric`: centreline stoichiometric distance, m;
- `lowest_flammable_height`: lowest flammable-gas height, m;
- `regime`: `grounded`, `low`, or `aloft`;
- `trajectory`: NumPy columns `[distance, centre height, mole fraction]`;
- `warnings`: validation-range and applicability warnings. Do not discard them.

Specify exactly one source geometry: `orifice` for a pressurised jet or
`pool_diameter` for a pool/evaporation source. `storage_pressure` is bar(a),
while `ambient_pressure` is Pa. All temperatures are K.

### Save the LH2 trajectory

```python
import numpy as np

np.savetxt(
    "lh2_trajectory.csv",
    result.trajectory,
    delimiter=",",
    header="distance_m,height_m,centreline_mole_fraction",
    comments="",
)
```

### Existing DEGADIS input decks

```python
from degali import run_steady, run_transient, run_jet_to_ground

profile, source = run_steady("CASE.INP", backend="legacy")
print(profile.distance_to(0.05))

transient = run_transient("TRANSIENT.INP")
profile, jet, source = run_jet_to_ground("JET.INO", "GROUND.IN")
```

The equivalent command-line workflows are:

```bash
degali steady CASE.INP
degali transient TRANSIENT.INP --snapshot 60 --snapshot 120
degali dose TRANSIENT.INP --at 50 --at 100 --at 200
degali jet JET.INO --bridge GROUND.IN
```

DEGADIS decks are positional files: a missing value changes the meaning of
every field that follows. Their pressure convention also differs from
`assess()`. Preserve the original deck, verify units, and change one field at
a time. Input decks and licensed third-party data are not bundled.

For complete installation, input, output, plotting, parameter-study, deck,
validation-range, and troubleshooting instructions, see the root-level
**[English user guide](USER_GUIDE.md)** or
**[한국어 사용자 가이드](USER_GUIDE_KO.md)**.

## Verification

Fast development checks:

```bash
python -m pytest -m "not slow" -q
```

The original Fortran oracle and raw REDIPHEM, SMEDIS and PRESLHY files are not
redistributed. Tests that need
them skip with an explicit message unless the documented environment variables
are configured. Their provenance, reduction method and aggregate validation
results are documented. See [publication scope](docs/publication-scope.md).

The consolidated [technical reference](docs/technical-reference.md) describes
the governing physics, LH2 extensions, validation process, quantitative
results and remaining model-form limitations.

## Repository layout

```text
src/degali/  model and validation code
tests/       self-contained numerical, physical and regression tests
tools/       research and release verification utilities
docs/        derivations, provenance, audits, limitations and results
```

## Scope and safety

The current validated domain does not cover arbitrary equipment conditions,
obstacles, indoor releases, pool spreading, downward/strongly wind-steered
jets, or general transient source behaviour. Always compare safety-critical
results with independent experiments and an accepted consequence-analysis
workflow. See [security and safety reporting](SECURITY.md).

## Citation and license

Citation metadata, including the version DOI, are provided in
[`CITATION.cff`](CITATION.cff). The README DOI badge uses the concept DOI so it
continues to resolve to the latest archived DEGALI release.

The complete GitHub, Zenodo DOI and PyPI release sequence is documented in the
[publication guide](docs/publication-guide.md).

The DEGALI Python implementation is available under the MIT License. Original
Fortran implementations and external experimental datasets are not bundled or
redistributed.

The chronological development log that previously occupied this front page is
preserved in [docs/README-development-log-2026-09-06.md](docs/README-development-log-2026-09-06.md)
and [CHANGELOG.md](CHANGELOG.md).
