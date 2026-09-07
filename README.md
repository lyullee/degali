# DEGALI — Dense Gas Dispersion for Liquid Hydrogen

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

Create a Python 3.10 or newer environment and install from a clone:

```bash
python -m pip install -e ".[test]"
```

CoolProp is optional for the DEGADIS 2.1 compatibility path and required for
most cryogenic-hydrogen calculations:

```bash
python -m pip install -e ".[coolprop]"
```

## Quick start

```python
from degali import run_steady, run_transient, run_jet_to_ground

profile, source = run_steady("B9.INP")
print(profile.distance_to(0.05))

transient = run_transient("B9T.INP")
profile, jet, source = run_jet_to_ground("EX2.INO", "EX2.IN")
```

The command-line interface exposes the same main workflows:

```bash
degali steady B9.INP
degali transient B9T.INP --snapshot 60 --snapshot 120
degali dose B9T.INP --at 200 --at 400 --at 800
degali jet EX2.INO --bridge EX2.IN
```

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

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Before the
first archived release, the maintainer must add the reserved DOI and release
date.

The complete GitHub, Zenodo DOI and PyPI release sequence is documented in the
[publication guide](docs/publication-guide.md).

The DEGALI Python implementation is available under the MIT License. Original
Fortran implementations and external experimental datasets are not bundled or
redistributed.

The chronological development log that previously occupied this front page is
preserved in [docs/README-development-log-2026-09-06.md](docs/README-development-log-2026-09-06.md)
and [CHANGELOG.md](CHANGELOG.md).
