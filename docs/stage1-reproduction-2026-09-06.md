# 1차 실행 방법

이 명령은 **이번 고정 7개 시험의 재현용**이다. 새로운 설비 입력을 자동으로
보정하거나 검증해 주는 설비 설계 프로그램은 아니다. 결과 해석은
`docs/stage1-results-2026-09-06.md`를 먼저 읽는다.

## 현재 PC에서 실행

PowerShell에서 프로젝트 폴더로 이동한다.

```powershell
Set-Location -LiteralPath '<project-root>\work\degali'
```

1. **파일과 실행 환경 확인** — 새 모델 계산이나 파일 변경 없이 검증한다.

```powershell
& '.\.venv\Scripts\python.exe' tools/stage1.py verify
```

`runtime_matches: true`이고 오류 없이 끝나야 한다. 이 검사는 납품 스냅샷의
일치 확인이지 물리식 재검증이나 모든 후보의 새 계산은 아니다. 과거 체크포인트의
알려진 실행 도구 해시 차이는 최종 보고서에 별도로 남아 있다.

2. **제공 기준 경로를 실제로 다시 계산** — 지면 반사 중복을 제거한
`BASE_SINGLE_IMAGE`로 7개 시험을 새로 적분한다.

```powershell
& '.\.venv\Scripts\python.exe' tools/stage1.py baseline --output-dir reference/preslhy/stage1_runs/my_baseline_01
```

결과는 지정한 새 폴더의 `baseline_single_image.json`이다. 7개 전체 12열 출력,
실행 입력, 38개 농도 비교, 17개 높이/폭 비교, 41개 온도 비교와 오차를 포함한다.
관측값이 들어 있는 재현 결과이며, 순수 무관측 예측 산출물로 혼동하지 않는다.
`maximum_saved_correction_reproduction_error`는 기준 감사 결과와의 차이다.

3. **연구 후보 대표 두 사례를 다시 계산** — 기본판 승격이 아니다.

```powershell
& '.\.venv\Scripts\python.exe' tools/stage1.py candidate --output-dir reference/preslhy/stage1_runs/my_candidate_01
```

시험10·23의 공급원·인터페이스·하류 궤적을 순서대로 다시 계산해 각각
`candidate_10.json`, `candidate_23.json`으로 저장한다. 전체 7개 신규 후보
실행이 아니다. 기준 경로보다 오래 걸릴 수 있다. 촘촘한 두 사례 진단은
다른 새 폴더를 지정하고 `--maximum-step 0.01`을 덧붙인다.

이미 존재하는 출력 폴더/결과를 덮어쓰지 않는다. 다시 실행할 때는
`my_baseline_02`처럼 새 이름을 쓴다. 중단했다면 부분 결과를 보존하고 새
폴더로 다시 시작한다. 다운로드·설치·삭제·외부 전송은 자동으로 하지 않는다.

## 고정 설정과 단위

- 전체 입력은 최종 assessment의 `executed_source_conditions`에 있다.
  BASE의 유량 kg/s, 지름/높이 m, 온도 K, 상대습도 %, 압력 barg를 구분한다.
  `storage_temperature: null`은 누락된 측정값을 0K로 쓰는 것이 아니라
  기존 함수가 탱크 압력의 포화액 온도를 계산하도록 맡기는 설정이다.
- BASE: 기존 보정 옵션 및 공급원 표/운동량 일관성 사용, 출력 간격 .2m,
  arclength limit 40m, 기존 허용오차 1e-4. .2m를 내부 최대 적분 간격이라고
  부르지 않는다. 실제 최종 출력 x는 적분 출력 방식상 40m를 조금 넘을 수 있다.
- 새 어댑터: 양의 높이·정상 폭·유한한 12열 JetPlume 보고 행만 허용.
  저장 노드에서 반사 전 농도로 복원 후 보간한다. 지면 도달이나 범위 밖
  센서를 조용히 대기온도로 대체하지 않는다.
- CANDIDATE: `full` 측정 배관 공급원, collective equilibrium bound,
  normal hydrogen, explicit ambient species, enthalpy profile,
  `polar_square_all_flux_v1`, 속도폭 비 1.16, 최대 간격 .02m, 하류 허용오차 1e-5.
- Python 실행 내부에서 현재 자식 프로세스의 경로와 수치 라이브러리 스레드만
  설정한다. PC 전역 환경변수/전원 정책/PowerShell 실행 정책은 바꾸지 않는다.

## 저장된 결과를 직접 확인하려면

모두 `reference/preslhy` 아래다.

| 파일 | 실제 역할 |
|---|---|
| `stage1_assessment_2026-09-06.json` | 최종 판정, 네 구성 점수, 공급원 설정, 전체 열 정의, 재현 검사 |
| `stage1_existing_evidence_2026-09-06.json` | 과거 분할 결과·입력·핵심 코드·원시 온도 파일 해시와 38/17 집계 |
| `stage1_temperature_2026-09-06.json` | 새 BASE 7개와 같은 41개 센서의 BASE/CONTROL/CANDIDATE 비교 |
| `stage1_reproduction_10_2026-09-06.json`, `stage1_reproduction_23_2026-09-06.json` | 후보 대표 두 사례 .02m 새 전체 궤적 |
| `stage1_resolution_diagnostic_2026-09-06.json` | BASE .2/.1m, 후보 .02/.01m 대표 두 사례 비교 |
| `stage1_single_image_2026-09-06.json` | 중복 반사 보정과 같은 38/17/41 평가 |
| `stage1_delivery_smoke_2026-09-06/baseline_single_image.json` | 납품 명령으로 실행한 보정 BASE 7개 새 적분 결과 |
| `stage1_delivery_tests_2026-09-06.xml` | 납품용 수정·도구 18개 검사 결과 |
| `stage1_delivery_manifest_2026-09-06.json` | 최종 파일 SHA-256 목록 |

`stage1_temperature`의 기존 `output_columns`는 처음 5열 이름만 담고 있지만
실제 행은 12열이다. 전체 정의는 assessment의 `reporting_schema_correction`
또는 새 실행 결과를 따른다. 기존 파일을 임의로 고치면 해시 검사가 실패한다.

## PC 이동과 재현 한계

현재 확인 환경: Windows 11, Python **3.12.14**, NumPy **2.3.5**,
SciPy **1.18.1**, CoolProp **8.0.0**, pytest **9.1.1**.
프로젝트 내부 경로는 실행 도구 위치를 기준으로 찾으므로 전체 프로젝트를
옮겨도 적용할 수 있다. 단 Windows 가상환경 폴더 자체는 이동만으로 재사용이
보장되지 않는다. 새 PC에서는 해당 버전의 별도 Python 환경을 구성해야 한다.
다른 라이브러리/OS의 숫자 일치까지 이번 결과가 보증하지 않는다.

`src`, `tools`, `docs`, `tests`, `reference/preslhy`의 납품 manifest 대상 파일과
관련 원시 자료를 함께 보존한다. 현재 프로젝트에는 원시 시험 자료가 들어 있으므로
외부 공유·논문 부록 배포 전에는 원자료의 이용/배포 조건을 별도로 확인한다.
오류가 나면 명령 전체와 오류 문구, `verify` 결과를 알려주면 된다.
