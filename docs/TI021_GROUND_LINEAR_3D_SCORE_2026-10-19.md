# TI021 3-D ground-linear 전체 검증 — 2026-10-19

## 실행

```powershell
.\.venv\Scripts\python.exe tools/run_thorney_ti021_inner_transient.py `
  --workbook ti21_input_ascii.xls --duration-s 600 --cell-size-m 5 `
  --inner-domain-m 450 640 40 --time-step-s 0.125 `
  --receptor-vertical-sampling ground_linear `
  --output-directory outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18

.\.venv\Scripts\python.exe tools/score_thorney_ti021_nested_145_sensors.py `
  --workbook ti21_input_ascii.xls `
  --inner-run-directory outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18 `
  --outer-run-directory outputs/thorney-ti021-wide-outer-domain-600s-x450y640-centered-2026-10-13-v1/outer `
  --output-directory outputs/thorney-ti021-groundlinear-nested-score-600s-2026-10-19
```

## 결과

- 600 s, 5 m cell, 90×128×8 grid, 4,800 steps
- solver runtime: **3044.86 s** (약 50.7분)
- 145-sensor nested dose MAPE: **84.45845%**
- observed-positive sensors: 84; predicted detections at 0.1 vol%: 53
- normalized dose bias: **−0.88678**
- selected-source breakdown: inner 77 positive sensors, MAPE 85.05%; outer 7 positive sensors, MAPE 78.00%
- overlap coupling gate: **fail** (peak/dose change >10%)
- decision: **withheld**, not promoted

## 보존성 및 수치 게이트

- species mass residual: `4.15e-15` relative
- global mass residual: `4.57e-9 kg/s`
- continuity residual: `3.08e-12 kg m-3 s-1`
- momentum balance residual: `1.92e-8 kg m/s`
- thermodynamic corrector non-converged steps: `0`
- all solver gates: passed

## 해석

`ground_linear`은 첫 5 m 셀 중심에 갇히는 저고도 관측을 완화했지만, 전체 dose를 회복하지 못했다. 3-D inner와 outer 모두 관측 dose를 크게 과소예측하므로 문제는 단순한 수직 관측 연산자보다 **근거리 source/transport 질량이 센서 영역으로 전달되지 않는 것**에 있다. 센서별 보정계수나 관측값 피팅은 사용하지 않았다.

다음 후보는 동일 600 s 재계산 전에 짧은 보존성 진단으로 초기 방출 운동량·근거리 dense-gas transport 및 inner-to-outer continuation의 기여를 확인해야 한다. 이 결과는 10% 목표를 충족하지 않으며 승격하지 않는다.
