# DEGALI 개발·검증 업데이트 요약 — 2026-10-19

이 문서는 최근 TI021 장애물·후류·3-D handoff 개발과 검증 결과를 한곳에 고정한
재현 기록이다. 실행 산출물 전체는 저장소에 넣지 않고, 소스·CLI·검증 문서와
정량 결과만 공개한다.

## 이번 기간에 반영된 내용

- 장애물 주변 dense-core 선택적 보존형 meander closure와 파라미터 스윕
- 풍향·장애물 높이별 wall-flux·수직 후류·원거리 entrainment 민감도 후보
- 보존형 3-D inner/outer handoff 및 nested 145-sensor scoring 경로
- `ground_linear` 수직 관측 연산자(저고도 센서 보간)
- 600 s 후보별 농도 dose MAPE, 실행시간, 질량·종·운동량 보존성 및 승격 게이트
- 실패 후보를 자동 승격하지 않고 거절 사유와 재현 명령을 문서화

## 정량 결과

| 경로 | 600 s dose MAPE | 실행시간 | 판정 |
|---|---:|---:|---|
| 최선의 빠른 hybrid (`core-meander-s005-l150-w50`) | **26.14946%** | 약 **33.65 s** | 10% 목표 미달, 보존성 통과 |
| 3-D `ground_linear` nested 145 sensors | **84.45845%** | 약 **3044.64 s** (50.7분) | 정량 검증 보류 |

최선 hybrid의 종 보존 오차는 약 `3.66e-16`이며, 센서별 dose 보정계수는 사용하지
않았다. 3-D run은 4,801개 프레임(0–600 s)을 완료했고, 관측 양성 84개 중
0.1 vol% 기준 예측 검출은 53개였다. normalized dose bias는 `-0.88678`이다.

## 3-D 수치·보존성 게이트

- species mass relative residual: `4.15e-15`
- global mass residual: `4.57e-9 kg/s`
- continuity residual: `3.08e-12 kg/m3/s`
- momentum residual: `1.92e-8 kg m/s`
- thermodynamic corrector non-converged steps: `0`
- solver 내부 게이트: 모두 통과
- overlap coupling gate: 실패(peak/dose 변화가 10% 이내 조건을 충족하지 못함)

따라서 3-D 결과는 보존성은 입증하지만 TI021 센서 정량 검증을 통과한 모델로
승격하지 않았다.

## 현재 확인된 핵심 원인 후보

handoff HDF5에서 양의 species flux를 시간 적분하면 약 **3292.60 kg**으로,
workbook 질량 proxy 4962 kg의 약 **66.3%**만 handoff plane을 통과한다. 이 값은
질량이 inner domain에 잔류하거나 outer domain으로 전달되지 않는 현상을 가리키는
진단 신호이지, 단독으로 원인을 확정하는 증거는 아니다. 다음 작업은 짧은 보존성
진단으로 명시적 finite-source mass, 초기 방출 운동량, 근거리 dense-gas transport와
inner-to-outer continuation의 기여를 분리하는 것이다.

## 재현 명령

```powershell
.\.venv\Scripts\python.exe tools\run_thorney_ti021_inner_transient.py `
  --workbook ti21_input_ascii.xls --duration-s 600 --cell-size-m 5 `
  --inner-domain-m 450 640 40 --time-step-s 0.125 `
  --receptor-vertical-sampling ground_linear `
  --output-directory outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18

.\.venv\Scripts\python.exe tools\score_thorney_ti021_nested_145_sensors.py `
  --workbook ti21_input_ascii.xls `
  --inner-run-directory outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18 `
  --outer-run-directory outputs/thorney-ti021-wide-outer-domain-600s-x450y640-centered-2026-10-13-v1/outer `
  --output-directory outputs/thorney-ti021-groundlinear-nested-score-600s-2026-10-19
```

## 사용 범위

현재 공개본은 연구용 시나리오 비교, 물리 민감도 분석, 재현 가능한 개발 검증에
사용할 수 있다. 10% 이하 정량 검증, 설비 인허가, 안전거리의 단독 근거로는 아직
사용할 수 없다. `outputs/` 실행 산출물과 개인 임시 로그는 저장소 정책상 제외한다.

상세 수치는 다음 문서에 기록되어 있다.

- [`TI021_GROUND_LINEAR_3D_SCORE_2026-10-19.md`](TI021_GROUND_LINEAR_3D_SCORE_2026-10-19.md)
- [`TI021_3D_UNDERPREDICTION_AUDIT_2026-10-18.md`](TI021_3D_UNDERPREDICTION_AUDIT_2026-10-18.md)
- [`TI021_MAPE_REDUCTION_SESSION_PLAN_2026-10-17.md`](TI021_MAPE_REDUCTION_SESSION_PLAN_2026-10-17.md)
