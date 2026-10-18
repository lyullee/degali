# DEGALI 업데이트 요약 — 2026-10-18

이번 업데이트는 TI021 장애물·후류·3-D handoff 검증을 포함한 현재 개발 상태를
GitHub 공개 저장소에 고정한 기록이다. 실행 산출물 전체를 저장소에 넣지 않고,
재현 가능한 소스·CLI·검증 문서만 포함한다.

## 포함된 변경

- 장애물 주변 dense-core의 선택적 보존형 meander closure
  (`core_meander_strength`, `core_meander_length_m`, `core_meander_width_m`)
- 600 s 전체 이벤트를 대상으로 한 TI021 후보별 MAPE·실행시간·보존성·높이/거리
  그룹 판정 기록
- 보존형 3-D handoff 인터페이스 점검 및 named packet/wide-inner packet 변환 감사
- 실패한 후보를 명시적으로 거절한 이유와 다음 물리 모델 작업 항목
- focused test 및 repository-wide test 결과 기록

## 현재 정량 상태

- 현재 최선의 개별 후보: `ti021-core-meander-s005-l150-w50-600s-2026-10-18.json`,
  dose MAPE **26.14946%**, 실행시간 약 **33.65 s**
- 10% 목표: **미달성**
- 종 보존: 통과(최선 후보의 `unmapped_species_fraction` 약 `3.66e-16`)
- 3-D named packet 결합 후보: x=100/160/220 m에서 각각 44.3493/37.7435/35.1191%
  로 악화되어 모두 거절
- wide-inner packet의 임의 schema 변환: 50.4307% 및 56.2863%로 악화되어 거절
- 강화 meander 후보(`s020/l250/w100`)는 **26.24680%**로 악화되어 거절했고,
  풍향 3° 후보도 **26.31148%**로 악화되어 거절했다. 센서별 dose 보정계수는
  사용하지 않았다.
- 풍향 9° 후보도 **26.23179%**로 기준보다 악화되어 거절했다.
- 장애물 높이별 wall-flux 분기(`height_resolved_wall_flux=True`)도 **26.22907%**로
  기준보다 악화되어 거절했다. 보존성은 동일하게 통과했다.
- 외부 원거리 entrainment을 제거한 후보(`exterior_entrainment_far_factor=0`)는
  **31.45541%**로 악화되어, 단순 원거리 희석이 주된 잔차 원인이 아님을 확인했다.
- 후속 600 s 물리 후보도 기록했다: 수직 diffusive-tail 확산계수 `0.20 m²/s`
  는 **30.98524%**, geometry-only 후류 상부분기(`wake_vertical_elevated_fraction=0.10`)
  는 **26.26299%**로 기준보다 악화되어 모두 거절했다. 두 후보 모두
  `unmapped_species_fraction ≈ 3.66e-16`으로 보존성은 통과했지만 정량 게이트는
  통과하지 못했다.
- 최신 전체 테스트 기준: **2,430 passed, 145 skipped**

## 최신 3-D nested 경로 감사

보존형 wide-inner 3-D 경로를 별도로 점검한 결과, 보존성은 통과하지만 센서
dose를 크게 과소예측했다. inner-only 양의 관측 dose MAPE는 **85.5286%**,
예측/관측 총 dose 비율은 약 **0.1033**이며, 기존 nested score MAPE는
**84.9011%**였다. 이는 10% 목표 후보가 아니므로 승격하지 않았다.

현재 우선순위는 동일 조건의 `ground_linear` 600 s 전체 재계산을 완료한 뒤
외부 센서 score를 재계산하고, 개선되지 않으면 해당 후보를 거절하는 것이다.
현재 실행은 `outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18/`에
체크포인트를 기록 중이며, 완료 전에는 결과를 승격하지 않는다. 상세 감사는
[`TI021_3D_UNDERPREDICTION_AUDIT_2026-10-18.md`](TI021_3D_UNDERPREDICTION_AUDIT_2026-10-18.md)
에 기록했다.

## 해석 및 사용 범위

현재 공개본은 연구용 시나리오 비교와 물리 모델 민감도 분석에 사용할 수 있다.
10% 이하 정량 검증과 설비 인허가·안전거리 단독 근거로의 승격은 아직 허용하지
않는다. 센서별 보정계수나 숨은 질량을 사용한 후보는 채택하지 않았다.

상세한 후보별 수치와 재현 경로는
[`TI021_MAPE_REDUCTION_SESSION_PLAN_2026-10-17.md`](TI021_MAPE_REDUCTION_SESSION_PLAN_2026-10-17.md)
에 기록되어 있다.

## 재현에 필요한 핵심 파일

- `src/degali/addons/dense_layer_patch.py`
- `tools/score_ti021_dense_layer_hybrid.py`
- `tools/run_thorney_ti021_outer_transient_checkpointed.py`
- `src/degali/addons/conservative_nested_handoff_3d.py`
- `src/degali/addons/binary_mixture_handoff_plane_3d.py`

실행 산출물(`outputs/`)과 개인별 임시 로그는 저장소 정책에 따라 제외한다.

