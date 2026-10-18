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

- 현재 최선의 개별 후보: dose MAPE **26.15868%**, 실행시간 **36.88 s**
- 10% 목표: **미달성**
- 종 보존: 통과(최선 후보의 `unmapped_species_fraction` 약 `3.66e-16`)
- 3-D named packet 결합 후보: x=100/160/220 m에서 각각 44.3493/37.7435/35.1191%
  로 악화되어 모두 거절
- wide-inner packet의 임의 schema 변환: 50.4307% 및 56.2863%로 악화되어 거절
- 최신 전체 테스트: **2,430 passed, 145 skipped**

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

