# DEGALI 현장 목표 충족 감사 (2026-10-07)

이 문서는 현재 작업트리의 현장 반해상 LH₂ 확장 목표를 요구사항별로
대조한 감사 기록이다. `pass`는 코드와 회귀가 해당 계약을 검증한다는
뜻이고, `conditional`/`withheld`는 모델의 적용범위나 외부 증거가 아직
제한된다는 뜻이다. 이 문서는 설계기준 승인이나 현장 validation 승격이
아니다.

현재 작업트리의 최신 전체 저장소 회귀는 **1,659 passed, 145 skipped,
0 failures**이며, 아래 field 표적 수치와 별도로 자유장·검증·현장 계층을
함께 포함한다.

## 요구사항 대조

| 요구사항 | 현재 구현·권위 증거 | 현재 판정 |
|---|---|---|
| 기존 자유장 물리 보존 | 기존 자유장·LH₂ 핵심 subset 회귀: **247 passed, 139 skipped, 0 failures**. 전체 저장소의 기존 baseline은 별도 보존 | 검증됨 (선택적 skip 존재) |
| 범용 누출원 | `field_lh2.py`의 explicit LH₂ flash 및 pressure-driven throat adapter, `field_source_io.py`의 fingerprinted atmospheric schedule, measured-history/blowdown/pressure-driven-history CSV/CLI/operational handoff, phase-routing handoff | 코드·field 표적 회귀 통과; 외부 현장 source는 미확보 |
| 불확실성 전파 | source/weather/surface/sensor/ambient/location/direction/duration/obstacle/distributed schedule deterministic corner envelope; pressure-driven opening-area/Cd corners are recomputed through the throat before flash, and measured/pressure-history source location/direction corners are fixed before field projection | **490 passed**; 확률 해석을 하지 않는 deterministic bound |
| 모델형 비교 | `field_comparison.py`, strict CSV SHA/row/operator contract, Test4 v3 replay, FFI Test 6 time-history observation replay | 재계산·무결성 검증 통과; native source/time operator 불일치로 `conditional`, model-selection `withheld`; owner-exported history 없이는 replay artifact 미생성 |
| 보존형 장애물 해석 | `semi_fv_obstacle.py` local x-z finite-volume routing, source별 mass ledger, final inventory identity, obstacle geometry projection | **490 passed**; reduced-order obstacle closure이며 CFD/wake validation 아님 |
| 적용범위 자동 판정 | `FieldApplicability`, operational decision, gate code, refinement/uncertainty/conditional-review gates | blocked/conditional/allowed를 fail-closed로 분리 |
| 실패안전 출력 | `field-verify`, `field-audit-verify`, nested provenance/inventory/decision checks, exclusive-create artifact writers | 변조·재계산 불일치 회귀 통과; 최신 v14 감사 3종과 Test 4 comparison execution v3를 실제 CLI로 재검증(exit 0) |
| 외부 채널 결합 경계 | `FieldEvidenceManifest`의 명시적 5채널 선택·SHA·event/clock/operator ID와 파일 재검증, strict validation case linkage | manifest/case drift 차단; `promotion_allowed=false` 유지 |
| 외부 현장 validation 승격 | 최신 SLABx 감사 v14 (`outputs/current-slabx-audit-2026-10-07-v14.json`; 로컬 생성 산출물) | **partial**: 16,139개 구조화 파일 중 source-boundary 후보 15개와 수용기 후보 33개를 노출하지만 weather/obstacle/common-clock가 없음. source-rate/wind 및 운전 historian near-miss는 유형별 집계와 `collection_requirements`로 다음 수집 작업만 안내 |

## 외부 evidence 판정

최신 SLABx 감사 v14은 `promotion_allowed=false`를 유지한다. 명시적
`source_bound`/기준풍속/평균농도 alias가 후보 triage에 포함됐지만, 이는
event join이나 validation 승격이 아니다. 후보 메모의 `diagnostic_counts`는
고정 수용기 농도 누락 3,268개, weather identity 누락
638개, source identity 누락 262개, 운전 historian 누락 560개, 센서 clock
누락 3개를 별도 집계한다. 각 채널의 `collection_requirements`는 후보가
발견된 경우에도 수동 event/clock/geometry/calibration qualification을
요구하며, 누락 채널에는 다음 수집 조치를 기록한다.
검사센터 폴더는
20개 구조화 파일 중 17개가 읽기 가능한 workbook candidate로 inventory되고
14개는 모호한 보조 sheet 진단을 보존하지만, 필수 채널이 없어 `withheld`다.
관련 자료 폴더는 구조화 파일이 없어 `withheld`다. 수용기 후보가 발견됐다는 사실만으로
source boundary, 기상 동시성, 장애물 geometry, common clock를 추론하지
않는다.

검사센터 및 관련 자료의 v14 감사 artifact는 로컬 생성 산출물인
`outputs/current-support-center-audit-2026-10-07-v14.json`과
`outputs/current-related-audit-2026-10-07-v14.json`이다. 이 파일들은
재생성 가능한 감사 결과이므로 공개 스냅샷에는 포함하지 않는다.
별도로 확인한 2톤 탱크 GA/P&ID와 TK-1101/TK-1102 계기 도면, 설명자료는
정적 geometry·벤트 topology·6.5→6.0 bar 운전 시나리오만 제공하며, 특정
사건의 source-rate·동시 기상·수용기 농도·common clock를 제공하지 않는다.
따라서 이 문서들은 향후 geometry/source export의 준비 자료로만 남기고
validation 채널로 승격하지 않는다.
Test4 matched-sensor 비교의 최신 실행 artifact는
`../outputs/applied-energy-evidence-2026-10-03/test4-external/field-model-comparison-execution-v3.json`이다.
`field-verify` 재계산은 통과했지만 native source boundary와 temporal
operator가 달라 수치 비교는 `conditional`, 모델 선택 영향은 `withheld`다.

## 다음 승격 조건

1. 동일 사건에 대한 명시적 source boundary와 phase/source-rate 근거
2. 동시 측정된 wind speed/direction/stability와 common clock
3. obstacle geometry/coordinate reference 및 적용 representation
4. 고정 수용기 H₂ 농도·위치·높이·교정·응답시간·평균연산자
5. 위 자료를 같은 사건/시계로 묶는 외부 manifest와 SHA-256

위 조건이 충족되기 전까지 현장 경로는 screening/sensitivity 도구로만
사용하며, 설계기준·인허가·모델 우열 결론으로 승격하지 않는다.
