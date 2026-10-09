# DEGALI 현재 인수인계 문서

최종 갱신: 2026-10-09
프로젝트: **DEGALI — Dense Gas Dispersion for Liquid Hydrogen**
Git remote: `https://github.com/lyullee/degali.git`
현재 브랜치: `main`
마지막 공개 태그(작업 시작 기준): `v0.2.0` (`b0e3c34`)
현재 릴리스 후보: `v0.3.0` (Zenodo DOI `10.5281/zenodo.23256538`)

이 문서는 다른 계정/세션에서 작업을 이어가기 위한 최신 기준점이다. 먼저 이
파일을 읽고, 그 다음 [docs/model-hardening-2026-09-21.md](docs/model-hardening-2026-09-21.md),
[docs/ffi-test6-decomposition.md](docs/ffi-test6-decomposition.md),
[docs/HANDOVER.md](docs/HANDOVER.md) 순서로 읽는다.

## 0B. 2026-10-05 현장 반해상 배포 경로

- 현장 증거 수집의 첫 입력 경계는 `FieldEvidenceManifest`다. 기존 5개
  채널(source boundary, weather, obstacle geometry, receptor observations,
  common clock)의 선택 파일과 SHA-256은 유지하면서, 이제
  `sensor_set_id`/`operator_id`를 기록한다. 센서 registry·교정 파일은
  `sensor_registry_artifact`로 별도 SHA-256 고정할 수 있다. 두 ID 또는
  registry가 없으면 manifest 자체는 추적 가능하되 `evidence_readiness`가
  `conditional`이며, 모델 qualification이나 설계기준 승격은 허용하지
  않는다. `docs/field-evidence-intake-template.md`의 패키지 구조와
  `--sensor-registry-path` CLI 옵션을 사용한다. 채널 누락은 기존처럼
  `field-audit`의 `partial`/`withheld`, 해시 변조는 hard failure이다.
- 새 현장 경로는 `src/degali/addons/field_*.py` 및
  `semi_fv_obstacle.py`에 있다. 이는 기존 자유장 DEGADIS 호환 경로를
  대체하지 않는, 국소 2-D wind-plane 보존형 scalar transport 확장이다.
  ReleaseSource, 기상, 표면, 센서, 장애물, direct-vapour/풀/액적 handoff와
  각 적용범위가 명시적으로 기록된다.
- `FieldSemiFVRequest`는 선택적인
  `ambient_temperature_uncertainty_k`, `ambient_pressure_uncertainty_pa`,
  `ambient_air_density_uncertainty_kg_m3` bounded boundary를 받는다. 각
  ambient corner는 flash/transport/sensor 변환까지 실제로 실행되고,
  report selection/provenance에 남는다. phase-routing/pool handoff도 같은
  corner를 phase ledger와 field request에 전파한다. nominal ambient bound를
  그대로 둔 operational screening은 `withheld`이며, measured-history joint
  adapter는 시간정렬 이력과 ambient corner가 모두 전파될 때까지 fail-safe
  보류한다. ambient pressure bound가 선언되면 reusable LH2 saturation table도
  양 끝점의 saturation temperature를 포함해 nominal-only table domain으로
  pressure corner가 거부되지 않도록 한다.
- `FieldValidationEvidence`는 외부 dataset identity/공통시각/기하를
  fingerprint하는 경계일 뿐 실제 observed-vs-model score가 아니다. 따라서
  evidence만 붙은 field applicability는 이제 `conditional`로 남고, 별도
  matched validation score가 생성되기 전에는 모델 qualification으로
  해석되지 않는다.
- field LH2 flash handoff는 phase mass partition, momentum 및
  total-specific-energy residual을 report에 보존한다. 비유한값 또는
  tolerance 초과 residual은 semi-FV 주입 전에 `blocked`가 되며,
  direct scalar branch가 잔여 post-flash liquid를 조용히 기화시키지 않는다.
  measured-history 및 cryogenic-blowdown schedule adapter도 동일한 closure
  gate를 적용한 뒤에만 시간변화 atmospheric source를 방출한다. report에는
  mass/phase/momentum/energy gate tolerance도 함께 남는다. gate 실패 시에도
  JSON-safe residual 진단은 보존하지만 unclosed flash plane은 transport에
  노출하지 않는다.
- semi-FV transport 진단은 전체 주입량과 별도로 `primary` 및 각
  `distributed:<label>` source schedule의 실제 주입 질량 원장을 JSON으로
  보존한다. 따라서 다중 source case에서 총 질량 보존뿐 아니라 source별
  bookkeeping도 재현 가능하다.
  typed diagnostics는 비유한/음수 질량·Courant, [0, 1] 밖 diverted fraction,
  중복 source-ledger label과 알 수 없는 applicability를 report 생성 전에
  거부한다. obstacle rectangle도 local inlet/outlet에 닿거나 domain 밖으로
  돌출되면 mesh clipping 없이 fail-safe로 차단한다.
- 다중 cuboid 장애물은 풍향별 local plane에 투영된다. off-plane source,
  sensor 또는 geometry와 검증되지 않은 lateral wake/3-D 난류는 계산값으로
  꾸미지 않고 conditional/blocked로 남긴다. 단일 및 복수 장애물은 서로
  섞을 수 없다.
- 계측 P/T/유량 이력은 alignment/응답시간/clock offset/interval 오차의
  evidence-based quality gate를 통과해야 transient source로 승격된다.
  `read_measured_history_csv()`는 명시적 SI column map, 채널별 교정/응답
  evidence, 별도 phase/flash-state evidence와 input SHA-256을 요구하는 유일한
  historian import boundary다.
  paired lower/upper calibration column도 각 행에서
  `lower <= nominal <= upper` 및 `lower <= upper`를 만족해야 하며, 뒤집힌
  bound나 nominal을 제외하는 bound는 flash/uncertainty 전파 전에 fail-safe로
  거부한다. 정렬·clipping으로 교정하지 않는다.
  event·채널 provenance는 field report에도 보존한다.
  긴 이력은 온도 범위 saturation table을 자동 사용하고, 짧은 이력은 table
  생성비용 때문에 direct reference path를 유지한다.
- 개발 PC에서 `tools/benchmark_lh2_flash_table.py --count 500 --nodes 161`은
  table build 0.051 s, direct/table 0.263/0.054 s (약 4.87x), quality 상대차
  6.2e-6을 기록했다. 이는 환경 의존 performance 진단이며 정확도 검증을
  대체하지 않는다.
  pressure-only workbook이나 nominal flow는 여기서 source로 추론되지 않는다.
  `run_field_measured_history_envelope()`는 source-history와 declared source
  location/direction corner를,
  `run_field_joint_measured_history_envelope()`는 history 재-flash와 source
  location/direction, 풍향·풍속·모든 in-plane detector calibration corner를
  함께 계산한다. detector corner는
  같은 true transport trace를 재사용하므로 transport를 다시 풀지 않는다. 이미
  measured flow로 대체된 static P/T/유량/opening/Cd는 두 번 변동시키지 않는다.
  measured-history envelope는 이력 온도 모든 lower/upper 범위를 덮는 LH2
  saturation table 하나를 자동 재사용하며 report에 사용 여부를 남긴다.
  direct vapour에서 unrouted liquid의
  표면 민감도는 주장하지 않으며 phase-routing/pool branch가 필요하다.
- 압력/온도만 있는 명시적 이력에는 `PressureDrivenMeasuredHistory` 경로를
  사용한다. 이 경로는 absolute pressure, declared opening area/Cd, ambient
  pressure와 source ID를 사용해 interval별 throat flow를 재계산한 뒤 flash하며,
  압력 추세를 누출률로 해석하지 않는다. `direct_vapour_schedule_envelope_from_pressure_driven_history()`와
  `run_field_pressure_driven_history_envelope()`는 P/T/liquid-fraction bound를
  coherent corner로 묶어 semi-FV까지 전파한다. ambient temperature/pressure,
  weather, stability, detector corner까지 필요하면
  `run_field_joint_pressure_driven_history_envelope()`를 사용해 throat 재계산
  전후를 함께 보존한다. source-only 경로는 unresolved ambient-pressure bound를
  별도 joint 경로 없이 fail-safe 거부하며, geometry/Cd uncertainty도 joint source
  envelope에서 throat 재계산과 함께 전파한다. source location/direction bound도
  각 field-history corner의 wind-plane projection 전에 고정한다. typed object만으로 외부 historian SHA/event/clock/
  calibration evidence가 충족되었다고 보지 않는다. operational refinement와
  aggregate gate까지 필요하면 `run_field_operational_joint_pressure_driven_history_envelope()`를
  사용하며, 하나라도 withheld corner이면 aggregate도 withheld다. 이 경로 추가 후 field/
  semi-FV 표적은 **490 passed**다.
- strict pressure-driven history CSV/CLI/replay와 operational wrapper를 포함한
  전체 저장소 회귀는 **1,659 passed, 145 skipped, 0 failures**다. skip은
  선택적 CoolProp/reference 환경 의존 항목이며, 현장 외부 evidence 승격을
  의미하지 않는다.
- 최신 v14 SLABx/support-center/related audit artifact 3종과 Test 4
  comparison execution v3를 `field-audit-verify`/`field-verify`로 다시 실행해
  모두 exit 0을 확인했다. 이는 파일·provenance 무결성 통과이며, 세 감사의
  `promotion_allowed=false`와 comparison의 현장 승격 보류를 변경하지 않는다.
- 연결된 PDF/PPTX도 별도 검토했다. 2톤 탱크 GA/P&ID와 TK-1101/TK-1102
  계기 도면은 정적 치수·노즐·벤트 topology를, 시나리오 설명자료는
  6.5→6.0 bar 자동 BOG 벤트 운전을 보여준다. 그러나 사건별 atmospheric
  source-rate, 동시 weather, 고정 수용기 H2 관측·교정, common clock가 없어
  `FieldEvidenceManifest` 채널로 승격하지 않고 context-only로 유지한다.
- integrity-only `field-verify`도 이제 serialized field applicability와
  operational allowance의 상태를 교차검증한다. 단일 screening과 nested
  envelope case 모두 `blocked` 물리 결과를 허용 decision으로 감싼 변조를
  거부하며, 이 보강을 포함한 field 표적은 **490 passed**다.
- common-clock field-evidence triage는 명시적 `synchronization_id`, `sync_id`,
  `time_sync_id`, `clock_sync_id` alias와 timestamp에 붙은 synchronization
  metadata를 인식하지만, 빈 값은 계속 증거로 세지 않는다. 이 보강은
  collection plan에도 새 explicit alias를 반영하며, 기존 v1 historical alias
  plan도 읽을 수 있게 유지한다. `promotion_allowed=false` 경계를 바꾸지 않는다.
- `field_*_report`와 batch manifest는 입력, quality evidence, 질량분할,
  numerical refinement, applicability, withheld sensor를 JSON으로 남긴다.
  현재 모델은 연구/운영 screening 보조용이고 독립 LH2 obstacle validation
  전에는 설계기준·인허가 단독 판단으로 승격하지 않는다.
- `evaluate_field_operational_screening()`은 screening report와 별도로
  transport 완료, in-plane sensor, grid/time refinement와 conditional review
  opt-in을 검사한다. source/weather/sensor의 남은 bounded uncertainty를
  propagation하지 않은 nominal run도 refinement 후 screening을 보류한다.
  batch manifest v3는 각 case의 이 decision, 필요 조치와 conditional review
  authorization을 기록한다. 어떤
  경우에도 design basis/approval을 허용하지 않는다.
- typed `gate_codes`는 이제 개별 decision뿐 아니라 uncertainty envelope,
  measured-history, sensor-array, joint source-sensor, phase-routing aggregate의
  최상위 decision에도 corner별 union으로 전파된다. aggregate conditional
  holdback은 `conditional_review_required`를 보존하고, legacy/manual withheld
  corner만 남은 경우에는 `corner_withheld`를 사용한다. 따라서 자동화가
  human-readable reasons를 파싱하지 않고도 최상위 보류 원인을 분기할 수 있다.
- fixed-sensor model comparison도 `gate_codes`를 내보낸다. 센서 집합/형상
  불일치와 execution block은 차단 코드로, basis·평균연산자·시간모드
  불일치와 conditional/unresolved execution은 조건부 코드로 구분한다.
  matched basis는 `matched_basis`로 표시하고, decision impact에는
  `classification_aligned`/`classification_disagreement` 및 필요 시
  `model_selection_withheld`를 추가한다. JSON report의 human-readable
  reasons는 그대로 유지되며, 이 코드는 batch/검증 자동화의 안정 키다.
  threshold는 물리적 mole-fraction 범위 `[0, 1]`을 벗어나면 case/function/
  typed result 생성 단계에서 거부해, 모든 수용기를 잘못된 below-limit으로
  분류하는 경로를 차단한다.
- operational uncertainty/sensor/history/phase/source-sensor envelope case도
  이제 attached physical applicability와 decision status를 교차검증한다.
  `screening_allowed`는 `accepted`, `conditional_allowed`는 `conditional`만
  허용하며, blocked 또는 상태가 어긋난 물리 결과를 수동 decision으로
  감싸 aggregate에 넣는 경계를 차단한다.
- matched field-validation score는 comparison의 gate code를 잃지 않고
  그대로 전파한다. 여기에 obstacle evidence 미지원과 lower-bound
  관측의 존재/위반/미해결 상태를 별도 코드로 붙이며, 실제 보류가 없는
  경우에만 `validation_qualified`를 기록한다. 따라서 validation/batch
  소비자는 점수의 human-readable reason을 재파싱할 필요가 없다.
- read-only field-evidence audit도 `gate_codes`를 내보낸다. 채널 누락/부분
  탐지/scan limit, coherent 또는 cross-file candidate 형태를 구분하고,
  inventory만으로는 절대 승격하지 않는 `promotion_not_allowed`를 항상
  보존한다. 구형 audit JSON은 optional code 없이도 계속 읽힌다. 추가로
  저장된 `channel_paths`, coherent 후보, status/scan completeness를 실제
  candidate의 `detected_channels`와 재구성해 비교하므로, 파일과 SHA가
  그대로여도 채널 인덱스만 바꾼 audit artifact는 `field-audit-verify`에서
  보류된다.
- phase/pool operational transport에서 실행 불능 corner는 일반
  `transport_incomplete`와 함께 `phase_routing_incomplete` 또는
  `pool_transport_withheld`를 보존한다. aggregate gate union이 ledger
  실패와 pool-to-field handoff 실패를 구분할 수 있다.
- semi-FV 진단과 operational decision 사이에도 원인 코드가 연결된다.
  질량 잔차 초과는 `transport_mass_residual_exceeded`, source별 ledger
  불일치는 `source_ledger_mismatch`, 장애물 포함 reduced-order 범위는
  `obstacle_transport_conditional`로 보존되며 상세 residual/warning은
  그대로 남는다.
- `SourceRateSchedule`은 `has_zero_endpoint`와
  `require_zero_endpoint()`를 공통으로 제공한다. 저수준 기록은 비영
  retained endpoint를 검사할 수 있지만, 실제 semi-FV primary/distributed
  주입은 명시적 zero endpoint 없이는 시작되지 않는다. strict CSV와 field
  request도 같은 검사를 사용해 종료 시각 이후의 묵시적 source 연장을 막는다.
- semi-FV diagnostics에는 이제 source별 선언 schedule 적분질량과 실제
  solver 주입량의 독립 residual이 들어간다. `source_mass_schedule_residual_kg`
  와 최대 residual은 report에 남고, 불일치 시
  `source_schedule_mass_mismatch` gate code로 operational screening을
  fail-closed한다. 기존 aggregate ledger residual만으로는 내부적으로
  일관된 잘못된 schedule bookkeeping을 잡지 못하던 경계를 보완한다.
  이 보강 뒤 field/semi-FV 전체 표적은 **398 passed**이며 `compileall`과
  `git diff --check`도 통과했다.
- measured-history와 pool-vapour schedule wrapper도 typed
  `SourceRateSchedule`, zero endpoint, warning/alignment metadata를 생성
  단계에서 검증한다. measured-flash direct-vapour mass와 pool evaporation
  ledger mass가 schedule 적분과 다르면 field source boundary로 승격하지
  않는다. 관련 history/pool 회귀는 **22 passed**다.
- joint measured-history envelope가 ambient temperature/pressure bound를
  history 재-flash에, ambient air-density bound를 detector 변환에 실제로
  전파한다. 각 field case는 평가한 ambient 값으로 uncertainty를
  해소하고, aggregate report selection에 ambient key를 보존한다. 관련
  field 전체 표적은 **401 passed**다.
- `field-verify`는 measured-history nominal/refined 및 joint-history envelope를
  현재 historian CSV에서 다시 계산하고, phase-routing·standalone
  sensor-array·joint source-sensor envelope도 저장된 옵션으로 재계산한다.
  source-sensor의 optional atmospheric schedule은 경로·SHA·schedule record를
  다시 확인한다. 여전히 replay할 수 없는 option-incomplete envelope와
  batch report 안의 nested field transport diagnostics는 semantic하게
  검사한다. 일반 bounded `nominal_field` uncertainty envelope도 저장된
  corner option으로 재계산한다. source schedule residual의 finite/unique
  형식과 per-source
  maximum 일치성을 확인하고, 수치 gate 초과가 있으면 동일 원인의 typed
  blocked applicability reason이 없을 때 fail-closed한다. 이 보강을 포함한
  execution 표적은 **6 passed**다.
- 같은 integrity-only 경계에서 compact transport concentration과 sensor
  extrema/time/density도 검사한다. withheld detector는 warning과 `result=null`
  을 함께 가져야 하며, true mole-fraction 범위·시간순서·ambient density가
  변조되면 `field-verify`가 통과시키지 않는다.
- strict pressure-derived source를 사용하는 경우에는 재실행 불가능한
  report와 nested sensor/envelope report에서도 typed ambient-pressure
  boundary, source ID, uncertainty record와 derivation metadata가 입력 case와
  일치하는지 추가로 확인한다.
- `field-verify`는 serialized operational decision의 status/allowance,
  reasons/actions, unique typed `gate_codes`와 approval 금지도 재검증한다.
  따라서 report와 decision을 각각 그럴듯하게 변조해도 integrity-only
  artifact가 허용 상태로 통과하지 않는다.
- transport diagnostics의 `source_mass_injected_kg` label 집합·합계도
  ledger residual 및 schedule residual label 집합과 교차검증한다. 따라서
  source별 원장과 총 주입량을 따로 바꾼 artifact도 fail-closed한다.
- 같은 integrity-only 경계에서 report의 `transport_input`에 기록된 direct
  vapour/distributed source schedule 질량도 `primary` 및
  `distributed:<label>` injected ledger와 대조한다. 선언 질량·label 집합이
  다르면 수치 재실행 없이도 verifier가 보류하여, 총 잔차가 작아 보이는
  source provenance 변조를 통과시키지 않는다. 이 보강을 포함한 field
  execution 표적은 **7 passed**다.
- 추가로 `mass_injected_kg = mass_domain_kg + mass_outflow_kg` 최종 inventory
  identity를 동일한 절대/상대 질량 gate로 재계산한다. 최대 잔차 필드만
  0으로 바꿔 domain/outflow 총량을 변조하거나 음수 질량을 넣은 execution
  artifact도 integrity-only 경계에서 보류된다. 관련 execution 표적은
  **8 passed**다.
- semi-FV solver 자체도 `final_mass_residual_kg`를 typed diagnostics에
  기록하고 동일한 최종 inventory gate를 직접 field screening에 적용한다.
  따라서 실행 JSON을 만들기 전의 Python 경로에서도 domain/outflow 폐쇄
  실패가 `blocked`로 남는다.
- `FieldAtmosphericSourceSchedule`/`field-source`는 이제 provenance가 있는
  `declared_atmospheric_vapour`도 primary atmospheric boundary로 받을 수
  있다. 이는 scenario의 명시적 위치·flash contract를 유지한 채 source
  schedule만 대체하며, 위치 없는 pool/droplet schedule은 여전히
  distributed/phase-routing handoff 없이는 거부한다.
- 이 변경 뒤 field/semi-FV 표적은 **452 passed**다(감사 artifact binding,
  다국어 헤더 inventory 및
  observed-volume-percent readiness 회귀 포함). pressure-driven adapter
  직전의 마지막 완전한 저장소 회귀 기준은 **1,601 passed, 145 skipped,
  0 failures**였고, 이번 adapter 자체는 `test_field_lh2.py`의 21개
  focused 회귀와 strict-case parser 회귀로 확인했다.
  기존 자유장·FFI 경로와 새 현장 경로의 보존 상태는 표적 회귀에서
  확인됐다.
- 현재 자유장 보존 subset(axisymmetric jet, cryogenic-air/blowdown, LH2
  scope/yaw/rainout, reference parity, model-comparison, screening gate)은
  **247 passed, 139 skipped, 0 failures**였다. skip은 선택적 CoolProp 또는
  reference 환경 의존 case이며, 이 수치는 전체 저장소 회귀를 대체하지
  않는다.
- 직접 생성하는 `FieldSemiFVRequest`의 ambient uncertainty도 strict case와
  동일한 `K`/`Pa`/`kg/m3` 단위를 요구한다. 잘못 표시된 ambient corner가
  flash·transport·sensor 변환에 들어가기 전에 거부된다.
- 직접 생성하는 `ReleaseSource`/`WeatherState`/`SensorModel`/
  `SurfaceBoundary`도 비어 있지 않은 단위 라벨이 물리량의 기대 SI 단위와
  다르면 생성 단계에서 fail-closed 한다. 기존 `BoundedValue` 직접 호출의
  빈 라벨은 호환성을 위해 허용하고, strict JSON은 여전히 정확한 단위를
  요구한다.
- `SourceRateSchedule`은 기존 `piecewise_constant`를 기본으로 유지하면서
  명시적 `linear` time-node operator를 지원한다. 두 모드 모두 zero endpoint,
  clipped mass integration, source별 ledger gate를 공유하며, strict
  JSON/CSV provenance와 lower/nominal/upper corner에 operator를 보존한다.
- 기존 v1 serialized artifact에 operator가 없으면 verifier가 이를 명확한
  legacy `piecewise_constant`로 메모리상 보정한다. 명시된 operator 값은
  그대로 provenance 비교에 사용한다.
- `pressure_driven_lh2_mass_flow()`는 기존 HEM throat closure를 이용해
  upstream pressure/temperature, opening area, Cd와 선택적 ambient-pressure
  endpoint를 deterministic corner로 교차하여 `kg/s` bounded rate를 만든다.
  `release_with_pressure_driven_lh2_mass_flow()`는 이를 `ReleaseSource`에
  붙이고 derivation/source ID metadata를 남기지만, 이미 선언된 non-zero
  measured rate는 덮어쓰지 않는다. corner가 해석되지 않으면 fail-closed한다.
  strict field JSON도 `pressure_driven_mass_flow` object를 명시하면 같은
  adapter를 사용하고 ambient bound/source ID를 provenance에 남기며, rate와
  opt-in object를 모두 생략하거나 동시에 선언하는 모호한 입력은 거부한다.
- pressure-derived rate는 typed derivation boundary로 보존된다. 따라서
  uncertainty envelope가 이미 계산된 `mass_flow_kg_s`를 pressure/temperature/
  opening/Cd와 독립적으로 다시 곱하지 않고, 선택된 ambient-pressure corner에서
  throat rate를 재계산한다. `ReleaseSource.as_dict()`에도 이 구조와 source ID가
  남는다.
- matched field validation의 observed CSV는 선택적으로
  `observation_kind_column`을 선언할 수 있다. `exact`와 `lower_bound`만
  허용하며, 후자는 검출/검열 하한으로서 MAE/RMSE/bias에서 제외된다.
  하한 만족도와 deficit을 별도로 기록하고, 전부 만족해도 `conditional`,
  위반 또는 matched comparison 불완전 시 `withheld`로 fail-closed한다.
  score JSON에는 exact 행과 lower-bound 행의 개수, one-sided constraint,
  exact 행에만 제한한 classification disagreement도 함께 남긴다. 관련
  field-validation 회귀는 현재 `20 passed`다.
- strict case의 optional `conditional_review`는 review ID, reviewer ID/role,
  UTC timestamp, evidence ID와 정확한 `conditional_screening_only` scope를
  요구한다. `--allow-conditional`은 이 기록 없이는 transport 전에 거부되고,
  execution JSON은 적용 여부를 보존한다. `export_field_screening_batch()`도
  conditional opt-in이면 모든 named case별 authorization을 요구하며 누락 시
  파일 생성 전에 실패한다. 이 기록은 refinement, uncertainty, in-plane
  receptor, obstacle representation 등 비협상 gate를 우회하지 못한다.
- `field_case_io.py`와 `degali field-screen`은 strict JSON nominal field case
  format을 제공한다. value마다 단위와 evidence source가 필수이며, case SHA,
  numerical report, operational disposition을 새 출력 파일 하나에 남긴다.
  base v1은 임의 post-flash schedule이나 jet/phase inference를 허용하지
  않지만, 위치·수직폭·zero-endpoint schedule과 evidence가 명시된
  `distributed_vapour_sources`는 strict 경계에서 받는다. 해당 source-shape
  corner는 envelope로 전파되고 `field_semi_fv_envelope_report()`는
  declared source/weather/primary sensor corner 모두를 JSON으로 보존하지만
  자체로 operational decision은 아니다.
- Python `FieldSemiFVRequest`에 직접 붙이는 post-flash `SourceRateSchedule`도
  이제 기본 `source_id="declared"`를 거부하고 명시적 source provenance를
  요구한다. measured-history 및 pool launch adapter는 각자의 evidence와
  함께 explicit ID를 전달한다.
- `run_field_operational_uncertainty_envelope()` 및 `degali field-screen
  --uncertainty-envelope`는 모든 nominal source/weather/primary sensor
  corner의 refinement와 decision을 합친다. 모든 corner가 동일한 gate를
  통과해야만 aggregate screening을 허용하며 `--no-refinement`는 결과를
  diagnostic/withheld로 남긴다.
- `export_field_screening_batch()`는 `atmospheric_source_schedules`에
  case별로 fingerprint된 post-flash atmospheric schedule을 받으면 해당
  case를 complete operational uncertainty envelope로 실행하고 source
  provenance와 모든 deterministic corner를 case JSON에 기록한다. 이 경로는
  refinement·in-plane sensor·resolved uncertainty gate를 항상 요구한다.
  conditional opt-in이면 named case별 authorization을 manifest에 보존하고,
  그 authorization은 모든 deterministic corner를 통과한 aggregate에만
  `conditional_allowed`를 부여한다. authorization이 비협상 gate를 우회하지는
  않는다.
- `field_batch_io.py`와 `degali field-batch`는 named batch JSON 경계를
  제공한다. 각 strict field/source JSON의 절대경로와 SHA-256은
  `batch-execution.json`에 남고, historian/phase-routing case를 이 direct
  batch 경로에 섞는 입력은 거부한다. `--require-screening`은 보류 case가
  있으면 exit 2를 반환하면서 보고서와 provenance는 보존한다.
- transient `ReleaseSource`의 `duration_s`도 bounded object 또는 별도
  `duration_uncertainty`로 명시할 수 있다. nominal duration과 lower/upper
  timing corner가 실제 field envelope case로 전파되며, transport domain은
  upper duration을 덮어야 한다. timing probability나 source schedule은
  추론하지 않는다.
- `ReleaseSource.location_m`도 세 좌표를 각각 bounded `m` object로 선언할
  수 있다. 모든 Cartesian lower/upper 조합이 실제 envelope case가 되어
  wind-plane projection·장애물·receptor gate를 다시 계산한다. 숫자와
  bounded 좌표를 섞은 배열은 거부한다. `direction_m`도 세 개의 bounded
  dimensionless component를 받을 수 있고 zero-vector corner는 거부된다.
  이 direction corner는 jet/phase handoff에 보존되며, direct scalar closure가
  운동량을 해석하지 않는 경우 그 효과를 발명하지 않는다.
- `field-screen` case JSON은 선택한 SI historian CSV를 `measured_history`
  object로 붙일 수 있다. 이 경로도 event/phase evidence, channel calibration,
  t90/clock/interval criteria와 CSV SHA를 강제하고 report에 보존한다. 그러나
  non-exact source history는 nominal schedule만으로 operational 허용되지 않으며
  `run_field_operational_joint_measured_history_envelope()` 또는
  `field-screen --uncertainty-envelope`로 모든 source/weather/detector
  corner 및 refinement를 전파해야 한다. 사용 형식은
  `docs/field-measured-history-input.md`에 있다.
- `tools/audit_test4_slabx_degalix.py`는 conditional external SLABx Test 4
  comparison에 `field_model_decision_impact` 기록도 출력한다. 2026-10-05
  인접 SLABx 결과 재현에서 30개 sensor 중 LFL 분류 차이 3개가 보였지만,
  native source mapping과 275 s temporal operator가 불일치하므로 decision
  상태는 `withheld`였다. 이 결과는 모델 선택/위험거리 판정 근거가 아니다.
- `degali field-compare`와 `field_comparison_io.py`는 이 Test 4 전용 도구를
  대체하지 않으면서, 일반적인 DEGALI–SLABx/CFD 및 steady–transient fixed-sensor
  비교를 strict JSON case와 SHA-256 CSV provenance로 재현한다. source, weather,
  sensor geometry, temporal operator/averaging ID, 또는 obstacle mask/wake
  representation ID 중 하나라도 다르거나
  steady/transient mode가 다르면 numerical row는 남지만 model-selection impact는
  `withheld`이고 `--require-comparable`은 exit 2를 반환한다. matched case도
  monitored receptor 결과일 뿐 연속 위험거리 또는 정확도 검증 주장이 아니다.
  CSV 밖에 있는 native source/flash 또는 평균시간 근거는 optional
  `comparison_evidence`의 manifest path/SHA-256/qualification으로 고정한다.
  입력 parser는 manifest가 바뀌면 거부하고 실행 report에도 그 제한을 그대로
  남긴다. 현재 verifier로 재계산한 Test 4 실행 기록은
  `outputs/applied-energy-evidence-2026-10-03/test4-external/field-model-comparison-execution-v3.json`이며,
  `report_recomputed=true`, `promotion_allowed=false`, comparison은
  `conditional`, selection impact는 `withheld`로 유지된다.
  30개 receptor 중 LFL 분류 차이 3개, status `conditional`, selection impact
  `withheld`이며, 이 결과는 native source 및 time operator가 다르다는 사실을
  더 명확히 보존한 것이지 모델 우열이나 위험거리 결론이 아니다.
- 풍속·풍향 lower/upper와 달리 categorical stability는 `FieldStabilityAlternatives`
  및 strict JSON `stability_alternatives`로 명시한다. 모든 nominal+alternative
  class에 evidence-backed `StabilityScalarMixingClosure` diffusivity가 있을 때만
  `run_field_semi_fv_envelope()`와 historian joint envelope가 각 class를 실제로
  실행한다. nominal result는 대안이 미전파된 동안 withheld이고, complete
  operational envelope의 각 refined case에만 resolved 표시가 붙는다. 확률가중치,
  자동 stability correlation, RANS/LES/wake closure를 도입한 것은 아니다.
- strict field case의 `scenario.weather.direction_deg`는 이제 원형 구간이다.
  따라서 `lower=350`, `nominal=0`, `upper=10`은 북쪽을 가로지르는 20°
  풍향 불확실성으로 audit/endpoint envelope에 남고, 340° 반대 방향 구간으로
  해석되거나 거부되지 않는다. 이 기능은 static wind-plane sensitivity일 뿐
  풍향 사행을 평균내는 transient wind solver가 아니다.
- `run_field_phase_routing_transport_envelope()`는 phase/rainout/dynamic-pool
  corner를 실제 local semi-FV source로 연결한다. direct flash는 유지하고,
  보존·시간분해된 pool evaporation만 independent distributed vapour source로
  추가한다. in-flight droplet vapour, pool launch temperature/momentum/footprint와
  3-D wake는 채우지 않으며, pool ledger가 없거나 보존/시간 조건이 불충분하면 그
  corner의 transport를 `None`/withheld로 남긴다.
  `run_field_operational_phase_routing_transport_envelope()`는 이 physical corner와
  모든 detector calibration corner를 다시 refinement/decision해 한 건이라도
  보류되면 aggregate를 `withheld`로 둔다. pool source가 들어간 single corner는
  `phase_routing_uncertainty_resolved=False`로 남아 단독 operational use가 불가하며,
  complete aggregate의 refined case에서만 true가 된다.
- 전용 strict JSON `degali.field-phase-routing-screening-input.v1`과
  `docs/field-phase-routing-case.example.json`이 위 phase/pool 경로를
  `degali field-screen --uncertainty-envelope`에 연결한다. historian 및 임의
  subordinate-model option은 함께 넣을 수 없고, explicit ground-level scalar
  launch/evidence, evidence-backed droplet diameter/mass classes와 unit-labelled
  gas/droplet numerical-domain whitelist만 받는다. pool footprint와 d-squared
  evaporation 가정도 별도 phase-routing evidence ID를 요구한다.
  phase/pool scalar에 evidence-backed `nominal/lower/upper`를 선언하면
  `phase_routing.uncertainty`가 pool area, d-squared coefficient, time-window
  및 pool timestep corner를 phase envelope/transport까지 전파하고, upper
  post-release를 넘지 못하는 transport는 즉시 보류한다. droplet class
  분포 범위는 아직 추정하지 않고 exact 선언으로만 유지한다.
  pool launch vertical scalar width도 root의 `pool_vertical_sigma_m` bounded
  object 또는 `pool_vertical_sigma_uncertainty`로 명시할 수 있으며, 각
  lower/upper source-shape corner가 별도 transport case로 남는다.
  droplet diameter/mass-fraction uncertainty도 `droplet_population.corners`
  와 별도 evidence ID로 complete simplex population을 선언하면 실제
  phase_model_options와 transport case까지 전파된다. nominal population이
  corner에 포함되지 않거나 fraction 합이 1이 아니면 입력을 거부한다.
  CLI execution record는 이 설정과 input SHA를 보존한다. 실제 예제 diagnostic
  실행은 physical 1/1 corner 완료, 0 withheld, conservative pool source 1건
  (0.2176433778 kg/10 s)이었고 `--no-refinement` 때문에 aggregate decision은
  의도대로 `withheld`/exit 2였다.
  synthetic short-release 실행에서 nominal/대체 droplet population simplex
  2 corner도 실제 phase ledger 2/2 완료로 확인했다.
- 운영 decision은 local wind plane에 실제 삽입되지 않은 declared cuboid를
  이제 무조건 보류한다. off-plane 및 완전히 upwind인 구조물은 2-D local plane에
  없더라도 3-D wake가 무해하다는 근거가 없기 때문이다. 반대로 plane을 가로지르는
  cuboid는 full-plane mask라는 conditional 근사로만 계산되며, validation된 wake
  correction으로 해석하면 안 된다.
- 장애물 입력은 global axis-aligned box뿐 아니라 회전된 직사각 footprint도
  지원한다. `OrientedCuboid`/strict JSON form은 centre, length, width, z bounds와
  global `+x` 기준 반시계 `math_to` long-axis bearing을 명시하며, 기상
  wind-from 방위와 혼용하지 않는다. 매 wind corner에서 실제 footprint를
  투영하고 report는 global geometry와 local mask를 모두 보존한다.
  입력 형식은 `docs/field-obstacle-input.md`에 있다. 이는 기하 정확도 보강일
  뿐 3-D wake/recirculation closure나 validation을 추가한 것이 아니다.
  축·치수·방위에 명시적 bounded geometry와 evidence ID를 붙이면 각
  obstacle-mask corner가 source/weather corner와 함께 재투영·재계산되며,
  nominal obstacle screen은 bounds가 해소되기 전 `withheld`로 남는다.
- strict `field-screen` case는 `coordinate_reference`를 반드시 요구한다.
  여기에 도면/grid ID, origin ID, vertical datum, evidence ID와 true east에서
  반시계로 잰 site `+x` 축 방위를 남긴다. code는 meteorological wind를 먼저 이
  site grid로 회전한 뒤 센서·장애물 plane을 계산한다. 따라서 engineering grid를
  동/북 grid로 조용히 가정하지 않는다. 기존 Python API에서 coordinate reference를
  생략한 연구 계산은 호환을 위해 east/north site grid를 전제로 유지된다.
- `lh2_validation_available`는 더 이상 독립적인 신뢰 플래그가 아니다.
  strict case에서 `true`를 쓰려면 SHA-256, 공통 clock, source/weather/
  obstacle/receptor geometry, temporal operator를 포함한
  `validation_evidence`가 필요하다. free-field evidence는 obstacle transport
  검증으로 승격되지 않으며, 입력과 field report 모두 fingerprint를 보존한다.
- `field_validation.py`와 `field-validate`는 별도 관측 CSV를 실제 고정
  수용점 검증 입력으로 받는다. 센서별 한 행, SHA/row count, common clock,
  obstacle ID와 averaging operator를 확인한 뒤 MAE/RMSE/bias만 계산하며,
  관측하지 않은 hazard distance나 설계기준으로 외삽하지 않는다.
  validation evidence의 source/weather/receptor-geometry/temporal-operator
  ID와 비교 basis ID가 다르면 fingerprint가 맞아도 dataset 생성 단계에서
  fail-closed로 거부한다. 물리 obstacle geometry ID와 모델 mask/wake
  representation ID는 서로 다른 namespace로 유지해 scoring gate에서
  별도로 확인한다.
- `field_source_io.py`와 `degali field-source`는 이미 대기 중인 post-flash/
  pool/droplet vapour rate CSV만 strict하게 받아 SHA, source boundary, common
  clock을 묶는다. 마지막 endpoint rate는 0이어야 하며, 압력·level·파일명에서
  phase나 누출률을 추론하지 않는다. 양의 integrated release mass도 요구한다.
  `rate_lower_kg_s`/
  `rate_upper_kg_s` 쌍이 있으면 같은 clock의 deterministic source corners로
  보존하며 확률분포로 해석하지 않는다. Python `run_field_semi_fv_envelope()`의
  `atmospheric_source_schedule=` 경로는 이 corner들을 weather/surface/sensor
  corner와 결합하고, 대체된 upstream source 값은 nominal로 고정해 이중 변동을
  피하며 envelope report에 provenance를 남긴다. 같은 인자를
  `run_field_operational_uncertainty_envelope()`에도 전달하면 source corner마다
  numerical refinement와 fail-safe operational decision을 다시 적용한다.
  `degali field-screen --atmospheric-source-case`는 별도 source JSON의
  경로/SHA와 schedule provenance를 실행 기록에 보존하면서 이 operational
  envelope에만 연결하며, measured-history/phase-routing 입력과 조합할 수 없다.
  direct envelope은 post-flash vapour만 받으며, 위치·수직폭이 없는
  pool/droplet schedule은 nozzle에 배치하지 않고 distributed/phase-routing
  handoff로 보류한다.
- nominal `direct_vapour_schedule`를 Python request에 직접 붙인 경우에도
  operational gate는 upstream thermodynamic/rate 경계만 대체된 것으로
  해석한다. source location/direction bound는 여전히 local plane 입력이므로
  matching envelope 없이 `uncertainty_resolved`가 되지 않으며, nominal
  schedule만으로 conditional screening을 허용하지 않는다.
  이 회귀를 포함한 decision/workflow/operational-envelope/source-sensor/batch
  표적 묶음은 **70 passed**로 재확인했다.
- atmospheric source schedule이 source rate를 대체하더라도 선언된 source
  location/direction corner는 유지한다. 고정 schedule이 표현하지 못하는
  non-exact release-duration uncertainty는 time-aligned schedule corner 없이
  거부한다.
- generic `FieldDistributedVapourSource`에도 evidence-backed bounded
  `position_uncertainty_m`와 `vertical_sigma_uncertainty`를 선언할 수 있다.
  semi-FV envelope는 각 Cartesian source-shape corner를 별도 계산하고,
  operational gate는 이를 전파하지 않은 nominal distributed-source 결과를
  보류한다. envelope report는 원래 bound와 source를 보존한다.
- `field_evidence_audit.py`와 `degali field-audit`는 외부 CSV/JSON/XLSX를 읽기
  전용으로 훑어 source boundary·weather·obstacle geometry·고정 수용점·
  common clock 다섯 채널의 준비도만 분류한다. `candidate_complete`도
  파일 간 사건/clock/기하 일치를 뜻하지 않으며 `promotion_allowed=false`를
  유지하고 `FieldValidationEvidence`를 자동 생성하지 않는다. 2026-10-05
  연결된 `SLABx_LH2` 폴더 전체를 점검한 결과 15,481개 구조화 파일 중 필요한
  채널을 모두 갖춘 후보는 없었고, 동적 풀 반경/CFD 시나리오 인벤토리는
  고정 H₂ 수용점 검증으로 분류되지 않았다. 큰 입력 디렉터리에서 scan limit에
  걸리면 `scan_complete=false`와 withheld/partial 사유를 남기며 완전 후보로
  승격하지 않는다. 센서 위치·4% 통과시각만 있는 Tsinghua 요약표도
  `receptor_observations`로 세지 않고 비정량 note로만 보존한다. XLSX는
  worksheet 앞부분만 bounded sample하고 row count·sheet identity를 승격하지
  않는다. 지원센터 폴더 재감사에서는 CSV/XLSX 20개가 모두 required channel
  없이 `withheld`였다. 2026-10-06 read-only 재감사에서는 명시적 weather
  identity가 없는 wind/temperature 단서를 weather 채널로 승격하지 않도록
  한 현재 판정 기준에 따라 `SLABx_LH2`의 15,604개 구조화 파일도
  `withheld`로 분류되었다. source_boundary, weather, obstacle_geometry,
  receptor_observations, common_clock가 모두 승격되지 않았고
  `promotion_allowed=false`였다. 지원센터 20개와 관련 자료 0개도 각각
  계속 `withheld`였다. 이 재감사 결과를 근거로 validation evidence나
  site acceptance를 생성하지 않았다.
  최신 실행 산출물은 `outputs/current-slabx-audit.json`,
  `outputs/current-support-center-audit.json`,
  `outputs/current-related-audit.json`에 보존했다.
- `obstacle_wake_validation.py`의 AIJ Case-H/SMEDIS 관측 행은 이제 좌표의
  유한성, 누락 sentinel의 명시성, 농도·RMS·난류운동에너지의 비음수/유한성
  경계를 생성 시점에 검사한다. 따라서 내부적으로 손상된 행은 neutral/dense
  gas 관측 triage에도 들어오지 않으며, 기존의 `quantitative_*_allowed=False`
  및 LH₂ obstacle-wake 승격 보류 경계는 그대로 유지된다. 해당 표적 회귀는
  **6 passed**다.
- `field_evidence_audit.py`는 header-only CSV와 빈 JSON row array를 실제
  source/weather/obstacle/receptor/clock 채널로 승격하지 않는다. 해당 파일은
  명시적 no-data note가 있는 diagnostic candidate로만 남아
  `candidate_complete`를 만들 수 없으며, 이 readiness 보강 표적 회귀는
  **16 passed**다. bounded XLSX sample도 첫 non-empty row 뒤의 실제 행이
  없으면 같은 방식으로 보류하며, XLSX row count 불확실성은 유지한다.
  값이 모두 빈 CSV 행과 JSON record도 같은 no-data 경계로 보류한다.
  `N/A`, `unknown`, `unspecified`, `not-recorded` 같은 placeholder도
  identifier/measurement로 세지 않는다.
  한 source/receptor record 내부의 필드 co-occurrence도 확인해 서로 다른
  행의 ID·rate·좌표·농도를 합성하지 않는다.
  audit record는 이제 다섯 채널이 한 파일 안에서 함께 확인된
  `coherent_candidate_paths`를 별도로 보존한다. 빈 목록이면 채널이 파일
  간에 분산된 것이며, 목록이 있어도 promotion은 허용하지 않는다.
- 이 시점에 현장/반해상 선택 회귀 286개(최근 XLSX 감사기·
  atmospheric source adapter, operational source-corner, named batch
  source-envelope, source geometry-corner/operational source-envelope,
  distributed-source-shape, strict distributed-case parser 및 validation
  basis/evidence identity gate 테스트 포함)를 통과했다. 전체 suite의 마지막
  complete run은 `1453 passed, 145 skipped`였고, 그 뒤 envelope/refinement
  contract 테스트를 추가해 current field subset을 `290 passed`로 확인했다.
  path-provenance 회귀를
  추가하기 직전 historical field/site subset은 `286 passed, 8 skipped`였다.
  이번 추가 변경의
  validation/comparison/audit/workflow/batch 집중 회귀는 `67 passed`였고,
  datum projection 및 수직-domain gate, applicability contract 보강 후
  모든 `test_field_*.py`와 semi-FV obstacle 회귀는 `264 passed`였다.
  operational decision/status contract 보강 후 decision/batch/operational
  envelope/history/phase 회귀는 `34 passed`였다.
  applicability contract 형식 검증까지 포함한 contracts/decision/workflow/report
  회귀는 `70 passed`였다.
  duration-corner preflight를 추가한 현재 field workflow/operational-envelope/
  case-IO 집중 회귀는 `78 passed`였고, 해당 실행에서 `compileall`과
  `git diff --check`도 함께 통과했다.
  모델 비교 불변식 보강 후 `test_field_*.py` 및 semi-FV obstacle 회귀는
  `267 passed`였고, comparison/comparison-IO 집중 회귀는 `19 passed`였다.
  distributed source rate-corner 추가 후 같은 field/semi-FV 묶음은
  `269 passed`로 다시 통과했고, `compileall`과 `git diff --check`도 통과했다.
  이후 rate-bound zero-endpoint 보강의 targeted workflow/case-IO/operational
  회귀도 `3 passed`로 확인했다.
  allowed decision의 refinement/uncertainty 불변식 보강 후 decision/operational/
  history/phase 회귀는 `31 passed`였다.
  현재 모든 `test_field_*.py` 및 semi-FV obstacle 회귀 기준은 `307 passed`이며,
  같은 실행에서 `compileall`과 `git diff --check`도 통과했다.
  이후 공통 `BoundedValue`/`CircularBoundedValue`와 저수준 source schedule/
  obstacle/semi-FV 수치의 boolean 강제변환, source/weather/surface/sensor/
  scenario/request의 untyped nested boundary를 생성 시 거부하도록 보강했으며,
  최신 field/semi-FV 회귀는 **312 passed**다.
  이후 blowdown/measured-history direct-vapour 장부의 음수 phase flow 및
  direct-vapour+unrouted-liquid 불일치 fail-safe, 센서 gain/averaging/좌표와
  semi-FV source-height의 boolean/비유한 입력 거부를 추가했으며, 관련
  집중 회귀는 **66 passed**였다.
  flash result/source preparation의 nested type·phase flow·residual·closure
  diagnostic 생성 검증과 blocked 결과의 flash-plane 비노출을 추가했으며,
  LH2/workflow/phase-routing 집중 회귀는 **67 passed**였다.
  이어서 `ReleaseSource`, `WeatherState`, `FieldSemiFVRequest` 및
  distributed-source 위치/폭의 boolean numeric coercion을 생성 시 거부하도록
  보강했다. 관련 계약·워크플로 회귀는 `63 passed`였고, 전체
  `test_field_*.py` 및 semi-FV obstacle 회귀는 **319 passed**로 통과했다.
  이번에는 완료된 field trace를 `FieldModelSensorSet`으로 내보내는
  `field_model_sensor_set_from_screening()` 어댑터를 추가했다. peak/final/
  full-trace average와 true/indicated 채널을 반드시 명시하게 하고,
  withheld sensor·blocked 실행은 거부하며, source applicability와 trace
  operator provenance를 비교 결과에 보존한다. 비교/validation/workflow
  집중 회귀는 `81 passed`, 전체 field/semi-FV 회귀는 **322 passed**였고,
  `compileall` 및 `git diff --check`도 통과했다.
  이어서 direct/pool/droplet sensor superposition도 같은 비교 계약으로
  내보내도록 `field_model_sensor_set_from_superposition()`을 추가했고,
  superposition이 steady/transient branch를 섞지 않도록 보강했다. 관련
  comparison/superposition 회귀는 `23 passed`였다.
  이어서 superposition 각 branch의 source/sensor/ambient/obstacle 및
  distributed-source bound와 stability/phase/history resolution flag를
  검사해 unresolved corner를 `uncertainty_complete=False`로 보존하도록
  보강했다. 관련 집중 회귀는 `24 passed`였다.
  또한 `FieldSemiFVEnvelope`가 duplicate corner selection을 생성 시 거부하도록
  보강해 상위 operational envelope의 set-based completeness 검사가 중복과
  누락을 상쇄하지 못하게 했다. 관련 workflow/operational 회귀는 `53 passed`다.
  measured-history flash/source-history/joint-history envelope도 typed case,
  bound-selection key, warning 및 unique corner을 생성 시 검증하도록 맞췄고,
  historian 집중 회귀는 `26 passed`였다.
  phase-routing result/transport corner, droplet-handoff refinement, historian
  import 및 batch-export 결과 컨테이너에도 동일한 construction-time
  type/finite/duplicate/complete-case gate를 추가했다.
  fingerprinted atmospheric source schedule의 직접 lower/upper bound 생성도
  CSV 경로와 동일한 zero-endpoint 및 bound-column 검증을 적용한다.
  field screening result/base envelope의 sensor·projection·corner collection도
  tuple을 요구해 mutable list가 deterministic report coverage를 우회하지
  못하게 했다.
  validation/comparison 집중 회귀는 `38 passed`로 통과했다.
  sensor observation/superposition targeted 회귀는 `5 passed`였고, 해당 실행에서
  `compileall`과 `git diff --check`도 통과했다.
  operational source-geometry,
  distributed-source-shape, bounded obstacle-geometry/phase propagation 및
  strict distributed-case parser 회귀를 포함한
  상태에서 `compileall`과
  `git diff --check`도 통과했다.
  추가로 semi-FV worst-step 질량 잔차가 절대/상대 수치 gate를 넘으면
  transport 진단은 보존하되 field screening을 `blocked`로 만드는 보존성
  gate와 회귀 테스트를 반영했다. 같은 진단에 `primary` 및 각
  `distributed:<label>`별 실제 주입 질량 원장도 보존해 다중 source
  bookkeeping을 총량과 분리해 감사할 수 있다. 원장 합계와 총 주입량의
  별도 residual gate도 추가되어, 격자 inventory가 우연히 닫혀 보여도
  source bookkeeping이 어긋난 경우에는 field 결과를 보류한다.
  read-only field-audit도 명시적 weather identity 없는 wind/temperature
  자료를 weather channel로 승격하지 않도록 보강했다.
  모델 비교 CSV provenance에는 averaging-time 열의 실제 존재 여부도
  남기며, basis의 averaging time을 CSV에서 검증할 수 없으면 수치 행은
  유지하되 model-selection impact를 withheld로 둔다.
  분산 atmospheric source도 지원된 `source_kind` 집합 밖의 값은 transport
  전에 거부하도록 보강했다.
  모델 비교 입력도 CSV provenance 타입과 row count를 실제 prediction 집합과
  대조해 metadata-only mismatch를 차단한다.
  matched validation CSV도 evidence digest/row count뿐 아니라 evidence path가
  실제로 읽은 CSV와 같은 파일을 가리키는지 확인하여, 동일 바이트 복사본에
  다른 dataset provenance가 붙는 경로를 차단한다.
  보존형 semi-FV obstacle도 상단이 선언된 수직 domain 밖으로 나가면
  solid mask가 조용히 잘라지지 않고 입력 단계에서 거부하도록 보강했다.
  전역 obstacle이 local z=0 datum을 가로지르는 경우도 동일하게 명시적
  `blocked` projection으로 반환하고, datum 아래에 완전히 놓인 geometry는
  삽입하지 않았다는 진단을 보존하도록 보강했다.
  operational screening decision record도 status/boolean gate flag와
  withheld 사유의 형식을 생성 시 검증해, 잘못된 상태 객체가 안전한
  screening 결과처럼 보고서에 직렬화되지 않도록 보강했다.
  동일한 fail-safe 형식 검사를 `FieldApplicability`의 uncertainty flag와
  reasons/warnings에도 적용해, malformed applicability가 report/gate로
  전파되지 않게 했다. 추가로 `blocked + uncertainty_complete=True`와
  `non-blocked + reasons` 조합은 생성 시 거부하지만, obstacle/validation
  scope 경고를 가진 accepted physical result는 기존 conditional-warning
  설계를 보존한다.
  field envelope도 distributed source schedule이 primary release-duration
  lower/upper corner를 넘어가면 corner 생성 중 일반 예외로 중단하지 않고,
  time-aligned source envelope 또는 명시적 post-release continuation을
  요구하는 명시적 입력 오류로 거부하도록 보강했다.
  모델 비교의 basis/prediction row/comparison result도 숫자·boolean 타입,
  ratio/difference/geometry 일관성, matched sensor coverage 및 disagreement/
  mean-difference 불변식을 생성 시 검증해 malformed direct API 객체가
  비교 report로 직렬화되지 않게 보강했다.
  generic distributed atmospheric source도 동일 source ID/time axis의
  lower/nominal/upper rate schedule을 받을 수 있고, 위치·수직폭 corner와
  함께 실제 transport case로 전파된다. nominal operational screening은
  이 rate bound가 미전파된 동안 보류하며, report는 원래 schedule bounds를
  보존한다.
  operational screening decision도 required refinement/uncertainty flag가
  미충족인데 `screening_allowed` 또는 `conditional_allowed`로 직접 조작된
  객체를 생성 단계에서 거부한다.
  matched validation score도 dataset/model/comparison identity와 finite
  error metric, row-derived MAE/RMSE/bias/max 및 qualified-status 일관성을
  생성 시 검증해 fabricated validation report를 차단한다.
  추가로 validation dataset의 관측 tuple·evidence row count/path와
  comparison right-side observed model을 직접 대조하고, conditional score는
  반드시 reason을 갖게 했다. model-comparison 직접 객체도 mutable list 입력,
  untyped basis/prediction, blocked 상태의 fabricated rows/disagreement를
  생성 단계에서 거부한다.
  범용 distributed vapour source와 `FieldSemiFVRequest`도 position/warning,
  obstacle/sensor/source-corner collection의 list 입력을 거부해 생성 후
  deterministic source/geometry coverage가 변하지 않도록 했다.
  operational screening decision도 reason/action을 immutable tuple로 고정하고,
  직접 만든 `conditional_allowed` 상태에는 설명 사유와 review action을
  요구한다.
  operational uncertainty/history/sensor-array/source-sensor/phase aggregate도
  모든 corner decision의 status, refinement availability, uncertainty
  resolution을 대조해 fabricated aggregate allowed status를 차단한다.
  개별 allowed corner도 완료된 screening/refinement 결과를 직접 가져야 한다.
  phase aggregate의 sensor/case collection도 tuple과 non-empty corner를
  요구하고, 빈 sensor-array aggregate는 withheld/resolution 미완료만 허용한다.
  `FieldBatchCase`도 physical applicability/status와 refinement
  completion/availability를 대조해 manifest 상태 조작을 차단한다.
  `FieldBatchExport`는 manifest/report 파일의 실제 존재·output directory
  포함·case label 파일명 일치를 확인한다.
  sensor trace와 multi-source superposition도 finite/strictly ordered time,
  equal array length, non-negative concentration 및 true mole-fraction
  consistency를 생성 시 검증해 malformed combined field report를 막는다.
  field screening/sensor-result record도 typed subrecord, 중복 sensor label,
  obstacle projection, transport 없는 trace 및 withheld 사유를 생성 시
  검증한다. 다만 불완전한 sensor_results를 operational gate가 별도로
  `withheld` 처리하는 방어 경로는 유지한다.
  Python `FieldSemiFVRequest.direct_vapour_schedule`도 strict CSV/source
  boundary와 동일하게 마지막 rate 0 endpoint를 요구해, 선언된 finite
  duration 뒤의 source continuation을 암묵적으로 허용하지 않는다.
  plain field uncertainty envelope, sensor-calibration envelope 및 refinement
  result도 empty corner set, malformed/duplicate/non-finite selection,
  untyped nested record를 report 직렬화 전에 거부한다. 이 추가 계약 테스트
  반영 후 field/semi-FV 회귀는 `290 passed`였다.
  public model-comparison function도 non-set 입력과 boolean threshold/
  position-tolerance를 숫자로 암묵 변환하지 않고 즉시 거부한다.
  strict/direct `FieldModelComparisonCase` 계약도 같은 boolean threshold/
  position-tolerance 거부를 적용해 parser와 Python API의 수치 경계를
  일치시킨다.
  matched field-validation scoring과 `FieldValidationCase`도 같은 boolean
  threshold/tolerance 거부를 적용한다. comparison/comparison-IO/validation
  표적 회귀는 현재 **38 passed**이며, direct case의 model/dataset nested
  type도 생성 시 검증한다.
  `FieldModelDecisionImpact`도 status, alert sensor ID, monitored-distance,
  runtime ratio 및 matched comparison과의 일관성을 생성 시 검증해 malformed
  model-selection 결과가 report로 전파되지 않게 한다.
  operational uncertainty aggregate도 각 corner의 `uncertainty_resolved`를
  논리곱으로 보존해, unresolved corner가 있는 withheld aggregate를
  uncertainty-complete로 표시하지 않는다.
  phase-routing/pool operational aggregate도 같은 규칙을 적용해 refinement
  부재와 physical uncertainty 미해소를 report에서 구분한다.
  measured-history operational aggregate도 withheld/conditional/allowed
  경로 모두에서 corner별 `uncertainty_resolved` 논리곱을 보존한다.
  phase-routing operational corner record도 phase/sensor selection의 finite·
  unique 조건과 screening/refinement/decision nested type을 검증하며,
  field result 없는 allowed decision을 거부한다.
  general 및 measured-history operational corner record도 selection의
  finite·unique·JSON-safe 조건과 envelope case type을 검증하고, required
  refinement 없는 allowed corner를 거부한다.
  phase-routing operational corner도 screening object가 존재하더라도
  required refinement가 없으면 allowed 상태를 거부한다.
  phase-routing operational envelope의 sensor selection 선언도 finite
  numeric value, unique key/corner 계약을 통과해야 physical/calibration
  corner 비교로 진행된다.
  sensor-calibration uncertainty envelope도 malformed pair, 중복 corner와
  completed screening의 empty case set을 거부한다.
  `run_field_operational_sensor_array_envelope()`는 shared scalar transport를
  재사용하면서 detector calibration corner별 exact trace, refinement와
  operational decision을 집계한다. 이는 sensor calibration만 해소하며
  source/weather/geometry 불확실성은 별도 gate로 계속 보류한다.
  사용 예와 범위는 `docs/field-operational-sensor-array.md`에 기록했다.
  strict `degali field-screen --sensor-array-envelope`도 이 경로를 사용하며,
  source/history/phase envelope와의 혼용을 거부하고 execution JSON에
  `sensor_calibration` kind를 남긴다.
  source/weather/geometry와 supplemental detector bounds를 함께 다루는
  `run_field_operational_source_sensor_envelope()` 및
  `--joint-source-sensor-envelope`도 추가했다. physical corner transport는
  nominal calibration으로 한 번만 풀고, 각 calibration corner trace에
  refinement/decision을 붙인다. 실행/범위 설명은
  `docs/field-operational-source-sensor.md`에 있다.
  공개 FFI/Spadeadam 수평 방출에 대해서는
  `degali.validation.ffi_source_state`의 source-rate/orifice/height/pressure,
  wind speed/direction, ambient 및 wind-reference bound를 Cartesian corner로
  전파하고 각 보고 센서 좌표의 lower-bound-aware residual map을 만드는
  `run_ffi_source_state_envelope()`도 추가했다. 이는 기존 자유장
  `hydrogen_jet`를 변경하지 않는 deterministic sensitivity이며, 센서 간
  plume miss·upwind·trajectory 미도달은 0으로 대체하지 않고 withheld로
  남긴다. 사용 범위는 `docs/ffi-source-state-envelope.md`에 있다.
  새 FFI source-state 경로의 API/CLI/JSON-safe/실패안전 테스트는 현재
  `6 passed`이며, 새 관측 연산자 집계 테스트는 `3 passed`다. 전체 suite
  최신 complete run은 `1474 passed, 145 skipped`였다. source-state report에는
  이제 각 corner의 `arc_max_table`/`sensor_height_table`도 포함된다.
  Test 6에 rate 0.82/0.833/0.85 kg/s와 wind 2.2/2.3/2.4 m/s를 넣은
  4-corner CLI 재감사는 arc 4행/height 12행을 생성했고, 첫 corner는
  lower-bound constraint `violated`로 남았지만 operational/validation
  qualification은 모두 `false`였다.
  실행 JSON은 이제 `reference/spadeadam/conditions.csv`와 `sensors.csv`의
  SHA-256을 보존하고, 계산 `complete`와 별도로 operational/validation
  qualification을 항상 `false`로 명시한다.
  lower-bound 관측 특성에 맞춰 센서별 one-sided 만족률과 worst deficit도
  보존하며, 이는 대칭 정확도 점수나 site acceptance가 아니다.
  전체 파일 상태는 여전히 의도적으로
  dirty/untracked이므로,
  공개 commit 전 `reference/`, 원본 FORTRAN, 외부 raw data가 제외됐는지
  반드시 별도 점검한다.

- 현장 batch 산출물의 교차 계층 불변식도 회귀로 고정했다. off-plane 장애물은
  case report의 `obstacle_projections[].local_obstacle=null`로 남고, 같은
  case의 off-plane detector는 `sensor_results[].withheld=true`로 남으며,
  manifest operational decision은 두 식별자를 이유에 포함한 `withheld`를
  유지한다. `tests/test_field_batch.py`와 `tests/test_field_report.py` 표적
  회귀는 현재 `19 passed`이며, 이는 obstacle/wake 정량 검증을 추가한 것이
  아니라 report와 fail-safe gate가 서로 다른 상태를 내보내지 않는지 확인하는
  출력 계약 검사다.

다음 우선순위는 site-specific historian 기준의 calibration evidence와
관련 LH2 obstacle 데이터로 source/phase 및 obstacle closure를 독립 검증하는
것이다. field flash handoff의 mass/momentum/energy residual gate는 구현됐지만,
사이트 대표성은 외부 증거 없이는 승격하지 않는다. 데이터 없이 wake 계수나
3-D 난류 보정값을 추가하지 않는다.

연결 폴더의 SLABx-LH2·저장탱크 운전 자료 적용성은
`docs/field-evidence-readiness-2026-10-05.md`에 정리했다. 현재 1초 벤트
이력은 sampling/운전 plausibility에는 유효하지만, 독립 mass flow·phase·풍장·
하류 H2 receptor가 함께 확보되지 않아 누출 확산 검교정으로 승격하지 않는다.

## 0A. 2026-10-02 SLABx/SLABx-LH2 물리 공백 보완

- `run_lh2_finite_release_research()`가 보존 근접장 → yawed 횡풍 플룸 →
  H2 종질량 clock의 정확한 plume/puff 전환 → native 3-D Gaussian puff를
  한 번에 실행한다. puff는 총질량/H2/3축 운동량/상대 총에너지,
  entrainment, 부력, form drag, 지면 접촉, 수평 중력확산 및 고정 수용점
  시간 이력을 계산한다.
- 전환 이후 puff는 `WindHistory`의 기상학적 풍향/풍속을 동·북 벡터로
  보간하며 관측 시각 경계를 정확히 분할한다. 기록이 방출 시작부터 전환까지
  존재하면 전환 전 풍속/풍향 변동 적용성도 자동 판정한다. 한계 초과 시
  조건부 표시하며 strict 모드에서는 즉시 거부한다. 전환 전 plume 자체는
  아직 하나의 선언된 정상풍 벡터를 사용한다.
- `run_lh2_rainout_pool_research()`가 post-flash 입력 → 크기군별 액적의
  횡풍/중력/d² 증발 → rainout → 방출 중 축대칭 얕은층 풀 확산/국부 접촉시간
  기판전도 증발을 한 번에 실행한다. 직접 기체, 비행 증발, 공중 액체,
  지면 액체, 풀 증발, 잔류 질량 및 계산영역 유출을 하나의 H2 ledger로
  닫는다. 고정면적 경로는 `pool_model='fixed'`로만 선택한다.
- `YawedCrosswind.solve()`는 유한 방출 material clock에서 조기 종료하며,
  고수준 함수는 보존상태 역산 실패 시 공간 간격을 자동으로 절반씩 줄여
  재시도한다.
- 5 s, 약 0.253 kg/s 연구 예제에서 plume/puff 전환은 약 21.74 m였고,
  puff H2 잔차는 0, 총질량 상대잔차는 약 1.6e-15였다. 이는 실행성/보존
  확인이며 독립 현장농도 검증을 뜻하지 않는다.
- 동적 풀 1차 통합·전환 전 풍장 gate까지 반영한 전체 회귀는 2026-10-02에
  **1158 passed, 145 skipped, 0 failed**로 완료했다(1270.93 s). 이후 검증
  기본 증발운동량 closure 정렬과 물 일정열유속 경계 추가 후 관련 회귀도
  **36 passed, 0 failed**로 다시 확인했다. 대표 5 s
  플래시/rainout 사례는 풀 유입 0.8254028 kg, 증발 0.8254028 kg,
  계산영역 유출 0 kg, 풀 원장 최대잔차 7.11e-15 kg, 전체 H2 잔차
  -1.04e-14 kg로 닫혔다.
- SLABx-LH2가 사용한 DynamicLH2PoolX와 동일한 물성·유입·격자 및 검증 기본
  `zero_radial_momentum_vapor` 조건으로 별도 대조했으며, 0--5 s 전 시각의
  액체질량과 누적증발량 최대 차이는 각각 **0 kg**, 최종 전면반경 차이도
  **0 m**였다. 이는 구성요소 수치동등성 확인이며 현장 반경 검증은 아니다.
- 로컬 공개 JUEL-3155 판독자료로 DEGALI 계산기를 직접 재실행해 동결된
  보류시험 지표를 재현했다: 물 Trial 4 반경 RMSE **0.13744 m**,
  알루미늄 Trial 6 **0.10159 m**, 60 s 반경 0.40/0.46 m, 최대 질량원장
  잔차 4.81e-13/3.49e-13 kg. 따라서 적용범위는 “선언 지면유입 + 평탄한
  매끄러운 축대칭 풀”의 제한 구성요소 검증 상태다.
- 상세 범위는
  [docs/slabx-physics-gap-closure-2026-10-02.md](docs/slabx-physics-gap-closure-2026-10-02.md)를
  우선한다. native puff 농도와 동적 pool 확산은 아직 독립 검교정 전 연구
  결과이며 운영/인허가 기본값으로 승격하지 않는다.

## 0. 2026-09-21 Test 6 후속 보완

- 공개 DNV 보고서에는 Test 6의 10 Hz 풍속·풍향·H2 **그림**이 있다. “시간열이
  없다”는 기존 표현은 “기계 판독 가능한 원시 채널은 공개되지 않았다”로
  정정했다.
- `src/degali/addons/transient_receptor.py`에 고정 센서 재생 연산자를 추가했다.
  풍향에 따라 센서를 순간 downwind/crosswind 좌표로 회전하고, 미리 계산한
  풍속별 정상해 테이블을 선형 보간하며, 선언된 `t90` 검출기 지연을 적용한다.
- `tools/audit_ffi_test6_transient_receptors.py`는 사용자가 제공한
  `time_s,wind_speed_ms,wind_direction_from_deg` CSV를 읽고 입력 SHA-256과 함께
  Test 6 센서별 평균/최댓값을 기록한다. 보고서의 min/mean/max로 가짜 시간열을
  만들지 않는다. 현재 도구는 `degali.ffi-test6-transient-receptor-execution.v1`
  compact artifact를 생성하며, over-range 행은 `null`로 보존하고
  `promotion_allowed=false`, `validation_qualified=false`를 고정한다. 실제
  owner-exported trace가 없으면 저장소는 여전히 합성 시간열을 만들지 않는다.
  `--verify` 경로는 plume을 재실행하지 않고 이력/reference SHA와 센서
  identity·window invariant만 integrity-check한다(`report_recomputed=false`).
  observation-operator discriminator, no-storage/source-history flag, response
  time, window sample count, sensor coverage 및 summary peak도 다시 계산해
  compact row 변조를 fail-closed한다.
  추가로 history time range, interval range, uniform-sampling diagnostic와
  wind-speed range를 저장·재검증하지만, irregular/sparse history를 자동
  quality 승인하지 않는다.
- 독립에너지 횡풍 모델에는 기본 OFF인 고정 지면 열경계를 추가했다. 계수 0은
  단열 경계이고, 양의 `h`와 표면온도를 모두 명시하고
  `include_ground_heat_transfer=True`일 때만 접촉폭에 비례한 열이 들어간다.
- 기존 finite wall-jet 단독 보정은 30 m 잔차를 해결하지 못해 계속
  비권장이다. 새 보완도 Test 6 피팅이나 기본 모델 승격이 아니다.
- 전체 회귀: **1132 passed, 145 skipped**, 실패 없음(2026-09-21).

## 1. 작업 디렉터리와 보존 원칙

현재 작업 디렉터리:

```text
<repository-root>
```

중요 원칙:

- 현재 worktree는 의도적으로 매우 dirty하다. 기존 변경사항을 삭제하거나
  `git reset --hard`, `git checkout --`를 실행하지 않는다.
- 원본 DEGADIS Fortran, 외부 실험 원자료, PDF/Excel/압축파일은 GitHub·PyPI·DOI
  배포물에 넣지 않는다. `reference/`는 로컬 검증용이며 `.gitignore` 대상이다.
- 외부 자료는 출처와 축약된 검증 결과만 문서화한다. 원자료를 `data/`나
  `tests/`에 복사하지 않는다.
- 새 물리 계수는 FFI/E3.4/E3.5 관측값에 맞춰 피팅하지 않는다. 물리 범위, 공개
  문헌, 독립 검증을 먼저 선언하고 opt-in 연구 경로로 둔다.
- 아직 commit/push/publish하지 않은 변경이 많다. 인수받은 계정에서 먼저 diff와
  패키지 제외 목록을 검토한 뒤 별도 커밋한다.

현재 `git status`에는 약 100개 이상의 수정/미추적 파일이 있다. 이는 이전
개발 작업 전체가 아직 하나의 공개 커밋으로 정리되지 않았다는 뜻이며, 작업을
잃어버렸다는 뜻이 아니다.

최근 `slow` marker를 제외한 저장소 회귀는 **1,556 passed, 134 skipped,
63 deselected, 0 failures**로 완료됐다. 63개는 명시적으로 표시된 장시간
통합시험이므로 전체 저장소 완전 회귀 기준을 대체하지 않는다.

## 2. 현재 모델 상태

### 기본 경로

- 기존 DEGADIS 2.1 재구현을 기반으로 한 LH2 source/jet/pool 경로가 있다.
- CoolProp 기반 실제 수소 물성, flashing source, jet plume, pool evaporation,
  quasi-steady plume 연결을 사용한다.
- `degali.lh2.assess()`는 결과와 함께 `warnings`, `screening_scope`를 반환한다.
- `strict_scope=True`이면 적용범위 경고를 `ApplicabilityError`로 차단한다.
- `assess_envelope()`는 사용자가 선언한 source rate·풍속 격자를 모두 계산한다.

### 새로 추가된 연구/감사 경로

`src/degali/lh2.py`:

- `run_lh2_yawed_crosswind_research()`
  - 임의 global wind/release bearing의 수평 횡풍 연구 경로
  - reverse-axial branch 거부
  - `validated=False`인 연구 경로이며 기본 `assess()`에 섞지 않는다.
- `assess_pool_history(..., response_time_s=...)`
  - 시간변화 pool source를 quasi-steady plume snapshot으로 연결
  - 명시된 1차 응답시간만 적용하며 transient puff/meander solver가 아니다.
- `assess_observation_envelope()`
  - source rate·풍속 가설을 전부 계산하는 no-fit admissible-set 진단
  - 여러 조합이 관측 배수 안에 들면 `non-identifiable`로 남긴다.
  - 기본 admissible-set 요약은 centreline 기준이며, FFI 감사 경로는
    `project_lh2_jet_to_sensors()`와 source-state corner별 arc/height 표를
    별도로 사용한다.
- `project_lh2_jet_to_sensors()`
  - `(x, y, z)` 실제 센서 좌표에서 Gaussian 수직·횡방향 프로파일을 직접 계산
  - 궤적 밖 점은 외삽하지 않고 `nan`으로 반환한다.

`src/degali/validation/model_comparison.py`:

- `ComparisonCase`, `ModelPrediction`, `compare_models()`
- source, 풍속·풍향, 방출 높이, 센서 연산자, 평균시간, phase closure, geometry가
  동일할 때만 모델 우열 ranking을 허용한다.
- 조건이 다르거나 관측값에 피팅된 모델은 통계는 계산하지만 ranking에서 제외한다.
- HyRAM/PHAST/EFFECTS를 패키지에 포함하지 않으며, 외부 실행 결과를 사용자가
  별도로 제공할 때만 비교한다.

`src/degali/validation/screening_gate.py`:

- `evaluate_screening()`과 `ScreeningDecision`
- 범위 밖 결과와 장애물 접촉 결과를 screening에서 차단한다.
- 조건부 결과는 `allow_conditional=True`를 명시해야 한다.
- 항상 `design_basis_allowed=False`, `approval_allowed=False`이다.

## 3. 현재 검증 기준점

문서화된 기준값:

- PRESLHY E3.5: 62 downstream arc maxima, MG **1.047**, VG **1.425**,
  FAC2 **0.839**.
- FFI/DNV six-arc screen: MG **1.245**, VG **1.373**, FAC2 **0.833**.
- FFI Test 6 30 m 보고 arc 최대: **21.0 vol%**.
- 보고 source 0.833 kg/s, low/mean/high wind 2.3/2.5/2.7 m/s.
- 측정 source·풍속 격자의 centreline 예측은 약 **6.08–7.87 vol%**이며
  21 vol%와 factor-2 이내가 아니다.
- 로컬 공개 extract의 30 m 센서 15개 좌표에 직접 투영한 최대값은
  **7.50 vol%**, 관측/예측 **2.80배**이다.
- source rate를 1.5–3.0 kg/s까지 확장하면 factor-2 조합이 나오지만,
  이는 현재 검증된 jet source 범위 0.084–0.285 kg/s 밖이고 여러 조합이
  동시에 맞아 source 식별이 불가능하다. 기본 source multiplier로 채택하지 않는다.

전체 회귀 기준:

- 최신 전체 테스트 집합(ambient/phase-routing 경계, pressure-bound table,
  validation-evidence conditional gate 및 FFI corner 관측 연산자 표 포함,
  2026-10-07)은 장시간 parity 프로세스를 분리해 검증했다:
  **1438 passed, 6 skipped** (parity 제외; non-field 3개 chunk 10분 41초,
  field/obstacle 표적 2분 32초) +
  **109 passed, 139 skipped** (reference parity, 7분 02초) =
  **1547 passed, 145 skipped, 0 failures**.
- FFI 관측 연산자·source-state 표적 테스트: **9 passed**; Test 6 30 m
  aggregate와 높이별 0.1/1.0/1.8 m 결과, source-state envelope를 함께
  검증했다.
- 최신 직접 validation/comparison 타깃 테스트:
  `tests/test_field_validation.py tests/test_field_comparison.py` **38 passed**.
- 최신 field/semi-FV obstacle 타깃 회귀: **348 passed** (2026-10-06,
  phase-routing result/transport corner, droplet-handoff refinement, historian
  import/batch-export construction gates, direct atmospheric source schedule
  zero-endpoint gate와 direct-schedule source-geometry gate,
  low-level/nested-boundary type gate,
  positive-downward settling convention, phase mass-ledger, sensor-boundary,
  flash-preparation, boolean numeric boundary, screening/superposition-to-model-
  comparison provenance 및 temporal-mode fail-safe 포함).
- 현재 field/evidence-audit/semi-FV/obstacle-wake 타깃 묶음은
  **384 passed**로 재확인했다. 여기에는 header-only CSV/JSON 및 XLSX
  readiness 보류, CSV duplicate/empty header·extra-value·bounded streaming,
  JSON duplicate-key fail-safe와 AIJ/SMEDIS 관측 행 물리 경계가 포함된다.
- typed field model sensor-set CSV writer의 fingerprint/round-trip 회귀도
  위 field 묶음에 포함된다. basis ID는 CSV에서 추론하지 않고 비교 case에
  명시하며, 기존 파일 overwrite는 fail-safe로 거부한다.
- strict comparison case의 CSV SHA-256/행 수 pinning과 post-authoring
  mutation 거부 회귀도 위 field 묶음에 포함된다.
- comparison report의 `csv_provenance.integrity_pinned`가 case의 기대 해시
  검증 여부를 명시하도록 보강했다. 단순 관찰 fingerprint와 고정 검증을
  후속 자동화가 구분할 수 있다.
- fixed-sensor comparison CSV의 중복/빈 header와 header보다 많은 row 값을
  parser가 거부하도록 해 `DictReader`의 silent column overwrite/loss를 막았다.
- sensor ID/좌표/prediction/averaging-time logical column role 충돌도
  입력 시 거부해 한 열의 이중 물리 해석을 차단했다.
- typed `FieldModelComparisonCase`를 strict JSON으로 내보내는 artifact bridge도
  추가했다. custom CSV column mapping, relative path, basis, fresh hash/row
  count를 기록하고 emitted JSON을 다시 읽어 round-trip을 검증한다.
- 범용 atmospheric source schedule CSV도 comparison CSV와 같은 duplicate/
  empty header, extra-value, column-role collision 경계를 적용하며, bound
  schedule column이 schedule object 없이 선언되는 malformed 상태를 거부한다.
- typed `FieldAtmosphericSourceSchedule`를 strict relative-path JSON case로
  내보내는 writer도 추가했다. source CSV hash/evidence를 재검증하고 emitted
  JSON을 다시 읽어 round-trip을 확인하며 기존 case overwrite는 거부한다.
- matched `FieldValidationCase`도 model/observed CSV column mapping, evidence
  identity, relative path와 hash를 보존하는 strict JSON writer를 갖게 됐다.
  양쪽 typed record와 CSV를 재검증한 뒤 emitted case를 다시 읽는다.
- matched validation CSV도 duplicate/empty header, extra-value와 sensor/
  geometry/time/concentration/clock/obstacle role collision을 score 전에
  거부하도록 강화했다.
- field-evidence CSV inventory도 duplicate/empty normalized header와 header
  열 수가 다른 행을 unreadable diagnostic candidate로 보류한다. `DictReader`
  의 silent column overwrite/extra-value loss가 readiness channel을 만들지
  못하도록 하는 감사 경계다. 대형 CSV는 모든 행의 열 수와 전체 row count를
  streaming으로 확인하되 첫 200개 data row만 값/co-occurrence 판정에 사용해
  감사 메모리 사용을 bounded하게 유지한다. JSON도 duplicate 또는
  normalization-colliding object key를 parser 단계에서 거부해 last-value-wins
  overwrite가 readiness channel을 바꾸지 못하게 한다. 각 candidate에는
  스캔한 파일 바이트의 lowercase SHA-256도 기록해 오래된 audit 출력과
  변경된 입력을 구분하지만, 이를 FieldValidationEvidence로 승격하지 않는다.
  `read_field_evidence_audit_json()`/`write_field_evidence_audit_json()` strict
  bridge도 추가되어 CLI execution schema를 round-trip하고 모든 candidate
  digest를 다시 확인하며 overwrite를 거부한다. candidate path가 audit root
  밖으로 escape하는 artifact도 파일을 열기 전에 거부한다.
  execution record에는 실제 `max_files` scan limit도 보존·교차검증해,
  다른 제한으로 만든 partial scan을 complete scan처럼 재사용하지 못한다.
  `degali field-audit-verify`도 추가해 저장 artifact의 candidate 파일을
  다시 검증하고 `--require-complete`에서 partial/withheld를 exit 2로 차단한다.
- measured-history historian CSV도 duplicate/empty header와 header 열 수가
  다른 행을 source flash import 전에 거부해 `DictReader`의 silent overwrite/
  extra-value loss가 P/T/flow/phase quality gate를 우회하지 못하게 한다.
- 범용 atmospheric source schedule, fixed-sensor model, matched-validation CSV도
  동일한 exact row-width 경계를 사용하도록 보강했다. truncated/extra row는
  행 번호와 함께 거부하고, stripped header 이름의 충돌과 header-only source
  schedule도 거부한다. 따라서 `DictReader`의 implicit missing/extra field가
  time/rate·prediction·observation mapping을 바꾸지 못한다.
- strict JSON reader 공통 guard(`field_json.strict_json_loads`)를 추가해
  screening case, source schedule, model comparison, matched validation,
  batch input 및 batch artifact의 exact/whitespace/case-normalized duplicate
  key를 schema mapping 전에 거부한다. 관련 case/source/batch/validation 및
  field regression은 **96 passed**로 확인했다.
- field-evidence XLSX bounded sampler도 CSV와 같은 empty/normalized-duplicate
  header fail-safe를 적용한다. 모호한 보조 sheet는 sheet-level diagnostic으로
  남기고 유효한 sheet는 계속 inventory하며, 모든 sheet가 모호한 workbook만
  unreadable diagnostic candidate로 남긴다. 건너뛴 sheet에서 readiness channel을
  만들지 않으며, `test_field_evidence_audit.py`는 **30 passed**다.
- `field_execution.py`와 `degali field-verify`를 추가했다. source schedule,
  model comparison, matched validation, field screening, batch execution JSON의
  input SHA와 참조 case/CSV/manifest/report를 다시 strict하게 읽는다.
  nominal/direct-source screening은 기록된 numerical option으로 report를
  재계산하고, 당시 옵션이 완전하지 않은 phase/history/sensor-array 및 batch
  record는 `report_recomputed=false`로 명시했다. 현재는 complete
  phase/history/sensor-array/source-sensor 옵션을 저장한 새 실행은 fresh
  replay하며, 여전히 불완전한 record만 `report_recomputed=false`로 남는다.
  입력/참조 artifact가 바뀌면 실패한다.
  verification record의 `promotion_allowed`는 항상 `false`다. 전용 및
  audit/batch 회귀는 **40 passed**로 확인했다.
- direct atmospheric-source `field-screen --uncertainty-envelope` fixture도
  실행 artifact를 실제 envelope solver로 재계산해 `report_recomputed=true`를
  확인했다. nominal/refined replay와 함께 execution 회귀 범위에 포함된다.
- `degali field-verify --require-recomputed`를 추가했다. 복합/구형
  screening 또는 batch처럼 report 재계산 없이 artifact 무결성만 확인된
  경우 exit 2를 반환해 자동화가 fresh replay 없는 결과를 허용하지 않게
  한다.
- batch verifier는 이제 각 report의 SHA만 확인하지 않고 JSON을 다시 읽어
  typed screening/operational-envelope schema와 manifest의 execution kind가
  일치하는지도 확인한다. 따라서 manifest와 execution의 digest를 함께
  바꾸더라도 임의 JSON을 유효한 batch case report로 위장할 수 없다.
- 이 배치 schema hardening 이후 실행/입력/범위 표적 회귀는 **84 passed**,
  field/evidence-audit/semi-FV/obstacle-wake 전체 표적 묶음은 다시
  **396 passed** (0 failures)였다. `compileall`과 `git diff --check`도
  통과했다.
- 이후 저장소 전체 회귀(`pytest -q`)도 **1,561 passed, 145 skipped,
  0 failed** (20분 51초)로 완료했다. 기존 자유장·Fortran parity와 현장
  field 계층을 함께 포함한 현재 상태의 기준 수치다.
- 2026-10-07 최신 연결 폴더를 read-only `field-audit`로 재감사했다.
  다국어 헤더와 source-rate/wind near-miss 진단을 보존하도록 감사기를
  보강한 뒤 최신 v14 `SLABx_LH2` 감사는 **16,139**개
  구조화 파일을 스캔했고 그중 source-boundary 후보 15개와 33개가
  `observed_time_mean_vol_pct` 등 명시적 고정 수용기 후보 농도 열을
  노출했다.
  따라서 상태는 `partial`이지만 source/weather/obstacle/common-clock가
  없어 새 validation evidence로 승격하지 않았으며
  `promotion_allowed=false`를 유지한다. 지원센터 데이터는 20개 중 17개
  workbook candidate가 읽기 가능한 inventory로 남고 14개는 모호한 보조
  sheet 진단을 보존하지만, 필수 채널이 없어 `withheld`였다. 관련 자료
  0개도 `withheld`였다. 최신 v14 실행 artifact는
  `outputs/current-slabx-audit-2026-10-07-v14.json`,
  `outputs/current-support-center-audit-2026-10-07-v14.json`,
  `outputs/current-related-audit-2026-10-07-v14.json`이다. v2/v3/v10/v11/v12/v13 결과는
  이전 감사 로직·시점의 재현용 historical artifact로 보존한다.
- 요구사항별 현재 증거와 외부 승격 조건은
  [docs/field-goal-audit-2026-10-07.md](docs/field-goal-audit-2026-10-07.md)에
  표로 정리했다. 이는 완료/승격 선언이 아니라, 자유장 보존·범용 source·
  uncertainty·model comparison·보존형 obstacle·fail-safe gate와 현재
  `partial` external evidence를 동일 기준으로 재검토하기 위한 감사 기준점이다.
- 수동 사건/시계 reconciliation을 위해
  `FieldEvidenceManifest.from_audit()`와
  `write_field_evidence_manifest_json()`을 추가했다. 다섯 채널의 선택 경로,
  SHA-256, event/source/weather/geometry/operator ID를 모두 요구하고,
  읽기 시 파일 digest를 다시 확인한다. 기존 `FieldValidationEvidence`로
  변환할 수 있지만 `promotion_allowed=false`를 고정해 자동 승격은 하지 않는다.
  `field-validation-input.v1`도 선택적 `validation_manifest` 경로를 보존하고
  manifest evidence와 case evidence의 ID/path/SHA가 다르면 거부한다.
- 같은 manifest 경계를 `degali field-evidence-manifest-create` CLI로도
  사용할 수 있다. 검증된 `field-audit` execution artifact, 다섯 개의 명시적
  `--selected-path channel=path`, 사건/시계/식별자 metadata를 요구하며,
  불완전 감사·중복 channel·기존 출력 파일은 fail-closed한다. 생성 결과는
  기존 verify CLI와 동일하게 `promotion_allowed=false`를 유지한다.
- 위 manifest CLI와 linkage 회귀를 포함한 field/semi-FV 표적은
  **450 passed**다. CLI 추가 직전 마지막 완전 저장소 `pytest -q`는
  **1,619 passed, 145 skipped, 0 failed** (약 16분 15초)였으며, 새 CLI
  경로 자체는 전용 회귀에서 검증했다.
- Test 6 transient replay 도구와 기존 관측 연산자 회귀는 **7 passed**다.
  이 도구는 steady free-field 계산을 재사용하는 관측 연산자일 뿐이며,
  plume storage·gust deformation·time-varying source를 추가하지 않는다.
- manifest CLI와 Test 6 replay를 포함한 최신 완전 저장소 회귀는
  **1,622 passed, 145 skipped, 0 failed** (약 16분 26초)로 완료했다.
- batch execution verifier는 report schema뿐 아니라 manifest의 typed
  summary와 각 compact operational status/`screening_allowed` field도
  교차검증한다. 구형 screening artifact를 replay하지 못하는 경우에도
  `field_report`와 `operational_screening`의 허용 schema discriminator를
  확인해 임의 JSON을 integrity-only 결과로 수용하지 않는다.
- Python `write_field_screening_report()`와
  `write_field_refinement_report()`도 exclusive-create 경계를 사용한다.
  기존 report evidence를 암묵적으로 덮어쓰지 않으며, 같은 경로에 대한
  동시 writer race도 실패안전하게 거부한다. operational-envelope report,
  `manifest.json`, CLI `batch-execution.json`도 payload 선계산 후 같은
  exclusive-create 경계를 사용한다. batch/report/execution 표적 회귀는
  이 overwrite 방어를 포함해 **36 passed**이며, 이어진 전체
  field/evidence-audit/semi-FV/obstacle-wake 묶음도 **397 passed**였다.
- 나머지 field artifact writer도 같은 경쟁조건을 닫았다. model-sensor CSV,
  comparison/validation/source case JSON, evidence-audit JSON은 payload 또는
  fingerprint 검증 뒤 exclusive-create를 사용하고, FFI source-state,
  field-screen, comparison, validation, verification, field-source CLI 출력도
  동일한 exclusive text 경계를 사용한다. strict IO/CLI 집중 회귀는 이번
  보강 후 **171 passed**다.
- operational envelope는 이제 두 독립 distributed source의
  lower/nominal/upper schedule 9개 cross-product와 source별 injected-mass
  ledger residual을 함께 회귀한다. 이 추가 검증 뒤 전체
  field/evidence-audit/semi-FV/obstacle-wake 묶음은 **397 passed**다.
- batch export는 이제 case report와 manifest JSON을 모두 선직렬화한 뒤
  output directory를 만든다. 늦은 JSON/finite-value 오류가 앞선 report를
  남기는 partial batch를 만들지 않으며, 이 경계를 확인하는 render-failure
  회귀를 포함한 batch/report/execution 표적은 **37 passed**, 전체 field 묶음은
  **398 passed**다.
- operational screening decision record에 optional typed `gate_codes`를
  추가했다. evaluator는 physical blocked/conditional scope, transport
  incomplete, refinement missing/failed, sensor missing/trace missing,
  obstacle unrepresented, uncertainty unresolved, conditional review를
  안정적인 코드로 보존하며, 기존 직접 생성 decision은 빈 tuple 기본값으로
  호환된다.
- 같은 코드를 모든 operational aggregate가 corner union으로 보존하도록
  확장했다. conditional aggregate holdback에는 review code를 추가하고,
  구형 직접 생성 withheld corner의 빈 code는 `corner_withheld`로 표시한다.
  관련 operational/field 표적은 **398 passed**다.
- 이 aggregate 전파와 회귀 assertion을 포함한 저장소 전체 `pytest -q`는
  **1,563 passed, 145 skipped, 0 failed** (803.70초)로 완료했다. 기존
  자유장·Fortran parity와 현장 field 계층 모두 영향 없이 통과했다.
- `batch-execution.json`의 compact case 요약도 manifest와 같은 typed
  `gate_codes`를 반복하고, `field-verify`가 status/allowance와 함께 그
  목록을 교차검증한다. 따라서 자동화가 withheld 원인을 다시 prose에서
  추출하거나 manifest를 별도 재해석할 필요가 없다.
- `degali.field-batch-summary.v1`에도 case 전반의 결정론적
  `gate_code_counts`를 추가해, batch-level 자동화가 각 report를 다시 읽지
  않고도 보류 원인별 집계를 얻도록 했다.
- `field-verify`는 manifest case decision에서 summary를 다시 계산해
  status/allowance/gate-code count 위조를 거부한다. typed code 이전의
  legacy v3 manifest는 코드가 없는 경우에 한해 계속 integrity-check할 수
  있다.
- verifier의 package-level export 이후 reference parity를 다시 실행해
  **109 passed, 139 skipped** (0 failures)를 확인했다. 기존 자유장/Fortran
  parity 결과는 변경되지 않았다.
- 이 변경 후 `test_field_source_io.py`, `test_field_comparison.py`,
  `test_field_validation.py`, `test_field_historian_io.py`는 **62 passed**,
  field/evidence-audit/semi-FV/obstacle-wake 표적 묶음은 **397 passed**다.
- ambient-pressure saturation-table 범위 보강 targeted 테스트: **2 passed**.
- fingerprinted validation evidence conditional-gate targeted 테스트: **20 passed**.
- 현장 closure gate 이후 자유장 scope/model-comparison/screening 회귀:
  **23 passed** (2026-10-06).
- field batch/report/decision/workflow 교차 계층 회귀:
  **70 passed** (2026-10-06). off-plane obstacle projection과 withheld
  detector의 식별자가 case report와 batch manifest decision에서 일치하는지
  포함한다.
- operational uncertainty envelope/report/batch 회귀: **36 passed** (2026-10-06).
  새 `deterministic_sensor_envelope`는 각 detector의 peak/final/time-average
  true·indicated extrema를 completed corner에서만 계산하고, 누락/off-plane
  trace는 값 없이 withheld count와 사유를 보존한다. 실제 CLI 실행 결과는
  `outputs/current-field-screen-sensor-envelope.json`에 기록했다. 이 경로와 아래
  `outputs/` 실행 결과는 재생성 가능한 로컬 감사 산출물이며 공개 스냅샷에는
  포함하지 않는다.
- measured-history/source-sensor/phase-pool operational report targeted 회귀:
  **35 passed** (2026-10-06). 세 경로도 같은 deterministic sensor-extrema
  요약을 내보내며, phase corner에 screening trace가 없으면 aggregate에서
  제거하지 않고 `withheld`로 유지한다.
- standalone sensor-calibration operational report targeted 회귀:
  **17 passed** (2026-10-06). `--sensor-array-envelope`도 동일한 detector
  extrema/withheld summary를 보존한다. 실제 CLI 확인 결과는
  `outputs/current-field-screen-sensor-array.json`에 기록했다.
- 최신 nominal uncertainty-envelope CLI 재실행도
  `deterministic_sensor_envelope.status=complete`, `case_count=4`,
  `sensor_count=2`, `operational_screening=withheld`를 반환했다. 결과는
  `outputs/current-field-screen-envelope-latest.json`에 저장했다.
- atmospheric-source batch conditional-review 경로를 보강했다. 이제 named-case
  authorization이 source-rate/weather/sensor deterministic envelope에도 적용되고,
  모든 corner의 refinement·in-plane receptor·resolved-uncertainty gate를 통과한
  경우에만 manifest decision이 `conditional_allowed`가 된다. Python/strict CLI
  batch/strict-input 회귀는 **17 passed**로 확인했다.
- 같은 batch의 `manifest.json`과 `batch-execution.json`에
  `degali.field-batch-summary.v1` typed `summary`를
  추가했다. 이 값은 physical applicability가 아니라 실제 operational decision
  에서만 집계되며, case/status별 count와 `all_screening_allowed`를 보존해
  후속 자동화가 conditional/withheld를 허용 결과로 오인하지 않게 한다.
  직접 만든 `FieldBatchCase`도 `conditional_allowed` 상태라면 review record와
  `conditional_review_applied=True`를 함께 가져야 한다.
- `FieldBatchCase`의 applicability/status/refinement 메타데이터 일치성과
  `FieldBatchExport`의 manifest/report 파일명·경로·존재성·분리 조건도 생성 시
  검증한다. 따라서 허용 상태를 물리적 applicability만으로 위조하거나, 다른
  디렉터리/manifest를 가리키는 batch 산출물을 만들 수 없다.
- `FieldBatchExport`는 이제 실제 manifest/report JSON을 다시 읽어 schema,
  case 순서/이름, applicability/refinement/review, operational decision 및
  summary가 typed 결과와 일치하는지 교차검증한다. 변조된 기존 파일은
  경로가 유효해도 export 객체로 수용하지 않는다.
- manifest에는 각 case report의 SHA-256이 추가되고,
  `batch-execution.json`에는 manifest와 report hashes가 함께 기록된다. 후속
  자동화가 export 이후 파일 변조를 다시 검출할 수 있다.
- strict batch input도 직접 생성 시 case collection을 immutable tuple로
  요구하고, 파일을 읽기 전에 safe report-name label을 검증한다. 이 경계의
  회귀는 batch/IO 합계 **17 passed**로 확인했다.
- `ReleaseSource.metadata`와 stability scalar closure mapping도 immutable
  snapshot으로 고정한다. source-contract/meteorology/report 집중 검증은
  **82 passed**였다.
- phase/pool model-option mapping도 immutable snapshot으로 고정하며,
  phase-routing 집중 회귀에 이 경계를 포함한다.
- 직접 `ReleaseSource`의 위치/방향 불확실성도 각각 `m`/무차원 `1` 단위를
  생성 시 강제해 strict case-file 경계와 일치시킨다.
- read-only `field-audit`의 candidate/channel/reason record도 tuple·중복·상태
  일관성을 생성 시 검증한다. `partial`인데 탐지 채널이 없거나, complete scan이
  모든 채널을 찾았는데 `partial`로 남는 malformed readiness record는 거부된다.
  관련 audit 회귀는 **9 passed**다.
- 이 보강 후 지원센터 폴더를 실제로 재감사한 결과도
  `scanned_files=20`, `scan_complete=true`, `status=withheld`,
  `promotion_allowed=false`로 유지됐다. 결과는
  `outputs/current-support-center-audit-latest.json`에 저장했다.
- 마지막 확인 명령:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe -m pytest -q tests/test_lh2_scope_guard.py tests/test_lh2_yawed.py tests/test_model_comparison.py tests/test_screening_gate.py tests/test_validation_integration.py
.venv\Scripts\python.exe -m pytest -q tests/test_ffi_sensor_operator.py tests/test_ffi_source_state.py
.venv\Scripts\python.exe -m py_compile src\degali\lh2.py src\degali\__init__.py tools\audit_ffi_sensor_operator.py tests\test_ffi_sensor_operator.py
git diff --check
```

FFI 센서 연산자 감사:

```powershell
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe tools\audit_ffi_sensor_operator.py
```

이 도구는 로컬 `reference/spadeadam`을 읽지만 원자료를 출력물/배포물에 복사하지
않는다. `arc_max_table`과 `sensor_height_table`을 분리해 내보내며 raw sensor
row는 포함하지 않는다. 기대되는 핵심 결과는 `sensor_count=15`,
`projected_sensor_max_vol_pct≈7.50`, `observed_over_projected≈2.80`이고,
0.1/1.0/1.8 m 높이군 비는 각각 약 2.919/2.817/2.481이다.

## 4. “불가” 항목의 현재 처리 방식

| 원래 우려 | 현재 처리 |
|---|---|
| 모든 LH2 조건에서 보편적으로 정확 | `screening_scope`, strict guard, source/wind envelope, no-fit observation envelope로 조건·식별불능 구간을 수치화 |
| HyRAM/PHAST보다 일반적으로 정확 | 동일조건 비교 프로토콜을 코드화. 조건이 다르면 ranking 차단 |
| 센서 높이/arc 연산자 불일치 | `project_lh2_jet_to_sensors()`로 실제 좌표 투영. 시간동기/사행은 별도 입력 |
| 완전 비정상 plume | source snapshot + advection + 명시적 1차 응답 surrogate. transient puff라고 부르지 않음 |
| 장애물/벽/후류 | `site_geometry`가 접촉을 검출하고 free-plume 결과를 차단. 정량 wake 계수는 데이터 없이 발명하지 않음 |
| 설계·인허가 단독 판단 | screening gate는 조건부 screening만 허용하고 design/approval을 항상 거부 |
| N₂/O₂ 혼상 EOS | 순수성분 경계·opt-in 연구 경로로 제한. 검증 전 기본 경로 승격 금지 |

핵심은 “불가”를 숨기는 것이 아니라, 계산 가능 범위·연구 범위·사용 금지 범위를
각각 실행 가능한 API와 보고서로 바꾼 것이다.

## 5. 다음 작업 순서

1. **공개 repo용 diff 정리**: 원본 Fortran·reference·외부 raw data가 staged 되지
   않았는지 확인한다.
2. **관측 연산자 확장 완료**: 로컬 공개 FFI sensor 좌표를 사용해
   `arc_max_table`과 `sensor_height_table`을 별도 집계하고, raw table은
   커밋하지 않는다. source-state envelope에도 같은 집계를 corner별로
   연결했으며, lower-bound 위반/withheld 상태를 함께 보존한다.
3. **FFI source-state envelope 후속 검증**: source, wind, ambient bound를
   분리한 residual-map API는 구현됐다. 다음에는 공개 raw time-history 또는
   독립 source-state evidence가 확보될 때만 같은 경계에 추가하고, 현재
   deterministic corner 결과를 confidence interval이나 site acceptance로
   승격하지 않는다.
4. **E3.4/E3.5 회귀**: E3.5 baseline audit는 2026-10-06 재실행해
   197 concentration rows, 40 thermal rows, 6 FFI arcs와 MG 1.244826 /
   VG 1.373001 / FAC2 0.833333을 재현했다. E3.4 pool heat-balance는
   손상된/incomplete archive에서 원본 workbook을 안전하게 승격할 수 없어
   기존 fail-safe 경계를 유지하며, intact workbook과 독립 substrate evidence가
   확보될 때만 추가한다.
5. **장애물·혼상·난류**: 공개 데이터가 추가될 때만 정량 closure를 opt-in으로
   구현한다. 검증 전 기본값 변경 금지.
6. **문서/논문**: [MANUSCRIPT_LH2_VALIDATION.md](MANUSCRIPT_LH2_VALIDATION.md)에
   현재 수치와 적용범위를 반영한다.
7. **커밋/배포**: [docs/distribution-commit-plan-2026-09-19.md](docs/distribution-commit-plan-2026-09-19.md)
   절차대로 별도 검토 후 commit → push → GitHub release → Zenodo/PyPI 순서로 진행한다.

## 6. 다른 계정에서 바로 시작하는 방법

새 계정에서 이 폴더를 열고 다음 순서로 실행한다.

```powershell
Set-Location "<repository-root>"
Get-Content HANDOVER_CURRENT.md
git status --short
$env:PYTHONPATH = "src"
.venv\Scripts\python.exe -m pytest -q tests/test_lh2_scope_guard.py tests/test_model_comparison.py tests/test_screening_gate.py
```

새 계정에 `.venv`가 없으면 프로젝트 루트에서 다음을 실행한다.

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[test]"
```

CoolProp이 없는 환경에서는 LH2 물성 테스트가 skip될 수 있다. 임의로 외부
데이터를 다운로드해 저장하지 말고, 필요하면 공개 URL·해시·재현 명령만 문서에
추가한다.

## 7. 주의할 점

- `strict_scope=False` 기본값은 연구/민감도 계산용이다. 운영 screening은 strict
  guard와 `evaluate_screening()`을 함께 사용한다.
- FFI Test 6을 맞추기 위해 질량률 배수, flashing fraction, 난류계수, wake factor를
  임의로 넣지 않는다.
- 기존 `docs/HANDOVER.md`에는 과거 세션의 상세 로그가 많다. 최신 판단은 이 문서와
  `docs/model-hardening-2026-09-21.md`를 우선한다.
- GitHub/PyPI/Zenodo에 올릴 때 원본 Fortran과 외부 실험 원자료가 포함되지 않았는지
  wheel/sdist archive 목록을 확인한다.
