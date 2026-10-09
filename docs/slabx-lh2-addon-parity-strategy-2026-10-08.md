# SLABx–LH2 addon parity와 DEGALI 승부 전략

이번 검토는 SLABx 공개 addon과 FFI `model-comparison/addons`를 DEGALI의
현재 모듈과 대조한 결과다. 목표는 SLABx 코드를 흉내 내는 것이 아니라,
SLABx가 실제로 사용하는 source/impact/profile/operator 가설 중 독립적으로
검증 가능한 부분만 DEGALI의 보존식·evidence gate 안으로 가져오는 것이다.

## 추가한 모듈

### 1. Source-state ledger

`src/degali/addons/source_state_ledger.py`는 source-plane의 H2 질량유량,
carrier 질량분율, 온도, 밀도, 면적, 속도, 액체분율, 높이, 방향을 공통
시간축에 보존한다. `SourceStateLedger.to_source_history()`는 명시적인
zero endpoint와 gas-only 상태를 요구한다. 따라서 SLABx의 `source_ledger.py`와
같은 handoff 정보를 보존하면서도, 한 행짜리 P04 source-state나 액체 상태를
자동으로 q(t) 또는 대기 source로 승격하지 않는다.
해결된 gas-only ledger는 `to_source_rate_schedule()`로 기존
`semi_fv_obstacle.SourceRateSchedule`에 변환할 수 있어, source ID·시간축·질량
보존을 잃지 않고 장애물 수송 경로에 연결된다.

### 2. Post-impact wall-jet

`src/degali/addons/post_impact_walljet.py`는 SLABx의
`radial_walljet_field`, `bifurcated_walljet`, `radial_walljet_mixing`,
`source_state_vertical_profile`에 해당하는 보존형 operator를 제공한다.

- radial/bifurcated angular redistribution
- source/vector 방향을 반영한 von-Mises bias
- roughness·wind·jet speed에서 계산한 lateral wall-jet mixing
- ground-reflected vertical finite-volume profile
- 모든 연산의 scalar inventory 보존 검사

이 모듈은 FFI Test4의 수평 단일 jet에 자동 적용하지 않는다. downward/
impinging release에서만 별도 branch로 실행하고, Test4·Test6 또는 외부
impingement 자료에서 독립 검증한다.

## SLABx addon별 처리 결정

| SLABx 기능 | DEGALI 상태 | 전략 결정 |
|---|---|---|
| `source_ledger`, `physical_transition_cloud` | field source/phase/transport 계약이 이미 존재 | 이번 ledger로 공통 source-plane 기록을 보강 |
| `pool`, `lh2pool`, ground heat balance | `dynamic_pool`, `pool_evaporation`, `field_phase_routing` 존재 | 중복 구현하지 않고 source ledger와 연결 |
| `air_condensation`, `water_ice` | axisymmetric/cryogenic-air 경로에 이미 존재 | 유지; flammable range 밖 condensation은 별도 sensitivity로 보고 |
| `diagnostics`, `lfl` | field decision·validation gate가 부분적으로 담당 | 향후 applicability/LFL report로 통합 |
| `plume_width` | SLABx 전역 monkey-patch이며 NASA lofting 보정용 | DEGALI에 복사하지 않음; source/vertical profile 검증으로 대체 |
| `vertical_drag` | SLABx 자체 closure patch, 채택되지 않은 exploratory 경로 | 기본 모델에 추가하지 않음 |
| `depth_cap` | SLABx 내부 검증에서 두 번 거부됨 | 금지; dense-gas depth를 수동 cap하지 않음 |
| prediluted Test6 source | 다른 시험에서 옮긴 조건부 source hypothesis | FFI Test4 기본 경로에 사용하지 않음 |
| receptor calibration factor | 관측으로부터 source/scalar를 조정할 위험 | 입력 계약과 gate 밖에서 금지 |

## 새 경쟁 전략

1. **Source 경쟁**: 동일 source-state ledger에서 H2·carrier mass,
   momentum, energy, phase residual을 먼저 닫는다.
2. **Impact 경쟁**: horizontal FFI Test4와 downward/impinging FFI 시험을
   분리한다. 후자에서만 radial/bifurcated wall-jet branch를 허용한다.
3. **Profile 경쟁**: point sensor 비교가 아니라 수직 profile, ground
   contact, vertical finite-volume inventory를 함께 점수화한다.
4. **시간 경쟁**: source q(t), wind(t), receptor time series가 모두 있을
   때만 causal packet/time-response 점수를 발행한다. aggregate 275 s 행은
   conditional 진단으로만 남긴다.
5. **장애물 경쟁**: wall-jet/obstacle branch와 wind-plane semi-FV를
   결합하되, 3-D wake가 검증되지 않은 경우 `conditional` 또는 `withheld`로
   유지한다.

## 사전 채택/중지 규칙

- retained SLABx–LH2 research route와 비교할 때 MAE·RMSE·2/4 vol% FN/FP를
  동시에 보고한다. 하나의 지표만 낮아진 경우 우위로 부르지 않는다.
- Test4의 native SLABx 근소 우위는 개발진단일 뿐, 독립 검증 우위가 아니다.
- wall-jet operator가 총 scalar inventory, source momentum, vertical profile
  residual을 통과하지 못하면 센서 MAE 개선 여부와 관계없이 중지한다.
- FFI 자료에서 source q(t) 또는 receptor time series가 없으면 결과는
  `withheld`이며, SLABx와 DEGALI 어느 쪽의 시간응답 우위도 선언하지 않는다.
