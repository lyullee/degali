# TI021 10% 이하 오차 탐색 결과 (2026-10-15)

## 결론

현재 저장소의 TI021 자료와 관측값 비사용 물리 탐색만으로는 600 s 누적 dose
MAPE 10% 이하를 입증하지 못했다. 기존 최저 재현 결과는 42.8053%였고,
이번 추가 탐색에서 벽 차폐 운동량 손실을 포함한 물리 후보의 최저는
**40.8064%**(벽 높이 3 m, 층 항력 0.2, 층 혼입 0.09, 차폐 강도 0.02)였다.
풍속을 0.1 m/s까지 낮춘 비현실적 민감도는 40.0445%였지만, workbook의
3.4 m/s를 위반하므로 생산 후보가 아니다.
이는 10%와 여전히 큰 차이이며, 다음 설정을 함께 사용한 진단 결과다.

- 5 m 수평 격자, sensor-envelope domain
- 40 m 근거리 hand-off 반경
- 5 m 수직 혼합, `D_z=2 m²/s`
- diffusive vertical observation tail (10 m, 2 m²/s)
- 수평 확산 `D_h=1 m²/s`
- wake deficit 0.5, wake length 50 m, base width 10 m, slope 0.2,
  recirculation 0.5
- drag 0.05, entrainment 0.1
- wind 1.5 m/s, horizontal source-spread sensitivity 1.5

이 중 wind 1.5 m/s는 TI021 workbook의 3.4 m/s와 다르고, source-spread 1.5는
실제 방출 형상 측정값이 아니므로 **생산 기본값이나 논문 검증값으로 승격하지
않는다**. 입력 풍속을 workbook 값으로 유지한 같은 계열 결과는 약 45.23%였다.

## 원인으로 확인된 것

기존 50 m hand-off는 14--40 m 센서가 누출 초기 고농도 구간을 잃게 만들었다.
40 m로 당기면 47.67%까지 개선되므로, 근거리 source-to-layer overlap은 실제
지배 오차원인 중 하나다. 그러나 hand-off를 10--30 m로 더 당기면 오히려
53.48--48.80%로 악화되어, 전환 반경만 줄이는 것으로는 해결되지 않는다.

풍속 sweep은 1.0--4.5 m/s에서 43.71--46.75% 범위였고, 항력 sweep은
`C_d=0.05`에서 49.51%까지 개선했지만 0.5 이상에서는 악화됐다. 수평 확산은
`D_h=1`에서 약간의 개선만 있었고, 수직층 수를 늘리면 49--51%대로 악화됐다.
따라서 잔여 오차는 단일 계수보다 **방출 직후의 3-D gravity-current head,
벽 전면/상단/측면 분기, 그리고 벽 뒤 재순환의 시간 이력**이 빠진 문제로
해석해야 한다.

추가로 확인한 사항은 다음과 같다.

- 소스 전단·난류 혼입계수를 0--0.05로 낮추면 42.65% 아래로 내려가지 않았다.
  소스 질량·풍속을 임의로 바꾸는 방식도 40%대에서 정체했다.
- 풍향 부호를 workbook 상대풍향의 반대쪽으로 바꾸면 72.74%로 악화되어,
  현재 +6.1° 좌표 변환이 맞다는 것을 확인했다.
- 센서 축 교환·반전 진단도 수행했다. 공식 변환은 40.8842%였고, 축 교환은
  89.95%, 횡축 반전은 77.05%, 종축 반전은 95.96%로 악화됐다. 좌표계가
  현재 잔차의 주원인은 아니다.
- 벽 내부 차폐 운동량 손실은 질량 및 운동량 잔차를 각각
  `8.6e-12 kg`, `7.0e-10 kg m/s` 수준으로 유지하면서 약 1.9 percentage
  points만 개선했다. 수직 handoff-Gaussian 관측자는 83.76%로 악화되어
  채택하지 않았다.
- 공동 접선 와류(`vortex rate=0.03 s⁻¹`)는 41.92--46.67%였고, 보존형
  공동 저장·반환(`capture=0.01 s⁻¹, return=0.005 s⁻¹`)은 41.79%였다.
  둘 다 단순 차폐보다 개선되지 않았다. 저층 농축 `dense_power` 관측자는
  80.25% 이상으로 악화되어 질량보존형 후보로 채택하지 않았다.
- 초기 dose를 보존하기 위해 공동 저장 시작을 30 s로 지연하고 중심 코어로
  반환한 경우도 41.7166%에 그쳤다. 저장소 최종 잔류 종질량은 123.06 kg로
  보존됐지만, 관측 잔차를 줄이는 구조적 해법은 아니었다.

## 10%를 정직하게 달성하는 방법

다음 두 가지가 필요하다.

1. TI021 workbook의 요약 dose를 보정하지 않고, 센서별 0--600 s 농도 시계열과
   교정·시간동기 정보를 확보한다. 현재 자료는 양의 dose 84개와 arrival/departure
   요약만 있어 순간 농도 MAPE를 검증할 수 없다.
2. 방출원에서 0--50 m까지를 3-D 또는 검증된 moving-source/overlap 모델로
   계산하고, 그 출구 질량·운동량·높이 분포를 현재 2-D 장애물 패치에 보존적으로
   연결한다. 특히 벽 전면 충돌, overtopping, 측면 bypass, lee recirculation을
   별도 flux로 계측·검증해야 한다.

관측 dose를 센서별로 곱해 맞추거나, workbook 풍속과 다른 값을 기본값으로
고정하면 수치상 10%에 가까워질 수 있어도 검증 성능이 아니라 calibration
artifact이므로 채택하지 않는다.

## 재현 산출물

- `outputs/integral-local-patch-2026-10-15/dose-envelope5m-best-spread1p5.json`
  : 42.8053% 진단 최저값
- `outputs/integral-local-patch-2026-10-15/dose-envelope5m-h40-diff-tail-wind3p4.json`
  : workbook 풍속 3.4 m/s 계열 결과
- `tools/score_ti021_dense_layer_hybrid.py`
  : `handoff_radius_m`, `wind_speed_m_s`, `wind_angle_deg`를 명시적으로 기록
- `outputs/integral-local-patch-2026-10-15/wh-3.json`
  : 현재 물리 후보 최저 40.8763%; 보존 게이트 통과
- `outputs/integral-local-patch-2026-10-15/current-best.json`
  : 현재 코드 기준 재현 결과 40.8064%; 보존 게이트 통과
- `outputs/integral-local-patch-2026-10-15/current-workbook-wall6.json`
  : workbook 벽 높이 6 m를 유지한 재현 결과 40.8842%; 보존 게이트 통과
- `outputs/integral-local-patch-2026-10-15/shelter-0p02.json`
  : 차폐 민감도 40.9679%; 보존 게이트 통과
- `outputs/integral-local-patch-2026-10-15/wind0p1.json`
  : 40.0445% 비현실적 풍속 민감도(승격 금지)
- `outputs/integral-local-patch-2026-10-15/res01.json`, `v03s05.json`,
  `dense1.json`
  : 공동 저장·와류·저층 프로파일의 재현 결과
- `outputs/integral-local-patch-2026-10-15/layerwind05.json`
  : 방출 풍속 3.4 m/s, 층 운반풍속 0.5 m/s 분리 진단 42.58%; 개선 없음
- `outputs/integral-local-patch-2026-10-15/res-start30.json`
  : 30 s 지연·중심 반환 공동 저장 진단 41.7166%; 보존 게이트 통과

## 후속 재검증(2026-10-15)

상부 확산 tail을 유지하면서 지면 쪽을 질량보존형 power profile로 재분배하는
`grounded_diffusive_tail` 후보도 동일한 600 s, 5 m, 40 m hand-off 조건에서
검토했다. 지수 0, 0.25, 0.5, 1, 2의 MAPE는 각각 82.20, 81.72, 81.22,
80.26, 79.30%였고, tail 연장 40--125 m에서도 86.47--86.99%였다.
이는 저층 농축을 단순 관측연산자로 추가하면 상부 센서와의 균형이 더 나빠진다는
뜻이므로 구현을 폐기했다. 현재 저장소에는 이 폐기 후보를 생산 경로로 남기지
않았다.

따라서 현재 물리적·보존형 후보의 최저는 여전히 `current-best.json`의
40.8064% (벽 3 m 진단)이며, workbook 벽 6 m를 유지한 결과는
`current-workbook-wall6.json`의 40.8842%다. 10% 이하를 주장하려면 센서별
교정계수나 비현실적 풍속을 쓰지 않고, 방출 직후 3-D head/벽 충돌/후류 재순환의
시간이력과 동기화된 원시 농도 시계열을 추가해야 한다.

## 유한 prism 방출원 연결 재검증

Gaussian source 대신 저장소의 보존형 `integrate_inertial_prism`을 40 m
hand-off까지 연결하고, `prism_layer_inventory`로 질량·수평운동량을 패치에
전달했다. 동일한 workbook 벽 6 m와 600 s 평가에서 `grounded_head` source는
MAPE **40.9728%**, hand-off 5.35 s, species residual `4.55e-12 kg`,
momentum residual `6.98e-10 kg m/s`였다. Gaussian 재검증 40.8842%보다
개선되지 않았으므로 생산 후보로 승격하지 않는다. 다만 source-to-layer
인터페이스가 보존형으로 실행되는 것은 확인했으며, 다음 반복은 이 경로에
검증된 벽 충돌/후류 flux를 연결하는 단계다.

## 장애물 공동 flux 반복

finite-prism을 15 m에서 hand-off하고, 벽 후류의 일반적인 저감 길이와
상부 확산계수를 함께 스크리닝했다. 가장 낮은 재현값은
`prism-best-s015-r05.json`의 **38.0557%**였다. 이 후보는 600 s, 5 m
격자, workbook 풍속 3.4 m/s를 사용하고 보존 게이트를 통과했지만 10%와는
여전히 큰 차이가 있다. 공동 저장소를 중심 Gaussian 대신 반경 10 m ring으로
반환하는 보존형 flux(`prism-ring-r10.json`)는 **39.4468%**로 악화되어
폐기했다. 따라서 이번 반복에서 실제로 확인된 개선은 근거리 finite-prism
중첩과 짧은 후류 감쇠뿐이며, 목표 달성에는 아직 3-D 벽 충돌/재순환 시간장이
필요하다.

모든 후보는 관측 dose를 계산 중 보정하지 않았고, 질량·운동량 보존 게이트를
통과한 것만 비교했다. 따라서 현재 목표 상태는 “10% 이하 달성”이 아니라
“10% 이하를 위해 어떤 물리와 데이터가 추가로 필요한지 정량적으로 특정”한
상태다.

## 2026-10-15 후속 민감도 감사

저층 센서 과소예측의 원인을 분리하기 위해 공동 내부에만 적용되는 보존형
관측 수직구배(`cavity_vertical_exponent`)를 추가했다. 공동 외부의 농도와
보존변수는 변경하지 않고, 공동 내부의 column-mean을 지수형 ground-hugging
분포로 재분배하는 폐쇄다. 계수 0.25–6의 범위에서 MAPE는 37.04–36.71%로
소폭만 개선되었다. 수직 다층 transport(3·6 layers)는 65.13–141.11%로
악화되어 생산 후보에서 제외했다.

공동 차폐·후류 길이·재순환을 함께 조정한 최저 재현값은
`outputs/integral-local-patch-2026-10-15/prism-bestm1p05.json`의
**34.4623%**이다. 조건은 600 s, 5 m, finite-prism, hand-off 15 m,
source mass sensitivity 1.05, cavity shelter 0.10, cavity vertical exponent
4, cavity entrainment factor 0.70, wake length 20 m, wake recirculation 0.70이며,
실행시간 16.28 s, species residual `6.82e-12 kg`, momentum residual
`6.98e-10 kg m/s`, momentum gate `true`다. source mass 1.05는 검증된
방출량이 아니므로 field-qualified 설정으로 승격하지 않고 민감도 최저값으로만
기록한다.

따라서 10% 이하 목표는 아직 달성되지 않았다. 2.5 m 격자는 239.43 s에
39.8762%로 악화되었고, 초기 형상·수직 프로파일·다층 transport·reservoir
반환·vortex·풍향 변형은 10%에 접근시키지 못했다. 잔여 오차는 근거리 공동
센서(반경 0–60 m)에서 약 49–55%로 집중되며, 현재 단일층 FV가 재현하지 못하는
3-D fence-top separation, 재순환 체류시간, 시간해상도 원시 농도자료의 부재가
핵심 차단점이다. 센서별 보정계수나 관측 dose를 입력으로 쓰는 방식은 목표에서
제외한다.

## 2026-10-15 반복 루프: 후류·수직·확산 후보

동일한 600 s 전체 누적 dose 평가로 다음 후보를 연속 검증했다.

| 후보 | MAPE | 판정 |
|---|---:|---|
| backflow carrier, cavity entrainment 0.80, tail diffusivity 0.60 | 32.8688% | 기준 개선, 보존 게이트 통과 |
| tail diffusivity 0.28, cavity exponent 3, wake length 27.5 m, wake recirculation 0.40 | **31.20%** | 현재 민감도 최저, field-qualified 아님 |
| conserved cavity momentum backflow rate 0.02--1.0 s⁻¹ | 33.06--39.52% | 악화, 폐기 |
| mass-conserving diffusive tail | 71.3569% | 악화, 폐기 |
| cavity stagnation carrier 0.10--0.90 | 33.06--36.56% | 악화, 폐기 |
| anisotropic crosswind diffusion factor 0.5--4 | 32.96--42.35% | 악화 또는 positivity 한계, 폐기 |

현재 최저 산출물은 `outputs/integral-local-patch-2026-10-15/prism-finDp275.json`이며,
source mass scale 1.05를 사용한 진단 민감도 결과다. 질량 잔차와 운동량 예산 게이트는
통과했지만, source mass scale 자체가 독립 측정값이 아니므로 생산 설정으로 승격하지
않는다. 동일 물리에서 strict source scale 1.0을 고정한
`outputs/integral-local-patch-2026-10-15/prism-finDp275-strict.json`은 MAPE
31.8574%, runtime 14.24 s, species residual `5e-12 kg`, momentum residual
`9.31e-10 kg m/s`, gate `true`다.

이 반복으로 확인된 결론은 계수 하나의 조정으로 10%에 도달할 수 없다는 것이다. 특히
전면 정체점과 상부/측면 bypass를 하나의 원형 cavity로 대체하면 S004·S005·S006
수직열과 75--100 m 상부열의 오차 방향이 서로 달라진다. 다음 채택 후보는 센서별
보정이 아니라, 독립적인 벽 전면/상단/측면 flux와 시간분해 후류 체류자료를 가진
3-D 중첩·출구 인터페이스여야 한다.

## 수직영역·방출 중첩 재검증

다층 solver가 센서 최고 높이 6.4 m를 수직 상한으로 사용하던 문제를 수정해,
finite-prism 방출원의 실제 상단 13 m를 수직영역에 포함했다. 보존 게이트는 모두
통과했지만 2--6층 후보의 최저 MAPE는 45.18--50.95%로 단일층 후보보다 높아
채택하지 않았다. 이는 기존 다층 악화가 단순히 수직영역 상한 하나의 문제만은
아님을 뜻한다.

finite-prism 수평 중첩 폭도 질량·종질량·운동량을 각 성분별로 재정규화하는 보존형
Gaussian overlap으로 시험했다. strict source scale 1.0에서 spread 1.10--3.0의
MAPE는 31.71--35.09%로 기준 31.86%보다 개선되지 않았다. 따라서 방출 폭을
임의로 넓히는 것은 승격하지 않았다.

## 월류 질량분포 보완 및 실행 인터페이스 재검증

단일층 벽 flux의 상단 통과율을 균일 column으로만 계산하던 항을 보완했다.
`wall_overtop_exponent=0`은 기존식이며, 양의 값은 지면에 집중된 수직 질량분포를
벽 상단 통과율에 적용한다. 동일한 600 s, 5 m, finite-prism, hand-off 15 m,
workbook 풍속 3.4 m/s 및 strict source scale 1.0 조건에서 결과는 다음과 같다.

| wall overtopping exponent | dose MAPE | 보존 게이트 |
|---:|---:|---|
| 0 | 32.03% | pass |
| 0.5 | 31.52% | pass |
| 3 | 31.13% | pass |
| 10 | 31.03% | pass |
| 30 | **31.02%** | pass |

최저 산출물은 `outputs/integral-local-patch-2026-10-15/prism-overtopx30.json`이다.
그러나 지수 10--30은 현재 단일층 식에서 사실상 벽 상단 월류를 거의 차단하는
진단 한계값이며, 등록된 6 m 장애물의 실제 3-D 상단/측면 분리를 검증한 값이
아니다. 따라서 코드의 기본값은 0으로 유지하고, 이 항은 물리 민감도 증거로만
보관한다. `vertical_grounding_exponent`와 함께 최적화한 조합의 최저값은
`prism-combov0p1.json`의 31.01%였으며 10% 목표와는 여전히 21.01%p 차이가
있다.

finite-prism handoff의 방출 폭을 풍향축과 횡풍축으로 분리하는 보존형
재매핑도 추가했다. strict source scale 1.0에서 횡풍 폭 1.5배·풍향 폭 1.0배의
`prism-anisoL1C1p5.json`이 **30.86%**, 횡풍 1.8배·풍향 1.0배가 31.19%였다.
풍향 폭 0.8배·횡풍 1.5배는 30.91%로, 모두 mass/momentum gate를 통과했지만
10% 이하와는 큰 차이가 있다. 이는 일부 저층·측면 잔차가 초기 source footprint
폭에 민감하다는 원인은 지지하지만, 손실의 주원인이 근거리 폭 하나는 아님을
확인한다.

후류 체류시간을 방출 횡풍폭 1.5배와 결합해 재탐색했다. `wake_length=32 m`,
`wake_recirculation=0.3`, `wake_deficit=0.8`의 결과는 wall overtopping 지수
3에서 31.02%, 지수 30에서 30.61%였다. 지수 30은 사실상 완전 차단 한계값이므로
진단 최저일 뿐 채택하지 않는다. 지수 3의 보존형 조합도 10% 목표와 21.02%p
차이가 남아 후류 단일 매개변수로는 목표를 달성할 수 없음을 확인했다.

지면층과 상부 월류층을 혼합하는 2-성분 수직 관측 폐쇄도 시험했다. 상부 모드
fraction 0.02--0.50의 MAPE는 30.80--35.94%였고, 기준 조합(0 fraction)의
30.75%보다 낮아지지 않았다. 상부 성분을 임의로 늘리는 방식은 폐기하고, 실제
상단 flux·수직운동량을 제공하는 3-D 인터페이스가 필요하다는 결론을 유지한다.

wall overtopping exponent를 수직 cavity 지수와 같은 3으로 제한한
비등방성 조합(`prism-aniso-wall3.json`)은 31.02%였다. 반면 기본 wall
closure(지수 0)의 같은 조합은 32.26%로, 상단 월류 폐쇄가 점수에 영향을 주는
것은 확인되지만 10%에 접근시키는 지배항은 아니다.

풍향·handoff·반경을 독립적으로 재검증한 결과도 31%대에서 멈췄다. 특히 CLI가
새 수직지수 인자를 실제 `run()`에 전달하는지 확인했고, focused test 50개와
Python compile, `git diff --check`를 통과시켰다. 따라서 이후 반복은 동일한
600 s 비교와 물리 보존 게이트를 유지하며, 센서별 보정이 아닌 벽 전면/상단/측면
flux와 후류 체류시간을 해상하는 3-D 인터페이스 추가가 되어야 10% 목표를
정직하게 재시도할 수 있다.

## 3-D solid-mask interface witness

기존 3-D dense-gas solver에 내부 solid-mask face blocking을 연결했다. 인접한
solid cell을 공유하는 면의 scalar flux를 0으로 만들고, 확산도 solid 경계에서
Neumann 처리하여 질량 손실 없이 전면/상단/측면 우회가 가능한 기반을 마련했다.
보존 검증은 focused test 60개 전체 통과와 함께 확인했다.

TI021 600 s diagnostic replay는 5 m 수평 격자, 2 m 수직 격자, 102×124×10
cells로 실행했다.

| 경로 | MAPE | 최대 질량 잔차 | 실행 결과 |
|---|---:|---:|---|
| 3-D scalar R12, 장애물 없음 | 94.891% | `7.28e-12 kg` | 기준 witness |
| 3-D scalar R12, 반원 solid mask | 94.160% | `4.09e-12 kg` | obstacle delta만 확인 |

산출물은 `outputs/integral-local-patch-2026-10-15/ti021-3d-no-obstacle-v2.json`와
`ti021-3d-obstacle-v3.json`이다. 장애물 경로는 28.38 s에 실행되며, 양의 센서
dose 합은 1,935.73 대 관측 17,317.13 vol%·s로 크게 과소예측된다. 장애물 효과는
보존되지만 scalar R12 소스와
layer solver의 재결합·농도 변환이 없어 MAPE가 10% 목표와 크게 어긋난다. 따라서
이 경로는 최종 모델로 승격하지 않고, 다음 구현의 필수 인터페이스(3-D solid
flux → layer/outer handoff)와 검증 증거로만 유지한다.

## 공동(cavity) 형상 분기와 체류 민감도

반원형 장애물의 공동을 전체 원판으로 처리하던 항을 분리해 `full_disk`,
`downwind_semidisk`, `upwind_semidisk` 세 형상을 동일한 600 s 조건으로
비교했다. 조건은 5 m 격자, 0.25 s 층 시간간격, finite-prism 방출,
횡풍 source 폭 1.5배, 후류 길이 32 m, 재순환 0.3, wall overtopping 지수
3이었다.

| cavity geometry | dose MAPE | 실행시간 | 보존 게이트 |
|---|---:|---:|---|
| `full_disk` | **30.75%** | 17.30 s | pass |
| `downwind_semidisk` | 31.35% | 16.42 s | pass |
| `upwind_semidisk` | 40.98% | 16.75 s | pass |

하류·상류 반원 모두 전체 원판보다 악화됐으므로 어느 반원도 승격하지 않았다.
결과 파일은 `outputs/integral-local-patch-2026-10-15/` 아래의
`prism-cavity-*.json`이다.

같은 조건의 빠른 600 s 민감도(0.5 s 시간간격)에서 수직 ground exponent
0→1→3→10→30은 MAPE 31.73→41.55→67.05→83.39→65.08%로 악화했고,
entrainment coefficient 0.05/0.09/0.15는 79.83/31.73/51.20%였다.
수평 확산계수 0/0.25/0.5/1/2 m²/s는 33.98/33.08/32.46/31.73/35.75%였다.
따라서 근거리 저층 부족을 단순 수직 집중이나 장기 저장으로 보정할 수 없으며,
현재 잔차의 주원인은 센서별 보정이 아니라 3-D obstacle separation과 지속
노출을 함께 전달하는 방출-후류 인터페이스로 남는다. 모든 후보는 질량 및
운동량 게이트를 통과했지만 10% 목표는 아직 미달이다.

풍향 부호도 동일 조건에서 확인했다. `+6.1°`는 31.73%였지만 `0°`, `+15.8°`,
`-6.1°`, `-15.8°`는 각각 48.46%, 44.39%, 74.63%, 99.67%였다. 따라서
문서의 상대각 부호를 기계적으로 뒤집는 것은 오차 원인이 아니며, 현재 좌표계의
`+6.1°`를 유지한다.

공동 내부 유속 저감(`cavity_shelter_strength`)도 0.1/0.3/0.5/0.7/0.9를
검사했다. MAPE는 31.73/34.91/39.72/45.93/74.41%로 단조 악화했다.
따라서 관측된 지속 노출을 단순히 공동 유속을 줄이는 방식으로 재현하는 것은
채택하지 않았다. source mass scale 1.0--1.5도 최저 31.72%에서 멈췄다.
현재 30.75% 결과는 0.25 s 정밀 시간간격 결과이며, 빠른 0.5 s 스캔값과
혼동하지 않는다.

후류 역류 속도를 보존변수 운동량에 직접 완화하는
`cavity_backflow_rate_s_inv`를 처음으로 활성화했다. 빠른 스캔에서
0/0.005/0.01/0.02/0.05/0.1 s⁻¹은 31.73/31.68/31.47/31.81/33.20/34.86%였다.
0.01 s⁻¹을 동일한 600 s, 5 m, 0.25 s 조건으로 재계산한 결과는
31.01%, species residual `4.09e-12 kg`, momentum residual
`6.98e-10 kg m/s`로 게이트를 통과했다. 다만 벽 월류 지수 30을 함께 쓴
진단 한계값은 30.59%에 그쳤고, 지수 30은 사실상 완전 차단이므로 승격하지
않았다. 역류 운동량 항은 실제 개선 신호가 있지만 10% 목표를 해결하는
지배 물리는 아니다.

공동 geometry 선택과 보존변수 역류율을 고정하는 회귀 테스트 2개를 추가했으며,
focused test 묶음은 **62 passed**로 갱신됐다.

방출 직후의 상부·측면 연행을 분리하는 인터페이스도 추가해 CLI로 재현 가능하게
했다. 기준 `.1/.1` 대비 상부 0·측면 .1은 31.59%, 상부 .1·측면 0은
31.69%, 양쪽 .05는 31.68%였다(빠른 600 s scan). 따라서 source 초기
연행의 면별 분리는 잔차를 설명하지 못해 기본값(None=기존 공통계수)을 유지했다.

## 2026-10-15 방출 지속시간 및 층풍 재검증

입력 workbook의 `release type=puff`, `release rate=4962 kg/s`, `release
duration=1 s`를 확인하고, finite-prism 방출 적분기에 `source_release_duration_s`
를 연결했다. 관측 dose를 사용한 보정은 하지 않고, 나머지 조건을 고정한 600 s
전체 비교를 수행했다.

| 후보 | coarse 600 s dose MAPE | hand-off | 판정 |
|---|---:|---:|---|
| 즉시 방출(None) | 31.47% | 2.60 s | 기준 |
| workbook 1 s 방출 | 32.99% | 3.60 s | 악화·폐기 |
| 2 s 방출 | 35.37% | 4.60 s | 악화·폐기 |
| 5 s 방출 | 47.35% | 6.75 s | 악화·폐기 |

따라서 방출시간 누락은 현재 30%대 잔차의 주원인이 아니다. 입력 조건을
그대로 반영한 strict 재계산에서 층풍 0.5 m/s, source scale 1.0 후보는
`outputs/strict-layer05-m1.json`의 **30.74%**였고, 종질량 잔차
`5.00e-12 kg`로 보존 게이트를 통과했다. source scale 1.1은 31.72%로
악화되어 field-qualified 설정으로 사용하지 않는다.

이번 반복에서 확인된 것은 방출시간·초기 방출속도·초기 면별 연행이 아니라,
센서 위치별 공간 분포를 결정하는 장애물 전면/상단/측면 분리와 시간분해
후류 체류장이 여전히 지배적인 미해결 항이라는 점이다. 따라서 10% 이하 목표는
아직 달성되지 않았으며, 다음 후보는 단일층 계수 조정이 아닌 보존형 3-D
obstacle-to-layer flux 인터페이스여야 한다.

## 추가 원인탐색: hand-off·운반풍·연행 경로

동일한 coarse 600 s 평가에서 hand-off 반경 5/10/15/20/30/40/50 m의 MAPE는
32.45/32.09/31.47/32.07/33.51/36.92/41.03%였다. source wind scale
0.1--1.0은 31.87--31.47% 범위였고, 층 운반풍 0.5--4.0 m/s도 최저
31.32% 수준에 머물렀다. 즉 시간 원점을 조금 이동하거나 source/층 평균풍을
재지정하는 것만으로는 공간별 dose 잔차가 해소되지 않는다.

상수 연행 대신 `parker_bulk_suppression`을 적용한 후보는 농도장이 과도하게
보존되어 MAPE가 451--10,568%까지 발산했다. 이는 현재 단일층 보존량에
Parker 억제를 바로 삽입할 수 없다는 구현·물리 불일치의 증거이며 생산 경로에
승격하지 않았다. 다음 구현은 계수 재조정이 아니라, 장애물 전면·상단·측면
flux를 높이별로 분리하여 후류/상부 bypass에 전달하는 보존형 3-D coupling으로
제한한다.

## 3-D fast obstacle witness 재연결

기존 fast-3-D potential-flow/scalar solver를 TI021 workbook 센서 145개, R12
species inventory `3290.12 kg`, 1 s source schedule, 180° zero-thickness arc에
연결한 독립 witness를 실행했다. 5×5×2 m, 600 s, `dt=0.1 s` 조건에서
실행시간은 **100.53 s**, 최대 질량 잔차는 `1.96e-11 kg`였지만 dose MAPE는
**86.98%**, 예측 양성 dose 합은 `1197.63` 대 관측 `17317.13 vol%·s`였다.

따라서 장애물 mask와 측면 bypass만 추가한 passive 3-D scalar route도 정확도
해결책이 아니다. 원인은 source-to-obstacle 근거리 binary mixture/부력과
분리 후류의 부재로 좁혀졌으며, 이 witness는 보존·geometry 경로 증거로만
남기고 생산 후보로 승격하지 않는다. 결과는
초기 prism에서 계산한 총 질량 `4964.972 kg`, R12 질량 `3290.116 kg`을
사용했으며, 결과는 `outputs/ti021-fast-3d-witness-v2.json`에 저장했다.

기존 variable-density 3-D diagnostic도 source duration 1 s, source sigma 7 m로
맞춰 재실행했다. `outputs/ti021-3d-obstacle-release1-sigma7.json`은 MAPE
**95.49%**, 예측 양성 dose 합 `2302.14`, 질량 잔차 `7.73e-12 kg`, runtime
`33.54 s`를 기록했다. 따라서 source duration과 buoyancy를 켠 것만으로는
해결되지 않으며, 이 경로 역시 binary-mixture/obstacle-to-layer 재결합 없이
생산 후보가 될 수 없다.

## 2026-10-15 높이분해 벽면 flux 및 접선 우회 재검증

단일 높이 평균으로 장애물 월류를 계산하는 것이 잔차의 원인인지 확인하기 위해
보존변수 수직층 5개에서 벽면 위 질량을 직접 집계하는
`height_resolved_wall_flux`를 추가했다. 600 s 전체 계산에서 초기 Gaussian 분포는
MAPE **58.72%**, 균일 수직 분포는 **66.19%**였다. 두 경우 모두 종질량 잔차는
각각 `4.55e-12 kg`, `4.55e-12 kg`이고 운동량 예산 게이트도 통과했지만, 기준
30.74%보다 악화되어 생산 경로에서 제외했다.

벽면 충돌 질량의 일부를 국소 fence tangent 방향으로 재분배하는
`wall_tangential_bypass_fraction`도 추가했다. strict `dt=0.25 s`에서 fraction
0.5는 MAPE **31.23%**, 종질량 잔차 `1.82e-12 kg`, 운동량 잔차
`9.31e-10 kg m/s`였다. coarse `dt=0.5 s` scan의 fraction 0.1/0.2/0.3도
각각 **31.62/31.70/31.56%**로 개선되지 않았다. 이 항은 보존성은 확인했지만
센서별 공간분포를 교정하지 못했으므로 기본값 0으로 두고 후보를 폐기했다.

이번 반복으로 확인된 결론은 높이 평균을 높이분해하거나 벽면을 따라 우회시키는
국소 closure만으로는 10% 목표에 접근하지 못한다는 것이다. 다음 후보는 3-D
장애물의 전면·상단·측면 경계 flux를 시간·높이별 packet으로 내보내고, 이를
동일한 질량·종·운동량으로 dense-layer 셀에 투영하는 실제
`obstacle-to-layer` 인터페이스로 제한한다. 새 인터페이스도 반드시 600 s 전체
MAPE와 두 보존 게이트를 함께 재계산해 채택 여부를 판정한다.

## 수평 확산 민감도 재확인

기준 strict 설정을 유지하고 `dt=0.5 s`의 빠른 600 s 전체 스크리닝으로
`horizontal_diffusivity_m2_s`를 0/0.25/0.5/1/2/4로 바꿨다. MAPE는
**33.79/32.57/32.04/31.37/34.55/42.46%**였고, 각 실행의 종질량 잔차는
`6.4e-12 kg` 이하, 운동량 게이트는 모두 통과했다. 수평 확산 1 부근이 가장
낮지만 strict 기준 30.74%보다 좋아지지 않았으므로 새 물리 후보로 승격하지
않았다. 이 반복은 잔차가 단순한 균일 수평 확산계수 부족이 아니라 장애물의
높이·면별 분리와 후류 재진입 경로에 있다는 판단을 지지한다.

## cavity 재순환 운동량 closure 스캔

장애물 후면 체류시간을 직접 나타내는 cavity 접선속도 목표장을 추가하고,
`cavity_vortex_rate_s_inv`/`cavity_vortex_strength`를 `(0.01,0.2)`,
`(0.01,0.5)`, `(0.05,0.2)`, `(0.05,0.5)`로 바꿔 600 s 전체 coarse scan을
수행했다. MAPE는 각각 **31.53/31.52/32.05/32.27%**였고 모든 경우
종질량·운동량 게이트를 통과했다. strict 기준 30.74%보다 낮아지지 않아
재순환 closure도 채택하지 않았다. 운동량을 보존하면서 회전시키는 것만으로는
센서별 dose 공간분포를 복원할 수 없으며, 실제 분리 shear layer와 높이별
경계 flux가 필요하다는 원인 판단을 유지한다.

## 관측 window 연산자 점검

워크북에 있는 센서별 arrival/departure 라벨을 observation operator로만
적용하는 후보도 별도로 계산했다. 두 라벨이 모두 있는 67개 센서는 해당
구간만 적분하고, 나머지는 기존 0--600 s 적분을 유지했다. 기준 trace의
전체 적분 MAPE는 **30.74%**였지만 window 혼합 적분은 **41.17%**로
악화되었고, 라벨이 모두 있는 subset에서도 44.86%였다(동일 subset의
전체 적분은 31.79%). 따라서 현재의 30.74% 잔차는 단순히 보고된 센서
window를 잘못 적분해서 생긴 것이 아니다. arrival/departure 값은 향후
동기화된 `FieldEvidenceManifest`의 시간축 검증에는 사용하되, 이번 MAPE
보정에는 채택하지 않았다.

## 보존형 obstacle-to-layer 계약 검증

`src/degali/addons/obstacle_layer_coupling.py`에 장애물 경계의 높이·횡방향
flux를 시간 packet으로 보존하는 `HeightResolvedBoundaryPacket`과 계층 수신용
collapse 연산을 추가했다. 종질량과 운동량을 packet 단위로 적분한 뒤 원래
rate history와 비교하는 테스트를 통과했으며, 전체 핵심 테스트는 **81 passed**다.
이 인터페이스를 추가한 상태에서 기존 strict 경로를 600 s, 5 m, `dt=0.25 s`로
재실행한 `outputs/strict-layer05-m1-interface-baseline.json`은 MAPE
**30.7435%**, 종질량 잔차 `5.00e-12 kg`, 운동량 예산 잔차
`9.31e-10 kg m/s`로 기준과 동일했다. 즉 packet 계약은 기존 결과를 조용히
바꾸지 않으며, 실제 3-D obstacle solver의 높이별 total/species/momentum
history를 연결할 준비만 된 상태다. 아직 이 계약만으로 10% 달성을 주장하지
않으며, 다음 반복에서 3-D 경계 history를 실제 dense-layer 수신면에 연결한다.

## 3-D 적분 patch의 전면·측면·상단 packet 외부 전달 재검증

`solve_integral_obstacle_patch`가 내보내는 `x_min/x_max/y_min/y_max/z_min/z_max`
종 flux를 모두 높이별로 외부 센서까지 전달하는 `all_boundary_trace`를 연결했다.
지면 반사 수직 kernel과 기존 x_max 전용 연산자의 일치 회귀도 추가했다. TI021
600 s, 5 m, flow 1000 조건에서 실행시간은 **26.87 s**로 30 s 목표를 통과했고
patch 질량 잔차는 `1.05e-11 kg`였다. 그러나 MAPE는 **97.94%**(기존 해당
passive integral 경로 94%대)로 악화되었고 carrier도 `flow_converged=false`였다.
따라서 모든 면을 외부로 계속 보내는 것만으로는 dense-layer의 30.74% 기준이나
10% 목표에 접근하지 못한다. 이 후보는 보존·실행시간은 확인했지만 정량 경로로
승격하지 않고 폐기했다.

## 실제 반원 fence cavity 형상 승격

TI021 장애물은 180° 후방 반원 fence이므로 `full_disk` cavity는 upstream
영역까지 재순환 영역으로 취급하는 형상 오류가 있었다. `downwind_semidisk`로
바꾸어 동일 strict 조건(5 m, `dt=0.25 s`, 600 s, source scale 1.0)을 다시
계산한 `outputs/ti021-cavity-geometry-downwind-semidisk-strict.json`은 MAPE
**30.5011%**를 기록했다. 기준 full-disk 30.7435% 대비 0.2424%p 개선이며,
종질량 잔차 `5.91e-12 kg`, 운동량 잔차 `6.98e-10 kg m/s`로 두 보존 게이트를
통과했다. 따라서 scorer 기본 형상을 downwind semidisk로 변경했다. 이 개선은
형상 오류를 제거한 결과이지 관측값 보정이 아니지만, 10% 목표와는 아직
20.50%p 차이가 남아 있다.

## 반원 cavity 역류 closure strict 승격

반원 형상을 고정한 뒤 `cavity_backflow_strength`/`rate` coarse scan에서
최저였던 `(0.9, 0.01 s^-1)`을 strict로 재계산했다. 결과
`outputs/ti021-semidisk-backflow-s0p9-r0p01-strict.json`은 MAPE
**30.2846%**, 종질량 잔차 `8.19e-12 kg`, 운동량 잔차
`9.31e-10 kg m/s`였다. downwind semidisk 기준 30.5011%보다 0.2164%p
추가 개선되어 scorer 기본값을 이 조합으로 승격했다. 이는 후류 역류 속도의
보존형 완화 closure이며 센서 dose를 직접 보정하지 않는다. 현재 best도
10% 목표에는 아직 20.28%p가 남아 있다.

## 후류 cavity 유효 반경 민감도

실제 fence 반경 50 m와 별도로 closure의 후류 유효 반경을 30/40/60/75 m로
coarse 600 s 스캔했다. MAPE는 **37.68/33.41/31.24/44.08%**로 모두 현재
strict best 30.2846%보다 높았다. 따라서 cavity 반경을 관측값에 맞춰 넓히거나
줄이는 방식은 채택하지 않았고, 실제 geometry 50 m와 분리 후류의 시간·높이
분포를 직접 해상하는 다음 단계가 필요하다.

## 반복 원인 탐색: 수직·질량·핸드오프·후류 후보 재검증

현재 strict best 조합을 기준으로 600 s 전체를 `dt=0.5 s`로 빠르게 재계산했다.
수직 cavity 분포 지수 1/2/4/6의 MAPE는 **31.87/30.60/30.52/31.68%**였고,
기준 지수 3보다 낮아지지 않았다. 방출 질량 스케일 0.85/1.00/1.15는
**32.36/30.05/32.90%**였으므로, 총량 보정은 공통 저편향의 원인이 아니며
채택하지 않았다.

핸드오프 반경 5/10/20/30 m는 **31.71/30.92/30.96/33.21%**로 모두 악화했다.
도착시각 감사에서 관측 유효 이송속도는 2.86 m/s였지만, 수송속도 2.4/2.8 m/s
민감도도 **30.22/30.15%**에 그쳐 기상풍을 단순히 낮추는 해법은 아니었다.
상층 월류 성분 0.05/0.10/0.20은 **30.26/30.58/31.30%**, 월류 개방도
지수 0.5/1/2는 **30.46/30.26/30.11%**로 개선되지 않았다.
따라서 현재 30.28% 잔차의 지배 원인은 단일 질량·시간·수직 보정이 아니라
장애물 전면 분리 shear layer와 후류 재진입의 3-D 결합 누락으로 판단한다.

추가로 관측 질량 proxy 방향의 source scale 1.3/1.5/1.75/2.0을 검사했을 때 MAPE는
**37.62/46.72/59.89/74.58%**로 급격히 악화했다. 따라서 현재의 28%대 저편향을
source inventory를 키워 해소하는 것은 잘못된 방향이다. 보고서의 `4962 kg/s × 1 s`
필드는 검증된 방출 질량 이력이 아니므로 생산 입력으로 승격하지 않는다.

벽 높이 불확실성 4/5/6/7 m도 **30.19/30.10/30.05/30.02%**로 변화가 작았다.
즉 30.28% 잔차는 단순한 fence 높이 입력 오차로 설명되지 않는다.

## 반원 fence 끝단 분리 우회 closure

crosswind 외곽 센서에서 상대오차가 커서, 기존 전 벽면 tangent bypass와 분리해
반원 양끝단에만 접선 우회를 허용하는 `wall_tip_bypass_fraction`을 추가했다.
coarse fraction 0.03/0.06/0.10은 **30.20/30.04/30.20%**였으나, 개선처럼 보인
0.06을 동일 600 s strict(`dt=0.25 s`)로 재계산한 결과 MAPE **30.3784%**로
기준 30.2846%보다 악화했다. 종질량 잔차 `2.73e-12 kg`, 운동량 잔차
`6.98e-10 kg m/s`는 통과했지만 정량 개선이 아니므로 기본값 0으로 유지하고
후보를 폐기했다. 끝단 우회 항 자체는 보존성 회귀 테스트를 추가해 검증했다.

## 후류 질량 체류 reservoir 승격

운동량만 역류시키는 기존 closure의 한계를 확인하기 위해, cavity 내부의
conserved total/species/momentum을 보존형 reservoir에 포획한 뒤 지연 반환하는
항을 현재 반원·역류 조합에 연결했다. 포획/반환률 `0.005/0.005 s^-1`, 반환
분포 `sigma=16.7 m`, 중심에서 10 m 떨어진 ring 반환의 strict 600 s 결과는
`outputs/ti021-semidisk-backflow-reservoir-c005-r005-radius10-strict.json`에서
MAPE **29.5441%**였다. 기존 strict best 30.2846%보다 **0.7405%p 개선**했고,
종질량 잔차 `1.00e-11 kg`, 운동량 잔차 `9.31e-10 kg m/s`, 실행시간
15.89 s로 보존·실행 게이트를 통과했다.

반환률과 반경 주변 후보는 strict 또는 coarse에서 모두 더 높았다. 이 항은
센서 dose를 직접 맞춘 보정계수가 아니라 후류 체류시간을 나타내는 내부
conserved state rearrangement이므로 현재 best로 승격했다. 그래도 10% 목표까지
**19.54%p**가 남아 있으며, 다음 후보는 reservoir가 실제로 어느 높이·측면으로
재진입해야 하는지와 상단·측면 boundary packet을 결합하는 것이다.

reservoir 반환 중심을 downwind으로 5/10 m 이동한 후보와, cavity 외곽만 선택적으로
포획하는 annulus 후보도 검사했다. 반환 중심 이동은 coarse **30.00/30.05%**로
기준보다 높았고, 포획 반경 20/30 m도 **30.15/30.02%**로 개선되지 않았다.
따라서 현재 승격 조합은 대칭 포획·10 m 반환 ring이며, 다음 개선은 단순한 반경
튜닝이 아니라 높이별·면별 packet 재진입으로 제한한다.

## 추가 후보 감사: cavity 수직 전단 및 면별 reservoir 반환

남은 19.5441 percentage points를 줄이기 위해 두 가지 구조 후보를 coarse
600 s로 검사했다. 첫 후보는 장애물 후류 내부의 전단 혼합을 geometry-wide하게
강화하는 `cavity_vertical_mixing_factor`였다. 기준 0에서 0.05/0.10/0.50/1.0의
MAPE는 각각 **29.9579/29.9740/30.1060/30.2659%**로 개선되지 않았다.
따라서 cavity 상층 혼합을 관측점별 gain으로 승격하지 않는다.

두 번째 후보는 저장고를 하나의 scalar로 섞지 않고 4개 방위 sector별로 보존한
뒤 같은 sector ring으로 반환하는 `cavity_reservoir_face_return_fraction`이다.
반환 fraction 0.5/1.0의 MAPE는 각각 **30.0400/30.0657%**였고, 종질량 잔차는
`4.55e-12/8.19e-12 kg`, 운동량 잔차는 `6.98e-10 kg m/s`로 보존성은 통과했다.
그러나 정량적으로 기준보다 악화했으므로 기본값 0으로 유지하고 후보를 폐기한다.

이 감사는 단순히 ``height`` 또는 ``face`` 계수를 추가하는 것만으로는 10% 목표에
도달하지 못한다는 것을 확인한다. 다음 구현 후보는 scalar/sector reservoir가
아니라 실제 finite-volume obstacle 전면의 분리 shear-layer와 상단·측면 flux
history를 시간 지연과 함께 재진입시키는 conservative packet coupling이다.

reservoir 반환 중심을 downwind 5/10 m 이동해도 coarse **30.00/30.05%**로
악화했고, cavity 외곽 선택 포획 반경 20/30 m도 **30.15/30.02%**로 개선되지
않았다. 반환 시 수평 확산비 1.25/1.5/2.0도 **29.85/30.76/33.34%**로
현재 strict best를 넘지 못했다. 따라서 다음 구현은 scalar 계수 조정이 아니라
높이·면별 packet 재진입으로 한정한다.

## 근거리 방출 형상 교차검증

reservoir best 조합에서 Gaussian source와 finite-prism 초기 반경을 교차검증했다.
Gaussian source spread 0.75/1.0/1.25/1.5의 coarse MAPE는 **30.72/30.82/31.11/31.31%**였고,
finite-prism 반경 8/10/14/16/18/20 m는 **69.21/54.17/29.96/35.06/51.83/75.40%**였다.
따라서 현재의 14 m finite-prism 근거리 방출을 유지한다. 소스 형상·반경만으로
10%에 접근할 수 없으며, 높이별 obstacle packet 결합이 남은 핵심 물리다.

## 방출 지속시간 민감도

workbook의 `4962 kg/s × 1 s` proxy를 임의로 장기 방출로 해석하지 않기 위해,
현재 reservoir best에서 source release duration 5/10/30/60 s를 비교했다.
MAPE는 **46.45/72.35/96.39/99.81%**로 급격히 악화했다. 따라서 TI021 dose
잔차는 단순한 장기 방출 지속시간 누락으로 설명되지 않으며, 1 s proxy/초기
finite-prism 방출을 유지한다.

reservoir 반환 ring의 crosswind 비등방성 scale 0.75/1.0/1.25는 coarse
**29.99/29.96/29.96%**로 유의미한 개선을 만들지 못했다. 따라서 측면 분기는
단순 반환 폭 조정이 아니라 실제 fence side/top boundary flux의 방향·높이
history를 연결해야 한다.

reservoir 반환 반경 15 m도 strict에서 **29.6029%**로 확인되어, 10 m 반환
ring의 **29.5441%**를 유지한다. 종질량 잔차 `5.91e-12 kg`, 운동량 잔차
`6.98e-10 kg m/s`로 보존성은 통과했지만 정량 승격하지 않았다.

후류 역류와 reservoir 시간·공간 파라미터를 동시에 바꾼 joint coarse scan도
수행했다. 대표 조합 3개는 **30.76/30.38/30.10%**로 현재 best보다 높았고,
상호작용 최적점이 29.54% 아래로 내려가는 증거는 없었다. 따라서 다음 변경은
동일한 scalar parameter joint tuning이 아니라 높이별·면별 boundary packet
재진입 구현으로 제한한다.

## 반복 목표 재확인 및 상류 재순환 오프셋 감사

세션 목표를 `원인 탐색 → 물리 보완 → 동일 600 s 재계산 → 보존·경로 게이트 확인`
반복으로 고정하고, 최종 종료 조건을 dose MAPE ≤10%로 둔다. 현재 기준은
`outputs/ti021-semidisk-backflow-reservoir-c005-r005-radius10-strict.json`의
29.5441%이며, 목표까지 19.5441 percentage points가 남아 있다.

후류 저장고가 정체점보다 상류로 재진입하는 후보를 허용하기 위해 반환 오프셋의
부호를 물리적으로 열고, 동일한 strict 설정으로 −30/−20/−10/−5/0/+5 m를
재계산했다. MAPE는 각각 **29.6831/29.6361/29.6336/29.6084/29.5441/29.5874%**였다.
모든 후보의 종질량 잔차는 `1.0e-11 kg` 이하였지만 0 m 기준을 넘지 못했으므로
상류 재진입 오프셋은 채택하지 않는다. 이 결과는 남은 잔차가 단순한 저장고 중심
이동이 아니라 높이·면별 경계 패킷의 재진입 구조에 있음을 재확인한다.

## 장애물 전면 정체 후보 감사

장애물 전면 중심부의 유효 이송속도를 낮추는 `cavity_stagnation_strength`를
0.05/0.10/0.20/0.35/0.50으로 600 s coarse 재계산했다. MAPE는 각각
**32.7071/32.4462/32.0350/32.0562/31.9631%**였다. 모든 실행에서 종질량
잔차는 `6.4e-12 kg` 이하, 운동량 보존 게이트는 통과했지만 현재 promoted
strict best **29.5441%**보다 낮아지지 않았다. 따라서 정체항은 채택하지 않고,
장애물 전면 flux의 높이·면별 재진입 구현으로 다음 반복을 진행한다.

## 내부 3-D handoff 패킷 인터페이스 후보 감사

중복 방출을 피하는 내부 인터페이스를 구현하고 resolved HDF5 plane을
`x=160 m`에서 layer source-rate로 연결했다. 600 s coarse 결과는 MAPE
**86.8543%**였고, 내부면 역류량 `1.97e9 kg`, 운동량 보존 게이트 실패가
확인됐다. 따라서 이 후보는 strict 계산으로 승격하지 않고 폐기했다. 원인은
resolved plane의 좌표·시간 원점 또는 총질량 anomaly를 layer 총질량으로
해석하는 방식이 아직 맞지 않는 것이다. adapter와 보존 회귀 테스트는 유지하되,
다음 재사용 전 좌표·시간·총질량/종질량 분해를 다시 확정한다.

## 2026-10-15 후속 후보 스크리닝: 소스 전면과 외부 수송

좌표 진단에서 실제 resolved-inner 센서 변환
`solver_x=60+site_local_y`, `solver_y=320-site_local_x`가 저장된 센서 메타데이터와
일치함을 확인했다. 따라서 좌표계 표기만 수정해 MAPE가 급락할 가능성은 배제했다.

그 다음 promoted reservoir 조건을 고정하고, 동일한 600 s dose 계산을 `dt=0.5 s`로
빠르게 스크리닝했다. 시간 간격 차이로 절대값은 strict 결과보다 약간 높지만,
후보 간 판정에는 충분하다.

| 후보 | 600 s coarse dose MAPE (%) | 판정 |
|---|---:|---|
| grounded head (기준, `dt=0.5 s`) | 29.98 | 기준 |
| entrained edge | 30.17 | 폐기 |
| diffuse edge | 29.98 | 개선 없음 |
| smooth edge | 29.98 | 개선 없음 |
| wake tail | 29.98 | 개선 없음 |
| density-current 외부 continuation (`strength=1`) | 29.98 | 개선 없음 |
| density-current 외부 continuation (`strength=0.5`) | 29.98 | 개선 없음 |

소스 초기 운동량을 0.5, 1, 2배로 바꾼 진단도 농도 예측에 유의한 변화를 만들지
않았다. 이는 finite-prism 경로에서 초기 운동량이 현재 600 s 잔차를 지배하지 않는다는
뜻이다. 모든 후보는 관측 dose를 계산 중 사용하지 않았고, 질량·운동량 보존 게이트를
유지했다. 따라서 promoted strict best는 계속
`outputs/ti021-semidisk-backflow-reservoir-c005-r005-radius10-strict.json`의
**29.5441%**이며, 10% 목표까지 **19.5441 percentage points**가 남아 있다.

다음 반복은 소스 계수 조정이 아니라 resolved HDF5 handoff를 장애물 전면의
상·측면·후면 flux history로 분해하는 단계다. 현재 내부 인터페이스 후보는 signed
total-anomaly를 one-way layer packet으로 해석할 때 역류/중복 방출이 발생해
MAPE 86.85%, momentum gate 실패로 폐기된 상태다. 이 분해가 보존적으로 닫히기
전에는 해당 패킷을 생산 경로에 승격하지 않는다.

추가로 HDF5 plane의 면적·시간 적분을 직접 감사했다. marked species flux는
`3292.60 kg`이고 음의 값이 없지만, `total_mass_flux_anomaly`는 양의
`3.350e6 kg`와 음의 `-3.348e6 kg`가 거의 상쇄되어 signed 적분이 `2501.95 kg`에
불과하다. 즉 total anomaly는 대기 carrier의 왕복 배경 flux이지 dense-layer에
그대로 주입할 총질량 packet이 아니다. 이 항을 one-way total source로 사용하면
내부면이 불안정해지는 원인이 정량적으로 확인됐다. resolved packet의 momentum도
횡방향 signed 적분이 `-2.119e7 kg m/s`로 큰 왕복 성분을 포함한다.

안정성 확인을 위해 interface에서 transverse/vertical momentum을 일단 0으로
두고 species packet만 주입한 진단은 MAPE **40.7408%**, negative interface mass
`406.70 kg`, species exported/injected `3036.70/3289.65 kg`, momentum gate 통과를
기록했다. 보존·안정성은 회복됐지만 promoted strict best보다 악화했으므로 생산
후보로 채택하지 않는다. 다음 보완은 signed total을 source로 쓰지 않고, resolved
plane의 **정방향 normal flux와 height/crosswind별 packet**을 분리해 carrier와
species를 각각 보존적으로 재결합하는 것이다.

## 2026-10-15 반복 루프: 방출 높이·질량·풍향 민감도

동일한 600 s 계산 루프에서 관측값을 파라미터에 직접 피팅하지 않고, 초기
방출 높이를 바꾸는 대신 `source_mass_scale=13 m/H`로 초기 종질량을 보존해
방출 형상만 시험했다. `H=4/6/8/10/13/16/18/19/20/21/22/24/26 m`의
스크리닝 중 `H=20 m`가 coarse 기준에서 가장 낮았고, 동일한 `dt=0.25 s`로
재계산한 결과는 다음과 같다.

| 후보 | 600 s dose MAPE (%) | 종질량 잔차 (kg) | 판정 |
|---|---:|---:|---|
| 기준 `H=13 m` | 29.5441 | `1.0e-11` 이하 | promoted baseline |
| `H=18 m` | 29.6115 | `7.7e-12` | 폐기 |
| `H=19 m` | 29.5947 | `6.4e-12` | 폐기 |
| `H=20 m` | **29.4881** | `4.1e-12` | 실측 근거 부족으로 보류 |
| `H=21 m` | 29.5369 | `4.5e-12` | 폐기 |
| `H=22 m` | 29.6421 | `3.6e-12` | 폐기 |
| `H=24 m` | 29.5794 | `4.1e-12` | 폐기 |

20 m 후보는 수치상 0.056 percentage point만 개선해 10% 목표와 거리가 크고,
TI021의 실제 초기 충전 높이 13 m를 대체할 독립 측정 근거가 없으므로 승격하지
않았다. 높이 변경은 잔차의 주원인이 아니다.

전체 방출 질량 민감도도 `source_mass_scale=0.70/0.85/1.00/1.15/1.30/1.50`로
검사했다. coarse MAPE는 **40.88/33.22/29.96/30.72/34.74/42.71%**였고,
질량 1.0 부근이 이미 최소였다. 따라서 단순한 공급량 부족으로 10%까지 내려갈
수 없으며, source mass를 센서에 맞춰 조정하는 보정은 채택하지 않는다.

좌표·풍향 원인 분리를 위해 풍향각도 스캔을 수행했다. `dt=0.5 s`에서
−20/−10/0/5/6.1/10/15/20°의 MAPE는 **97.61/87.23/46.76/29.68/29.96/
34.21/43.29/52.71%**였다. 좁은 재현성 범위에서 `6.0°`를 `dt=0.25 s`로
재계산하면 **29.5128%**로 기준보다 0.031 percentage point 낮지만, 10%와는
여전히 19.5128 percentage points 차이가 난다. 풍향 표기 오류만으로는 잔차가
설명되지 않는다.

따라서 이번 반복의 결론은 `방출 높이·전체 질량·풍향`이 단독 원인이 아니라는
것이다. 20 m 후보와 6.0° 후보는 보존성은 통과하지만 외부 실측 근거가 없어
생산 경로로 승격하지 않았고, 다음 반복은 계획대로 resolved plane의 정방향
normal flux를 장애물 전면·상단·측면 packet으로 분해하는 데 집중한다.

## 추가 장애물·수직 물리 감사

측면·끝단 우회가 누락 원인인지 확인하기 위해 `wall_tangential_bypass_fraction`
과 `wall_tip_bypass_fraction`을 각각 0–0.75로 스캔했다. coarse MAPE 범위는
각각 **29.96–30.37%**와 **29.96–30.20%**였으며, 모든 후보가 질량·운동량
보존 게이트는 통과했지만 기준보다 개선되지 않았다.

수직 다층 수송(`vertical_layers=2/3/4/6/8`)도 실제 수직 관측 연산자 후보로
시험했다. 가장 나은 설정도 MAPE **47.61%**였고, 나머지는 84.60–112.44%로
악화했다. 따라서 현재 TI021 잔차를 단순히 층상 수직 혼합 누락으로 설명할 수
없으며, 다층 경로는 생산 후보에서 제외했다.

resolved HDF5 interface의 추가 진단으로 총질량을 종 flux의 일정한 물리적
질량분율 환산(`total = alpha × species`)으로 닫는 후보를 시험했다. `alpha`
1.0/1.2/1.5/2.0/2.5의 coarse MAPE는 **72.88/68.55/63.85/58.24/54.23%**로,
species-only 후보(40.74%)와 기준 hybrid(29.54%)보다 모두 나빴다. HDF5 면의
정방향 종 flux는 보존되지만, 이 면을 layer patch에 그대로 접속하는 구조 자체가
근거리 센서와 downstream 관측을 동시에 재현하지 못한다. 따라서 일정 환산계수
또한 채택하지 않으며, 향후에는 면의 시간 원점과 receptor 영역 분리를 먼저
재구성해야 한다.

## 수평 확산·중력수송 민감도 감사

수평 확산이 농도를 과도하게 퍼뜨리는지 확인하기 위해
`horizontal_diffusivity_m2_s=0/0.1/0.25/0.5/1/2/4/8`을 검사했다. MAPE는
**32.20/31.27/30.81/30.31/29.96/34.08/43.32/54.21%**로, 현재 1.0 부근이
최적이며 확산을 줄이면 오히려 오차가 커졌다. 횡방향 확산비도 0.25–3 범위에서
최선이 29.85% coarse에 머물렀고, `dt=0.25 s` 재계산에서는 기준보다 높았다.

중력수송 스케일은 `gravity_scale=0/0.25/0.5/0.75/1/1.25/1.5/2`에서
**43.30/37.67/35.12/31.03/29.96/29.82/30.63/32.02%**였다. 1.1을 엄밀한
`dt=0.25 s`로 재계산해도 **29.4985%**에 그쳤다. 중력수송은 잔차에 기여하지만,
단순 스케일 조정만으로 10%에 접근하지 못하며 현장 근거 없는 1.1–1.25 배는
승격하지 않는다.

## bounded joint interaction scan

cavity shelter/entrainment·backflow·reservoir capture/return·wake deficit·중력
수송의 상호작용을 물리적 범위 안에서 12개 조합으로 coarse 재계산했다. MAPE
범위는 **29.70–35.61%**였고, 가장 낮은 조합도 기준 29.96% coarse보다
개선되지 않았다. 종질량 잔차는 모두 `1.1e-11 kg` 이하, 운동량 잔차는
`9.4e-10 kg m/s` 이하였지만, 상호작용만으로 10%에 접근하는 근거는 없었다.

## cavity 형상·근거리 source transport 추가 감사

장애물 cavity 선택 자체를 바꾼 control도 확인했다. 물리적으로 기대되는
`downwind_semidisk`는 coarse **29.96%**, 전체 원형 cavity는 **31.93%**,
upwind 반원은 **40.62%**였다. 따라서 현재 downwind cavity 선택은 관측 방향과
일관되며 형상 변경으로 10%에 접근하지 못한다.

초기 source wind-scale은 0.25/0.5/0.75/1/1.25/1.5/2에서
**30.24/29.99/29.91/29.96/30.17/30.84/33.89%**였고, radial damping도
0–2 범위에서 **29.96–31.52%**였다. 이는 source transport의 단순 감쇠나
초기 이송속도 조정이 잔차의 주원인이 아님을 확인한다.

HDF5 수직 flux 분포를 반영하는 질량보존형 상부 가중 관측 연산자도 별도
후보로 구현해 `top_loading_exponent=0.05–2`를 coarse 시험했다. 기준 29.96%에서
MAPE가 **30.26–53.00%**로 악화했으므로 후보를 폐기하고 구현은 되돌렸다.
상부 packet이 존재한다는 사실만으로 균일 column을 상부 가중으로 바꾸는 것은
TI021 dose 잔차를 설명하지 못한다.

## 좌표 매핑 독립 감사

실제 입력 매핑 `solver_x=60+local_y`, `solver_y=320-local_x`가 맞는지 확인하기
위해 진단용으로 센서 좌표를 좌우반전·상하반전·축교환·180° 회전해 재계산했다.
기준 identity는 **29.96%**였지만 mirror-x **69.34%**, mirror-y **96.74%**,
swap **86.14%**, rot180 **96.72%**로 모두 악화했다. 따라서 현재 좌표 변환은
잔차 원인이 아니며, 변환을 변경하지 않는다.

## 세션 목표 재설정: 원인 탐색 → 물리 보완 → 600 s MAPE 반복

이번 세션의 종료 조건은 후보 파일을 많이 만드는 것이 아니라, 동일한 TI021
600 s 센서 dose 평가에서 **MAPE ≤ 10.0%**를 달성하고 재현 가능한 후보를
승격하는 것이다. 후보마다 다음 루프를 반복한다.

1. **원인 탐색**: 센서별 잔차를 저층/상층, 전면/상단/측면, 근거리/원거리와
   시간 도착·누적 dose로 분해하고, 누락된 물리를 하나만 지목한다.
2. **물리 보완**: 지목된 항만 질량·운동량 보존식과 경계조건을 통해 구현한다.
   센서별 보정계수, 관측값 직접 피팅, 임의의 source mass 증폭은 금지한다.
3. **동일 조건 재계산**: 600 s 전체를 재계산하고 `dose MAPE`, 정규화 편향,
   종질량 잔차, 운동량 잔차, 실행시간을 함께 기록한다.
4. **승격/폐기**: MAPE가 감소하고 모든 보존 게이트와 입력 근거가 통과할
   때만 승격한다. 감소하지 않으면 원인 가설을 폐기하고 다음 가설로 간다.

현재 승격 기준선은 `outputs/ti021-semidisk-backflow-reservoir-c005-r005-radius10-strict.json`의
**29.5441%**이다. 따라서 남은 정량 목표는 **19.5441 percentage points**이며,
10% 미만이 확인되기 전까지 세션 목표는 완료로 표시하지 않는다.

### 2026-10-15 수직 분해능 재검증

5 m 단일 수직층이 0.4 m와 2.4 m 센서를 같은 예측으로 만드는 결함을 확인한
뒤, 동일한 600 s 전체 계산으로 수직 보완 후보를 재검증했다. 원시 결과와
보존 게이트는 `outputs/ti021-vertical-resolution-search-2026-10-15.json`에
고정했다.

| 후보 | 600 s dose MAPE | 판정 |
|---|---:|---|
| promoted 1-layer diffusive tail | **29.5441%** | 현재 기준 |
| top-hat + ground exponent 2 | 95.1998% | 폐기 |
| dense-power exponent 2 | 95.1998% | 폐기 |
| diffusive-tail extension 5 m | 29.8941% | 폐기 |
| conservative 3 layers | 86.9729% | 폐기 |
| conservative 7 layers | 93.3317% | 폐기 |
| conservative 10 layers | 92.6086% | 폐기 |

모든 후보는 종질량 및 운동량 보존 게이트를 통과했지만 MAPE가 기준보다
낮아진 후보는 없었다. 따라서 수직층 수나 센서 관측 shape만 늘리는 것으로는
10% 목표를 달성할 수 없다. 다음 구현 후보는 관측값에 맞춘 gain이 아니라
장애물 전면·상단·측면 분리 shear layer의 시간·높이별 flux history와 후류
재진입을 보존적으로 연결하는 것이다. 해당 TI021 면별 실측 history가 없는
상태에서 합성 packet으로 10% 달성이라고 주장하지 않는다.

### 장애물 면별 데이터 게이트

로컬 `outputs`의 HDF5 handoff 파일을 전수 점검한 결과, 제공된 데이터셋은
`species_mass_flux_kg_m2_s`, `advective_species_mass_flux_kg_m2_s`,
`total_mass_flux_anomaly_kg_m2_s`, `lateral_momentum_flux_anomaly_pa`,
`vertical_momentum_flux_anomaly_pa` 등 **단일 downstream 평면의 시간·높이·횡방향
flux**만 포함한다. 장애물 전면·상단·측면의 독립 wall-face flux history 또는
face별 공통 시간축은 존재하지 않는다. 따라서 면별 재결합은 아직 실행 데이터가
없는 가설이며, 새 입력이 확보되기 전에는 임의 packet 분해로 10% 목표를
주장하지 않는다.

## 외부 경계 재유입 가설의 600 s 검증

`remaining_physics`에 남아 있던 patch 외부 재순환을 검증하기 위해, 실제
downstream 경계에서 방출된 4-성분 보존 packet의 일부를 지연시켜 cavity로
되돌리는 선택적 경로를 추가했다. 센서별 보정이나 질량 생성은 사용하지 않았다.
그러나 기준선보다 모두 악화했고 운동량 예산도 통과하지 못했다.

| 재유입 fraction | delay (s) | sigma (m) | MAPE (%) | 종질량 잔차 (kg) | 운동량 잔차 (kg m/s) | 판정 |
|---:|---:|---:|---:|---:|---:|---|
| 0.02 | 20 | 15 | 36.3740 | `6.4e-12` | 99.66 | 폐기 |
| 0.05 | 20 | 15 | 36.0420 | `5.0e-12` | 249.43 | 폐기 |
| 0.10 | 20 | 15 | 35.8313 | `5.0e-12` | 499.79 | 폐기 |
| 0.05 | 40 | 20 | 36.2221 | `4.0e-12` | 249.12 | 폐기 |

이 결과는 단순한 downstream 경계 재유입이 근거리 저층 under-prediction의
원인이 아님을 보여준다. 해당 경로는 기본 설정에서 비활성으로 유지하며,
면별 CFD flux가 확보되기 전에는 승격하지 않는다.

## 근거리 3-D/층 모델 인계 반경 감사

근거리 방출·중첩이 부족한지 확인하기 위해 2-D dense-layer로 넘기는
`handoff_radius_m`을 15 m 기준에서 20/30/40 m로 늘려 동일한 600 s를 다시
계산했다. 인계 반경을 늘리면 근거리 3-D 대체 구간이 넓어지지만 MAPE는 오히려
악화했다.

| 인계 반경 (m) | MAPE (%) | 종질량 잔차 (kg) | 판정 |
|---:|---:|---:|---|
| 15 (기준) | 29.5441 | `1.0e-11` 이하 | 기준 |
| 20 | 36.7609 | `1.1e-11` | 폐기 |
| 30 | 37.9660 | `6.0e-12` | 폐기 |
| 40 | 39.5359 | `7.0e-12` | 폐기 |

따라서 단순히 3-D 영역을 넓히는 것이 저층 센서 under-prediction을 해결하지
못한다. 다음 후보는 인계 반경 조정이 아니라, 실제 장애물 면별/높이별 유동과
센서 높이에서의 수직 분리를 제공하는 입력이어야 한다.

## 방출 시간 이력의 인계 처리 감사

초기 방출이 실제로 수 초에 걸쳐 형성되는 경우를 확인하기 위해
`source_release_duration_s`를 켜고, 방출 완료 시점 이후에 dense-layer로
인계되도록 반경을 함께 늘렸다. 단순 staged release는 인계 시점 이전의 일부
질량만 layer에 전달하고 이후 방출을 버리는 구조이므로, 물리적 source history
연결로 승격하지 않았다.

| 방출 기간 (s) | 인계 반경 (m) | 인계 시각 (s) | MAPE (%) | 판정 |
|---:|---:|---:|---:|---|
| 5 | 30 | 9.55 | 47.6816 | 폐기 |
| 10 | 40 | 15.18 | 79.6808 | 폐기 |

두 후보의 종질량 잔차는 `2e-12 kg` 이하였지만 기준선 29.5441%보다 크게
악화했다. 이 결과는 방출 기간 자체가 불필요하다는 뜻이 아니라, 방출 이력을
초기 cloud 상태에서 잘라 전달하는 현재 인터페이스가 부적합하다는 뜻이다.
다음 source-history 후보는 방출 packet을 dense-layer 시간축에 직접 연결하고
동일 질량·운동량을 보존하는 별도 handoff가 필요하다.

### 시간축 source-injection prototype

위 결함을 보완하기 위해 인계 시점 이후의 잔여 source inventory를 layer
시간축에 직접 주입하는 보존 callback을 구현했다. 5 s release, 10 m handoff
후보의 결과는 MAPE **42.9466%**였고, 종질량 잔차는 `2.4e-12 kg`였지만
운동량 ledger 수정 후 잔차는 `6.98e-10 kg m/s`로 게이트를 통과했다. 그러나
MAPE가 기준선보다 악화했으므로 이 callback은 구현 prototype으로만 남기고
승격하지 않는다. 현재 결과는 source history를
연결하는 것만으로 충분하지 않으며, 주입 packet의 위치·운동량과 layer 경계
flux를 하나의 동일한 보존 ledger로 결합해야 함을 보여준다.
