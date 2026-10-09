# 현장 이벤트 증거 패키지 입력 템플릿

실제 이벤트를 등록할 때는 아래 디렉터리를 하나의 읽기 전용 패키지로
만든다. 파일을 수정하거나 교체하면 기존 SHA-256 manifest는 폐기하고
`field-audit`부터 다시 실행한다.

```text
<event-id>/
├── source_boundary.csv          # 누출 위치·형상·상태·대기 질량유량/시간축
├── weather.csv                  # 풍속·풍향·측정높이·안정도·시간축
├── obstacle_geometry.json       # 장애물 좌표·높이·폭·길이·좌표계/수직기준
├── receptor_observations.csv    # sensor_id, x, y, z, timestamp, H2, unit
├── sensor_registry.csv          # sensor_id, calibration, response, validity
├── common_clock.json            # 시간원·동기화 방법·offset·적용 구간
└── package-notes.md             # 이벤트/운영자/수집 절차 메모
```

## 최소 필드

| 채널 | 반드시 확인할 항목 |
|---|---|
| `source_boundary` | `event_id`, `source_boundary_id`, timestamp, 위치/형상, phase, `mass_flow_kg_s` 또는 승인된 pressure-derived provenance |
| `weather` | `event_id`, `weather_id`, timestamp, `wind_speed_m_s`, wind direction, reference height, stability, units |
| `obstacle_geometry` | `event_id`, `obstacle_geometry_id`, coordinate reference, vertical datum, vertices 또는 dimensions/height |
| `receptor_observations` | `event_id`, `sensor_id`, x/y/z, timestamp, H2 concentration, unit, quality flag |
| `sensor_registry` | `sensor_set_id`, sensor IDs, calibration ID/date, gain/bias, response time, validity interval |
| `common_clock` | `common_clock_id`, source clock, synchronization method, offset/uncertainty, covered channels |

`temporal_operator_id`는 평균·최대·dose 등 비교에 사용한 시간 연산자를
고정하고, `operator_id`는 파일을 수집·대조한 책임 주체를 기록한다. 여러
센서가 있으면 manifest에는 하나의 `sensor_set_id`를 기록하되, 실제 센서
식별자는 관측 파일의 각 행에 남긴다.

## 생성 순서

```console
degali field-audit <event-dir> --output audit-execution.json
degali field-audit-verify audit-execution.json --require-complete
degali field-evidence-manifest-create audit-execution.json \
  --selected-path source_boundary=<event-dir>/source_boundary.csv \
  --selected-path weather=<event-dir>/weather.csv \
  --selected-path obstacle_geometry=<event-dir>/obstacle_geometry.json \
  --selected-path receptor_observations=<event-dir>/receptor_observations.csv \
  --selected-path common_clock=<event-dir>/common_clock.json \
  --manifest-id manifest-<event-id> \
  --event-id <event-id> \
  --dataset-id receptors-<event-id> \
  --observed-row-count <positive-row-count> \
  --source-boundary-id <source-id> \
  --weather-id <weather-id> \
  --obstacle-geometry-id <obstacle-id> \
  --receptor-geometry-id <receptor-id> \
  --temporal-operator-id <operator-id> \
  --common-clock-id <clock-id> \
  --sensor-set-id <sensor-set-id> \
  --operator-id <accountable-operator-id> \
  --sensor-registry-path <event-dir>/sensor_registry.csv \
  --sensor-calibration-status certified \
  --output manifest-execution.json
degali field-evidence-manifest-verify manifest-execution.json
```

현재 `field-audit`가 5개 채널 중 하나라도 찾지 못하면 manifest를 만들지
않고 `partial` 또는 `withheld`로 남긴다. 모든 채널과 해시는 맞지만 센서
registry 또는 책임 운영자 ID가 없으면 manifest는 `conditional`로 생성될 수
있다. registry에 nominal specification만 있거나 사건별 교정성적서가 없으면
`--sensor-calibration-status specification_only` 또는 `missing`으로 기록하며
역시 `conditional`이다. `certified`를 명시한 경우에도
`promotion_allowed=false`는 유지되며, 해당 메타데이터가 보완되기 전에는
검증용 qualification 입력으로 승격하지 않는다.
