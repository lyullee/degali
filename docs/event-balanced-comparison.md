# Event-balanced fixed-receptor comparison

센서 행은 하나의 방출 사건 안에서 반복 측정되는 관측이다. 따라서
`N=210` 센서 행을 `N=210`개의 독립 시험으로 간주하면 센서가 많이 배치된
한 사건이 전체 점수를 지배할 수 있다. `degali.event-balanced-metrics.v1`은
기존 pooled 점수를 없애지 않고, 다음 두 점수를 함께 고정한다.

- `pooled`: 모든 `exact` 센서 행을 한 번에 계산한 추적용 MAE/RMSE/bias/max error.
- `macro_equal_event_weight`: 방출 사건마다 먼저 MAE/RMSE/bias를 계산한 뒤 사건별
  동일 가중치로 평균한 점수. 사건 수가 성능의 통계적 단위다.

`lower_bound` 관측은 대칭 오차에 넣지 않는다. 대신 모델 예측이 관측 하한
이상인지 `satisfied_count`와 `satisfaction_fraction`으로 보고한다. 이 구분은
검출한계, 포화, 또는 공개 데이터의 일방 정보가 좋은 점수로 변환되는 것을 막는다.

## Python contract

```python
from degali.addons.event_balanced_metrics import (
    EventMetricObservation,
    score_event_balanced_metrics,
)

score = score_event_balanced_metrics(
    [
        EventMetricObservation("FFI-Test4", "OC_01", 0.02, 0.018),
        EventMetricObservation("FFI-Test6", "OC_01", 0.04, 0.039),
        EventMetricObservation(
            "FFI-Test6", "OC_02", 0.04, 0.035, observation_kind="lower_bound"
        ),
    ],
    threshold_mole_fraction=0.04,
)
report = score.as_dict()
```

입력은 `(event_id, sensor_id)`가 유일해야 하고, 농도는 mole fraction `[0, 1]`
범위로 전달한다. `as_dict()` 결과에는 스키마 버전, 사건별 요약, pooled/macro
점수, 2/4 vol% 판정에 해당하는 threshold count, 하한 제약이 함께 남는다.

## FFI strategy integration

`tools/analyse_slabx_competitive_strategy.py`는 현재 Test4 진단에서도 이 계약을
생성한다. 단일 사건에서는 pooled와 macro가 같지만, 7개 FFI 방출을 결합할 때는
각 사건을 하나의 독립 단위로 유지해야 한다. 비교 입력은 먼저
`FieldEvidenceManifest`의 event/source/weather/receptor/sensor 시간축과 연결하고,
source·weather·센서 형상·시간 연산자가 비교 가능한 경우에만 경쟁 점수를 발행한다.
그 외에는 수치와 원인을 남기되 `conditional` 또는 `withheld`로 유지한다.

이 점수는 모델을 자동 선택하거나 승인하는 규칙이 아니다. 동일 사건군, 동일
관측 연산자, 동일 threshold와 결측·포화 정책을 잠근 뒤, 정확도·FN/FP·경보시간·
보존량 잔차를 함께 검토하기 위한 비교 계층이다.
