# DEGALI 상세 사용 가이드

이 문서는 공개 배포판 DEGALI를 처음 설치하는 사용자부터 기존 DEGADIS
입력 덱을 사용하는 연구자까지를 대상으로 한다. DEGALI는 **Dense Gas
Dispersion for Liquid Hydrogen**의 약자이며, DEGADIS 2.1의 검증된 Python
재구현 위에 액화수소(LH₂) 물리를 확장한 연구용 확산 모델이다.

> **중요:** DEGALI는 연구, 시나리오 비교 및 민감도 분석용 소프트웨어다.
> 단독으로 인허가 거리, 설비 배치 또는 안전 필수 설계를 확정하는 도구가
> 아니다. 계산 결과와 함께 `warnings`를 반드시 검토하고, 안전에 중요한
> 판단은 독립 실험 또는 승인된 해석 절차와 교차 확인해야 한다.

## 1. 어떤 실행 방법을 선택할까

| 목적 | 권장 인터페이스 | 입력 |
|---|---|---|
| LH₂ 누출을 빠르게 평가 | `degali.lh2.assess()` | 물리량을 Python 인자로 입력 |
| 기존 DEGADIS 정상상태 계산 | `run_steady()` 또는 `degali steady` | `.INP` 덱 |
| 기존 DEGADIS 비정상 계산 | `run_transient()` 또는 `degali transient` | `.INP` 덱 |
| 고정 수용점 농도·노출량 | `TransientOutput.dose()` 또는 `degali dose` | `.INP` 덱과 거리 |
| 기존 JETPLU 제트 계산 | `run_jet()` 또는 `degali jet` | `.INO` 덱 |
| 제트 착지 후 지표 확산 | `run_jet_to_ground()` | `.INO`와 `.IN` 덱 |

처음 사용하는 경우에는 입력 덱이 필요 없는 `assess()`부터 시작하는 것이
가장 간단하다. 기존 DEGADIS 호환 경로는 사용자가 적법하게 확보한 입력
덱을 가지고 있을 때 사용한다. 원본 FORTRAN과 제3자 실험 데이터는 배포본에
포함되지 않는다.

## 2. 설치

### 2.1 권장 환경

- Python 3.10 이상
- 64비트 Python 권장
- 새 가상환경 사용 권장

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install "degali[coolprop]"
```

Linux 또는 macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "degali[coolprop]"
```

설치 확인:

```bash
python -c "import degali; print(degali.__version__)"
degali --help
```

LH₂ 계산에는 CoolProp이 필요하다. `pip install degali`만 설치하면 전통
DEGADIS 호환 계산은 가능하지만 `assess()` 호출 시 CoolProp 가져오기 오류가
발생할 수 있다.

개발·테스트 환경이 필요한 경우:

```bash
git clone https://github.com/lyullee/degali.git
cd degali
python -m pip install -e ".[coolprop,test]"
python -m pytest -m "not slow" -q
```

## 3. 가장 간단한 LH₂ 계산

다음 예제는 6 bar(a)에 저장된 포화 액체수소가 높이 0.5 m의 10 mm
오리피스에서 0.1 kg/s로 수평 누출되고, 누출 높이의 풍속이 2.0 m/s인
경우를 계산한다.

```python
from degali.lh2 import assess

result = assess(
    rate=0.10,                 # 수소 질량유량, kg/s
    wind=2.0,                  # 누출 높이 풍속, m/s
    height=0.50,               # 누출구 지상 높이, m
    orifice=0.010,             # 오리피스 직경, m
    storage_pressure=6.0,      # 저장 절대압력, bar(a)
    ambient_temperature=288.15,# 대기온도, K
    relative_humidity=65.0,    # 상대습도, %
    ambient_pressure=101325.0, # 대기압, Pa
    max_distance=30.0,         # 계산 최대거리, m
    at_distance=10.0,          # 최저 가연 높이를 확인할 거리, m
)

print(result.report())
```

`storage_pressure`는 **게이지압이 아니라 절대압 bar(a)** 이다. 예를 들어
5 bar(g)는 대기압이 약 1 bar일 때 약 6 bar(a)로 입력한다.

## 4. `assess()` 입력값

| 인자 | 단위 | 의미 | 기본값 |
|---|---:|---|---:|
| `rate` | kg/s | 수소 질량 누출률 | 필수 |
| `wind` | m/s | 누출 높이에서의 풍속 | 필수 |
| `height` | m | 누출구 높이; 풀은 보통 0 | `0.0` |
| `orifice` | m | 가압 누출 오리피스 직경 | 둘 중 하나 필수 |
| `pool_diameter` | m | 풀 또는 증발원 직경 | 둘 중 하나 필수 |
| `storage_pressure` | bar(a) | 가압 저장 절대압력 | `1.013` |
| `ambient_temperature` | K | 대기온도 | `288.15` |
| `relative_humidity` | % | 상대습도 | `65.0` |
| `ambient_pressure` | Pa | 대기 절대압력 | `101325.0` |
| `max_distance` | m | 수치 적분 최대거리 | `100.0` |
| `at_distance` | m | 최저 가연 높이 평가 위치 | 선택 |

`orifice`와 `pool_diameter`는 정확히 하나만 지정해야 한다.

### 4.1 가압 제트

```python
from degali.lh2 import assess

jet = assess(
    rate=0.20,
    wind=1.5,
    height=0.50,
    orifice=0.0254,
    storage_pressure=3.0,
    max_distance=40.0,
)
print(jet.report())
```

가압 제트 경로는 플래싱 등가원과 제트 궤적을 계산한다. `height`는 0보다
커야 하며, 모델은 현재 수평·평균풍 정렬 누출을 중심으로 검증됐다.

### 4.2 풀 또는 저운동량 증발원

```python
from degali.lh2 import assess

pool = assess(
    rate=9.5,
    wind=3.0,
    pool_diameter=9.1,
    ambient_temperature=288.15,
    relative_humidity=65.0,
    max_distance=35.0,
)
print(pool.report())
```

현재 풀 경로는 주어진 수소 증발률과 직경에서 시작한다. 액체 유출량으로부터
시간에 따른 풀 크기와 증발률을 산정하는 독립 풀 확장 모델은 포함하지 않는다.
따라서 `rate`와 `pool_diameter`는 외부 누출원 모델 또는 보수적 공학 가정으로
먼저 산정해야 한다.

증발 이력이 있으면 `assess_pool_history()`가 준정상 스냅샷을 유지한다.
수신점과 사용자가 선언한 응답시간을 알고 있는 경우
`response_time_s=...`를 추가하면 인과적 1차 관측 응답 커널을 적용한
`receptor_transient_mole_fraction`도 계산한다. 이는 저차 연구용 근사이며
비정상 CFD나 검증된 puff 모델이 아니고, 응답시간을 데이터에서 자동 피팅하지
않는다.

## 5. 결과 해석

`Assessment`의 주요 속성은 다음과 같다.

| 속성 | 의미 |
|---|---|
| `distance_to_lfl` | 중심선 농도가 수소 LFL 4 mol%로 감소하는 거리(m) |
| `distance_to_stoichiometric` | 중심선 농도가 화학양론 농도로 감소하는 거리(m) |
| `lowest_flammable_height` | 지정 위치 또는 가연영역 말단에서 가장 낮은 가연성 가스 높이(m) |
| `regime` | `grounded`, `low`, `aloft` 중 하나 |
| `neutral_buoyancy` | 혼합물이 공기보다 가벼워지는 수소 몰분율 |
| `trajectory` | `[x, z, 중심선 몰분율]` 형식의 NumPy 배열 |
| `warnings` | 검증범위 이탈 또는 물리적 적용성 경고 목록 |
| `notes` | 사용한 원천·모델 경로·플래싱 후 수소 질량분율 |

`result.screening_scope`는 `qualified`, `conditional`, `out_of_scope` 중
하나다. 외삽값을 자동으로 차단해야 하는 운영 경로에서는
`strict_scope=True`를 전달한다. 검증범위 경고가 발생하면
`ApplicabilityError`가 발생한다. 기본 경고 모드는 연구·민감도 분석용이며
인허가나 설계 승인 판단을 의미하지 않는다.

### 5.1.1 다른 모델과의 비교

HyRAM·PHAST·EFFECTS 또는 CFD의 독립 계산값이 있을 때는
`degali.validation.compare_models`를 사용한다. `ComparisonCase`가 source,
풍속·풍향, 방출 높이, 센서 연산자, 평균시간, 상평형 폐쇄, 형상을 함께 기록한다.
조건이 다르거나 관측값에 맞춰 피팅된 결과는 진단 통계는 계산하지만 우열 순위에서는
자동 제외된다. 서로 다른 실행조건의 수치로 “일반적으로 더 정확하다”는 주장을 만드는
것을 막는 비교 프로토콜이다.

단일 중심선 관측점의 원인 분해에는 `assess_observation_envelope()`를 사용한다.
사용자가 선언한 source 질량률·풍속 조합을 모두 계산하고, 관측값의 지정 배수 안에
들어오는 모든 조합을 반환한다. 여러 조합이 남으면 결과를 식별불능으로 표시하며,
보정계수를 피팅하거나 하나를 임의로 선택하지 않는다.
센서 높이값이나 arc 최대값은 별도의 검증된 관측 연산자가 필요하며 중심선 값으로
자동 대체하지 않는다.
지정한 마스트 좌표 `(x, y, z)`에서 실제 연산자를 적용하려면
`project_lh2_jet_to_sensors(points_m=[...])`를 사용한다. Gaussian 수직·횡방향
프로파일로 농도와 온도를 직접 계산하며, 계산 궤적 밖 점은 `nan`으로 남긴다.

`nan`은 보통 계산 실패가 아니라, 지정한 적분 범위 안에서 해당 농도 교차가
발견되지 않았다는 뜻이다. 이때 무조건 `max_distance`만 늘리지 말고 궤적,
경고 및 입력값 단위를 먼저 확인한다.

구조화된 결과 사용 예:

```python
import math

if math.isnan(result.distance_to_lfl):
    print("계산 범위에서 LFL 교차가 확인되지 않았습니다.")
else:
    print(f"LFL 중심선 거리: {result.distance_to_lfl:.2f} m")

print("플룸 상태:", result.regime)
print("모델 경로:", result.notes["model path"])

for warning in result.warnings:
    print("주의:", warning)
```

### 5.1 궤적을 CSV로 저장

```python
import numpy as np

np.savetxt(
    "lh2_trajectory.csv",
    result.trajectory,
    delimiter=",",
    header="distance_m,height_m,centreline_mole_fraction",
    comments="",
)
```

### 5.2 궤적과 농도 그리기

Matplotlib은 DEGALI의 필수 의존성이 아니므로 별도로 설치한다.

```bash
python -m pip install matplotlib
```

```python
import matplotlib.pyplot as plt

x = result.trajectory[:, 0]
z = result.trajectory[:, 1]
c = result.trajectory[:, 2]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(x, z)
ax1.set_ylabel("Plume centre height (m)")
ax1.grid(True)

ax2.semilogy(x, c)
ax2.axhline(0.04, color="red", linestyle="--", label="H2 LFL")
ax2.set_xlabel("Downwind distance (m)")
ax2.set_ylabel("Centreline mole fraction")
ax2.grid(True)
ax2.legend()

fig.tight_layout()
plt.show()
```

## 6. 여러 조건 비교

동일한 모델 설정에서 풍속 민감도를 비교하는 예다. 경고가 있는 계산값은
표에 그대로 기록해 검증범위 밖의 외삽을 숨기지 않는다.

```python
import csv
from degali.lh2 import assess

rows = []
for wind in (0.6, 1.0, 2.0, 3.0, 4.2):
    r = assess(
        rate=0.10,
        wind=wind,
        height=0.50,
        orifice=0.010,
        storage_pressure=6.0,
        max_distance=30.0,
    )
    rows.append({
        "wind_m_s": wind,
        "lfl_distance_m": r.distance_to_lfl,
        "regime": r.regime,
        "warnings": " | ".join(r.warnings),
    })

with open("wind_sensitivity.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
```

시나리오 비교에서는 한 번에 하나의 입력만 바꾸고, 모든 입력·버전·경고를
결과 파일에 함께 보존하는 것이 좋다.

원천 유량과 풍속의 신뢰 가능한 범위를 직접 넣어 입력 민감도 봉투를
계산할 수도 있다. DEGALI가 임의의 오차율이나 보정계수를 만들지는 않는다.

```python
from degali.lh2 import assess_envelope

envelope = assess_envelope(
    rates=[0.25, 0.285, 0.32],
    winds=[2.3, 2.5, 2.7],
    height=0.50, orifice=0.0254, storage_pressure=6.0,
    max_distance=6.0,
)
print(envelope.distance_to_lfl_range)
print(envelope.screening_scope, envelope.warnings)
```

이는 FFI Test 6에 맞춘 피팅이나 신뢰구간이 아니라, 사용자가 선언한
입력 범위의 결정론적 민감도 계산이다.

## 7. 현재 검증범위와 경고

공개 `assess()` 경로가 코드에서 확인하는 범위는 다음과 같다.

| 경로 | 질량유량 | 풍속 | 거리 | 기하 |
|---|---:|---:|---:|---|
| 수평 제트 | 0.084–0.285 kg/s | 0.6–4.2 m/s | 0.79–6.0 m | 오리피스 0.006–0.0254 m, 높이 0.5 또는 1.5 m |
| 풀 | 9.2–10.3 kg/s | 1.55–6.30 m/s | 0–33.8 m | 9.1 m 단일 풀 규모 |

이는 사용 가능 영역을 강제로 자르는 수치 한계가 아니라, 현재 직접 비교된
실험영역이다. 범위 밖에서도 계산은 수행될 수 있지만 `warnings`가 발생하며
불확실성이 커진다. 특히 출구속도/풍속 비가 10보다 낮은 제트는 바람에 강하게
휘어질 수 있어 정상·평균풍 정렬 제트 농도를 방어하기 어렵다.

다음 조건은 현재 단독 설계 판단 범위가 아니다.

- 장애물, 건물 후류 및 복잡 지형(아래 기하 화면은 중심 궤적 접촉만
  검출하며, 후류 모델은 아니다)
- 실내 또는 부분 밀폐 공간
- 기본 `assess()` 경로에서의 임의 방향·강한 횡풍·하향 충돌 제트. 수평
  방출의 yaw와 풍향을 계산하는 별도
  `run_lh2_yawed_crosswind_research()` 경로는 추가되었지만, 연구 전용이며
  설계 정량값으로 검증된 경로가 아니다.
- 독립 풀 확장·기화율 계산
- 임의 설비 조건 전체에 대한 인증된 사고결과 해석
- 점화, 화염, 복사열, 폭발 과압 계산

### 7.1 연구용 시설 기하 화면

별도로 계산된 3차원 중심 궤적에 대해, 선택형 기하 화면은 **자유 플룸**
적분 결과가 더 이상 적용되지 않는 지점을 찾는다. 좌표는 전역 `(x, y, z)`
미터다. 풍향 좌표계의 각도는 전역 `+x`축에서 반시계 방향으로 잰 라디안
**향하는(to)** 방향이며, 기상학적 불어오는(from) 방향이 아니다.

```python
from degali.addons import (
    AxisAlignedCuboid, TransverseWall, WindFrame, screen_trajectory,
)

path = [(0.0, 0.0, 1.0), (4.0, 0.0, 1.2), (8.0, 0.2, 1.8)]
obstacles = [
    AxisAlignedCuboid(3.0, 4.5, -0.8, 0.8, 0.0, 3.0, label="compressor"),
    TransverseWall(WindFrame(direction_rad=0.0), downwind_m=6.0,
                   base_height_m=0.0, height_m=2.5, half_width_m=None),
]
screen = screen_trajectory(path, obstacles)
if not screen.free_plume_prediction_applicable:
    first = screen.encounters[0]
    print(first.obstacle_label, first.entry_point_m)
    # 이 지점 이후 자유 플룸 농도장은 보고하지 않는다.
    assert screen.free_plume_applicable_through(first.entry_arc_length_m - 1e-6)
```

`half_width_m=None`은 횡방향 무한 방벽의 이상화다. 이때 반환되는
`required_overflight_height_m`는 연속적인 중심 경로가 넘기 위해 필요한 최소
기하학 높이일 뿐, 플룸 상승 예측값이 아니다. 유한 폭 방벽은 위로 넘거나
양옆으로 우회할 수 있으므로 화면은 어느 경로도 선택하지 않는다. 선택형
`YawedCrosswind` 계산 뒤에는 `YawedTrajectory.obstacle_screen(obstacles)`로
동일한 검사를 할 수 있다.

이 화면은 후류 계수, 우회 유량 비율, 벽 열전달, 반사 농도장 또는
장애물 보정 농도를 절대 적용하지 않는다. 중심선이 비어 있어도 유한한
가우시안 플룸 외곽이 구조물을 피한다는 뜻은 아니다. 건물·방벽 하류 농도는
장애물 해상 및 별도 검증을 갖춘 해석으로 평가해야 한다.

## 8. 기존 DEGADIS 입력 덱 사용

### 8.1 정상상태

명령행:

```bash
degali steady CASE.INP
degali steady CASE.INP --backend coolprop
degali steady CASE.INP --er1 CUSTOM.ER1 --er2 CUSTOM.ER2
```

Python:

```python
from degali import run_steady

profile, source = run_steady("CASE.INP", backend="legacy")
print("5 mol% 도달거리:", profile.distance_to(0.05))
print("LLC 이상 질량:", profile.mass_above_lfl)

distance = profile.column("dist")
mole_fraction = profile.column("yc")
temperature = profile.column("temp")
```

정상상태 `Profile` 열은 `dist`, `yc`, `cc`, `rho`, `gamma`, `temp`, `b`,
`sz`, `sy` 순서다. 가능한 경우 `column("이름")`을 사용하고 배열 열 번호를
직접 하드코딩하지 않는다.

### 8.2 비정상상태와 순간별 구름

명령행:

```bash
degali transient CASE.INP --snapshot 60 --snapshot 120 --snapshot 180
```

Python:

```python
import numpy as np
from degali import run_transient

run = run_transient(
    "CASE.INP",
    times=np.array([60.0, 120.0, 180.0]),
    backend="legacy",
)

for snapshot in run.snapshots:
    print(snapshot.time, snapshot.distance_to(0.05), snapshot.mass_above_llc)
```

순간장 열은 `dist`, `dist0`, `yc`, `cc`, `ccstr`, `rho`, `gamma`, `temp`,
`sz`, `sy`, `b`다.

### 8.3 고정 수용점 농도와 노출량

명령행:

```bash
degali dose CASE.INP --at 50 --at 100 --at 200
```

Python에서는 횡방향·높이 방향 오프셋도 지정할 수 있다.

```python
from degali import Receptor, run_transient

run = run_transient("CASE.INP")
receptors = [
    Receptor(x=100.0),
    Receptor(x=200.0, offsets=[(0.0, 1.5), (5.0, 1.5)]),
]

for history in run.dose(receptors):
    peak_fraction, peak_time = history.peak
    print(history.receptor.x, peak_fraction, peak_time)
    print("1차 농도 노출량:", history.dose(exponent=1.0))
```

`DoseHistory.dose()`의 단위는 `(몰분율)^n·s`이다. 독성 또는 위해성 기준의
ppm·min 등과 혼용하지 말고 필요한 단위 변환과 물질별 지수를 별도로 적용한다.

### 8.4 제트와 지표 확산 연결

```bash
degali jet JET.INO
degali jet JET.INO --bridge GROUND.IN
```

```python
from degali import run_jet, run_jet_to_ground

jet, deck = run_jet("JET.INO")
print(jet.touchdown, jet.distance, jet.halfwidth)

profile, jet, source = run_jet_to_ground("JET.INO", "GROUND.IN")
```

플룸이 관심 농도에서 지면에 도달하지 않으면 지표 확산으로 연결되지 않는다.
이는 반드시 오류가 아니라 해당 모델 경로에서의 물리적 결과일 수 있다.

### 8.5 입력 덱 주의사항

`.INP`, `.INO`, `.IN`, `.ER1`, `.ER2`는 자유 형식처럼 보이지만 **값의 순서에
의미가 있는 위치 기반 파일**이다. 단위를 포함한 핵심 원칙은 다음과 같다.

- 풍속 m/s, 길이 m, 시간 s, 온도 K
- 대기압은 전통 덱에서 atm, `assess()`에서는 Pa
- 농도 한계는 몰분율(예: 4%는 `0.04`)
- `CHECK4`가 참이면 정상상태, 거짓이면 비정상상태
- `.INO`는 제트 조건과 자체 물성 표를 포함
- 사용자 정의 `.ER1`·`.ER2`가 없으면 EPA 예제 기본 수치계수가 사용됨

덱 전체 필드 순서는 [`src/degali/io/inp.py`](src/degali/io/inp.py)와
[`src/degali/io/jetdeck.py`](src/degali/io/jetdeck.py)의 모듈 설명에 기록돼
있다. 기존 덱을 수정할 때는 원본을 보존하고 한 항목씩 변경한다.

## 9. `legacy`와 `coolprop` 선택

- `legacy`: DEGADIS 2.1 재현과 기존 결과 비교에 사용한다.
- `coolprop`: 지원되는 물질의 실제 유체 물성과 현대적 열역학 역산을 사용한다.
- LH₂ `assess()`는 내부적으로 CoolProp을 사용한다.

백엔드 변경은 단순한 수치 정밀도 옵션이 아니라 물성 모델 변경이다. 두 결과가
다를 때 하나를 임의로 정답으로 간주하지 말고, 입력 온도·압력·상태와 적용
목적을 확인한다.

## 10. 오류와 문제 해결

### `ModuleNotFoundError: No module named 'CoolProp'`

```bash
python -m pip install "degali[coolprop]"
```

### `give exactly one of an orifice diameter or a pool diameter`

`orifice`와 `pool_diameter` 중 하나만 지정한다.

### `a pressurised jet requires a positive release elevation`

제트 계산에서 `height`를 0보다 크게 입력한다.

### 결과가 `nan`

해당 농도 교차가 적분 구간에 없을 수 있다. 단위, `trajectory`의 마지막 농도,
`max_distance` 및 `warnings`를 순서대로 확인한다.

### 실행이 오래 걸림

LH₂ 실제 물성, 플래싱 및 제트 적분은 단순 가우시안 식보다 계산량이 크다.
민감도 분석에서는 먼저 적은 조건과 짧은 `max_distance`로 시험한 뒤 범위를
늘린다. 같은 입력의 결과를 반복 계산하지 말고 파일로 저장한다.

### 입력 덱 파싱 실패

덱은 위치 기반이다. 누락된 한 값이 이후 모든 필드의 의미를 바꿀 수 있다.
헤더 4줄, 행 개수, 물성표 행 수 및 논리값 위치를 확인한다.

## 11. 결과 기록 권장사항

재현 가능한 보고서에는 최소한 다음을 기록한다.

1. DEGALI 버전과 DOI
2. 모든 입력값과 단위
3. 사용한 경로(제트/풀/기존 덱)와 물성 백엔드
4. `warnings` 전체
5. 출력 CSV와 실행 코드
6. 결과를 비교한 독립 자료 또는 해석 도구
7. 검증범위 밖 외삽 여부와 적용한 안전여유

버전 확인:

```python
import degali
print(degali.__version__)
```

DEGALI 0.2.0 버전 DOI는
[`10.5281/zenodo.23105451`](https://doi.org/10.5281/zenodo.23105451)이며,
모든 버전을 묶는 개념 DOI는
[`10.5281/zenodo.22646258`](https://doi.org/10.5281/zenodo.22646258)이다.

## 12. 추가 기술문서

- [물리이론·검증·결과 통합 기술문서](docs/technical-reference.md)
- [현재 모델 상태와 한계](docs/status.md)
- [검증 과정](docs/validation.md)
- [정량 현장시험 검증](docs/field-validation.md)
- [검증 주장 등급](docs/claim-grading.md)
- [외부 자료 출처와 재현 절차](docs/DATA_AND_REPRODUCTION.md)
- [공개 배포 범위](docs/publication-scope.md)
- [안전 및 보안 보고](SECURITY.md)

오류 보고 시 운영체제, Python 버전, DEGALI 버전, 최소 재현 코드, 전체 오류
메시지를 포함하되 제3자 원시 실험 데이터는 GitHub 이슈에 첨부하지 않는다.
# 유한 방출 제트–플룸–퍼프 실행

유한 시간 동안 방출되는 수평 LH2 제트는 연구 API 한 번으로 근접장,
횡풍 플룸, 정확한 방출 종료 전환, 3차원 퍼프 및 고정 수용점 이력까지
계산할 수 있습니다. 모든 입력은 SI 단위입니다.

```python
from degali.addons.axisymmetric_jet import AxisymmetricJetSource
from degali.lh2 import run_lh2_finite_release_research

source = AxisymmetricJetSource(
    diameter=0.001,
    velocity=500.0,
    density=0.5,
    temperature=45.0,
    theta=0.0,
    y=0.5,
)

result = run_lh2_finite_release_research(
    source,
    source_duration_s=5.0,
    puff_duration_s=60.0,
    wind=2.5,
    wind_angle=0.0,       # +x 방향으로 부는 바람, rad
    release_angle=0.0,    # +x 방향 수평 방출, rad
)

print(result.report())
trace = result.receptor_trace((30.0, 0.0, 0.5))
print(trace.peak_mole_fraction, trace.peak_time_s)
```

`trace.time_s`는 방출 시작을 0초로 하는 절대 시각입니다. 정상 플룸은
물질 이동 시간이 `source_duration_s`에 도달하는 단면까지만 계산되고,
그 이후에는 정상 플룸을 외삽하지 않고 native puff가 계산합니다.

현재 puff 경로는 질량·H2 보존 및 수치 회귀시험을 통과했지만 독립적인
유한 LH2 현장 농도 검증은 아직 완료되지 않았으므로 결과에는 연구 경고가
유지됩니다. 측정 풍장이 방출 시작부터 전환 시점까지 포함되면 고수준 함수가
`assess_steady_wind_applicability`를 자동 실행합니다. 한계를 넘으면 결과는
조건부가 되고, `strict_scope=True`이면 정상 플룸 계산을 시작하기 전에
`ApplicabilityError`를 냅니다.

전환 이후 puff에는 측정 풍장을 직접 줄 수 있습니다. `WindHistory`의 시간은
방출 시작 기준 절대시각이며, 기록이 `source_duration_s`부터 puff 종료시각까지
이어져야 합니다.

```python
from degali.addons import WindHistory

measured_wind = WindHistory(
    time_s=[0.0, 5.0, 10.0, 20.0, 65.0],
    speed_m_s=[2.5, 2.5, 3.0, 2.0, 2.0],
    direction_from_deg=[270.0, 270.0, 250.0, 220.0, 220.0],
)

result = run_lh2_finite_release_research(
    source,
    source_duration_s=5.0,
    puff_duration_s=60.0,
    wind=2.5,
    puff_wind_history=measured_wind,
)
```

풍향은 기상학적 `from` 방위각으로 입력하고 내부에서는 동·북 벡터로
변환해 보간합니다. 관측 시각 경계에서 적분 step을 정확히 분할합니다.
전환 전 plume은 아직 `wind`와 `wind_angle`로 선언한 정상풍을 사용하므로,
그 구간의 풍향 변화가 크면 결과는 조건부 연구 결과입니다.

## 플래시 액적–rainout–풀 실행

`flashing_hydrogen_droplet_source()`로 만든 post-flash 입력은
`run_lh2_rainout_pool_research()`에 바로 전달할 수 있습니다. 이 함수는
잔류 액적을 크기군별로 횡풍·중력·증발 계산하고, 지면 도달분을 방출과
동시에 선언된 침적 footprint에 공급합니다. 기본값은 축대칭 얕은층
방정식으로 풀의 중력 확산과 각 반경 셀의 기판 접촉시간별 증발을 함께
계산합니다.

```python
from degali.addons import (
    DropletClass,
    SolidSubstrate,
    flashing_hydrogen_droplet_source,
)
from degali.lh2 import run_lh2_rainout_pool_research

flash = flashing_hydrogen_droplet_source(
    mass_flow=0.285,
    orifice_diameter=0.006,
    upstream_temperature=28.2550342766,  # 600 kPa 포화온도
    upstream_pressure=600000.0,
    upstream_quality=0.0,
)
ground = SolidSubstrate(
    conductivity_w_m_k=1.4,
    density_kg_m3=2200.0,
    heat_capacity_j_kg_k=850.0,
    initial_temperature_k=293.15,
    depth_m=0.5,
)

result = run_lh2_rainout_pool_research(
    flash,
    release_duration_s=5.0,
    post_release_duration_s=60.0,
    release_position_m=(0.0, 0.0, 1.5),
    release_azimuth_rad=0.0,
    release_elevation_rad=0.0,
    wind_speed_m_s=2.5,
    wind_to_angle_rad=0.0,
    evaporation_coefficient_m2_s=1.0e-9,
    pool_area_m2=1.0,
    pool_time_step_s=0.1,
    substrate=ground,
    droplet_classes=(
        DropletClass(diameter_m=1.0e-4, mass_fraction=0.25),
        DropletClass(diameter_m=1.0e-3, mass_fraction=0.75),
    ),
)
print(result.report())
```

`evaporation_coefficient_m2_s`와 `pool_area_m2`는 누출 조건만으로 유일하게
정해지지 않으므로 측정값 또는 선택한 상관식 결과를 명시해야 합니다.
여기서 `pool_area_m2`는 동적 풀의 최종 면적이 아니라 초기 지면 침적
footprint 면적입니다. 모델이 impact 점 몇 개에서 이 면적을 임의로
만들어내지는 않습니다. 계산영역 반경·격자·유입 환형을 직접 통제하려면
`DynamicPoolNumerics`를 `dynamic_pool_numerics=`로 전달하고, 물리적으로
고정된 방유제 면적을 계산할 때만 `pool_model="fixed"`를 명시합니다.
짧은 물/얼음 표면 시험처럼 독립적으로 열유속이 알려진 경우에는
`ConstantHeatFluxSurface`를 사용할 수 있습니다. 기본 증발 운동량 closure는
제한 검증에 사용된 `zero_radial_momentum_vapor`이며,
`liquid_velocity_carryoff`는 명시적 구조 민감도 옵션입니다.
