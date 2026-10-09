# 공개·로컬 실측 검증 데이터 인벤토리

작성일: 2026-10-08
목적: 누출·장애물/후류·센서 데이터의 실제 사용 가능성을 확인하고, 현재 `FieldEvidenceManifest`에 승격할 수 있는 범위를 고정한다.

## 이번 감사에서 확보한 데이터

| 데이터 | 원본 위치 | 실제 감사 결과 | 현재 사용 범위 |
|---|---|---|---|
| ELVHYS WP4.2 HSE | `SLABx_LH2/.../10.18710-JXJP0H/data` | 대표 9개 시험(3, 10, 24, 29, 30, 31, 36, 43, 48)을 concentration/temperature/pressure/misc까지 원본 CSV로 읽음. 시험 30은 100,001행, 시험 31은 273,649행. 100 kHz 점화 압력 파일의 예약 빈 행은 시간축 끝에서만 잘라냄. | 밀폐 transfer-connection 공간의 극저온 H₂ 관측. `FLMT`는 환기 유량이며 H₂ 질량유량이 아니다. outdoor LH₂ pool 정량 점수에는 사용하지 않음. |
| PRESLHY E3.4 pool | `SLABx_LH2/.../10.35097-1319/data/dataset` | 10개 XLSX 중 9개가 보정 `m(LH2) [g]` 및 열전대 기록으로 읽힘. Concrete01은 해당 보정 질량열이 `No Data!`라 입력을 승격하지 않음. 예: Gravel02의 호출자 선언 구간 4000–8000 s에서 질량감소율 `0.00102506 kg/s`, `R²=0.9878`이 계산됨(자동 구간선택 아님). | 액체 풀 질량손실·기판 열응답 관측. 풍속·농도·장애물과 결합한 야외 plume 점수에는 별도 입력이 필요. |
| PRESLHY E3.1 raw | `SLABx_LH2/.../10.35097-1317/raw` | 5개 `*-Final.xlsx`가 확보되어 있으나 현재 E3.1 reader가 기대하는 `*-Press`/`*-Temp` sheet 명명과 실제 `*-VD10`/`*-TE`가 다름. 원본은 보존하고 adapter 작업 대상으로 표시함. | 압력·온도 원본 후보. 현 상태에서는 자동 정량 결과를 만들지 않고 `withheld`. |
| AIJ UWE Case-H | `reference/obstacle_wake/aij_case_h` | `RS_caseH.csv` 876개 위치 레코드, `AF_caseH.csv` 17개 평균/난류 레코드, `LF_caseH.xlsx`(inflow/후류 단면 시트), 이미지와 README를 Zenodo API로 확보. 평균 농도의 배경보정 음수 96건을 허용하고 RMS/TKE만 음수 금지로 검증. | 단일 cuboid 후류·중성 ethylene plume의 풍동 관측. signed mean은 보존하지만 LH₂ 또는 현장 공통시간축으로 간주하지 않음. |
| USN open-channel H₂ | `reference/hydrogen/open_channel_usn` | Dataverse 공개 ZIP 22개와 ReadMe를 확보. 22/22 MD5 일치, 전부 29 센서열을 읽었고 flow meter 1/2 변형을 모두 처리. 각 시험은 flow clock과 sensor clock을 별도로 보존. 4% 센서 도달/피크/이탈 스크리닝을 완료. 해시 고정 timing audit는 22/22 파일, 638 센서, 455 threshold crossings를 기록했다. | 실제 H₂ 방출과 센서 시계열. 장애물·기상·현장 좌표계가 없어 standalone source/receptor 검증이며 후류 점수에는 사용하지 않음. `boundary_only`, `promotion_allowed=false`. |
| EPA irregular urban array | `reference/obstacle_wake/epa_irregular_array` | 공개 CSV 5개를 확보하여 7,272행, 결측 0, 좌표·속도/농도 필드를 확인. Base/Tank/TB65 농도장과 JRII Tank 속도장을 분리 보존. | 장애물 배열/후류의 tracer·풍동 검증. H₂/LH₂ 물성 검증이 아니므로 H₂ source score에는 직접 사용하지 않음. |

## 코드 보완 및 회귀검증

- ELVHYS 파일 매칭은 `HSE004` prefix glob을 제거하고, stream 앞의 전체 숫자 토큰을 정수로 해석한다. 따라서 `HSE0043`을 test 4로 오인하지 않고 `HSE0024`도 test 24로 읽는다.
- ELVHYS reader는 100 kHz 파일의 **끝부분 예약 빈 행만** 제거한다. 기록 중간의 빈 시간행은 계속 오류로 남겨 시간축 단절을 숨기지 않는다.
- 선택적으로 요구하지 않은 `FLMT`가 `Time`만 가진 placeholder이면 누출 데이터에서 환기 유량 부재로 처리한다. `FLMT`를 required로 선언하면 여전히 실패한다.
- AIJ reader는 풍동의 배경보정 평균 농도 음수를 유효한 signed mean으로 보존하고, `concentration_rms_ppm`만 비음수로 제한한다.
- USN reader는 `mass flow meter 1 [g/s]`와 `mass flow meter 2 [g/s]`를 모두 지원하며 실제 선택 열 이름을 결과에 남긴다.

실행한 집중 회귀검증:

```text
tests/test_elvhys.py                         6 passed
tests/test_obstacle_wake_validation.py       7 passed
tests/test_open_channel_h2.py                3 passed
USN extracted experiments                   22/22 parsed, 0 errors
전체 pytest 회귀                         1668 passed, 145 skipped
```

## `FieldEvidenceManifest` 승격 판정

이번에 확보한 공개·로컬 데이터는 모두 유용하지만, 하나의 운영용 `FieldEvidenceManifest`로 자동 승격할 수 있는 데이터는 아직 없다.

- **conditional**: ELVHYS, PRESLHY E3.4, USN — 실측 H₂ 또는 LH₂ 관측은 있으나 목적별로 질량유량/기상/장애물/공통시간축/운영자·교정 registry가 동시에 완비되지 않았다.
- **conditional (obstacle-only)**: AIJ, EPA — 장애물·후류 검증에는 유효하지만 중성가스/tracer 풍동 자료이며 LH₂ 누출 source 증거가 아니다.
- **withheld**: PRESLHY E3.1 raw — 원본은 확보했지만 현재 sheet schema adapter가 없어 자동 수치 결과를 발행하지 않는다.

따라서 이 자료들을 실제 현장 사건의 manifest로 만들 때는 원본을 재가공해 빈 칸을 채우지 말고 다음을 추가로 받아야 한다.

1. 사건 ID와 모든 장비의 공통 시간축/오프셋·동기화 근거
2. 누출 위치·형상·상태·질량유량(또는 승인된 압력 기반 입력)
3. 풍속·풍향·측정 높이·대기 안정도와 기상 장비 교정정보
4. 장애물 좌표계·위치·높이·형상 및 후류 측정 조건
5. 센서 좌표·높이·시계열·교정/검출한계·operator/curator ID
6. 원본 파일 SHA-256과 처리/변환 이력

위 항목 중 하나라도 누락되면 manifest는 `conditional` 또는 `withheld`로 유지하는 것이 현재 게이트 정책이다.
