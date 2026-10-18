# TI021 수직 관측 연산자·풍향 진단 — 2026-10-18

## 짧은 물리 진단

동일한 10 s, 5 m inner 조건에서 60개 센서 trace를 비교했다. 아래 dose는
10 s 부분 적분이며, 600 s 관측 dose와 직접 동일한 성능 지표로 사용하지
않고 물리 후보의 방향성만 판단한다.

| 조건 | 최대 receptor 농도 (vol%) | 10 s 예측 dose 합 / 600 s 관측 dose 합 |
|---|---:|---:|
| `cell_center` | 9.3257 | 0.02645 |
| `ground_linear` | 13.2425 | 0.03077 |
| `ground_linear` v2 | 12.6933 | 0.02994 |
| finite-source 10 s | 6.2824 | 0.01245 |
| finite-source 4962 kg + `ground_linear` | 11.2154 | 0.02133 |

`ground_linear`은 저고도 관측 연산자로서 near-ground peak을 일부 회복하지만,
10 s에서 관측 전체 dose를 설명하지 못한다. finite-source 후보는 오히려
receptor dose를 낮췄으므로 현재 prior에서는 채택하지 않았다.

## 600 s 검증 실행

수직 연산자 효과를 최종 판단하기 위해 5 m, 450×640×40 m, Δt=0.125 s,
600 s의 `ground_linear` inner 실행을 시작했다.

```powershell
.\.venv\Scripts\python.exe tools/run_thorney_ti021_inner_transient.py `
  --workbook ti21_input_ascii.xls --duration-s 600 --cell-size-m 5 `
  --inner-domain-m 450 640 40 --time-step-s 0.125 `
  --receptor-vertical-sampling ground_linear `
  --output-directory outputs/thorney-ti021-inner-600s-groundlinear-2026-10-18
```

실행 중에는 결과를 승격하지 않는다. 완료 후 nested outer score를 동일한
600 s 조건으로 계산하고, 종·질량·운동량 보존 및 높이/거리 그룹 비악화를
모두 확인한다.

## 풍향 민감도

현재 최선 설정에서 풍향을 6.1°에서 3°로 바꾼 600 s 후보는 dose MAPE
**26.3115%**로 기준 26.1495%보다 악화되어 거절했다. 이는 단순 전역 풍향
회전만으로 10% 목표를 달성할 수 없음을 보여준다.
