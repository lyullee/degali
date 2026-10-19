# TI021 현재 코드 재현 및 물리 sweep — 2026-10-19

## 목적

과거 산출물의 구현 SHA가 현재 checkout과 달라지는 문제를 확인한 뒤, 현재
코드에서 같은 후보를 다시 실행하고 전역 후류 체류 파라미터만 바꿔 재검증했다.
센서별 보정계수와 관측 dose 피팅은 사용하지 않았다. 각 항목은 동일한 600 s
전체 이벤트를 다시 계산했다.

## 구현 재현성 감사

기존 기록의 `ti021-core-meander-s005-l150-w50`는 당시
`score_ti021_dense_layer_hybrid.py` SHA-256이
`6d731aa8238810a88b9422cba91a1a4c4193a7bf6986d0072ebf6380a83414d9`였지만,
현재 checkout은 `6bd7ccbc9805f9b2259ec5fea95de9633460dc85255a6c77d3eaa19c80c77e4c`다.
현재 코드로 동일 설정을 재실행한 결과는 **26.22907%**였으며, 종 보존 및
운동량 게이트는 통과했다. 따라서 이전 26.14946%를 현재 코드의 재현값으로
간주하지 않는다.

## 현재 코드 600 s 후보

| 후보 | dose MAPE | species residual | momentum gate |
|---|---:|---:|---|
| 재부착 0.40 + wake capture 0.025 s⁻¹ | **25.86961%** | `3.29e-15` | pass |
| 재부착 0.50 + wake capture 0.025 s⁻¹ | 25.89101% | `3.47e-15` | pass |
| 재부착 0.50 + wake capture 0.020 s⁻¹ | 25.92528% | `4.75e-15` | pass |

상단 이송 branch를 결합한 후속 600 s 재계산에서는
`wall_impulse_elevated_transport_fraction=0.20`이 **25.79593%**로 더 낮았다.
species residual은 `4.02e-15`, momentum gate는 pass였다. 횡방향 속도·중심
위치 변경은 각각 25.91228%, 26.31736%로 악화되어 채택하지 않았다.

최선 후보는 현재 checkout 기준의 새 최선이지만 10% 목표는 아직 미달성이다.
이 sweep은 장애물 후면 재순환과 재부착 체류라는 전역 물리를 변경했으며,
센서별 gain을 적용하지 않았다.

## 재현 명령

```powershell
.\.venv\Scripts\python.exe tools\sweep_ti021_current_physics.py `
  --base outputs/ti021-core-meander-s005-l150-w50-600s-2026-10-18.json `
  --output-directory outputs/ti021-current-physics-sweep
```

상세 실행 결과는 `outputs/ti021-physical-sweep2-2026-10-19/`에 보존되어
있고, 저장소에는 재현 스크립트와 이 요약만 커밋한다.
