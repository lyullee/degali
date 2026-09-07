# 서로 다른 열/성분 확산률의 엔탈피 수송 —2026-09-06

현재 연구 도구는 rho*D_H*grad(h)로 열량 확산을 표현하며 열/수소 확산률
비를 명시적으로 받는다. 실제 최신 약형/TKE 시도는 비율1을 사용했으므로
아래 보완이 그 결과를 바꾸거나 기존 잔차의 원인을 확정한다는 뜻은 아니다.
고정1차 모델도 변경하지 않는다.

## 보완할 구조

고정 압력에서 유효 이성분 h=h(T,Y)를 쓰면
grad(h)=h_T grad(T)+h_Y grad(Y)이다. 성분 교환은
j_Y=-rho*D_Y*grad(Y), 반대 성분의 유량은 -j_Y다. 열만의 수송을
-rho*D_T*h_T*grad(T)로 표현하려면 총 확산 엔탈피 유량은

```
J_H = -rho*D_T*(grad(h)-h_Y*grad(Y)) + h_Y*j_Y
    = -rho*D_T*grad(h) - rho*(D_Y-D_T)*h_Y*grad(Y).
```

gas ideal binary에서 h_Y=h_H2-h_air이며 성분이 지닌 엔탈피의 차이에
해당한다. 이는 추가 외부 열원이 아니고, 두 종의 순질량 확산 합은0이다.
같은 확산률에서는 기존 식과 정확히 같아야 한다.

이 식의 필수 검사는 성분 기준 엔탈피 변경이다. h'=h+aY+b일 때
grad(h')=grad(h)+a grad(Y), h'_Y=h_Y+a이므로
J'_H=J_H+a*j_Y여야 한다. h'_t=h_t+aY_t와 함께 쓰면 예측T는 변하지 않는다.
단순 -rho*D_T*grad(h)는 D_T!=D_Y일 때 이 조건을 만족하지 못한다.

## 구현 범위와 검사

독립 선택형 순수 함수로 구현한다. rho,D_Y,D_T,grad(h),grad(Y),h_Y는 모두
명시적 입력이며 출처 없는 난류 계수를 기본값으로 만들지 않는다. 성분
질량유량, 열만의 유량, 성분이 나르는 유량 및 합계를 분리해 반환한다.
단위, 같은확산률/등온/균일조성/0확산, 기준변경 불변성 및 닫힌2셀의
질량/에너지 보존과 T의 기준 불변성을 독립 검사한다.

혼상에서는 응축 입자가 어떤 속도로 움직이고 어떤 조성으로 교환되는지에
따라 부분엔탈피와 유효Y의 의미가 달라진다. 기체식 h_H2-h_air를59K 혼상에
그대로 적용하지 않는다. 새 함수는 사용자가 제공한 단일온도 이성분 접선
구성식에 조건부이며, 다성분 Maxwell--Stefan/Soret/Dufour 또는 별도 상의
완전한 수송 구현이라고 부르지 않는다. 압력 구배·상별 슬립도 별도다.

NASA TMR은 총에너지의 난류 열유량과 응력일을 구분한다. Cantera의
에너지/성분 방정식은 확산 성분과 엔탈피의 결합을 포함한다. 위 유효
이성분의 두 확산률 식 및 기준 불변성은 그 구조에 대한 본 프로젝트의
유도이며 두 문서의 LH2용 검증된 폐쇄상수라는 주장은 아니다.

- [NASA TMR](https://tmbwg.github.io/turbmodels/implementrans.html)
- [Cantera 에너지·성분 방정식](https://www.cantera.org/stable/reference/onedim/governing-equations.html)
- [Cantera 에너지 식 구현 설명](https://www.cantera.org/3.2/reference/onedim/discretization.html)

결과를 보고 확산률이나 h_Y를 조정하지 않는다. 이후 실제 하류 연결에는
변화된 총열 유량의 약형/외부 경계/성분 수송을 함께 재검증해야 한다.
