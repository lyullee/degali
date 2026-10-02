# Pressure--strain / normal-stress status

Date: 2026-09-12

## What has been added

`reynolds_stress_redistribution` adds an opt-in, algebraically transparent
Rotta-type slow pressure--strain comparison term:

```text
Pi_ij = -C_phi (epsilon/k) [R_ij - 2 k delta_ij / 3].
```

It verifies `trace(R)=2k`, makes `trace(Pi)=0` exactly, and returns a maximum
forward time step for which `R + dt Pi` remains positive semidefinite. It
never clips or projects a non-realizable tensor. This is important because a
normal stress affects both axial momentum (`rho R_ss`) and stress work
(`rho U R_ss`), whereas pressure--strain redistributes existing fluctuation
energy rather than creating heat or TKE.

## Why it is not a default closure

The exact Reynolds-stress equation needs production, pressure--strain,
dissipation and turbulent/pressure transport. The classical
[Launder--Reece--Rodi paper](https://doi.org/10.1017/S0022112075001814)
explicitly treats pressure--strain as a closure problem and states its
high-Reynolds-number applicability. Publicly available cryogenic-hydrogen
material located so far provides no full normal-stress tensor, no pressure
correlation and no matching epsilon profile. A 2026 hydrogen crossflow DNS
comparison also reports that RANS can underpredict Reynolds stresses even
where mean flow is satisfactory; its geometry is not transferable to the
Sandia cryogenic free jet.

Therefore `C_phi`, `epsilon`, `k` and the full covariance are mandatory
inputs. The module is a conservation/realizability building block, not an
adopted model, a calibration, or evidence that a classical incompressible
pressure--strain law transfers through a cryogenic shock train.

## Next evidence gate

Before coupling this branch into `FiniteTkeModalTransport`, obtain or recover
at least matched axial and transverse RMS/covariance information and a
dissipation or length-scale observation. Until then, finite-TKE candidates
remain limited to the public lower-bound screen and declared end-member
sensitivity protocol.
