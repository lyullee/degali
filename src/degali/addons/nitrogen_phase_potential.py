"""Bounded beta-N2 Gibbs potential using Ziegler & Mullins(1963) calorics.

This is a pure solid property module, not a N2/O2 solid-solution model.
The source's mean solid volume28.67cm3/mol neglects thermal expansion and
compressibility. Beta domain35.62--63.151K,0<P<=200kPa only. No core changes.
"""
from __future__ import annotations

from functools import lru_cache
import math

from scipy.optimize import brentq

from .oxygen_phase_potential import PhasePotential

MIN_T = 35.62
MAX_T = 63.151
MAX_P = 200_000.0
REFERENCE_GAS_PRESSURE = 100_000.0
FUSION_J_MOL = 172.3 * 4.184
VOLUME_M3_MOL = 28.67e-6
# TableIV headers multiply A1,A2,A3,A4 by10^1,10^2,10^3,10^5.
CP_COEFFICIENTS = tuple(4.184*x for x in (
    24.326576,-16.223561e-1,5.7461088e-2,-0.84664594e-3,0.46237120e-5))


def _check(temperature: float, pressure: float | None = None) -> None:
    if not math.isfinite(temperature) or not MIN_T <= temperature <= MAX_T:
        raise ValueError(f"beta-N2 requires {MIN_T}--{MAX_T} K")
    if pressure is not None and (not math.isfinite(pressure) or not 0 < pressure <= MAX_P):
        raise ValueError(f"beta-N2 requires 0<P<={MAX_P} Pa")


@lru_cache(maxsize=1)
def _triple_reference() -> dict[str,float]:
    import CoolProp as CP
    state = CP.AbstractState('HEOS','Nitrogen')
    state.update(CP.QT_INPUTS,0.,MAX_T)
    return dict(temperature_K=MAX_T,pressure_Pa=state.p(),
        liquid_h_J_mol=state.hmolar(),liquid_s_J_mol_K=state.smolar(),
        liquid_g_J_mol=state.gibbsmolar(),molecular_weight=state.molar_mass(),
        gas_constant=state.gas_constant(),
        solid_h_J_mol=state.hmolar()-FUSION_J_MOL,
        solid_s_J_mol_K=state.smolar()-FUSION_J_MOL/MAX_T)


def triple_reference() -> dict[str,float]:
    return dict(_triple_reference())


def _power_difference(t: float, t0: float, exponent: int) -> float:
    # Exact factorization avoids cancellation for T extremely near Tt.
    return (t-t0)*math.fsum(t**j*t0**(exponent-1-j) for j in range(exponent))


def solid(temperature: float, pressure: float) -> PhasePotential:
    _check(temperature,pressure)
    t,t0 = temperature,MAX_T
    ref = _triple_reference()
    h0 = ref['solid_h_J_mol']+math.fsum(a*_power_difference(t,t0,j+1)/(j+1)
                                       for j,a in enumerate(CP_COEFFICIENTS))
    s0 = ref['solid_s_J_mol_K']+CP_COEFFICIENTS[0]*math.log1p((t-t0)/t0)+math.fsum(
        a*_power_difference(t,t0,j)/j for j,a in enumerate(CP_COEFFICIENTS) if j)
    h = h0+(pressure-ref['pressure_Pa'])*VOLUME_M3_MOL
    cp = math.fsum(a*t**j for j,a in enumerate(CP_COEFFICIENTS))
    return PhasePotential(t,pressure,h-t*s0,h,s0,VOLUME_M3_MOL,cp)


@lru_cache(maxsize=8192)
def _ideal_reference(temperature: float) -> tuple[float,float,float,float]:
    import CoolProp as CP
    _check(temperature)
    state = CP.AbstractState('HEOS','Nitrogen')
    state.update(CP.DmolarT_INPUTS,1.,temperature)
    r = state.gas_constant()
    s = state.smolar_idealgas()+r*math.log(r*temperature/REFERENCE_GAS_PRESSURE)
    cp = r*(1-state.tau()**2*state.d2alpha0_dTau2())
    return state.hmolar_idealgas(),s,cp,r


def ideal_gas(temperature: float, partial_pressure: float) -> PhasePotential:
    _check(temperature,partial_pressure)
    h,s0,cp,r = _ideal_reference(temperature)
    s = s0-r*math.log(partial_pressure/REFERENCE_GAS_PRESSURE)
    return PhasePotential(temperature,partial_pressure,h-temperature*s,h,s,
                          r*temperature/partial_pressure,cp)


def equilibrium_partial_pressure(temperature: float, total_pressure: float) -> float:
    """Effective solid fugacity in an IDEAL mixture at mechanical totalP."""
    s = solid(temperature,total_pressure)
    h,s0,_,r = _ideal_reference(temperature)
    return REFERENCE_GAS_PRESSURE*math.exp((s.gibbs_J_mol-h+temperature*s0)/(r*temperature))


def pure_vapor_pressure(temperature: float) -> float:
    _check(temperature)
    root = brentq(lambda logp: math.log(equilibrium_partial_pressure(temperature,math.exp(logp)))-logp,
                  math.log(1e-8),math.nextafter(math.log(MAX_P),-math.inf),xtol=2e-13)
    return math.exp(root)


def ideal_gas_to_solid_enthalpy(temperature: float, total_pressure: float) -> float:
    return _ideal_reference(temperature)[0]-solid(temperature,total_pressure).enthalpy_J_mol
