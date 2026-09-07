"""Opt-in gamma-solid O2 Gibbs potential; NOT a mixed-air phase model.

Caloric/volume data: Roder, NBSIR77-859 (1977), sections6.1--6.5.
Only43.801--54.361K and0<P<=200kPa. Saturation heat capacity is approximated
as reference isobaric Cp; the audit bounds that approximation. Incompressible
at fixed T. One potential generates h,s,Cp and equilibrium gas fugacity.
The frozen phase model and defaults are untouched.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

from scipy.optimize import brentq

MIN_T = 43.801
MAX_T = 54.361
MAX_P = 200_000.0
REFERENCE_GAS_PRESSURE = 100_000.0
FUSION_J_MOL = 106.3 * 4.184
CP_COEFFICIENTS = (16.908081 * 4.184, -0.24181777 * 4.184, 0.0024809089 * 4.184)
V_COEFFICIENTS = (23.2808187e-6, -0.06772142868e-6, 0.001339285715e-6)


def _check(temperature: float, pressure: float | None = None) -> None:
    if not math.isfinite(temperature) or not MIN_T <= temperature <= MAX_T:
        raise ValueError(f"gamma-O2 requires {MIN_T}--{MAX_T} K")
    if pressure is not None and (not math.isfinite(pressure) or not 0 < pressure <= MAX_P):
        raise ValueError(f"gamma-O2 requires 0<P<={MAX_P} Pa")


@dataclass(frozen=True)
class PhasePotential:
    temperature_K: float
    pressure_Pa: float
    gibbs_J_mol: float
    enthalpy_J_mol: float
    entropy_J_mol_K: float
    volume_m3_mol: float
    heat_capacity_J_mol_K: float


@lru_cache(maxsize=1)
def _triple_reference() -> dict[str, float]:
    """Current HEOS liquid triple reference, plus measured fusion heat."""
    import CoolProp as CP
    state = CP.AbstractState("HEOS", "Oxygen")
    state.update(CP.QT_INPUTS, 0.0, MAX_T)
    return dict(temperature_K=MAX_T, pressure_Pa=state.p(),
                liquid_h_J_mol=state.hmolar(), liquid_s_J_mol_K=state.smolar(),
                liquid_g_J_mol=state.gibbsmolar(), molecular_weight=state.molar_mass(),
                gas_constant=state.gas_constant(),
                solid_h_J_mol=state.hmolar()-FUSION_J_MOL,
                solid_s_J_mol_K=state.smolar()-FUSION_J_MOL/MAX_T)


def triple_reference() -> dict[str, float]:
    """Return a defensive copy; callers cannot mutate the reference cache."""
    return dict(_triple_reference())


def volume_and_derivatives(temperature: float) -> tuple[float, float, float]:
    _check(temperature)
    a, b, c = V_COEFFICIENTS
    return a+b*temperature+c*temperature**2, b+2*c*temperature, 2*c


def solid(temperature: float, pressure: float) -> PhasePotential:
    """Gamma-solid state from a single G(T,P), with molar SI units."""
    _check(temperature, pressure)
    ref = triple_reference()
    a, b, c = CP_COEFFICIENTS
    t, t0 = temperature, MAX_T
    dt = t-t0
    # Factor differences to preserve accuracy near the triple point.
    h0 = ref['solid_h_J_mol']+dt*(a+b*(t+t0)/2+c*(t*t+t*t0+t0*t0)/3)
    s0 = ref['solid_s_J_mol_K']+a*math.log1p(dt/t0)+b*dt+c*dt*(t+t0)/2
    v, dv, ddv = volume_and_derivatives(t)
    dp = pressure-ref['pressure_Pa']
    g = h0-t*s0+dp*v
    entropy = s0-dp*dv
    h = h0+dp*(v-t*dv)
    cp = a+b*t+c*t*t-t*dp*ddv
    return PhasePotential(t, pressure, g, h, entropy, v, cp)


@lru_cache(maxsize=8192)
def _ideal_reference(temperature: float) -> tuple[float, float, float, float]:
    import CoolProp as CP
    _check(temperature)
    state = CP.AbstractState("HEOS", "Oxygen")
    # h0 and s0 expose only the Helmholtz IDEAL contribution, even though
    # the disposable state itself includes residual terms.
    state.update(CP.DmolarT_INPUTS, 1.0, temperature)
    r = state.gas_constant()
    h = state.hmolar_idealgas()
    s = state.smolar_idealgas()+r*math.log(r*temperature/REFERENCE_GAS_PRESSURE)
    # cp0 from the ideal Helmholtz term, not actual Cp at density1mol/m3.
    cp = r*(1-state.tau()**2*state.d2alpha0_dTau2())
    return h, s, cp, r


def ideal_gas(temperature: float, partial_pressure: float) -> PhasePotential:
    """Matching ideal gas; p here is O2 partial pressure, not total pressure."""
    _check(temperature, partial_pressure)
    h, s0, cp, r = _ideal_reference(temperature)
    s = s0-r*math.log(partial_pressure/REFERENCE_GAS_PRESSURE)
    return PhasePotential(temperature, partial_pressure, h-temperature*s, h,
                          s, r*temperature/partial_pressure, cp)


def equilibrium_partial_pressure(temperature: float, total_pressure: float) -> float:
    """Gas partial pressure at equality to pure solid at mechanical totalP.

    This is the effective pure-solid fugacity in an ideal gas MIXTURE.
    It is allowed to exceed totalP (then that solid cannot be saturated).
    It is not an assertion of zero O2/N2 mutual solid or liquid solubility.
    """
    condensed = solid(temperature, total_pressure)
    h, s0, _, r = _ideal_reference(temperature)
    return REFERENCE_GAS_PRESSURE*math.exp((condensed.gibbs_J_mol-(h-temperature*s0))/(r*temperature))


def pure_vapor_pressure(temperature: float) -> float:
    """Pure ideal vapor pressure with the same solid mechanical pressure."""
    _check(temperature)
    root = brentq(lambda logp: math.log(equilibrium_partial_pressure(temperature, math.exp(logp)))-logp,
                  math.log(1e-8), math.nextafter(math.log(MAX_P), -math.inf), xtol=2e-13)
    return math.exp(root)


def ideal_gas_to_solid_enthalpy(temperature: float, total_pressure: float) -> float:
    """Molar enthalpy decrement for an IDEAL gas ledger, including reference."""
    return _ideal_reference(temperature)[0]-solid(temperature,total_pressure).enthalpy_J_mol


def nbs1977_vapor_pressure(temperature: float) -> float:
    """Published comparator only; natural logarithm, P initially in mmHg.

    Not used to construct this module's equilibrium or enthalpy.
    """
    if not math.isfinite(temperature) or not MIN_T <= temperature <= 54.359:
        raise ValueError('published gamma-O2 pressure comparator is43.801--54.359K')
    return 133.322368*math.exp(-1096.562485/temperature-2.025578307*math.log(temperature)+28.35976524)
