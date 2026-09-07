"""Common-reference pure/ideal-mixed liquid potential, an opt-in diagnostic.

First-order pressure expansion from HEOS stable saturated liquid. This is
NOT a complete air phase diagram: no subtriple liquid, mixed solids, excess
mixing enthalpy or dissolved H2. No frozen field/model defaults are changed.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

from .oxygen_phase_potential import PhasePotential
from .mixed_air_liquid import ideal_flash, AirLiquidFlash

R = 8.31446261815324
P0 = 100_000.0
MAX_P = 200_000.0
MAX_T = 100.0
TRIPLE_T = {'Nitrogen': 63.151, 'Oxygen': 54.361}
GAS_SPECIES = ('Hydrogen', 'Nitrogen', 'Oxygen')


def _pressure(pressure: float):
    if not math.isfinite(pressure) or not 0 < pressure <= MAX_P:
        raise ValueError('potential requires0<P<=200000Pa')


def _liquid_domain(species: str, temperature: float):
    if species not in TRIPLE_T:
        raise ValueError('liquid species must be Nitrogen or Oxygen')
    if not math.isfinite(temperature) or not TRIPLE_T[species] <= temperature <= MAX_T:
        raise ValueError(f'{species} liquid reference requires{TRIPLE_T[species]}--{MAX_T}K')


@dataclass(frozen=True)
class SaturatedReference:
    temperature_K: float
    pressure_Pa: float
    gibbs_J_mol: float
    enthalpy_J_mol: float
    entropy_J_mol_K: float
    volume_m3_mol: float
    volume_T_m3_mol_K: float
    volume_TT_m3_mol_K2: float
    pressure_T_Pa_K: float
    enthalpy_T_J_mol_K: float
    molecular_weight_kg_mol: float
    backend_gas_constant: float


@lru_cache(maxsize=8192)
def saturated_reference(species: str, temperature: float) -> SaturatedReference:
    import CoolProp as CP
    _liquid_domain(species, temperature)
    st = CP.AbstractState('HEOS', species)
    st.update(CP.QT_INPUTS, 0., temperature)
    rho = st.rhomolar()
    # CoolProp8 only implements second saturation derivatives w.r.t. P,P.
    tp = st.first_saturation_deriv(CP.iT, CP.iP)
    tpp = st.second_saturation_deriv(CP.iT, CP.iP, CP.iP)
    rp = st.first_saturation_deriv(CP.iDmolar, CP.iP)
    rpp = st.second_saturation_deriv(CP.iDmolar, CP.iP, CP.iP)
    pt = 1/tp
    ptt = -tpp/tp**3
    rt = rp*pt
    rtt = rpp*pt**2+rp*ptt
    vt = -rt/rho**2
    vtt = 2*rt**2/rho**3-rtt/rho**2
    return SaturatedReference(temperature, st.p(), st.gibbsmolar(), st.hmolar(),
        st.smolar(), 1/rho, vt, vtt, pt,
        st.first_saturation_deriv(CP.iHmolar, CP.iT), st.molar_mass(), st.gas_constant())


def liquid(species: str, temperature: float, pressure: float) -> PhasePotential:
    """One G generates every molar property; total mechanical P, not partialP."""
    _pressure(pressure)
    ref = saturated_reference(species, temperature)
    t, dp = temperature, pressure-ref.pressure_Pa
    v, vt, vtt = ref.volume_m3_mol, ref.volume_T_m3_mol_K, ref.volume_TT_m3_mol_K2
    g = ref.gibbs_J_mol+dp*v
    s = ref.entropy_J_mol_K-dp*vt
    h = ref.enthalpy_J_mol+dp*(v-t*vt)
    cp = ref.enthalpy_T_J_mol_K-ref.pressure_T_Pa_K*(v-t*vt)-t*dp*vtt
    if not all(math.isfinite(z) for z in (g, s, h, cp, v)) or cp <= 0 or v <= 0:
        raise ValueError('liquid pressure expansion is not thermally admissible here')
    return PhasePotential(t, pressure, g, h, s, v, cp)


@lru_cache(maxsize=8192)
def ideal_standard(species: str, temperature: float) -> tuple[float, float, float, float]:
    """Ideal h,s atP0 and Cp, molecular weight; no mutable cached objects."""
    import CoolProp as CP
    if species not in GAS_SPECIES:
        raise ValueError('ideal species must be Hydrogen, Nitrogen or Oxygen')
    if not math.isfinite(temperature) or not 35.62 <= temperature <= 300.:
        raise ValueError('ideal reference requires35.62--300K')
    st = CP.AbstractState('HEOS', species)
    st.update(CP.DmolarT_INPUTS, 1., temperature)
    rb = st.gas_constant()
    return (st.hmolar_idealgas(), st.smolar_idealgas()+rb*math.log(rb*temperature/P0),
            rb*(1-st.tau()**2*st.d2alpha0_dTau2()), st.molar_mass())


def ideal_gas(species: str, temperature: float, partial_pressure: float) -> PhasePotential:
    _pressure(partial_pressure)
    h, s0, cp, _ = ideal_standard(species, temperature)
    s = s0-R*math.log(partial_pressure/P0)
    return PhasePotential(temperature, partial_pressure, h-temperature*s, h, s,
                          R*temperature/partial_pressure, cp)


def equilibrium_partial_pressure(species: str, temperature: float, total_pressure: float) -> float:
    gl = liquid(species, temperature, total_pressure).gibbs_J_mol
    h, s0, _, _ = ideal_standard(species, temperature)
    return P0*math.exp((gl-(h-temperature*s0))/(R*temperature))


def ideal_gas_to_liquid_enthalpy(species: str, temperature: float, total_pressure: float) -> float:
    return ideal_standard(species, temperature)[0]-liquid(species, temperature, total_pressure).enthalpy_J_mol


def _composition(x_nitrogen: float):
    if not math.isfinite(x_nitrogen) or not 0 <= x_nitrogen <= 1:
        raise ValueError('liquid mole fraction must be in[0,1]')


def ideal_liquid_mixture(temperature: float, pressure: float, x_nitrogen: float) -> PhasePotential:
    """Fixed-composition molar G; ideal mixing has zero excess h,v,Cp."""
    _composition(x_nitrogen)
    states = (liquid('Nitrogen', temperature, pressure), liquid('Oxygen', temperature, pressure))
    x = (x_nitrogen, 1-x_nitrogen)
    mixing = math.fsum(z*math.log(z) for z in x if z > 0)
    average = lambda name: math.fsum(z*getattr(st, name) for z, st in zip(x, states))
    g = average('gibbs_J_mol')+R*temperature*mixing
    s = average('entropy_J_mol_K')-R*mixing
    return PhasePotential(temperature, pressure, g, average('enthalpy_J_mol'), s,
                          average('volume_m3_mol'), average('heat_capacity_J_mol_K'))


def liquid_chemical_potentials(temperature: float, pressure: float, x_nitrogen: float) -> tuple[float,float]:
    _composition(x_nitrogen)
    return tuple(liquid(sp, temperature, pressure).gibbs_J_mol+R*temperature*math.log(x)
                 if x > 0 else -math.inf
                 for sp, x in zip(('Nitrogen','Oxygen'), (x_nitrogen,1-x_nitrogen)))


def flash(temperature: float, pressure: float, hydrogen_mol: float,
          nitrogen_mol: float, oxygen_mol: float, *, reservoir_vapor_fraction: float = 0.) -> AirLiquidFlash:
    """Ideal gas/liquid amounts with chemical-potential-derived K ratios.

    Reservoir handling is inherited: callers must verify finite water inventory.
    No stability assertion against omitted mixed-solid phases is made.
    """
    kn, ko = (equilibrium_partial_pressure(sp, temperature, pressure)/pressure
              for sp in ('Nitrogen','Oxygen'))
    return ideal_flash(hydrogen_mol, nitrogen_mol, oxygen_mol, kn, ko,
                       reservoir_vapor_fraction=reservoir_vapor_fraction)
