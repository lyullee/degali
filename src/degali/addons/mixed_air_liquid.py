"""Bounded ideal N2/O2 solution flash, separate from the frozen phase EOS.

This is a liquid-phase diagnostic, not an all-temperature plume closure.
Hydrogen is an insoluble inert gas; liquid mixing is ideal. No solid-liquid
extrapolation, dissolved hydrogen, excess enthalpy or fitted field parameter.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq

R = 8.31446261815324
MW_N = .0280134
MW_O = .0319988
MIN_T = 63.151
MAX_T = 120.


@dataclass(frozen=True)
class AirLiquidFlash:
    gas_inert_mol: float
    gas_nitrogen_mol: float
    gas_oxygen_mol: float
    reservoir_vapor_mol: float
    liquid_nitrogen_mol: float
    liquid_oxygen_mol: float
    vapor_fraction_without_reservoir: float
    dew_index: float
    maximum_inventory_error: float
    equilibrium_error: float

    @property
    def gas_mol(self):
        return self.gas_inert_mol+self.gas_nitrogen_mol+self.gas_oxygen_mol+self.reservoir_vapor_mol

    @property
    def liquid_mol(self):
        return self.liquid_nitrogen_mol+self.liquid_oxygen_mol


def ideal_flash(inert_mol: float, nitrogen_mol: float, oxygen_mol: float,
                nitrogen_k: float, oxygen_k: float, *, reservoir_vapor_fraction: float = 0.) -> AirLiquidFlash:
    """Flash given n_i and K_i=p_i_sat/P, with optional saturated ice vapor.

    A reservoir vapor fraction q uses V'=V*(1-q) and K'=K/(1-q).
    Reservoir inventory must be checked by the caller; no material is created
    by treating the returned vapor amount as a freely available water source.
    All mole inputs can use any common extensive basis; outputs use that basis.
    """
    values = (inert_mol, nitrogen_mol, oxygen_mol, nitrogen_k, oxygen_k, reservoir_vapor_fraction)
    if not all(math.isfinite(v) for v in values):
        raise ValueError('flash inputs must be finite')
    if min(inert_mol, nitrogen_mol, oxygen_mol) < 0 or min(nitrogen_k, oxygen_k) <= 0:
        raise ValueError('non-negative moles and positive K ratios required')
    q = reservoir_vapor_fraction
    if not 0 <= q < 1:
        raise ValueError('reservoir vapor fraction must lie in [0,1)')
    total = inert_mol+nitrogen_mol+oxygen_mol
    if total <= 0:
        raise ValueError('positive total inventory required')
    z = np.array([nitrogen_mol, oxygen_mol])/total
    zi = inert_mol/total
    k = np.array([nitrogen_k, oxygen_k])/(1-q)
    dew = float(np.sum(z/k))
    def rr(beta):
        return zi/beta + float(np.sum(z*(k-1)/(1+beta*(k-1))))
    if dew <= 1:
        beta = 1.
    elif zi == 0 and float(np.sum(z*(k-1))) <= 0:
        beta = 0.
    else:
        # With a noncondensable, beta >= zi. No arbitrary trace threshold.
        lower = zi if zi > 0 else 0.
        def objective(beta):
            if beta == 0:
                return float(np.sum(z*(k-1)))
            return rr(beta)
        beta = float(brentq(objective, lower, 1., xtol=5e-324, rtol=1e-14))
    liquid = total*(1-beta)*z/(1+beta*(k-1))
    # Compute both phases directly: subtracting two nearly equal inventories
    # loses the gas composition when almost all air is liquid.
    gas = total*beta*z*k/(1+beta*(k-1))
    gas_no_reservoir = inert_mol+float(np.sum(gas))
    vapor_reservoir = gas_no_reservoir*q/(1-q)
    inventory_error = float(max(abs(gas+liquid-total*z)))/total
    equilibrium_error = 0.
    if 0 < beta < 1:
        actual_y = gas/(gas_no_reservoir+vapor_reservoir)
        target_y = np.array([nitrogen_k, oxygen_k])*liquid/float(np.sum(liquid))
        equilibrium_error = float(max(abs(actual_y-target_y)))
    if min(*gas, *liquid) < -1e-12*total or inventory_error > 1e-12 or equilibrium_error > 1e-10:
        raise RuntimeError('air-liquid flash failed independent amount/equilibrium checks')
    return AirLiquidFlash(inert_mol, float(gas[0]), float(gas[1]), vapor_reservoir,
        float(liquid[0]), float(liquid[1]), beta, dew, inventory_error, equilibrium_error)


@lru_cache(maxsize=4096)
def liquid_properties(temperature: float):
    """Pure liquid references only where BOTH N2 and O2 have stable liquids."""
    from CoolProp.CoolProp import PropsSI
    if not math.isfinite(temperature) or not MIN_T <= temperature <= MAX_T:
        raise ValueError(f'liquid screen requires {MIN_T} <= T <= {MAX_T} K')
    result = {}
    for name in ('Nitrogen', 'Oxygen'):
        result[name] = dict(pressure=PropsSI('P', 'T', temperature, 'Q', 0, name),
            density=PropsSI('D', 'T', temperature, 'Q', 0, name),
            latent=PropsSI('H', 'T', temperature, 'Q', 1, name)-PropsSI('H', 'T', temperature, 'Q', 0, name))
    return result


def flash_at_temperature(temperature: float, pressure: float, inert_mol: float,
                         nitrogen_mol: float, oxygen_mol: float, *, reservoir_vapor_fraction: float = 0.):
    if not math.isfinite(pressure) or pressure <= 0:
        raise ValueError('positive finite pressure required')
    properties = liquid_properties(float(temperature))
    return ideal_flash(inert_mol, nitrogen_mol, oxygen_mol,
        properties['Nitrogen']['pressure']/pressure, properties['Oxygen']['pressure']/pressure,
        reservoir_vapor_fraction=reservoir_vapor_fraction)


def nbs_activity_coefficients(temperature: float, liquid_nitrogen_fraction: float):
    """Reproduce NBS3921 eq1a/1b with Tables4/6, NOT a new validation fit.

    Independent of this module's ideal flash. No temperature extrapolation.
    Molar volumes cm3/mol and A in cal/cm3 give RT in cal/mol.
    """
    if not math.isfinite(temperature) or not 65 <= temperature <= 77.5:
        raise ValueError('NBS correlation interpolation only covers65-77.5K')
    x = liquid_nitrogen_fraction
    if not math.isfinite(x) or not 0 <= x <= 1:
        raise ValueError('liquid mole fraction must lie in [0,1]')
    grid = [65., 70., 77.5]
    vn = float(np.interp(temperature, grid, [32.58, 33.40, 34.73]))
    vo = float(np.interp(temperature, grid, [25.32, 25.82, 26.58]))
    a = float(np.interp(temperature, grid, [1.47, 1.38, 1.22]))
    phi_n = x*vn/(x*vn+(1-x)*vo)
    return math.exp(vn*a*(1-phi_n)**2/(R/4.184*temperature)), math.exp(vo*a*phi_n**2/(R/4.184*temperature))
