"""Bounded NBS3921 excess-G extension; no field-fitted coefficient.

The temperature interpolation is explicit and caloric predictions from it
are NOT independently measured heat-of-mixing data. H2 is insoluble. Gas
remains ideal; mixed solids and excess volume are omitted. No frozen edits.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from . import liquid_phase_potential as lp
from .mixed_air_liquid import ideal_flash, AirLiquidFlash
from .oxygen_phase_potential import PhasePotential

MIN_T, MAX_T = 65., 77.5
GRID = (65., 70., 77.5)
TABLE = ((1.47, 1.38, 1.22), (32.58, 33.40, 34.73), (25.32, 25.82, 26.58))


def _domain(t, x, interpolation):
    if not math.isfinite(t) or not MIN_T <= t <= MAX_T:
        raise ValueError('excess liquid potential requires65--77.5K; no extrapolation')
    if not math.isfinite(x) or not 0 <= x <= 1:
        raise ValueError('liquid nitrogen mole fraction must lie in[0,1]')
    if interpolation not in ('quadratic', 'pchip'):
        raise ValueError('explicit quadratic or pchip temperature interpolation required')


@lru_cache(maxsize=2)
def _curves(interpolation):
    if interpolation == 'quadratic':
        return tuple(Polynomial.fit(GRID, v, 2) for v in TABLE)
    if interpolation == 'pchip':
        return tuple(PchipInterpolator(GRID, v, extrapolate=False) for v in TABLE)
    raise ValueError('unsupported interpolation')


def _parameter_derivatives(t, interpolation):
    result = []
    for curve in _curves(interpolation):
        derivative = curve.deriv if interpolation == 'quadratic' else curve.derivative
        result.append(tuple(float(curve(t) if k == 0 else derivative(k)(t)) for k in range(3)))
    return result


@dataclass(frozen=True)
class ExcessPotential:
    gibbs_J_mol: float
    enthalpy_J_mol: float
    entropy_J_mol_K: float
    heat_capacity_J_mol_K: float
    chemical_N2_J_mol: float
    chemical_O2_J_mol: float
    composition_curvature_J_mol: float


def excess(temperature, x_nitrogen, *, interpolation) -> ExcessPotential:
    _domain(temperature, x_nitrogen, interpolation)
    t, x = temperature, x_nitrogen
    (a, at, att), (vn, vnt, vntt), (vo, vot, vott) = _parameter_derivatives(t, interpolation)
    d, dt, dtt = x*vn+(1-x)*vo, x*vnt+(1-x)*vot, x*vntt+(1-x)*vott
    n = a*vn*vo
    nt = at*vn*vo+a*vnt*vo+a*vn*vot
    ntt = att*vn*vo+a*vntt*vo+a*vn*vott+2*(at*vnt*vo+at*vn*vot+a*vnt*vot)
    factor = 4.184*x*(1-x)
    g = factor*n/d
    gt = factor*(nt/d-n*dt/d**2)
    gtt = factor*(ntt/d-2*nt*dt/d**2-n*dtt/d**2+2*n*dt**2/d**3)
    pn, po = x*vn/d, (1-x)*vo/d
    return ExcessPotential(g, g-t*gt, -gt, -t*gtt,
        4.184*a*vn*po**2, 4.184*a*vo*pn**2, -2*4.184*a*vn**2*vo**2/d**3)


def activity_coefficients(temperature, x_nitrogen, *, interpolation):
    e = excess(temperature, x_nitrogen, interpolation=interpolation)
    return tuple(math.exp(mu/(lp.R*temperature)) for mu in (e.chemical_N2_J_mol,e.chemical_O2_J_mol))


def liquid_mixture(temperature, pressure, x_nitrogen, *, interpolation) -> PhasePotential:
    e = excess(temperature, x_nitrogen, interpolation=interpolation)
    old = lp.ideal_liquid_mixture(temperature, pressure, x_nitrogen)
    cp = old.heat_capacity_J_mol_K+e.heat_capacity_J_mol_K
    if cp <= 0: raise ValueError('nonpositive liquid mixture heat capacity')
    return PhasePotential(temperature, pressure, old.gibbs_J_mol+e.gibbs_J_mol,
        old.enthalpy_J_mol+e.enthalpy_J_mol, old.entropy_J_mol_K+e.entropy_J_mol_K,
        old.volume_m3_mol, cp)


def chemical_potentials(temperature, pressure, x_nitrogen, *, interpolation):
    e = excess(temperature, x_nitrogen, interpolation=interpolation)
    mu = lp.liquid_chemical_potentials(temperature, pressure, x_nitrogen)
    return mu[0]+e.chemical_N2_J_mol, mu[1]+e.chemical_O2_J_mol


def composition_curvature(temperature, x_nitrogen, *, interpolation):
    e = excess(temperature, x_nitrogen, interpolation=interpolation)
    x = x_nitrogen
    return math.inf if x in (0.,1.) else lp.R*temperature/(x*(1-x))+e.composition_curvature_J_mol


@dataclass(frozen=True)
class NonidealFlash:
    amounts: AirLiquidFlash
    incipient_liquid_nitrogen_fraction: float | None
    minimum_liquid_tangent_distance_RT: float | None
    chemical_equilibrium_error_RT: float
    interpolation: str


def flash(temperature, pressure, inert_mol, nitrogen_mol, oxygen_mol, *,
          interpolation, reservoir_vapor_fraction=0.) -> NonidealFlash:
    """Positive-inert equilibrium; activity and heat use the SAME excessG.

    The inert can include available water in the preliminary no-ice trial.
    Reservoir handling is inherited and its finite inventory is caller-owned.
    A stable all-gas state's hypothetical liquid composition is not inventory.
    The amounts.dew_index is a frozen-K RR diagnostic, NOT nonideal stability.
    """
    _domain(temperature, .5, interpolation)
    lp._pressure(pressure)
    values = (inert_mol,nitrogen_mol,oxygen_mol,reservoir_vapor_fraction)
    if (not all(math.isfinite(v) for v in values) or inert_mol <= 0
            or min(nitrogen_mol,oxygen_mol) < 0 or not 0 <= reservoir_vapor_fraction < 1):
        raise ValueError('positive inert, nonnegative air and reservoir fraction[0,1) required')
    kn, ko = (lp.equilibrium_partial_pressure(sp,temperature,pressure)/pressure
              for sp in ('Nitrogen','Oxygen'))
    kwargs = dict(reservoir_vapor_fraction=reservoir_vapor_fraction)
    if nitrogen_mol == 0 or oxygen_mol == 0:
        amounts = ideal_flash(inert_mol,nitrogen_mol,oxygen_mol,kn,ko,**kwargs)
        return NonidealFlash(amounts, None, None, 0., interpolation)
    total = inert_mol+nitrogen_mol+oxygen_mol
    zn, zo = nitrogen_mol/total, oxygen_mol/total
    q = reservoir_vapor_fraction
    lower, upper = np.nextafter(0.,1.), np.nextafter(1.,0.)

    def incipient_gradient(x):
        e = excess(temperature,x,interpolation=interpolation)
        return math.log(x)-math.log1p(-x)+math.log(kn/zn)-math.log(ko/zo)+(e.chemical_N2_J_mol-e.chemical_O2_J_mol)/(lp.R*temperature)
    xi = float(brentq(incipient_gradient,lower,upper,xtol=5e-324,rtol=1e-14))
    e = excess(temperature,xi,interpolation=interpolation)
    tpd = (xi*math.log(xi*kn/(zn*(1-q)))+(1-xi)*math.log((1-xi)*ko/(zo*(1-q)))+e.gibbs_J_mol/(lp.R*temperature))

    def at_x(x):
        gammas = activity_coefficients(temperature,x,interpolation=interpolation)
        ratios = (kn*gammas[0], ko*gammas[1])
        f = ideal_flash(inert_mol,nitrogen_mol,oxygen_mol,*ratios,**kwargs)
        beta = f.vapor_fraction_without_reservoir
        dn, do = (1+beta*(k/(1-q)-1) for k in ratios)
        residual = math.log(x)-math.log1p(-x)-math.log(zn/dn)+math.log(zo/do)
        return f, ratios, residual
    if tpd >= 0:
        f, _, _ = at_x(xi)
        if f.liquid_mol != 0: raise RuntimeError('stable gas and RR amount classifications disagree')
        return NonidealFlash(f,xi,tpd,0.,interpolation)
    x = float(brentq(lambda x:at_x(x)[2],lower,upper,xtol=5e-324,rtol=1e-14))
    f, ratios, _ = at_x(x)
    if f.liquid_mol <= 0: raise RuntimeError('unstable gas failed to form liquid')
    xl = f.liquid_nitrogen_mol/f.liquid_mol
    actual = (f.gas_nitrogen_mol/f.gas_mol, f.gas_oxygen_mol/f.gas_mol)
    expected = (xl*ratios[0], (1-xl)*ratios[1])
    error = max(abs(math.log(a/b)) for a,b in zip(actual,expected))
    if abs(x-xl)>1e-10 or error>1e-9:
        raise RuntimeError('nonideal chemical-equilibrium check failed')
    return NonidealFlash(f,xi,tpd,error,interpolation)
