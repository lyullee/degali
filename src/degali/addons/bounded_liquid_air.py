"""Bounded liquid-only air property adapter for a partial downstream contrast.

N2/O2 common liquid G, legacy H2/gas h interpolation and water phase ledger.
No mixed solids or sub64K continuation; no complete upstream-EOS claim.
"""
from __future__ import annotations

from functools import lru_cache
import numpy as np
from scipy.interpolate import CubicHermiteSpline

from . import liquid_phase_potential as lp
from .axisymmetric_jet import _air_phase_property_table

MIN_T = 64.0
MAX_T = 300.0
MW = np.array([.0280134,.0319988])


def binary_flash(inert, nitrogen, oxygen, kn, ko, reservoir=0.):
    """Vectorized exact quadratic form of inert/N2/O2 Rachford-Rice.

    Returns gasN,gasO,liquidN,liquidO,reservoir vapor. No fixed trace cutoff.
    Inputs use any common mole scale. Reservoir inventory remains caller-owned.
    """
    inert,nitrogen,oxygen,kn,ko,q = np.broadcast_arrays(*[
        np.asarray(x,dtype=float) for x in (inert,nitrogen,oxygen,kn,ko,reservoir)])
    if (not all(np.all(np.isfinite(a)) for a in (inert,nitrogen,oxygen,kn,ko,q))
            or np.any(inert<0) or np.any(nitrogen<0) or np.any(oxygen<0)
            or np.any(kn<=0) or np.any(ko<=0) or np.any(q<0) or np.any(q>=1)):
        raise ValueError('invalid binary flash inputs')
    total = inert+nitrogen+oxygen
    if np.any(total<=0): raise ValueError('positive mole inventory required')
    zi,zn,zo = inert/total,nitrogen/total,oxygen/total
    k1,k2 = kn/(1-q),ko/(1-q)
    a,b = k1-1,k2-1
    dew = zn/k1+zo/k2
    liquid_only = (zi==0) & (zn*a+zo*b<=0) & (dew>1)
    two = (dew>1) & ~liquid_only
    beta = np.where(liquid_only,0.,1.)
    if np.any(two):
        aa,bb,ii = (a*b)[two],(zi*(a+b)+zn*a+zo*b)[two],zi[two]
        disc = bb*bb-4*aa*ii
        if np.any(disc < -1e-13*np.maximum(bb*bb,1)):
            raise RuntimeError('negative RR discriminant')
        radical = np.sqrt(np.maximum(disc,0))
        stable_q = -.5*(bb+np.copysign(radical,bb))
        root = np.empty_like(ii)
        use_small = bb<0
        # In the two-phase domain, the sign selects the root inside[zi,1].
        root[use_small] = ii[use_small]/stable_q[use_small]
        root[~use_small] = stable_q[~use_small]/aa[~use_small]
        # Noninert two-phase mixtures have a zero extraneous polynomial root.
        no_inert = ii==0
        root[no_inert] = -bb[no_inert]/aa[no_inert]
        if np.any(root < zi[two]-2e-13) or np.any(root>1+2e-13):
            raise RuntimeError('quadratic flash root outside physical interval')
        beta[two] = np.clip(root,zi[two],1.)
    gn,go = np.where(liquid_only,0.,nitrogen),np.where(liquid_only,0.,oxygen)
    ln,lo = np.where(liquid_only,nitrogen,0.),np.where(liquid_only,oxygen,0.)
    if np.any(two):
        # Avoid0/0 for an absent species whose tiny K rounds K-1 to-1
        # in a gas-only endpoint. Evaluate denominators only in two phases.
        den_n,den_o = 1+beta[two]*a[two],1+beta[two]*b[two]
        gn[two] = total[two]*beta[two]*zn[two]*k1[two]/den_n
        go[two] = total[two]*beta[two]*zo[two]*k2[two]/den_o
        ln[two] = total[two]*(1-beta[two])*zn[two]/den_n
        lo[two] = total[two]*(1-beta[two])*zo[two]/den_o
    reservoir_vapor = (inert+gn+go)*q/(1-q)
    if (np.any(abs(gn+ln-nitrogen)>2e-12*total)
            or np.any(abs(go+lo-oxygen)>2e-12*total)):
        raise RuntimeError('flash component conservation failed')
    return gn,go,ln,lo,reservoir_vapor


@lru_cache(maxsize=8)
def _potential_tables(pressure: float):
    if not 0<pressure<=lp.MAX_P: raise ValueError('unsupported liquid adapter pressure')
    grid = np.linspace(MIN_T,lp.MAX_T,1001)
    result = []
    for sp in ('Nitrogen','Oxygen'):
        l = [lp.liquid(sp,float(t),pressure) for t in grid]
        g = [lp.ideal_standard(sp,float(t)) for t in grid]
        vol = [lp.saturated_reference(sp,float(t)) for t in grid]
        gl = CubicHermiteSpline(grid,[s.gibbs_J_mol for s in l],[-s.entropy_J_mol_K for s in l],extrapolate=False)
        gg = CubicHermiteSpline(grid,[s[0]-t*s[1] for s,t in zip(g,grid)],[-s[1] for s in g],extrapolate=False)
        v = CubicHermiteSpline(grid,[s.volume_m3_mol for s in vol],[s.volume_T_m3_mol_K for s in vol],extrapolate=False)
        if lp.equilibrium_partial_pressure(sp,lp.MAX_T,pressure)<=pressure:
            raise ValueError('100K gas-only branch is not allowed at this pressure')
        result.append((gl,gg,v))
    return tuple(result)


class BoundedLiquidAir:
    """Delegating wrapper, never mutates the original atmosphere or model."""
    def __init__(self, original):
        if (not original.consistent_phase_ambient or original._dry_argon_mass_fraction != 0
                or original._phase_enthalpy_tables is None or not original.equilibrium_air_condensation):
            raise ValueError('requires explicit N2/O2/H2O temperature-dependent atmosphere')
        self.original = original
        self.tables = _potential_tables(original.ambient_pressure)
        self.phase_table = _air_phase_property_table()
        self.domain_failures = 0
        self.last_domain_failure = None
        hot = self.phase_table['temperature']>=lp.MAX_T
        if any(np.any(self.phase_table[f'{sp}_saturation'][hot]<=original.ambient_pressure)
               for sp in ('nitrogen','oxygen')):
            raise ValueError('legacy warm table permits dry air condensation')

    def __getattr__(self,name):
        return getattr(self.original,name)

    def air_properties(self,t):
        t = np.asarray(t,dtype=float)
        if np.any(t<MIN_T) or np.any(t>lp.MAX_T):
            raise ValueError('air liquid potential interpolation outside64--100K')
        k,h,v = [],[],[]
        for gl,gg,vs in self.tables:
            dg = gl(t)-gg(t)
            k.append(lp.P0*np.exp(dg/(lp.R*t))/self.ambient_pressure)
            h.append(-dg+t*(gl(t,1)-gg(t,1)))
            v.append(vs(t))
        return np.stack(k),np.stack(h),np.stack(v)

    def forward(self,temperature,fraction):
        t,y = np.broadcast_arrays(np.asarray(temperature,dtype=float),np.asarray(fraction,dtype=float))
        if (not np.all(np.isfinite(t)) or not np.all(np.isfinite(y))
                or np.any(t<MIN_T) or np.any(t>MAX_T) or np.any(y<0) or np.any(y>1)):
            raise ValueError('bounded liquid forward requires64--300K and0<=Y<=1')
        shape = t.shape
        t,y = t.ravel(),y.ravel()
        dry = (1-y)/(1+self.ambient_absolute_humidity)
        n = dry*self._dry_nitrogen_mass_fraction/MW[0]
        o = dry*self._dry_oxygen_mass_fraction/MW[1]
        water = dry*self.ambient_absolute_humidity/self.water_molecular_weight
        inert = y/self.fuel_molecular_weight
        pt = self.phase_table
        qw = np.interp(t,pt['temperature'],pt['water_saturation'])/self.ambient_pressure
        gn,go,gw = n.copy(),o.copy(),water.copy()
        ln,lo = np.zeros_like(t),np.zeros_like(t)
        latent = np.zeros((2,len(t)))
        volumes = np.zeros((2,len(t)))
        cold = t<=lp.MAX_T
        if np.any(cold):
            k,h,v = self.air_properties(t[cold])
            latent[:,cold],volumes[:,cold] = h,v
            cn,co,cln,clo,_ = binary_flash(inert[cold]+water[cold],n[cold],o[cold],*k)
            gas_total = inert[cold]+water[cold]+cn+co
            ice = water[cold]>qw[cold]*gas_total
            if np.any(ice):
                inn,io,iln,ilo,iw = binary_flash(inert[cold][ice],n[cold][ice],o[cold][ice],
                                               k[0,ice],k[1,ice],qw[cold][ice])
                cn[ice],co[ice],cln[ice],clo[ice] = inn,io,iln,ilo
                wsub = water[cold].copy()
                wsub[ice] = iw
                gw[cold] = wsub
            gn[cold],go[cold],ln[cold],lo[cold] = cn,co,cln,clo
        hot_ice = ~cold & (qw<1) & (water>qw*(inert+n+o+water))
        gw[hot_ice] = qw[hot_ice]*(inert[hot_ice]+n[hot_ice]+o[hot_ice])/(1-qw[hot_ice])
        if np.any(gw>water+1e-12*np.maximum(inert+n+o+water,1)):
            raise RuntimeError('finite water inventory exceeded')
        ice_mass = (water-gw)*self.water_molecular_weight
        water_v = ice_mass/np.interp(t,pt['temperature'],pt['water_density'])
        gas_mol = inert+gn+go+gw
        volume = gas_mol*lp.R*t/self.ambient_pressure+ln*volumes[0]+lo*volumes[1]+water_v
        h = self._mixture_enthalpy(t,y)-ln*latent[0]-lo*latent[1]-ice_mass*np.interp(t,pt['temperature'],pt['water_latent'])
        return dict(density=(1/volume).reshape(shape),enthalpy=h.reshape(shape),
                    liquid_N2=(ln*MW[0]).reshape(shape),liquid_O2=(lo*MW[1]).reshape(shape),
                    ice_water=ice_mass.reshape(shape),gas_moles=gas_mol.reshape(shape))

    def _condensed_air_state(self,density,fuel_fraction):
        rho,y = np.broadcast_arrays(np.asarray(density,dtype=float),np.asarray(fuel_fraction,dtype=float))
        if np.any(~np.isfinite(rho)) or np.any(rho<=0): raise ValueError('positive finite density required')
        low,high = np.full_like(rho,MIN_T),np.full_like(rho,MAX_T)
        rlow,rhigh = self.forward(low,y)['density'],self.forward(high,y)['density']
        eps = 8*np.finfo(float).eps
        invalid = (rho>rlow*(1+eps)) | (rho<rhigh*(1-eps))
        if np.any(invalid):
            self.domain_failures += 1
            ix = int(np.flatnonzero(invalid)[0])
            self.last_domain_failure = dict(density=float(rho.flat[ix]),fraction=float(y.flat[ix]),
                maximum_density_at_64K=float(rlow.flat[ix]),minimum_density_at_300K=float(rhigh.flat[ix]))
            raise ValueError(f'bounded liquid density outside64--300K:{self.last_domain_failure}')
        # Near an air-rich vaporization transition, tiny T errors amplify in
        # h.42steps retain the original enthalpy round-trip test tolerance.
        for _ in range(42):
            mid = .5*(low+high)
            dense = self.forward(mid,y)['density']>rho
            low,high = np.where(dense,mid,low),np.where(dense,high,mid)
        t = .5*(low+high)
        h = self.forward(t,y)['enthalpy']
        ambient = (y==0) & np.isclose(rho,self.ambient_density,rtol=eps,atol=0)
        return np.where(ambient,self.ambient_temperature,t),np.where(ambient,0.,rho*h)

    _condensed_air_state_exact = _condensed_air_state
