"""Opt-in six-flux yaw extension; bounded coflow, not a full RANS closure."""
from __future__ import annotations

from copy import copy
import math

import numpy as np
from scipy.optimize import brentq

from .energy_crosswind import IndependentEnergyCrosswind
from ..core.constants import VKC
from ..core.entrainment import phi

# rho, H2 mass fraction, sigma_y*sigma_z, pitch, yaw, excess axial velocity,
# horizontal travelled distance, Cartesian X, Y(right), Z(up).
RHO, YH2, AREA, PITCH, YAW, UC, LENGTH, X, Y, Z = range(10)


def direction(pitch, yaw):
    ct = math.cos(pitch)
    return np.array([ct*math.cos(yaw), ct*math.sin(yaw), math.sin(pitch)])


def normal_drag(wind, tangent, coefficient):
    wind, tangent = np.asarray(wind, float), np.asarray(tangent, float)
    if (wind.shape != (3,) or tangent.shape != (3,)
            or not np.all(np.isfinite([wind, tangent]))
            or abs(float(tangent @ tangent)-1) > 1e-10
            or not math.isfinite(coefficient) or coefficient < 0):
        raise ValueError('finite wind, unit tangent and nonnegative drag coefficient required')
    normal = wind-float(wind @ tangent)*tangent
    return coefficient*np.linalg.norm(normal)*normal


class YawedCrosswind:
    """Retain original thermochemistry/coefficients and add conserved Py/yaw.

    Ambient width growth uses horizontal travelled distance, not global X.
    No reverse axial ambient branch, resolved TKE or mixed-solid EOS is added.
    """

    def __init__(self, base, wind_angle):
        if type(base) is not IndependentEnergyCrosswind:
            raise ValueError('extension requires the explicit original density model')
        if not math.isfinite(wind_angle):
            raise ValueError('finite wind angle required')
        self.base = base
        self.wind_angle = float(wind_angle)

    @staticmethod
    def lift(old, yaw=0.0):
        old = np.asarray(old, float)
        IndependentEnergyCrosswind._physical(old)
        if not math.isfinite(yaw):
            raise ValueError('finite yaw required')
        rho, fraction, area, pitch, uc, length, z = old
        return np.array([rho, fraction, area, pitch, yaw, uc, length,
                         length*math.cos(yaw), length*math.sin(yaw), z])

    @staticmethod
    def proxy(state):
        state = np.asarray(state, float)
        if state.shape != (10,) or not np.all(np.isfinite(state)) or state[LENGTH] < 0:
            raise ValueError('ten finite states with nonnegative travelled distance required')
        old = state[[RHO, YH2, AREA, PITCH, UC, LENGTH, Z]]
        IndependentEnergyCrosswind._physical(old)
        return old

    def view(self, state):
        old = self.proxy(state)
        cosine = math.cos(float(state[YAW])-self.wind_angle)
        if cosine < -1e-14 and self.base._wind(old) > 0:
            raise ValueError('reverse ambient axial branch unsupported; do not clip wind')
        # A fresh shallow view changes no base attributes or thermodynamics.
        projected = copy(self.base)
        projected._wind = lambda probe: cosine*self.base._wind(probe)
        return projected, old

    def fluxes(self, state, *, quadrature_points=None):
        view, old = self.view(state)
        f = view.integral_fluxes(old, quadrature_points=quadrature_points)
        p = math.hypot(f.momentum_x, f.momentum_z)*direction(state[PITCH], state[YAW])
        return np.array([f.total_mass, f.contaminant_mass, *p, f.energy])

    def match(self, target, initial, *, tolerance=1e-9):
        target, guess = np.asarray(target, float), np.asarray(initial, float).copy()
        self.proxy(guess)
        if target.shape != (6,) or not np.all(np.isfinite(target)):
            raise ValueError('six finite target fluxes required')
        horizontal = math.hypot(target[2], target[3])
        if horizontal <= 0:
            raise ValueError('vertical or zero horizontal momentum is unsupported')
        guess[PITCH] = math.atan2(target[4], horizontal)
        principal_yaw = math.atan2(target[3], target[2])
        guess[YAW] += math.remainder(principal_yaw-guess[YAW], 2*math.pi)
        view, old = self.view(guess)
        target5 = target[[0, 1, 2, 4, 5]].copy()
        target5[2] = horizontal
        fitted = view._match_flux_array(target5, old, relative_tolerance=tolerance)
        guess[[RHO, YH2, AREA, PITCH, UC, LENGTH, Z]] = fitted
        scale = np.maximum(abs(target), [1e-12, 1e-12, 1, 1, 1, 1])
        if np.max(abs(self.fluxes(guess)-target)/scale) > tolerance:
            raise RuntimeError('six-flux projection failed')
        return guess

    def sources(self, state):
        view, old = self.view(state)
        rho, _fraction, area, pitch, uc, _length, z = old
        (_sy, sz, ua, perimeter, ground_width, cleared,
         alpha_wind, sya, sza, dsya, dsza) = self.base._geometry(old)
        k, th, rhoa = self.base.k, self.base.thermodynamics, self.base.rhoa
        st, ct = math.sin(pitch), math.cos(pitch)
        n = direction(pitch, state[YAW])
        wind = ua*np.array([math.cos(self.wind_angle), math.sin(self.wind_angle), 0.])
        parallel = float(wind @ n)
        perpendicular = np.linalg.norm(wind-parallel*n)
        ri = 9.81*abs(rho-rhoa)*math.sqrt(area)/max(rhoa*max(uc,1e-8)**2,1e-12)
        alpha = k.alfa1
        if k.plume_transition:
            alpha = k.ALPHA_JET+(k.ALPHA_PLUME-k.ALPHA_JET)*min(abs(ri)/k.RI_PLUME,1)**2
        alpha /= math.sqrt(max(rho/rhoa,1e-12))
        shear = alpha*uc*perimeter
        width = self.base.houf_velocity_width(old)
        centre_velocity = max(uc+parallel,1e-9)
        contrast = abs(rhoa-rho)
        buoyant = 0.
        if contrast > 1e-12 and st > 0:
            froude = centre_velocity**2*max(rho,1e-12)/(9.81*width*contrast)
            buoyant = th._buoyancy_coefficient/max(froude,1e-30)*2*math.pi*centre_velocity*width*st
        local = min(shear+buoyant, th.plume_entrainment_limit*2*math.pi*width*centre_velocity)
        cross = k.alfa2*parallel*perpendicular/ua*perimeter if ua > 0 else 0.
        ambient = k.profile_integral(0.)*ua*(sza*dsya+sya*dsza)
        ground = 0.
        if self.base.ground_interaction == 'surface_layer' and ground_width > 0:
            depth = max(z+k.delta*sz*ct,1e-9)
            ri_layer = 9.81*max(rho-rhoa,0)*depth/(rhoa*max(self.base.jetplume.ustar,1e-6)**2*k.delta)
            ground = VKC*self.base.jetplume.ustar*(1+alpha_wind)/phi(ri_layer,0.,3)*ground_width
        mass = rhoa*(local+cross+ambient+ground)
        momentum = mass*wind+normal_drag(wind,n,k.cd*perimeter*rhoa/2)
        momentum[2] += self.base.buoyancy_force(old)
        if self.base.ground_interaction != 'free' and cleared < 1 and momentum[2] > 0:
            momentum[2] *= cleared
        energy = .5*mass*ua**2 if self.base.energy_transport == 'total' else 0.
        result = np.array([mass,0.,*momentum,energy])
        if not np.all(np.isfinite(result)) or mass < 0:
            raise RuntimeError('nonphysical source vector')
        return result

    @staticmethod
    def path_rate(state):
        return np.r_[math.cos(state[PITCH]), direction(state[PITCH], state[YAW])]

    def solve(self, initial, *, distance, step=.02, checkpoint=None):
        """Conservative RK4 with strict flux inversion at every stage."""
        self.view(initial)
        if not all(math.isfinite(v) and v > 0 for v in (distance,step)):
            raise ValueError('positive finite distance and step required')
        count = math.ceil(distance/step)
        arc = np.linspace(0.,distance,count+1)
        states, fluxes, sources, cumulative = [np.array(initial)], [], [], [np.zeros(6)]
        fluxes.append(self.fluxes(initial))
        sources.append(self.sources(initial))

        def stage(target, guess, position):
            guess = guess.copy()
            guess[LENGTH:Z+1] = position
            return self.match(target,guess)

        def result():
            actual = np.asarray(fluxes)
            total = np.asarray(cumulative)
            scale = np.maximum.reduce([abs(actual),abs(total),np.ones_like(actual)])
            return dict(arc_length=arc[:len(states)],states=np.asarray(states),fluxes=actual,
                        sources=np.asarray(sources),cumulative_sources=total,
                        maximum_relative_balance_residual=float(np.max(abs(actual-actual[0]-total)/scale)))

        try:
            for i in range(count):
                h = arc[i+1]-arc[i]
                current, f = states[-1], fluxes[-1]
                s1, p1 = sources[-1], self.path_rate(current)
                b = stage(f+h/2*s1,current,current[6:10]+h/2*p1)
                s2, p2 = self.sources(b),self.path_rate(b)
                c = stage(f+h/2*s2,b,current[6:10]+h/2*p2)
                s3, p3 = self.sources(c),self.path_rate(c)
                d = stage(f+h*s3,c,current[6:10]+h*p3)
                s4, p4 = self.sources(d),self.path_rate(d)
                df = h/6*(s1+2*s2+2*s3+s4)
                dp = h/6*(p1+2*p2+2*p3+p4)
                new = stage(f+df,d,current[6:10]+dp)
                states.append(new)
                fluxes.append(self.fluxes(new))
                sources.append(self.sources(new))
                cumulative.append(cumulative[-1]+df)
                if checkpoint is not None and ((i+1)%10 == 0 or i+1 == count):
                    checkpoint(i+1,result(),None)
        except Exception as error:
            if checkpoint is not None:
                checkpoint(len(states)-1,result(),repr(error))
            raise
        return result()


class YawedTrajectory:
    """Original vertical profile at a uniquely bracketed horizontal-normal cut."""
    def __init__(self, model, states):
        self.model = model
        self.states = np.asarray(states,float)
        if (self.states.ndim != 2 or self.states.shape[1] != 10
                or len(self.states) < 2 or not np.all(np.isfinite(self.states))):
            raise ValueError('at least two finite ten-state sections required')

    def section(self,x,y):
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError('finite receptor coordinates required')
        states = self.states
        def residual(state):
            return (x-state[X])*math.cos(state[YAW])+(y-state[Y])*math.sin(state[YAW])
        values = np.array([residual(s) for s in states])
        exact = np.flatnonzero(values == 0.)
        brackets = np.flatnonzero(values[:-1]*values[1:] < 0.)
        if len(exact)+len(brackets) != 1:
            raise ValueError('receptor has no unique bracketed section; no extrapolation')
        if len(exact):
            state = states[exact[0]].copy()
        else:
            i = brackets[0]
            weight = brentq(lambda a: residual((1-a)*states[i]+a*states[i+1]),0.,1.,xtol=1e-13)
            state = (1-weight)*states[i]+weight*states[i+1]
        lateral = -(x-state[X])*math.sin(state[YAW])+(y-state[Y])*math.cos(state[YAW])
        return state,lateral

    def temperature_at(self,x,y,z):
        state,lateral = self.section(x,y)
        return self.model.base.point_temperature(self.model.proxy(state),lateral,z)

    def concentration_at(self,x,y,z):
        """Return mole fraction (0--1), not volume percent."""
        state,lateral = self.section(x,y)
        return self.model.base.point_mole_fraction(self.model.proxy(state),lateral,z)
