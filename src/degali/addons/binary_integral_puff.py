"""Generic binary-gas bridge to the existing conservative integral puff.

Uses the actual finite-puff mass/momentum/energy integration, not an advected
Gaussian with prescribed widths. Species names are generic here even though
the older puff dataclasses call their conserved species field hydrogen_mass.
This bridge does not add obstacle physics to a free-field integral puff.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.special import ndtr

from .binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from .finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from .finite_release import FiniteReleasePuffHandoff
from .finite_puff import GaussianPuffConfig, GaussianPuffResult, integrate_finite_gaussian_puff


@dataclass(frozen=True)
class BinaryPuffThermodynamics:
    """Calorically-perfect ideal-gas adapter, with enthalpy zero at ambient T."""
    closure: BinaryIdealGasMixtureClosure

    @property
    def ambient_density(self):
        return self.closure.ambient_density_kg_m3

    @property
    def ambient_temperature(self):
        return self.closure.ambient_temperature_k

    @property
    def fuel_molecular_weight(self):
        return self.closure.species_molar_mass_kg_mol

    @property
    def _humid_ambient_molecular_weight(self):
        return self.closure.carrier_molar_mass_kg_mol

    @property
    def _ambient_enthalpy(self):
        return 0.

    def _temperature_from_enthalpy(self, enthalpy, fraction):
        cp = (fraction*self.closure.species_heat_capacity_j_kg_k
              +(1-fraction)*self.closure.carrier_heat_capacity_j_kg_k)
        temperature = self.ambient_temperature + enthalpy/cp
        if not self.closure.minimum_temperature_k <= temperature <= self.closure.maximum_temperature_k:
            raise RuntimeError('binary integral puff temperature outside closure bounds')
        return temperature

    def _density_from_temperature(self, temperature, fraction):
        reciprocal = (fraction/self.fuel_molecular_weight
                      +(1-fraction)/self._humid_ambient_molecular_weight)
        return self.closure.ambient_pressure_pa/(8.314462618*temperature*reciprocal)


def initial_prism_puff(cloud: RegularPolygonPrismInitialCloud,
                      closure: BinaryIdealGasMixtureClosure,
                      config: GaussianPuffConfig,
                      initial_velocity_m_s: tuple[float, float, float] = (0., 0., 0.)) -> GaussianPuffResult:
    """Match exact inventory, horizontal second moments and vertical bounds.

    A uniform prism has sigma_z=H/sqrt(12). Using cutoff=sqrt(3) places the
    initial grounded centre at H/2 without the old default 3-sigma geometric
    clamp shifting it upwards. This is a source-geometric choice, not a fit.
    Gaussian tails are an approximation to the initial top-hat holder.
    """
    if cloud.z_min_m != 0 or not math.isclose(config.width_cutoff_sigma,math.sqrt(3)):
        raise ValueError('grounded uniform prism requires z_min=0 and cutoff=sqrt(3)')
    thermo = BinaryPuffThermodynamics(closure)
    fraction = float(closure.mass_fraction_from_mole_fraction(cloud.species_mole_fraction))
    rho = thermo._density_from_temperature(cloud.temperature_k,fraction)
    mass = cloud.exact_volume_m3*rho
    species = mass*fraction
    # Polygon second moments by exact boundary integral (isotropic for a
    # regular n>=4 polygon), independent of source azimuth.
    p = cloud.vertices_xy_m()-np.array(cloud.center_xy_m)
    q = np.roll(p,-1,axis=0)
    cross = p[:,0]*q[:,1]-q[:,0]*p[:,1]
    sigma_xy = math.sqrt(float(np.sum(cross*(p[:,0]**2+p[:,0]*q[:,0]+q[:,0]**2))
                               /(12*cloud.horizontal_area_m2)))
    height = cloud.z_max_m
    cp = fraction*closure.species_heat_capacity_j_kg_k+(1-fraction)*closure.carrier_heat_capacity_j_kg_k
    initial_velocity = np.asarray(initial_velocity_m_s, dtype=float)
    if initial_velocity.shape != (3,) or not np.all(np.isfinite(initial_velocity)):
        raise ValueError('initial_velocity_m_s must contain three finite values')
    handoff = FiniteReleasePuffHandoff(
        status='transition_ready',source_duration_s=0.,transition_time_s=0.,
        transition_arc_length_m=0.,centre_position_m=(*cloud.center_xy_m,height/2),
        transverse_widths_m=(sigma_xy,height/math.sqrt(12)),
        longitudinal_half_length_m=config.width_cutoff_sigma*sigma_xy,
        bulk_velocity_m_s=tuple(float(v) for v in initial_velocity),total_mass_kg=mass,hydrogen_mass_kg=species,
        momentum_kg_m_s=tuple(float(v) for v in mass*initial_velocity),
        relative_energy_j=(mass*cp*(cloud.temperature_k-thermo.ambient_temperature)
                           + 0.5*mass*float(initial_velocity@initial_velocity)),
        hydrogen_flux_residual_kg_s=0.,hydrogen_inventory_residual_kg=0.,
        requires_transient_puff_continuation=True,
        qualification='uniform binary prism mapped to Gaussian integral moments; obstacles not resolved',
    )
    return integrate_finite_gaussian_puff(handoff,thermo,config)


def binary_puff_traces(puff: GaussianPuffResult, points_m: np.ndarray, *, state_count: int | None = None) -> np.ndarray:
    """Vectorized mol% observation, equivalent to the native point observer.

    Conserved species mass density converts directly to ideal-gas mole
    fraction. Reject values above pure species rather than silently clipping.
    Return shape (number of receptors, number of integrated states).
    """
    points = np.asarray(points_m,float)
    if points.ndim!=2 or points.shape[1]!=3 or np.any(~np.isfinite(points)) or np.any(points[:,2]<0):
        raise ValueError('points require finite (n,3) coordinates above ground')
    states=puff.states if state_count is None else puff.states[:state_count]
    result = np.empty((len(points),len(states)))
    cutoff = puff.config.width_cutoff_sigma
    for i,state in enumerate(states):
        centre = np.asarray(state.centre_position_m)
        velocity = np.asarray(state.bulk_velocity_m_s[:2])
        if np.linalg.norm(velocity)<=1e-12:
            velocity = np.asarray(puff.config.wind_velocity_m_s[:2])
        along = velocity/np.linalg.norm(velocity) if np.linalg.norm(velocity)>1e-12 else np.array([1.,0.])
        cross = np.array([-along[1],along[0]])
        offset = points[:,:2]-centre[:2]
        sx,sy,sz = np.array([state.longitudinal_half_length_m,state.lateral_half_width_m,
                            state.vertical_half_width_m])/cutoff
        horizontal = np.exp(-.5*((offset@along/sx)**2+(offset@cross/sy)**2))
        shape = horizontal*(np.exp(-.5*((points[:,2]-centre[2])/sz)**2)
                            +np.exp(-.5*((points[:,2]+centre[2])/sz)**2))
        density = state.hydrogen_mass_kg/((2*np.pi)**1.5*sx*sy*sz)*shape
        temperature = (puff.thermodynamics.ambient_temperature
            +(state.bulk_temperature_k-puff.thermodynamics.ambient_temperature)*np.minimum(shape,1.))
        mole = density*8.314462618*temperature/(puff.thermodynamics.closure.ambient_pressure_pa
                                               *puff.thermodynamics.fuel_molecular_weight)
        if np.any(mole>1+1e-12):
            raise RuntimeError('Gaussian source mapping exceeds pure-species concentration')
        result[:,i] = 100*mole
    return result


def axisymmetric_puff_inventory(puff: GaussianPuffResult, state_index: int,
                               x_edges_m, y_edges_m, z_edges_m):
    """Exact cell-integrated species handoff, without truncation renormalizing.

    Only horizontal axisymmetry is accepted: a yawed anisotropic puff needs
    bivariate rectangle integration. Returned exterior inventory is explicit;
    it must not be silently injected or discarded by the receiving solver.
    This map transfers species only, NOT density/momentum/enthalpy feedback.
    """
    state=puff.states[state_index]
    sx,sy,sz=np.array([state.longitudinal_half_length_m,state.lateral_half_width_m,
                       state.vertical_half_width_m])/puff.config.width_cutoff_sigma
    if not np.isclose(sx,sy,rtol=1e-12,atol=1e-12):
        raise ValueError('inventory mapping requires horizontally axisymmetric puff')
    edges=[np.asarray(a,float) for a in (x_edges_m,y_edges_m,z_edges_m)]
    if any(a.ndim!=1 or len(a)<2 or np.any(~np.isfinite(a)) or np.any(np.diff(a)<=0) for a in edges):
        raise ValueError('handoff edges must be strictly increasing and finite')
    if edges[2][0]<0:
        raise ValueError('ground-reflected handoff requires z>=0')
    cx,cy,cz=state.centre_position_m
    px=np.diff(ndtr((edges[0]-cx)/sx))
    py=np.diff(ndtr((edges[1]-cy)/sy))
    pz=np.diff(ndtr((edges[2]-cz)/sz))+np.diff(ndtr((edges[2]+cz)/sz))
    mass=state.hydrogen_mass_kg*pz[:,None,None]*py[None,:,None]*px[None,None,:]
    volume=np.diff(edges[2])[:,None,None]*np.diff(edges[1])[None,:,None]*np.diff(edges[0])[None,None,:]
    return mass/volume, float(state.hydrogen_mass_kg-mass.sum())
