"""Compact uniform-prism integral research model and conservative interface.

Unlike a Gaussian observation overlay, every observer and cell handoff uses
the SAME finite-support prism and the same bulk EOS. Initial concentration,
inventory, geometry and momentum therefore describe the same cloud.

The optional variational spreading model is an isothermal affine-motion
hypothesis, not a validated gravity-current head model. It retains vertical
inertia and starts with zero radial velocity; no fitted time delay is used.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.special import erfc

from .binary_integral_puff import BinaryPuffThermodynamics
from .finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud, _clip_polygon_axis


@dataclass(frozen=True)
class PrismState:
    elapsed_s: float
    centre_xy_m: tuple[float,float]
    apothem_m: float
    height_m: float
    total_mass_kg: float
    species_mass_kg: float
    bulk_velocity_xy_m_s: tuple[float,float]
    radial_velocity_m_s: float
    density_kg_m3: float
    temperature_k: float
    entrained_air_kg: float


@dataclass
class PrismIntegralResult:
    states: tuple[PrismState,...]
    cloud: RegularPolygonPrismInitialCloud
    thermodynamics: BinaryPuffThermodynamics
    diagnostics: dict
    front_profile: str = 'uniform'
    handoff_profile: str = 'uniform'
    grounded_head_fraction: float = 0.38
    edge_extension_fraction: float = 0.5
    edge_smoothing_fraction: float = 0.25
    longitudinal_ramp_fraction: float = 0.5
    crosswind_scale: float = 1.0
    crosswind_observation_tail_sigma: float = 0.0
    vertical_profile: str = 'uniform'
    vertical_profile_amplitude: float = 0.0
    wake_extension_fraction: float = 0.5
    vertical_tail_scale_m: float = 1.0
    vertical_tail_extension_m: float = 4.0
    leading_edge_scale: float = 1.0
    leading_edge_onset_fraction: float = 0.0
    crosswind_drift_ratio: float = 0.0
    vertical_observation_profile: str = 'top_hat'
    vertical_observation_tail_scale_m: float = 1.0
    vertical_observation_tail_extension_m: float = 0.0
    vertical_observation_tail_diffusivity_m2_s: float = 1.0
    vertical_observation_tail_until_s: float | None = None
    vertical_crosswind_spread_m_per_m: float = 0.0
    source_wind_scale: float = 1.0
    release_duration_s: float | None = None
    release_injection_velocity_xy_m_s: tuple[float,float] | None = None

    def vertices(self,index):
        state=self.states[index]
        base=self.cloud.vertices_xy_m()-np.array(self.cloud.center_xy_m)
        return base*(state.apothem_m/self.cloud.apothem_m)+state.centre_xy_m

    def traces(self,points_m):
        points=np.asarray(points_m,float)
        if points.ndim!=2 or points.shape[1]!=3 or np.any(~np.isfinite(points)) or np.any(points[:,2]<0):
            raise ValueError('receptors require finite (n,3) points above ground')
        angles=math.radians(self.cloud.rotation_deg)+np.arange(self.cloud.sides)*2*np.pi/self.cloud.sides
        normals=np.column_stack((np.cos(angles),np.sin(angles)))
        traces=np.empty((len(points),len(self.states)))
        for i,state in enumerate(self.states):
            offset=points[:,:2]-state.centre_xy_m
            inside=(np.max(offset@normals.T,axis=1)<=state.apothem_m+1e-12)&(points[:,2]<=state.height_m+1e-12)
            fraction=state.species_mass_kg/state.total_mass_kg
            if self.vertical_profile == 'cosine_dense' and state.height_m > 1.0e-30:
                # Smooth, mass-mean-preserving dense-gas stratification.  The
                # cosine has zero vertical mean; it changes only the local
                # observation profile until a matching vertical handoff exists.
                normalized_z=np.clip(points[:,2]/state.height_m,0.,1.)
                factor=1.0+self.vertical_profile_amplitude*np.cos(np.pi*normalized_z)
                local_fraction=np.clip(fraction*factor,0.,1.)
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(inside,concentration,0.)
                continue
            if self.front_profile == 'wake_tail' and i > 0:
                velocity=np.asarray(state.bulk_velocity_xy_m_s,dtype=float)
                speed=float(np.linalg.norm(velocity))
                along=velocity/speed if speed>1.0e-12 else np.array([1.,0.])
                cross=np.array([-along[1],along[0]])
                axial=offset@along
                lateral=offset@cross
                half_width=state.apothem_m
                tail_length=self.wake_extension_fraction*half_width
                tail_inside=(axial < -half_width)
                tail_inside &= axial >= -(half_width+tail_length)-1.0e-12
                tail_inside &= np.abs(lateral) <= half_width+1.0e-12
                tail_factor=np.clip(
                    1.0+(axial+half_width)/max(tail_length,1.0e-30),0.,1.)
                tail_concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(
                        fraction*tail_factor), dtype=float)
                core_concentration=100.*float(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(fraction))
                traces[:,i]=np.where(inside,core_concentration,
                    np.where(tail_inside,tail_concentration,0.))
                continue
            if self.front_profile == 'vertical_tail' and i > 0:
                horizontal_inside=np.max(offset@normals.T,axis=1)<=state.apothem_m+1.0e-12
                above=points[:,2]>state.height_m
                tail_inside=above & (
                    points[:,2] <= state.height_m+self.vertical_tail_extension_m+1.0e-12)
                factor=np.where(
                    above,
                    np.exp(-(points[:,2]-state.height_m)/max(self.vertical_tail_scale_m,1.0e-30)),
                    1.0)
                factor=np.where(tail_inside|~above,factor,0.0)
                local_fraction=np.clip(fraction*factor,0.,1.)
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(horizontal_inside,concentration,0.)
                continue
            if self.front_profile in {
                'wind_anisotropic_leading','longitudinal_ramp_anisotropic',
                'longitudinal_ramp_gaussian'} and i > 0:
                velocity=np.asarray(state.bulk_velocity_xy_m_s,dtype=float)
                speed=float(np.linalg.norm(velocity))
                along=velocity/speed if speed>1.0e-12 else np.array([1.,0.])
                cross=np.array([-along[1],along[0]])
                axial=offset@along
                lateral=offset@cross
                if self.front_profile == 'wind_anisotropic_leading':
                    onset=max(self.leading_edge_onset_fraction,0.)*state.apothem_m
                    scaled_axial=np.where(axial>onset,
                        onset+(axial-onset)/max(self.leading_edge_scale,1.0e-30),axial)
                else:
                    scaled_axial=axial
                if self.front_profile == 'longitudinal_ramp_gaussian':
                    core_inside=(np.max(offset@normals.T,axis=1)<=state.apothem_m+1.0e-12)&(
                        points[:,2]<=state.height_m+1.0e-12)
                    inside=core_inside.copy()
                    vertical_factor=np.ones(len(points))
                    vertical_tail=np.zeros(len(points),dtype=bool)
                    if self.vertical_observation_profile in {
                        'exponential_tail','diffusive_tail'}:
                        horizontal_inside=np.max(offset@normals.T,axis=1)<=state.apothem_m+1.0e-12
                        above=points[:,2]>state.height_m
                        tail_active=(self.vertical_observation_tail_until_s is None or
                            state.elapsed_s <= self.vertical_observation_tail_until_s+1.0e-12)
                        vertical_tail=(above & tail_active &
                            (points[:,2]<=state.height_m+
                            self.vertical_observation_tail_extension_m+1.0e-12))
                        if self.vertical_observation_profile == 'diffusive_tail':
                            penetration=np.sqrt(4.0*max(
                                self.vertical_observation_tail_diffusivity_m2_s,
                                1.0e-30)*max(state.elapsed_s,1.0e-12))
                            vertical_factor=np.where(
                                above,
                                erfc((points[:,2]-state.height_m)/max(
                                    penetration,1.0e-30)),1.0)
                        else:
                            vertical_factor=np.where(
                                above,
                                np.exp(-(points[:,2]-state.height_m)/max(
                                    self.vertical_observation_tail_scale_m,1.0e-30)),1.0)
                        vertical_factor=np.where(vertical_tail|~above,vertical_factor,0.)
                        inside=horizontal_inside & ((points[:,2]<=state.height_m+1.0e-12)|vertical_tail)
                else:
                    if self.front_profile == 'longitudinal_ramp_anisotropic':
                        lateral=lateral/max(self.crosswind_scale,1.0e-30)
                    transformed=np.outer(scaled_axial,along)+np.outer(lateral,cross)
                    support=np.max(transformed@normals.T,axis=1)/state.apothem_m
                    inside=(support<=1.+1.0e-12)&(points[:,2]<=state.height_m+1.0e-12)
                if self.front_profile in {'longitudinal_ramp_anisotropic','longitudinal_ramp_gaussian'}:
                    ramp=self.longitudinal_ramp_fraction*state.apothem_m
                    u=np.clip((axial-(state.apothem_m-ramp))/max(ramp,1e-30),0.,1.)
                    factor=.5*(1.+np.cos(np.pi*u))
                    if self.front_profile == 'longitudinal_ramp_gaussian':
                        sigma=max(self.crosswind_scale*state.apothem_m,1.0e-30)
                        # Optional 3-D observation diagnostic: unresolved
                        # overturning can broaden the lateral footprint above
                        # the collapsing dense layer.  The conserved prism
                        # remains unchanged and the default ratio is zero.
                        vertical_excess=np.maximum(points[:,2]-state.height_m,0.0)
                        sigma_eff=np.sqrt(sigma*sigma+
                            (self.vertical_crosswind_spread_m_per_m*vertical_excess)**2)
                        if self.crosswind_observation_tail_sigma > 0.0:
                            tail_extent=self.crosswind_observation_tail_sigma*sigma
                            lateral_tail=(np.abs(lateral)<=state.apothem_m+tail_extent+1.0e-12)
                            axial_tail=np.abs(axial)<=state.apothem_m+1.0e-12
                            inside=(core_inside|(
                                lateral_tail & axial_tail &
                                ((points[:,2]<=state.height_m+1.0e-12)|vertical_tail)))
                        drifted_lateral=lateral-self.crosswind_drift_ratio*np.maximum(axial,0.)
                        factor*=np.exp(-.5*(drifted_lateral/np.maximum(
                            sigma_eff,1.0e-30))**2)
                        factor*=vertical_factor
                    local_fraction=fraction*factor
                    concentration=100.*np.asarray(
                        self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                        dtype=float)
                else:
                    concentration=100.*float(
                        self.thermodynamics.closure.mole_fraction_from_mass_fraction(fraction))
                traces[:,i]=np.where(inside,concentration,0.)
                continue
            if self.front_profile == 'grounded_head' and i > 0:
                beta=self.grounded_head_fraction
                inside &= points[:,2] <= beta*state.height_m+1e-12
                local_mole_fraction=min(self.cloud.species_mole_fraction/beta,1.)
                concentration=100.*local_mole_fraction
                traces[:,i]=np.where(inside,concentration,0.)
            elif self.front_profile == 'entrained_edge' and state.entrained_air_kg > 0.:
                normalized=self._entrained_edge_factor(offset,normals,state)
                local_fraction=fraction*normalized
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(inside,concentration,0.)
            elif self.front_profile == 'diffuse_edge' and i > 0:
                extension=self.edge_extension_fraction
                support=np.max(offset@normals.T,axis=1)/state.apothem_m
                inside=(support<=1.+extension+1e-12)&(points[:,2]<=state.height_m+1e-12)
                factor=self._diffuse_edge_factor(support)
                # The normalized diagnostic tail can exceed the binary
                # mixture's admissible mass fraction in the core.  Cap it at
                # pure species rather than allowing an impossible state; this
                # branch is explicitly diagnostic and is not a conserved
                # handoff profile.
                local_fraction=np.clip(fraction*factor,0.,1.)
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(inside,concentration,0.)
            elif self.front_profile == 'smooth_edge' and i > 0:
                width=self.edge_smoothing_fraction
                support=np.max(offset@normals.T,axis=1)/state.apothem_m
                core=1.-width
                u=np.clip((support-core)/max(width,1e-30),0.,1.)
                factor=np.where(support<=core,1.,.5*(1.+np.cos(np.pi*u)))
                local_fraction=fraction*factor
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(inside,concentration,0.)
            elif self.front_profile == 'longitudinal_ramp' and i > 0:
                ramp=self.longitudinal_ramp_fraction*state.apothem_m
                velocity=np.asarray(state.bulk_velocity_xy_m_s,dtype=float)
                speed=float(np.linalg.norm(velocity))
                along=velocity/speed if speed>1e-12 else np.array([1.,0.])
                axial=offset@along
                u=np.clip((axial-(state.apothem_m-ramp))/max(ramp,1e-30),0.,1.)
                factor=.5*(1.+np.cos(np.pi*u))
                local_fraction=fraction*factor
                concentration=100.*np.asarray(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(local_fraction),
                    dtype=float)
                traces[:,i]=np.where(inside,concentration,0.)
            elif self.front_profile == 'anisotropic_support' and i > 0:
                velocity=np.asarray(state.bulk_velocity_xy_m_s,dtype=float)
                speed=float(np.linalg.norm(velocity))
                along=velocity/speed if speed>1e-12 else np.array([1.,0.])
                cross=np.array([-along[1],along[0]])
                axial=offset@along
                lateral=offset@cross/max(self.crosswind_scale,1e-30)
                transformed=np.outer(axial,along)+np.outer(lateral,cross)
                support=np.max(transformed@normals.T,axis=1)/state.apothem_m
                inside=(support<=1.+1e-12)&(points[:,2]<=state.height_m+1e-12)
                concentration=100.*float(
                    self.thermodynamics.closure.mole_fraction_from_mass_fraction(fraction))
                traces[:,i]=np.where(inside,concentration,0.)
            else:
                concentration=100*float(self.thermodynamics.closure.mole_fraction_from_mass_fraction(fraction))
                traces[:,i]=np.where(inside,concentration,0.)
        return traces

    def _entrained_edge_factor(self,offset,normals,state):
        """Mass-conserving radial H2 profile for the diagnostic edge model.

        Entrainment is represented as a homothetic outer transition shell.
        The shell factor is zero at the front and one in the core; its
        polygon-area mean is normalized analytically so the integrated species
        mass remains the state value. This is a profile hypothesis, not a
        resolved head-vortex solution.
        """
        if state.entrained_air_kg <= 0.:
            return np.ones(len(offset))
        volume=state.total_mass_kg/max(state.density_kg_m3,1e-30)
        air_volume=state.entrained_air_kg/max(self.thermodynamics.ambient_density,1e-30)
        volume_fraction=float(np.clip(air_volume/max(volume,1e-30),0.,1.))
        width_ratio=1.-math.sqrt(max(1.-volume_fraction,0.))
        if width_ratio<=1e-12:
            return np.ones(len(offset))
        normalized_area_mean=1.-width_ratio+2.*width_ratio**2/3.
        support=np.max(offset@normals.T,axis=1)/state.apothem_m
        core=1.-width_ratio
        taper=np.where(support<=core,1.,np.clip((1.-support)/width_ratio,0.,1.))
        return taper/max(normalized_area_mean,1e-30)

    def _diffuse_edge_factor(self,support):
        """Outer concentration taper used only as a front-shape diagnostic.

        The original polygon remains the conserved state geometry.  This
        profile extends the observed support by a declared fraction and
        normalizes a plateau-plus-linear-tail profile by its homothetic area
        mean.  It is intentionally not transferred to the dense layer until a
        matching expanded-volume handoff is implemented.
        """
        extension=self.edge_extension_fraction
        if extension <= 1e-12:
            return np.ones_like(support)
        radius=1.+extension
        mean=(radius**3/3.-radius+5./3.)/(radius**2)
        taper=np.where(support<=1.,1.,np.clip((radius-support)/extension,0.,1.))
        return taper/max(mean,1e-30)


def prism_from_gaussian_trajectory(puff,cloud):
    """Use old bulk/width dynamics, but replace Gaussian profile consistently.

    This isolates the geometry/interface change. The legacy width trajectory
    still has instantaneous spreading and is NOT the repaired inertia model.
    """
    first=puff.states[0]
    if any(not np.isclose(s.longitudinal_half_length_m,s.lateral_half_width_m,rtol=1e-12) for s in puff.states):
        raise ValueError('compact prism bridge requires horizontal isotropy')
    times=np.array([s.elapsed_s for s in puff.states])
    radii=cloud.apothem_m*np.array([s.longitudinal_half_length_m/first.longitudinal_half_length_m for s in puff.states])
    speeds=np.gradient(radii,times)
    states=[]
    for s,radius,speed in zip(puff.states,radii,speeds):
        area=cloud.horizontal_area_m2*(radius/cloud.apothem_m)**2
        depth=s.total_mass_kg/(s.bulk_density_kg_m3*area)
        states.append(PrismState(s.elapsed_s,s.centre_position_m[:2],float(radius),depth,
            s.total_mass_kg,s.hydrogen_mass_kg,s.bulk_velocity_m_s[:2],float(speed),
            s.bulk_density_kg_m3,s.bulk_temperature_k,s.cumulative_entrained_air_kg))
    return PrismIntegralResult(tuple(states),cloud,puff.thermodynamics,{
        'model':'compact prism on legacy integral width trajectory',
        'radial_inertia':False,'density_feedback':True,
        'maximum_relative_species_residual':puff.maximum_relative_hydrogen_residual,
        'temperature_model':'legacy puff enthalpy',
        'initial_radial_velocity_m_s':float(speeds[0])})


def _radial_inertia_factor(cloud):
    p=cloud.vertices_xy_m()-np.array(cloud.center_xy_m);q=np.roll(p,-1,axis=0)
    cross=p[:,0]*q[:,1]-q[:,0]*p[:,1]
    polar=float(np.sum(cross*(np.sum(p*p,axis=1)+np.sum(p*q,axis=1)+np.sum(q*q,axis=1)))/(12*cloud.horizontal_area_m2))
    return polar/cloud.apothem_m**2


def integrate_inertial_prism(cloud,closure,*,duration_s=10.,dt_s=.0125,
                             wind_xy_m_s=(3.9,0.),entrainment_shear=.1,
                             source_wind_scale=1.0,
                             entrainment_ambient=.1,drag_coefficient=.2,
                             vertical_inertia=True,
                             entrainment_top_shear=None,
                             entrainment_top_ambient=None,
                             entrainment_side_shear=None,
                             entrainment_side_ambient=None,
                             front_profile='uniform',
                             handoff_profile='uniform',
                             grounded_head_fraction=.38,
                             edge_extension_fraction=.5,
                             edge_smoothing_fraction=.25,
                             longitudinal_ramp_fraction=.5,
                             crosswind_scale=1.0,
                             crosswind_observation_tail_sigma=0.0,
                             radial_damping_coefficient=0.,
                             radial_damping_aspect_ratio=0.,
                             radial_damping_aspect_threshold=0.,
                             minimum_height_m=0.,
                             radial_froude_limit=None,
                             vertical_profile='uniform',
                             vertical_profile_amplitude=0.0,
                             dense_entrainment_model='constant',
                             wake_extension_fraction=0.5,
                             vertical_tail_scale_m=1.0,
                             vertical_tail_extension_m=4.0,
                             leading_edge_scale=1.0,
                             leading_edge_onset_fraction=0.0,
                             crosswind_drift_ratio=0.0,
                             vertical_observation_profile='top_hat',
                             vertical_observation_tail_scale_m=1.0,
                             vertical_observation_tail_extension_m=0.0,
                             vertical_observation_tail_diffusivity_m2_s=1.0,
                             vertical_observation_tail_until_s=None,
                             vertical_crosswind_spread_m_per_m=0.0,
                             release_duration_s=None,
                             release_injection_velocity_xy_m_s=None):
    """Affine isothermal potential-to-kinetic spreading, RK4 radial dynamics.

    For a prism a(t), H=V/(n tan(pi/n) a^2), velocity fields are radial
    u_r=(a_dot/a)r and w=-(2 a_dot/a)z in the fixed-volume gravity substep.
    K=.5 M [C_r+4 H^2/(3 a^2)] a_dot^2; U=(M-rho_air V)gH/2.
    The Euler-Lagrange equation retains the changing inertia metric. Ambient
    entrainment is split and its energy change is explicitly ledgered, not
    called a closed full thermodynamic energy budget.
    """
    # Legacy coefficients apply to both exposed surfaces. Optional
    # face-specific coefficients expose upper and lateral exchange separately
    # while preserving old calls exactly when all four are omitted.
    face_coefficients=(
        entrainment_shear if entrainment_top_shear is None else entrainment_top_shear,
        entrainment_ambient if entrainment_top_ambient is None else entrainment_top_ambient,
        entrainment_shear if entrainment_side_shear is None else entrainment_side_shear,
        entrainment_ambient if entrainment_side_ambient is None else entrainment_side_ambient,
    )
    controls=np.array([duration_s,dt_s,entrainment_shear,entrainment_ambient,
                       drag_coefficient,*face_coefficients],float)
    if np.any(~np.isfinite(controls)) or duration_s<=0 or dt_s<=0 or np.any(controls[2:]<0):
        raise ValueError('invalid prism integration controls')
    if not math.isfinite(source_wind_scale) or source_wind_scale <= 0.0:
        raise ValueError('source_wind_scale must be positive and finite')
    if front_profile not in {'uniform','entrained_edge','grounded_head','diffuse_edge','smooth_edge','longitudinal_ramp','longitudinal_ramp_anisotropic','longitudinal_ramp_gaussian','anisotropic_support','wind_anisotropic_leading','wake_tail','vertical_tail'}:
        raise ValueError('unsupported front_profile')
    if handoff_profile not in {'uniform','match_front_profile'}:
        raise ValueError('unsupported handoff_profile')
    if (not math.isfinite(grounded_head_fraction) or
            not 0.0 < grounded_head_fraction <= 1.0):
        raise ValueError('grounded_head_fraction must lie in (0, 1]')
    if (not math.isfinite(edge_extension_fraction) or edge_extension_fraction < 0.):
        raise ValueError('edge_extension_fraction must be finite and non-negative')
    if (not math.isfinite(edge_smoothing_fraction) or
            not 0.0 < edge_smoothing_fraction <= 1.0):
        raise ValueError('edge_smoothing_fraction must lie in (0, 1]')
    if (not math.isfinite(longitudinal_ramp_fraction) or
            not 0.0 < longitudinal_ramp_fraction <= 1.0):
        raise ValueError('longitudinal_ramp_fraction must lie in (0, 1]')
    if (not math.isfinite(crosswind_scale) or not 0.0 < crosswind_scale <= 1.0):
        raise ValueError('crosswind_scale must lie in (0, 1]')
    if (not math.isfinite(crosswind_observation_tail_sigma)
            or crosswind_observation_tail_sigma < 0.0):
        raise ValueError('crosswind_observation_tail_sigma must be finite and non-negative')
    if (not math.isfinite(leading_edge_scale) or not 0.0 < leading_edge_scale <= 1.0):
        raise ValueError('leading_edge_scale must lie in (0, 1]')
    if (not math.isfinite(leading_edge_onset_fraction) or
            not 0.0 <= leading_edge_onset_fraction <= 1.0):
        raise ValueError('leading_edge_onset_fraction must lie in [0, 1]')
    if not math.isfinite(crosswind_drift_ratio) or crosswind_drift_ratio < 0.0:
        raise ValueError('crosswind_drift_ratio must be finite and non-negative')
    if vertical_observation_profile not in {
            'top_hat','exponential_tail','diffusive_tail'}:
        raise ValueError('unsupported vertical observation profile')
    if (not math.isfinite(vertical_observation_tail_scale_m)
            or vertical_observation_tail_scale_m<=0
            or not math.isfinite(vertical_observation_tail_extension_m)
            or vertical_observation_tail_extension_m<0):
        raise ValueError('invalid vertical observation tail controls')
    if (not math.isfinite(vertical_observation_tail_diffusivity_m2_s)
            or vertical_observation_tail_diffusivity_m2_s<=0):
        raise ValueError('vertical_observation_tail_diffusivity_m2_s must be positive and finite')
    if (vertical_observation_tail_until_s is not None and
            (not math.isfinite(vertical_observation_tail_until_s)
             or vertical_observation_tail_until_s < 0.0)):
        raise ValueError('vertical_observation_tail_until_s must be non-negative and finite')
    if (not math.isfinite(vertical_crosswind_spread_m_per_m)
            or vertical_crosswind_spread_m_per_m<0):
        raise ValueError('vertical_crosswind_spread_m_per_m must be finite and non-negative')
    if release_duration_s is not None and (
            not math.isfinite(release_duration_s) or release_duration_s<=0):
        raise ValueError('release_duration_s must be positive and finite when supplied')
    if release_injection_velocity_xy_m_s is not None:
        release_velocity=np.asarray(release_injection_velocity_xy_m_s,float)
        if release_velocity.shape!=(2,) or np.any(~np.isfinite(release_velocity)):
            raise ValueError('release_injection_velocity_xy_m_s must be a finite 2-vector')
    else:
        release_velocity=None
    if (not math.isfinite(radial_damping_coefficient) or radial_damping_coefficient < 0.):
        raise ValueError('radial_damping_coefficient must be finite and non-negative')
    if (not math.isfinite(radial_damping_aspect_ratio) or
            radial_damping_aspect_ratio < 0.):
        raise ValueError('radial_damping_aspect_ratio must be finite and non-negative')
    if (not math.isfinite(radial_damping_aspect_threshold) or
            radial_damping_aspect_threshold < 0.):
        raise ValueError('radial_damping_aspect_threshold must be finite and non-negative')
    if (not math.isfinite(minimum_height_m) or minimum_height_m < 0.):
        raise ValueError('minimum_height_m must be finite and non-negative')
    if minimum_height_m > cloud.z_max_m-cloud.z_min_m:
        raise ValueError('minimum_height_m cannot exceed initial cloud height')
    if radial_froude_limit is not None and (
            not math.isfinite(radial_froude_limit) or radial_froude_limit <= 0.):
        raise ValueError('radial_froude_limit must be positive when supplied')
    if vertical_profile not in {'uniform','cosine_dense'}:
        raise ValueError('unsupported vertical_profile')
    if (not math.isfinite(vertical_profile_amplitude)
            or not 0.0 <= vertical_profile_amplitude < 1.0):
        raise ValueError('vertical_profile_amplitude must lie in [0, 1)')
    if not math.isfinite(wake_extension_fraction) or wake_extension_fraction < 0.0:
        raise ValueError('wake_extension_fraction must be finite and non-negative')
    if not math.isfinite(vertical_tail_scale_m) or vertical_tail_scale_m <= 0.0:
        raise ValueError('vertical_tail_scale_m must be positive and finite')
    if not math.isfinite(vertical_tail_extension_m) or vertical_tail_extension_m < 0.0:
        raise ValueError('vertical_tail_extension_m must be finite and non-negative')
    if dense_entrainment_model not in {
        'constant','parker_bulk_suppression',
        'parker_top_suppression','parker_side_suppression',
    }:
        raise ValueError('unsupported dense_entrainment_model')
    if cloud.z_min_m!=0 or closure.species_molar_mass_kg_mol<=closure.carrier_molar_mass_kg_mol:
        raise ValueError('inertial prism requires grounded dense species')
    if not math.isclose(cloud.temperature_k,closure.ambient_temperature_k):
        raise ValueError('isothermal prism currently requires ambient source temperature')
    thermo=BinaryPuffThermodynamics(closure)
    ambient_wind=np.asarray(wind_xy_m_s,float)
    if ambient_wind.shape!=(2,) or np.any(~np.isfinite(ambient_wind)):
        raise ValueError('finite wind vector required')
    # A grounded dense current can sample a lower wind than the reference
    # meteorological vector.  Keep this explicit and diagnostic: the receiving
    # layer still uses the measured ambient wind, while source transport and
    # default release momentum use the declared source-scale vector.
    wind=ambient_wind*float(source_wind_scale)
    if release_velocity is None:
        release_velocity=wind.copy()
    y=float(closure.mass_fraction_from_mole_fraction(cloud.species_mole_fraction))
    if vertical_profile == 'cosine_dense' and y*(1.0+vertical_profile_amplitude) > 1.0+1.0e-12:
        raise ValueError('vertical cosine profile exceeds admissible mass fraction')
    rho=thermo._density_from_temperature(cloud.temperature_k,y)
    full_mass=cloud.exact_volume_m3*rho;full_species=full_mass*y
    initial_mass=full_mass;initial_species=full_species
    if release_duration_s is None:
        mass=full_mass;species=full_species
    else:
        # Stage the marked source from a tiny positive seed so the primitive
        # variables remain defined.  The seed carries the same composition;
        # subsequent source increments are conservative and wind-momentum
        # consistent.  This is an opt-in source-history diagnostic, not a
        # fitted default.
        seed_fraction=1.0e-12
        mass=full_mass*seed_fraction;species=mass*y
    radius=cloud.apothem_m;radial=0.;centre=np.array(cloud.center_xy_m,float);momentum=np.zeros(2)
    shape=cloud.sides*np.tan(np.pi/cloud.sides);c_r=_radial_inertia_factor(cloud)
    vertical=1. if vertical_inertia else 0.
    def geometry(a,m):
        density=thermo._density_from_temperature(cloud.temperature_k,species/m)
        volume=m/density;height=volume/(shape*a*a)
        excess=m-thermo.ambient_density*volume
        metric=c_r+vertical*(4/3)*(height/a)**2
        return height,density,excess,metric
    def energy(a,v,m):
        height,_,excess,metric=geometry(a,m)
        return .5*m*metric*v*v+excess*9.80665*height/2
    def rhs(z,m):
        a,v=z;height,_,excess,metric=geometry(a,m)
        acceleration=((excess/m)*9.80665*height/a+vertical*4*height**2/a**3*v*v)/metric
        aspect_excess=max(a/max(height,1.0e-30)-1.0-
                          radial_damping_aspect_threshold,0.0)
        effective_damping=radial_damping_coefficient*(
            1.0+radial_damping_aspect_ratio*aspect_excess)
        acceleration-=effective_damping*v*abs(v)
        return np.array([v,acceleration])
    time=0.;air=0.;top_air=0.;side_air=0.;states=[];energy0=energy(radius,radial,mass);mixwork=0.;release_work=0.;radial_dissipation=0.;maxerr=0.;maxspecies=0.;released_source_mass=mass
    def save():
        height,density,_,_=geometry(radius,mass)
        states.append(PrismState(time,tuple(centre),radius,height,mass,species,
            tuple(momentum/mass),radial,density,cloud.temperature_k,air))
    save()
    while time<duration_s-1e-12:
        dt=min(dt_s,duration_s-time)
        if release_duration_s is not None and time < release_duration_s-1.0e-12:
            before_release=energy(radius,radial,mass)
            target=full_mass*min((time+dt)/release_duration_s,1.0)
            dm_release=max(target-mass,0.0)
            mass+=dm_release;species+=dm_release*y;released_source_mass+=dm_release
            momentum+=dm_release*release_velocity
            release_work+=energy(radius,radial,mass)-before_release
        height,density,_,metric=geometry(radius,mass)
        slip=wind-momentum/mass
        area=shape*radius**2
        perimeter=2*cloud.sides*radius*np.tan(np.pi/cloud.sides)
        # Upper and lateral interfaces can now be ablated independently. The
        # slip/ambient pieces are retained on each face so the historical
        # closure remains recoverable when all face-specific values are None.
        top_velocity=(face_coefficients[0]*np.linalg.norm(slip)
                      +face_coefficients[1]*np.linalg.norm(wind))
        side_velocity=(face_coefficients[2]*np.linalg.norm(slip)
                       +face_coefficients[3]*np.linalg.norm(wind))
        entrainment_suppression=1.0
        top_entrainment_suppression=1.0
        side_entrainment_suppression=1.0
        if dense_entrainment_model != 'constant':
            # Parker et al. (1987) bulk-Richardson suppression, shared with
            # the finite-puff closure.  Apply it to both exposed surfaces;
            # this is a selectable physics closure, not an observation fit.
            g_prime=9.80665*max(density/thermo.ambient_density-1.0,0.0)
            relative_speed=max(float(np.linalg.norm(slip)),1.0e-15)
            ri=g_prime*height/(relative_speed**2)
            entrainment_suppression=(1.0+718.0*ri**2.4)**-0.5
            if dense_entrainment_model in {
                'parker_bulk_suppression','parker_top_suppression',
            }:
                top_entrainment_suppression=entrainment_suppression
                top_velocity*=top_entrainment_suppression
            if dense_entrainment_model in {
                'parker_bulk_suppression','parker_side_suppression',
            }:
                side_entrainment_suppression=entrainment_suppression
                side_velocity*=side_entrainment_suppression
        dm_top=thermo.ambient_density*area*top_velocity*dt
        dm_side=thermo.ambient_density*perimeter*height*side_velocity*dt
        dm=dm_top+dm_side
        top_air+=dm_top;side_air+=dm_side
        before=energy(radius,radial,mass)
        radial_p=mass*metric*radial
        old_velocity=momentum/mass
        mass+=dm;air+=dm;momentum+=dm*wind
        height,density,_,newmetric=geometry(radius,mass)
        radial=radial_p/(mass*newmetric)
        mixwork+=energy(radius,radial,mass)-before
        # Ambient drag changes only cloud translation here; it is not an
        # unrecorded sink in the radial gravity energy balance.
        projected=2*radius*height
        factor=1/(1+.5*drag_coefficient*thermo.ambient_density*projected*np.linalg.norm(slip)*dt/mass)
        momentum=mass*wind+factor*(momentum-mass*wind)
        release_active=(release_duration_s is not None and
                        time+dt < release_duration_s-1.0e-12)
        z=np.array([radius,radial])
        if release_active:
            # The source is still being formed; do not spend the unfinished
            # cloud's gravity potential as radial spread before its declared
            # release history has supplied the bulk volume.
            k1=k2=k3=k4=np.array([0.0,0.0])
        else:
            k1=rhs(z,mass);k2=rhs(z+.5*dt*k1,mass)
            k3=rhs(z+.5*dt*k2,mass);k4=rhs(z+dt*k3,mass)
        if radial_damping_coefficient > 0. and not release_active:
            def damping_rate(stage):
                stage_a,stage_v=stage
                stage_height,_,_,stage_metric=geometry(stage_a,mass)
                stage_aspect_excess=max(stage_a/max(stage_height,1.0e-30)-1.0-
                                        radial_damping_aspect_threshold,0.0)
                effective_damping=radial_damping_coefficient*(
                    1.0+radial_damping_aspect_ratio*stage_aspect_excess)
                return mass*stage_metric*effective_damping*abs(stage_v)**3
            radial_dissipation+=dt*(damping_rate(z)+2.*damping_rate(z+.5*dt*k1)
                +2.*damping_rate(z+.5*dt*k2)+damping_rate(z+dt*k3))/6.
        radius,radial=z+dt*(k1+2*k2+2*k3+k4)/6
        if minimum_height_m > 0.:
            limited_height,limited_density,_,_=geometry(radius,mass)
            if limited_height < minimum_height_m:
                before_limit=energy(radius,radial,mass)
                volume=mass/max(limited_density,1e-30)
                radius=math.sqrt(volume/(shape*minimum_height_m))
                radial=0.
                radial_dissipation+=max(before_limit-energy(radius,radial,mass),0.)
        centre+=.5*(old_velocity+momentum/mass)*dt
        time+=dt
        if not np.all(np.isfinite([radius,radial,mass])) or radius<=0 or radial<0:
            raise RuntimeError('prism positivity/finite gate failed')
        if radial_froude_limit is not None:
            height,_,excess,_=geometry(radius,mass)
            reduced_gravity=max(excess/mass,0.)*9.80665
            cap=radial_froude_limit*math.sqrt(max(reduced_gravity*height,0.))
            if radial>cap:
                before_cap=energy(radius,radial,mass)
                radial=cap
                radial_dissipation+=max(before_cap-energy(radius,radial,mass),0.)
        maxerr=max(maxerr,abs(energy(radius,radial,mass)-energy0-mixwork-release_work+
            radial_dissipation)/max(abs(energy0)+abs(release_work),1e-30))
        save();maxspecies=max(maxspecies,abs(states[-1].species_mass_kg-species))
    return PrismIntegralResult(tuple(states),cloud,thermo,{
        'model':'isothermal affine prism with radial and vertical inertia',
        'radial_inertia':True,'vertical_inertia':bool(vertical_inertia),
        'initial_radial_velocity_m_s':0.,'initial_available_potential_energy_j':energy0,
        'release_duration_s':release_duration_s,
        'source_release_staged':release_duration_s is not None,
        'release_target_mass_kg':full_mass,
        'release_target_species_mass_kg':full_species,
        'released_source_mass_kg':released_source_mass,
        'release_energy_input_j':release_work,
        'release_injection_velocity_xy_m_s':None if release_injection_velocity_xy_m_s is None
            else tuple(float(v) for v in release_velocity),
        'source_wind_scale':float(source_wind_scale),
        'source_wind_xy_m_s':tuple(float(v) for v in wind),
        'ambient_wind_xy_m_s':tuple(float(v) for v in ambient_wind),
        'entrainment_radial_energy_exchange_j':mixwork,
        'radial_dissipation_j':radial_dissipation,
        'radial_damping_coefficient':radial_damping_coefficient,
        'radial_damping_aspect_ratio':radial_damping_aspect_ratio,
        'radial_damping_aspect_threshold':radial_damping_aspect_threshold,
        'minimum_height_m':minimum_height_m,
        'radial_froude_limit':radial_froude_limit,
        'maximum_relative_radial_energy_ledger_residual':maxerr,
        'maximum_relative_species_residual':maxspecies/max(species,1e-30),
        'mixture_mass_residual_kg':mass-initial_mass-air,
        'coefficients':{'entrainment_shear':entrainment_shear,'entrainment_ambient':entrainment_ambient,
            'entrainment_top_shear':face_coefficients[0],
            'entrainment_top_ambient':face_coefficients[1],
            'entrainment_side_shear':face_coefficients[2],
            'entrainment_side_ambient':face_coefficients[3],
            'drag':drag_coefficient},
        'entrained_air_by_surface_kg':{'top':top_air,'side':side_air},
        'face_specific_entrainment':any(value is not None for value in (
            entrainment_top_shear,entrainment_top_ambient,
            entrainment_side_shear,entrainment_side_ambient)),
        'temperature_model':'isothermal ambient; full thermal energy not advanced',
        'operating_envelope':'pre-obstacle dense prism, no head vortices or 3D pressure field',
        'front_profile':front_profile,
        'handoff_profile':handoff_profile,
        'vertical_profile':vertical_profile,
        'vertical_profile_amplitude':vertical_profile_amplitude,
        'wake_extension_fraction':wake_extension_fraction,
        'vertical_tail_scale_m':vertical_tail_scale_m,
        'vertical_tail_extension_m':vertical_tail_extension_m,
        'dense_entrainment_model':dense_entrainment_model,
        'last_entrainment_suppression':entrainment_suppression,
        'last_top_entrainment_suppression':top_entrainment_suppression,
        'last_side_entrainment_suppression':side_entrainment_suppression,
        'grounded_head_fraction':grounded_head_fraction,
        'edge_extension_fraction':edge_extension_fraction,
        'edge_smoothing_fraction':edge_smoothing_fraction,
        'longitudinal_ramp_fraction':longitudinal_ramp_fraction,
        'crosswind_scale':crosswind_scale,
        'crosswind_observation_tail_sigma':crosswind_observation_tail_sigma,
        'leading_edge_scale':leading_edge_scale,
        'leading_edge_onset_fraction':leading_edge_onset_fraction,
        'crosswind_drift_ratio':crosswind_drift_ratio,
        'vertical_observation_profile':vertical_observation_profile,
        'vertical_observation_tail_scale_m':vertical_observation_tail_scale_m,
        'vertical_observation_tail_extension_m':vertical_observation_tail_extension_m,
        'vertical_observation_tail_diffusivity_m2_s':vertical_observation_tail_diffusivity_m2_s,
        'vertical_observation_tail_until_s':vertical_observation_tail_until_s,
        'vertical_crosswind_spread_m_per_m':vertical_crosswind_spread_m_per_m,
        'no_observed_arrivals_or_concentrations_used_for_coefficients':True},
        front_profile=front_profile,handoff_profile=handoff_profile,
        grounded_head_fraction=grounded_head_fraction,
        edge_extension_fraction=edge_extension_fraction,
        edge_smoothing_fraction=edge_smoothing_fraction,
        longitudinal_ramp_fraction=longitudinal_ramp_fraction,
        crosswind_scale=crosswind_scale,
        crosswind_observation_tail_sigma=crosswind_observation_tail_sigma,
        vertical_profile=vertical_profile,
        vertical_profile_amplitude=vertical_profile_amplitude,
        wake_extension_fraction=wake_extension_fraction,
        vertical_tail_scale_m=vertical_tail_scale_m,
        vertical_tail_extension_m=vertical_tail_extension_m,
        leading_edge_scale=leading_edge_scale,
        leading_edge_onset_fraction=leading_edge_onset_fraction,
        crosswind_drift_ratio=crosswind_drift_ratio,
        vertical_observation_profile=vertical_observation_profile,
        vertical_observation_tail_scale_m=vertical_observation_tail_scale_m,
        vertical_observation_tail_extension_m=vertical_observation_tail_extension_m,
        vertical_observation_tail_diffusivity_m2_s=vertical_observation_tail_diffusivity_m2_s,
        vertical_observation_tail_until_s=vertical_observation_tail_until_s,
        vertical_crosswind_spread_m_per_m=vertical_crosswind_spread_m_per_m,
        source_wind_scale=source_wind_scale,
        release_duration_s=release_duration_s,
        release_injection_velocity_xy_m_s=None if release_injection_velocity_xy_m_s is None
            else tuple(float(v) for v in release_velocity))


def prism_layer_inventory(result,index,x_edges,y_edges,*,radial_momentum=True):
    """Exact polygon/cell mass AND affine momentum integrals per unit area.

    Exterior inventory is explicit; no tail renormalization or cell-centre
    injection. z momentum/enthalpy are NOT consumed by the hydrostatic layer.
    """
    xe,ye=np.asarray(x_edges,float),np.asarray(y_edges,float)
    if any(a.ndim!=1 or len(a)<2 or np.any(~np.isfinite(a)) or np.any(np.diff(a)<=0) for a in (xe,ye)):
        raise ValueError('finite strictly increasing face edges required')
    state=result.states[index]
    offset=np.array(state.centre_xy_m)
    # Work in source-relative coordinates to avoid cancellation at site coords.
    vertices=result.vertices(index)-offset
    nx,ny=len(xe)-1,len(ye)-1
    q=np.zeros((4,ny,nx));covered=0.;profile_species_integral=0.
    profile_enabled=result.front_profile=='entrained_edge' and state.entrained_air_kg>0.
    handoff_front_enabled=(result.handoff_profile=='match_front_profile' and
        result.front_profile in {'longitudinal_ramp','longitudinal_ramp_anisotropic',
                                 'longitudinal_ramp_gaussian'})
    profile_width_ratio=0.;profile_norm=1.
    if profile_enabled:
        volume=state.total_mass_kg/max(state.density_kg_m3,1e-30)
        air_volume=state.entrained_air_kg/max(result.thermodynamics.ambient_density,1e-30)
        volume_fraction=float(np.clip(air_volume/max(volume,1e-30),0.,1.))
        profile_width_ratio=1.-math.sqrt(max(1.-volume_fraction,0.))
        profile_norm=1.-profile_width_ratio+2.*profile_width_ratio**2/3.

    def inside_convex(points,polygon):
        p=np.asarray(polygon,float);nxt=np.roll(p,-1,axis=0)
        cross=(nxt[:,0]-p[:,0])[:,None]*(points[:,1]-p[:,1,None])-(nxt[:,1]-p[:,1])[:,None]*(points[:,0]-p[:,0,None])
        return np.all(cross>=-1e-10,axis=0) | np.all(cross<=1e-10,axis=0)

    def mean_profile_factor(polygon,cell_x,cell_y):
        if not profile_enabled and not handoff_front_enabled:
            return 1.
        # Tensor Gauss sampling is used only for the optional front profile;
        # q[0] continues to use the exact clipped area. The ledger records the
        # quadrature residual and full-domain species is normalized below.
        nodes=np.array([-0.8611363116,-0.3399810436,0.3399810436,0.8611363116])
        weights=np.array([0.3478548451,0.6521451549,0.6521451549,0.3478548451])
        xs=cell_x[0]-offset[0]+(nodes+1.)*(cell_x[1]-cell_x[0])/2.
        ys=cell_y[0]-offset[1]+(nodes+1.)*(cell_y[1]-cell_y[0])/2.
        samples=np.array([(x,y) for y in ys for x in xs])
        mask=inside_convex(samples,np.asarray(polygon))
        if not np.any(mask):
            return 0.
        normals=np.column_stack((
            np.cos(math.radians(result.cloud.rotation_deg)+np.arange(result.cloud.sides)*2*np.pi/result.cloud.sides),
            np.sin(math.radians(result.cloud.rotation_deg)+np.arange(result.cloud.sides)*2*np.pi/result.cloud.sides)))
        if profile_enabled:
            support=np.max(samples@normals.T,axis=1)/state.apothem_m
            core=1.-profile_width_ratio
            taper=np.where(support<=core,1.,np.clip((1.-support)/profile_width_ratio,0.,1.))/max(profile_norm,1e-30)
        else:
            velocity=np.asarray(state.bulk_velocity_xy_m_s,dtype=float)
            speed=float(np.linalg.norm(velocity))
            along=velocity/speed if speed>1e-12 else np.array([1.,0.])
            cross=np.array([-along[1],along[0]])
            axial=samples@along
            lateral=samples@cross
            ramp=max(result.longitudinal_ramp_fraction*state.apothem_m,1e-30)
            u=np.clip((axial-(state.apothem_m-ramp))/ramp,0.,1.)
            taper=.5*(1.+np.cos(np.pi*u))
            if result.front_profile=='longitudinal_ramp_anisotropic':
                lateral=lateral/max(result.crosswind_scale,1e-30)
            if result.front_profile=='longitudinal_ramp_gaussian':
                sigma=max(result.crosswind_scale*state.apothem_m,1e-30)
                taper*=np.exp(-.5*(lateral/sigma)**2)
        ww=np.array([wy*wx for wy in weights for wx in weights])
        return float(np.sum(ww*mask*taper)/np.sum(ww*mask))
    candidate_x=np.flatnonzero((xe[:-1]<state.centre_xy_m[0]+vertices[:,0].max())&(xe[1:]>state.centre_xy_m[0]+vertices[:,0].min()))
    candidate_y=np.flatnonzero((ye[:-1]<state.centre_xy_m[1]+vertices[:,1].max())&(ye[1:]>state.centre_xy_m[1]+vertices[:,1].min()))
    for j in candidate_y:
        for i in candidate_x:
            polygon=[tuple(row) for row in vertices]
            for axis,bound,keep in ((0,xe[i]-offset[0],True),(0,xe[i+1]-offset[0],False),
                                    (1,ye[j]-offset[1],True),(1,ye[j+1]-offset[1],False)):
                polygon=_clip_polygon_axis(polygon,axis=axis,bound=bound,keep_greater=keep)
            if len(polygon)<3:continue
            p=np.array(polygon);following=np.roll(p,-1,axis=0)
            cross=p[:,0]*following[:,1]-following[:,0]*p[:,1]
            area=float(cross.sum()/2)
            if area<0:raise RuntimeError('polygon orientation reversed')
            first=np.sum((p+following)*cross[:,None],axis=0)/6
            volume=area*state.height_m;cellarea=(xe[i+1]-xe[i])*(ye[j+1]-ye[j])
            mixture=state.density_kg_m3*volume
            profile_factor=mean_profile_factor(p, (xe[i],xe[i+1]), (ye[j],ye[j+1]))
            if handoff_front_enabled:
                # Shape the complete marked mixture, then normalize the whole
                # column below. This preserves the local species fraction and
                # avoids an impossible q_species>q_mixture after normalization.
                mixture*=profile_factor
                species=mixture*state.species_mass_kg/state.total_mass_kg
            else:
                species=mixture*state.species_mass_kg/state.total_mass_kg
            if profile_enabled:
                species=state.density_kg_m3*state.height_m*area*state.species_mass_kg/(
                    state.total_mass_kg*max(1e-30,1.))
                species*=profile_factor
            momentum=mixture*np.array(state.bulk_velocity_xy_m_s)
            if radial_momentum:
                momentum+=state.density_kg_m3*state.height_m*state.radial_velocity_m_s/state.apothem_m*first
            q[:,j,i]=np.array([mixture,species,*momentum])/cellarea
            covered+=volume
            profile_species_integral+=species
    cellareas=np.diff(ye)[:,None]*np.diff(xe)[None,:]
    mapped=np.sum(q*cellareas,axis=(1,2))
    total=np.array([state.total_mass_kg,state.species_mass_kg,
                    *(state.total_mass_kg*np.array(state.bulk_velocity_xy_m_s))])
    grid_kinetic=float(np.sum(np.divide(q[2]**2+q[3]**2,2*q[0],
        out=np.zeros_like(q[0]),where=q[0]>0)*cellareas))
    # Exact radial horizontal KE is available for a full-domain mapping. The
    # receiving layer keeps only cell-average velocity; expose the lost
    # subcell KE and the vertical affine KE, rather than calling this a full
    # energy-conserving interface merely because mass and momentum balance.
    full_coverage=abs(total[0]-mapped[0])<=1e-10*max(total[0],1.)
    geometric_coverage=abs(covered-state.total_mass_kg/max(state.density_kg_m3,1e-30)) <= 1e-10*max(
        state.total_mass_kg/max(state.density_kg_m3,1e-30),1.)
    profile_scale=1.
    if handoff_front_enabled and geometric_coverage and mapped[0]>0.:
        profile_scale=state.total_mass_kg/mapped[0]
        q*=profile_scale
        mapped=np.sum(q*cellareas,axis=(1,2))
    elif profile_enabled and full_coverage and mapped[1]>0.:
        profile_scale=state.species_mass_kg/mapped[1]
        q[1]*=profile_scale
        if np.any(q[1]>q[0]+1e-10):
            raise RuntimeError('handoff profile exceeds admissible species mass')
        mapped[1]=float(np.sum(q[1]*cellareas))
    full_horizontal_ke=.5*state.total_mass_kg*(
        sum(v*v for v in state.bulk_velocity_xy_m_s)
        +(_radial_inertia_factor(result.cloud)*state.radial_velocity_m_s**2 if radial_momentum else 0.))
    vertical_spreading_ke=state.total_mass_kg*2/3*(state.height_m/state.apothem_m)**2*state.radial_velocity_m_s**2
    return q,{'mapped_conserved':mapped.tolist(),'exterior_conserved':(total-mapped).tolist(),
        'covered_volume_m3':covered,'radial_momentum_transferred':bool(radial_momentum),
        'vertical_momentum_transferred':False,'enthalpy_transferred':False,
        'cell_average_horizontal_kinetic_energy_j':grid_kinetic,
        'subcell_horizontal_kinetic_energy_not_transferred_j':full_horizontal_ke-grid_kinetic if full_coverage else None,
        'full_source_vertical_affine_spreading_kinetic_energy_not_transferred_j':vertical_spreading_ke,
        'full_energy_conserving_interface':False,
        'front_profile':result.front_profile,
        'handoff_profile':result.handoff_profile,
        'grounded_head_fraction':result.grounded_head_fraction,
        'vertical_profile':result.vertical_profile,
        'vertical_profile_amplitude':result.vertical_profile_amplitude,
        'vertical_profile_transferred':result.vertical_profile=='uniform',
        'front_profile_species_normalization':profile_scale,
        'front_profile_species_quadrature_integral_kg':profile_species_integral if (profile_enabled or handoff_front_enabled) else None,
        'front_profile_species_residual_kg':(mapped[1]-state.species_mass_kg) if (profile_enabled or handoff_front_enabled) else None}
