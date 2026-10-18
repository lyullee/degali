"""Experimental hydrostatic dense layer on a LOCAL horizontal FV patch.

Conserved variables: marked mixture mass, species mass, x/y momentum per area.
Height is recovered from binary EOS, not fixed to a vertical grid cell. A
hydrostatic density-current pressure drives gravity spreading. This is an
isothermal single-layer research closure, NOT a 3-D CFD/wake or light-gas model.
Finite-height fences reflect the sub-fence fraction and permit overtopping.
"""
from dataclasses import dataclass
import math
import time
from typing import Callable
import numpy as np
from scipy.special import erfc,ndtr,erf
from .binary_integral_puff import BinaryPuffThermodynamics
from .conservative_column_observation import reconstruct_column


@dataclass
class DenseLayerResult:
    time_s: np.ndarray
    species_traces_vol_percent: np.ndarray
    boundary_species_rate_kg_s: np.ndarray
    boundary_current_speed_m_s: np.ndarray
    boundary_height_m: np.ndarray
    boundary_positions_m: np.ndarray
    mixture_mass_kg: np.ndarray
    species_mass_kg: np.ndarray
    final_conserved: np.ndarray
    diagnostics: dict
    observation_component_traces_vol_percent: dict[str,np.ndarray] | None = None


def _primitive(q,thermo,gravity_scale=1.):
    m,s,px,py=q
    fraction=np.divide(s,m,out=np.zeros_like(m),where=m>1e-30)
    rho=thermo._density_from_temperature(thermo.ambient_temperature,fraction)
    height=m/rho
    vx=np.divide(px,m,out=np.zeros_like(m),where=m>1e-30)
    vy=np.divide(py,m,out=np.zeros_like(m),where=m>1e-30)
    reduced=gravity_scale*9.80665*np.maximum(1-thermo.ambient_density/rho,0.)
    pressure=.5*reduced*m*height
    return height,vx,vy,np.sqrt(reduced*height),pressure,fraction


def _flux_pair(left,right,axis,thermo,wall_height,gravity_scale=1.,wall_overtop_exponent=0.0,
               openness_override=None,wall_blockage_fraction=None,
               return_shared=False):
    hl,ul,vl,cl,pl,_=_primitive(left,thermo,gravity_scale)
    hr,ur,vr,cr,pr,_=_primitive(right,thermo,gravity_scale)
    nl,nr=(ul,ur) if axis==0 else (vl,vr)
    fl=left*nl; fr=right*nr
    fl[2+axis]+=pl; fr[2+axis]+=pr
    speed=np.maximum(np.abs(nl)+cl,np.abs(nr)+cr)
    flux=.5*(fl+fr)-.5*speed*(right-left)
    # Species/mixture use one shared flux, so arbitrary walls cannot create
    # species. Momentum reaction belongs to the wall, not to the other fluid.
    raw_openness=np.divide(np.maximum(np.maximum(hl,hr)-wall_height,0.),np.maximum(hl,hr),
                           out=np.ones_like(hl),where=np.maximum(hl,hr)>1e-30)
    # A single hydrostatic column cannot resolve the vertical mass profile at
    # a finite-height wall.  ``wall_overtop_exponent`` applies the declared
    # normalized profile to the above-wall fraction: exponent 0 is the
    # historical uniform-column closure, while positive values represent a
    # ground-loaded dense layer and reduce overtopping without changing total
    # mass or wall impulse accounting.
    if openness_override is None:
        # If c(z) is proportional to (1-z/H)**p, its *integrated* mass above
        # a wall is (1-wall/H)**(p+1).  The previous p exponent made the
        # uniform p=0 column completely open as soon as H exceeded the wall
        # and over-predicted every positive-p overtopping flux.
        openness=np.power(
            np.clip(raw_openness,0.0,1.0),wall_overtop_exponent+1.0)
        openness=np.where(raw_openness<=0.0,0.0,openness)
        openness=np.where(wall_height>0,openness,1.)
    else:
        openness=np.asarray(openness_override,float)
        if openness.shape != wall_height.shape or np.any(~np.isfinite(openness)):
            raise ValueError('openness_override must match finite wall heights')
        openness=np.clip(openness,0.0,1.0)
        openness=np.where(wall_height>0,openness,1.)
    if wall_blockage_fraction is not None:
        blockage=np.asarray(wall_blockage_fraction,float)
        if blockage.shape != wall_height.shape or np.any(~np.isfinite(blockage)):
            raise ValueError('wall_blockage_fraction must match finite wall heights')
        blockage=np.clip(blockage,0.0,1.0)
        # A thin obstacle can cover only part of a coarse face/cell.  The
        # uncovered part remains open while the covered part uses the same
        # finite-height overtopping/reflection closure.
        openness=1.0-blockage*(1.0-openness)
        openness=np.where(wall_height>0,openness,1.)
    shared=flux*openness
    wall_left=np.zeros_like(fl);wall_right=np.zeros_like(fr)
    wall_left[2+axis]=pl+left[0]*nl**2+speed*left[0]*nl
    wall_right[2+axis]=pr+right[0]*nr**2-speed*right[0]*nr
    left_flux=shared+(1-openness)*wall_left
    right_flux=shared+(1-openness)*wall_right
    if return_shared:
        return left_flux,right_flux,shared
    return left_flux,right_flux


def _diffuse_horizontal_through_open_faces(q, height, *, dxx, dyy, dt_s,
                                          cell_size_m, wall_x_height_m,
                                          wall_y_height_m,
                                          wall_blockage_fraction=1.0):
    """Diffuse conserved columns through the geometrically open face area.

    A common nonnegative conductance multiplies all four conserved fluxes.
    This preserves mixture/species bounds and both signed momenta.  The
    above-wall fraction uses the taller neighbouring column, as in the
    hydrostatic face reconstruction, so a tall cloud can enter a dry cell.
    Outer faces have zero flux.  The returned blocked species throughput is
    an accumulated absolute transfer diagnostic, not a unique mass inventory.
    """
    if min(dxx,dyy,dt_s)<0.0 or cell_size_m<=0.0:
        raise ValueError('diffusion coefficients and step must be nonnegative')
    if (dxx+dyy)*dt_s/(cell_size_m*cell_size_m)>.5+1.e-12:
        raise ValueError('horizontal face diffusion violates positivity CFL')
    if not 0.0<=wall_blockage_fraction<=1.0:
        raise ValueError('wall_blockage_fraction must lie in [0,1]')
    def openness(left_height,right_height,wall):
        face_height=np.maximum(left_height,right_height)
        above=np.divide(np.maximum(face_height-wall,0.0),face_height,
                        out=np.zeros_like(face_height),where=face_height>1.e-30)
        return np.where(wall>0.0,
                        1.0-wall_blockage_fraction*(1.0-above),1.0)
    open_x=openness(height[:,:-1],height[:,1:],wall_x_height_m[:,1:-1])
    open_y=openness(height[:-1,:],height[1:,:],wall_y_height_m[1:-1,:])
    tx=dxx*dt_s/(cell_size_m*cell_size_m)*(q[:,:,:-1]-q[:,:,1:])
    ty=dyy*dt_s/(cell_size_m*cell_size_m)*(q[:,:-1,:]-q[:,1:,:])
    blocked_species=float((np.sum(np.abs(tx[1])*(1.0-open_x))
                           +np.sum(np.abs(ty[1])*(1.0-open_y)))
                          *cell_size_m*cell_size_m)
    tx*=open_x[None,:,:];ty*=open_y[None,:,:]
    result=q.copy()
    result[:,:,:-1]-=tx;result[:,:,1:]+=tx
    result[:,:-1,:]-=ty;result[:,1:,:]+=ty
    return result,blocked_species


def _apply_tangential_wall_bypass(q, *, origin_xy_m, cell_size_m,
                                  wall_x_height_m, wall_y_height_m,
                                  centre_xy_m, wind_xy_m_s, dt_s,
                                  fraction, thermo, gravity_scale=1.0,
                                  wall_overtop_exponent=0.0,
                                  tip_only=False, tip_power=2.0,
                                  tip_radius_m=1.0):
    """Redirect a bounded blocked incident flux along a curved fence tangent.

    The operation transfers mixture, species, and a new tangential momentum
    state together.  The returned impulse is the wall reaction required by
    the momentum ledger.  No receptor values enter this operator.
    """
    fraction=float(fraction)
    if fraction <= 0.0:
        return q, np.zeros(2), 0.0
    h=float(cell_size_m); dt=float(dt_s)
    ny,nx=q.shape[1:]
    wx=np.asarray(wall_x_height_m,float); wy=np.asarray(wall_y_height_m,float)
    wind=np.asarray(wind_xy_m_s,float); cx,cy=(float(v) for v in centre_xy_m)
    height,u,v,_,_,_=_primitive(q,thermo,gravity_scale)
    available=q[0].copy(); delta=np.zeros_like(q); impulse=np.zeros(2); moved=0.0
    def transfer(donor,target,desired,tangent):
        nonlocal moved
        if donor==target or not (0<=donor[0]<ny and 0<=donor[1]<nx
                                 and 0<=target[0]<ny and 0<=target[1]<nx):
            return
        amount=min(float(desired),float(available[donor]))
        if amount<=0.0:return
        dm=float(q[0,donor[0],donor[1]])
        old_velocity=np.array([
            q[2,donor[0],donor[1]]/dm if dm>1e-30 else 0.0,
            q[3,donor[0],donor[1]]/dm if dm>1e-30 else 0.0])
        tangent=np.asarray(tangent,float)
        tangent/=max(float(np.linalg.norm(tangent)),1e-30)
        new_velocity=tangent*float(np.linalg.norm(wind))
        species=q[1,donor[0],donor[1]]/dm if dm>1e-30 else 0.0
        delta[0,donor[0],donor[1]]-=amount;delta[1,donor[0],donor[1]]-=amount*species
        delta[2,donor[0],donor[1]]-=amount*old_velocity[0];delta[3,donor[0],donor[1]]-=amount*old_velocity[1]
        delta[0,target[0],target[1]]+=amount;delta[1,target[0],target[1]]+=amount*species
        delta[2,target[0],target[1]]+=amount*new_velocity[0];delta[3,target[0],target[1]]+=amount*new_velocity[1]
        available[donor]-=amount; moved+=amount*h*h
        impulse[:]+=amount*h*h*(new_velocity-old_velocity)
    # x-normal faces
    for j,i in zip(*np.nonzero(wx[:,1:-1]>0.0)):
        donor=(j,i if wind[0]>=0.0 else i+1)
        H=max(float(height[j,i]),float(height[j,i+1]))
        wall=float(wx[j,i+1]); open_fraction=0.0 if H<=1e-30 else max(H-wall,0.0)/H
        open_fraction=open_fraction**(float(wall_overtop_exponent)+1.0)
        px=float(origin_xy_m[0]+(i+1)*h); py=float(origin_xy_m[1]+(j+.5)*h)
        radial=np.array([px-cx,py-cy]); tangent=np.array([-radial[1],radial[0]])
        if np.dot(tangent,wind)<0.0:tangent=-tangent
        side=1 if tangent[1]>=0.0 else -1
        target=(donor[0]+side,donor[1])
        desired=dt/h*max(float(q[0,donor[0],donor[1]])*abs(float(u[donor[0],donor[1]])),0.0)
        if tip_only:
            tip_weight=np.clip(abs(float(radial[1]))/max(float(tip_radius_m),1e-30),0.0,1.0)
            desired*=tip_weight**float(tip_power)
        transfer(donor,target,desired*fraction*(1.0-open_fraction),tangent)
    # y-normal faces
    for j,i in zip(*np.nonzero(wy[1:-1,:]>0.0)):
        donor=(j if wind[1]>=0.0 else j+1,i)
        H=max(float(height[j,i]),float(height[j+1,i]))
        wall=float(wy[j+1,i]); open_fraction=0.0 if H<=1e-30 else max(H-wall,0.0)/H
        open_fraction=open_fraction**(float(wall_overtop_exponent)+1.0)
        px=float(origin_xy_m[0]+(i+.5)*h); py=float(origin_xy_m[1]+(j+1)*h)
        radial=np.array([px-cx,py-cy]); tangent=np.array([-radial[1],radial[0]])
        if np.dot(tangent,wind)<0.0:tangent=-tangent
        side=1 if tangent[0]>=0.0 else -1
        target=(donor[0],donor[1]+side)
        desired=dt/h*max(float(q[0,donor[0],donor[1]])*abs(float(v[donor[0],donor[1]])),0.0)
        if tip_only:
            tip_weight=np.clip(abs(float(radial[0]))/max(float(tip_radius_m),1e-30),0.0,1.0)
            desired*=tip_weight**float(tip_power)
        transfer(donor,target,desired*fraction*(1.0-open_fraction),tangent)
    q=q+delta
    if np.any(q[:2]<-1e-10) or np.any(q[1]>q[0]+1e-10):
        raise RuntimeError('tangential wall bypass violated conserved positivity')
    return q,impulse,moved


def _transport_multilayer_relative(vertical_q, relative_u, relative_v, *,
                                   dt_s, cell_size_m, wall_x_height_m,
                                   wall_y_height_m, z_edges_m):
    """Move fixed-height layers relative to the column carrier conservatively.

    The hydrostatic column solve remains the mean carrier.  This correction
    transports every conserved component of each fixed-height layer with its
    *relative* horizontal velocity.  Internal transfers cancel exactly and
    the outer boundary is closed, so this operator cannot create mixture,
    species, or momentum.  A finite wall admits only the geometric fraction
    of a layer above its height.
    """
    state=np.asarray(vertical_q,float)
    ru=np.asarray(relative_u,float);rv=np.asarray(relative_v,float)
    if state.ndim!=4 or state.shape[0]!=4:
        raise ValueError('vertical_q must have shape (4,layers,ny,nx)')
    layers,ny,nx=state.shape[1:]
    if ru.shape!=(layers,ny,nx) or rv.shape!=(layers,ny,nx):
        raise ValueError('relative layer velocities must have shape (layers,ny,nx)')
    wx=np.asarray(wall_x_height_m,float);wy=np.asarray(wall_y_height_m,float)
    if wx.shape!=(ny,nx+1) or wy.shape!=(ny+1,nx):
        raise ValueError('wall grids do not match multilayer state')
    z_edges=np.asarray(z_edges_m,float)
    if z_edges.shape!=(layers+1,) or np.any(np.diff(z_edges)<=0.0):
        raise ValueError('z_edges_m must be increasing and match layer count')
    if dt_s<=0.0 or cell_size_m<=0.0:
        raise ValueError('positive dt_s and cell_size_m are required')
    if np.any(~np.isfinite(state)) or np.any(~np.isfinite(ru)) or np.any(~np.isfinite(rv)):
        raise ValueError('multilayer relative transport inputs must be finite')

    # Use a positivity-preserving donor-cell update.  Subcycle only this
    # inexpensive correction when vertical shear is stronger than the column
    # CFL; the main hydrostatic step need not be shortened.
    rate=(float(np.max(np.abs(ru)))+float(np.max(np.abs(rv))))/cell_size_m
    substeps=max(1,int(math.ceil(float(dt_s)*rate/.42)))
    sub_dt=float(dt_s)/substeps
    moved=0.0
    for _ in range(substeps):
        flux_x=np.zeros((4,layers,ny,nx+1),dtype=float)
        flux_y=np.zeros((4,layers,ny+1,nx),dtype=float)
        ux=.5*(ru[:,:,:-1]+ru[:,:,1:])
        vy=.5*(rv[:,:-1,:]+rv[:,1:,:])
        donor_x=np.where(ux[None,:,:,:]>=0.0,
                         state[:,:,:,:-1],state[:,:,:,1:])
        donor_y=np.where(vy[None,:,:,:]>=0.0,
                         state[:,:,:-1,:],state[:,:,1:,:])
        dz=(z_edges[1:]-z_edges[:-1])[:,None,None]
        open_x=np.clip((z_edges[1:,None,None]
                        -np.maximum(z_edges[:-1,None,None],wx[None,:,1:nx]))
                       /dz,0.0,1.0)
        open_y=np.clip((z_edges[1:,None,None]
                        -np.maximum(z_edges[:-1,None,None],wy[None,1:ny,:]))
                       /dz,0.0,1.0)
        open_x=np.where(wx[None,:,1:nx]>0.0,open_x,1.0)
        open_y=np.where(wy[None,1:ny,:]>0.0,open_y,1.0)
        flux_x[:,:,:,1:nx]=ux[None,:,:,:]*donor_x*open_x[None,:,:,:]
        flux_y[:,:,1:ny,:]=vy[None,:,:,:]*donor_y*open_y[None,:,:,:]
        state-=sub_dt/cell_size_m*(
            flux_x[:,:,:,1:]-flux_x[:,:,:,:-1]
            +flux_y[:,:,1:,:]-flux_y[:,:,:-1,:])
        moved+=float((np.abs(flux_x[0]).sum()+np.abs(flux_y[0]).sum())
                     *sub_dt*cell_size_m)
        if (np.any(~np.isfinite(state)) or np.min(state[:2])<-1.e-10
                or np.max(state[1]-state[0])>1.e-10):
            raise RuntimeError(
                'multilayer relative transport positivity gate failed: '
                f'min_mass={float(np.min(state[0])):.17g}, '
                f'min_species={float(np.min(state[1])):.17g}, '
                f'max_species_minus_mass={float(np.max(state[1]-state[0])):.17g}, '
                f'sub_dt_s={sub_dt:.17g}, rate_s_inv={rate:.17g}')
    return state,moved,substeps


def _repair_multilayer_species_bounds(vertical_q):
    """Project layer species into ``0<=s<=m`` without changing a column sum.

    Independent Rusanov mixture/species diffusion can leave an admissible
    column with a tiny-to-moderate layer-local excess.  Clipping that excess
    would destroy species.  This bounded water-fill moves only the invalid
    amount to available mass in another height layer and records the moved
    inventory for diagnostics.
    """
    state=np.asarray(vertical_q,float)
    mass=np.maximum(state[0],0.0)
    species=state[1]
    target=species.sum(axis=0)
    bounded=np.clip(species,0.0,mass)
    difference=target-bounded.sum(axis=0)
    add=np.maximum(difference,0.0)
    capacity=np.maximum(mass-bounded,0.0)
    capacity_total=capacity.sum(axis=0)
    bounded+=capacity*np.divide(
        add,capacity_total,out=np.zeros_like(add),where=capacity_total>1.e-30)[None,:,:]
    remove=np.maximum(-difference,0.0)
    available=bounded.sum(axis=0)
    bounded-=bounded*np.divide(
        remove,available,out=np.zeros_like(remove),where=available>1.e-30)[None,:,:]
    residual=target-bounded.sum(axis=0)
    if np.max(np.abs(residual))>1.e-10:
        raise RuntimeError('multilayer species bound repair could not close column inventory')
    moved=.5*float(np.abs(bounded-species).sum())
    state[1]=bounded
    return state,moved


def solve_dense_layer_patch(*,initial_conserved,thermo:BinaryPuffThermodynamics,
                            origin_xy_m,cell_size_m,duration_s,wind_xy_m_s,
                            wall_x_height_m,wall_y_height_m,receptor_points_m,
                            max_dt_s=.25,entrainment_coefficient=.1,
                            entrainment_far_coefficient=None,
                            entrainment_relaxation_time_s=60.0,
                            drag_coefficient=.2,
                            friction_velocity_m_s=None,roughness_length_m=.005,
                            dense_entrainment_model='constant',gravity_scale=1.,
                            vertical_observation_profile='top_hat',
                            dense_vertical_exponent=0.0,
                            vertical_grounding_exponent=0.0,
                            vertical_elevated_fraction=0.0,
                            wall_overtop_exponent=0.0,
                            vertical_tail_scale_m=1.0,
                            vertical_tail_extension_m=0.0,
                            vertical_tail_diffusivity_m2_s=1.0,
                            mass_conserving_vertical_profile=False,
                            vertical_observation_mode='legacy_mass_fraction',
                            vertical_observation_bins=64,
                            vertical_observation_component_diagnostics=False,
                            vertical_layers=1,
                            vertical_mixing_diffusivity_m2_s=0.0,
                             multilayer_entrainment_placement='top',
                            multilayer_horizontal_shear_fraction=0.0,
                            multilayer_obstacle_crosswind_shear_rate_s_inv=0.0,
                            multilayer_obstacle_crosswind_decay_length_m=70.0,
                            multilayer_relative_transport_interval_s=0.0,
                            obstacle_vertical_lift_fraction=0.0,
                            obstacle_vertical_lift_decay_m=30.0,
                            wall_impulse_vertical_lift_fraction=0.0,
                            wall_impulse_elevated_transfer_fraction=0.0,
                             wall_impulse_elevated_state_mode='species_marker',
                             wall_impulse_elevated_conserved_fraction=1.0,
                            wall_flux_elevated_transfer_fraction=0.0,
                            wall_impulse_elevated_transport_fraction=0.0,
        wall_impulse_elevated_crosswind_velocity_m_s=0.0,
        wall_impulse_elevated_crosswind_diffusivity_m2_s=0.0,
        wall_impulse_elevated_crosswind_outward_velocity_m_s=0.0,
        wall_impulse_elevated_crosswind_outward_scale_m=50.0,
        wall_impulse_elevated_crosswind_outward_decay_s=0.0,
        wall_impulse_elevated_crosswind_spinup_s=0.0,
        wall_impulse_elevated_transport_interval_s=0.0,
        wall_impulse_elevated_advection_scale=1.0,
        wall_impulse_elevated_advection_duration_s=0.0,
                            wall_impulse_elevated_return_rate_s_inv=0.0,
                            wall_impulse_elevated_center_m=6.0,
                            wall_impulse_elevated_sigma_m=2.0,
                            wall_impulse_elevated_observation_ramp_m=0.0,
                             wall_impulse_elevated_settling_velocity_m_s=0.0,
                             wall_impulse_elevated_vertical_diffusivity_m2_s=0.0,
                             wall_impulse_elevated_min_center_m=0.0,
                             wall_impulse_elevated_rise_velocity_m_s=0.0,
                             wall_impulse_elevated_rise_time_s=1.0,
                             wall_impulse_elevated_conserved_center_m=None,
                             wall_impulse_elevated_conserved_sigma_m=None,
                             wall_impulse_elevated_conserved_settling_velocity_m_s=None,
                             wall_impulse_elevated_conserved_vertical_diffusivity_m2_s=None,
                             wall_impulse_elevated_conserved_min_center_m=None,
                             wall_impulse_elevated_conserved_rise_velocity_m_s=0.0,
                             wall_impulse_elevated_conserved_rise_time_s=1.0,
                             wall_impulse_elevated_conserved_transport_fraction=None,
                             wall_impulse_elevated_conserved_crosswind_velocity_m_s=None,
                             wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s=None,
                             wall_impulse_elevated_conserved_crosswind_outward_scale_m=None,
                             wall_impulse_elevated_conserved_center_mode='fixed',
                             wall_impulse_reattached_fraction=0.0,
                             wall_impulse_reattached_transport_fraction=0.0,
                             wall_impulse_reattached_transport_height_m=None,
                             wall_impulse_reattached_crosswind_velocity_m_s=0.0,
                             wall_impulse_reattached_crosswind_outward_velocity_m_s=0.0,
                             wall_impulse_reattached_crosswind_outward_scale_m=50.0,
                             wall_impulse_reattached_return_rate_s_inv=0.0,
                             wall_impulse_reattached_center_m=0.7,
                             wall_impulse_reattached_sigma_m=0.35,
                             wall_impulse_reattached_settling_velocity_m_s=0.002,
                             wall_impulse_reattached_vertical_diffusivity_m2_s=0.0,
                             wall_impulse_reattached_min_center_m=0.0,
                            wall_impulse_wake_capture_rate_s_inv=0.0,
                            wall_impulse_wake_capture_mode='bulk',
                            wall_impulse_wake_return_rate_s_inv=0.0,
                            wall_impulse_wake_transport_fraction=0.0,
                            wall_impulse_wake_crosswind_velocity_m_s=0.0,
                            wall_impulse_wake_crosswind_outward_velocity_m_s=0.0,
                            wall_impulse_wake_crosswind_outward_scale_m=50.0,
                            wall_impulse_wake_observation_fraction=1.0,
                            wall_impulse_wake_observation_center_m=None,
                            wall_impulse_wake_observation_sigma_m=2.0,
                             wall_impulse_wake_settling_velocity_m_s=0.0,
                             wall_impulse_wake_vertical_diffusivity_m2_s=0.0,
                             wall_impulse_wake_min_center_m=0.0,
                            initial_elevated_species_fraction=0.0,
                            dense_settling_velocity_m_s=0.0,
                            height_resolved_wall_flux=False,
                            receptor_interpolation='nearest',
                            initial_vertical_gaussian=None,
                            vertical_domain_top_m=None,
                            horizontal_diffusivity_m2_s=0.0,
                            horizontal_diffusion_respect_walls=False,
                            crosswind_diffusivity_factor=1.0,
                            wake_deficit_strength=0.0,wake_length_m=1.0,
                            wake_base_width_m=1.0,wake_spread_slope=0.0,
                            wake_center_xy_m=None,wake_radius_m=0.0,
                            wake_recirculation_strength=0.0,
                            wake_meander_strength=0.0,
                            wake_shedding_strouhal=0.2,
                            wake_meander_phase_deg=0.0,
                            cavity_shelter_strength=0.0,
                            cavity_entrainment_factor=1.0,
                            cavity_core_entrainment_reduction=0.0,
                            cavity_vortex_rate_s_inv=0.0,
                            cavity_vortex_strength=0.0,
                            cavity_reservoir_capture_rate_s_inv=0.0,
                            cavity_reservoir_return_rate_s_inv=0.0,
                            cavity_reservoir_return_target='main',
                            cavity_reservoir_reattached_fraction=1.0,
                            cavity_reservoir_capture_mode='cavity_rate',
                            cavity_reservoir_vertical_capture_mode='column',
                            cavity_reservoir_spatial_mode='well_mixed',
                            cavity_reservoir_content_mode='conserved',
                            cavity_reservoir_transport_fraction=0.0,
                            cavity_reservoir_transport_height_m=None,
                             cavity_reservoir_transport_mask='cavity',
                             cavity_reservoir_wake_length_m=0.0,
                             cavity_reservoir_wake_half_width_m=None,
                             cavity_reservoir_crosswind_convergence_rate_s_inv=0.0,
                             cavity_reservoir_crosswind_target_offset_m=0.0,
                            local_reservoir_observation_fraction=0.0,
                            local_reservoir_observation_center_m=None,
                            local_reservoir_observation_sigma_m=2.0,
                            cavity_reservoir_return_sigma_m=None,
                            cavity_reservoir_return_radius_m=0.0,
                            cavity_reservoir_return_offset_m=0.0,
                            cavity_reservoir_return_crosswind_scale=1.0,
                            cavity_reservoir_face_return_fraction=0.0,
                            cavity_reservoir_capture_radius_m=0.0,
                            cavity_reservoir_start_time_s=0.0,
                            cavity_reservoir_return_delay_s=0.0,
                            cavity_return_streamwise_bias=0.0,
                            wall_blockage_fraction=1.0,
                            cavity_vertical_exponent=0.0,
                            cavity_vertical_radial_decay_power=0.0,
                            cavity_vertical_mixing_factor=0.0,
                            exterior_entrainment_far_factor=1.0,
                            exterior_entrainment_decay_length_m=100.0,
                            cavity_backflow_strength=0.0,
                            cavity_backflow_rate_s_inv=0.0,
                            cavity_backflow_frame='cavity',
                            interface_x_m=None,
                            interface_source_callback: Callable | None = None,
                            interface_time_offset_s=0.0,
                            source_injection_callback: Callable | None = None,
                            cavity_stagnation_strength=0.0,
                            exterior_reentry_fraction=0.0,
                            exterior_reentry_delay_s=0.0,
                            exterior_reentry_sigma_m=15.0,
                            wall_tangential_bypass_fraction=0.0,
                            wall_tip_bypass_fraction=0.0,
                            wall_tip_bypass_power=2.0,
                            wall_overflow_capture_fraction=0.0,
                            wall_overflow_return_rate_s_inv=0.0,
                            wall_overflow_side_gate_length_m=0.0,
                            cavity_geometry='full_disk',
                            cavity_center_xy_m=None,cavity_radius_m=0.0,
                            cavity_align_with_wind=False,
                            wind_history_time_s=None,
                            wind_history_xy_m_s=None,
                            wind_history_time_offset_s=0.0,
                            wind_history_advection_reference_x_m=None,
                            wind_history_advection_speed_m_s=None,
                            wind_history_scale_friction_velocity=True):
    """Rusanov FV with dry ambient boundaries and explicit wall reactions.

    Ambient entrainment enters through layer top only, with ambient momentum.
    Existing integral .1 shear and .2 form-drag priors are explicit hypotheses
    here, not fitted experimental coefficients. No gravity force is clipped.
    """
    started=time.perf_counter()
    if dense_entrainment_model not in {'constant','parker_bulk_suppression'}:
        raise ValueError('unsupported dense_entrainment_model')
    if (not math.isfinite(entrainment_coefficient)
            or entrainment_coefficient < 0.0):
        raise ValueError('entrainment_coefficient must be finite and non-negative')
    if entrainment_far_coefficient is not None:
        if (not math.isfinite(entrainment_far_coefficient)
                or entrainment_far_coefficient < 0.0):
            raise ValueError(
                'entrainment_far_coefficient must be finite and non-negative')
        if (not math.isfinite(entrainment_relaxation_time_s)
                or entrainment_relaxation_time_s <= 0.0):
            raise ValueError(
                'entrainment_relaxation_time_s must be positive and finite')
    if not math.isfinite(gravity_scale) or gravity_scale<0:
        raise ValueError('gravity_scale must be finite and nonnegative')
    if vertical_observation_profile not in {'top_hat','dense_power','exponential_tail','diffusive_tail','handoff_gaussian','multilayer'}:
        raise ValueError('unsupported vertical observation profile')
    if not math.isfinite(dense_vertical_exponent) or dense_vertical_exponent<0:
        raise ValueError('dense_vertical_exponent must be finite and non-negative')
    if not math.isfinite(vertical_grounding_exponent) or vertical_grounding_exponent<0:
        raise ValueError('vertical_grounding_exponent must be finite and non-negative')
    if not math.isfinite(vertical_elevated_fraction) or not 0.0 <= vertical_elevated_fraction <= 1.0:
        raise ValueError('vertical_elevated_fraction must lie in [0,1]')
    if not math.isfinite(wall_overtop_exponent) or wall_overtop_exponent < 0.0:
        raise ValueError('wall_overtop_exponent must be finite and non-negative')
    for value,name in (
            (wall_impulse_elevated_settling_velocity_m_s,
             'wall_impulse_elevated_settling_velocity_m_s'),
            (wall_impulse_elevated_vertical_diffusivity_m2_s,
             'wall_impulse_elevated_vertical_diffusivity_m2_s'),
            (wall_impulse_elevated_min_center_m,
             'wall_impulse_elevated_min_center_m'),
            (wall_impulse_elevated_rise_velocity_m_s,
             'wall_impulse_elevated_rise_velocity_m_s')):
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f'{name} must be finite and non-negative')
    if (not math.isfinite(wall_impulse_elevated_rise_time_s)
            or wall_impulse_elevated_rise_time_s <= 0.0):
        raise ValueError(
            'wall_impulse_elevated_rise_time_s must be finite and positive')
    conserved_vertical_values=(
        (wall_impulse_elevated_conserved_center_m,
         'wall_impulse_elevated_conserved_center_m',True),
        (wall_impulse_elevated_conserved_sigma_m,
         'wall_impulse_elevated_conserved_sigma_m',True),
        (wall_impulse_elevated_conserved_settling_velocity_m_s,
         'wall_impulse_elevated_conserved_settling_velocity_m_s',False),
        (wall_impulse_elevated_conserved_vertical_diffusivity_m2_s,
         'wall_impulse_elevated_conserved_vertical_diffusivity_m2_s',False),
        (wall_impulse_elevated_conserved_min_center_m,
         'wall_impulse_elevated_conserved_min_center_m',False),
        (wall_impulse_elevated_conserved_rise_velocity_m_s,
         'wall_impulse_elevated_conserved_rise_velocity_m_s',False),
        (wall_impulse_elevated_conserved_rise_time_s,
         'wall_impulse_elevated_conserved_rise_time_s',True),
    )
    for value,name,strictly_positive in conserved_vertical_values:
        if value is None:
            continue
        if (not math.isfinite(value) or value < 0.0
                or (strictly_positive and value <= 0.0)):
            qualifier='positive' if strictly_positive else 'non-negative'
            raise ValueError(f'{name} must be finite and {qualifier}')
    conserved_prompt_center=(
        float(wall_impulse_elevated_center_m)
        if wall_impulse_elevated_conserved_center_m is None
        else float(wall_impulse_elevated_conserved_center_m))
    conserved_prompt_sigma=(
        float(wall_impulse_elevated_sigma_m)
        if wall_impulse_elevated_conserved_sigma_m is None
        else float(wall_impulse_elevated_conserved_sigma_m))
    conserved_prompt_settling=(
        float(wall_impulse_elevated_settling_velocity_m_s)
        if wall_impulse_elevated_conserved_settling_velocity_m_s is None
        else float(wall_impulse_elevated_conserved_settling_velocity_m_s))
    conserved_prompt_diffusivity=(
        float(wall_impulse_elevated_vertical_diffusivity_m2_s)
        if wall_impulse_elevated_conserved_vertical_diffusivity_m2_s is None
        else float(wall_impulse_elevated_conserved_vertical_diffusivity_m2_s))
    conserved_prompt_min_center=(
        float(wall_impulse_elevated_min_center_m)
        if wall_impulse_elevated_conserved_min_center_m is None
        else float(wall_impulse_elevated_conserved_min_center_m))
    conserved_prompt_transport_fraction=(
        float(wall_impulse_elevated_transport_fraction)
        if wall_impulse_elevated_conserved_transport_fraction is None
        else float(wall_impulse_elevated_conserved_transport_fraction))
    conserved_prompt_crosswind_velocity=(
        float(wall_impulse_elevated_crosswind_velocity_m_s)
        if wall_impulse_elevated_conserved_crosswind_velocity_m_s is None
        else float(wall_impulse_elevated_conserved_crosswind_velocity_m_s))
    conserved_prompt_outward_velocity=(
        float(wall_impulse_elevated_crosswind_outward_velocity_m_s)
        if wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s is None
        else float(
            wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s))
    conserved_prompt_outward_scale=(
        float(wall_impulse_elevated_crosswind_outward_scale_m)
        if wall_impulse_elevated_conserved_crosswind_outward_scale_m is None
        else float(wall_impulse_elevated_conserved_crosswind_outward_scale_m))
    if wall_impulse_elevated_conserved_center_mode not in {'fixed','impact_head'}:
        raise ValueError(
            'wall_impulse_elevated_conserved_center_mode must be fixed or '
            'impact_head')
    if (not 0.0<=conserved_prompt_transport_fraction<=1.0
            or not math.isfinite(conserved_prompt_transport_fraction)):
        raise ValueError(
            'wall_impulse_elevated_conserved_transport_fraction must lie in [0,1]')
    if not math.isfinite(conserved_prompt_crosswind_velocity):
        raise ValueError(
            'wall_impulse_elevated_conserved_crosswind_velocity_m_s must be finite')
    if (not math.isfinite(conserved_prompt_outward_velocity)
            or conserved_prompt_outward_velocity<0.0):
        raise ValueError(
            'wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s '
            'must be finite and non-negative')
    if (not math.isfinite(conserved_prompt_outward_scale)
            or conserved_prompt_outward_scale<=0.0):
        raise ValueError(
            'wall_impulse_elevated_conserved_crosswind_outward_scale_m must '
            'be positive and finite')
    for value,name in (
            (wall_impulse_reattached_fraction,
             'wall_impulse_reattached_fraction'),
            (wall_impulse_reattached_transport_fraction,
             'wall_impulse_reattached_transport_fraction')):
        if not math.isfinite(value) or not 0.0<=value<=1.0:
            raise ValueError(f'{name} must lie in [0,1]')
    if wall_impulse_reattached_transport_height_m is not None:
        if (not math.isfinite(wall_impulse_reattached_transport_height_m)
                or wall_impulse_reattached_transport_height_m<=0.0):
            raise ValueError(
                'wall_impulse_reattached_transport_height_m must be positive '
                'and finite')
        if friction_velocity_m_s is None:
            raise ValueError(
                'height-based reattached transport requires friction velocity')
    for value,name,strictly_positive in (
            (wall_impulse_reattached_crosswind_outward_velocity_m_s,
             'wall_impulse_reattached_crosswind_outward_velocity_m_s',False),
            (wall_impulse_reattached_crosswind_outward_scale_m,
             'wall_impulse_reattached_crosswind_outward_scale_m',True),
            (wall_impulse_reattached_return_rate_s_inv,
             'wall_impulse_reattached_return_rate_s_inv',False),
            (wall_impulse_reattached_center_m,
             'wall_impulse_reattached_center_m',True),
            (wall_impulse_reattached_sigma_m,
             'wall_impulse_reattached_sigma_m',True),
            (wall_impulse_reattached_settling_velocity_m_s,
             'wall_impulse_reattached_settling_velocity_m_s',False),
            (wall_impulse_reattached_vertical_diffusivity_m2_s,
             'wall_impulse_reattached_vertical_diffusivity_m2_s',False),
            (wall_impulse_reattached_min_center_m,
             'wall_impulse_reattached_min_center_m',False)):
        if (not math.isfinite(value) or value<0.0
                or (strictly_positive and value<=0.0)):
            qualifier='positive' if strictly_positive else 'non-negative'
            raise ValueError(f'{name} must be finite and {qualifier}')
    if not math.isfinite(wall_impulse_reattached_crosswind_velocity_m_s):
        raise ValueError(
            'wall_impulse_reattached_crosswind_velocity_m_s must be finite')
    if (not math.isfinite(wall_blockage_fraction)
            or not 0.0 <= wall_blockage_fraction <= 1.0):
        raise ValueError('wall_blockage_fraction must lie in [0,1]')
    if not math.isfinite(cavity_vertical_exponent) or cavity_vertical_exponent<0:
        raise ValueError('cavity_vertical_exponent must be finite and non-negative')
    if (not math.isfinite(cavity_vertical_radial_decay_power)
            or cavity_vertical_radial_decay_power < 0.0):
        raise ValueError(
            'cavity_vertical_radial_decay_power must be finite and non-negative')
    if cavity_reservoir_spatial_mode not in {'well_mixed','local'}:
        raise ValueError(
            "cavity_reservoir_spatial_mode must be 'well_mixed' or 'local'")
    if cavity_reservoir_content_mode not in {'conserved','species_only'}:
        raise ValueError(
            "cavity_reservoir_content_mode must be 'conserved' or 'species_only'")
    if cavity_reservoir_return_target not in {'main','reattached'}:
        raise ValueError(
            "cavity_reservoir_return_target must be 'main' or 'reattached'")
    if (not math.isfinite(cavity_reservoir_reattached_fraction)
            or not 0.0<=cavity_reservoir_reattached_fraction<=1.0):
        raise ValueError(
            'cavity_reservoir_reattached_fraction must lie in [0,1]')
    if cavity_reservoir_return_target=='reattached':
        if cavity_reservoir_spatial_mode!='local':
            raise ValueError('reattached reservoir return requires local spatial mode')
        if cavity_reservoir_content_mode!='conserved':
            raise ValueError('reattached reservoir return requires conserved content')
        if (wall_impulse_elevated_state_mode not in {'conserved','hybrid'}
                or wall_impulse_elevated_transfer_fraction<=0.0
                or wall_impulse_reattached_fraction<=0.0):
            raise ValueError(
                'reattached reservoir return requires an active conserved '
                'reattached branch')
    if cavity_reservoir_capture_mode not in {'cavity_rate','wall_impact'}:
        raise ValueError(
            "cavity_reservoir_capture_mode must be 'cavity_rate' or 'wall_impact'")
    if cavity_reservoir_vertical_capture_mode not in {'column','below_wall_core'}:
        raise ValueError(
            "cavity_reservoir_vertical_capture_mode must be 'column' or "
            "'below_wall_core'")
    if (cavity_reservoir_vertical_capture_mode=='below_wall_core'
            and cavity_reservoir_capture_mode!='wall_impact'):
        raise ValueError(
            "below_wall_core vertical capture requires wall_impact capture")
    if (cavity_reservoir_capture_mode=='wall_impact'
            and cavity_reservoir_spatial_mode!='local'):
        raise ValueError('wall-impact reservoir capture requires local spatial mode')
    if (cavity_reservoir_content_mode=='species_only'
            and cavity_reservoir_spatial_mode!='local'):
        raise ValueError('species-only reservoir requires local spatial mode')
    if (not math.isfinite(cavity_reservoir_transport_fraction)
            or not 0.0<=cavity_reservoir_transport_fraction<=1.0):
        raise ValueError('cavity_reservoir_transport_fraction must lie in [0,1]')
    if (cavity_reservoir_transport_fraction>0.0
            and cavity_reservoir_spatial_mode!='local'):
        raise ValueError('reservoir transport requires local spatial mode')
    if cavity_reservoir_transport_height_m is not None:
        if (not math.isfinite(cavity_reservoir_transport_height_m)
                or cavity_reservoir_transport_height_m<=0.0):
            raise ValueError(
                'cavity_reservoir_transport_height_m must be positive and finite')
        if cavity_reservoir_spatial_mode!='local':
            raise ValueError(
                'height-based reservoir transport requires local spatial mode')
        if friction_velocity_m_s is None:
            raise ValueError(
                'height-based reservoir transport requires friction velocity')
    if (not math.isfinite(cavity_reservoir_crosswind_convergence_rate_s_inv)
            or cavity_reservoir_crosswind_convergence_rate_s_inv<0.0):
        raise ValueError(
            'cavity_reservoir_crosswind_convergence_rate_s_inv must be finite and non-negative')
    if not math.isfinite(cavity_reservoir_crosswind_target_offset_m):
        raise ValueError(
            'cavity_reservoir_crosswind_target_offset_m must be finite')
    if (cavity_reservoir_crosswind_convergence_rate_s_inv>0.0
            and cavity_reservoir_spatial_mode!='local'):
        raise ValueError('reservoir crosswind convergence requires local spatial mode')
    if cavity_reservoir_transport_mask not in {'cavity','wake'}:
        raise ValueError("cavity_reservoir_transport_mask must be 'cavity' or 'wake'")
    if (not math.isfinite(cavity_reservoir_wake_length_m)
            or cavity_reservoir_wake_length_m<0.0):
        raise ValueError('cavity_reservoir_wake_length_m must be finite and non-negative')
    if cavity_reservoir_wake_half_width_m is not None and (
            not math.isfinite(cavity_reservoir_wake_half_width_m)
            or cavity_reservoir_wake_half_width_m<=0.0):
        raise ValueError(
            'cavity_reservoir_wake_half_width_m must be positive and finite')
    if (cavity_reservoir_transport_mask=='wake'
            and cavity_reservoir_wake_length_m<=0.0):
        raise ValueError('wake reservoir transport requires positive wake length')
    if (cavity_reservoir_transport_mask=='wake'
            and cavity_center_xy_m is None):
        raise ValueError('wake reservoir transport requires a cavity centre')
    if (not math.isfinite(local_reservoir_observation_fraction)
            or not 0.0<=local_reservoir_observation_fraction<=1.0):
        raise ValueError('local_reservoir_observation_fraction must lie in [0,1]')
    if (local_reservoir_observation_fraction>0.0
            and cavity_reservoir_spatial_mode!='local'):
        raise ValueError('local reservoir observation requires local spatial mode')
    if local_reservoir_observation_center_m is not None:
        if (not math.isfinite(local_reservoir_observation_center_m)
                or local_reservoir_observation_center_m<0.0
                or not math.isfinite(local_reservoir_observation_sigma_m)
                or local_reservoir_observation_sigma_m<=0.0):
            raise ValueError('local reservoir observation centre/sigma are invalid')
        if cavity_reservoir_spatial_mode!='local':
            raise ValueError('local reservoir vertical observation requires local mode')
    if not isinstance(cavity_align_with_wind,bool):
        raise TypeError('cavity_align_with_wind must be boolean')
    if not isinstance(wind_history_scale_friction_velocity,bool):
        raise TypeError('wind_history_scale_friction_velocity must be boolean')
    if not math.isfinite(float(wind_history_time_offset_s)):
        raise ValueError('wind_history_time_offset_s must be finite')
    if (wind_history_time_s is None) != (wind_history_xy_m_s is None):
        raise ValueError('wind history times and vectors must be provided together')
    if (not math.isfinite(cavity_return_streamwise_bias)
            or not -1.0 <= cavity_return_streamwise_bias <= 1.0):
        raise ValueError('cavity_return_streamwise_bias must lie in [-1,1]')
    if (not math.isfinite(cavity_vertical_mixing_factor)
            or cavity_vertical_mixing_factor < 0.0):
        raise ValueError('cavity_vertical_mixing_factor must be finite and non-negative')
    if (not math.isfinite(exterior_entrainment_far_factor)
            or not 0.0 <= exterior_entrainment_far_factor <= 1.0):
        raise ValueError('exterior_entrainment_far_factor must lie in [0,1]')
    if (not math.isfinite(exterior_entrainment_decay_length_m)
            or exterior_entrainment_decay_length_m <= 0.0):
        raise ValueError(
            'exterior_entrainment_decay_length_m must be positive and finite')
    if not math.isfinite(cavity_backflow_strength) or not 0.0 <= cavity_backflow_strength <= 1.0:
        raise ValueError('cavity_backflow_strength must lie in [0,1]')
    if not math.isfinite(cavity_backflow_rate_s_inv) or cavity_backflow_rate_s_inv < 0.0:
        raise ValueError('cavity_backflow_rate_s_inv must be finite and non-negative')
    if cavity_backflow_frame not in {'cavity','legacy_x'}:
        raise ValueError("cavity_backflow_frame must be 'cavity' or 'legacy_x'")
    if interface_x_m is not None and not math.isfinite(float(interface_x_m)):
        raise ValueError('interface_x_m must be finite when provided')
    if interface_source_callback is not None and not callable(interface_source_callback):
        raise TypeError('interface_source_callback must be callable')
    if not math.isfinite(float(interface_time_offset_s)):
        raise ValueError('interface_time_offset_s must be finite')
    if not math.isfinite(cavity_stagnation_strength) or not 0.0 <= cavity_stagnation_strength <= 1.0:
        raise ValueError('cavity_stagnation_strength must lie in [0,1]')
    if (not math.isfinite(wall_tangential_bypass_fraction)
            or not 0.0 <= wall_tangential_bypass_fraction <= 1.0):
        raise ValueError('wall_tangential_bypass_fraction must lie in [0,1]')
    if (not math.isfinite(wall_tip_bypass_fraction)
            or not 0.0 <= wall_tip_bypass_fraction <= 1.0):
        raise ValueError('wall_tip_bypass_fraction must lie in [0,1]')
    if not math.isfinite(wall_tip_bypass_power) or wall_tip_bypass_power < 0.0:
        raise ValueError('wall_tip_bypass_power must be finite and non-negative')
    if (not math.isfinite(wall_overflow_capture_fraction)
            or not 0.0 <= wall_overflow_capture_fraction <= 1.0):
        raise ValueError('wall_overflow_capture_fraction must lie in [0,1]')
    if (not math.isfinite(wall_overflow_return_rate_s_inv)
            or wall_overflow_return_rate_s_inv < 0.0):
        raise ValueError('wall_overflow_return_rate_s_inv must be finite and non-negative')
    if (wall_overflow_capture_fraction > 0.0
            and wall_overflow_return_rate_s_inv <= 0.0):
        raise ValueError('captured wall overflow requires a positive return rate')
    if (not math.isfinite(wall_overflow_side_gate_length_m)
            or wall_overflow_side_gate_length_m < 0.0):
        raise ValueError('wall_overflow_side_gate_length_m must be finite and non-negative')
    if (wall_overflow_side_gate_length_m > 0.0
            and (cavity_center_xy_m is None or cavity_radius_m <= 0.0)):
        raise ValueError('side-gate overflow storage requires cavity geometry')
    if cavity_geometry not in {'full_disk','downwind_semidisk','upwind_semidisk'}:
        raise ValueError('unsupported cavity_geometry')
    if (not math.isfinite(vertical_tail_scale_m) or vertical_tail_scale_m<=0
            or not math.isfinite(vertical_tail_extension_m)
            or vertical_tail_extension_m<0):
        raise ValueError('invalid vertical observation tail controls')
    if (not math.isfinite(vertical_tail_diffusivity_m2_s)
            or vertical_tail_diffusivity_m2_s<=0):
        raise ValueError('vertical_tail_diffusivity_m2_s must be positive and finite')
    if not isinstance(mass_conserving_vertical_profile,bool):
        raise TypeError('mass_conserving_vertical_profile must be boolean')
    if vertical_observation_mode not in {'legacy_mass_fraction','conservative_volume'}:
        raise ValueError('unsupported vertical_observation_mode')
    if (isinstance(vertical_observation_bins,bool)
            or not isinstance(vertical_observation_bins,int) or vertical_observation_bins<8):
        raise ValueError('vertical_observation_bins must be an integer >=8')
    if not isinstance(vertical_observation_component_diagnostics,bool):
        raise TypeError('vertical_observation_component_diagnostics must be boolean')
    if (vertical_observation_mode=='conservative_volume'
            and vertical_layers==1
            and vertical_observation_profile=='handoff_gaussian'):
        raise ValueError('conservative volume observer requires a supported single-layer profile')
    if (isinstance(vertical_layers,bool) or not isinstance(vertical_layers,int)
            or vertical_layers<1):
        raise ValueError('vertical_layers must be a positive integer')
    if (not math.isfinite(vertical_mixing_diffusivity_m2_s)
            or vertical_mixing_diffusivity_m2_s<0):
        raise ValueError('vertical_mixing_diffusivity_m2_s must be finite and non-negative')
    if (not math.isfinite(multilayer_horizontal_shear_fraction)
            or not 0.0<=multilayer_horizontal_shear_fraction<=1.0):
        raise ValueError('multilayer_horizontal_shear_fraction must lie in [0,1]')
    if not math.isfinite(multilayer_obstacle_crosswind_shear_rate_s_inv):
        raise ValueError(
            'multilayer_obstacle_crosswind_shear_rate_s_inv must be finite')
    if (not math.isfinite(multilayer_obstacle_crosswind_decay_length_m)
            or multilayer_obstacle_crosswind_decay_length_m<=0.0):
        raise ValueError(
            'multilayer_obstacle_crosswind_decay_length_m must be positive and finite')
    if (not math.isfinite(multilayer_relative_transport_interval_s)
            or multilayer_relative_transport_interval_s<0.0):
        raise ValueError(
            'multilayer_relative_transport_interval_s must be finite and non-negative')
    if (not math.isfinite(obstacle_vertical_lift_fraction)
            or obstacle_vertical_lift_fraction<0.0):
        raise ValueError('obstacle_vertical_lift_fraction must be finite and non-negative')
    if (not math.isfinite(obstacle_vertical_lift_decay_m)
            or obstacle_vertical_lift_decay_m<=0.0):
        raise ValueError('obstacle_vertical_lift_decay_m must be positive and finite')
    if (not math.isfinite(wall_impulse_vertical_lift_fraction)
            or not 0.0<=wall_impulse_vertical_lift_fraction<=1.0):
        raise ValueError('wall_impulse_vertical_lift_fraction must lie in [0,1]')
    if wall_impulse_vertical_lift_fraction>0.0 and vertical_layers<=1:
        raise ValueError('wall-impulse vertical lift requires multilayer transport')
    if (not math.isfinite(wall_impulse_elevated_transfer_fraction)
            or not 0.0<=wall_impulse_elevated_transfer_fraction<=1.0):
        raise ValueError('wall_impulse_elevated_transfer_fraction must lie in [0,1]')
    if wall_impulse_elevated_state_mode not in {
            'species_marker','conserved','hybrid'}:
        raise ValueError(
            'wall_impulse_elevated_state_mode must be species_marker, conserved, '
            'or hybrid')
    if (not math.isfinite(wall_impulse_elevated_conserved_fraction)
            or not 0.0<=wall_impulse_elevated_conserved_fraction<=1.0):
        raise ValueError(
            'wall_impulse_elevated_conserved_fraction must lie in [0,1]')
    if (not math.isfinite(wall_flux_elevated_transfer_fraction)
            or not 0.0<=wall_flux_elevated_transfer_fraction<=1.0):
        raise ValueError('wall_flux_elevated_transfer_fraction must lie in [0,1]')
    if (not math.isfinite(wall_impulse_elevated_transport_fraction)
            or not 0.0<=wall_impulse_elevated_transport_fraction<=1.0):
        raise ValueError('wall_impulse_elevated_transport_fraction must lie in [0,1]')
    if not math.isfinite(wall_impulse_elevated_crosswind_velocity_m_s):
        raise ValueError(
            'wall_impulse_elevated_crosswind_velocity_m_s must be finite')
    if (not math.isfinite(wall_impulse_elevated_crosswind_diffusivity_m2_s)
            or wall_impulse_elevated_crosswind_diffusivity_m2_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_crosswind_diffusivity_m2_s must be finite '
            'and non-negative')
    if (not math.isfinite(wall_impulse_elevated_crosswind_outward_velocity_m_s)
            or wall_impulse_elevated_crosswind_outward_velocity_m_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_crosswind_outward_velocity_m_s must be '
            'finite and non-negative')
    if (not math.isfinite(wall_impulse_elevated_crosswind_outward_scale_m)
            or wall_impulse_elevated_crosswind_outward_scale_m<=0.0):
        raise ValueError(
            'wall_impulse_elevated_crosswind_outward_scale_m must be positive '
            'and finite')
    if (not math.isfinite(wall_impulse_elevated_crosswind_outward_decay_s)
            or wall_impulse_elevated_crosswind_outward_decay_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_crosswind_outward_decay_s must be finite '
            'and non-negative')
    if (not math.isfinite(wall_impulse_elevated_crosswind_spinup_s)
            or wall_impulse_elevated_crosswind_spinup_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_crosswind_spinup_s must be finite and non-negative')
    if (not math.isfinite(wall_impulse_elevated_transport_interval_s)
            or wall_impulse_elevated_transport_interval_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_transport_interval_s must be finite and non-negative')
    if (not math.isfinite(wall_impulse_elevated_advection_scale)
            or wall_impulse_elevated_advection_scale<0.0):
        raise ValueError(
            'wall_impulse_elevated_advection_scale must be finite and non-negative')
    if (not math.isfinite(wall_impulse_elevated_advection_duration_s)
            or wall_impulse_elevated_advection_duration_s<0.0):
        raise ValueError(
            'wall_impulse_elevated_advection_duration_s must be finite and non-negative')
    if (wall_impulse_elevated_advection_scale!=1.0
            and wall_impulse_elevated_advection_duration_s<=0.0):
        raise ValueError(
            'non-unit elevated advection scale requires a positive duration')
    if (not math.isfinite(wall_impulse_elevated_return_rate_s_inv)
            or wall_impulse_elevated_return_rate_s_inv<0.0):
        raise ValueError('wall_impulse_elevated_return_rate_s_inv must be non-negative')
    for value,name in (
            (wall_impulse_wake_capture_rate_s_inv,
             'wall_impulse_wake_capture_rate_s_inv'),
            (wall_impulse_wake_return_rate_s_inv,
             'wall_impulse_wake_return_rate_s_inv')):
        if not math.isfinite(value) or value<0.0:
            raise ValueError(f'{name} must be finite and non-negative')
    if wall_impulse_wake_capture_mode not in {'bulk','species_weighted'}:
        raise ValueError(
            'wall_impulse_wake_capture_mode must be bulk or species_weighted')
    for value,name in (
            (wall_impulse_wake_transport_fraction,
             'wall_impulse_wake_transport_fraction'),
            (wall_impulse_wake_observation_fraction,
             'wall_impulse_wake_observation_fraction')):
        if not math.isfinite(value) or not 0.0<=value<=1.0:
            raise ValueError(f'{name} must lie in [0,1]')
    if not math.isfinite(wall_impulse_wake_crosswind_velocity_m_s):
        raise ValueError(
            'wall_impulse_wake_crosswind_velocity_m_s must be finite')
    if (not math.isfinite(wall_impulse_wake_crosswind_outward_velocity_m_s)
            or wall_impulse_wake_crosswind_outward_velocity_m_s<0.0):
        raise ValueError(
            'wall_impulse_wake_crosswind_outward_velocity_m_s must be finite '
            'and non-negative')
    if (not math.isfinite(wall_impulse_wake_crosswind_outward_scale_m)
            or wall_impulse_wake_crosswind_outward_scale_m<=0.0):
        raise ValueError(
            'wall_impulse_wake_crosswind_outward_scale_m must be positive '
            'and finite')
    if wall_impulse_wake_observation_center_m is not None:
        if (not math.isfinite(wall_impulse_wake_observation_center_m)
                or wall_impulse_wake_observation_center_m<0.0
                or not math.isfinite(wall_impulse_wake_observation_sigma_m)
                or wall_impulse_wake_observation_sigma_m<=0.0):
            raise ValueError('wall-impact wake observation centre/sigma are invalid')
    for value,name in (
            (wall_impulse_wake_settling_velocity_m_s,
             'wall_impulse_wake_settling_velocity_m_s'),
            (wall_impulse_wake_vertical_diffusivity_m2_s,
             'wall_impulse_wake_vertical_diffusivity_m2_s'),
            (wall_impulse_wake_min_center_m,
             'wall_impulse_wake_min_center_m')):
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f'{name} must be finite and non-negative')
    if (wall_impulse_wake_capture_rate_s_inv>0.0
            and wall_impulse_wake_return_rate_s_inv<=0.0):
        raise ValueError('wall-impact wake capture requires a positive return rate')
    if (not math.isfinite(wall_impulse_elevated_center_m)
            or wall_impulse_elevated_center_m<=0.0
            or not math.isfinite(wall_impulse_elevated_sigma_m)
            or wall_impulse_elevated_sigma_m<=0.0):
        raise ValueError('elevated branch centre and sigma must be positive')
    if (not math.isfinite(wall_impulse_elevated_observation_ramp_m)
            or wall_impulse_elevated_observation_ramp_m<0.0):
        raise ValueError('wall elevated observation ramp must be non-negative')
    if (not math.isfinite(initial_elevated_species_fraction)
            or not 0.0<=initial_elevated_species_fraction<=1.0):
        raise ValueError('initial_elevated_species_fraction must lie in [0,1]')
    if ((wall_impulse_elevated_transfer_fraction>0.0
         or wall_flux_elevated_transfer_fraction>0.0
         or initial_elevated_species_fraction>0.0) and vertical_layers>1):
        raise ValueError('elevated two-branch state and multilayer state are mutually exclusive')
    if ((wall_impulse_elevated_transfer_fraction>0.0
         or wall_flux_elevated_transfer_fraction>0.0
         or initial_elevated_species_fraction>0.0)
            and (exterior_reentry_fraction>0.0
                 or cavity_reservoir_return_delay_s>0.0
                 or cavity_reservoir_face_return_fraction>0.0)):
        raise ValueError('elevated branch does not support delayed/sector reservoirs')
    if (not math.isfinite(dense_settling_velocity_m_s)
            or dense_settling_velocity_m_s<0.0):
        raise ValueError('dense_settling_velocity_m_s must be finite and non-negative')
    if vertical_layers>1 and vertical_mixing_diffusivity_m2_s<=0:
        raise ValueError('multilayer transport requires positive vertical diffusivity')
    if vertical_layers==1 and vertical_observation_profile=='multilayer':
        raise ValueError('multilayer profile requires vertical_layers>1')
    if vertical_layers>1 and vertical_observation_profile!='multilayer':
        raise ValueError('vertical_layers>1 requires vertical_observation_profile=multilayer')
    if multilayer_entrainment_placement not in {'top','legacy_ground'}:
        raise ValueError(
            "multilayer_entrainment_placement must be 'top' or 'legacy_ground'")
    if ((multilayer_horizontal_shear_fraction>0.0
         or multilayer_obstacle_crosswind_shear_rate_s_inv!=0.0)
            and vertical_layers<=1):
        raise ValueError('height-dependent horizontal transport requires multilayer state')
    if (multilayer_obstacle_crosswind_shear_rate_s_inv!=0.0
            and (cavity_center_xy_m is None or cavity_radius_m<=0.0)):
        raise ValueError('obstacle crosswind shear requires positive cavity geometry')
    if vertical_layers>1 and wall_overflow_capture_fraction>0.0:
        raise ValueError('wall overflow storage is not yet layer-resolved')
    if vertical_layers>1 and wall_impulse_wake_capture_rate_s_inv>0.0:
        raise ValueError('wall-impact wake reservoir is not yet layer-resolved')
    if not isinstance(height_resolved_wall_flux,bool):
        raise TypeError('height_resolved_wall_flux must be boolean')
    if receptor_interpolation not in {'nearest','bilinear'}:
        raise ValueError("receptor_interpolation must be 'nearest' or 'bilinear'")
    if initial_vertical_gaussian is not None:
        if len(initial_vertical_gaussian)!=2 or any(not math.isfinite(float(v)) for v in initial_vertical_gaussian):
            raise ValueError('initial_vertical_gaussian must contain finite centre and sigma')
        if float(initial_vertical_gaussian[1])<=0.0 or float(initial_vertical_gaussian[0])<0.0:
            raise ValueError('initial_vertical_gaussian requires non-negative centre and positive sigma')
    if not math.isfinite(horizontal_diffusivity_m2_s) or horizontal_diffusivity_m2_s<0:
        raise ValueError('horizontal_diffusivity_m2_s must be finite and non-negative')
    if not isinstance(horizontal_diffusion_respect_walls,bool):
        raise TypeError('horizontal_diffusion_respect_walls must be boolean')
    if not math.isfinite(crosswind_diffusivity_factor) or crosswind_diffusivity_factor <= 0.0:
        raise ValueError('crosswind_diffusivity_factor must be positive and finite')
    if vertical_domain_top_m is not None and (
            not math.isfinite(vertical_domain_top_m) or vertical_domain_top_m <= 0.0):
        raise ValueError('vertical_domain_top_m must be positive and finite')
    for value,name in ((wake_deficit_strength,'wake_deficit_strength'),
                       (wake_length_m,'wake_length_m'),
                       (wake_base_width_m,'wake_base_width_m'),
                       (wake_spread_slope,'wake_spread_slope'),
                       (wake_radius_m,'wake_radius_m'),
                       (wake_recirculation_strength,'wake_recirculation_strength'),
                       (wake_meander_strength,'wake_meander_strength'),
                       (wake_shedding_strouhal,'wake_shedding_strouhal'),
                       (cavity_shelter_strength,'cavity_shelter_strength'),
                       (cavity_backflow_rate_s_inv,'cavity_backflow_rate_s_inv'),
                       (cavity_stagnation_strength,'cavity_stagnation_strength'),
                       (exterior_reentry_fraction,'exterior_reentry_fraction'),
                       (exterior_reentry_delay_s,'exterior_reentry_delay_s'),
                       (exterior_reentry_sigma_m,'exterior_reentry_sigma_m'),
                       (cavity_entrainment_factor,'cavity_entrainment_factor'),
                       (cavity_vortex_rate_s_inv,'cavity_vortex_rate_s_inv'),
                       (cavity_vortex_strength,'cavity_vortex_strength'),
                       (cavity_reservoir_capture_rate_s_inv,'cavity_reservoir_capture_rate_s_inv'),
                       (cavity_reservoir_return_rate_s_inv,'cavity_reservoir_return_rate_s_inv'),
                       (cavity_reservoir_return_radius_m,'cavity_reservoir_return_radius_m'),
                       (cavity_reservoir_return_crosswind_scale,'cavity_reservoir_return_crosswind_scale'),
                       (cavity_reservoir_face_return_fraction,'cavity_reservoir_face_return_fraction'),
                       (cavity_reservoir_capture_radius_m,'cavity_reservoir_capture_radius_m'),
                       (cavity_reservoir_start_time_s,'cavity_reservoir_start_time_s'),
                       (cavity_reservoir_return_delay_s,'cavity_reservoir_return_delay_s'),
                       (cavity_radius_m,'cavity_radius_m')):
        if not math.isfinite(value) or value<0:
            raise ValueError(f'{name} must be finite and non-negative')
    # A recirculation return point may lie upstream of the cavity centre;
    # unlike rates and widths, this is a signed along-wind displacement.
    if not math.isfinite(cavity_reservoir_return_offset_m):
        raise ValueError('cavity_reservoir_return_offset_m must be finite')
    if cavity_reservoir_return_delay_s < 0.0:
        raise ValueError('cavity_reservoir_return_delay_s must be non-negative')
    if cavity_reservoir_return_delay_s > 0.0 and cavity_reservoir_face_return_fraction > 0.0:
        raise ValueError('delayed reservoir return does not support sector history yet')
    if (cavity_reservoir_spatial_mode=='local'
            and cavity_reservoir_return_delay_s>0.0):
        raise ValueError('local cavity reservoir does not support a separate return delay')
    if wake_deficit_strength>1 or wake_recirculation_strength>1:
        raise ValueError('wake strengths must be <=1')
    if wake_meander_strength>1:
        raise ValueError('wake_meander_strength must be <=1')
    if not math.isfinite(wake_meander_phase_deg):
        raise ValueError('wake_meander_phase_deg must be finite')
    if cavity_shelter_strength>1:
        raise ValueError('cavity_shelter_strength must be <=1')
    if cavity_entrainment_factor>1:
        raise ValueError('cavity_entrainment_factor must be <=1')
    if (not math.isfinite(cavity_core_entrainment_reduction)
            or not 0.0 <= cavity_core_entrainment_reduction <= 1.0):
        raise ValueError('cavity_core_entrainment_reduction must lie in [0,1]')
    if cavity_vortex_strength>1:
        raise ValueError('cavity_vortex_strength must be <=1')
    if cavity_vortex_rate_s_inv>0 and (cavity_center_xy_m is None or cavity_radius_m<=0):
        raise ValueError('cavity geometry is required for positive vortex rate')
    if ((cavity_reservoir_capture_rate_s_inv>0 or cavity_reservoir_return_rate_s_inv>0)
            and (cavity_center_xy_m is None or cavity_radius_m<=0)):
        raise ValueError('cavity geometry is required for positive reservoir rates')
    if cavity_reservoir_return_sigma_m is None:
        cavity_reservoir_return_sigma_m=(
            float(cavity_radius_m)/3.0 if cavity_radius_m>0 else 1.0
        )
    if ((cavity_reservoir_capture_rate_s_inv>0 or cavity_reservoir_return_rate_s_inv>0)
            and (not math.isfinite(cavity_reservoir_return_sigma_m)
                 or cavity_reservoir_return_sigma_m<=0)):
        raise ValueError('cavity_reservoir_return_sigma_m must be positive and finite')
    if cavity_reservoir_return_radius_m > cavity_radius_m:
        raise ValueError('cavity_reservoir_return_radius_m cannot exceed cavity radius')
    if cavity_reservoir_capture_radius_m > cavity_radius_m:
        raise ValueError('cavity_reservoir_capture_radius_m cannot exceed cavity radius')
    if abs(cavity_reservoir_return_offset_m) > cavity_radius_m:
        raise ValueError('absolute cavity_reservoir_return_offset_m cannot exceed cavity radius')
    if not math.isfinite(cavity_reservoir_return_crosswind_scale) or cavity_reservoir_return_crosswind_scale <= 0.0:
        raise ValueError('cavity_reservoir_return_crosswind_scale must be finite and positive')
    if (not math.isfinite(cavity_reservoir_face_return_fraction)
            or not 0.0 <= cavity_reservoir_face_return_fraction <= 1.0):
        raise ValueError('cavity_reservoir_face_return_fraction must lie in [0,1]')
    if cavity_shelter_strength>0 and (cavity_center_xy_m is None or cavity_radius_m<=0):
        raise ValueError('cavity geometry is required for positive shelter strength')
    if ((cavity_backflow_rate_s_inv>0 or cavity_backflow_strength>0)
            and (cavity_center_xy_m is None or cavity_radius_m<=0)):
        raise ValueError('cavity geometry is required for backflow')
    if cavity_stagnation_strength>0 and (cavity_center_xy_m is None or cavity_radius_m<=0):
        raise ValueError('cavity geometry is required for positive stagnation strength')
    if exterior_reentry_fraction>1.0:
        raise ValueError('exterior_reentry_fraction must be <=1')
    if exterior_reentry_fraction>0.0 and (cavity_center_xy_m is None or cavity_radius_m<=0):
        raise ValueError('cavity geometry is required for positive exterior re-entry')
    if exterior_reentry_fraction>0.0 and exterior_reentry_sigma_m<=0.0:
        raise ValueError('exterior_reentry_sigma_m must be positive for re-entry')
    if ((wake_deficit_strength>0 or wake_meander_strength>0)
            and (wake_length_m<=0 or wake_base_width_m<=0
                 or wake_center_xy_m is None or wake_radius_m<=0)):
        raise ValueError('wake geometry and positive scales are required')
    if wake_center_xy_m is not None and len(wake_center_xy_m)!=2:
        raise ValueError('wake_center_xy_m must contain two coordinates')
    q=np.asarray(initial_conserved,float).copy()
    if q.ndim!=3 or q.shape[0]!=4 or np.any(~np.isfinite(q)) or np.any(q[:2]<0) or np.any(q[1]>q[0]):
        raise ValueError('initial q requires finite (4,ny,nx), 0<=species<=mixture')
    if duration_s<=0 or max_dt_s<=0 or cell_size_m<=0 or entrainment_coefficient<0 or drag_coefficient<0:
        raise ValueError('invalid layer controls')
    ny,nx=q.shape[1:];h=cell_size_m;area=h*h
    interface_face_index=None
    if interface_x_m is not None:
        interface_face_index=int(round((float(interface_x_m)-float(origin_xy_m[0]))/h))
        if interface_face_index < 1 or interface_face_index >= nx:
            raise ValueError('interface_x_m must lie on an interior cell face')
    wx,wy=np.asarray(wall_x_height_m,float),np.asarray(wall_y_height_m,float)
    if wx.shape!=(ny,nx+1) or wy.shape!=(ny+1,nx) or np.any(wx<0) or np.any(wy<0):
        raise ValueError('wall heights must match face grids and be nonnegative')
    if np.any(~np.isfinite(wx)) or np.any(~np.isfinite(wy)):
        raise ValueError('wall heights must be finite')
    reference_wind=np.asarray(wind_xy_m_s,float)
    if friction_velocity_m_s is not None and (friction_velocity_m_s<=0 or roughness_length_m<=0):
        raise ValueError('positive friction velocity and roughness required')
    if friction_velocity_m_s is not None and np.linalg.norm(reference_wind)==0:
        raise ValueError('log wind needs a nonzero wind direction')
    points=np.asarray(receptor_points_m,float)
    if reference_wind.shape!=(2,) or points.ndim!=2 or points.shape[1]!=3 or np.any(~np.isfinite(reference_wind)) or np.any(~np.isfinite(points)) or np.any(points[:,2]<0):
        raise ValueError('invalid wind/receptors')
    x=origin_xy_m[0]+(np.arange(nx)+.5)*h
    y=origin_xy_m[1]+(np.arange(ny)+.5)*h
    spatial_wind_history=(wind_history_advection_reference_x_m is not None
                          or wind_history_advection_speed_m_s is not None)
    if spatial_wind_history:
        if (wind_history_advection_reference_x_m is None
                or wind_history_advection_speed_m_s is None
                or not math.isfinite(float(wind_history_advection_reference_x_m))
                or not math.isfinite(float(wind_history_advection_speed_m_s))
                or float(wind_history_advection_speed_m_s)<=0.0):
            raise ValueError('spatial wind history requires a finite reference x and positive advection speed')
        if wind_history_time_s is None:
            raise ValueError('spatial wind history requires a wind time series')
    history_times=history_vectors=None
    if wind_history_time_s is not None:
        history_times=np.asarray(wind_history_time_s,dtype=float)
        history_vectors=np.asarray(wind_history_xy_m_s,dtype=float)
        if (history_times.ndim!=1 or history_times.size<2
                or history_vectors.shape!=(history_times.size,2)
                or np.any(~np.isfinite(history_times))
                or np.any(~np.isfinite(history_vectors))
                or np.any(np.diff(history_times)<=0.0)
                or np.any(np.linalg.norm(history_vectors,axis=1)<=0.0)):
            raise ValueError('wind history requires increasing finite times and nonzero (n,2) vectors')
        spatial_shift=(
            (float(wind_history_advection_reference_x_m)-x)
            /float(wind_history_advection_speed_m_s)
            if spatial_wind_history else np.zeros_like(x)
        )
        required_start=float(wind_history_time_offset_s)+float(np.min(spatial_shift))
        required_stop=(float(wind_history_time_offset_s)+float(duration_s)
                       +float(np.max(spatial_shift)))
        if (history_times[0]>required_start+1.0e-12
                or history_times[-1]<required_stop-1.0e-12):
            raise ValueError('wind history does not cover the complete layer interval')

    representative_x=(float(cavity_center_xy_m[0])
                      if cavity_center_xy_m is not None else float(np.mean(x)))

    def wind_at_time(model_time_s):
        if history_times is None:
            return reference_wind.copy()
        event_time=float(model_time_s)+float(wind_history_time_offset_s)
        if spatial_wind_history:
            event_time+=(float(wind_history_advection_reference_x_m)-representative_x)/float(
                wind_history_advection_speed_m_s)
        return np.array([
            np.interp(event_time,history_times,history_vectors[:,0]),
            np.interp(event_time,history_times,history_vectors[:,1]),
        ],dtype=float)

    def wind_field_at_time(model_time_s):
        if history_times is None or not spatial_wind_history:
            return np.broadcast_to(
                wind_at_time(model_time_s)[:,None,None],(2,ny,nx)).copy()
        event_time=(float(model_time_s)+float(wind_history_time_offset_s)
                    +(float(wind_history_advection_reference_x_m)-x)
                    /float(wind_history_advection_speed_m_s))
        columns=np.stack((
            np.interp(event_time,history_times,history_vectors[:,0]),
            np.interp(event_time,history_times,history_vectors[:,1])))
        return np.broadcast_to(columns[:,None,:],(2,ny,nx)).copy()

    wind=wind_at_time(0.0)
    reference_wind_speed=max(float(np.linalg.norm(reference_wind)),1.0e-30)
    # Receiver must preserve the dense-gas envelope: reject buoyant initial
    # mixtures instead of using gravity-current physics for H2 by accident.
    if np.any(q[1]>0) and thermo.fuel_molecular_weight<thermo._humid_ambient_molecular_weight:
        raise ValueError('single ground layer does not model positively buoyant species')
    xx,yy=np.meshgrid(x,y)
    cavity_streamwise_unit=np.array([1.0,0.0],dtype=float)
    if cavity_align_with_wind:
        cavity_reference_speed=float(np.linalg.norm(reference_wind))
        if cavity_reference_speed<=1.e-30:
            raise ValueError('wind-aligned cavity requires nonzero reference wind')
        cavity_streamwise_unit=reference_wind/cavity_reference_speed
    cavity_crosswind_unit=np.array(
        [-cavity_streamwise_unit[1],cavity_streamwise_unit[0]],dtype=float)
    cavity_streamwise_coordinate=np.zeros_like(xx)
    cavity_crosswind_coordinate=np.zeros_like(xx)
    if cavity_center_xy_m is not None:
        cx,cy=(float(v) for v in cavity_center_xy_m)
        cavity_dx=xx-cx
        cavity_dy=yy-cy
        cavity_streamwise_coordinate=(
            cavity_dx*cavity_streamwise_unit[0]
            +cavity_dy*cavity_streamwise_unit[1])
        cavity_crosswind_coordinate=(
            cavity_dx*cavity_crosswind_unit[0]
            +cavity_dy*cavity_crosswind_unit[1])
        circle=(cavity_dx**2+cavity_dy**2 <= float(cavity_radius_m)**2)
        if cavity_geometry=='downwind_semidisk':
            cavity_mask=circle & (cavity_streamwise_coordinate>=0.0)
        elif cavity_geometry=='upwind_semidisk':
            cavity_mask=circle & (cavity_streamwise_coordinate<=0.0)
        else:
            cavity_mask=circle
    else:
        cavity_mask=np.zeros((ny,nx),dtype=bool)
    wall_overflow_side_mask_x=np.zeros((ny,nx+1),dtype=bool)
    wall_overflow_side_mask_y=np.zeros((ny+1,nx),dtype=bool)
    wall_overflow_side_outward_sign_x=np.zeros((ny,nx+1),dtype=np.int8)
    wall_overflow_side_outward_sign_y=np.zeros((ny+1,nx),dtype=np.int8)
    if wall_overflow_side_gate_length_m>0.0:
        # Close the otherwise open, upwind half of the semicircular obstacle
        # with a diagnostic control surface.  Material leaving this surface
        # is the physical side-bypass branch: it escaped the cavity upstream
        # and can subsequently travel around either fence tip.  This replaces
        # the earlier pair of single-row tip gates, which missed most of that
        # curved path on a Cartesian grid.
        x_faces=float(origin_xy_m[0])+np.arange(nx+1,dtype=float)*h
        x_centres=float(origin_xy_m[0])+(np.arange(nx,dtype=float)+.5)*h
        upwind_x_face=((x_faces>=float(cx)-float(wall_overflow_side_gate_length_m)-.75*h)
                       &(x_faces<=float(cx)+.25*h))
        upwind_x_cell=((x_centres>=float(cx)-float(wall_overflow_side_gate_length_m)-.75*h)
                       &(x_centres<=float(cx)+.25*h))
        left_inside=circle[:,:-1]
        right_inside=circle[:,1:]
        crosses_x=left_inside^right_inside
        wall_overflow_side_mask_x[:,1:nx]=(
            crosses_x&upwind_x_face[None,1:nx]&(wx[:,1:nx]<=0.0))
        wall_overflow_side_outward_sign_x[:,1:nx]=np.where(
            left_inside&~right_inside,1,
            np.where(~left_inside&right_inside,-1,0)).astype(np.int8)
        below_inside=circle[:-1,:]
        above_inside=circle[1:,:]
        crosses_y=below_inside^above_inside
        wall_overflow_side_mask_y[1:ny,:]=(
            crosses_y&upwind_x_cell[None,:]&(wy[1:ny,:]<=0.0))
        wall_overflow_side_outward_sign_y[1:ny,:]=np.where(
            below_inside&~above_inside,1,
            np.where(~below_inside&above_inside,-1,0)).astype(np.int8)
    cavity_reservoir=np.zeros(4,dtype=float)
    local_cavity_reservoir=np.zeros_like(q)
    wall_impulse_wake_reservoir=np.zeros_like(q)
    # First temporal moment of species held by the independently transported
    # slow wake.  Its ratio to stored species is the local residence age used
    # by the normalized vertical observer.  It follows exactly the same
    # horizontal transport and first-order return as the wake inventory.
    wall_impulse_wake_species_age_moment=np.zeros_like(q[1])
    wall_impulse_wake_captured=np.zeros(4,dtype=float)
    wall_impulse_wake_returned=np.zeros(4,dtype=float)
    wall_impulse_wake_transport_moved=0.0
    wall_impulse_elevated_reservoir=np.zeros_like(q)
    wall_impulse_elevated_reservoir_age_moment=np.zeros_like(q[1])
    wall_impulse_elevated_reservoir_center_moment=np.zeros_like(q[1])
    wall_impulse_reattached_reservoir=np.zeros_like(q)
    wall_impulse_reattached_reservoir_age_moment=np.zeros_like(q[1])
    wall_impulse_reattached_captured=np.zeros(4,dtype=float)
    wall_impulse_elevated_captured=np.zeros(4,dtype=float)
    wall_impulse_elevated_returned=np.zeros(4,dtype=float)
    wall_impulse_elevated_reservoir_transport_moved=0.0
    wall_impulse_reattached_reservoir_transport_moved=0.0
    cell_wall_height=np.maximum.reduce((
        wx[:,:-1],wx[:,1:],wy[:-1,:],wy[1:,:]))
    # Directional face stores for material that has physically crossed the
    # finite-height part of an obstacle face.  Direction 0 is positive x/y;
    # direction 1 is negative x/y.  Values are conserved amounts (kg or
    # kg m/s), not per-area states.  A first-order return represents the
    # unresolved spread of overflow residence times without moving material
    # at wall-free faces or applying receptor-specific corrections.
    wall_overflow_x=np.zeros((2,4,ny,nx+1),dtype=float)
    wall_overflow_y=np.zeros((2,4,ny+1,nx),dtype=float)
    # Subset markers follow the same face-local storage path as the conserved
    # overflow packet.  These arrays contain absolute kg (or kg s for the age
    # moments), matching ``wall_overflow_x/y`` rather than cell-area density.
    # Keeping them direction resolved makes release to the receiver side
    # exact and lets the elevated observer coexist with delayed wall overflow.
    wall_overflow_elevated_x=np.zeros((2,ny,nx+1),dtype=float)
    wall_overflow_elevated_y=np.zeros((2,ny+1,nx),dtype=float)
    wall_overflow_elevated_age_x=np.zeros_like(wall_overflow_elevated_x)
    wall_overflow_elevated_age_y=np.zeros_like(wall_overflow_elevated_y)
    wall_overflow_flux_marker_x=np.zeros_like(wall_overflow_elevated_x)
    wall_overflow_flux_marker_y=np.zeros_like(wall_overflow_elevated_y)
    wall_overflow_flux_age_x=np.zeros_like(wall_overflow_elevated_x)
    wall_overflow_flux_age_y=np.zeros_like(wall_overflow_elevated_y)
    wall_overflow_captured=np.zeros(4,dtype=float)
    wall_overflow_returned=np.zeros(4,dtype=float)
    wall_overflow_top_captured=np.zeros(4,dtype=float)
    wall_overflow_side_captured=np.zeros(4,dtype=float)
    def wall_overflow_total():
        return (wall_overflow_x.sum(axis=(0,2,3))
                +wall_overflow_y.sum(axis=(0,2,3)))
    # Conserved cavity material can spend a declared residence delay in the
    # sub-grid storage before becoming eligible for the existing return
    # kernel.  Keeping it in a separate bucket makes the delay explicit while
    # preserving the same total mass/momentum ledger at every step.
    cavity_reservoir_pending=[]
    cavity_reservoir_sector=np.zeros((4,4),dtype=float)
    def cavity_reservoir_total():
        return (cavity_reservoir
                +local_cavity_reservoir.sum(axis=(1,2))*area
                +wall_impulse_wake_reservoir.sum(axis=(1,2))*area
                +wall_impulse_elevated_reservoir.sum(axis=(1,2))*area
                +wall_impulse_reattached_reservoir.sum(axis=(1,2))*area)
    # Optional finite-domain exterior recirculation.  A fraction of the
    # actually exported downstream packet is delayed and returned to the
    # registered cavity.  The packet carries all four conserved quantities;
    # this is a domain-reentry hypothesis, not a receptor correction.
    exterior_reentry_queue=[]
    exterior_reentry_returned=np.zeros(4,dtype=float)
    exterior_reentry_captured=np.zeros(4,dtype=float)
    exterior_reentry_momentum_impulse=np.zeros(2,dtype=float)
    source_injected=np.zeros(4,dtype=float)
    largest_transport_budget_detail=None
    cavity_sector_index=None
    if cavity_mask.any():
        theta=np.mod(np.arctan2(
            cavity_crosswind_coordinate[cavity_mask],
            cavity_streamwise_coordinate[cavity_mask])+2.0*np.pi,2.0*np.pi)
        cavity_sector_index=np.minimum((theta/(0.5*np.pi)).astype(int),3)
    reservoir_transport_mask=cavity_mask
    if cavity_reservoir_transport_mask=='wake':
        transport_direction=(reference_wind
                             /max(float(np.linalg.norm(reference_wind)),1.e-30))
        transport_dx=xx-cx;transport_dy=yy-cy
        transport_downstream=(transport_dx*transport_direction[0]
                              +transport_dy*transport_direction[1])
        transport_lateral=(-transport_dx*transport_direction[1]
                           +transport_dy*transport_direction[0])
        transport_half_width=(
            float(cavity_reservoir_wake_half_width_m)
            if cavity_reservoir_wake_half_width_m is not None
            else float(cavity_radius_m))
        wake_corridor=(
            (transport_downstream>=0.0)
            &(transport_downstream<=float(cavity_reservoir_wake_length_m))
            &(np.abs(transport_lateral)<=transport_half_width
              +float(wake_spread_slope)*transport_downstream))
        reservoir_transport_mask=cavity_mask|wake_corridor
    inside=(points[:,0]>=x[0]-h/2)&(points[:,0]<=x[-1]+h/2)&(points[:,1]>=y[0]-h/2)&(points[:,1]<=y[-1]+h/2)
    ix=np.clip(np.floor((points[:,0]-origin_xy_m[0])/h).astype(int),0,nx-1)
    iy=np.clip(np.floor((points[:,1]-origin_xy_m[1])/h).astype(int),0,ny-1)
    if receptor_interpolation=='bilinear':
        # Cell-centred interpolation is an observation operator only.  The
        # prognostic FV state and all conservation ledgers remain unchanged.
        gx=(points[:,0]-float(origin_xy_m[0]))/h-.5
        gy=(points[:,1]-float(origin_xy_m[1]))/h-.5
        ix0=np.clip(np.floor(gx).astype(int),0,max(nx-2,0))
        iy0=np.clip(np.floor(gy).astype(int),0,max(ny-2,0))
        tx=np.clip(gx-ix0,0.0,1.0)
        ty=np.clip(gy-iy0,0.0,1.0)
    else:
        ix0=iy0=tx=ty=None
    def sample_receptors(field):
        arr=np.asarray(field)
        if receptor_interpolation=='nearest':
            return arr[...,iy,ix]
        f00=arr[...,iy0,ix0]; f10=arr[...,iy0,ix0+1]
        f01=arr[...,iy0+1,ix0]; f11=arr[...,iy0+1,ix0+1]
        shape=(1,)*(arr.ndim-2)+(points.shape[0],)
        return (f00*((1.0-tx)*(1.0-ty)).reshape(shape)
                +f10*(tx*(1.0-ty)).reshape(shape)
                +f01*((1.0-tx)*ty).reshape(shape)
                +f11*(tx*ty).reshape(shape))
    vertical_q=None;vertical_dz=None
    elevated_species=(
        q[1].copy()*float(initial_elevated_species_fraction)
        if (wall_impulse_elevated_state_mode in {'species_marker','hybrid'}
            and (wall_impulse_elevated_transfer_fraction>0.0
            or wall_flux_elevated_transfer_fraction>0.0
            or initial_elevated_species_fraction>0.0)) else None
    )
    # First temporal moment of the independently transported elevated marker.
    # Dividing this conservative moment by ``elevated_species`` gives the
    # local time since wall transfer.  That age drives a distance-independent
    # settling/diffusion law in the observation operator, replacing the old
    # assumption that every downwind cell has one fixed vertical Gaussian.
    elevated_species_age_moment=np.zeros_like(q[1])
    # Provenance subset for species that crossed the geometrically open part
    # of a finite wall face.  It remains part of ``elevated_species`` and q;
    # this extra label lets the observer preserve the column's original
    # vertical kernel instead of assigning true over-wall flux the same
    # fixed Gaussian used for unresolved wall-impact lift.
    wall_flux_elevated_species=np.zeros_like(q[1])
    wall_flux_elevated_age_moment=np.zeros_like(q[1])
    local_cavity_reservoir_elevated_species=np.zeros_like(q[1])
    local_cavity_reservoir_elevated_age_moment=np.zeros_like(q[1])
    local_cavity_reservoir_wall_flux_species=np.zeros_like(q[1])
    local_cavity_reservoir_wall_flux_age_moment=np.zeros_like(q[1])
    cavity_reservoir_elevated_species=0.0
    cavity_reservoir_elevated_age_moment=0.0
    cavity_reservoir_wall_flux_species=0.0
    cavity_reservoir_wall_flux_age_moment=0.0
    elevated_species_transferred=0.0
    wall_flux_elevated_species_marked=0.0
    wall_flux_elevated_relative_transport=0.0
    elevated_species_relative_transport=0.0
    elevated_species_crosswind_diffusion=0.0
    elevated_transport_pending_s=0.0
    conserved_transport_pending_s=0.0
    reattached_transport_pending_s=0.0
    elevated_transport_update_count=0
    elevated_species_returned=0.0
    elevated_species_subset_correction=0.0
    if vertical_layers>1:
        vertical_top=max(float(np.max(points[:,2])),float(np.max(wx)),
                         float(np.max(wy)),1.0)
        if vertical_domain_top_m is not None:
            vertical_top=max(vertical_top,float(vertical_domain_top_m))
        elif initial_vertical_gaussian is not None:
            # Keep the declared source head inside the prognostic vertical
            # domain.  Previously a 0--13 m prism was projected into a grid
            # whose top was only the highest sensor (6.4 m), silently folding
            # upper mass into the lower layers before transport.
            centre_z,sigma_z=(float(v) for v in initial_vertical_gaussian)
            vertical_top=max(vertical_top,centre_z+4.0*sigma_z)
        vertical_dz=vertical_top/vertical_layers
        vertical_q=np.zeros((4,vertical_layers,ny,nx),dtype=float)
        initial_height,*_=_primitive(q,thermo)
        z_edges=np.arange(vertical_layers+1,dtype=float)*vertical_dz
        if initial_vertical_gaussian is None:
            overlap=np.maximum(np.minimum(initial_height[None,:,:],z_edges[1:,None,None])
                               -z_edges[:-1,None,None],0.0)
            fraction=np.divide(overlap,initial_height[None,:,:],
                               out=np.zeros_like(overlap),where=initial_height[None,:,:]>1e-30)
        else:
            centre_z,sigma_z=(float(v) for v in initial_vertical_gaussian)
            # A declared 3-D source distribution is defined on the fixed
            # vertical domain.  Do not truncate it by the shallower
            # hydrostatic height reconstructed from the column-mean q; that
            # would collapse a 0--13 m source into its ground layer before
            # the first transport step and make this option ineffective.
            upper=np.broadcast_to(
                z_edges[1:,None,None],(vertical_layers,ny,nx))
            lower=np.broadcast_to(
                z_edges[:-1,None,None],(vertical_layers,ny,nx))
            fraction=(ndtr((upper-centre_z)/sigma_z)-ndtr((lower-centre_z)/sigma_z)
                      +ndtr((upper+centre_z)/sigma_z)-ndtr((lower+centre_z)/sigma_z))
            fraction=np.maximum(fraction,0.0)
            total=fraction.sum(axis=0,keepdims=True)
            fraction=np.divide(fraction,total,out=np.zeros_like(fraction),where=total>1e-30)
        vertical_q=q[:,None,:,:]*fraction
    overflow_total=wall_overflow_total()
    reservoir_total=cavity_reservoir_total()
    times=[0.];mass=[float(q[0].sum()*area+reservoir_total[0]+overflow_total[0])];species=[float(q[1].sum()*area+reservoir_total[1]+overflow_total[1])]
    momentum_initial=(q[2:4].sum(axis=(1,2))*area+reservoir_total[2:4]
                      +overflow_total[2:4])
    maximum_momentum_budget_residual=0.0
    maximum_vertical_momentum_residual=0.0
    maximum_multilayer_horizontal_closure_residual=0.0
    maximum_multilayer_relative_velocity=0.0
    multilayer_relative_transport_moved=0.0
    multilayer_relative_transport_pending_s=0.0
    multilayer_relative_transport_update_count=0
    multilayer_relative_transport_substeps=0
    multilayer_species_bound_redistributed=0.0
    maximum_obstacle_vertical_velocity=0.0
    maximum_wall_impulse_vertical_velocity=0.0
    wall_impulse_lifted_species_kg=0.0
    height_resolved_wall_steps=0
    maximum_wall_bypass_mass_kg=0.0
    horizontal_diffusion_blocked_species_throughput_kg=0.0
    vertical_observation_maximum_relative_volume_residual=0.0
    vertical_observation_saturated_column_samples=0
    vertical_observation_component_traces={}
    traces=[];packets=[];heights=[];current_speeds=[];exported=0.;entrained=0.;mixture_exported=0.;maxres=0.;minq=float(q[:2].min())
    interface_exported=np.zeros(4,dtype=float)
    interface_injected=np.zeros(4,dtype=float)
    interface_negative_mass=0.0
    bx=np.concatenate((np.full(ny,x[0]-h/2),np.full(ny,x[-1]+h/2),x,x))
    by=np.concatenate((y,y,np.full(nx,y[0]-h/2),np.full(nx,y[-1]+h/2)))
    boundary_positions=np.column_stack((bx,by))
    _,observation_column_indices,observation_column_inverse=np.unique(
        points[:,:2],axis=0,return_index=True,return_inverse=True)
    observation_column_points=points[observation_column_indices]
    def observe_conservative(local_q):
        nonlocal vertical_observation_maximum_relative_volume_residual
        nonlocal vertical_observation_saturated_column_samples
        local_q=local_q[:,observation_column_indices]
        def sample_columns(field):
            return sample_receptors(field)[...,observation_column_indices]
        pure_density=float(thermo._density_from_temperature(thermo.ambient_temperature,1.))
        height,*_=_primitive(local_q,thermo)
        core_volume=local_q[1]/pure_density
        exponent=np.full_like(height,float(vertical_grounding_exponent))
        if vertical_observation_profile=='dense_power':
            exponent+=float(dense_vertical_exponent)
        cavity_weight=sample_columns(cavity_mask.astype(float))
        if cavity_vertical_exponent>0.:
            radial_weight=np.ones_like(height)
            if cavity_vertical_radial_decay_power>0.:
                radius=np.hypot(observation_column_points[:,0]-float(cx),observation_column_points[:,1]-float(cy))
                radial_weight=np.maximum(1.-radius/float(cavity_radius_m),0.)**float(cavity_vertical_radial_decay_power)
            exponent+=float(cavity_vertical_exponent)*cavity_weight*radial_weight
        development=np.ones_like(height)
        if wall_impulse_elevated_observation_ramp_m>0.:
            speed=max(float(np.linalg.norm(wind)),1.e-30)
            along=((observation_column_points[:,0]-float(cx))*wind[0]+(observation_column_points[:,1]-float(cy))*wind[1])/speed
            development=np.clip((along-float(cavity_radius_m))
                /float(wall_impulse_elevated_observation_ramp_m),0.,1.)
        branches=[];core_components=[]
        dynamic=np.zeros_like(height)
        if elevated_species is not None:
            marked=sample_columns(elevated_species)
            flux_marked=np.minimum(
                sample_columns(wall_flux_elevated_species),marked)
            impact_marked=np.maximum(marked-flux_marked,0.0)
            impact_age_moment=np.maximum(
                sample_columns(elevated_species_age_moment)
                -sample_columns(wall_flux_elevated_age_moment),0.0)
            age=np.divide(impact_age_moment,impact_marked,
                out=np.zeros_like(impact_marked),where=impact_marked>1.e-30)
            flux_share=np.clip(np.divide(
                flux_marked,local_q[1],out=np.zeros_like(flux_marked),
                where=local_q[1]>1.e-30),0.,1.)
            flux_volume=core_volume*flux_share
            if np.any(flux_volume>0.0):
                core_components.append(('wall_flux_core',flux_volume))
                core_volume=core_volume-flux_volume
            remaining_species=np.maximum(local_q[1]-flux_marked,0.0)
            dynamic=np.clip(np.divide(
                impact_marked,remaining_species,
                out=np.zeros_like(impact_marked),where=remaining_species>1.e-30),
                0.,1.)*development
            centre=np.maximum(float(wall_impulse_elevated_min_center_m),
                float(wall_impulse_elevated_center_m)
                +float(wall_impulse_elevated_rise_velocity_m_s)*float(wall_impulse_elevated_rise_time_s)
                *(-np.expm1(-age/float(wall_impulse_elevated_rise_time_s)))
                -float(wall_impulse_elevated_settling_velocity_m_s)*age)
            sigma=np.sqrt(float(wall_impulse_elevated_sigma_m)**2
                +2.*float(wall_impulse_elevated_vertical_diffusivity_m2_s)*age)
        else:
            centre=.65*height
            sigma=np.maximum(.22*height,1.e-20)
        share=1.-(1.-float(vertical_elevated_fraction))*(1.-dynamic)
        if np.any(share>0.):
            branches.append(('marker_upper',core_volume*share,centre,sigma))
            core_volume=core_volume*(1.-share)

        def add_reservoir(label,state,age_moment,centre,sigma,settling,kz,min_center,
                          visibility=1.,ramp=None,rise=0.,rise_time=1.,center_moment=None):
            nonlocal core_volume,height
            local=sample_columns(state)
            volume=local[1]/pure_density*visibility
            if not np.any(volume>0.):
                return
            age=np.divide(sample_columns(age_moment),local[1],
                out=np.zeros_like(volume),where=local[1]>1.e-30)
            if center_moment is not None:
                centre=np.divide(sample_columns(center_moment),local[1],
                    out=np.full_like(volume,centre),where=local[1]>1.e-30)
            centre=np.maximum(min_center,centre+rise*rise_time*(-np.expm1(-age/rise_time))-settling*age)
            sigma=np.sqrt(sigma*sigma+2.*kz*age)
            if ramp is not None:
                # Before the upper branch develops, its gas remains visible
                # in the lower column; an observation ramp cannot delete it.
                core_components.append((label+'_lower',volume*(1.-ramp)))
                extra_height,*_=_primitive(local,thermo)
                height+=extra_height*visibility*(1.-ramp)
                volume=volume*ramp
            branches.append((label+'_upper',volume,centre,sigma))

        if wall_impulse_elevated_state_mode in {'conserved','hybrid'} and wall_impulse_elevated_transfer_fraction>0.:
            add_reservoir('conserved_overflow',wall_impulse_elevated_reservoir,
                wall_impulse_elevated_reservoir_age_moment,conserved_prompt_center,
                conserved_prompt_sigma,conserved_prompt_settling,conserved_prompt_diffusivity,
                conserved_prompt_min_center,ramp=development,
                rise=float(wall_impulse_elevated_conserved_rise_velocity_m_s),
                rise_time=float(wall_impulse_elevated_conserved_rise_time_s),
                center_moment=(wall_impulse_elevated_reservoir_center_moment
                    if wall_impulse_elevated_conserved_center_mode=='impact_head' else None))
            if wall_impulse_reattached_fraction>0.:
                add_reservoir('reattached',wall_impulse_reattached_reservoir,
                    wall_impulse_reattached_reservoir_age_moment,float(wall_impulse_reattached_center_m),
                    float(wall_impulse_reattached_sigma_m),float(wall_impulse_reattached_settling_velocity_m_s),
                    float(wall_impulse_reattached_vertical_diffusivity_m2_s),
                    float(wall_impulse_reattached_min_center_m),ramp=development)
        if local_reservoir_observation_center_m is not None and local_reservoir_observation_fraction>0.:
            # This store has no age-based observation trajectory.
            add_reservoir('local_reservoir',local_cavity_reservoir,np.zeros_like(q[1]),
                float(local_reservoir_observation_center_m),float(local_reservoir_observation_sigma_m),
                0.,0.,0.,visibility=float(local_reservoir_observation_fraction))
        if wall_impulse_wake_observation_center_m is not None and wall_impulse_wake_observation_fraction>0.:
            add_reservoir('slow_wake',wall_impulse_wake_reservoir,wall_impulse_wake_species_age_moment,
                float(wall_impulse_wake_observation_center_m),float(wall_impulse_wake_observation_sigma_m),
                float(wall_impulse_wake_settling_velocity_m_s),float(wall_impulse_wake_vertical_diffusivity_m2_s),
                float(wall_impulse_wake_min_center_m),visibility=float(wall_impulse_wake_observation_fraction))
        tail_kind='none';tail_scale=1.
        if vertical_tail_extension_m>0.:
            if vertical_observation_profile=='exponential_tail':
                tail_kind='exponential';tail_scale=float(vertical_tail_scale_m)
            elif vertical_observation_profile=='diffusive_tail':
                tail_kind='diffusive'
                tail_scale=np.sqrt(4.*float(vertical_tail_diffusivity_m2_s)
                    *(1.+float(cavity_vertical_mixing_factor)*cavity_weight)*max(times[-1],1.e-12))
        profile=reconstruct_column(height_m=height,core_species_volume_m=core_volume,
            exponent=exponent,tail_kind=tail_kind,tail_scale_m=tail_scale,
            tail_extent_m=float(vertical_tail_extension_m),gaussian_branches=branches,
            core_components=core_components,
            bins=vertical_observation_bins)
        vertical_observation_maximum_relative_volume_residual=max(
            vertical_observation_maximum_relative_volume_residual,profile.maximum_relative_volume_residual)
        vertical_observation_saturated_column_samples+=profile.saturated_column_count
        values=np.where(inside,100.*profile.sample(points[:,2],observation_column_inverse),0.)
        if vertical_observation_component_diagnostics:
            attributed=100.*profile.sample_components(points[:,2],observation_column_inverse)
            attributed=np.where(inside[None,:],attributed,0.)
            current=dict(zip(profile.component_labels,attributed))
            prior_count=len(traces)
            for label in current:
                if label not in vertical_observation_component_traces:
                    vertical_observation_component_traces[label]=[
                        np.zeros(points.shape[0]) for _ in range(prior_count)]
            for label,history in vertical_observation_component_traces.items():
                history.append(current.get(label,np.zeros(points.shape[0])))
            if not np.allclose(sum(current.values()),values,rtol=1.e-12,atol=1.e-12):
                raise RuntimeError('vertical observation component attribution does not close')
        return values

    def observe():
        if vertical_q is not None:
            k=np.clip(np.floor(points[:,2]/vertical_dz).astype(int),0,vertical_layers-1)
            sampled_vertical=sample_receptors(vertical_q)
            receptor_index=np.arange(points.shape[0])
            local_m=sampled_vertical[0,k,receptor_index]/vertical_dz
            local_s=sampled_vertical[1,k,receptor_index]/vertical_dz
            local_fraction=np.divide(local_s,local_m,out=np.zeros_like(local_s),where=local_m>1e-30)
            mole=thermo.closure.mole_fraction_from_mass_fraction(np.clip(local_fraction,0.,1.))
            # A nonzero fixed-height layer is not necessarily filled by the
            # marked cloud. Recover its volume occupancy from the same binary
            # EOS used by the hydrostatic column. Without this factor, an
            # arbitrarily small diffusive tail is reported as a fully
            # occupied gas slab and grossly overstates concentration.
            marked_density=thermo._density_from_temperature(
                thermo.ambient_temperature,np.clip(local_fraction,0.,1.))
            occupancy=np.clip(np.divide(
                local_m,marked_density,out=np.zeros_like(local_m),
                where=marked_density>1e-30),0.0,1.0)
            values=np.where(inside&(local_m>1e-30),100*mole*occupancy,0.)
            traces.append(values)
            return
        observation_q=(q if local_reservoir_observation_center_m is not None
                       else q+float(local_reservoir_observation_fraction)
                       *local_cavity_reservoir)
        if wall_impulse_wake_observation_center_m is None:
            observation_q=(observation_q
                           +float(wall_impulse_wake_observation_fraction)
                           *wall_impulse_wake_reservoir)
        local_q=sample_receptors(observation_q)
        if vertical_observation_mode=='conservative_volume':
            traces.append(observe_conservative(local_q))
            return
        height,*_,fraction=_primitive(local_q,thermo)
        # Uniform hydrostatic layer; height below/above sensor is resolved
        # continuously. The optional exponential tail is a diagnostic
        # observation operator for unresolved turbulent vertical exchange; it
        # does not add mass to q and must not be interpreted as a conserved
        # 3-D plume closure.
        local_height=np.asarray(height).reshape(-1)
        local_fraction=np.asarray(fraction).reshape(-1)
        below=points[:,2]<=local_height+1.0e-12
        if cavity_vertical_exponent > 0.0:
            # A fence cavity can retain a ground-hugging recirculation layer
            # even when the outer patch is represented by a column mean.
            # This is an observation closure only: q remains unchanged and
            # the normalized shape integrates to the original column mass.
            in_cavity=sample_receptors(cavity_mask.astype(float))>0.0
            in_cavity &= below
            a=np.maximum(local_height,0.0)
            cavity_exponent=np.full_like(a,float(cavity_vertical_exponent))
            if cavity_vertical_radial_decay_power > 0.0:
                receptor_radius=np.sqrt(
                    (points[:,0]-float(cx))**2
                    +(points[:,1]-float(cy))**2
                )
                radial_weight=np.power(
                    np.clip(1.0-receptor_radius/max(float(cavity_radius_m),1.e-30),
                            0.0,1.0),
                    float(cavity_vertical_radial_decay_power),
                )
                cavity_exponent*=radial_weight
            shape=(cavity_exponent+1.0)*np.power(
                np.maximum(1.0-np.divide(points[:,2],np.maximum(a,1e-30)),0.0),
                cavity_exponent)
            local_fraction=np.where(in_cavity,local_fraction*shape,local_fraction)
        profile_base_fraction=local_fraction.copy()
        if vertical_grounding_exponent > 0.0:
            # Declared ground-current observation closure for the resolved
            # column.  It changes only the receptor projection, not q or its
            # mass ledger, and is intentionally separate from the cavity
            # wall-recirc shape.
            a=np.maximum(local_height,0.0)
            shape=(vertical_grounding_exponent+1.0)*np.power(
                np.maximum(1.0-np.divide(points[:,2],np.maximum(a,1e-30)),0.0),
                vertical_grounding_exponent)
            local_fraction=np.where(below,profile_base_fraction*shape,local_fraction)
        dynamic_elevated_fraction=np.zeros_like(local_fraction)
        if elevated_species is not None:
            local_elevated=sample_receptors(elevated_species)
            local_elevated_age=sample_receptors(elevated_species_age_moment)
            local_species=np.asarray(local_q[1]).reshape(-1)
            dynamic_elevated_fraction=np.clip(np.divide(
                local_elevated,local_species,
                out=np.zeros_like(local_elevated),where=local_species>1.e-30),0.0,1.0)
            if wall_impulse_elevated_observation_ramp_m>0.0:
                wind_norm=max(float(np.linalg.norm(wind)),1.e-30)
                streamwise=((points[:,0]-float(cx))*wind[0]
                            +(points[:,1]-float(cy))*wind[1])/wind_norm
                development=np.clip(
                    (streamwise-float(cavity_radius_m))
                    /float(wall_impulse_elevated_observation_ramp_m),0.0,1.0)
                dynamic_elevated_fraction*=development
        effective_elevated_fraction=(1.0-(1.0-float(vertical_elevated_fraction))
                                     *(1.0-dynamic_elevated_fraction))
        dynamic_elevated_visible=np.zeros_like(below,dtype=bool)
        elevated=np.zeros_like(local_fraction)
        if np.any(effective_elevated_fraction > 0.0):
            # Two-component unresolved vertical closure: retain the declared
            # ground-loaded profile while allowing a bounded overtopping mode
            # centred in the upper part of the resolved column.  Both shapes
            # integrate to one over 0<=z<=H, so this is an observation
            # projection only and cannot create or remove conserved mass.
            a=np.maximum(local_height,0.0)
            base_exp=float(vertical_grounding_exponent)
            ground_shape=(base_exp+1.0)*np.power(
                np.maximum(1.0-np.divide(points[:,2],np.maximum(a,1e-30)),0.0),
                base_exp)
            if elevated_species is None:
                centre=.65*a
                sigma=np.maximum(.22*a,1e-30)
            else:
                mean_age=np.divide(
                    local_elevated_age,local_elevated,
                    out=np.zeros_like(local_elevated_age),
                    where=local_elevated>1.e-30)
                centre=np.maximum(
                    float(wall_impulse_elevated_min_center_m),
                    float(wall_impulse_elevated_center_m)
                    +float(wall_impulse_elevated_rise_velocity_m_s)
                    *float(wall_impulse_elevated_rise_time_s)
                    *(1.0-np.exp(-mean_age/float(
                        wall_impulse_elevated_rise_time_s)))
                    -float(wall_impulse_elevated_settling_velocity_m_s)
                    *mean_age)
                sigma=np.sqrt(
                    float(wall_impulse_elevated_sigma_m)**2
                    +2.0*float(
                        wall_impulse_elevated_vertical_diffusivity_m2_s)
                    *mean_age)
                sigma=np.maximum(sigma,1.e-30)
            root=math.sqrt(2.0)*sigma
            upper=((a-centre)/root if elevated_species is None
                   else np.full_like(a,np.inf))
            lower=-centre/root
            integral=sigma*math.sqrt(math.pi/2.0)*(erf(upper)-erf(lower))
            elevated=np.divide(
                a*np.exp(-.5*((points[:,2]-centre)/sigma)**2),
                integral,out=np.zeros_like(a),where=integral>1e-30)
            shape=((1.0-effective_elevated_fraction)*ground_shape
                   +effective_elevated_fraction*elevated)
            local_fraction=np.where(below,profile_base_fraction*shape,local_fraction)
            if elevated_species is not None:
                above_elevated=~below
                local_fraction=np.where(
                    above_elevated,
                    profile_base_fraction*effective_elevated_fraction*elevated,
                    local_fraction)
                dynamic_elevated_visible=(
                    (effective_elevated_fraction>1.e-12)
                    &(points[:,2]<=centre+4.0*sigma))
        if vertical_observation_profile=='dense_power':
            # A dense-ground-current observation operator.  The conserved
            # layer still stores only the column mean, but the assumed
            # hydrostatic vertical shape is normalized analytically so its
            # integral over 0<=z<=H remains one.  This concentrates the
            # marked species near the ground without creating mass.
            a=np.maximum(local_height,0.0)
            exponent=float(dense_vertical_exponent)
            shape=(exponent+1.0)*np.power(
                np.maximum(1.0-np.divide(points[:,2],np.maximum(a,1e-30)),0.0),
                exponent)
            local_fraction=np.where(below,local_fraction*shape,local_fraction)
        elif vertical_observation_profile=='handoff_gaussian':
            if initial_vertical_gaussian is None:
                raise ValueError('handoff_gaussian requires initial_vertical_gaussian')
            centre_z,sigma_z=(float(v) for v in initial_vertical_gaussian)
            root=math.sqrt(2.0)*sigma_z
            shape=np.exp(-0.5*((points[:,2]-centre_z)/sigma_z)**2)
            shape+=np.exp(-0.5*((points[:,2]+centre_z)/sigma_z)**2)
            a=np.maximum(local_height,0.0)
            integral=sigma_z*math.sqrt(math.pi/2.0)*(
                erf((a-centre_z)/root)-erf(-centre_z/root)
                +erf((a+centre_z)/root)-erf(centre_z/root))
            factor=np.divide(a*shape,integral,out=np.zeros_like(a),where=integral>1e-30)
            local_fraction=np.where(below,local_fraction*factor,local_fraction)
        if (vertical_observation_profile in {'exponential_tail','diffusive_tail'}
                and vertical_tail_extension_m>0):
            above=points[:,2]>local_height
            tail=above & (points[:,2]<=local_height+vertical_tail_extension_m+1.0e-12)
            if vertical_observation_profile=='diffusive_tail':
                elapsed=max(times[-1],1.0e-12)
                # Obstacle-separated shear layers enhance vertical exchange
                # inside the registered cavity. This factor is geometry-wide
                # (never sensor-specific), while q remains conservative.
                cavity_factor=1.0 + float(cavity_vertical_mixing_factor) * (
                (sample_receptors(cavity_mask.astype(float))
                 if cavity_mask.size else 0.0)
                )
                penetration=np.sqrt(4.0*vertical_tail_diffusivity_m2_s
                                    *cavity_factor*elapsed)
                factor=erfc(np.divide(points[:,2]-local_height,
                                      np.maximum(penetration,1.0e-30)))
                if mass_conserving_vertical_profile:
                    # Normalize the resolved top-hat plus unresolved tail so
                    # its vertical integral remains the original layer mass.
                    # Integral_0^a erfc(u) du = a*erfc(a) +
                    # (1-exp(-a^2))/sqrt(pi).
                    ratio=np.divide(vertical_tail_extension_m,
                                    np.maximum(penetration,1.0e-30))
                    tail_integral=penetration*(
                        ratio*erfc(ratio)
                        +(1.0-np.exp(-ratio*ratio))/np.sqrt(np.pi))
                    norm=np.divide(local_height,
                                   local_height+tail_integral,
                                   out=np.ones_like(local_height),
                                   where=local_height>1.0e-30)
                    local_fraction=np.where(below,local_fraction*norm,local_fraction)
                    factor=factor*norm
            else:
                factor=np.exp(-(points[:,2]-local_height)/vertical_tail_scale_m)
                if mass_conserving_vertical_profile:
                    tail_integral=vertical_tail_scale_m*(
                        1.0-np.exp(-vertical_tail_extension_m/
                                   vertical_tail_scale_m))
                    norm=np.divide(local_height,
                                   local_height+tail_integral,
                                   out=np.ones_like(local_height),
                                   where=local_height>1.0e-30)
                    local_fraction=np.where(below,local_fraction*norm,local_fraction)
                    factor=factor*norm
            if elevated_species is not None:
                ground_tail=(profile_base_fraction
                             *(1.0-effective_elevated_fraction)*factor)
                elevated_tail=(profile_base_fraction
                               *effective_elevated_fraction*elevated)
                local_fraction=np.where(
                    tail,ground_tail+elevated_tail,local_fraction)
            else:
                local_fraction=np.where(tail,local_fraction*factor,local_fraction)
            visible=below|tail|dynamic_elevated_visible
        else:
            visible=below|dynamic_elevated_visible
        mole=thermo.closure.mole_fraction_from_mass_fraction(np.clip(local_fraction,0.,1.))
        values=np.where(inside&visible,100*mole,0.)
        if (wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                and wall_impulse_elevated_transfer_fraction>0.0):
            prompt_local=sample_receptors(wall_impulse_elevated_reservoir)
            prompt_mass=np.asarray(prompt_local[0]).reshape(-1)
            prompt_species=np.asarray(prompt_local[1]).reshape(-1)
            prompt_fraction=np.divide(
                prompt_species,prompt_mass,out=np.zeros_like(prompt_species),
                where=prompt_mass>1.e-30)
            prompt_density=thermo._density_from_temperature(
                thermo.ambient_temperature,np.clip(prompt_fraction,0.0,1.0))
            prompt_height=np.divide(
                prompt_mass,prompt_density,out=np.zeros_like(prompt_mass),
                where=prompt_density>1.e-30)
            prompt_mole=thermo.closure.mole_fraction_from_mass_fraction(
                np.clip(prompt_fraction,0.0,1.0))
            prompt_age_moment=np.asarray(sample_receptors(
                wall_impulse_elevated_reservoir_age_moment)).reshape(-1)
            prompt_mean_age=np.divide(
                prompt_age_moment,prompt_species,
                out=np.zeros_like(prompt_species),where=prompt_species>1.e-30)
            if wall_impulse_elevated_conserved_center_mode=='impact_head':
                prompt_center_moment=np.asarray(sample_receptors(
                    wall_impulse_elevated_reservoir_center_moment)).reshape(-1)
                prompt_initial_center=np.divide(
                    prompt_center_moment,prompt_species,
                    out=np.full_like(prompt_species,conserved_prompt_center),
                    where=prompt_species>1.e-30)
            else:
                prompt_initial_center=conserved_prompt_center
            prompt_centre=np.maximum(
                conserved_prompt_min_center,
                prompt_initial_center
                +float(wall_impulse_elevated_conserved_rise_velocity_m_s)
                *float(wall_impulse_elevated_conserved_rise_time_s)
                *(1.0-np.exp(-prompt_mean_age/float(
                    wall_impulse_elevated_conserved_rise_time_s)))
                -conserved_prompt_settling
                *prompt_mean_age)
            prompt_sigma=np.sqrt(
                conserved_prompt_sigma**2
                +2.0*conserved_prompt_diffusivity
                *prompt_mean_age)
            prompt_sigma=np.maximum(prompt_sigma,1.e-30)
            prompt_root=math.sqrt(2.0)*prompt_sigma
            prompt_norm=(prompt_sigma*math.sqrt(math.pi/2.0)
                         *(1.0+erf(prompt_centre/prompt_root)))
            prompt_profile=(np.exp(
                -0.5*((points[:,2]-prompt_centre)/prompt_sigma)**2)
                /prompt_norm)
            prompt_development=np.ones_like(prompt_species)
            if wall_impulse_elevated_observation_ramp_m>0.0:
                wind_norm=max(float(np.linalg.norm(wind)),1.e-30)
                streamwise=((points[:,0]-float(cx))*wind[0]
                            +(points[:,1]-float(cy))*wind[1])/wind_norm
                prompt_development=np.clip(
                    (streamwise-float(cavity_radius_m))
                    /float(wall_impulse_elevated_observation_ramp_m),0.0,1.0)
            prompt_vol_percent=(100.0*prompt_mole*prompt_height
                                *prompt_profile*prompt_development)
            values=np.where(
                inside,np.clip(values+prompt_vol_percent,0.0,100.0),0.0)
            if wall_impulse_reattached_fraction>0.0:
                lower_local=sample_receptors(wall_impulse_reattached_reservoir)
                lower_mass=np.asarray(lower_local[0]).reshape(-1)
                lower_species=np.asarray(lower_local[1]).reshape(-1)
                lower_fraction=np.divide(
                    lower_species,lower_mass,out=np.zeros_like(lower_species),
                    where=lower_mass>1.e-30)
                lower_density=thermo._density_from_temperature(
                    thermo.ambient_temperature,np.clip(lower_fraction,0.0,1.0))
                lower_height=np.divide(
                    lower_mass,lower_density,out=np.zeros_like(lower_mass),
                    where=lower_density>1.e-30)
                lower_mole=thermo.closure.mole_fraction_from_mass_fraction(
                    np.clip(lower_fraction,0.0,1.0))
                lower_age_moment=np.asarray(sample_receptors(
                    wall_impulse_reattached_reservoir_age_moment)).reshape(-1)
                lower_mean_age=np.divide(
                    lower_age_moment,lower_species,
                    out=np.zeros_like(lower_species),where=lower_species>1.e-30)
                lower_centre=np.maximum(
                    float(wall_impulse_reattached_min_center_m),
                    float(wall_impulse_reattached_center_m)
                    -float(wall_impulse_reattached_settling_velocity_m_s)
                    *lower_mean_age)
                lower_sigma=np.sqrt(
                    float(wall_impulse_reattached_sigma_m)**2
                    +2.0*float(
                        wall_impulse_reattached_vertical_diffusivity_m2_s)
                    *lower_mean_age)
                lower_sigma=np.maximum(lower_sigma,1.e-30)
                lower_root=math.sqrt(2.0)*lower_sigma
                lower_norm=(lower_sigma*math.sqrt(math.pi/2.0)
                            *(1.0+erf(lower_centre/lower_root)))
                lower_profile=(np.exp(
                    -0.5*((points[:,2]-lower_centre)/lower_sigma)**2)
                    /lower_norm)
                lower_vol_percent=(100.0*lower_mole*lower_height
                                   *lower_profile*prompt_development)
                values=np.where(
                    inside,np.clip(values+lower_vol_percent,0.0,100.0),0.0)
        if (local_reservoir_observation_center_m is not None
                and local_reservoir_observation_fraction>0.0):
            stored_species=sample_receptors(local_cavity_reservoir[1])
            pure_density=float(thermo._density_from_temperature(
                thermo.ambient_temperature,1.0))
            equivalent_height=stored_species/max(pure_density,1.e-30)
            centre=float(local_reservoir_observation_center_m)
            sigma=float(local_reservoir_observation_sigma_m)
            root=math.sqrt(2.0)*sigma
            normalizer=(sigma*math.sqrt(math.pi/2.0)
                        *(1.0+math.erf(centre/root)))
            profile=np.exp(-0.5*((points[:,2]-centre)/sigma)**2)/normalizer
            reservoir_vol_percent=(100.0
                *float(local_reservoir_observation_fraction)
                *equivalent_height*profile)
            values=np.where(inside,np.clip(values+reservoir_vol_percent,0.0,100.0),0.0)
        if (wall_impulse_wake_observation_center_m is not None
                and wall_impulse_wake_observation_fraction>0.0):
            wake_local=sample_receptors(wall_impulse_wake_reservoir)
            wake_mass=np.asarray(wake_local[0]).reshape(-1)
            wake_species=np.asarray(wake_local[1]).reshape(-1)
            wake_age_moment=np.asarray(sample_receptors(
                wall_impulse_wake_species_age_moment)).reshape(-1)
            wake_mean_age=np.divide(
                wake_age_moment,wake_species,
                out=np.zeros_like(wake_species),where=wake_species>1.e-30)
            wake_fraction=np.divide(
                wake_species,wake_mass,out=np.zeros_like(wake_species),
                where=wake_mass>1.e-30)
            wake_density=thermo._density_from_temperature(
                thermo.ambient_temperature,np.clip(wake_fraction,0.0,1.0))
            wake_height=np.divide(
                wake_mass,wake_density,out=np.zeros_like(wake_mass),
                where=wake_density>1.e-30)
            wake_mole=thermo.closure.mole_fraction_from_mass_fraction(
                np.clip(wake_fraction,0.0,1.0))
            centre=np.maximum(
                float(wall_impulse_wake_min_center_m),
                float(wall_impulse_wake_observation_center_m)
                -float(wall_impulse_wake_settling_velocity_m_s)
                *wake_mean_age)
            sigma=np.sqrt(
                float(wall_impulse_wake_observation_sigma_m)**2
                +2.0*float(wall_impulse_wake_vertical_diffusivity_m2_s)
                *wake_mean_age)
            sigma=np.maximum(sigma,1.e-30)
            root=math.sqrt(2.0)*sigma
            normalizer=(sigma*math.sqrt(math.pi/2.0)
                        *(1.0+erf(centre/root)))
            profile=np.exp(-0.5*((points[:,2]-centre)/sigma)**2)/normalizer
            wake_vol_percent=(100.0
                *float(wall_impulse_wake_observation_fraction)
                *wake_mole*wake_height*profile)
            values=np.where(
                inside,np.clip(values+wake_vol_percent,0.0,100.0),0.0)
        traces.append(values)
    def split_multilayer_face_flux(face_flux,component,axis,wall_height):
        """Allocate a solved column flux to donor-side fixed vertical layers.

        The aggregate Rusanov solve remains authoritative.  This operator
        transports its mixture/species flux with the donor column's vertical
        distribution.  At a finite wall only the geometrical part of a layer
        above the wall is eligible, and eligible weights are renormalized to
        the already-computed aggregate overtopping flux.
        """
        layer_state=vertical_q[component]
        fraction=np.divide(
            layer_state,q[component][None,:,:],out=np.zeros_like(layer_state),
            where=q[component][None,:,:]>1.e-30)
        if axis==0:
            dry=np.zeros((vertical_layers,ny,1),dtype=float)
            left=np.concatenate((dry,fraction),axis=2)
            right=np.concatenate((fraction,dry),axis=2)
        else:
            dry=np.zeros((vertical_layers,1,nx),dtype=float)
            left=np.concatenate((dry,fraction),axis=1)
            right=np.concatenate((fraction,dry),axis=1)
        donor=np.where(face_flux[None,:,:]>=0.0,left,right)
        if height_resolved_wall_flux:
            layer_above=np.clip(
                (z_edges[1:,None,None]
                 -np.maximum(z_edges[:-1,None,None],wall_height[None,:,:]))
                /vertical_dz,0.0,1.0)
            eligible=np.where(wall_height[None,:,:]>0.0,layer_above,1.0)
            weights=donor*eligible
        else:
            weights=donor
        # Do not renormalize the above-wall weights.  For a wall face the
        # caller supplies the unrestricted scalar flux, so their sum is the
        # donor-specific open fraction.  Renormalizing here would force the
        # full column flux through a small elevated inventory and can remove
        # more mass from that layer than it contains.
        return face_flux[None,:,:]*weights
    def split_elevated_species_flux(face_flux,axis,field=None):
        if field is None:
            field=elevated_species
        fraction=np.divide(
            field,q[1],out=np.zeros_like(field),
            where=q[1]>1.e-30)
        if axis==0:
            left=np.concatenate((np.zeros((ny,1)),fraction),axis=1)
            right=np.concatenate((fraction,np.zeros((ny,1))),axis=1)
        else:
            left=np.concatenate((np.zeros((1,nx)),fraction),axis=0)
            right=np.concatenate((fraction,np.zeros((1,nx))),axis=0)
        donor=np.where(face_flux>=0.0,left,right)
        return face_flux*donor
    def transport_elevated_relative(relative_u,relative_v,step_dt):
        """Conservatively slip the wall-elevated marker through the carrier.

        The ground branch has already followed the hydrostatic layer flux.
        This bounded correction moves only the elevated marked-gas subset at
        a fraction of its ambient-minus-layer slip velocity.  Face transfers
        are limited first by donor inventory and then by receiver mixture
        capacity, preserving both global species mass and ``0<=q1<=q0``.
        """
        nonlocal elevated_species,elevated_species_age_moment
        nonlocal wall_flux_elevated_species,wall_flux_elevated_age_moment
        # The wall-flux provenance subset has its own upper-layer carrier.
        # Move only the unresolved impact marker here so the two physical
        # paths do not inherit the same empirical crosswind velocity.
        e=np.maximum(elevated_species-wall_flux_elevated_species,0.0)
        age_density=np.maximum(
            elevated_species_age_moment-wall_flux_elevated_age_moment,0.0)
        age_per_species=np.divide(
            age_density,e,out=np.zeros_like(age_density),where=e>1.e-30)
        flux_marker=np.zeros_like(e)
        flux_fraction=np.divide(
            flux_marker,e,out=np.zeros_like(flux_marker),where=e>1.e-30)
        flux_age_per_species=np.divide(
            wall_flux_elevated_age_moment,flux_marker,
            out=np.zeros_like(flux_marker),where=flux_marker>1.e-30)
        ux=0.5*(relative_u[:,:-1]+relative_u[:,1:])
        vy=0.5*(relative_v[:-1,:]+relative_v[1:,:])
        wall_x=wx[:,1:nx]
        wall_y=wy[1:ny,:]
        root=math.sqrt(2.0)*float(wall_impulse_elevated_sigma_m)
        open_x=np.where(
            wall_x>0.0,
            0.5*erfc((wall_x-float(wall_impulse_elevated_center_m))/root),
            1.0)
        open_y=np.where(
            wall_y>0.0,
            0.5*erfc((wall_y-float(wall_impulse_elevated_center_m))/root),
            1.0)
        scale=float(step_dt)/h
        xp=np.maximum(ux,0.0)*scale*e[:,:-1]*open_x
        xn=np.maximum(-ux,0.0)*scale*e[:,1:]*open_x
        yp=np.maximum(vy,0.0)*scale*e[:-1,:]*open_y
        yn=np.maximum(-vy,0.0)*scale*e[1:,:]*open_y
        outgoing=np.zeros_like(e)
        outgoing[:,:-1]+=xp;outgoing[:,1:]+=xn
        outgoing[:-1,:]+=yp;outgoing[1:,:]+=yn
        donor_scale=np.minimum(
            1.0,np.divide(e,outgoing,out=np.ones_like(e),where=outgoing>1.e-30))
        xp*=donor_scale[:,:-1];xn*=donor_scale[:,1:]
        yp*=donor_scale[:-1,:];yn*=donor_scale[1:,:]
        outgoing.fill(0.0)
        outgoing[:,:-1]+=xp;outgoing[:,1:]+=xn
        outgoing[:-1,:]+=yp;outgoing[1:,:]+=yn
        incoming=np.zeros_like(e)
        incoming[:,1:]+=xp;incoming[:,:-1]+=xn
        incoming[1:,:]+=yp;incoming[:-1,:]+=yn
        capacity=np.maximum(q[0]-q[1]+outgoing,0.0)
        receiver_scale=np.minimum(
            1.0,np.divide(capacity,incoming,out=np.ones_like(e),where=incoming>1.e-30))
        xp*=receiver_scale[:,1:];xn*=receiver_scale[:,:-1]
        yp*=receiver_scale[1:,:];yn*=receiver_scale[:-1,:]
        age_xp=xp*age_per_species[:,:-1]
        age_xn=xn*age_per_species[:,1:]
        age_yp=yp*age_per_species[:-1,:]
        age_yn=yn*age_per_species[1:,:]
        flux_xp=xp*flux_fraction[:,:-1]
        flux_xn=xn*flux_fraction[:,1:]
        flux_yp=yp*flux_fraction[:-1,:]
        flux_yn=yn*flux_fraction[1:,:]
        flux_age_xp=flux_xp*flux_age_per_species[:,:-1]
        flux_age_xn=flux_xn*flux_age_per_species[:,1:]
        flux_age_yp=flux_yp*flux_age_per_species[:-1,:]
        flux_age_yn=flux_yn*flux_age_per_species[1:,:]
        delta=np.zeros_like(e)
        delta[:,:-1]-=xp;delta[:,1:]+=xp
        delta[:,1:]-=xn;delta[:,:-1]+=xn
        delta[:-1,:]-=yp;delta[1:,:]+=yp
        delta[1:,:]-=yn;delta[:-1,:]+=yn
        age_delta=np.zeros_like(age_density)
        age_delta[:,:-1]-=age_xp;age_delta[:,1:]+=age_xp
        age_delta[:,1:]-=age_xn;age_delta[:,:-1]+=age_xn
        age_delta[:-1,:]-=age_yp;age_delta[1:,:]+=age_yp
        age_delta[1:,:]-=age_yn;age_delta[:-1,:]+=age_yn
        flux_delta=np.zeros_like(flux_marker)
        flux_delta[:,:-1]-=flux_xp;flux_delta[:,1:]+=flux_xp
        flux_delta[:,1:]-=flux_xn;flux_delta[:,:-1]+=flux_xn
        flux_delta[:-1,:]-=flux_yp;flux_delta[1:,:]+=flux_yp
        flux_delta[1:,:]-=flux_yn;flux_delta[:-1,:]+=flux_yn
        flux_age_delta=np.zeros_like(flux_marker)
        flux_age_delta[:,:-1]-=flux_age_xp;flux_age_delta[:,1:]+=flux_age_xp
        flux_age_delta[:,1:]-=flux_age_xn;flux_age_delta[:,:-1]+=flux_age_xn
        flux_age_delta[:-1,:]-=flux_age_yp;flux_age_delta[1:,:]+=flux_age_yp
        flux_age_delta[1:,:]-=flux_age_yn;flux_age_delta[:-1,:]+=flux_age_yn
        if abs(float(delta.sum()))>1.e-10:
            raise RuntimeError('elevated relative transport lost species mass')
        if abs(float(age_delta.sum()))>1.e-8:
            raise RuntimeError('elevated relative transport lost age moment')
        elevated_species+=delta
        elevated_species_age_moment+=age_delta
        wall_flux_elevated_species+=flux_delta
        wall_flux_elevated_age_moment+=flux_age_delta
        q[1]+=delta
        if (np.min(elevated_species)<-1.e-10
                or np.max(elevated_species-q[1])>1.e-10
                or np.min(wall_flux_elevated_species)<-1.e-10
                or np.max(wall_flux_elevated_species-elevated_species)>1.e-10
                or np.min(q[1])<-1.e-10
                or np.max(q[1]-q[0])>1.e-10):
            raise RuntimeError('elevated relative transport violated subset bounds')
        return float((xp.sum()+xn.sum()+yp.sum()+yn.sum())*area)
    def diffuse_elevated_crosswind(step_dt,carrier_wind):
        """Conservatively disperse only the wall-elevated marked subset.

        This is a sub-grid closure for the two lateral shear layers shed by
        the obstacle.  It deliberately leaves the ground-loaded carrier
        untouched.  Pairwise face transfers, donor limiting, receiver
        capacity limiting and zero boundary flux preserve marked-species
        mass and ``0 <= elevated_species <= q[1] <= q[0]``.
        """
        nonlocal elevated_species,elevated_species_age_moment,q
        nonlocal wall_flux_elevated_species,wall_flux_elevated_age_moment
        diffusivity=float(wall_impulse_elevated_crosswind_diffusivity_m2_s)
        if diffusivity<=0.0 or step_dt<=0.0:
            return 0.0
        substeps=max(1,int(math.ceil(diffusivity*float(step_dt)/(.24*h*h))))
        local_dt=float(step_dt)/substeps
        moved_total=0.0
        for _ in range(substeps):
            e=np.maximum(
                elevated_species-wall_flux_elevated_species,0.0)
            age_density=np.maximum(
                elevated_species_age_moment-wall_flux_elevated_age_moment,0.0)
            age_per_species=np.divide(
                age_density,e,out=np.zeros_like(age_density),where=e>1.e-30)
            flux_marker=np.zeros_like(e)
            flux_fraction=np.divide(
                flux_marker,e,out=np.zeros_like(flux_marker),where=e>1.e-30)
            flux_age_per_species=np.divide(
                wall_flux_elevated_age_moment,flux_marker,
                out=np.zeros_like(flux_marker),where=flux_marker>1.e-30)
            speed=np.hypot(carrier_wind[0],carrier_wind[1])
            cross_x=np.divide(
                -carrier_wind[1],speed,out=np.zeros_like(speed),
                where=speed>1.e-30)
            cross_y=np.divide(
                carrier_wind[0],speed,out=np.zeros_like(speed),
                where=speed>1.e-30)
            dxx=diffusivity*cross_x*cross_x
            dyy=diffusivity*cross_y*cross_y
            root=math.sqrt(2.0)*float(wall_impulse_elevated_sigma_m)
            open_x=np.where(
                wx[:,1:nx]>0.0,
                0.5*erfc((wx[:,1:nx]-float(
                    wall_impulse_elevated_center_m))/root),1.0)
            open_y=np.where(
                wy[1:ny,:]>0.0,
                0.5*erfc((wy[1:ny,:]-float(
                    wall_impulse_elevated_center_m))/root),1.0)
            flux_x=(local_dt/(h*h)*0.5*(dxx[:,:-1]+dxx[:,1:])
                    *(e[:,:-1]-e[:,1:])*open_x)
            flux_y=(local_dt/(h*h)*0.5*(dyy[:-1,:]+dyy[1:,:])
                    *(e[:-1,:]-e[1:,:])*open_y)
            xp=np.maximum(flux_x,0.0);xn=np.maximum(-flux_x,0.0)
            yp=np.maximum(flux_y,0.0);yn=np.maximum(-flux_y,0.0)
            outgoing=np.zeros_like(e)
            outgoing[:,:-1]+=xp;outgoing[:,1:]+=xn
            outgoing[:-1,:]+=yp;outgoing[1:,:]+=yn
            donor_scale=np.minimum(
                1.0,np.divide(e,outgoing,out=np.ones_like(e),
                              where=outgoing>1.e-30))
            xp*=donor_scale[:,:-1];xn*=donor_scale[:,1:]
            yp*=donor_scale[:-1,:];yn*=donor_scale[1:,:]
            outgoing.fill(0.0)
            outgoing[:,:-1]+=xp;outgoing[:,1:]+=xn
            outgoing[:-1,:]+=yp;outgoing[1:,:]+=yn
            incoming=np.zeros_like(e)
            incoming[:,1:]+=xp;incoming[:,:-1]+=xn
            incoming[1:,:]+=yp;incoming[:-1,:]+=yn
            capacity=np.maximum(q[0]-q[1]+outgoing,0.0)
            receiver_scale=np.minimum(
                1.0,np.divide(capacity,incoming,out=np.ones_like(e),
                              where=incoming>1.e-30))
            xp*=receiver_scale[:,1:];xn*=receiver_scale[:,:-1]
            yp*=receiver_scale[1:,:];yn*=receiver_scale[:-1,:]
            age_xp=xp*age_per_species[:,:-1]
            age_xn=xn*age_per_species[:,1:]
            age_yp=yp*age_per_species[:-1,:]
            age_yn=yn*age_per_species[1:,:]
            flux_xp=xp*flux_fraction[:,:-1]
            flux_xn=xn*flux_fraction[:,1:]
            flux_yp=yp*flux_fraction[:-1,:]
            flux_yn=yn*flux_fraction[1:,:]
            flux_age_xp=flux_xp*flux_age_per_species[:,:-1]
            flux_age_xn=flux_xn*flux_age_per_species[:,1:]
            flux_age_yp=flux_yp*flux_age_per_species[:-1,:]
            flux_age_yn=flux_yn*flux_age_per_species[1:,:]
            delta=np.zeros_like(e)
            delta[:,:-1]-=xp;delta[:,1:]+=xp
            delta[:,1:]-=xn;delta[:,:-1]+=xn
            delta[:-1,:]-=yp;delta[1:,:]+=yp
            delta[1:,:]-=yn;delta[:-1,:]+=yn
            age_delta=np.zeros_like(age_density)
            age_delta[:,:-1]-=age_xp;age_delta[:,1:]+=age_xp
            age_delta[:,1:]-=age_xn;age_delta[:,:-1]+=age_xn
            age_delta[:-1,:]-=age_yp;age_delta[1:,:]+=age_yp
            age_delta[1:,:]-=age_yn;age_delta[:-1,:]+=age_yn
            flux_delta=np.zeros_like(flux_marker)
            flux_delta[:,:-1]-=flux_xp;flux_delta[:,1:]+=flux_xp
            flux_delta[:,1:]-=flux_xn;flux_delta[:,:-1]+=flux_xn
            flux_delta[:-1,:]-=flux_yp;flux_delta[1:,:]+=flux_yp
            flux_delta[1:,:]-=flux_yn;flux_delta[:-1,:]+=flux_yn
            flux_age_delta=np.zeros_like(flux_marker)
            flux_age_delta[:,:-1]-=flux_age_xp;flux_age_delta[:,1:]+=flux_age_xp
            flux_age_delta[:,1:]-=flux_age_xn;flux_age_delta[:,:-1]+=flux_age_xn
            flux_age_delta[:-1,:]-=flux_age_yp;flux_age_delta[1:,:]+=flux_age_yp
            flux_age_delta[1:,:]-=flux_age_yn;flux_age_delta[:-1,:]+=flux_age_yn
            if abs(float(delta.sum()))>1.e-10:
                raise RuntimeError('elevated crosswind diffusion lost species mass')
            if abs(float(age_delta.sum()))>1.e-8:
                raise RuntimeError('elevated crosswind diffusion lost age moment')
            elevated_species+=delta
            elevated_species_age_moment+=age_delta
            wall_flux_elevated_species+=flux_delta
            wall_flux_elevated_age_moment+=flux_age_delta
            q[1]+=delta
            if (np.min(elevated_species)<-1.e-10
                    or np.max(elevated_species-q[1])>1.e-10
                    or np.min(wall_flux_elevated_species)<-1.e-10
                    or np.max(wall_flux_elevated_species-elevated_species)>1.e-10
                    or np.min(q[1])<-1.e-10
                    or np.max(q[1]-q[0])>1.e-10):
                raise RuntimeError(
                    'elevated crosswind diffusion violated subset bounds')
            moved_total+=float(
                (xp.sum()+xn.sum()+yp.sum()+yn.sum())*area)
        return moved_total
    def transport_wall_flux_relative(relative_u,relative_v,step_dt):
        """Move true over-wall species with its height-resolved carrier.

        ``wall_flux_elevated_species`` is a provenance subset of the ordinary
        elevated marker.  Reusing the bounded marker transport in an isolated
        swap keeps its species and age moments conservative while preventing
        the wall-impact marker's empirical crosswind drift from being applied
        to material whose velocity follows the measured upper-level wind.
        """
        nonlocal elevated_species,elevated_species_age_moment
        nonlocal wall_flux_elevated_species,wall_flux_elevated_age_moment
        if step_dt<=0.0 or not np.any(wall_flux_elevated_species>0.0):
            return 0.0
        impact_species=np.maximum(
            elevated_species-wall_flux_elevated_species,0.0)
        impact_age=np.maximum(
            elevated_species_age_moment-wall_flux_elevated_age_moment,0.0)
        flux_species=wall_flux_elevated_species.copy()
        flux_age=wall_flux_elevated_age_moment.copy()
        elevated_species=flux_species
        elevated_species_age_moment=flux_age
        wall_flux_elevated_species=np.zeros_like(flux_species)
        wall_flux_elevated_age_moment=np.zeros_like(flux_age)
        moved=transport_elevated_relative(relative_u,relative_v,step_dt)
        flux_species=elevated_species
        flux_age=elevated_species_age_moment
        wall_flux_elevated_species=flux_species
        wall_flux_elevated_age_moment=flux_age
        elevated_species=impact_species+flux_species
        elevated_species_age_moment=impact_age+flux_age
        if np.max(wall_flux_elevated_species-elevated_species)>1.e-10:
            raise RuntimeError('wall-flux marker ceased to be an elevated subset')
        return moved
    observe()
    while times[-1]<duration_s-1e-12:
        pending_total=(np.sum(np.stack([packet for _,packet in cavity_reservoir_pending],axis=0),axis=0)
                       if cavity_reservoir_pending else np.zeros(4))
        overflow_total=wall_overflow_total()
        momentum_before_transport=(q[2:4].sum(axis=(1,2))*area
                                   +cavity_reservoir_total()[2:4]+pending_total[2:4]
                                   +overflow_total[2:4])
        height,u,v,c,_,_=_primitive(q,thermo,gravity_scale)
        rate=(float(np.max(np.abs(u)+c))+float(np.max(np.abs(v)+c)))/h
        dt=min(max_dt_s,.42/max(rate,1e-30),duration_s-times[-1])
        # Use the midpoint ambient vector for every forcing term in this
        # interval.  The resolved conserved state still determines the FV
        # wave speed and face flux; this history changes only the measured
        # ambient carrier, entrainment, drag and declared obstacle closures.
        wind=wind_at_time(times[-1]+0.5*dt)
        if vertical_q is not None:
            if vertical_mixing_diffusivity_m2_s>0.0:
                dt=min(dt,.45*vertical_dz*vertical_dz/vertical_mixing_diffusivity_m2_s)
            vertical_speed_bound=(obstacle_vertical_lift_fraction
                                  *float(np.linalg.norm(wind))
                                  +dense_settling_velocity_m_s)
            if wall_impulse_vertical_lift_fraction>0.0:
                vertical_speed_bound+=float(wall_impulse_vertical_lift_fraction)*max(
                    float(np.max(np.abs(u))),float(np.max(np.abs(v))))
            if vertical_speed_bound>0.0:
                dt=min(dt,.45*vertical_dz/vertical_speed_bound)
        openness_x = openness_y = None
        if height_resolved_wall_flux and vertical_q is not None:
            # Resolve the fraction of each wall face occupied by gas above the
            # registered wall height.  This is a conservative flux modifier:
            # it only changes which already-conserved vertical mass crosses a
            # face; it cannot create species or carrier mass.  At a boundary
            # with no wall the fraction is exactly one.
            z_edges=np.arange(vertical_layers+1,dtype=float)*vertical_dz
            def face_openness(mass_left,mass_right,wall):
                wall=np.asarray(wall,float)
                layer_above=np.clip(
                    (z_edges[1:,None,None]-np.maximum(z_edges[:-1,None,None],wall[None,:,:]))
                    /vertical_dz,0.0,1.0)
                total_left=mass_left.sum(axis=0)
                total_right=mass_right.sum(axis=0)
                above_left=(mass_left*layer_above).sum(axis=0)
                above_right=(mass_right*layer_above).sum(axis=0)
                total=np.maximum(total_left,total_right)
                above=np.maximum(above_left,above_right)
                result=np.divide(above,total,out=np.ones_like(above),where=total>1e-30)
                return np.clip(result,0.0,1.0)
            # Pad with dry exterior columns.  A wall-free outer face remains
            # open; a finite-height outer wall is closed unless the adjacent
            # column has mass above it.
            zero=np.zeros((vertical_layers,ny,1),dtype=float)
            # The staggered strips leave the exterior face dry and pair each
            # interior face with its two adjacent cell columns.
            mass_x_l=np.concatenate((zero,vertical_q[0]),axis=2)
            mass_x_r=np.concatenate((vertical_q[0],zero),axis=2)
            openness_x=face_openness(mass_x_l,mass_x_r,wx)
            zero_y=np.zeros((vertical_layers,1,nx),dtype=float)
            mass_y_l=np.concatenate((zero_y,vertical_q[0]),axis=1)
            mass_y_r=np.concatenate((vertical_q[0],zero_y),axis=1)
            openness_y=face_openness(mass_y_l,mass_y_r,wy)
            height_resolved_wall_steps+=1
        if horizontal_diffusivity_m2_s>0.:
            dt=min(dt,.20*h*h/(horizontal_diffusivity_m2_s
                                *(1.0+crosswind_diffusivity_factor)))
        ex=np.pad(q,((0,0),(0,0),(1,1)))
        ey=np.pad(q,((0,0),(1,1),(0,0)))
        fx_left,fx_right,shared_fx=_flux_pair(
            ex[:,:,:-1],ex[:,:,1:],0,thermo,wx,gravity_scale,
            wall_overtop_exponent,openness_override=openness_x,
            wall_blockage_fraction=np.where(
                wx>0.0,wall_blockage_fraction,0.0),return_shared=True)
        fy_left,fy_right,shared_fy=_flux_pair(
            ey[:,:-1,:],ey[:,1:,:],1,thermo,wy,gravity_scale,
            wall_overtop_exponent,openness_override=openness_y,
            wall_blockage_fraction=np.where(
                wy>0.0,wall_blockage_fraction,0.0),return_shared=True)
        raw_fx_left=raw_fx_right=raw_fy_left=raw_fy_right=None
        raw_shared_fx=raw_shared_fy=None
        if ((vertical_q is not None and (
                height_resolved_wall_flux
                or wall_impulse_vertical_lift_fraction>0.0))
                or elevated_species is not None
                or (wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                    and wall_impulse_elevated_transfer_fraction>0.0)
                or cavity_reservoir_capture_mode=='wall_impact'
                or wall_impulse_wake_capture_rate_s_inv>0.0):
            # Scalar layer transport applies donor-specific geometrical wall
            # eligibility itself. Obtain the corresponding unrestricted
            # column flux here; momentum keeps the hydrostatic wall reaction
            # from the primary solve above.
            raw_fx_left,raw_fx_right,raw_shared_fx=_flux_pair(
                ex[:,:,:-1],ex[:,:,1:],0,thermo,np.zeros_like(wx),
                gravity_scale,wall_overtop_exponent,return_shared=True)
            raw_fy_left,raw_fy_right,raw_shared_fy=_flux_pair(
                ey[:,:-1,:],ey[:,1:,:],1,thermo,np.zeros_like(wy),
                gravity_scale,wall_overtop_exponent,return_shared=True)
        overflow_advective_fx=overflow_advective_fy=None
        if wall_overflow_capture_fraction>0.0:
            def donor_advective_flux(shared,left_state,right_state,wall_mask):
                mass_flux=shared[0]
                donor=np.where((mass_flux>=0.0)[None,:,:],left_state,right_state)
                donor_mass=donor[0]
                species_fraction=np.divide(
                    donor[1],donor_mass,out=np.zeros_like(donor_mass),
                    where=donor_mass>1.e-30)
                velocity_x=np.divide(
                    donor[2],donor_mass,out=np.zeros_like(donor_mass),
                    where=donor_mass>1.e-30)
                velocity_y=np.divide(
                    donor[3],donor_mass,out=np.zeros_like(donor_mass),
                    where=donor_mass>1.e-30)
                advective=np.stack((
                    mass_flux,mass_flux*species_fraction,
                    mass_flux*velocity_x,mass_flux*velocity_y))
                return advective*wall_mask[None,:,:]
            internal_wall_x=wx>0.0
            internal_wall_x[:,0]=False;internal_wall_x[:,-1]=False
            internal_wall_y=wy>0.0
            internal_wall_y[0,:]=False;internal_wall_y[-1,:]=False
            storage_face_x=internal_wall_x|wall_overflow_side_mask_x
            storage_face_y=internal_wall_y|wall_overflow_side_mask_y
            overflow_advective_fx=donor_advective_flux(
                shared_fx,ex[:,:,:-1],ex[:,:,1:],storage_face_x)
            overflow_advective_fy=donor_advective_flux(
                shared_fy,ey[:,:-1,:],ey[:,1:,:],storage_face_y)
            # At storage-enabled wall faces use the donor's bounded species
            # fraction for the advective part.  Momentum pressure/Rusanov
            # residuals remain in the instantaneous shared flux and are not
            # misrepresented as material residing in the wake store.
            corrected_fx=shared_fx.copy()
            corrected_fy=shared_fy.copy()
            corrected_fx[1]=np.where(
                storage_face_x,overflow_advective_fx[1],shared_fx[1])
            corrected_fy[1]=np.where(
                storage_face_y,overflow_advective_fy[1],shared_fy[1])
            correction_fx=corrected_fx-shared_fx
            correction_fy=corrected_fy-shared_fy
            fx_left+=correction_fx;fx_right+=correction_fx
            fy_left+=correction_fy;fy_right+=correction_fy
            shared_fx=corrected_fx;shared_fy=corrected_fy
        # Record the physical wall reaction before any crossing material is
        # diverted into the internal overflow store.  Storage is an internal
        # transfer, not an additional obstacle impulse.
        wall_x_mask=wx[:,1:-1]>0.
        wall_y_mask=wy[1:-1,:]>0.
        wall_momentum_rate=np.array([
            -h*np.sum((fx_left[2,:,1:-1]-fx_right[2,:,1:-1])*wall_x_mask),
            -h*np.sum((fy_left[3,1:-1,:]-fy_right[3,1:-1,:])*wall_y_mask)])
        wall_impulse_approach_speed=np.zeros((ny,nx),dtype=float)
        wall_impulse_lift_velocity=np.zeros((ny,nx),dtype=float)
        if (wall_impulse_vertical_lift_fraction>0.0
                or elevated_species is not None
                or (wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                    and wall_impulse_elevated_transfer_fraction>0.0)
                or cavity_reservoir_capture_mode=='wall_impact'
                or wall_impulse_wake_capture_rate_s_inv>0.0):
            # Convert only the actually blocked normal approach speed into an
            # upward wall-impingement branch.  The unrestricted and finite-
            # wall mass fluxes establish the blocked fraction; the donor
            # momentum establishes its physical speed.  Numerical diffusion
            # with zero donor normal velocity therefore cannot create lift.
            def add_wall_lift(raw_shared,wall_shared,left_state,right_state,
                              wall_faces,axis):
                raw_mass=raw_shared[0]
                donor=np.where((raw_mass>=0.0)[None,:,:],left_state,right_state)
                donor_mass=donor[0]
                normal_speed=np.abs(np.divide(
                    donor[2+axis],donor_mass,out=np.zeros_like(donor_mass),
                    where=donor_mass>1.e-30))
                blocked_fraction=np.clip(np.divide(
                    np.maximum(np.abs(raw_mass)-np.abs(wall_shared[0]),0.0),
                    np.abs(raw_mass),out=np.zeros_like(raw_mass),
                    where=np.abs(raw_mass)>1.e-30),0.0,1.0)
                candidate=normal_speed*blocked_fraction*wall_faces
                if axis==0:
                    positive=raw_mass[:,1:nx]>=0.0
                    face=candidate[:,1:nx]
                    wall_impulse_approach_speed[:,:-1]=np.maximum(
                        wall_impulse_approach_speed[:,:-1],
                        np.where(positive,face,0.0))
                    wall_impulse_approach_speed[:,1:]=np.maximum(
                        wall_impulse_approach_speed[:,1:],
                        np.where(~positive,face,0.0))
                else:
                    positive=raw_mass[1:ny,:]>=0.0
                    face=candidate[1:ny,:]
                    wall_impulse_approach_speed[:-1,:]=np.maximum(
                        wall_impulse_approach_speed[:-1,:],
                        np.where(positive,face,0.0))
                    wall_impulse_approach_speed[1:,:]=np.maximum(
                        wall_impulse_approach_speed[1:,:],
                        np.where(~positive,face,0.0))
            internal_wall_x=wx>0.0
            internal_wall_x[:,0]=False;internal_wall_x[:,-1]=False
            internal_wall_y=wy>0.0
            internal_wall_y[0,:]=False;internal_wall_y[-1,:]=False
            add_wall_lift(raw_shared_fx,shared_fx,ex[:,:,:-1],ex[:,:,1:],
                          internal_wall_x,0)
            add_wall_lift(raw_shared_fy,shared_fy,ey[:,:-1,:],ey[:,1:,:],
                          internal_wall_y,1)
            wall_impulse_lift_velocity=(
                float(wall_impulse_vertical_lift_fraction)
                *wall_impulse_approach_speed)
            maximum_wall_impulse_vertical_velocity=max(
                maximum_wall_impulse_vertical_velocity,
                float(np.max(wall_impulse_lift_velocity)))
        if wall_overflow_capture_fraction>0.0:
            wall_mask_x=wx>0.0
            wall_mask_x[:,0]=False;wall_mask_x[:,-1]=False
            wall_mask_y=wy>0.0
            wall_mask_y[0,:]=False;wall_mask_y[-1,:]=False
            storage_mask_x=wall_mask_x|wall_overflow_side_mask_x
            storage_mask_y=wall_mask_y|wall_overflow_side_mask_y
            captured_fx=(float(wall_overflow_capture_fraction)
                         *overflow_advective_fx*storage_mask_x[None,:,:])
            captured_fy=(float(wall_overflow_capture_fraction)
                         *overflow_advective_fy*storage_mask_y[None,:,:])
            if wall_overflow_side_gate_length_m>0.0:
                outward_x=((wall_overflow_side_outward_sign_x>0)&(captured_fx[0]>0.0)
                           |(wall_overflow_side_outward_sign_x<0)&(captured_fx[0]<0.0))
                outward_y=((wall_overflow_side_outward_sign_y>0)&(captured_fy[0]>0.0)
                           |(wall_overflow_side_outward_sign_y<0)&(captured_fy[0]<0.0))
                captured_fx=np.where(
                    (wall_mask_x|outward_x)[None,:,:],captured_fx,0.0)
                captured_fy=np.where(
                    (wall_mask_y|outward_y)[None,:,:],captured_fy,0.0)
            for captured_flux,reservoir in (
                    (captured_fx,wall_overflow_x),
                    (captured_fy,wall_overflow_y)):
                positive=captured_flux[0]>0.0
                negative=captured_flux[0]<0.0
                packet=captured_flux*dt*h
                positive_packet=np.where(positive[None,:,:],packet,0.0)
                negative_packet=np.where(negative[None,:,:],-packet,0.0)
                for candidate in (positive_packet,negative_packet):
                    if (np.any(candidate[0]<-1.e-12)
                            or np.any(candidate[1]<-1.e-12)
                            or np.any(candidate[1]>candidate[0]+1.e-10)):
                        raise RuntimeError(
                            'wall overflow flux produced an inadmissible stored packet')
                reservoir[0]+=positive_packet
                reservoir[1]+=negative_packet
                wall_overflow_captured+=positive_packet.sum(axis=(1,2))
                wall_overflow_captured+=negative_packet.sum(axis=(1,2))
            if elevated_species is not None:
                def store_captured_marker(
                        captured_flux,axis,marker,age_moment,
                        marker_store,age_store):
                    if axis==0:
                        padded_marker=np.pad(marker,((0,0),(1,1)))
                        padded_age=np.pad(age_moment,((0,0),(1,1)))
                        left_marker=padded_marker[:,:-1]
                        right_marker=padded_marker[:,1:]
                        left_age=padded_age[:,:-1]
                        right_age=padded_age[:,1:]
                        left_species=ex[1,:,:-1]
                        right_species=ex[1,:,1:]
                    else:
                        padded_marker=np.pad(marker,((1,1),(0,0)))
                        padded_age=np.pad(age_moment,((1,1),(0,0)))
                        left_marker=padded_marker[:-1,:]
                        right_marker=padded_marker[1:,:]
                        left_age=padded_age[:-1,:]
                        right_age=padded_age[1:,:]
                        left_species=ey[1,:-1,:]
                        right_species=ey[1,1:,:]
                    positive=captured_flux[1]>=0.0
                    donor_marker=np.where(positive,left_marker,right_marker)
                    donor_age=np.where(positive,left_age,right_age)
                    donor_species=np.where(positive,left_species,right_species)
                    marker_fraction=np.divide(
                        donor_marker,donor_species,
                        out=np.zeros_like(donor_marker),
                        where=donor_species>1.e-30)
                    marker_fraction=np.clip(marker_fraction,0.0,1.0)
                    marker_flux=captured_flux[1]*marker_fraction
                    age_per_marker=np.divide(
                        donor_age,donor_marker,
                        out=np.zeros_like(donor_age),
                        where=donor_marker>1.e-30)
                    age_flux=marker_flux*age_per_marker
                    positive_packet=np.where(
                        positive,marker_flux*dt*h,0.0)
                    negative_packet=np.where(
                        ~positive,-marker_flux*dt*h,0.0)
                    positive_age=np.where(
                        positive,age_flux*dt*h,0.0)
                    negative_age=np.where(
                        ~positive,-age_flux*dt*h,0.0)
                    marker_store[0]+=positive_packet
                    marker_store[1]+=negative_packet
                    age_store[0]+=positive_age
                    age_store[1]+=negative_age
                store_captured_marker(
                    captured_fx,0,elevated_species,
                    elevated_species_age_moment,
                    wall_overflow_elevated_x,
                    wall_overflow_elevated_age_x)
                store_captured_marker(
                    captured_fy,1,elevated_species,
                    elevated_species_age_moment,
                    wall_overflow_elevated_y,
                    wall_overflow_elevated_age_y)
                store_captured_marker(
                    captured_fx,0,wall_flux_elevated_species,
                    wall_flux_elevated_age_moment,
                    wall_overflow_flux_marker_x,
                    wall_overflow_flux_age_x)
                store_captured_marker(
                    captured_fy,1,wall_flux_elevated_species,
                    wall_flux_elevated_age_moment,
                    wall_overflow_flux_marker_y,
                    wall_overflow_flux_age_y)
            wall_overflow_top_captured+=np.abs(
                captured_fx*wall_mask_x[None,:,:]).sum(axis=(1,2))*dt*h
            wall_overflow_top_captured+=np.abs(
                captured_fy*wall_mask_y[None,:,:]).sum(axis=(1,2))*dt*h
            wall_overflow_side_captured+=np.abs(
                captured_fx*wall_overflow_side_mask_x[None,:,:]).sum(axis=(1,2))*dt*h
            wall_overflow_side_captured+=np.abs(
                captured_fy*wall_overflow_side_mask_y[None,:,:]).sum(axis=(1,2))*dt*h
            # The donor must still lose the complete crossing packet; only
            # the receiver-side flux is withheld.  Subtracting from both
            # sides would leave the captured amount in the donor and create
            # a duplicate when it is added to storage.
            pos_x=captured_fx[0]>0.0;neg_x=captured_fx[0]<0.0
            pos_y=captured_fy[0]>0.0;neg_y=captured_fy[0]<0.0
            fx_right-=np.where(pos_x[None,:,:],captured_fx,0.0)
            fx_left-=np.where(neg_x[None,:,:],captured_fx,0.0)
            fy_right-=np.where(pos_y[None,:,:],captured_fy,0.0)
            fy_left-=np.where(neg_y[None,:,:],captured_fy,0.0)
        rates=np.concatenate((-fx_right[1,:,0],fx_left[1,:,-1],-fy_right[1,0,:],fy_left[1,-1,:]))*h
        boundary_momentum_rate=np.array([
            (-fx_right[2,:,0].sum()+fx_left[2,:,-1].sum()
             -fy_right[2,0,:].sum()+fy_left[2,-1,:].sum())*h,
            (-fx_right[3,:,0].sum()+fx_left[3,:,-1].sum()
             -fy_right[3,0,:].sum()+fy_left[3,-1,:].sum())*h])
        # Internal finite-height faces exchange momentum with the obstacle.
        # Include that wall impulse in the ledger rather than misclassifying
        # reflection as non-conservation of the fluid state.
        totalout=float(rates.sum()*dt)
        mixtureout=float((-fx_right[0,:,0].sum()+fx_left[0,:,-1].sum()-fy_right[0,0,:].sum()+fy_left[0,-1,:].sum())*h*dt)
        interface_step_injected=np.zeros(4,dtype=float)
        if interface_face_index is not None:
            # The upstream side remains an open outlet; the downstream side
            # receives only the declared time-resolved packet.  This prevents
            # an interior handoff from counting the same marked mass twice.
            k=interface_face_index
            raw_interface_rate=fx_left[:,:,k].copy()
            interface_negative_mass += float(
                -np.minimum(raw_interface_rate[0],0.0).sum()*h*dt
            )
            # This is a one-way interface contract. Reverse carrier flow is
            # not interpreted as an incoming packet; it is clipped and
            # recorded so the candidate can fail closed if the plane is
            # materially bidirectional.
            interface_out_rate=raw_interface_rate.copy()
            interface_out_rate[0]=np.maximum(raw_interface_rate[0],0.0)
            interface_out_rate[1]=np.clip(
                raw_interface_rate[1],0.0,interface_out_rate[0]
            )
            reverse=interface_out_rate[0]<=0.0
            interface_out_rate[2:,reverse]=0.0
            fx_left[:,:,k]=interface_out_rate
            interface_exported += interface_out_rate.sum(axis=1)*h*dt
            fx_right[:,:,k]=0.0
            if interface_source_callback is not None:
                provided=interface_source_callback(
                    float(interface_time_offset_s+times[-1]), float(dt), int(ny), float(h),
                )
                provided=np.asarray(provided,dtype=float)
                if provided.shape != (4,ny):
                    raise ValueError('interface callback must return shape (4, ny)')
                if np.any(~np.isfinite(provided)) or np.any(provided[0] < -1.0e-12):
                    raise ValueError('interface total packet must be finite and non-negative')
                if np.any(provided[1] < -1.0e-12) or np.any(provided[1] > provided[0] + 1.0e-10):
                    raise ValueError('interface species packet must lie in [0, total]')
                interface_step_injected=provided.copy()*dt
                interface_injected += interface_step_injected.sum(axis=1)
        if vertical_q is not None:
            # Carry the actual donor-side vertical distribution through every
            # solved horizontal face.  The previous implementation updated
            # only the aggregate column and deposited every newly populated
            # column into layer zero, erasing obstacle-overtopping and source
            # elevation information after each horizontal step.
            for component in (0,1):
                scalar_fx_left=(np.where(wx>0.0,raw_fx_left[component],fx_left[component])
                                if raw_fx_left is not None else fx_left[component])
                scalar_fx_right=(np.where(wx>0.0,raw_fx_right[component],fx_right[component])
                                 if raw_fx_right is not None else fx_right[component])
                scalar_fy_left=(np.where(wy>0.0,raw_fy_left[component],fy_left[component])
                                if raw_fy_left is not None else fy_left[component])
                scalar_fy_right=(np.where(wy>0.0,raw_fy_right[component],fy_right[component])
                                 if raw_fy_right is not None else fy_right[component])
                vfx_left=split_multilayer_face_flux(
                    scalar_fx_left,component,0,wx)
                vfx_right=split_multilayer_face_flux(
                    scalar_fx_right,component,0,wx)
                vfy_left=split_multilayer_face_flux(
                    scalar_fy_left,component,1,wy)
                vfy_right=split_multilayer_face_flux(
                    scalar_fy_right,component,1,wy)
                fx_left[component]=vfx_left.sum(axis=0)
                fx_right[component]=vfx_right.sum(axis=0)
                fy_left[component]=vfy_left.sum(axis=0)
                fy_right[component]=vfy_right.sum(axis=0)
                vertical_q[component]-=dt/h*(
                    vfx_left[:,:,1:]-vfx_right[:,:,:-1]
                    +vfy_left[:,1:,:]-vfy_right[:,:-1,:])
            if np.min(vertical_q[:2]) < -1.e-10:
                raise RuntimeError(
                    'multilayer horizontal transport positivity gate failed: '
                    f'min={float(np.min(vertical_q[:2])):.17g}, '
                    f'time_s={float(times[-1]):.17g}, dt_s={float(dt):.17g}')
        if elevated_species is not None:
            efx_left=split_elevated_species_flux(fx_left[1],0)
            efx_right=split_elevated_species_flux(fx_right[1],0)
            efy_left=split_elevated_species_flux(fy_left[1],1)
            efy_right=split_elevated_species_flux(fy_right[1],1)
            elevated_species-=dt/h*(
                efx_left[:,1:]-efx_right[:,:-1]
                +efy_left[1:,:]-efy_right[:-1,:])
            age_fx_left=split_elevated_species_flux(
                fx_left[1],0,elevated_species_age_moment)
            age_fx_right=split_elevated_species_flux(
                fx_right[1],0,elevated_species_age_moment)
            age_fy_left=split_elevated_species_flux(
                fy_left[1],1,elevated_species_age_moment)
            age_fy_right=split_elevated_species_flux(
                fy_right[1],1,elevated_species_age_moment)
            elevated_species_age_moment-=dt/h*(
                age_fx_left[:,1:]-age_fx_right[:,:-1]
                +age_fy_left[1:,:]-age_fy_right[:-1,:])
            flux_marker_fx_left=split_elevated_species_flux(
                fx_left[1],0,wall_flux_elevated_species)
            flux_marker_fx_right=split_elevated_species_flux(
                fx_right[1],0,wall_flux_elevated_species)
            flux_marker_fy_left=split_elevated_species_flux(
                fy_left[1],1,wall_flux_elevated_species)
            flux_marker_fy_right=split_elevated_species_flux(
                fy_right[1],1,wall_flux_elevated_species)
            wall_flux_elevated_species-=dt/h*(
                flux_marker_fx_left[:,1:]-flux_marker_fx_right[:,:-1]
                +flux_marker_fy_left[1:,:]-flux_marker_fy_right[:-1,:])
            flux_age_fx_left=split_elevated_species_flux(
                fx_left[1],0,wall_flux_elevated_age_moment)
            flux_age_fx_right=split_elevated_species_flux(
                fx_right[1],0,wall_flux_elevated_age_moment)
            flux_age_fy_left=split_elevated_species_flux(
                fy_left[1],1,wall_flux_elevated_age_moment)
            flux_age_fy_right=split_elevated_species_flux(
                fy_right[1],1,wall_flux_elevated_age_moment)
            wall_flux_elevated_age_moment-=dt/h*(
                flux_age_fx_left[:,1:]-flux_age_fx_right[:,:-1]
                +flux_age_fy_left[1:,:]-flux_age_fy_right[:-1,:])
            if wall_flux_elevated_transfer_fraction>0.0:
                # Mark only the species packet that actually crosses a finite-
                # height obstacle face. The donor loses the ordinary species
                # flux through q below; the marker is created on the receiver
                # side, so ground/cavity material that never crossed the wall
                # is not spuriously removed from the low-level branch.
                fraction=float(wall_flux_elevated_transfer_fraction)
                wall_x=wx[:,1:nx]>0.0
                species_x=shared_fx[1,:,1:nx]
                marker_x=efx_left[:,1:nx]
                new_x=(fraction*np.maximum(
                    np.abs(species_x)-np.abs(marker_x),0.0)*wall_x*dt/h)
                positive_x=species_x>=0.0
                elevated_species[:,1:]+=np.where(positive_x,new_x,0.0)
                elevated_species[:,:-1]+=np.where(~positive_x,new_x,0.0)
                wall_flux_elevated_species[:,1:]+=(
                    np.where(positive_x,new_x,0.0))
                wall_flux_elevated_species[:,:-1]+=(
                    np.where(~positive_x,new_x,0.0))
                wall_y=wy[1:ny,:]>0.0
                species_y=shared_fy[1,1:ny,:]
                marker_y=efy_left[1:ny,:]
                new_y=(fraction*np.maximum(
                    np.abs(species_y)-np.abs(marker_y),0.0)*wall_y*dt/h)
                positive_y=species_y>=0.0
                elevated_species[1:,:]+=np.where(positive_y,new_y,0.0)
                elevated_species[:-1,:]+=np.where(~positive_y,new_y,0.0)
                wall_flux_elevated_species[1:,:]+=(
                    np.where(positive_y,new_y,0.0))
                wall_flux_elevated_species[:-1,:]+=(
                    np.where(~positive_y,new_y,0.0))
                wall_flux_elevated_species_marked+=float(
                    (new_x.sum()+new_y.sum())*area)
        packets.append(rates)
        heights.append(np.concatenate((height[:,0],height[:,-1],height[0,:],height[-1,:])))
        current_speeds.append(np.concatenate((c[:,0],c[:,-1],c[0,:],c[-1,:])))
        q-=dt/h*(fx_left[:,:,1:]-fx_right[:,:,:-1]+fy_left[:,1:,:]-fy_right[:,:-1,:])
        if elevated_species is not None:
            if (np.min(elevated_species)<-1.e-10
                    or np.max(elevated_species-q[1])>1.e-10
                    or np.min(wall_flux_elevated_species)<-1.e-10
                    or np.max(wall_flux_elevated_species-elevated_species)>1.e-10):
                raise RuntimeError('elevated species transport violated subset bounds')
            elevated_species=np.clip(elevated_species,0.0,np.maximum(q[1],0.0))
            wall_flux_elevated_species=np.clip(
                wall_flux_elevated_species,0.0,elevated_species)
        if wall_overflow_capture_fraction>0.0:
            return_factor=1.0-math.exp(-wall_overflow_return_rate_s_inv*dt)
            released_x=wall_overflow_x*return_factor
            released_y=wall_overflow_y*return_factor
            wall_overflow_x-=released_x
            wall_overflow_y-=released_y
            # Positive x/y packets enter the cell on the positive side of an
            # internal face; negative packets return to the negative side.
            q[:,:,1:]+=released_x[0,:,:,1:nx]/area
            q[:,:,:-1]+=released_x[1,:,:,1:nx]/area
            q[:,1:,:]+=released_y[0,:,1:ny,:]/area
            q[:,:-1,:]+=released_y[1,:,1:ny,:]/area
            if elevated_species is not None:
                def release_stored_marker(marker_store,age_store,axis):
                    released_marker=marker_store*return_factor
                    released_age=age_store*return_factor
                    marker_store-=released_marker
                    age_store-=released_age
                    target_marker=np.zeros_like(elevated_species)
                    target_age=np.zeros_like(elevated_species)
                    if axis==0:
                        target_marker[:,1:]+=released_marker[0,:,1:nx]/area
                        target_marker[:,:-1]+=released_marker[1,:,1:nx]/area
                        target_age[:,1:]+=released_age[0,:,1:nx]/area
                        target_age[:,:-1]+=released_age[1,:,1:nx]/area
                    else:
                        target_marker[1:,:]+=released_marker[0,1:ny,:]/area
                        target_marker[:-1,:]+=released_marker[1,1:ny,:]/area
                        target_age[1:,:]+=released_age[0,1:ny,:]/area
                        target_age[:-1,:]+=released_age[1,1:ny,:]/area
                    return target_marker,target_age
                returned_elevated_x,returned_elevated_age_x=(
                    release_stored_marker(
                        wall_overflow_elevated_x,
                        wall_overflow_elevated_age_x,0))
                returned_elevated_y,returned_elevated_age_y=(
                    release_stored_marker(
                        wall_overflow_elevated_y,
                        wall_overflow_elevated_age_y,1))
                returned_flux_x,returned_flux_age_x=release_stored_marker(
                    wall_overflow_flux_marker_x,wall_overflow_flux_age_x,0)
                returned_flux_y,returned_flux_age_y=release_stored_marker(
                    wall_overflow_flux_marker_y,wall_overflow_flux_age_y,1)
                elevated_species+=returned_elevated_x+returned_elevated_y
                elevated_species_age_moment+=(
                    returned_elevated_age_x+returned_elevated_age_y)
                wall_flux_elevated_species+=returned_flux_x+returned_flux_y
                wall_flux_elevated_age_moment+=(
                    returned_flux_age_x+returned_flux_age_y)
                if (np.min(elevated_species)<-1.e-10
                        or np.max(elevated_species-q[1])>1.e-10
                        or np.min(wall_flux_elevated_species)<-1.e-10
                        or np.max(wall_flux_elevated_species-elevated_species)>1.e-10):
                    raise RuntimeError(
                        'released wall-overflow marker violated subset bounds')
            wall_overflow_returned+=released_x.sum(axis=(0,2,3))
            wall_overflow_returned+=released_y.sum(axis=(0,2,3))
        if vertical_q is not None:
            horizontal_sum=vertical_q[:2].sum(axis=1)
            maximum_multilayer_horizontal_closure_residual=max(
                maximum_multilayer_horizontal_closure_residual,
                float(np.max(np.abs(horizontal_sum-q[:2]))))
            if maximum_multilayer_horizontal_closure_residual>1.e-8:
                raise RuntimeError('multilayer horizontal column closure gate failed')
            old_vertical=vertical_q.sum(axis=1)
        step_reentry_impulse=np.zeros(2,dtype=float)
        step_source_impulse=np.zeros(2,dtype=float)
        if source_injection_callback is not None:
            provided=np.asarray(source_injection_callback(float(times[-1]),float(dt)),dtype=float)
            if provided.shape != q.shape or np.any(~np.isfinite(provided)):
                raise ValueError('source injection callback must return finite q-shaped increments')
            if np.any(provided[:2]<-1.0e-12) or np.any(provided[1]>provided[0]+1.0e-10):
                raise ValueError('source injection callback must return non-negative admissible increments')
            q+=provided
            source_injected+=provided.sum(axis=(1,2))*area
            step_source_impulse+=provided[2:4].sum(axis=(1,2))*area
        if exterior_reentry_fraction>0.0:
            # Capture only outward packets at the wind-aligned outer x face.
            # The packet is delayed and returned to the registered cavity;
            # this is a finite-domain recirculation hypothesis, not a sensor
            # correction or a new source.
            raw_out=(fx_left[:,:,-1] if wind[0]>=0.0 else -fx_right[:,:,0]).copy()
            packet_rate=np.zeros_like(raw_out)
            packet_rate[0]=np.maximum(raw_out[0],0.0)
            packet_rate[1]=np.clip(raw_out[1],0.0,packet_rate[0])
            packet_rate[2:]=np.where(packet_rate[0][None,:]>0.0,raw_out[2:],0.0)
            packet=packet_rate*h*dt*float(exterior_reentry_fraction)
            if float(packet[0].sum())>0.0:
                exterior_reentry_queue.append((
                    float(times[-1]+dt+exterior_reentry_delay_s),packet))
                exterior_reentry_captured+=packet.sum(axis=1)
        if exterior_reentry_queue:
            due=[]
            while exterior_reentry_queue and float(exterior_reentry_queue[0][0])<=float(times[-1]+dt)+1.0e-12:
                due.append(exterior_reentry_queue.pop(0)[1])
            if due:
                returned=np.sum(np.stack(due,axis=0),axis=(0,2))
                if float(returned[0])>0.0 and cavity_mask.any():
                    rr=np.sqrt((xx[cavity_mask]-cx)**2+(yy[cavity_mask]-cy)**2)
                    weights=np.exp(-0.5*(rr/max(float(exterior_reentry_sigma_m),1.0e-30))**2)
                    wsum=float(weights.sum())
                    if wsum>0.0:
                        q[:,cavity_mask]+=returned[:,None]*weights[None,:]/(wsum*area)
                        exterior_reentry_returned+=returned
                        step_reentry_impulse+=returned[2:4]
                        exterior_reentry_momentum_impulse+=returned[2:4]
        if interface_face_index is not None and interface_source_callback is not None:
            q[:,:,interface_face_index] += interface_step_injected / area
        wall_bypass_impulse=np.zeros(2)
        if wall_tangential_bypass_fraction>0.0:
            q,wall_bypass_impulse,bypassed_mass=_apply_tangential_wall_bypass(
                q,origin_xy_m=origin_xy_m,cell_size_m=h,
                wall_x_height_m=wx,wall_y_height_m=wy,
                centre_xy_m=(cx,cy) if cavity_center_xy_m is not None else origin_xy_m,
                wind_xy_m_s=wind,dt_s=dt,
                fraction=wall_tangential_bypass_fraction,thermo=thermo,
                gravity_scale=gravity_scale,wall_overtop_exponent=wall_overtop_exponent)
            maximum_wall_bypass_mass_kg=max(maximum_wall_bypass_mass_kg,bypassed_mass)
        if wall_tip_bypass_fraction>0.0:
            q,tip_bypass_impulse,bypassed_mass=_apply_tangential_wall_bypass(
                q,origin_xy_m=origin_xy_m,cell_size_m=h,
                wall_x_height_m=wx,wall_y_height_m=wy,
                centre_xy_m=(cx,cy) if cavity_center_xy_m is not None else origin_xy_m,
                wind_xy_m_s=wind,dt_s=dt,
                fraction=wall_tip_bypass_fraction,thermo=thermo,
                gravity_scale=gravity_scale,wall_overtop_exponent=wall_overtop_exponent,
                tip_only=True,tip_power=wall_tip_bypass_power,
                tip_radius_m=cavity_radius_m)
            wall_bypass_impulse += tip_bypass_impulse
            maximum_wall_bypass_mass_kg=max(maximum_wall_bypass_mass_kg,bypassed_mass)
        reservoir_impulse=np.zeros(2)
        if cavity_reservoir_capture_rate_s_inv>0 or cavity_reservoir_return_rate_s_inv>0:
            # Sub-grid cavity storage: transfer the local conserved state to
            # a well-mixed recirculation reservoir, then return it uniformly
            # over the cavity.  The reservoir carries all four conserved
            # components, so this is an internal rearrangement rather than a
            # source or a sensor-specific correction.
            return_sigma=float(cavity_reservoir_return_sigma_m)
            if cavity_reservoir_pending:
                ready=[]
                waiting=[]
                cutoff=float(times[-1])+float(dt)+1.0e-12
                for due,packet in cavity_reservoir_pending:
                    (ready if due <= cutoff else waiting).append((due,packet))
                if ready:
                    cavity_reservoir += np.sum(
                        np.stack([packet for _,packet in ready],axis=0),axis=0)
                cavity_reservoir_pending=waiting
            if (cavity_reservoir_capture_rate_s_inv>0
                    and times[-1] >= cavity_reservoir_start_time_s):
                capture_factor=1.-math.exp(-cavity_reservoir_capture_rate_s_inv*dt)
                if cavity_reservoir_capture_mode=='wall_impact':
                    # Split only inventory that encounters a physically
                    # blocked finite-wall approach.  The normalized blocked
                    # approach speed is a geometry/flow weight, not a sensor
                    # correction.  All four conserved components follow the
                    # same transfer into the independently transported local
                    # wake state below.
                    reference_speed=max(
                        float(np.linalg.norm(reference_wind)),1.e-30)
                    capture_weight=np.clip(
                        wall_impulse_approach_speed[cavity_mask]
                        /reference_speed,0.0,1.0)
                elif cavity_reservoir_capture_radius_m>0.0:
                    capture_radius=np.sqrt((xx[cavity_mask]-cx)**2
                                           +(yy[cavity_mask]-cy)**2)
                    capture_weight=np.exp(-0.5*((capture_radius-
                        cavity_reservoir_capture_radius_m)/return_sigma)**2)
                else:
                    capture_weight=np.ones(int(np.count_nonzero(cavity_mask)))
                marker_capture_weight=None
                if cavity_reservoir_vertical_capture_mode=='below_wall_core':
                    # The fast outer carrier remains a single conservative
                    # column.  Only the geometrical portion below an adjacent
                    # finite wall, and only the non-elevated species share,
                    # enters the local slow recirculation state.  This is a
                    # height-resolved obstacle exchange without replacing the
                    # wall-free carrier by a global multilayer solve.
                    below_wall=np.clip(np.divide(
                        cell_wall_height[cavity_mask],
                        np.maximum(height[cavity_mask],1.e-30)),0.0,1.0)
                    core_share=np.ones_like(below_wall)
                    if elevated_species is not None:
                        core_share=np.clip(1.0-np.divide(
                            elevated_species[cavity_mask],q[1,cavity_mask],
                            out=np.zeros_like(below_wall),
                            where=q[1,cavity_mask]>1.e-30),0.0,1.0)
                    capture_weight*=below_wall*core_share
                    # The elevated marker is explicitly the prompt/over-wall
                    # complement and must remain outside the lower store.
                    marker_capture_weight=np.zeros_like(capture_weight)
                captured=q[:,cavity_mask]*(capture_factor*capture_weight)[None,:]
                if cavity_reservoir_content_mode=='species_only':
                    species_captured=np.zeros_like(captured)
                    species_captured[1]=captured[1]
                    captured=species_captured
                if elevated_species is not None:
                    marker_weight=(capture_weight if marker_capture_weight is None
                                   else marker_capture_weight)
                    captured_elevated=(elevated_species[cavity_mask]
                                       *capture_factor*marker_weight)
                    captured_elevated_age=(
                        elevated_species_age_moment[cavity_mask]
                        *capture_factor*marker_weight)
                    elevated_species[cavity_mask]-=captured_elevated
                    elevated_species_age_moment[cavity_mask]-=(
                        captured_elevated_age)
                    captured_wall_flux=(
                        wall_flux_elevated_species[cavity_mask]
                        *capture_factor*marker_weight)
                    captured_wall_flux_age=(
                        wall_flux_elevated_age_moment[cavity_mask]
                        *capture_factor*marker_weight)
                    wall_flux_elevated_species[cavity_mask]-=(
                        captured_wall_flux)
                    wall_flux_elevated_age_moment[cavity_mask]-=(
                        captured_wall_flux_age)
                    if cavity_reservoir_spatial_mode=='local':
                        local_cavity_reservoir_elevated_species[cavity_mask]+=(
                            captured_elevated)
                        local_cavity_reservoir_elevated_age_moment[cavity_mask]+=(
                            captured_elevated_age)
                        local_cavity_reservoir_wall_flux_species[cavity_mask]+=(
                            captured_wall_flux)
                        local_cavity_reservoir_wall_flux_age_moment[cavity_mask]+=(
                            captured_wall_flux_age)
                    else:
                        cavity_reservoir_elevated_species+=float(
                            captured_elevated.sum()*area)
                        cavity_reservoir_elevated_age_moment+=float(
                            captured_elevated_age.sum()*area)
                        cavity_reservoir_wall_flux_species+=float(
                            captured_wall_flux.sum()*area)
                        cavity_reservoir_wall_flux_age_moment+=float(
                            captured_wall_flux_age.sum()*area)
                q[:,cavity_mask]-=captured
                captured_total=captured.sum(axis=1)*area
                if cavity_reservoir_spatial_mode=='local':
                    local_cavity_reservoir[:,cavity_mask]+=captured
                elif cavity_reservoir_return_delay_s > 0.0:
                    cavity_reservoir_pending.append(
                        (float(times[-1])+float(dt)+float(cavity_reservoir_return_delay_s),
                         captured_total.copy()))
                else:
                    cavity_reservoir+=captured_total
                if (cavity_reservoir_spatial_mode=='well_mixed'
                        and cavity_sector_index is not None):
                    for sector in range(4):
                        mask=cavity_sector_index==sector
                        if np.any(mask):
                            cavity_reservoir_sector[:,sector]+=captured[:,mask].sum(axis=1)*area
            if (cavity_reservoir_spatial_mode=='local'
                    and np.any(local_cavity_reservoir>0)):
                if (cavity_reservoir_transport_fraction>0.0
                        or cavity_reservoir_crosswind_convergence_rate_s_inv>0.0):
                    if cavity_reservoir_transport_height_m is None:
                        transport_velocity=(
                            float(cavity_reservoir_transport_fraction)*wind)
                    else:
                        wind_speed=float(np.linalg.norm(wind))
                        if wind_speed<=1.e-30:
                            transport_velocity=np.zeros(2,dtype=float)
                        else:
                            transport_speed=(
                                float(friction_velocity_m_s)/0.41
                                *math.log1p(float(
                                    cavity_reservoir_transport_height_m)
                                    /float(roughness_length_m)))
                            transport_velocity=(
                                float(cavity_reservoir_transport_fraction)
                                *transport_speed*wind/wind_speed)
                    cfl_x=abs(float(transport_velocity[0]))*dt/h
                    cfl_y=abs(float(transport_velocity[1]))*dt/h
                    if cfl_x+cfl_y>1.0+1.e-12:
                        raise RuntimeError('local reservoir transport CFL exceeded')
                    def advect_masked(field):
                        updated=field.copy()
                        if transport_velocity[0]>=0.0:
                            transfer=cfl_x*field[...,:,:-1]
                            gate=(reservoir_transport_mask[:,:-1]
                                  &reservoir_transport_mask[:,1:])
                            transfer*=gate
                            updated[...,:,:-1]-=transfer
                            updated[..., :,1:]+=transfer
                        else:
                            transfer=cfl_x*field[..., :,1:]
                            gate=(reservoir_transport_mask[:,:-1]
                                  &reservoir_transport_mask[:,1:])
                            transfer*=gate
                            updated[..., :,1:]-=transfer
                            updated[...,:,:-1]+=transfer
                        field=updated;updated=field.copy()
                        if transport_velocity[1]>=0.0:
                            transfer=cfl_y*field[...,:-1,:]
                            gate=(reservoir_transport_mask[:-1,:]
                                  &reservoir_transport_mask[1:,:])
                            transfer*=gate
                            updated[...,:-1,:]-=transfer
                            updated[...,1:,:]+=transfer
                        else:
                            transfer=cfl_y*field[...,1:,:]
                            gate=(reservoir_transport_mask[:-1,:]
                                  &reservoir_transport_mask[1:,:])
                            transfer*=gate
                            updated[...,1:,:]-=transfer
                            updated[...,:-1,:]+=transfer
                        return updated
                    local_cavity_reservoir=advect_masked(
                        local_cavity_reservoir)
                    if elevated_species is not None:
                        local_cavity_reservoir_elevated_species=advect_masked(
                            local_cavity_reservoir_elevated_species)
                        local_cavity_reservoir_elevated_age_moment=advect_masked(
                            local_cavity_reservoir_elevated_age_moment)
                        local_cavity_reservoir_wall_flux_species=advect_masked(
                            local_cavity_reservoir_wall_flux_species)
                        local_cavity_reservoir_wall_flux_age_moment=advect_masked(
                            local_cavity_reservoir_wall_flux_age_moment)
                    if cavity_reservoir_crosswind_convergence_rate_s_inv>0.0:
                        target_y=(float(cy)
                                  +float(cavity_reservoir_crosswind_target_offset_m))
                        y_faces=float(origin_xy_m[1])+np.arange(1,ny)*h
                        face_velocity=(-float(
                            cavity_reservoir_crosswind_convergence_rate_s_inv)
                            *(y_faces-target_y))
                        face_cfl=np.abs(face_velocity)*dt/h
                        if np.max(face_cfl)>1.0+1.e-12:
                            raise RuntimeError(
                                'local reservoir convergence CFL exceeded')
                        gate=(reservoir_transport_mask[:-1,:]
                              &reservoir_transport_mask[1:,:])
                        def converge_crosswind(field):
                            lead=(1,)*(field.ndim-2)
                            positive=(face_velocity>=0.0).reshape(
                                lead+(ny-1,1))
                            cfl=face_cfl.reshape(lead+(ny-1,1))
                            donor=np.where(
                                positive,field[...,:-1,:],field[...,1:,:])
                            transfer=cfl*donor*gate
                            upward=np.where(positive,transfer,0.0)
                            downward=np.where(~positive,transfer,0.0)
                            updated=field.copy()
                            updated[...,:-1,:]-=upward
                            updated[...,1:,:]+=upward
                            updated[...,1:,:]-=downward
                            updated[...,:-1,:]+=downward
                            return updated
                        local_cavity_reservoir=converge_crosswind(
                            local_cavity_reservoir)
                        if elevated_species is not None:
                            local_cavity_reservoir_elevated_species=(
                                converge_crosswind(
                                    local_cavity_reservoir_elevated_species))
                            local_cavity_reservoir_elevated_age_moment=(
                                converge_crosswind(
                                    local_cavity_reservoir_elevated_age_moment))
                            local_cavity_reservoir_wall_flux_species=(
                                converge_crosswind(
                                    local_cavity_reservoir_wall_flux_species))
                            local_cavity_reservoir_wall_flux_age_moment=(
                                converge_crosswind(
                                    local_cavity_reservoir_wall_flux_age_moment))
                if cavity_reservoir_return_rate_s_inv>0.0:
                    return_factor=1.-math.exp(
                        -cavity_reservoir_return_rate_s_inv*dt)
                    returned=local_cavity_reservoir*return_factor
                    local_cavity_reservoir-=returned
                    if cavity_reservoir_return_target=='main':
                        q+=returned
                    else:
                        # A separated lower cavity first resides in the local
                        # wake and then reattaches into its own conservative
                        # carrier.  This sequential path avoids forcing the
                        # returned packet immediately back into the prompt
                        # hydrostatic column while preserving all four
                        # conserved quantities exactly.
                        reattached_return=(
                            returned*cavity_reservoir_reattached_fraction)
                        main_return=returned-reattached_return
                        q+=main_return
                        wall_impulse_reattached_reservoir+=reattached_return
                        wall_impulse_reattached_captured+=(
                            reattached_return.sum(axis=(1,2))*area)
                    if elevated_species is not None:
                        returned_elevated=(
                            local_cavity_reservoir_elevated_species*return_factor)
                        returned_elevated_age=(
                            local_cavity_reservoir_elevated_age_moment*return_factor)
                        local_cavity_reservoir_elevated_species-=returned_elevated
                        local_cavity_reservoir_elevated_age_moment-=(
                            returned_elevated_age)
                        elevated_species+=returned_elevated
                        elevated_species_age_moment+=returned_elevated_age
                        returned_wall_flux=(
                            local_cavity_reservoir_wall_flux_species*return_factor)
                        returned_wall_flux_age=(
                            local_cavity_reservoir_wall_flux_age_moment*return_factor)
                        local_cavity_reservoir_wall_flux_species-=(
                            returned_wall_flux)
                        local_cavity_reservoir_wall_flux_age_moment-=(
                            returned_wall_flux_age)
                        wall_flux_elevated_species+=returned_wall_flux
                        wall_flux_elevated_age_moment+=returned_wall_flux_age
            if (cavity_reservoir_spatial_mode=='well_mixed'
                    and cavity_reservoir_return_rate_s_inv>0
                    and np.any(cavity_reservoir>0)):
                return_factor=1.-math.exp(-cavity_reservoir_return_rate_s_inv*dt)
                return_radius=np.sqrt(
                    (cavity_streamwise_coordinate[cavity_mask]
                     -float(cavity_reservoir_return_offset_m))**2
                    +(cavity_crosswind_coordinate[cavity_mask]
                      /float(cavity_reservoir_return_crosswind_scale))**2)
                if cavity_reservoir_return_radius_m>0.0:
                    return_weight=np.exp(-0.5*((return_radius-
                        cavity_reservoir_return_radius_m)/return_sigma)**2)
                else:
                    return_weight=np.exp(-0.5*(return_radius/return_sigma)**2)
                if cavity_return_streamwise_bias != 0.0:
                    streamwise=np.clip(
                        cavity_streamwise_coordinate[cavity_mask]
                        /max(float(cavity_radius_m),1e-30),-1.0,1.0)
                    return_weight*=np.maximum(
                        1.0+float(cavity_return_streamwise_bias)*streamwise,0.0)
                if (cavity_sector_index is not None
                        and cavity_reservoir_face_return_fraction > 0.0):
                    # Preserve a bounded portion of the captured face/sector
                    # history. The remainder follows the existing isotropic
                    # ring closure; every sector is normalized independently
                    # so total/species/momentum remain exactly conserved.
                    returned_sector=cavity_reservoir_sector*return_factor
                    cavity_reservoir_sector-=returned_sector
                    returned=returned_sector.sum(axis=1)
                    cavity_reservoir-=returned
                    weight_sum=max(float(return_weight.sum()),1e-30)
                    q[:,cavity_mask]+=((1.0-cavity_reservoir_face_return_fraction)
                        *returned[:,None]*return_weight[None,:]/(weight_sum*area))
                    for sector in range(4):
                        mask=cavity_sector_index==sector
                        sector_weight=return_weight*mask
                        sector_sum=max(float(sector_weight.sum()),1e-30)
                        q[:,cavity_mask]+= (cavity_reservoir_face_return_fraction
                            *returned_sector[:,sector,None]
                            *sector_weight[None,:]/(sector_sum*area))
                else:
                    returned=cavity_reservoir*return_factor
                    cavity_reservoir-=returned
                    weight_sum=max(float(return_weight.sum()),1e-30)
                    q[:,cavity_mask]+=returned[:,None]*return_weight[None,:]/(weight_sum*area)
                    if elevated_species is not None:
                        returned_elevated=(cavity_reservoir_elevated_species
                                           *return_factor)
                        returned_elevated_age=(
                            cavity_reservoir_elevated_age_moment*return_factor)
                        cavity_reservoir_elevated_species-=returned_elevated
                        cavity_reservoir_elevated_age_moment-=(
                            returned_elevated_age)
                        elevated_species[cavity_mask]+=(
                            returned_elevated*return_weight/(weight_sum*area))
                        elevated_species_age_moment[cavity_mask]+=(
                            returned_elevated_age*return_weight/(weight_sum*area))
                        returned_wall_flux=(
                            cavity_reservoir_wall_flux_species*return_factor)
                        returned_wall_flux_age=(
                            cavity_reservoir_wall_flux_age_moment*return_factor)
                        cavity_reservoir_wall_flux_species-=(
                            returned_wall_flux)
                        cavity_reservoir_wall_flux_age_moment-=(
                            returned_wall_flux_age)
                        wall_flux_elevated_species[cavity_mask]+=(
                            returned_wall_flux*return_weight/(weight_sum*area))
                        wall_flux_elevated_age_moment[cavity_mask]+=(
                            returned_wall_flux_age*return_weight/(weight_sum*area))
        if wall_impulse_wake_capture_rate_s_inv>0.0:
            # A second conserved carrier for the separated, slow wall wake.
            # The ordinary cavity reservoir remains untouched.  Only cells
            # with a nonzero blocked-face approach speed can enter this
            # branch, so prompt over-wall material and slow wake material can
            # coexist instead of replacing one another.
            reference_speed=max(
                float(np.linalg.norm(reference_wind)),1.e-30)
            impact_weight=np.clip(
                wall_impulse_approach_speed/reference_speed,0.0,1.0)
            capture_factor=(1.0-np.exp(
                -float(wall_impulse_wake_capture_rate_s_inv)
                *impact_weight
                *(np.divide(q[1],q[0],out=np.zeros_like(q[1]),
                            where=q[0]>1.e-30)
                  if wall_impulse_wake_capture_mode=='species_weighted'
                  else 1.0)
                *dt))
            captured=q*capture_factor[None,:,:]
            q-=captured
            wall_impulse_wake_reservoir+=captured
            wall_impulse_wake_captured+=captured.sum(axis=(1,2))*area
            if elevated_species is not None:
                # The separated branch is the lower wake path, not the
                # prompt elevated path.  Remove the matching tagged subset;
                # physical species remains conserved in the wake reservoir.
                elevated_species*=(1.0-capture_factor)
                elevated_species_age_moment*=(1.0-capture_factor)
                wall_flux_elevated_species*=(1.0-capture_factor)
                wall_flux_elevated_age_moment*=(1.0-capture_factor)
            if (wall_impulse_wake_transport_fraction>0.0
                    or wall_impulse_wake_crosswind_velocity_m_s!=0.0
                    or wall_impulse_wake_crosswind_outward_velocity_m_s>0.0):
                wind_speed=max(float(np.linalg.norm(wind)),1.e-30)
                cross_x=-float(wind[1])/wind_speed
                cross_y=float(wind[0])/wind_speed
                cross_coordinate=(xx-float(cx))*cross_x+(yy-float(cy))*cross_y
                wake_crosswind=(
                    float(wall_impulse_wake_crosswind_velocity_m_s)
                    +float(wall_impulse_wake_crosswind_outward_velocity_m_s)
                    *np.tanh(cross_coordinate/float(
                        wall_impulse_wake_crosswind_outward_scale_m)))
                carrier_u=(np.full(
                    (1,ny,nx),float(wall_impulse_wake_transport_fraction)
                    *float(wind[0]))+wake_crosswind[None,:,:]*cross_x)
                carrier_v=(np.full(
                    (1,ny,nx),float(wall_impulse_wake_transport_fraction)
                    *float(wind[1]))+wake_crosswind[None,:,:]*cross_y)
                transported,moved,_=_transport_multilayer_relative(
                    wall_impulse_wake_reservoir[:,None,:,:],carrier_u,carrier_v,
                    dt_s=dt,cell_size_m=h,
                    wall_x_height_m=np.zeros_like(wx),
                    wall_y_height_m=np.zeros_like(wy),
                    z_edges_m=np.array([0.0,1.0]))
                wall_impulse_wake_reservoir=transported[:,0]
                # Reuse the conservative four-component transport kernel for
                # the scalar moment by duplicating it into the admissible
                # mass/species slots.  Both copies receive identical fluxes.
                age_transport_state=np.zeros_like(transported)
                age_transport_state[0,0]=wall_impulse_wake_species_age_moment
                age_transport_state[1,0]=wall_impulse_wake_species_age_moment
                age_transport_state,_,_=_transport_multilayer_relative(
                    age_transport_state,carrier_u,carrier_v,
                    dt_s=dt,cell_size_m=h,
                    wall_x_height_m=np.zeros_like(wx),
                    wall_y_height_m=np.zeros_like(wy),
                    z_edges_m=np.array([0.0,1.0]))
                wall_impulse_wake_species_age_moment=age_transport_state[0,0]
                wall_impulse_wake_transport_moved+=moved
            return_factor=(1.0-math.exp(
                -float(wall_impulse_wake_return_rate_s_inv)*dt))
            returned=wall_impulse_wake_reservoir*return_factor
            returned_age=wall_impulse_wake_species_age_moment*return_factor
            wall_impulse_wake_reservoir-=returned
            wall_impulse_wake_species_age_moment-=returned_age
            q+=returned
            wall_impulse_wake_returned+=returned.sum(axis=(1,2))*area
            # New captures enter with age zero; inventory remaining after the
            # return step accrues one time-step of residence.
            wall_impulse_wake_species_age_moment+=(
                wall_impulse_wake_reservoir[1]*dt)
        shelter_impulse=np.zeros(2)
        vortex_impulse=np.zeros(2)
        if cavity_vortex_rate_s_inv>0:
            # A resolved obstacle cavity is represented as a weakly forced
            # recirculation cell.  The target is counter-clockwise for a
            # positive-x inflow: return flow at the top of the cavity points
            # upstream, while the lower branch points downstream.  This is a
            # momentum exchange with the obstacle, not a sensor calibration.
            before_vortex=q[2:4,cavity_mask].copy()
            cxx=(xx[cavity_mask]-cx)/max(float(cavity_radius_m),1e-30)
            cyy=(yy[cavity_mask]-cy)/max(float(cavity_radius_m),1e-30)
            rr=np.sqrt(cxx*cxx+cyy*cyy)
            radial=np.clip(rr,0.,1.)
            tangent=np.column_stack((-cyy,cxx))
            tangent_norm=np.sqrt(np.sum(tangent*tangent,axis=1))
            tangent=np.divide(tangent,np.maximum(tangent_norm[:,None],1e-30))
            # Vanish at the cavity centre and wall; peak circulation is in
            # the shear layer where a separated wake is expected.
            envelope=4.*radial*(1.-radial)
            target=np.zeros((2,int(np.count_nonzero(cavity_mask))))
            target[:]=tangent.T*(float(np.linalg.norm(wind))*cavity_vortex_strength*envelope)[None,:]
            local_mass=q[0,cavity_mask]
            target_momentum=target*local_mass[None,:]
            factor=1.-np.exp(-cavity_vortex_rate_s_inv*dt)
            q[2:4,cavity_mask]=(1.-factor)*q[2:4,cavity_mask]+factor*target_momentum
            vortex_impulse=(q[2:4,cavity_mask]-before_vortex).sum(axis=1)*area
        if cavity_shelter_strength>0:
            before_shelter=q[2:4,cavity_mask].copy()
            factor=math.exp(-cavity_shelter_strength*dt)
            q[2:4,cavity_mask]*=factor
            shelter_impulse=(q[2:4,cavity_mask]-before_shelter).sum(axis=1)*area
        backflow_impulse=np.zeros(2)
        if cavity_backflow_rate_s_inv > 0.0:
            # Relax the conserved cavity momentum toward the same bounded
            # separated-bubble carrier used by entrainment.  This is an
            # obstacle reaction (recorded in the momentum ledger), not a
            # receptor-specific concentration correction.
            radius=max(float(cavity_radius_m),1e-30)
            if cavity_backflow_frame=='legacy_x':
                legacy_dx=xx-float(cavity_center_xy_m[0])
                legacy_dy=yy-float(cavity_center_xy_m[1])
                downstream=np.clip(legacy_dx/radius,0.0,1.0)
                radial=np.clip(
                    np.sqrt(legacy_dx**2+legacy_dy**2)/radius,0.0,1.0)
                bubble=cavity_backflow_strength*downstream*(0.5+0.5*radial)
                target_u=wind[0]*(1.0-2.0*bubble)
                target_v=(-wind[0]*0.5*bubble
                          *np.clip(legacy_dy/radius,-1.0,1.0))
            else:
                downstream=np.clip(
                    cavity_streamwise_coordinate/radius,0.0,1.0)
                radial=np.clip(np.sqrt(
                    cavity_streamwise_coordinate**2
                    +cavity_crosswind_coordinate**2)/radius,0.0,1.0)
                bubble=cavity_backflow_strength*downstream*(0.5+0.5*radial)
                along_speed=float(np.dot(wind,cavity_streamwise_unit))
                crosswind_speed=float(np.dot(wind,cavity_crosswind_unit))
                target_along=along_speed*(1.0-2.0*bubble)
                target_crosswind=(
                    crosswind_speed
                    -along_speed*0.5*bubble*np.clip(
                        cavity_crosswind_coordinate/radius,-1.0,1.0))
                target_u=(target_along*cavity_streamwise_unit[0]
                          +target_crosswind*cavity_crosswind_unit[0])
                target_v=(target_along*cavity_streamwise_unit[1]
                          +target_crosswind*cavity_crosswind_unit[1])
            target=np.stack((target_u[cavity_mask],target_v[cavity_mask]))
            before_backflow=q[2:4,cavity_mask].copy()
            factor=1.0-math.exp(-cavity_backflow_rate_s_inv*dt)
            q[2:4,cavity_mask]=(1.0-factor)*before_backflow + factor*(q[0,cavity_mask][None,:]*target)
            backflow_impulse=(q[2:4,cavity_mask]-before_backflow).sum(axis=1)*area
        pending_total=(np.sum(np.stack([packet for _,packet in cavity_reservoir_pending],axis=0),axis=0)
                       if cavity_reservoir_pending else np.zeros(4))
        overflow_total=wall_overflow_total()
        momentum_after_transport=(q[2:4].sum(axis=(1,2))*area
                                  +cavity_reservoir_total()[2:4]+pending_total[2:4]
                                  +overflow_total[2:4])
        transport_expected=(momentum_before_transport
                            -dt*boundary_momentum_rate
                            +dt*wall_momentum_rate+shelter_impulse+vortex_impulse
                            +backflow_impulse+wall_bypass_impulse
                            +step_reentry_impulse+step_source_impulse)
        if interface_face_index is not None:
            transport_expected -= interface_out_rate[2:4].sum(axis=1)*h*dt
            transport_expected += interface_step_injected[2:4].sum(axis=1)
        _transport_residual=float(np.max(np.abs(momentum_after_transport-transport_expected)))
        if _transport_residual>maximum_momentum_budget_residual:
            maximum_momentum_budget_residual=_transport_residual
            largest_transport_budget_detail={
                'time_s':float(times[-1]+dt),
                'actual':momentum_after_transport.tolist(),
                'expected':transport_expected.tolist(),
                'source_impulse':step_source_impulse.tolist(),
                'reentry_impulse':step_reentry_impulse.tolist(),
            }
        if horizontal_diffusivity_m2_s>0.:
            # Diagnostic unresolved turbulent mixing.  Edge padding imposes
            # zero normal diffusive flux at the patch boundary; the same
            # conservative operator is applied to mixture, species and
            # momentum so 0<=species<=mixture remains an invariant under the
            # explicit stability bound below.
            diffusion_number=(horizontal_diffusivity_m2_s
                              *(1.0+crosswind_diffusivity_factor)*dt/(h*h))
            if diffusion_number>0.25+1.0e-12:
                raise RuntimeError('horizontal diffusion timestep stability gate failed')
            direction=wind/np.linalg.norm(wind)
            dpar=horizontal_diffusivity_m2_s
            dcross=dpar*crosswind_diffusivity_factor
            dxx=dpar*direction[0]**2+dcross*direction[1]**2
            dyy=dpar*direction[1]**2+dcross*direction[0]**2
            # The mixed derivative is omitted deliberately: on this coarse
            # FV patch the nine-point cross stencil is not positivity
            # preserving.  The diagonal tensor is the conservative,
            # bounded approximation at the measured six-degree wind angle.
            if horizontal_diffusion_respect_walls:
                diffusion_height,*_=_primitive(q,thermo)
                q,blocked_species=_diffuse_horizontal_through_open_faces(
                    q,diffusion_height,dxx=dxx,dyy=dyy,dt_s=dt,
                    cell_size_m=h,wall_x_height_m=wx,wall_y_height_m=wy,
                    wall_blockage_fraction=wall_blockage_fraction)
                horizontal_diffusion_blocked_species_throughput_kg+=(
                    blocked_species)
            else:
                # Historical unrestricted stencil, retained explicitly for
                # exact replay of old research artifacts. It does not impose
                # an impermeable wall on turbulent diffusive flux.
                padded=np.pad(q,((0,0),(1,1),(1,1)),mode='edge')
                qxx=(padded[:,1:-1,:-2]-2.*padded[:,1:-1,1:-1]
                     +padded[:,1:-1,2:])/(h*h)
                qyy=(padded[:,:-2,1:-1]-2.*padded[:,1:-1,1:-1]
                     +padded[:,2:,1:-1])/(h*h)
                q+=dt*(dxx*qxx+dyy*qyy)
        if np.any(~np.isfinite(q)) or np.min(q[:2])<-1e-10 or np.max(q[1]-q[0])>1e-10:
            raise RuntimeError('layer positivity/finite gate failed; no clipping applied')
        if vertical_q is not None:
            for component in (0,1):
                ratio=np.divide(q[component],old_vertical[component],
                                out=np.zeros_like(q[component]),
                                where=old_vertical[component]>1e-30)
                vertical_q[component]*=ratio[None,:,:]
                # A Rusanov step can populate a previously dry column. Put
                # that newly transported conserved amount in the ground layer
                # before vertical mixing; never silently lose it.
                vertical_q[component,0]+=q[component]-vertical_q[component].sum(axis=0)
            column_mass=vertical_q[0].sum(axis=0,keepdims=True)
            weights=np.divide(vertical_q[0],column_mass,
                              out=np.zeros_like(vertical_q[0]),where=column_mass>1e-30)
            vertical_q[2]=q[2][None,:,:]*weights
            vertical_q[3]=q[3][None,:,:]*weights
        height,u,v,_,_,_=_primitive(q,thermo)
        # ``broadcast_to`` returns a read-only view; wake diagnostics apply a
        # cell-local deficit, so materialize a writable carrier field first.
        local_wind=wind_field_at_time(times[-1]+0.5*dt)
        if friction_velocity_m_s is not None:
            # Analytic depth average of neutral log wind. Numerically stable
            # at vanishing h using log1p; no minimum sensor height substituted.
            ratio=height/roughness_length_m
            mean=np.divide((1+ratio)*np.log1p(ratio)-ratio,ratio,
                           out=np.zeros_like(ratio),where=ratio>1e-20)*friction_velocity_m_s/.4
            if spatial_wind_history:
                local_speed=np.linalg.norm(local_wind,axis=0)
                direction=np.divide(
                    local_wind,local_speed[None,:,:],out=np.zeros_like(local_wind),
                    where=local_speed[None,:,:]>1.0e-30)
                friction_scale=(
                    local_speed/reference_wind_speed
                    if wind_history_scale_friction_velocity else 1.0
                )
                local_wind=direction*(mean*friction_scale)[None,:,:]
            else:
                direction=wind/np.linalg.norm(wind)
                friction_scale=(
                    float(np.linalg.norm(wind))/reference_wind_speed
                    if history_times is not None and wind_history_scale_friction_velocity
                    else 1.0
                )
                local_wind=direction[:,None,None]*(mean*friction_scale)
        if wake_deficit_strength>0:
            direction=wind/np.linalg.norm(wind)
            xx,yy=np.meshgrid(x,y)
            dx=xx-float(wake_center_xy_m[0]);dy=yy-float(wake_center_xy_m[1])
            downstream=dx*direction[0]+dy*direction[1]-wake_radius_m
            lateral=-dx*direction[1]+dy*direction[0]
            width=wake_base_width_m+wake_spread_slope*np.maximum(downstream,0.)
            shape=np.exp(-np.maximum(downstream,0.)/wake_length_m)
            shape*=np.exp(-0.5*(lateral/np.maximum(width,1e-30))**2)
            shape=np.where(downstream>0.,shape,0.)
            deficit=np.clip(wake_deficit_strength*shape,0.,1.)
            local_wind*=1.-deficit[None,:,:]
            if wake_recirculation_strength>0:
                local_wind-=direction[:,None,None]*(
                    wake_recirculation_strength*abs(float(np.linalg.norm(wind)))
                    *shape[None,:,:])
        if wake_meander_strength>0:
            # A finite obstacle sheds alternating lateral momentum.  The
            # phase is evaluated at retarded time (t - x/U), so the signal
            # is convected downstream instead of rotating the whole wake in
            # lockstep.  This is a carrier-velocity closure: mass/species are
            # still moved only by the conservative FV equations, and the
            # resulting drag/entrainment impulses stay in the momentum
            # ledger below.
            direction=wind/max(float(np.linalg.norm(wind)),1.e-30)
            cross_direction=np.array([-direction[1],direction[0]])
            xx,yy=np.meshgrid(x,y)
            dx=xx-float(wake_center_xy_m[0])
            dy=yy-float(wake_center_xy_m[1])
            downstream=dx*direction[0]+dy*direction[1]-wake_radius_m
            lateral=-dx*direction[1]+dy*direction[0]
            width=wake_base_width_m+wake_spread_slope*np.maximum(downstream,0.)
            envelope=np.exp(-np.maximum(downstream,0.)/wake_length_m)
            envelope*=np.exp(-0.5*(lateral/np.maximum(width,1.e-30))**2)
            envelope=np.where(downstream>0.,envelope,0.)
            obstacle_diameter=max(2.0*float(wake_radius_m),1.e-30)
            current_speed=abs(float(np.linalg.norm(wind)))
            convected_distance=(
                current_speed*float(times[-1]+0.5*dt)-downstream)
            shedding_phase=(
                2.0*math.pi*float(wake_shedding_strouhal)
                *convected_distance/obstacle_diameter
                +math.radians(float(wake_meander_phase_deg)))
            lateral_speed=(
                float(wake_meander_strength)*current_speed*envelope
                *np.sin(shedding_phase))
            local_wind+=cross_direction[:,None,None]*lateral_speed[None,:,:]
        if cavity_backflow_strength > 0.0:
            # A separated recirculation bubble occupies the downstream half
            # of the fence cavity.  Blend a bounded reverse carrier there;
            # the resulting drag and entrainment impulses remain in the
            # momentum ledger, so no receptor-specific concentration gain is
            # introduced.
            radius=max(float(cavity_radius_m),1e-30)
            if cavity_backflow_frame=='legacy_x':
                legacy_dx=xx-float(cavity_center_xy_m[0])
                legacy_dy=yy-float(cavity_center_xy_m[1])
                downstream=np.clip(legacy_dx/radius,0.0,1.0)
                radial=np.clip(
                    np.sqrt(legacy_dx**2+legacy_dy**2)/radius,0.0,1.0)
                bubble=cavity_backflow_strength*downstream*(0.5+0.5*radial)
                target_u=wind[0]*(1.0-2.0*bubble)
                target_v=(-wind[0]*0.5*bubble
                          *np.clip(legacy_dy/radius,-1.0,1.0))
            else:
                downstream=np.clip(
                    cavity_streamwise_coordinate/radius,0.0,1.0)
                radial=np.clip(np.sqrt(
                    cavity_streamwise_coordinate**2
                    +cavity_crosswind_coordinate**2)/radius,0.0,1.0)
                bubble=cavity_backflow_strength*downstream*(0.5+0.5*radial)
                along_speed=float(np.dot(wind,cavity_streamwise_unit))
                crosswind_speed=float(np.dot(wind,cavity_crosswind_unit))
                target_along=along_speed*(1.0-2.0*bubble)
                target_crosswind=(
                    crosswind_speed
                    -along_speed*0.5*bubble*np.clip(
                        cavity_crosswind_coordinate/radius,-1.0,1.0))
                target_u=(target_along*cavity_streamwise_unit[0]
                          +target_crosswind*cavity_crosswind_unit[0])
                target_v=(target_along*cavity_streamwise_unit[1]
                          +target_crosswind*cavity_crosswind_unit[1])
            alpha=np.where(cavity_mask,bubble,0.0)
            local_wind[0]=(1.0-alpha)*local_wind[0]+alpha*target_u
            local_wind[1]=(1.0-alpha)*local_wind[1]+alpha*target_v
        if cavity_stagnation_strength > 0.0:
            # Front-tip stagnation layer: the incoming carrier loses speed
            # before splitting around the finite-height fence.  The Gaussian
            # footprint is tied to the registered obstacle centre/radius and
            # is applied only inside the cavity; it is not a sensor gain.
            radius=max(float(cavity_radius_m),1e-30)
            stagnation=np.exp(
                -0.5*(cavity_streamwise_coordinate/(0.30*radius))**2
                -0.5*(cavity_crosswind_coordinate/(0.70*radius))**2)
            alpha=np.where(cavity_mask,cavity_stagnation_strength*stagnation,0.0)
            local_wind*=(1.0-alpha)[None,:,:]
        if cavity_shelter_strength>0:
            # The shelter parameter also reduces the local carrier wind used
            # by entrainment/drag.  Applying damping only to q would be
            # immediately undone by the ambient-wind relaxation below.
            local_wind[:,cavity_mask]*=(1.0-cavity_shelter_strength)
        slip=local_wind-np.stack((u,v))
        speed=np.sqrt(np.sum(slip**2,axis=0))
        # Only marked cloud surface entrains; never grow gas in a dry cell.
        suppression=1.
        if dense_entrainment_model=='parker_bulk_suppression':
            rho=np.divide(q[0],height,out=np.full_like(height,thermo.ambient_density),where=height>1e-30)
            ri=gravity_scale*9.80665*np.maximum(rho/thermo.ambient_density-1.,0.)*height/np.maximum(speed**2,1e-30)
            suppression=(1+718*ri**2.4)**-.5
        local_entrainment_coefficient=float(entrainment_coefficient)
        if entrainment_far_coefficient is not None:
            # Collapse-generated mixing can be strongest close to the source
            # and relax toward an atmospheric far-field value.  The closure
            # is global in time, uses no receptor information, and entrains
            # only ambient mixture and its momentum exactly as before.
            relax=math.exp(-float(times[-1])/entrainment_relaxation_time_s)
            local_entrainment_coefficient=(
                float(entrainment_far_coefficient)
                +(float(entrainment_coefficient)
                  -float(entrainment_far_coefficient))*relax
            )
        dm=(thermo.ambient_density*local_entrainment_coefficient*speed*dt
            *(q[0]>1e-20)*suppression)
        if cavity_entrainment_factor < 1.0 or cavity_core_entrainment_reduction > 0.0:
            if cavity_core_entrainment_reduction > 0.0 and cavity_mask.any():
                radial_fraction=np.clip(
                    np.sqrt((xx-cx)**2+(yy-cy)**2)
                    / max(float(cavity_radius_m),1.0e-30), 0.0, 1.0,
                )
                core_factor=(1.0-cavity_core_entrainment_reduction)
                local_cavity_factor=cavity_entrainment_factor * (
                    core_factor + (1.0-core_factor)*radial_fraction
                )
                dm=np.where(cavity_mask,dm*local_cavity_factor,dm)
            else:
                dm=np.where(cavity_mask,dm*cavity_entrainment_factor,dm)
        if exterior_entrainment_far_factor < 1.0:
            if cavity_center_xy_m is None or cavity_radius_m <= 0.0:
                raise ValueError(
                    'exterior entrainment relaxation requires cavity geometry')
            radial_distance=np.sqrt((xx-cx)**2+(yy-cy)**2)
            outside_distance=np.maximum(
                radial_distance-float(cavity_radius_m),0.0)
            exterior_factor=(
                float(exterior_entrainment_far_factor)
                +(1.0-float(exterior_entrainment_far_factor))*np.exp(
                    -outside_distance/float(exterior_entrainment_decay_length_m)
                )
            )
            dm=np.where(radial_distance>float(cavity_radius_m),
                        dm*exterior_factor,dm)
        entrainment_momentum_impulse=np.array([
            np.sum(dm*local_wind[0]),np.sum(dm*local_wind[1])])*area
        q[0]+=dm;q[2]+=dm*local_wind[0];q[3]+=dm*local_wind[1]
        # A bounded implicit drag relaxation avoids drag timestep instability.
        momentum_before_drag=q[2:4].sum(axis=(1,2))*area
        drag=.5*drag_coefficient*thermo.ambient_density*speed/np.maximum(q[0],1e-30)
        factor=1/(1+drag*dt)
        q[2]=q[0]*local_wind[0]+factor*(q[2]-q[0]*local_wind[0])
        q[3]=q[0]*local_wind[1]+factor*(q[3]-q[0]*local_wind[1])
        drag_momentum_impulse=q[2:4].sum(axis=(1,2))*area-momentum_before_drag
        if vertical_q is not None:
            if multilayer_entrainment_placement=='legacy_ground':
                # Historical replay only.  This placement contradicted the
                # documented top-entrainment boundary condition and
                # spuriously diluted the lowest sensor layer.
                vertical_q[0,0]+=dm
            else:
                # The bulk closure above represents ambient gas crossing the
                # cloud's free upper surface.  Deposit that mass in the
                # highest materially occupied fixed layer of each column,
                # rather than at the ground.  A relative threshold ignores
                # numerical Gaussian tails while remaining independent of
                # receptors and absolute inventory scale.
                layer_mass=np.maximum(vertical_q[0],0.0)
                column_mass_before=layer_mass.sum(axis=0)
                active=(layer_mass>
                        np.maximum(column_mass_before*1.e-12,1.e-30)[None,:,:])
                layer_number=np.arange(1,vertical_layers+1)[:,None,None]
                top_index=np.maximum(
                    np.max(np.where(active,layer_number,0),axis=0)-1,0)
                row_index=np.arange(ny)[:,None]
                column_index=np.arange(nx)[None,:]
                vertical_q[0,top_index,row_index,column_index]+=dm
            column_mass=vertical_q[0].sum(axis=0,keepdims=True)
            weights=np.divide(vertical_q[0],column_mass,
                              out=np.zeros_like(vertical_q[0]),where=column_mass>1e-30)
            vertical_q[2]=q[2][None,:,:]*weights
            vertical_q[3]=q[3][None,:,:]*weights
            vertical_velocity=np.full((ny,nx),-dense_settling_velocity_m_s,dtype=float)
            if obstacle_vertical_lift_fraction>0.0:
                if cavity_center_xy_m is None or cavity_radius_m<=0.0:
                    raise ValueError('obstacle vertical lift requires a positive cavity geometry')
                direction=wind/max(float(np.linalg.norm(wind)),1.e-30)
                dx=xx-cx;dy=yy-cy
                along=dx*direction[0]+dy*direction[1]
                cross=-dx*direction[1]+dy*direction[0]
                on_span=np.abs(cross)<=float(cavity_radius_m)
                wall_along=np.sqrt(np.maximum(float(cavity_radius_m)**2-cross*cross,0.0))
                downstream=along-wall_along
                # Begin lifting one horizontal cell before the curved wall,
                # then decay through the separated shear layer downstream.
                approach=np.clip((downstream+h)/h,0.0,1.0)
                wake=np.exp(-np.maximum(downstream,0.0)/obstacle_vertical_lift_decay_m)
                lift=(obstacle_vertical_lift_fraction*float(np.linalg.norm(wind))
                      *approach*wake*on_span)
                vertical_velocity+=lift
                maximum_obstacle_vertical_velocity=max(
                    maximum_obstacle_vertical_velocity,float(np.max(lift)))
            if wall_impulse_vertical_lift_fraction>0.0:
                vertical_velocity+=wall_impulse_lift_velocity
            if np.any(vertical_velocity!=0.0):
                advective_flux=np.zeros(
                    (4,vertical_layers+1,ny,nx),dtype=float)
                for level in range(1,vertical_layers):
                    advective_flux[:,level]=np.where(
                        vertical_velocity[None,:,:]>=0.0,
                        vertical_velocity[None,:,:]*vertical_q[:,level-1],
                        vertical_velocity[None,:,:]*vertical_q[:,level])
                if wall_impulse_vertical_lift_fraction>0.0:
                    wall_impulse_lifted_species_kg+=float(np.sum(
                        wall_impulse_lift_velocity*vertical_q[1,0]
                    )*dt*area/vertical_dz)
                vertical_q-=dt/vertical_dz*(
                    advective_flux[:,1:]-advective_flux[:,:-1])
            vertical_flux=np.zeros((4,vertical_layers+1,ny,nx),dtype=float)
            vertical_flux[:,1:-1]=(-vertical_mixing_diffusivity_m2_s
                *(vertical_q[:,1:]-vertical_q[:,:-1])/vertical_dz)
            vertical_q-=dt/vertical_dz*(vertical_flux[:,1:]-vertical_flux[:,:-1])
            if np.any(~np.isfinite(vertical_q)) or np.min(vertical_q[:2])<-1e-10:
                raise RuntimeError('multilayer vertical positivity gate failed')
            vertical_q,redistributed=_repair_multilayer_species_bounds(vertical_q)
            multilayer_species_bound_redistributed+=redistributed*area
            multilayer_relative_transport_pending_s+=dt
            relative_transport_due=(
                multilayer_relative_transport_interval_s<=0.0
                or multilayer_relative_transport_pending_s+1.e-12
                   >=multilayer_relative_transport_interval_s
                or times[-1]+dt>=duration_s-1.e-12)
            if ((multilayer_horizontal_shear_fraction>0.0
                    or multilayer_obstacle_crosswind_shear_rate_s_inv!=0.0)
                    and relative_transport_due):
                # Resolve the horizontal carrier about the column mean.  The
                # mass-weighted relative velocity is zero in every cell, so
                # this adds vertical shear without changing the declared
                # depth-mean wind or using receptor observations.
                layer_mass=np.maximum(vertical_q[0],0.0)
                total_layer_mass=layer_mass.sum(axis=0)
                z_centres=(.5*(z_edges[:-1]+z_edges[1:]))[:,None,None]
                log_shape=np.log1p(z_centres/max(float(roughness_length_m),1.e-30))
                mean_log_shape=np.divide(
                    np.sum(layer_mass*log_shape,axis=0),total_layer_mass,
                    out=np.zeros_like(total_layer_mass),where=total_layer_mass>1.e-30)
                local_speed=np.hypot(local_wind[0],local_wind[1])
                unit_x=np.divide(
                    local_wind[0],local_speed,out=np.zeros_like(local_speed),
                    where=local_speed>1.e-30)
                unit_y=np.divide(
                    local_wind[1],local_speed,out=np.zeros_like(local_speed),
                    where=local_speed>1.e-30)
                streamwise_relative=np.zeros_like(layer_mass)
                if multilayer_horizontal_shear_fraction>0.0:
                    layer_scale=np.divide(
                        log_shape,mean_log_shape[None,:,:],
                        out=np.ones_like(layer_mass),
                        where=mean_log_shape[None,:,:]>1.e-30)
                    streamwise_relative=(
                        float(multilayer_horizontal_shear_fraction)
                        *local_speed[None,:,:]*(layer_scale-1.0))
                    streamwise_relative=np.where(
                        total_layer_mass[None,:,:]>1.e-30,
                        streamwise_relative,0.0)
                cross_relative=np.zeros_like(layer_mass)
                if multilayer_obstacle_crosswind_shear_rate_s_inv!=0.0:
                    mean_z=np.divide(
                        np.sum(layer_mass*z_centres,axis=0),total_layer_mass,
                        out=np.zeros_like(total_layer_mass),
                        where=total_layer_mass>1.e-30)
                    wind_direction=wind/max(float(np.linalg.norm(wind)),1.e-30)
                    grid_x,grid_y=np.meshgrid(x,y)
                    dx=grid_x-float(cavity_center_xy_m[0])
                    dy=grid_y-float(cavity_center_xy_m[1])
                    along=dx*wind_direction[0]+dy*wind_direction[1]
                    cross=-dx*wind_direction[1]+dy*wind_direction[0]
                    wall_along=np.sqrt(np.maximum(
                        float(cavity_radius_m)**2-cross*cross,0.0))
                    downstream=along-wall_along
                    wake_weight=(
                        (np.abs(cross)<=float(cavity_radius_m))
                        &(downstream>=0.0))
                    wake_weight=(wake_weight*np.exp(
                        -np.maximum(downstream,0.0)
                        /float(multilayer_obstacle_crosswind_decay_length_m)))
                    cross_relative=(
                        float(multilayer_obstacle_crosswind_shear_rate_s_inv)
                        *(z_centres-mean_z[None,:,:])
                        *wake_weight[None,:,:])
                # Positive cross-relative speed follows the meteorological
                # left-normal of the local carrier direction.
                relative_u=(streamwise_relative*unit_x[None,:,:]
                            -cross_relative*unit_y[None,:,:])
                relative_v=(streamwise_relative*unit_y[None,:,:]
                            +cross_relative*unit_x[None,:,:])
                maximum_multilayer_relative_velocity=max(
                    maximum_multilayer_relative_velocity,
                    float(np.max(np.hypot(relative_u,relative_v))))
                vertical_q,moved,substeps=_transport_multilayer_relative(
                    vertical_q,relative_u,relative_v,
                    dt_s=multilayer_relative_transport_pending_s,
                    cell_size_m=h,wall_x_height_m=wx,
                    wall_y_height_m=wy,z_edges_m=z_edges)
                multilayer_relative_transport_moved+=moved
                multilayer_relative_transport_substeps+=substeps
                multilayer_relative_transport_update_count+=1
                multilayer_relative_transport_pending_s=0.0
            momentum_before_vertical=q[2:4].sum(axis=(1,2))*area
            q=vertical_q.sum(axis=1)
            maximum_vertical_momentum_residual=max(
                maximum_vertical_momentum_residual,
                float(np.max(np.abs(q[2:4].sum(axis=(1,2))*area-momentum_before_vertical))))
        if elevated_species is not None:
            elevated_transport_pending_s+=dt
            bounded=np.clip(elevated_species,0.0,np.maximum(q[1],0.0))
            bound_scale=np.divide(
                bounded,elevated_species,out=np.zeros_like(bounded),
                where=elevated_species>1.e-30)
            elevated_species_age_moment*=bound_scale
            wall_flux_elevated_species*=bound_scale
            wall_flux_elevated_age_moment*=bound_scale
            elevated_species_subset_correction+=float(
                np.abs(bounded-elevated_species).sum()*area)
            elevated_species=bounded
            transport_due=(
                wall_impulse_elevated_transport_interval_s<=0.0
                or elevated_transport_pending_s+1.e-12
                >=wall_impulse_elevated_transport_interval_s
                or times[-1]+dt>=duration_s-1.e-12)
            marker_transport_dt=elevated_transport_pending_s
            if ((wall_impulse_elevated_transport_fraction>0.0
                    or wall_impulse_elevated_crosswind_velocity_m_s!=0.0
                    or wall_impulse_elevated_crosswind_diffusivity_m2_s>0.0
                    or wall_impulse_elevated_crosswind_outward_velocity_m_s>0.0
                    or wall_impulse_elevated_advection_scale!=1.0)
                    and transport_due):
                _,branch_u,branch_v,_,_,_=_primitive(q,thermo,gravity_scale)
                relative_u=(float(wall_impulse_elevated_transport_fraction)
                            *(local_wind[0]-branch_u))
                relative_v=(float(wall_impulse_elevated_transport_fraction)
                            *(local_wind[1]-branch_v))
                local_age=np.divide(
                    elevated_species_age_moment,elevated_species,
                    out=np.zeros_like(elevated_species),
                    where=elevated_species>1.e-30)
                if wall_impulse_elevated_advection_scale!=1.0:
                    active=(local_age
                            <float(wall_impulse_elevated_advection_duration_s))
                    scale=np.where(
                        active,float(wall_impulse_elevated_advection_scale),1.0)
                    relative_u+=(scale-1.0)*branch_u
                    relative_v+=(scale-1.0)*branch_v
                local_speed=np.hypot(local_wind[0],local_wind[1])
                crosswind=float(wall_impulse_elevated_crosswind_velocity_m_s)
                if wall_impulse_elevated_crosswind_outward_velocity_m_s>0.0:
                    cross_x=np.divide(
                        -local_wind[1],local_speed,
                        out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                    cross_y=np.divide(
                        local_wind[0],local_speed,
                        out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                    cross_coordinate=(
                        (xx-float(cx))*cross_x+(yy-float(cy))*cross_y)
                    outward=(
                        float(wall_impulse_elevated_crosswind_outward_velocity_m_s)
                        *np.tanh(cross_coordinate/float(
                            wall_impulse_elevated_crosswind_outward_scale_m)))
                    if wall_impulse_elevated_crosswind_outward_decay_s>0.0:
                        outward*=np.exp(-local_age/float(
                            wall_impulse_elevated_crosswind_outward_decay_s))
                    crosswind=crosswind+outward
                if wall_impulse_elevated_crosswind_spinup_s>0.0:
                    crosswind=(crosswind*np.clip(
                        local_age/float(wall_impulse_elevated_crosswind_spinup_s),
                        0.0,1.0))
                relative_u+=np.divide(
                    -crosswind*local_wind[1],local_speed,
                    out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                relative_v+=np.divide(
                    crosswind*local_wind[0],local_speed,
                    out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                elevated_species_relative_transport+=transport_elevated_relative(
                    relative_u,relative_v,elevated_transport_pending_s)
                elevated_species_crosswind_diffusion+=diffuse_elevated_crosswind(
                    elevated_transport_pending_s,local_wind)
                elevated_transport_pending_s=0.0
                elevated_transport_update_count+=1
            if (wall_flux_elevated_transfer_fraction>0.0
                    and np.any(wall_flux_elevated_species>0.0)
                    and transport_due):
                _,branch_u,branch_v,_,_,_=_primitive(q,thermo,gravity_scale)
                local_speed=np.hypot(local_wind[0],local_wind[1])
                direction_x=np.divide(
                    local_wind[0],local_speed,out=np.zeros_like(local_speed),
                    where=local_speed>1.e-30)
                direction_y=np.divide(
                    local_wind[1],local_speed,out=np.zeros_like(local_speed),
                    where=local_speed>1.e-30)
                if friction_velocity_m_s is not None:
                    wall_reference_height=max(
                        float(np.max(wx)),float(np.max(wy)),
                        float(roughness_length_m))
                    upper_speed=(float(friction_velocity_m_s)/.4
                                 *math.log1p(
                                     wall_reference_height
                                     /float(roughness_length_m)))
                    depth_mean=np.maximum(mean,1.e-30)
                    history_scale=np.divide(
                        local_speed,depth_mean,out=np.ones_like(local_speed),
                        where=depth_mean>1.e-30)
                    upper_speed=upper_speed*history_scale
                else:
                    upper_speed=local_speed
                upper_u=direction_x*upper_speed
                upper_v=direction_y*upper_speed
                wall_flux_elevated_relative_transport+=(
                    transport_wall_flux_relative(
                        upper_u-branch_u,upper_v-branch_v,
                        marker_transport_dt))
                elevated_transport_pending_s=0.0
            if wall_impulse_elevated_return_rate_s_inv>0.0:
                return_factor=1.0-math.exp(
                    -float(wall_impulse_elevated_return_rate_s_inv)*dt)
                returned_elevated=elevated_species*return_factor
                returned_elevated_age=(
                    elevated_species_age_moment*return_factor)
                elevated_species-=returned_elevated
                elevated_species_age_moment-=returned_elevated_age
                wall_flux_elevated_species*=(1.0-return_factor)
                wall_flux_elevated_age_moment*=(1.0-return_factor)
                elevated_species_returned+=float(returned_elevated.sum()*area)
            # Age every member of the elevated branch, including material in
            # a sub-grid cavity store.  Newly transferred ground material is
            # added below with zero age, so age is measured from wall impact.
            elevated_species_age_moment+=elevated_species*dt
            wall_flux_elevated_age_moment+=(
                wall_flux_elevated_species*dt)
            local_cavity_reservoir_elevated_age_moment+=(
                local_cavity_reservoir_elevated_species*dt)
            local_cavity_reservoir_wall_flux_age_moment+=(
                local_cavity_reservoir_wall_flux_species*dt)
            cavity_reservoir_elevated_age_moment+=(
                cavity_reservoir_elevated_species*dt)
            cavity_reservoir_wall_flux_age_moment+=(
                cavity_reservoir_wall_flux_species*dt)
            wall_overflow_elevated_age_x+=wall_overflow_elevated_x*dt
            wall_overflow_elevated_age_y+=wall_overflow_elevated_y*dt
            wall_overflow_flux_age_x+=wall_overflow_flux_marker_x*dt
            wall_overflow_flux_age_y+=wall_overflow_flux_marker_y*dt
            transfer_factor=1.0-np.exp(
                -float(wall_impulse_elevated_transfer_fraction)
                *wall_impulse_approach_speed*dt/h)
            ground_species=np.maximum(q[1]-elevated_species,0.0)
            transferred=ground_species*transfer_factor
            elevated_species+=transferred
            elevated_species_transferred+=float(transferred.sum()*area)
            if (np.min(elevated_species)<-1.e-10
                    or np.max(elevated_species-q[1])>1.e-10):
                raise RuntimeError('wall elevated species transfer violated subset bounds')
        if (wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                and wall_impulse_elevated_transfer_fraction>0.0):
            if wall_impulse_reattached_fraction>0.0:
                conserved_transport_pending_s+=dt
                prompt_transport_pending_s=conserved_transport_pending_s
            else:
                elevated_transport_pending_s+=dt
                prompt_transport_pending_s=elevated_transport_pending_s
            transport_due=(
                wall_impulse_elevated_transport_interval_s<=0.0
                or prompt_transport_pending_s+1.e-12
                >=wall_impulse_elevated_transport_interval_s
                or times[-1]+dt>=duration_s-1.e-12)
            if np.any(wall_impulse_elevated_reservoir>0.0) and transport_due:
                prompt_transport_dt=prompt_transport_pending_s
                if (wall_impulse_reattached_fraction>0.0
                        and wall_impulse_elevated_transport_interval_s>0.0
                        and times[-1]+dt<duration_s-1.e-12):
                    prompt_transport_dt=float(
                        wall_impulse_elevated_transport_interval_s)
                prompt_mass=wall_impulse_elevated_reservoir[0]
                prompt_u=np.divide(
                    wall_impulse_elevated_reservoir[2],prompt_mass,
                    out=np.zeros_like(prompt_mass),where=prompt_mass>1.e-30)
                prompt_v=np.divide(
                    wall_impulse_elevated_reservoir[3],prompt_mass,
                    out=np.zeros_like(prompt_mass),where=prompt_mass>1.e-30)
                ambient_fraction=conserved_prompt_transport_fraction
                carrier_u=((1.0-ambient_fraction)*prompt_u
                           +ambient_fraction*local_wind[0])
                carrier_v=((1.0-ambient_fraction)*prompt_v
                           +ambient_fraction*local_wind[1])
                prompt_local_age=np.divide(
                    wall_impulse_elevated_reservoir_age_moment,
                    wall_impulse_elevated_reservoir[1],
                    out=np.zeros_like(wall_impulse_elevated_reservoir[1]),
                    where=wall_impulse_elevated_reservoir[1]>1.e-30)
                if wall_impulse_elevated_advection_scale!=1.0:
                    active=(prompt_local_age
                            <float(wall_impulse_elevated_advection_duration_s))
                    scale=np.where(
                        active,float(wall_impulse_elevated_advection_scale),1.0)
                    carrier_u*=scale;carrier_v*=scale
                local_speed=np.hypot(local_wind[0],local_wind[1])
                crosswind=conserved_prompt_crosswind_velocity
                if conserved_prompt_outward_velocity>0.0:
                    cross_x=np.divide(
                        -local_wind[1],local_speed,out=np.zeros_like(local_speed),
                        where=local_speed>1.e-30)
                    cross_y=np.divide(
                        local_wind[0],local_speed,out=np.zeros_like(local_speed),
                        where=local_speed>1.e-30)
                    cross_coordinate=(
                        (xx-float(cx))*cross_x+(yy-float(cy))*cross_y)
                    outward=(
                        conserved_prompt_outward_velocity
                        *np.tanh(cross_coordinate/conserved_prompt_outward_scale))
                    if wall_impulse_elevated_crosswind_outward_decay_s>0.0:
                        outward*=np.exp(-prompt_local_age/float(
                            wall_impulse_elevated_crosswind_outward_decay_s))
                    crosswind=crosswind+outward
                if wall_impulse_elevated_crosswind_spinup_s>0.0:
                    crosswind=(crosswind*np.clip(
                        prompt_local_age/float(
                            wall_impulse_elevated_crosswind_spinup_s),0.0,1.0))
                carrier_u+=np.divide(
                    -crosswind*local_wind[1],local_speed,
                    out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                carrier_v+=np.divide(
                    crosswind*local_wind[0],local_speed,
                    out=np.zeros_like(local_speed),where=local_speed>1.e-30)
                transported,moved,_=_transport_multilayer_relative(
                    wall_impulse_elevated_reservoir[:,None,:,:],
                    carrier_u[None,:,:],carrier_v[None,:,:],
                    dt_s=prompt_transport_dt,cell_size_m=h,
                    wall_x_height_m=np.zeros_like(wx),
                    wall_y_height_m=np.zeros_like(wy),
                    z_edges_m=np.array([0.0,1.0]))
                wall_impulse_elevated_reservoir=transported[:,0]
                age_state=np.zeros_like(transported)
                age_state[0,0]=wall_impulse_elevated_reservoir_age_moment
                age_state[1,0]=wall_impulse_elevated_reservoir_age_moment
                age_state,_,_=_transport_multilayer_relative(
                    age_state,carrier_u[None,:,:],carrier_v[None,:,:],
                    dt_s=prompt_transport_dt,cell_size_m=h,
                    wall_x_height_m=np.zeros_like(wx),
                    wall_y_height_m=np.zeros_like(wy),
                    z_edges_m=np.array([0.0,1.0]))
                wall_impulse_elevated_reservoir_age_moment=age_state[0,0]
                if wall_impulse_elevated_conserved_center_mode=='impact_head':
                    center_state=np.zeros_like(transported)
                    center_state[0,0]=(
                        wall_impulse_elevated_reservoir_center_moment)
                    center_state[1,0]=(
                        wall_impulse_elevated_reservoir_center_moment)
                    center_state,_,_=_transport_multilayer_relative(
                        center_state,carrier_u[None,:,:],carrier_v[None,:,:],
                        dt_s=prompt_transport_dt,cell_size_m=h,
                        wall_x_height_m=np.zeros_like(wx),
                        wall_y_height_m=np.zeros_like(wy),
                        z_edges_m=np.array([0.0,1.0]))
                    wall_impulse_elevated_reservoir_center_moment=(
                        center_state[0,0])
                wall_impulse_elevated_reservoir_transport_moved+=moved
                if wall_impulse_reattached_fraction>0.0:
                    conserved_transport_pending_s=max(
                        conserved_transport_pending_s-prompt_transport_dt,0.0)
                else:
                    elevated_transport_pending_s=0.0
                elevated_transport_update_count+=1
            if wall_impulse_elevated_return_rate_s_inv>0.0:
                return_factor=1.0-math.exp(
                    -float(wall_impulse_elevated_return_rate_s_inv)*dt)
                returned=wall_impulse_elevated_reservoir*return_factor
                returned_age=(wall_impulse_elevated_reservoir_age_moment
                              *return_factor)
                wall_impulse_elevated_reservoir-=returned
                wall_impulse_elevated_reservoir_age_moment-=returned_age
                if wall_impulse_elevated_conserved_center_mode=='impact_head':
                    wall_impulse_elevated_reservoir_center_moment*=(
                        1.0-return_factor)
                q+=returned
                wall_impulse_elevated_returned+=(
                    returned.sum(axis=(1,2))*area)
                elevated_species_returned+=float(returned[1].sum()*area)
            if wall_impulse_reattached_fraction>0.0:
                reattached_transport_pending_s+=dt
                reattached_due=(
                    wall_impulse_elevated_transport_interval_s<=0.0
                    or reattached_transport_pending_s+1.e-12
                    >=wall_impulse_elevated_transport_interval_s
                    or times[-1]+dt>=duration_s-1.e-12)
                if (np.any(wall_impulse_reattached_reservoir>0.0)
                        and reattached_due):
                    reattached_transport_dt=reattached_transport_pending_s
                    if (wall_impulse_elevated_transport_interval_s>0.0
                            and times[-1]+dt<duration_s-1.e-12):
                        reattached_transport_dt=float(
                            wall_impulse_elevated_transport_interval_s)
                    lower_mass=wall_impulse_reattached_reservoir[0]
                    lower_u=np.divide(
                        wall_impulse_reattached_reservoir[2],lower_mass,
                        out=np.zeros_like(lower_mass),where=lower_mass>1.e-30)
                    lower_v=np.divide(
                        wall_impulse_reattached_reservoir[3],lower_mass,
                        out=np.zeros_like(lower_mass),where=lower_mass>1.e-30)
                    ambient_fraction=float(
                        wall_impulse_reattached_transport_fraction)
                    reattached_ambient_u=local_wind[0]
                    reattached_ambient_v=local_wind[1]
                    if wall_impulse_reattached_transport_height_m is not None:
                        local_speed=np.hypot(local_wind[0],local_wind[1])
                        reference_speed=max(float(np.linalg.norm(wind)),1.e-30)
                        height_speed=(
                            float(friction_velocity_m_s)/0.41
                            *math.log1p(float(
                                wall_impulse_reattached_transport_height_m)
                                /float(roughness_length_m)))
                        scale=height_speed/reference_speed
                        reattached_ambient_u=local_wind[0]*scale
                        reattached_ambient_v=local_wind[1]*scale
                    carrier_u=((1.0-ambient_fraction)*lower_u
                               +ambient_fraction*reattached_ambient_u)
                    carrier_v=((1.0-ambient_fraction)*lower_v
                               +ambient_fraction*reattached_ambient_v)
                    local_speed=np.hypot(local_wind[0],local_wind[1])
                    crosswind=float(
                        wall_impulse_reattached_crosswind_velocity_m_s)
                    if (wall_impulse_reattached_crosswind_outward_velocity_m_s
                            >0.0):
                        cross_x=np.divide(
                            -local_wind[1],local_speed,
                            out=np.zeros_like(local_speed),
                            where=local_speed>1.e-30)
                        cross_y=np.divide(
                            local_wind[0],local_speed,
                            out=np.zeros_like(local_speed),
                            where=local_speed>1.e-30)
                        cross_coordinate=(
                            (xx-float(cx))*cross_x+(yy-float(cy))*cross_y)
                        crosswind=(crosswind
                            +float(
                                wall_impulse_reattached_crosswind_outward_velocity_m_s)
                            *np.tanh(cross_coordinate/float(
                                wall_impulse_reattached_crosswind_outward_scale_m)))
                    carrier_u+=np.divide(
                        -crosswind*local_wind[1],local_speed,
                        out=np.zeros_like(local_speed),
                        where=local_speed>1.e-30)
                    carrier_v+=np.divide(
                        crosswind*local_wind[0],local_speed,
                        out=np.zeros_like(local_speed),
                        where=local_speed>1.e-30)
                    lower_transport,moved,_=_transport_multilayer_relative(
                        wall_impulse_reattached_reservoir[:,None,:,:],
                        carrier_u[None,:,:],carrier_v[None,:,:],
                        dt_s=reattached_transport_dt,cell_size_m=h,
                        wall_x_height_m=np.zeros_like(wx),
                        wall_y_height_m=np.zeros_like(wy),
                        z_edges_m=np.array([0.0,1.0]))
                    wall_impulse_reattached_reservoir=lower_transport[:,0]
                    lower_age_state=np.zeros_like(lower_transport)
                    lower_age_state[0,0]=(
                        wall_impulse_reattached_reservoir_age_moment)
                    lower_age_state[1,0]=(
                        wall_impulse_reattached_reservoir_age_moment)
                    lower_age_state,_,_=_transport_multilayer_relative(
                        lower_age_state,carrier_u[None,:,:],carrier_v[None,:,:],
                        dt_s=reattached_transport_dt,cell_size_m=h,
                        wall_x_height_m=np.zeros_like(wx),
                        wall_y_height_m=np.zeros_like(wy),
                        z_edges_m=np.array([0.0,1.0]))
                    wall_impulse_reattached_reservoir_age_moment=(
                        lower_age_state[0,0])
                    wall_impulse_reattached_reservoir_transport_moved+=moved
                    reattached_transport_pending_s=max(
                        reattached_transport_pending_s-reattached_transport_dt,
                        0.0)
                if wall_impulse_reattached_return_rate_s_inv>0.0:
                    lower_return_factor=1.0-math.exp(
                        -float(wall_impulse_reattached_return_rate_s_inv)*dt)
                    lower_returned=(wall_impulse_reattached_reservoir
                                    *lower_return_factor)
                    lower_returned_age=(
                        wall_impulse_reattached_reservoir_age_moment
                        *lower_return_factor)
                    wall_impulse_reattached_reservoir-=lower_returned
                    wall_impulse_reattached_reservoir_age_moment-=(
                        lower_returned_age)
                    q+=lower_returned
                    wall_impulse_elevated_returned+=(
                        lower_returned.sum(axis=(1,2))*area)
                    elevated_species_returned+=float(
                        lower_returned[1].sum()*area)
            transfer_factor=1.0-np.exp(
                -float(wall_impulse_elevated_transfer_fraction)
                *wall_impulse_approach_speed*dt/h)
            if wall_impulse_elevated_state_mode=='hybrid':
                transfer_factor*=float(
                    wall_impulse_elevated_conserved_fraction)
            captured=q*transfer_factor[None,:,:]
            if elevated_species is not None:
                # The marker is a subset label on q, not extra mass.  Remove
                # precisely the labelled portion transferred to the conserved
                # branch so prompt-path species is never double counted.
                marker_captured=elevated_species*transfer_factor
                marker_age_captured=elevated_species_age_moment*transfer_factor
                elevated_species-=marker_captured
                elevated_species_age_moment-=marker_age_captured
                wall_flux_elevated_species*=(1.0-transfer_factor)
                wall_flux_elevated_age_moment*=(1.0-transfer_factor)
            q-=captured
            # A shallow layer blocked by a tall face primarily feeds the
            # separated lower reattachment path; material whose local layer
            # depth exceeds the adjacent face retains a larger upper-overflow
            # share.  This geometry/field ratio is common to all receptors
            # and partitions already captured conserved inventory only.
            reattached_partition=(
                float(wall_impulse_reattached_fraction)
                *np.clip(np.divide(
                    cell_wall_height,np.maximum(height,1.e-30),
                    out=np.zeros_like(height),where=height>1.e-30),0.0,1.0))
            reattached_captured=captured*reattached_partition[None,:,:]
            wall_impulse_reattached_reservoir+=reattached_captured
            upper_captured=captured-reattached_captured
            wall_impulse_elevated_reservoir+=upper_captured
            wall_impulse_reattached_captured+=(
                reattached_captured.sum(axis=(1,2))*area)
            wall_impulse_elevated_captured+=(
                captured.sum(axis=(1,2))*area)
            elevated_species_transferred+=float(captured[1].sum()*area)
            wall_impulse_elevated_reservoir_age_moment+=(
                wall_impulse_elevated_reservoir[1]*dt)
            if wall_impulse_elevated_conserved_center_mode=='impact_head':
                impact_ratio=np.divide(
                    wall_impulse_approach_speed,np.maximum(c,1.e-30),
                    out=np.zeros_like(c),where=c>1.e-30)
                impact_center=np.minimum(
                    cell_wall_height,
                    0.5*height*(1.0+impact_ratio*impact_ratio))
                impact_center=np.maximum(
                    impact_center,conserved_prompt_min_center)
                wall_impulse_elevated_reservoir_center_moment+=(
                    upper_captured[1]*impact_center)
            wall_impulse_reattached_reservoir_age_moment+=(
                wall_impulse_reattached_reservoir[1]*dt)
        pending_total=(np.sum(np.stack([packet for _,packet in cavity_reservoir_pending],axis=0),axis=0)
                       if cavity_reservoir_pending else np.zeros(4))
        overflow_total=wall_overflow_total()
        momentum_after_step=(q[2:4].sum(axis=(1,2))*area
                             +cavity_reservoir_total()[2:4]+pending_total[2:4]
                             +overflow_total[2:4])
        expected_after_step=(momentum_after_transport
                             +entrainment_momentum_impulse
                             +drag_momentum_impulse)
        maximum_momentum_budget_residual=max(
            maximum_momentum_budget_residual,
            float(np.max(np.abs(momentum_after_step-expected_after_step))))
        exported+=totalout;mixture_exported+=mixtureout;entrained+=float(dm.sum()*area)
        pending_total=(np.sum(np.stack([packet for _,packet in cavity_reservoir_pending],axis=0),axis=0)
                       if cavity_reservoir_pending else np.zeros(4))
        overflow_total=wall_overflow_total()
        reservoir_total=cavity_reservoir_total()
        times.append(times[-1]+dt);mass.append(float(q[0].sum()*area+reservoir_total[0]+pending_total[0]+overflow_total[0]));species.append(float(q[1].sum()*area+reservoir_total[1]+pending_total[1]+overflow_total[1]))
        maxres=max(maxres,abs(species[-1]+exported
                              -exterior_reentry_returned[1]
                              +interface_exported[1]-species[0]
                              -interface_injected[1]-source_injected[1]))
        minq=min(minq,float(q[:2].min()));observe()
    component_result=({label:np.array(history).T for label,history in
                       vertical_observation_component_traces.items()}
                      if vertical_observation_component_diagnostics else None)
    return DenseLayerResult(np.array(times),np.array(traces).T,np.array(packets),
        np.array(current_speeds),np.array(heights),boundary_positions,
        np.array(mass),np.array(species),q,{
            'runtime_s':time.perf_counter()-started,'steps':len(times)-1,
            'maximum_species_mass_residual_kg':maxres,
            'maximum_relative_species_mass_residual':maxres/max(species[0],1e-30),
            'momentum_initial_kg_m_s':momentum_initial.tolist(),
            'momentum_final_kg_m_s':(
                q[2:4].sum(axis=(1,2))*area+cavity_reservoir_total()[2:4]
                +pending_total[2:4]+overflow_total[2:4]).tolist(),
            'maximum_momentum_budget_residual_kg_m_s':maximum_momentum_budget_residual,
            'maximum_vertical_momentum_residual_kg_m_s':maximum_vertical_momentum_residual,
            'momentum_conservation_gate':bool(maximum_momentum_budget_residual<1e-7
                                               and maximum_vertical_momentum_residual<1e-7),
            'mixture_mass_residual_kg':(
                mass[-1]+mixture_exported-exterior_reentry_returned[0]
                +interface_exported[0]
                -mass[0]-entrained-interface_injected[0]-source_injected[0]
            ),
            'interface_x_m':None if interface_x_m is None else float(interface_x_m),
            'interface_time_offset_s':float(interface_time_offset_s),
            'interface_exported_total_kg':float(interface_exported[0]),
            'interface_exported_species_kg':float(interface_exported[1]),
            'interface_injected_total_kg':float(interface_injected[0]),
            'interface_injected_species_kg':float(interface_injected[1]),
            'interface_negative_mass_kg':float(interface_negative_mass),
            'source_injected_total_kg':float(source_injected[0]),
            'source_injected_species_kg':float(source_injected[1]),
            'source_injected_momentum_kg_m_s':source_injected[2:4].tolist(),
            'largest_transport_budget_detail':largest_transport_budget_detail,
            'entrained_air_kg':entrained,'exported_species_mass_kg':exported,
            'vertical_observation_profile':vertical_observation_profile,
            'dense_vertical_exponent':float(dense_vertical_exponent),
            'vertical_grounding_exponent':float(vertical_grounding_exponent),
            'vertical_elevated_fraction':float(vertical_elevated_fraction),
            'wall_overtop_exponent':float(wall_overtop_exponent),
            'wall_blockage_fraction':float(wall_blockage_fraction),
            'wall_overflow_capture_fraction':float(wall_overflow_capture_fraction),
            'wall_overflow_return_rate_s_inv':float(
                wall_overflow_return_rate_s_inv),
            'wall_overflow_side_gate_length_m':float(
                wall_overflow_side_gate_length_m),
            'wall_overflow_side_gate_face_count':int(
                np.count_nonzero(wall_overflow_side_mask_x)
                +np.count_nonzero(wall_overflow_side_mask_y)),
            'wall_overflow_captured_total_kg':float(wall_overflow_captured[0]),
            'wall_overflow_captured_species_kg':float(wall_overflow_captured[1]),
            'wall_overflow_top_captured_species_kg':float(
                wall_overflow_top_captured[1]),
            'wall_overflow_side_captured_species_kg':float(
                wall_overflow_side_captured[1]),
            'wall_overflow_returned_total_kg':float(wall_overflow_returned[0]),
            'wall_overflow_returned_species_kg':float(wall_overflow_returned[1]),
            'wall_overflow_stored_total_kg':float(overflow_total[0]),
            'wall_overflow_stored_species_kg':float(overflow_total[1]),
            'wall_overflow_elevated_stored_species_kg':float(
                wall_overflow_elevated_x.sum()
                +wall_overflow_elevated_y.sum()),
            'wall_overflow_elevated_age_moment_kg_s':float(
                wall_overflow_elevated_age_x.sum()
                +wall_overflow_elevated_age_y.sum()),
            'wall_overflow_flux_marker_stored_species_kg':float(
                wall_overflow_flux_marker_x.sum()
                +wall_overflow_flux_marker_y.sum()),
            'vertical_tail_scale_m':float(vertical_tail_scale_m),
            'vertical_tail_extension_m':float(vertical_tail_extension_m),
            'vertical_tail_diffusivity_m2_s':float(vertical_tail_diffusivity_m2_s),
            'mass_conserving_vertical_profile':mass_conserving_vertical_profile,
            'vertical_observation_mode':vertical_observation_mode,
            'vertical_observation_bins':vertical_observation_bins,
            'vertical_observation_component_diagnostics':vertical_observation_component_diagnostics,
            'vertical_observation_maximum_relative_volume_residual':(
                vertical_observation_maximum_relative_volume_residual),
            'vertical_observation_saturated_column_samples':(
                vertical_observation_saturated_column_samples),
            'vertical_layers':vertical_layers,
            'multilayer_entrainment_placement':multilayer_entrainment_placement,
            'vertical_mixing_diffusivity_m2_s':float(vertical_mixing_diffusivity_m2_s),
            'multilayer_horizontal_shear_fraction':float(
                multilayer_horizontal_shear_fraction),
            'multilayer_obstacle_crosswind_shear_rate_s_inv':float(
                multilayer_obstacle_crosswind_shear_rate_s_inv),
            'multilayer_obstacle_crosswind_decay_length_m':float(
                multilayer_obstacle_crosswind_decay_length_m),
            'multilayer_relative_transport_interval_s':float(
                multilayer_relative_transport_interval_s),
            'maximum_multilayer_relative_velocity_m_s':float(
                maximum_multilayer_relative_velocity),
            'multilayer_relative_transport_moved_kg':float(
                multilayer_relative_transport_moved),
            'multilayer_relative_transport_update_count':int(
                multilayer_relative_transport_update_count),
            'multilayer_relative_transport_substeps':int(
                multilayer_relative_transport_substeps),
            'multilayer_species_bound_redistributed_kg':float(
                multilayer_species_bound_redistributed),
            'obstacle_vertical_lift_fraction':float(obstacle_vertical_lift_fraction),
            'obstacle_vertical_lift_decay_m':float(obstacle_vertical_lift_decay_m),
            'wall_impulse_vertical_lift_fraction':float(
                wall_impulse_vertical_lift_fraction),
            'wall_impulse_elevated_transfer_fraction':float(
                wall_impulse_elevated_transfer_fraction),
            'wall_impulse_elevated_state_mode':wall_impulse_elevated_state_mode,
            'wall_impulse_elevated_conserved_fraction':float(
                wall_impulse_elevated_conserved_fraction),
            'wall_flux_elevated_transfer_fraction':float(
                wall_flux_elevated_transfer_fraction),
            'wall_impulse_elevated_transport_fraction':float(
                wall_impulse_elevated_transport_fraction),
            'wall_impulse_elevated_crosswind_velocity_m_s':float(
                wall_impulse_elevated_crosswind_velocity_m_s),
            'wall_impulse_elevated_crosswind_diffusivity_m2_s':float(
                wall_impulse_elevated_crosswind_diffusivity_m2_s),
            'wall_impulse_elevated_crosswind_outward_velocity_m_s':float(
                wall_impulse_elevated_crosswind_outward_velocity_m_s),
            'wall_impulse_elevated_crosswind_outward_scale_m':float(
                wall_impulse_elevated_crosswind_outward_scale_m),
            'wall_impulse_elevated_crosswind_outward_decay_s':float(
                wall_impulse_elevated_crosswind_outward_decay_s),
            'wall_impulse_elevated_crosswind_spinup_s':float(
                wall_impulse_elevated_crosswind_spinup_s),
            'wall_impulse_elevated_transport_interval_s':float(
                wall_impulse_elevated_transport_interval_s),
            'wall_impulse_elevated_advection_scale':float(
                wall_impulse_elevated_advection_scale),
            'wall_impulse_elevated_advection_duration_s':float(
                wall_impulse_elevated_advection_duration_s),
            'wall_impulse_elevated_transport_update_count':int(
                elevated_transport_update_count),
            'wall_impulse_elevated_return_rate_s_inv':float(
                wall_impulse_elevated_return_rate_s_inv),
            'wall_impulse_elevated_center_m':float(wall_impulse_elevated_center_m),
            'wall_impulse_elevated_sigma_m':float(wall_impulse_elevated_sigma_m),
            'wall_impulse_elevated_observation_ramp_m':float(
                wall_impulse_elevated_observation_ramp_m),
            'wall_impulse_elevated_settling_velocity_m_s':float(
                wall_impulse_elevated_settling_velocity_m_s),
            'wall_impulse_elevated_vertical_diffusivity_m2_s':float(
                wall_impulse_elevated_vertical_diffusivity_m2_s),
            'wall_impulse_elevated_min_center_m':float(
                wall_impulse_elevated_min_center_m),
            'wall_impulse_elevated_rise_velocity_m_s':float(
                wall_impulse_elevated_rise_velocity_m_s),
            'wall_impulse_elevated_rise_time_s':float(
                wall_impulse_elevated_rise_time_s),
            'wall_impulse_elevated_conserved_center_m':float(
                conserved_prompt_center),
            'wall_impulse_elevated_conserved_sigma_m':float(
                conserved_prompt_sigma),
            'wall_impulse_elevated_conserved_settling_velocity_m_s':float(
                conserved_prompt_settling),
            'wall_impulse_elevated_conserved_vertical_diffusivity_m2_s':float(
                conserved_prompt_diffusivity),
            'wall_impulse_elevated_conserved_min_center_m':float(
                conserved_prompt_min_center),
            'wall_impulse_elevated_conserved_rise_velocity_m_s':float(
                wall_impulse_elevated_conserved_rise_velocity_m_s),
            'wall_impulse_elevated_conserved_rise_time_s':float(
                wall_impulse_elevated_conserved_rise_time_s),
            'wall_impulse_elevated_conserved_transport_fraction':float(
                conserved_prompt_transport_fraction),
            'wall_impulse_elevated_conserved_crosswind_velocity_m_s':float(
                conserved_prompt_crosswind_velocity),
            'wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s':float(
                conserved_prompt_outward_velocity),
            'wall_impulse_elevated_conserved_crosswind_outward_scale_m':float(
                conserved_prompt_outward_scale),
            'wall_impulse_elevated_conserved_center_mode':(
                wall_impulse_elevated_conserved_center_mode),
            'wall_impulse_reattached_fraction':float(
                wall_impulse_reattached_fraction),
            'wall_impulse_reattached_transport_fraction':float(
                wall_impulse_reattached_transport_fraction),
            'wall_impulse_reattached_transport_height_m':(
                None if wall_impulse_reattached_transport_height_m is None
                else float(wall_impulse_reattached_transport_height_m)),
            'wall_impulse_reattached_crosswind_velocity_m_s':float(
                wall_impulse_reattached_crosswind_velocity_m_s),
            'wall_impulse_reattached_crosswind_outward_velocity_m_s':float(
                wall_impulse_reattached_crosswind_outward_velocity_m_s),
            'wall_impulse_reattached_crosswind_outward_scale_m':float(
                wall_impulse_reattached_crosswind_outward_scale_m),
            'wall_impulse_reattached_return_rate_s_inv':float(
                wall_impulse_reattached_return_rate_s_inv),
            'wall_impulse_reattached_center_m':float(
                wall_impulse_reattached_center_m),
            'wall_impulse_reattached_sigma_m':float(
                wall_impulse_reattached_sigma_m),
            'wall_impulse_reattached_settling_velocity_m_s':float(
                wall_impulse_reattached_settling_velocity_m_s),
            'wall_impulse_reattached_vertical_diffusivity_m2_s':float(
                wall_impulse_reattached_vertical_diffusivity_m2_s),
            'wall_impulse_reattached_min_center_m':float(
                wall_impulse_reattached_min_center_m),
            'wall_impulse_wake_capture_rate_s_inv':float(
                wall_impulse_wake_capture_rate_s_inv),
            'wall_impulse_wake_capture_mode':wall_impulse_wake_capture_mode,
            'wall_impulse_wake_return_rate_s_inv':float(
                wall_impulse_wake_return_rate_s_inv),
            'wall_impulse_wake_transport_fraction':float(
                wall_impulse_wake_transport_fraction),
            'wall_impulse_wake_crosswind_velocity_m_s':float(
                wall_impulse_wake_crosswind_velocity_m_s),
            'wall_impulse_wake_crosswind_outward_velocity_m_s':float(
                wall_impulse_wake_crosswind_outward_velocity_m_s),
            'wall_impulse_wake_crosswind_outward_scale_m':float(
                wall_impulse_wake_crosswind_outward_scale_m),
            'wall_impulse_wake_observation_fraction':float(
                wall_impulse_wake_observation_fraction),
            'wall_impulse_wake_observation_center_m':(
                None if wall_impulse_wake_observation_center_m is None
                else float(wall_impulse_wake_observation_center_m)),
            'wall_impulse_wake_observation_sigma_m':float(
                wall_impulse_wake_observation_sigma_m),
            'wall_impulse_wake_settling_velocity_m_s':float(
                wall_impulse_wake_settling_velocity_m_s),
            'wall_impulse_wake_vertical_diffusivity_m2_s':float(
                wall_impulse_wake_vertical_diffusivity_m2_s),
            'wall_impulse_wake_min_center_m':float(
                wall_impulse_wake_min_center_m),
            'wall_impulse_wake_species_age_moment_kg_s':float(
                wall_impulse_wake_species_age_moment.sum()*area),
            'wall_impulse_wake_mean_species_age_s':float(
                wall_impulse_wake_species_age_moment.sum()
                /max(float(wall_impulse_wake_reservoir[1].sum()),1.e-30)),
            'wall_impulse_wake_captured_total_kg':float(
                wall_impulse_wake_captured[0]),
            'wall_impulse_wake_captured_species_kg':float(
                wall_impulse_wake_captured[1]),
            'wall_impulse_wake_returned_total_kg':float(
                wall_impulse_wake_returned[0]),
            'wall_impulse_wake_returned_species_kg':float(
                wall_impulse_wake_returned[1]),
            'wall_impulse_wake_stored_total_kg':float(
                wall_impulse_wake_reservoir[0].sum()*area),
            'wall_impulse_wake_stored_species_kg':float(
                wall_impulse_wake_reservoir[1].sum()*area),
            'wall_impulse_wake_transport_moved_kg':float(
                wall_impulse_wake_transport_moved),
            'initial_elevated_species_fraction':float(
                initial_elevated_species_fraction),
            'dense_settling_velocity_m_s':float(dense_settling_velocity_m_s),
            'maximum_obstacle_vertical_velocity_m_s':float(
                maximum_obstacle_vertical_velocity),
            'maximum_wall_impulse_vertical_velocity_m_s':float(
                maximum_wall_impulse_vertical_velocity),
            'wall_impulse_lifted_species_kg':float(wall_impulse_lifted_species_kg),
            'wall_impulse_elevated_species_transferred_kg':float(
                elevated_species_transferred),
            'wall_flux_elevated_species_marked_kg':float(
                wall_flux_elevated_species_marked),
            'wall_flux_elevated_active_species_kg':float(
                wall_flux_elevated_species.sum()*area
                +local_cavity_reservoir_wall_flux_species.sum()*area
                +cavity_reservoir_wall_flux_species
                +wall_overflow_flux_marker_x.sum()
                +wall_overflow_flux_marker_y.sum()),
            'wall_flux_elevated_age_moment_kg_s':float(
                wall_flux_elevated_age_moment.sum()*area
                +local_cavity_reservoir_wall_flux_age_moment.sum()*area
                +cavity_reservoir_wall_flux_age_moment
                +wall_overflow_flux_age_x.sum()
                +wall_overflow_flux_age_y.sum()),
            'wall_flux_elevated_relative_transport_kg':float(
                wall_flux_elevated_relative_transport),
            'wall_impulse_elevated_relative_transport_kg':float(
                elevated_species_relative_transport),
            'wall_impulse_elevated_crosswind_diffusion_kg':float(
                elevated_species_crosswind_diffusion),
            'wall_impulse_elevated_species_returned_kg':float(
                elevated_species_returned),
            'wall_impulse_elevated_marker_active_species_kg':float(
                0.0 if elevated_species is None
                else elevated_species.sum()*area
                     +local_cavity_reservoir_elevated_species.sum()*area
                     +cavity_reservoir_elevated_species
                     +wall_overflow_elevated_x.sum()
                     +wall_overflow_elevated_y.sum()),
            'wall_impulse_elevated_conserved_active_species_kg':float(
                (wall_impulse_elevated_reservoir[1].sum()
                 +wall_impulse_reattached_reservoir[1].sum())*area),
            'wall_impulse_reattached_active_species_kg':float(
                wall_impulse_reattached_reservoir[1].sum()*area),
            'wall_impulse_elevated_species_active_kg':float(
                ((wall_impulse_elevated_reservoir[1].sum()
                  +wall_impulse_reattached_reservoir[1].sum())*area)
                +((0.0 if elevated_species is None
                   else elevated_species.sum()*area
                        +local_cavity_reservoir_elevated_species.sum()*area
                        +cavity_reservoir_elevated_species
                        +wall_overflow_elevated_x.sum()
                        +wall_overflow_elevated_y.sum())
                  if wall_impulse_elevated_state_mode=='hybrid' else 0.0)
                if wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                else (0.0 if elevated_species is None
                      else elevated_species.sum()*area
                           +local_cavity_reservoir_elevated_species.sum()*area
                           +cavity_reservoir_elevated_species
                           +wall_overflow_elevated_x.sum()
                           +wall_overflow_elevated_y.sum())),
            'wall_impulse_elevated_age_moment_kg_s':float(
                wall_impulse_elevated_reservoir_age_moment.sum()*area
                +wall_impulse_reattached_reservoir_age_moment.sum()*area
                +((elevated_species_age_moment.sum()*area
                   +local_cavity_reservoir_elevated_age_moment.sum()*area
                   +cavity_reservoir_elevated_age_moment
                   +wall_overflow_elevated_age_x.sum()
                   +wall_overflow_elevated_age_y.sum())
                  if wall_impulse_elevated_state_mode=='hybrid' else 0.0)
                if wall_impulse_elevated_state_mode in {'conserved','hybrid'}
                else (elevated_species_age_moment.sum()*area
                      +local_cavity_reservoir_elevated_age_moment.sum()*area
                      +cavity_reservoir_elevated_age_moment
                      +wall_overflow_elevated_age_x.sum()
                      +wall_overflow_elevated_age_y.sum())),
            'wall_impulse_elevated_conserved_mean_center_m':float(
                (wall_impulse_elevated_reservoir_center_moment.sum()
                 /max(float(wall_impulse_elevated_reservoir[1].sum()),1.e-30))
                if wall_impulse_elevated_conserved_center_mode=='impact_head'
                else conserved_prompt_center),
            'wall_impulse_elevated_captured_total_kg':float(
                wall_impulse_elevated_captured[0]),
            'wall_impulse_elevated_returned_total_kg':float(
                wall_impulse_elevated_returned[0]),
            'wall_impulse_elevated_stored_total_kg':float(
                (wall_impulse_elevated_reservoir[0].sum()
                 +wall_impulse_reattached_reservoir[0].sum())*area),
            'wall_impulse_reattached_stored_total_kg':float(
                wall_impulse_reattached_reservoir[0].sum()*area),
            'wall_impulse_reattached_captured_total_kg':float(
                wall_impulse_reattached_captured[0]),
            'wall_impulse_reattached_captured_species_kg':float(
                wall_impulse_reattached_captured[1]),
            'wall_impulse_reattached_effective_capture_fraction':float(
                wall_impulse_reattached_captured[0]
                /max(wall_impulse_elevated_captured[0],1.e-30)),
            'wall_impulse_elevated_reservoir_transport_moved_kg':float(
                wall_impulse_elevated_reservoir_transport_moved
                +wall_impulse_reattached_reservoir_transport_moved),
            'wall_impulse_reattached_reservoir_transport_moved_kg':float(
                wall_impulse_reattached_reservoir_transport_moved),
            'local_cavity_reservoir_elevated_species_kg':float(
                local_cavity_reservoir_elevated_species.sum()*area),
            'wall_impulse_elevated_subset_correction_kg':float(
                elevated_species_subset_correction),
            'multilayer_horizontal_flux_transport':bool(vertical_q is not None),
            'maximum_multilayer_horizontal_closure_residual_kg_m2':float(
                maximum_multilayer_horizontal_closure_residual),
            'vertical_domain_top_m':None if vertical_domain_top_m is None else float(vertical_domain_top_m),
            'receptor_interpolation':receptor_interpolation,
            'wake_deficit_strength':float(wake_deficit_strength),
            'wake_length_m':float(wake_length_m),
            'wake_base_width_m':float(wake_base_width_m),
            'wake_spread_slope':float(wake_spread_slope),
            'wake_recirculation_strength':float(wake_recirculation_strength),
            'wake_meander_strength':float(wake_meander_strength),
            'wake_shedding_strouhal':float(wake_shedding_strouhal),
            'wake_meander_phase_deg':float(wake_meander_phase_deg),
            'cavity_shelter_strength':float(cavity_shelter_strength),
            'cavity_entrainment_factor':float(cavity_entrainment_factor),
            'cavity_core_entrainment_reduction':float(cavity_core_entrainment_reduction),
            'cavity_vortex_rate_s_inv':float(cavity_vortex_rate_s_inv),
            'cavity_vortex_strength':float(cavity_vortex_strength),
            'cavity_reservoir_capture_rate_s_inv':float(cavity_reservoir_capture_rate_s_inv),
            'cavity_reservoir_return_rate_s_inv':float(cavity_reservoir_return_rate_s_inv),
            'cavity_reservoir_return_target':cavity_reservoir_return_target,
            'cavity_reservoir_reattached_fraction':float(
                cavity_reservoir_reattached_fraction),
            'cavity_reservoir_capture_mode':cavity_reservoir_capture_mode,
            'cavity_reservoir_vertical_capture_mode':(
                cavity_reservoir_vertical_capture_mode),
            'cavity_reservoir_spatial_mode':cavity_reservoir_spatial_mode,
            'cavity_reservoir_content_mode':cavity_reservoir_content_mode,
            'cavity_reservoir_transport_fraction':float(
                cavity_reservoir_transport_fraction),
            'cavity_reservoir_transport_height_m':(
                None if cavity_reservoir_transport_height_m is None
                else float(cavity_reservoir_transport_height_m)),
            'cavity_reservoir_transport_mask':cavity_reservoir_transport_mask,
            'cavity_reservoir_wake_length_m':float(
                cavity_reservoir_wake_length_m),
            'cavity_reservoir_wake_half_width_m':(
                None if cavity_reservoir_wake_half_width_m is None
                else float(cavity_reservoir_wake_half_width_m)),
            'cavity_reservoir_crosswind_convergence_rate_s_inv':float(
                cavity_reservoir_crosswind_convergence_rate_s_inv),
            'cavity_reservoir_crosswind_target_offset_m':float(
                cavity_reservoir_crosswind_target_offset_m),
            'local_reservoir_observation_fraction':float(
                local_reservoir_observation_fraction),
            'local_reservoir_observation_center_m':(
                None if local_reservoir_observation_center_m is None
                else float(local_reservoir_observation_center_m)),
            'local_reservoir_observation_sigma_m':float(
                local_reservoir_observation_sigma_m),
            'cavity_reservoir_return_sigma_m':float(cavity_reservoir_return_sigma_m),
            'cavity_reservoir_return_radius_m':float(cavity_reservoir_return_radius_m),
            'cavity_reservoir_return_offset_m':float(cavity_reservoir_return_offset_m),
            'cavity_reservoir_return_crosswind_scale':float(cavity_reservoir_return_crosswind_scale),
            'cavity_reservoir_face_return_fraction':float(cavity_reservoir_face_return_fraction),
            'cavity_reservoir_capture_radius_m':float(cavity_reservoir_capture_radius_m),
            'cavity_vertical_exponent':float(cavity_vertical_exponent),
            'cavity_vertical_radial_decay_power':float(
                cavity_vertical_radial_decay_power),
            'exterior_entrainment_far_factor':float(
                exterior_entrainment_far_factor),
            'exterior_entrainment_decay_length_m':float(
                exterior_entrainment_decay_length_m),
            'cavity_backflow_strength':float(cavity_backflow_strength),
            'cavity_backflow_rate_s_inv':float(cavity_backflow_rate_s_inv),
            'cavity_backflow_frame':cavity_backflow_frame,
            'cavity_stagnation_strength':float(cavity_stagnation_strength),
            'cavity_geometry':cavity_geometry,
            'cavity_align_with_wind':bool(cavity_align_with_wind),
            'cavity_reservoir_start_time_s':float(cavity_reservoir_start_time_s),
            'cavity_reservoir_return_delay_s':float(cavity_reservoir_return_delay_s),
            'cavity_return_streamwise_bias':float(cavity_return_streamwise_bias),
            'cavity_reservoir_mass_kg':float(cavity_reservoir_total()[0]),
            'cavity_reservoir_species_mass_kg':float(cavity_reservoir_total()[1]),
            'local_cavity_reservoir_mass_kg':float(
                local_cavity_reservoir[0].sum()*area),
            'local_cavity_reservoir_outside_cavity_mass_kg':float(
                local_cavity_reservoir[0,~cavity_mask].sum()*area),
            'local_cavity_reservoir_species_mass_kg':float(
                local_cavity_reservoir[1].sum()*area),
            'cavity_reservoir_pending_mass_kg':float(pending_total[0]),
            'cavity_reservoir_pending_species_mass_kg':float(pending_total[1]),
            'cavity_radius_m':float(cavity_radius_m),
            'horizontal_diffusivity_m2_s':float(horizontal_diffusivity_m2_s),
            'horizontal_diffusion_respect_walls':horizontal_diffusion_respect_walls,
            'horizontal_diffusion_blocked_species_throughput_kg':float(
                horizontal_diffusion_blocked_species_throughput_kg),
            'crosswind_diffusivity_factor':float(crosswind_diffusivity_factor),
            'minimum_conserved_mass_per_area_kg_m2':minq,
            'density_momentum_feedback':True,'vertical_ground_flux':0.,
            'gravity_scale':gravity_scale,
            'obstacle_model':'finite-height hydrostatic face reflection/overtopping; no 3-D wake turbulence',
            'entrainment_coefficient':entrainment_coefficient,
            'entrainment_far_coefficient':(
                None if entrainment_far_coefficient is None
                else float(entrainment_far_coefficient)),
            'entrainment_relaxation_time_s':float(
                entrainment_relaxation_time_s),
            'drag_coefficient':drag_coefficient,
            'dense_entrainment_model':dense_entrainment_model,
            'friction_velocity_m_s':friction_velocity_m_s,
            'wind_history_enabled':bool(history_times is not None),
            'wind_history_time_offset_s':float(wind_history_time_offset_s),
            'wind_history_sample_count':0 if history_times is None else int(history_times.size),
            'wind_history_spatial_advection_enabled':bool(spatial_wind_history),
            'wind_history_advection_reference_x_m':(
                None if wind_history_advection_reference_x_m is None
                else float(wind_history_advection_reference_x_m)),
            'wind_history_advection_speed_m_s':(
                None if wind_history_advection_speed_m_s is None
                else float(wind_history_advection_speed_m_s)),
            'wind_history_scale_friction_velocity':bool(wind_history_scale_friction_velocity),
            'wind_profile':'depth averaged neutral logarithmic' if friction_velocity_m_s is not None else 'uniform',
            'vertical_observer':('conserved multilayer state with EOS volume occupancy'
                if vertical_q is not None else
                'continuous uniform layer, not 5m cell-centre clamping'),
            'height_resolved_wall_flux':bool(height_resolved_wall_flux),
            'height_resolved_wall_steps':int(height_resolved_wall_steps),
            'exterior_reentry_fraction':float(exterior_reentry_fraction),
            'exterior_reentry_delay_s':float(exterior_reentry_delay_s),
            'exterior_reentry_sigma_m':float(exterior_reentry_sigma_m),
            'exterior_reentry_captured_total_kg':float(exterior_reentry_captured[0]),
            'exterior_reentry_captured_species_kg':float(exterior_reentry_captured[1]),
            'exterior_reentry_returned_total_kg':float(exterior_reentry_returned[0]),
            'exterior_reentry_returned_species_kg':float(exterior_reentry_returned[1]),
            'exterior_reentry_queue_remaining_kg':float(sum(p[1][0].sum() for p in exterior_reentry_queue)),
            'wall_tangential_bypass_fraction':float(wall_tangential_bypass_fraction),
            'wall_tip_bypass_fraction':float(wall_tip_bypass_fraction),
            'wall_tip_bypass_power':float(wall_tip_bypass_power),
            'maximum_wall_bypass_mass_kg':float(maximum_wall_bypass_mass_kg)},
        component_result)
