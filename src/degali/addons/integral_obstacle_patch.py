"""Grid-free exterior transport joined to a local resolved obstacle patch.

First operating envelope: steady wind, passive isothermal H2. The exterior is
the cross-section-integrated Gaussian advection/diffusion solution (not a
DEGADIS dense-gas closure). Only the patch is time-marched on a grid. A viscous
resolved carrier is computed once, preserving its bypass and reverse flow.
Incoming species mass is a prescribed face flux, never a concentration overwrite.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
import time
from typing import Callable, Sequence

import numpy as np
from scipy.special import ndtr

from .fast_obstacle_transport_3d import (
    FastObstacleReceptor3D, _receptor_stencil, _transport_step,
)
from .resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig, solve_resolved_obstacle_flow_3d,
    _circular_arc_blocked_faces, _projected_faces_from_cell_velocity,
)
from .semi_fv_obstacle import SourceRateSchedule


def schedule_mass(schedule: SourceRateSchedule, start: float, end: float) -> float:
    """Integrate only the intersection with a finite source event."""
    left, right = max(0.0, start), min(schedule.duration_s, end)
    return 0.0 if right <= left else schedule.mass_between(left, right)


@dataclass(frozen=True)
class GaussianIntegralInlet:
    """Exact advective delay and transverse Gaussian cell-area integration.

    Ground reflection preserves the z>=0 inventory. Tails outside the patch
    are retained as an exterior bypass; they are not renormalized into it.
    Longitudinal diffusion, density feedback and changing wind are omitted.
    """
    schedule: SourceRateSchedule
    distance_m: float
    wind_speed_m_s: float
    diffusivity_m2_s: float
    source_y_m: float
    source_z_m: float
    source_sigma_m: float

    def __post_init__(self):
        values = np.array([
            self.distance_m, self.wind_speed_m_s, self.diffusivity_m2_s,
            self.source_y_m, self.source_z_m, self.source_sigma_m,
        ])
        if not np.all(np.isfinite(values)) or self.distance_m < 0 or self.source_z_m < 0:
            raise ValueError("inlet values must be finite with nonnegative distance/height")
        if self.wind_speed_m_s <= 0 or self.diffusivity_m2_s < 0 or self.source_sigma_m <= 0:
            raise ValueError("positive wind/sigma and nonnegative diffusivity required")
        self.schedule.require_zero_endpoint("integral inlet")

    @property
    def delay_s(self):
        return self.distance_m / self.wind_speed_m_s

    def area_probabilities(self, y_edges, z_edges):
        y, z = np.asarray(y_edges, float), np.asarray(z_edges, float)
        if any(a.ndim != 1 or len(a) < 2 or np.any(~np.isfinite(a)) or np.any(np.diff(a) <= 0) for a in (y, z)):
            raise ValueError("face edges must be finite and strictly increasing")
        if z[0] < 0:
            raise ValueError("ground-reflected inlet requires nonnegative heights")
        sigma = math.sqrt(self.source_sigma_m**2 + 2*self.diffusivity_m2_s*self.delay_s)
        py = np.diff(ndtr((y-self.source_y_m)/sigma))
        pz = np.diff(ndtr((z-self.source_z_m)/sigma)) + np.diff(ndtr((z+self.source_z_m)/sigma))
        return pz[:, None]*py[None, :]

    def mass_between(self, start, end):
        return schedule_mass(self.schedule, start-self.delay_s, end-self.delay_s)


@dataclass
class IntegralObstaclePatchResult:
    time_s: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    solid_mask: np.ndarray
    final_kg_m3: np.ndarray
    maximum_kg_m3: np.ndarray
    domain_mass_kg: np.ndarray
    inlet_mass_kg: np.ndarray
    outlet_mass_kg: np.ndarray
    other_boundary_mass_kg: np.ndarray
    outlet_rate_kg_s: np.ndarray
    receptor_traces_kg_m3: dict[str, np.ndarray]
    diagnostics: dict
    boundary_rate_kg_s: dict[str, np.ndarray]

    @property
    def mass_residual_kg(self):
        return (self.domain_mass_kg+self.outlet_mass_kg+self.other_boundary_mass_kg
                -self.inlet_mass_kg-self.diagnostics.get('initial_species_mass_kg',0.))


def solve_integral_obstacle_patch(
    config: ResolvedObstacleFlow3DConfig, *,
    incoming_flux: Callable[[float, float, np.ndarray, np.ndarray], np.ndarray],
    receptors: Sequence[FastObstacleReceptor3D] = (),
    origin_m: tuple[float, float, float] = (0., 0., 0.),
    initial_species_mass_kg_m3: np.ndarray | None = None,
    additional_diffusivity_m2_s: np.ndarray | None = None,
    carrier_velocity_offset_m_s: tuple[float,float,float] = (0.,0.,0.),
) -> IntegralObstaclePatchResult:
    """Advance local scalar transport with an integral-model inlet callback.

    Callback returns interval-mean kg m^-2 s^-1 on global (z,y) faces.
    Per-step outlet face rates are exported for grid-free continuation. Ambient
    carrier momentum is resolved in the patch; scalar-density feedback is not.
    """
    started = time.perf_counter()
    if config.nonstreamwise_wall_boundary == 'all_pressure_open':
        raise ValueError('open ground carrier is incompatible with impermeable scalar ground; use ground_free_slip_lateral_top_pressure_open')
    if config.thermal_buoyancy_coupling or config.initial_binary_cloud is not None:
        raise ValueError("frozen-carrier patch currently supports passive isothermal transport only")
    if config.source_schedule.released_mass_kg != 0:
        raise ValueError("patch source must be zero; supply source through integral face flux")
    origin = np.asarray(origin_m, float)
    if origin.shape != (3,) or not np.all(np.isfinite(origin)):
        raise ValueError("origin must have three finite coordinates")
    if not callable(incoming_flux):
        raise TypeError("incoming_flux must be callable")
    labels = [r.label for r in receptors]
    if len(labels) != len(set(labels)):
        raise ValueError("receptor labels must be unique")
    # Resolve the carrier once. This call is part of the measured runtime and
    # is not an offline cache or a prepaid warm-start.
    h = config.domain_m[0]/config.cells[0]
    seed_dt = min(config.scalar_time_step_s, 1e-5)
    zero = SourceRateSchedule((0., seed_dt), (0., 0.))
    seed = solve_resolved_obstacle_flow_3d(replace(
        config, duration_s=seed_dt, scalar_time_step_s=seed_dt,
        source_schedule=zero, store_scalar_fields=False,
    ))
    carrier_wall = time.perf_counter()-started
    solid = seed.solid_mask
    u, v, w = seed.u_face_m_s, seed.v_face_m_s, seed.w_face_m_s
    blocked = _circular_arc_blocked_faces(config,seed.x_m,seed.y_m,seed.z_m)
    velocity_offset = np.asarray(carrier_velocity_offset_m_s,float)
    if velocity_offset.shape!=(3,) or np.any(~np.isfinite(velocity_offset)):
        raise ValueError('carrier velocity offset must have three finite values')
    if np.any(velocity_offset) or any(np.any(mask) for mask in blocked):
        u,v,w,*_ = _projected_faces_from_cell_velocity(
            u=seed.velocity_u_m_s+velocity_offset[0],
            v=seed.velocity_v_m_s+velocity_offset[1],
            w=seed.velocity_w_m_s+velocity_offset[2],solid=solid,
            inlet_physical=seed.u_face_m_s[:,:,0]+velocity_offset[0],
            spacing_m=h,blocked_faces=blocked,
            nonstreamwise_wall_boundary=config.nonstreamwise_wall_boundary,
        )
    diffusivity = np.full_like(seed.velocity_u_m_s, config.scalar_diffusivity_m2_s)
    if config.effective_schmidt_number is not None:
        diffusivity += seed.eddy_viscosity_m2_s/config.effective_schmidt_number
    if additional_diffusivity_m2_s is not None:
        extra = np.asarray(additional_diffusivity_m2_s,float)
        if extra.shape!=solid.shape or np.any(~np.isfinite(extra)) or np.any(extra<0):
            raise ValueError('additional diffusivity must be finite nonnegative with patch shape')
        diffusivity += extra
    # Conservative global CFL bound; do not relax it to hit a runtime target.
    rate = (np.max(np.abs(u))+np.max(np.abs(v))+np.max(np.abs(w)))/h+6*float(np.max(diffusivity))/h**2
    steps = max(1, math.ceil(config.duration_s/min(config.scalar_time_step_s, .8/max(rate, 1e-30))))
    dt = config.duration_s/steps
    times = np.linspace(0., config.duration_s, steps+1)
    y_edges = origin[1]+np.arange(config.cells[1]+1)*h
    z_edges = origin[2]+np.arange(config.cells[2]+1)*h
    x, y, z = seed.x_m+origin[0], seed.y_m+origin[1], seed.z_m+origin[2]
    stencils = [_receptor_stencil(x, y, z, solid, r) for r in receptors]
    c = np.zeros_like(solid, dtype=float)
    if initial_species_mass_kg_m3 is not None:
        c = np.asarray(initial_species_mass_kg_m3,float).copy()
        if c.shape!=solid.shape or np.any(~np.isfinite(c)) or np.any(c<0) or np.any(c[solid]!=0):
            raise ValueError('initial species field must be finite nonnegative fluid cell averages')
    initial_mass = float(c.sum()*h**3)
    maximum = c.copy()
    domain, incoming, outgoing, other = (np.zeros(steps+1) for _ in range(4))
    domain[0] = initial_mass
    outlet_rate = np.zeros((steps, config.cells[2], config.cells[1]))
    # Signed outward packets on every external face, including lateral/roof
    # exports. Retain their distribution rather than silently losing it when
    # continuing only x_max. Ground packets must remain identically zero.
    boundary_rates = {
        'x_min': np.zeros_like(outlet_rate), 'x_max': outlet_rate,
        'y_min': np.zeros((steps, config.cells[2], config.cells[0])),
        'y_max': np.zeros((steps, config.cells[2], config.cells[0])),
        'z_min': np.zeros((steps, config.cells[1], config.cells[0])),
        'z_max': np.zeros((steps, config.cells[1], config.cells[0])),
    }
    traces = {r.label: np.zeros(steps+1) for r in receptors}
    for r,stencil in zip(receptors,stencils):
        traces[r.label][0] = sum(weight*c[k,j,l] for k,j,l,weight in stencil)
    open_faces = (("x_min", "x_max", "y_min", "y_max", "z_max")
                  if config.nonstreamwise_wall_boundary == 'ground_free_slip_lateral_top_pressure_open'
                  else ("x_min", "x_max"))
    # Flux callback is invoked once for upwind/MUSCL-Hancock. Multi-stage RK
    # needs stage weights, so is intentionally not used in this executor.
    scheme = config.scalar_advection_scheme
    if scheme == "muscl_ssprk2":
        raise ValueError("stage-weighted RK flux export is not implemented")
    for i in range(steps):
        flux = np.asarray(incoming_flux(times[i], dt, y_edges, z_edges), float)
        if flux.shape != solid.shape[:2] or np.any(~np.isfinite(flux)) or np.any(flux < 0):
            raise ValueError("integral inlet must return finite nonnegative (nz,ny) flux")
        if np.any(flux[solid[:, :, 0]] != 0):
            raise ValueError("inlet flux intersects solid; enlarge patch fetch")
        inlet_step = float(flux.sum()*h*h*dt)
        captured = {}
        def capture(fx, fy, fz):
            if np.min(fx[:, :, -1]) < -1e-12:
                raise RuntimeError("reverse flow at outlet; extend the wake patch")
            outlet_rate[i] = fx[:, :, -1]*h*h
            boundary_rates['x_min'][i] = -fx[:, :, 0]*h*h
            boundary_rates['y_min'][i] = -fy[:, 0, :]*h*h
            boundary_rates['y_max'][i] = fy[:, -1, :]*h*h
            boundary_rates['z_min'][i] = -fz[0, :, :]*h*h
            boundary_rates['z_max'][i] = fz[-1, :, :]*h*h
            captured['input'] = float(fx[:, :, 0].sum()*h*h*dt)
            captured['other'] = float((fy[:, -1, :].sum()-fy[:, 0, :].sum()
                                      +fz[-1, :, :].sum()-fz[0, :, :].sum())*h*h*dt)
        c, boundary_net = _transport_step(
            c, u_face=u, v_face=v, w_face=w, solid=solid,
            diffusivity=diffusivity, dt=dt, spacing=np.full(3, h),
            advection_scheme=scheme, open_boundary_faces=open_faces,
            signed_boundary_flux=True,
            prescribed_x_min_inflow_flux=flux, face_flux_callback=capture,
            blocked_faces=blocked,
        )
        if np.min(c) < -1e-12 or np.any(~np.isfinite(c)):
            raise RuntimeError("patch concentration violates positivity/finite gate")
        # No clipping or post-hoc renormalization is used to close the ledger.
        np.maximum(maximum, c, out=maximum)
        exit_mass = float(outlet_rate[i].sum()*dt)
        if not math.isclose(captured['input'], inlet_step, rel_tol=1e-10, abs_tol=1e-15):
            raise RuntimeError("applied inlet flux differs from integral packet; enlarge upstream fetch")
        if not math.isclose(boundary_net, exit_mass+captured['other']-inlet_step,
                            rel_tol=1e-10, abs_tol=1e-15):
            raise RuntimeError("transport and independently integrated face ledgers disagree")
        incoming[i+1] = incoming[i]+inlet_step
        outgoing[i+1] = outgoing[i]+exit_mass
        other[i+1] = other[i]+captured['other']
        domain[i+1] = c.sum()*h**3
        for r, stencil in zip(receptors, stencils):
            traces[r.label][i+1] = sum(weight*c[k,j,l] for k,j,l,weight in stencil)
    scale = max(initial_mass+float(incoming[-1]), 1e-30)
    residual = domain+outgoing+other-incoming-initial_mass
    if np.max(np.abs(residual))/scale > 1e-8:
        raise RuntimeError("local patch relative species balance exceeds 1e-8")
    return IntegralObstaclePatchResult(
        times, x, y, z, solid, c, maximum, domain, incoming, outgoing,
        other, outlet_rate, traces, {
            "runtime_s": time.perf_counter()-started,
            "carrier_runtime_s": carrier_wall,
            "patch_cells": list(config.cells), "cell_size_m": h,
            "time_step_s": dt, "steps": steps,
            "exterior_time_marched_grid_cells": 0,
            "carrier_recomputed_each_scalar_step": False,
            "carrier_flow_steps": seed.diagnostics.get("flow_steps_completed"),
            "carrier_flow_residual": (float(seed.flow_residual_history[-1, 1])
                                      if seed.flow_residual_history.size else None),
            "carrier_final_statistical_residual": seed.diagnostics.get("final_flow_residual"),
            "carrier_flow_converged": bool(seed.diagnostics.get("flow_converged",False)),
            "carrier_flow_average_samples": config.flow_average_last_steps,
            "maximum_applied_scalar_diffusivity_m2_s": float(np.max(diffusivity)),
            "recirculating_lee_cell_count": seed.diagnostics.get("recirculating_lee_cell_count"),
            "maximum_absolute_divergence_s_1": seed.diagnostics.get("maximum_absolute_divergence_s_1"),
            "maximum_mass_residual_kg": float(np.max(np.abs(residual))),
            "maximum_relative_mass_residual": float(np.max(np.abs(residual))/scale),
            "initial_species_mass_kg": initial_mass,
            "embedded_surface_blocked_face_count": sum(int(np.count_nonzero(mask)) for mask in blocked),
            "operating_envelope": "steady-wind passive isothermal species on a resolved frozen carrier",
            "density_buoyancy_feedback": False,
            "wake_time_evolution": "scalar transient on fixed resolved mean flow",
            "ground_maximum_normal_velocity_m_s": float(np.max(np.abs(w[0]))),
            "ground_boundary": seed.diagnostics.get('ground_boundary'),
            "signed_outward_boundary_mass_kg": {
                name: float(rates.sum()*dt) for name,rates in boundary_rates.items()
            },
            "x_max_only_continuation_unrouted_fraction": float(
                sum(np.maximum(boundary_rates[name],0.).sum()*dt
                    for name in ('y_min','y_max','z_max'))/scale),
        },
        boundary_rates,
    )


def downstream_trace(result, *, point_m, wind_speed_m_s, diffusivity_m2_s,
                     crosswind_speed_m_s=0.):
    """Grid-free transverse Green-function continuation of outlet packets.

    Every outlet face is retained, rather than collapsed to one plume centre.
    Arrival uses the same exact characteristic delay as the inlet. Point
    sampling adds no exterior time-marched mesh.
    """
    if not np.all(np.isfinite(np.asarray((*point_m, wind_speed_m_s, diffusivity_m2_s),float))):
        raise ValueError("downstream arguments must be finite")
    h = result.diagnostics["cell_size_m"]
    outlet_x = result.x_m[-1]+h/2
    px, py, pz = point_m
    if px <= outlet_x or pz < 0 or wind_speed_m_s <= 0 or diffusivity_m2_s <= 0:
        raise ValueError("downstream point, positive wind and diffusivity required")
    delay = (px-outlet_x)/wind_speed_m_s
    py -= crosswind_speed_m_s*delay
    sigma = math.sqrt(2*diffusivity_m2_s*delay)
    # Uniform finite-cell outlet distribution, integrated analytically. No
    # cell-centre delta sources and no arbitrary numerical width are used.
    ky = (ndtr((result.y_m+h/2-py)/sigma)-ndtr((result.y_m-h/2-py)/sigma))/h
    kz = ((ndtr((result.z_m+h/2-pz)/sigma)-ndtr((result.z_m-h/2-pz)/sigma))
          +(ndtr((result.z_m+h/2+pz)/sigma)-ndtr((result.z_m-h/2+pz)/sigma)))/h
    rates = np.einsum("tzy,z,y->t", result.outlet_rate_kg_s, kz, ky)/wind_speed_m_s
    query = result.time_s-delay
    index = np.searchsorted(result.time_s, query, side="right")-1
    valid = (index >= 0)&(index < len(rates))
    trace = np.zeros_like(query)
    trace[valid] = rates[index[valid]]
    return trace


def all_boundary_trace(result, *, point_m, wind_xy_m_s,
                       horizontal_diffusivity_m2_s,
                       vertical_diffusivity_m2_s=None):
    """Continue positive species packets from every exported open face.

    This is a grid-free observation operator for the local 3-D patch.  It
    retains the packet height and face location for ``x``, ``y`` and ``z``
    exits, then applies a finite-cell Gaussian kernel along the crosswind and
    vertical directions.  It is intentionally an exterior continuation, not
    a recoupled outer CFD domain; packet mass is never multiplied or
    renormalized.
    """
    point = np.asarray(point_m, dtype=float)
    wind = np.asarray(wind_xy_m_s, dtype=float)
    if point.shape != (3,) or not np.all(np.isfinite(point)) or point[2] < 0.:
        raise ValueError("point_m must contain finite x/y and nonnegative z")
    if wind.shape != (2,) or not np.all(np.isfinite(wind)):
        raise ValueError("wind_xy_m_s must contain two finite values")
    speed = float(np.linalg.norm(wind))
    if speed <= 0.:
        raise ValueError("wind_xy_m_s must have positive magnitude")
    dh = float(horizontal_diffusivity_m2_s)
    dv = dh if vertical_diffusivity_m2_s is None else float(vertical_diffusivity_m2_s)
    if not math.isfinite(dh) or not math.isfinite(dv) or dh <= 0. or dv <= 0.:
        raise ValueError("exterior diffusivities must be positive and finite")
    h = float(result.diagnostics["cell_size_m"])
    unit = wind / speed
    t = np.asarray(result.time_s, dtype=float)
    rates_by_face = result.boundary_rate_kg_s
    x, y, z = result.x_m, result.y_m, result.z_m
    trace = np.zeros(len(t), dtype=float)

    def cdf_bin(offset, half_width, sigma):
        sigma = np.maximum(np.asarray(sigma, dtype=float), 1.e-12)
        return ndtr((half_width - offset) / sigma) - ndtr((-half_width - offset) / sigma)

    # Flatten each face into independent source cells.  This keeps the
    # staggered shapes explicit: x-faces are (z,y), y-faces (z,x), and
    # horizontal faces (y,x).
    face_sources = (
        ("x_min", np.repeat(x[0] - h / 2, len(z) * len(y)),
         np.tile(y, len(z)), np.repeat(z, len(y)),
         np.asarray(rates_by_face["x_min"], float)),
        ("x_max", np.repeat(x[-1] + h / 2, len(z) * len(y)),
         np.tile(y, len(z)), np.repeat(z, len(y)),
         np.asarray(rates_by_face["x_max"], float)),
        ("y_min", np.tile(np.repeat(x, len(z)), 1),
         np.repeat(y[0] - h / 2, len(z) * len(x)), np.repeat(z, len(x)),
         np.asarray(rates_by_face["y_min"], float)),
        ("y_max", np.tile(np.repeat(x, len(z)), 1),
         np.repeat(y[-1] + h / 2, len(z) * len(x)), np.repeat(z, len(x)),
         np.asarray(rates_by_face["y_max"], float)),
        ("z_min", np.tile(x, len(y)), np.repeat(y, len(x)),
         np.full(len(y) * len(x), z[0] - h / 2),
         np.asarray(rates_by_face["z_min"], float)),
        ("z_max", np.tile(x, len(y)), np.repeat(y, len(x)),
         np.full(len(y) * len(x), z[-1] + h / 2),
         np.asarray(rates_by_face["z_max"], float)),
    )
    for name, sx, sy, sz, raw_rates in face_sources:
        rates = np.maximum(raw_rates, 0.)
        n_steps = len(t) - 1
        expected = (n_steps, len(z), len(y)) if name.startswith("x_") else (
            (n_steps, len(z), len(x)) if name.startswith("y_") else (n_steps, len(y), len(x))
        )
        if rates.shape != expected:
            raise ValueError(f"{name} boundary rate shape is inconsistent with patch")
        cell_rates = rates.reshape(n_steps, -1)
        delta_x = point[0] - sx
        delta_y = point[1] - sy
        along = delta_x * unit[0] + delta_y * unit[1]
        valid = along > 0.
        if not np.any(valid):
            continue
        lag = along[valid] / speed
        cross = -unit[1] * delta_x[valid] + unit[0] * delta_y[valid]
        sigma_h = np.sqrt(np.maximum(2 * dh * lag, 1.e-24))
        sigma_v = np.sqrt(np.maximum(2 * dv * lag, 1.e-24))
        cross_weight = cdf_bin(cross, h / 2, sigma_h) / h
        # Reflect the vertical kernel at the impermeable ground.  Without
        # this term a packet leaving an x/y face loses the half of its
        # Gaussian tail below z=0, unlike the local patch's ground ledger.
        height_weight = (
            cdf_bin(sz[valid] - point[2], h / 2, sigma_v)
            + cdf_bin(-sz[valid] - point[2], h / 2, sigma_v)
        ) / h
        weights = cross_weight * height_weight / speed
        indices = np.flatnonzero(valid)
        delayed = np.searchsorted(t, t[:, None] - lag[None, :], side="right") - 1
        for source_index, cell_index in enumerate(indices):
            idx = delayed[:, source_index]
            ok = (idx >= 0) & (idx < len(t) - 1)
            if np.any(ok):
                trace[ok] += cell_rates[idx[ok], cell_index] * weights[source_index]
    return trace
