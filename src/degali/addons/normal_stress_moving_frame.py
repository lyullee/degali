"""Straight expanding-frame axial Reynolds stress/energy identities only.

The returned shear covariance includes the moving-coordinate normal-stress
correction. No curvature, pressure, transverse-stress or turbulence closure
is supplied. Negative production is retained, never clipped.
"""

import math
import numpy as np


def _checked(result):
    if any(isinstance(value, np.ndarray) and not np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('nonfinite transformed stress, flux or production')
    return result


def _vector(value, count, name):
    result = np.asarray(value, float)
    if result.shape != (count, 2) or not np.all(np.isfinite(result)):
        raise ValueError(f'{name} must be a finite (points,2) array')
    return result


def _scalar(value, count, name):
    result = np.asarray(value, float)
    if result.shape != (count,) or not np.all(np.isfinite(result)):
        raise ValueError(f'{name} must be a finite (points,) array')
    return result


def _inputs(density, velocity, axial_normal_stress, coordinates, area,
            transverse_lengths, log_width_rates, curvature):
    density = np.asarray(density, float)
    if density.ndim != 1 or not len(density) or not np.all(np.isfinite(density)) or np.any(density <= 0.):
        raise ValueError('positive finite density samples required')
    count = len(density)
    velocity = _scalar(velocity, count, 'velocity')
    normal = _scalar(axial_normal_stress, count, 'axial normal stress')
    coords = _vector(coordinates, count, 'coordinates')
    if np.any(normal < 0.):
        raise ValueError('nonnegative axial covariance stress required')
    area = float(area)
    lengths, beta = np.asarray(transverse_lengths, float), np.asarray(log_width_rates, float)
    if (not math.isfinite(area) or area <= 0. or lengths.shape != (2,) or beta.shape != (2,)
            or not np.all(np.isfinite(lengths)) or not np.all(np.isfinite(beta)) or np.any(lengths <= 0.)):
        raise ValueError('positive finite area/lengths and two finite log-width rates required')
    if not np.isscalar(curvature) or not np.isfinite(curvature) or curvature != 0.:
        raise ValueError('explicit zero curvature required by this straight-axis identity')
    return density, velocity, normal, coords, area, lengths, beta


def transform_axial_fluxes(*, density, velocity, transverse_velocity,
                          axial_normal_stress, physical_shear_stress,
                          coordinates, area, transverse_lengths, log_width_rates, curvature):
    rho, u, normal, xi, area, lengths, beta = _inputs(density, velocity, axial_normal_stress,
        coordinates, area, transverse_lengths, log_width_rates, curvature)
    transverse = _vector(transverse_velocity, len(rho), 'transverse velocity')
    shear = _vector(physical_shear_stress, len(rho), 'physical shear stress')
    mass = area*rho[:, None]*(transverse/lengths-u[:, None]*beta*xi)
    momentum = u[:, None]*mass+area*shear/lengths-area*normal[:, None]*beta*xi
    return _checked(dict(relative_mass=mass, relative_momentum=momentum,
        relative_kinetic_and_normal_work=u[:, None]*momentum-.5*u[:, None]**2*mass,
        axial_mass=area*rho*u, axial_momentum=area*(rho*u*u+normal),
        axial_kinetic_and_normal_work=area*(.5*rho*u**3+u*normal),
        straight_axis_only=True, physical_closure_passed=False, adopted=False))


def recover_physical_shear(*, relative_mass, relative_momentum, density, velocity,
                           axial_normal_stress, coordinates, area, transverse_lengths,
                           log_width_rates, curvature):
    rho, u, normal, xi, area, lengths, beta = _inputs(density, velocity, axial_normal_stress,
        coordinates, area, transverse_lengths, log_width_rates, curvature)
    mass = _vector(relative_mass, len(rho), 'relative mass flux')
    momentum = _vector(relative_momentum, len(rho), 'relative momentum flux')
    relative = momentum-u[:, None]*mass
    correction = area*normal[:, None]*beta*xi
    shear = lengths/area*(relative+correction)
    return _checked(dict(physical_shear_stress=shear, shear_covariance=shear/rho[:, None],
        relative_stress=relative, geometric_correction=correction,
        straight_axis_only=True, physical_closure_passed=False, adopted=False))


def axial_row_production(*, relative_mass, relative_momentum, density, velocity,
                         axial_normal_stress, coordinates, area, transverse_lengths,
                         log_width_rates, curvature, normalized_velocity_gradient,
                         axial_velocity_rate_at_fixed_coordinates):
    out = recover_physical_shear(relative_mass=relative_mass, relative_momentum=relative_momentum,
        density=density, velocity=velocity, axial_normal_stress=axial_normal_stress,
        coordinates=coordinates, area=area, transverse_lengths=transverse_lengths,
        log_width_rates=log_width_rates, curvature=curvature)
    count = len(out['relative_stress'])
    gradient = _vector(normalized_velocity_gradient, count, 'normalized velocity gradient')
    rate = _scalar(axial_velocity_rate_at_fixed_coordinates, count, 'fixed-coordinate velocity rate')
    xi, beta, lengths = np.asarray(coordinates), np.asarray(log_width_rates), np.asarray(transverse_lengths)
    physical_rate = rate-np.sum(beta*xi*gradient, axis=1)
    physical_gradient = gradient/lengths
    normal = np.asarray(axial_normal_stress)
    moving = -np.sum(out['relative_stress']*gradient, axis=1)-area*normal*rate
    physical = -area*(normal*physical_rate+np.sum(out['physical_shear_stress']*physical_gradient, axis=1))
    return _checked(dict(**out, axial_row_production=moving, physical_coordinate_production=physical,
        physical_axial_velocity_gradient=physical_rate, physical_transverse_velocity_gradient=physical_gradient,
        full_tensor_production=False, pressure_work_included=False, negative_production_clipped=False))
