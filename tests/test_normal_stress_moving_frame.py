import numpy as np
import pytest
from degali.addons.normal_stress_moving_frame import (transform_axial_fluxes,
    recover_physical_shear, axial_row_production)


def physical_fields(point):
    s, a, b = point
    lengths = np.array([.4*np.exp(.15*s), .3*np.exp(-.04*s)])
    beta = np.array([.15, -.04])
    y, z = lengths*np.array([a, b])
    rho = np.array([1.+.1*s+.02*y+.03*z])
    u = np.array([4.+.2*s-.3*y*y-.1*z*z])
    normal = np.array([.2+.01*s+.005*(y*y+z*z)])
    physical_gradient = np.array([[-.6*y, -.2*z]])
    common = dict(density=rho, velocity=u, axial_normal_stress=normal,
        coordinates=np.array([[a, b]]), area=np.prod(lengths)/4.,
        transverse_lengths=lengths, log_width_rates=beta, curvature=0.)
    shear = np.array([[.04*y, .02*z]])
    transverse = np.array([[.1*y+.03*s, -.07*z]])
    normalized_gradient = physical_gradient*lengths
    rate = np.array([.2+np.sum(beta*np.array([a, b])*normalized_gradient[0])])
    return common, shear, transverse, normalized_gradient, rate


def fluxes(point):
    common, shear, transverse, _, _ = physical_fields(point)
    result = transform_axial_fluxes(**common, physical_shear_stress=shear, transverse_velocity=transverse)
    # Independent actual values in divergence order (s,xi_y,xi_n).
    return np.array([[result['axial_mass'][0], *result['relative_mass'][0]],
                     [result['axial_momentum'][0], *result['relative_momentum'][0]],
                     [result['axial_kinetic_and_normal_work'][0], *result['relative_kinetic_and_normal_work'][0]]])


def test_inverse_shear_retains_moving_normal_stress_correction():
    common, shear, transverse, _, _ = physical_fields([.7, .8, .6])
    transformed = transform_axial_fluxes(**common, physical_shear_stress=shear, transverse_velocity=transverse)
    out = recover_physical_shear(**common, relative_mass=transformed['relative_mass'],
                                 relative_momentum=transformed['relative_momentum'])
    np.testing.assert_allclose(out['physical_shear_stress'], shear, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(out['shear_covariance'], shear/common['density'][:, None], rtol=1e-13)
    naive = common['transverse_lengths']/common['area']*out['relative_stress']
    assert max(abs((naive-shear).ravel())) > 1e-3
    assert not out['physical_closure_passed']


@pytest.mark.parametrize('zero', ['normal', 'width_rate'])
def test_old_shear_expression_is_recovered_only_in_declared_limits(zero):
    common, shear, transverse, _, _ = physical_fields([.7, .8, .6])
    common['axial_normal_stress' if zero == 'normal' else 'log_width_rates'] *= 0.
    transformed = transform_axial_fluxes(**common, physical_shear_stress=shear, transverse_velocity=transverse)
    out = recover_physical_shear(**common, relative_mass=transformed['relative_mass'],
                                 relative_momentum=transformed['relative_momentum'])
    np.testing.assert_array_equal(out['geometric_correction'], np.zeros((1, 2)))
    np.testing.assert_allclose(common['transverse_lengths']/common['area']*out['relative_stress'], shear, rtol=1e-12)


@pytest.mark.parametrize('step', [1e-4, 5e-5])
def test_actual_flux_divergences_obey_product_identity_even_off_solution(step):
    point = np.array([.7, .8, .6])
    common, shear, transverse, gradient, rate = physical_fields(point)
    transformed = transform_axial_fluxes(**common, physical_shear_stress=shear, transverse_velocity=transverse)
    result = axial_row_production(**common, relative_mass=transformed['relative_mass'],
        relative_momentum=transformed['relative_momentum'], normalized_velocity_gradient=gradient,
        axial_velocity_rate_at_fixed_coordinates=rate)
    divergence = sum((fluxes(point+step*direction)-fluxes(point-step*direction))[:, i]/(2*step)
                     for i, direction in enumerate(np.eye(3)))
    u = common['velocity'][0]
    left = divergence[2]+result['axial_row_production'][0]
    right = u*divergence[1]-.5*u*u*divergence[0]
    assert abs(left-right)/max(1., abs(left), abs(right)) < 1e-7
    assert min(abs(divergence[:2])) > 1e-3  # Do not quietly assume solved continuity/momentum.
    np.testing.assert_allclose(result['axial_row_production'], result['physical_coordinate_production'], rtol=1e-12, atol=1e-12)
    assert result['physical_axial_velocity_gradient'][0] == pytest.approx(.2, abs=1e-14)
    assert result['axial_row_production'][0] < 0.  # No sign clipping.
    assert not result['full_tensor_production'] and not result['pressure_work_included']


def test_curvature_and_invalid_geometry_are_not_silently_omitted():
    common, shear, transverse, _, _ = physical_fields([.7, .8, .6])
    for name, value, match in (('curvature', .001, 'zero curvature'), ('area', 0., 'area'),
                               ('axial_normal_stress', np.array([-.1]), 'nonnegative'),
                               ('density', np.array([0.]), 'density'),
                               ('transverse_lengths', np.array([.2]), 'lengths')):
        with pytest.raises(ValueError, match=match):
            transform_axial_fluxes(**dict(common, **{name: value}), physical_shear_stress=shear,
                                  transverse_velocity=transverse)


def test_finite_inputs_that_overflow_fluxes_are_rejected():
    common, shear, transverse, _, _ = physical_fields([.7, .8, .6])
    common['velocity'] = np.array([1e200])
    with np.errstate(over='ignore', invalid='ignore'), pytest.raises(ValueError, match='nonfinite transformed'):
        transform_axial_fluxes(**common, physical_shear_stress=shear, transverse_velocity=transverse)
