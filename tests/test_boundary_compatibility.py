import numpy as np
import pytest
from degali.addons.boundary_compatibility import boundary_compatibility


def sample(**changes):
    args=dict(y=.2,specific_h=-5000.,velocity=15.,ambient_speed=2.,theta=.4,
        diffusion_y=-.02,diffusion_h=480.,scales=[.02,500.,3.])
    args.update(changes)
    return boundary_compatibility(**args)


def test_analytic_minimax_is_a_sharp_bound_not_a_fitted_flux_law():
    out=sample()
    assert out['scaled_left_null']@out['scaled_matrix']==pytest.approx([0.,0.],abs=1e-10)
    assert max(abs(out['unconstrained_residual']))==pytest.approx(out['unavoidable_minimax'],abs=1e-12)
    assert out['physical_minimax']+1e-12>=out['unavoidable_minimax']
    assert max(abs(out['physical_residual']))==pytest.approx(out['physical_minimax'],abs=1e-10)
    assert not out['sufficient_for_conservative_flow']


def test_compatible_gradient_gives_zero_for_all_three_boundaries():
    out=sample()
    perfect=sample(diffusion_h=out['effective_specific_h']/.2*(-.02))
    assert perfect['physical_minimax']<1e-12
    mass,momentum=perfect['physical_optimal_flux']
    assert mass<0. and momentum-15.*mass>0.


def test_minimax_scales_with_uniform_diagnostic_denominator():
    a,b=sample(),sample(scales=[.04,1000.,6.])
    assert b['unavoidable_minimax']==pytest.approx(a['unavoidable_minimax']/2.)
    assert b['physical_minimax']==pytest.approx(a['physical_minimax']/2.)


@pytest.mark.parametrize('change',[dict(y=0.),dict(scales=[1.,0.,1.]),dict(diffusion_h=np.nan)])
def test_invalid_compatibility_inputs_rejected(change):
    with pytest.raises(ValueError):
        sample(**change)
