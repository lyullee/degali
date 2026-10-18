import numpy as np
import pytest

from degali.addons.conservative_column_observation import reconstruct_column


@pytest.mark.parametrize('power', [0., 2., 12.])
@pytest.mark.parametrize('volume', [.05, .5, 2.])
def test_column_species_volume_is_conserved_even_when_core_saturates(power, volume):
    r = reconstruct_column(height_m=[2.], core_species_volume_m=[volume], exponent=power)
    assert np.sum(r.mole_fraction*np.diff(r.edges_m)) == pytest.approx(volume, abs=1e-12)
    assert r.mole_fraction.min() >= 0.
    assert r.mole_fraction.max() <= 1.+1e-12
    # Independently sample the entire profile: do not just trust a ledger.
    z = (np.arange(64000)+.5)*2./64000
    cell = np.searchsorted(r.edges_m[0], z, side='right')-1
    assert np.mean(r.mole_fraction[0, cell])*2. == pytest.approx(volume, abs=2e-10)


@pytest.mark.parametrize('tail', ['exponential', 'diffusive'])
def test_tail_and_overlapping_upper_branches_share_one_bounded_inventory(tail):
    r = reconstruct_column(height_m=[2., .4], core_species_volume_m=[.5, .1],
        exponent=[12., 0.], tail_kind=tail, tail_extent_m=3., tail_scale_m=.5,
        gaussian_branches=[([.2, .05], [.3, .2], [.08, .1]),
                           ([.1, .04], [2.4, .1], [.5, .1])])
    np.testing.assert_allclose(np.sum(r.mole_fraction*np.diff(r.edges_m), axis=1), [.8, .19], atol=1e-12)
    assert r.mole_fraction.max() <= 1.+1e-12
    assert r.maximum_relative_volume_residual < 1e-12


def test_pure_source_reconstructs_one_hundred_percent_not_lost_clipped_mass():
    r = reconstruct_column(height_m=[13.], core_species_volume_m=[13.], exponent=12.)
    np.testing.assert_allclose(r.mole_fraction, 1.)
    assert r.saturated_column_count == 1


def test_empty_column_and_outside_support_are_zero():
    r = reconstruct_column(height_m=[0., 2.], core_species_volume_m=[0., .5])
    np.testing.assert_array_equal(r.sample([0., 3.]), [0., 0.])


def test_impossible_inventory_is_rejected_not_clipped_or_silently_expanded():
    with pytest.raises(ValueError, match='support capacity'):
        reconstruct_column(height_m=[1.], core_species_volume_m=[2.])


def test_bin_refinement_and_receptor_height_do_not_change_total_or_geometry():
    for bins in (32, 64, 128):
        r = reconstruct_column(height_m=[2.], core_species_volume_m=[.5], exponent=2., bins=bins)
        exact = .5*3./2.*(1.-.4/2.)**2
        assert abs(r.sample([.4])[0]-exact) < .04/bins*32
        assert np.sum(r.mole_fraction*np.diff(r.edges_m)) == pytest.approx(.5)


def test_linear_samples_preserve_integral_and_bounds_with_overlapping_kernels():
    r = reconstruct_column(height_m=[2.], core_species_volume_m=[.5], exponent=12.,
        gaussian_branches=[([.2], [2.4], [.5]), ([.1], [2.4+1e-15], [.5])])
    z = (np.arange(100000)+.5)*7./100000
    sampled = r.sample(z, np.zeros(len(z), int))
    assert sampled.min() >= 0.
    assert sampled.max() <= 1.+1e-12
    assert sampled.mean()*7. == pytest.approx(.8, abs=2e-5)
    near = r.sample([2.4-1e-12, 2.4, 2.4+1e-12], [0,0,0])
    assert np.ptp(near) < .002


def test_each_linear_bin_retains_its_exact_mean_and_capacity():
    r = reconstruct_column(height_m=[2.], core_species_volume_m=[.9], exponent=12.)
    width = np.diff(r.edges_m)
    low = r.mole_fraction-.5*width*r.slopes_per_m
    high = r.mole_fraction+.5*width*r.slopes_per_m
    assert min(low.min(), high.min()) >= -1e-12
    assert max(low.max(), high.max()) <= 1.+1e-12
    assert np.sum(.5*(low+high)*width) == pytest.approx(.9, abs=1e-12)
