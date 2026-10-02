import numpy as np
import pytest

from degali.validation.scalar_fluctuations import (
    conditional_active_scalar, intermittency_from_relative_rms,
)


def test_hecht_intermitttency_relation_is_preserved_without_clipping():
    out = intermittency_from_relative_rms([.22, .23], .049)
    np.testing.assert_allclose(out['occupancy'], (1.049)/(1.+np.array([.22, .23])**2))
    # Rounded public inputs straddle gamma=1; retain that information.
    np.testing.assert_array_equal(out['occupancy_physical'], [False, True])
    assert not out['clipped'] and not out['adopted']


def test_incompatible_occupancy_is_flagged_not_repaired():
    out = intermittency_from_relative_rms(.1, .049)
    assert out['occupancy'] > 1.
    assert not out['occupancy_physical']


def test_conditional_active_identity_and_input_rejection():
    out = conditional_active_scalar([2., 3.], [.5, 1.])
    np.testing.assert_allclose(out['conditional_active_mean'], [4., 3.])
    with pytest.raises(ValueError, match='occupancy'):
        conditional_active_scalar(1., 1.1)
