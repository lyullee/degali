import numpy as np
import pytest

from tools.audit_active_thermal_corrector_grid import (
    _conservative_remap_cell_averages,
)


def _coordinates(count, length):
    spacing = length / count
    return (np.arange(count) + 0.5) * spacing


def test_conservative_remap_preserves_constant_and_integral():
    domain = (4.8, 3.6, 2.4)
    source = (
        _coordinates(36, domain[0]),
        _coordinates(27, domain[1]),
        _coordinates(18, domain[2]),
    )
    target = (
        _coordinates(16, domain[0]),
        _coordinates(12, domain[1]),
        _coordinates(8, domain[2]),
    )
    rng = np.random.default_rng(18)
    field = rng.random((18, 27, 36))
    remapped = _conservative_remap_cell_averages(field, source, target)
    constant = _conservative_remap_cell_averages(
        np.ones_like(field), source, target,
    )
    source_volume = np.prod([
        domain[0] / 36, domain[1] / 27, domain[2] / 18,
    ])
    target_volume = np.prod([
        domain[0] / 16, domain[1] / 12, domain[2] / 8,
    ])

    assert np.all(constant == 1.0)
    assert np.sum(remapped) * target_volume == pytest.approx(
        np.sum(field) * source_volume,
    )


def test_conservative_remap_retains_an_aligned_solid_nan_region():
    domain = (1.2, 1.2, 1.2)
    source = tuple(_coordinates(12, length) for length in domain)
    target = tuple(_coordinates(6, length) for length in domain)
    field = np.ones((12, 12, 12))
    field[:4, 4:8, 4:8] = np.nan

    remapped = _conservative_remap_cell_averages(field, source, target)

    assert np.all(np.isnan(remapped[:2, 2:4, 2:4]))
    assert np.all(remapped[np.isfinite(remapped)] == 1.0)
