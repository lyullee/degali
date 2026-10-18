import numpy as np
import pytest

from degali.addons.dense_gas_stratification_3d import (
    diagnose_gradient_richardson_3d,
    finite_volume_vertical_gradient,
    summarize_gradient_richardson_regions,
    stratification_damped_diffusivity_tensor_3d,
)


def test_finite_volume_vertical_gradient_preserves_linear_interior():
    z = np.arange(5.0)[:, None, None]
    field = np.broadcast_to(2.0 * z + 3.0, (5, 2, 3)).copy()
    solid = np.zeros_like(field, dtype=bool)
    gradient = finite_volume_vertical_gradient(
        field, solid, grid_spacing_m=1.0,
    )
    assert gradient[1:-1] == pytest.approx(2.0)
    assert gradient[0] == pytest.approx(1.0)
    assert gradient[-1] == pytest.approx(1.0)


def test_gradient_richardson_detects_stable_dense_layer_without_flux_change():
    density = np.broadcast_to(
        np.array([1.4, 1.3, 1.2, 1.1])[:, None, None], (4, 2, 2),
    ).copy()
    u = np.zeros_like(density)
    v = np.zeros_like(density)
    solid = np.zeros_like(density, dtype=bool)
    result = diagnose_gradient_richardson_3d(
        density, u, v, solid,
        grid_spacing_m=1.0,
        shear_regularization_s_inv=0.1,
    )
    expected = 9.80665 * 0.1 / 1.3 / 0.1**2
    assert result.gradient_richardson_number[1] == pytest.approx(expected)
    record = result.record()
    assert record["fluxes_modified"] is False
    assert record["observation_derived_coefficients"] == []


def test_gradient_richardson_clips_unstable_density_gradient_and_solids():
    density = np.broadcast_to(
        np.array([1.0, 1.1, 1.2])[:, None, None], (3, 1, 2),
    ).copy()
    u = np.zeros_like(density)
    v = np.zeros_like(density)
    solid = np.zeros_like(density, dtype=bool)
    solid[1, 0, 1] = True
    result = diagnose_gradient_richardson_3d(
        density, u, v, solid,
        grid_spacing_m=1.0,
        shear_regularization_s_inv=0.1,
    )
    assert np.all(result.gradient_richardson_number == 0.0)
    assert result.gradient_richardson_number[1, 0, 1] == 0.0


def test_gradient_richardson_accepts_nan_obstacle_cells_and_excludes_them():
    density = np.broadcast_to(
        np.array([1.4, 1.3, 1.2])[:, None, None], (3, 2, 2),
    ).copy()
    u = np.zeros_like(density)
    v = np.zeros_like(density)
    solid = np.zeros_like(density, dtype=bool)
    solid[1, 0, 1] = True
    density[solid] = np.nan
    result = diagnose_gradient_richardson_3d(
        density, u, v, solid,
        grid_spacing_m=1.0,
        shear_regularization_s_inv=0.1,
    )
    assert np.all(np.isfinite(result.gradient_richardson_number))
    assert result.gradient_richardson_number[1, 0, 1] == 0.0


def test_region_summary_reports_undamped_neutral_exposure_separately():
    density = np.broadcast_to(
        np.array([1.4, 1.3, 1.2])[:, None, None], (3, 2, 2),
    ).copy()
    zeros = np.zeros_like(density)
    solid = np.zeros_like(density, dtype=bool)
    result = diagnose_gradient_richardson_3d(
        density, zeros, zeros, solid,
        grid_spacing_m=2.0,
        shear_regularization_s_inv=0.2,
    )
    summary = summarize_gradient_richardson_regions(
        result,
        {"wake": np.ones_like(solid, dtype=bool)},
        neutral_scalar_eddy_diffusivity_m2_s=np.full_like(density, 0.3),
        cell_volume_m3=8.0,
    )
    wake = summary["regions"]["wake"]
    assert wake["fluid_cell_count"] == 12
    assert wake["sampling_volume_m3"] == pytest.approx(96.0)
    assert wake["neutral_scalar_eddy_diffusivity_mean_m2_s"] == pytest.approx(0.3)
    assert wake["stable_neutral_diffusivity_volume_integral_m5_s"] == (
        pytest.approx(0.3 * 12 * 8.0)
    )
    assert summary["damping_applied"] is False


def test_stratification_candidate_damps_only_vertical_turbulent_part():
    fluid = np.ones((2, 1, 3), dtype=bool)
    fluid[:, :, 2] = False
    ri = np.broadcast_to(np.array([0.0, 1.0, 99.0]), fluid.shape).copy()
    turbulent = np.full(fluid.shape, 0.4)
    molecular = 0.02
    before = turbulent.copy()

    dx, dy, dz = stratification_damped_diffusivity_tensor_3d(
        molecular, turbulent, ri, fluid, alpha=5.0, exponent=1.0,
    )

    assert dx[0, 0, 0] == pytest.approx(molecular + 0.4)
    assert dy[0, 0, 1] == pytest.approx(molecular + 0.4)
    assert dz[0, 0, 0] == pytest.approx(molecular + 0.4)
    assert dz[0, 0, 1] == pytest.approx(molecular + 0.4 / 6.0)
    assert dz[0, 0, 1] > molecular  # no artificial critical-Ri cutoff
    assert np.all(dx[~fluid] == 0.0)
    assert np.all(dy[~fluid] == 0.0)
    assert np.all(dz[~fluid] == 0.0)
    assert turbulent == pytest.approx(before)
    assert not np.shares_memory(dx, dy)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"alpha": 0.0, "exponent": 1.0},
        {"alpha": 1.0, "exponent": -1.0},
    ],
)
def test_stratification_candidate_rejects_invalid_parameters(kwargs):
    fluid = np.ones((2, 2, 2), dtype=bool)
    with pytest.raises(ValueError):
        stratification_damped_diffusivity_tensor_3d(
            0.01, np.full(fluid.shape, 0.2), np.zeros(fluid.shape), fluid,
            **kwargs,
        )
