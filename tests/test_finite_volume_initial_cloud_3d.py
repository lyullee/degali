import numpy as np
import pytest

from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
    initialise_regular_polygon_prism_cloud,
)


def _r12_closure() -> BinaryIdealGasMixtureClosure:
    return BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )


def _centres(length: float, spacing: float) -> np.ndarray:
    count = int(round(length / spacing))
    return (np.arange(count) + 0.5) * spacing


@pytest.mark.parametrize("spacing", [1.0, 0.5, 0.25])
def test_ti021_dodecagonal_prism_has_grid_independent_exact_volume_and_mass(spacing):
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(20.0, 20.0),
        flat_to_flat_m=14.0,
        z_min_m=0.0,
        z_max_m=13.0,
        species_mole_fraction=0.32,
        temperature_k=291.15,
        sides=12,
        label="TI021 initial gas holder",
    )
    initial = initialise_regular_polygon_prism_cloud(
        _centres(40.0, spacing),
        _centres(40.0, spacing),
        _centres(16.0, spacing),
        cloud=cloud,
        closure=_r12_closure(),
    )

    assert cloud.exact_volume_m3 == pytest.approx(2048.20, rel=1.0e-4)
    assert initial.represented_geometry_volume_m3 == pytest.approx(
        cloud.exact_volume_m3, rel=2.0e-13,
    )
    assert initial.integrated_species_mass_kg == pytest.approx(
        initial.exact_species_mass_kg, rel=2.0e-13,
    )
    assert np.min(initial.volume_fraction) >= 0.0
    assert np.max(initial.volume_fraction) <= 1.0


def test_initial_cloud_fails_closed_when_it_overlaps_a_solid_cell():
    spacing = 1.0
    x = _centres(20.0, spacing)
    y = _centres(20.0, spacing)
    z = _centres(16.0, spacing)
    solid = np.zeros((len(z), len(y), len(x)), dtype=bool)
    solid[0, 10, 10] = True
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(10.0, 10.0), flat_to_flat_m=6.0,
        z_min_m=0.0, z_max_m=4.0,
        species_mole_fraction=0.32, temperature_k=291.15,
    )

    with pytest.raises(ValueError, match="intersects a solid"):
        initialise_regular_polygon_prism_cloud(
            x, y, z, cloud=cloud, closure=_r12_closure(), solid_mask=solid,
        )


def test_partial_cells_preserve_species_mass_without_boolean_mask_aliasing():
    spacing = 1.0
    cloud = RegularPolygonPrismInitialCloud(
        center_xy_m=(10.25, 9.75), flat_to_flat_m=5.4,
        z_min_m=0.2, z_max_m=3.7,
        species_mole_fraction=0.32, temperature_k=291.15,
        rotation_deg=17.0,
    )
    initial = initialise_regular_polygon_prism_cloud(
        _centres(20.0, spacing), _centres(20.0, spacing),
        _centres(8.0, spacing), cloud=cloud, closure=_r12_closure(),
    )

    assert np.any((initial.volume_fraction > 0.0) & (initial.volume_fraction < 1.0))
    assert initial.relative_geometry_volume_error == pytest.approx(0.0, abs=2e-14)
    assert initial.relative_species_mass_error == pytest.approx(0.0, abs=2e-14)
