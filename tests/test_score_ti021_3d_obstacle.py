import pytest

from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from tools.score_ti021_3d_obstacle import vol_percent_per_kg_m3


def test_r12_species_mass_concentration_roundtrips_known_32_vol_percent():
    closure = BinaryIdealGasMixtureClosure('R12', .120913, 606.,
        ambient_temperature_k=291.15, ambient_pressure_pa=100500.)
    state = closure.state_from_mole_fraction_and_temperature(.32, 291.15)
    observed = float(state.species_mass_concentration_kg_m3)*vol_percent_per_kg_m3(closure)
    assert observed == pytest.approx(32., rel=1e-10)
    # A mass-mixing-ratio conversion cannot stand in for this observation.
    legacy = float(state.species_mass_concentration_kg_m3)*100/closure.ambient_density_kg_m3
    assert legacy > 100.
