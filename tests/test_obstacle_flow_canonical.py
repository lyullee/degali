import numpy as np
import pytest

from degali.validation.obstacle_flow_canonical import (
    D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
    D3Q19PoiseuilleVerificationConfig,
    POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
    PoiseuilleChannelVerificationConfig,
    verify_d3q19_poiseuille_channel,
    verify_poiseuille_channel,
)


@pytest.fixture(scope="module")
def poiseuille_result():
    return verify_poiseuille_channel(
        PoiseuilleChannelVerificationConfig(
            nx=32, nz=20, steps=5_000, convergence_check_interval=500,
        )
    )


def test_d2q9_channel_matches_the_analytic_poiseuille_profile(poiseuille_result):
    result = poiseuille_result

    assert result.schema == POISEUILLE_CHANNEL_VERIFICATION_SCHEMA
    assert result.passed
    assert result.relative_l2_error < 0.003
    assert result.relative_linf_error < 0.002
    assert result.relative_mass_flux_error < 0.003
    assert result.maximum_streamwise_nonuniformity < 1.0e-12
    assert result.maximum_relative_density_variation < 1.0e-10
    assert result.runtime_s > 0.0


def test_channel_profile_is_symmetric_and_zero_at_the_hydrodynamic_walls(
    poiseuille_result,
):
    numerical = poiseuille_result.numerical_velocity_lattice

    assert np.allclose(numerical, numerical[::-1], rtol=0.0, atol=1.0e-14)
    assert numerical[0] < numerical[len(numerical) // 2]
    assert numerical[-1] < numerical[len(numerical) // 2]


def test_canonical_verifier_fails_closed_when_a_predeclared_gate_is_too_tight():
    result = verify_poiseuille_channel(
        PoiseuilleChannelVerificationConfig(
            nx=24, nz=16, steps=3_000, convergence_check_interval=300,
            relative_l2_tolerance=1.0e-5,
            relative_linf_tolerance=1.0e-5,
        )
    )

    assert not result.passed
    assert result.relative_l2_error > 1.0e-5


def test_canonical_config_rejects_unstable_or_nonphysical_inputs():
    with pytest.raises(ValueError, match="stability limit"):
        PoiseuilleChannelVerificationConfig(relaxation_time=0.5)
    with pytest.raises(ValueError, match="positive and finite"):
        PoiseuilleChannelVerificationConfig(body_acceleration_lattice=0.0)


def test_d3q19_poiseuille_matches_the_analytic_channel_profile():
    result = verify_d3q19_poiseuille_channel(
        D3Q19PoiseuilleVerificationConfig(
            nx=12,
            ny=4,
            nz=12,
            steps=6_000,
            convergence_check_interval=500,
        )
    )

    assert result.schema == D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA
    assert result.passed
    assert result.relative_l2_error < 0.01
    assert result.relative_linf_error < 0.01
    assert result.relative_mass_flux_error < 0.01
    assert result.maximum_horizontal_nonuniformity < 1.0e-12
    assert result.maximum_relative_density_variation < 1.0e-12
    assert result.diagnostics["flow_model"].startswith("periodic D3Q19")


def test_d3q19_poiseuille_rejects_the_bgk_stability_limit():
    with pytest.raises(ValueError, match="exceed 0.5"):
        D3Q19PoiseuilleVerificationConfig(relaxation_time=0.5)
