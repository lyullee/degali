from degali.addons.cryogenic_air import relax_two_velocity_drag
from degali.addons.tke_parameterization import dissipation_time_from_integral_scale
from degali.addons.transport_evidence import (
    PARTICLE_SLIP_BOUNDARY,
    TURBULENCE_BOUNDARY,
)


def test_turbulence_boundary_is_not_a_default_or_quantitative_lh2_closure():
    assert TURBULENCE_BOUNDARY.numerical_boundary_verified
    assert not TURBULENCE_BOUNDARY.physical_closure_validated
    assert not TURBULENCE_BOUNDARY.default_prediction_enabled
    assert not TURBULENCE_BOUNDARY.quantitative_lh2_prediction_allowed
    out = dissipation_time_from_integral_scale(.1, 1.0, .2)
    assert out["transport_evidence_boundary"] == TURBULENCE_BOUNDARY.identifier
    assert not out["quantitative_lh2_prediction_allowed"]


def test_particle_slip_kernel_carries_the_same_nonpromotion_boundary():
    assert PARTICLE_SLIP_BOUNDARY.numerical_boundary_verified
    assert not PARTICLE_SLIP_BOUNDARY.physical_closure_validated
    result = relax_two_velocity_drag(
        gas_mass_flow=1.0, particle_mass_flow=.1, gas_velocity=20.0,
        particle_velocity=0.0, duration=.01, particle_relaxation=.02,
    )
    assert result.transport_evidence_boundary == PARTICLE_SLIP_BOUNDARY.identifier
    assert not result.physical_closure_validated
    assert not result.quantitative_lh2_prediction_allowed
