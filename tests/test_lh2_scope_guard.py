import pytest
import numpy as np

from degali.lh2 import (
    ApplicabilityError,
    Assessment,
    AssessmentEnvelope,
)


def _assessment(*warnings):
    return Assessment(
        distance_to_lfl=float("nan"),
        distance_to_stoichiometric=float("nan"),
        lowest_flammable_height=0.0,
        regime="grounded",
        neutral_buoyancy=0.999,
        warnings=list(warnings),
    )


def test_scope_guard_distinguishes_qualified_conditional_and_steered():
    assert _assessment().screening_scope == "qualified"
    conditional = _assessment("pool diameter 2 is outside the pool range")
    assert conditional.screening_scope == "conditional"
    conditional.require_screening_scope(allow_conditional=True)
    with pytest.raises(ApplicabilityError):
        conditional.require_screening_scope()

    steered = _assessment(
        "exit velocity is only 3.0 times the wind; below 10 the wind steers "
        "the plume and the concentrations are not defensible"
    )
    assert steered.screening_scope == "out_of_scope"
    with pytest.raises(ApplicabilityError):
        steered.require_screening_scope(allow_conditional=True)


def test_envelope_reports_input_sensitivity_without_statistical_claims():
    envelope = AssessmentEnvelope(
        scenarios=(_assessment(), _assessment("distance outside range")),
        rates_kg_s=(1.0, 2.0),
        winds_m_s=(1.0, 2.0),
    )
    assert envelope.screening_scope == "conditional"
    assert envelope.warnings == ("distance outside range",)
    assert envelope.distance_to_lfl_range[0] != envelope.distance_to_lfl_range[0]


def test_assess_envelope_evaluates_every_rate_wind_corner(monkeypatch):
    import degali.lh2 as lh2

    calls = []

    def fake_assess(**kwargs):
        calls.append((kwargs["rate"], kwargs["wind"]))
        return _assessment()

    monkeypatch.setattr(lh2, "assess", fake_assess)
    result = lh2.assess_envelope(
        rates=[1.0, 2.0], winds=[3.0, 4.0], pool_diameter=1.0
    )
    assert calls == [(1.0, 3.0), (1.0, 4.0), (2.0, 3.0), (2.0, 4.0)]
    assert len(result.scenarios) == 4


def test_observation_envelope_keeps_all_admissible_source_wind_hypotheses(monkeypatch):
    import degali.lh2 as lh2

    calls = []

    def fake_assess(**kwargs):
        calls.append((kwargs["rate"], kwargs["wind"]))
        concentration = kwargs["rate"] / kwargs["wind"]
        return Assessment(
            distance_to_lfl=float("nan"), distance_to_stoichiometric=float("nan"),
            lowest_flammable_height=0.0, regime="grounded", neutral_buoyancy=0.999,
            trajectory=np.array([[0.0, 0.0, concentration], [10.0, 0.0, concentration]]),
        )

    monkeypatch.setattr(lh2, "assess", fake_assess)
    envelope = lh2.assess_observation_envelope(
        rates=[1.0, 2.0], winds=[1.0, 2.0], observed_mole_fraction=1.0,
        distance_m=10.0, pool_diameter=1.0, acceptance_factor=1.5,
    )
    assert calls == [(1.0, 1.0), (1.0, 2.0), (2.0, 1.0), (2.0, 2.0)]
    assert len(envelope.rows) == 4
    assert len(envelope.admissible_rows) == 2
    assert not envelope.identifiable
    assert envelope.source_rate_range == (1.0, 2.0)
    assert envelope.wind_range == (1.0, 2.0)


def test_observation_envelope_does_not_silently_use_sensor_height_operator(monkeypatch):
    import degali.lh2 as lh2

    with pytest.raises(ValueError, match="observation adapter"):
        lh2.assess_observation_envelope(
            rates=[1.0], winds=[1.0], observed_mole_fraction=0.1,
            distance_m=1.0, pool_diameter=1.0,
            observation_operator="arc_max_1s_peak",
        )


def test_sensor_projection_uses_exact_vertical_and_lateral_operator(monkeypatch):
    import degali.validation.nearfield as nearfield
    import degali.lh2 as lh2

    class FakeSource:
        velocity = 100.0

    class FakeJet:
        axisymmetric_source = FakeSource()
        th = type("Thermo", (), {"table": object()})()

        def run(self, *args, **kwargs):
            return type("Run", (), {"rows": [[0.0], [10.0], [20.0]]})()

    class FakeTrajectory:
        def __init__(self, table, rows):
            pass

        @property
        def ok(self):
            return True

        def concentration_at(self, x, y, z):
            return 100.0 * (1.0 + y + z) / 10.0

        def at(self, x):
            return object()

        def temperature_at(self, x, y, z):
            return 80.0 + z

    monkeypatch.setattr(nearfield, "hydrogen_jet", lambda **kwargs: (FakeJet(), object()))
    monkeypatch.setattr(nearfield, "Trajectory", FakeTrajectory)
    monkeypatch.setattr(nearfield, "STEP", 0.1)
    projection = lh2.project_lh2_jet_to_sensors(
        rate=0.2, wind=2.0, height=0.5, orifice=0.01,
        points_m=[(5.0, 0.0, 0.5), (5.0, 1.0, 0.5)], max_distance=5.0,
    )
    assert projection.complete
    assert projection.mole_fractions == (0.15, 0.25)
    assert projection.temperatures_k == (80.5, 80.5)
