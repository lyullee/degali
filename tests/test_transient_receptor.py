import math

import numpy as np
import pytest

from degali.addons.transient_receptor import (
    FixedReceptor,
    SteadyPlumeTable,
    SourceHistory,
    WindHistory,
    first_order_sensor_response,
    meteorological_from_to_math_radians,
    replay_fixed_receptors,
    replay_fixed_receptors_with_source_history,
)


def test_meteorological_west_wind_transports_toward_east():
    assert meteorological_from_to_math_radians(270.0) == pytest.approx(0.0)
    assert meteorological_from_to_math_radians(180.0) == pytest.approx(
        math.pi / 2.0
    )


def test_first_order_response_reaches_ninety_percent_at_t90():
    result = first_order_sensor_response(
        [0.0, 6.0], [1.0, 1.0], t90_s=6.0, initial_value=0.0
    )
    assert result == pytest.approx([0.0, 0.9])
    assert first_order_sensor_response(
        [0.0, 1.0], [0.1, 0.2], t90_s=None
    ) == pytest.approx([0.1, 0.2])


def test_steady_table_interpolates_wind_and_rejects_extrapolation():
    table = SteadyPlumeTable(
        [2.0, 4.0],
        [
            lambda x, y, z: 0.1 + 0.01 * x,
            lambda x, y, z: 0.3 + 0.01 * x,
        ],
    )
    assert table(3.0, 5.0, 0.0, 1.0) == pytest.approx(0.25)
    with pytest.raises(ValueError, match="outside the steady table"):
        table(1.0, 5.0, 0.0, 1.0)


def test_replay_rotates_fixed_receptor_and_applies_sensor_lag():
    history = WindHistory(
        time_s=[0.0, 3.0, 6.0],
        speed_m_s=[2.0, 2.0, 2.0],
        direction_from_deg=[270.0, 270.0, 180.0],
    )
    receptor = FixedReceptor("east", 30.0, 0.0, 1.0, response_t90_s=6.0)

    def plume(speed, alongwind, crosswind, height):
        assert speed == pytest.approx(2.0)
        return math.exp(-0.5 * (crosswind / 5.0) ** 2)

    trace = replay_fixed_receptors(history, [receptor], plume)[0]
    assert trace.alongwind_m == pytest.approx([30.0, 30.0, 0.0], abs=1.0e-12)
    assert trace.crosswind_m == pytest.approx([0.0, 0.0, -30.0])
    assert trace.true_mole_fraction == pytest.approx([1.0, 1.0, 0.0])
    assert trace.indicated_mole_fraction == pytest.approx(
        [0.0, 1.0 - math.sqrt(0.1), math.sqrt(0.1) - 0.1]
    )
    stats = trace.statistics(0.0, 3.0)
    assert stats.sample_count == 2
    assert stats.true_maximum == pytest.approx(1.0)


def test_replay_validates_histories_receptors_and_outputs():
    with pytest.raises(ValueError, match="strictly increasing"):
        WindHistory([0.0, 0.0], [1.0, 1.0], [0.0, 0.0]).arrays()
    bad = FixedReceptor("bad", 1.0, 0.0, -1.0)
    with pytest.raises(ValueError, match="height"):
        replay_fixed_receptors(
            WindHistory([0.0], [1.0], [270.0]), [bad],
            lambda speed, x, y, z: 0.0,
        )
    with pytest.raises(ValueError, match="invalid mole fraction"):
        replay_fixed_receptors(
            WindHistory([0.0], [1.0], [270.0]),
            [FixedReceptor("ok", 1.0, 0.0, 0.0)],
            lambda speed, x, y, z: np.nan,
        )


def test_source_history_packet_quadrature_conserves_mass_and_carries_state():
    history = SourceHistory(
        time_s=[0.0, 1.0, 2.0],
        mass_rate_kg_s=[2.0, 1.0, 0.0],
        pressure_pa=[100_000.0, 90_000.0, 80_000.0],
        temperature_k=[25.0, 26.0, 27.0],
        source_direction_to_deg=[350.0, 10.0, 30.0],
        rate_operator="linear",
    )
    packets = history.to_packets(subdivisions_per_interval=2)
    assert sum(packet.mass_kg for packet in packets) == pytest.approx(2.0)
    assert packets[0].release_time_s == pytest.approx(0.25)
    assert packets[0].pressure_pa == pytest.approx(97_500.0)
    assert packets[0].temperature_k == pytest.approx(25.25)
    assert packets[0].source_direction_to_deg == pytest.approx(355.0)
    assert history.released_mass_kg() == pytest.approx(2.0)


def test_source_history_replay_preserves_causal_travel_and_fixed_source_axis():
    source = SourceHistory(
        time_s=[0.0, 1.0, 2.0],
        mass_rate_kg_s=[1.0, 0.0, 0.0],
        source_direction_to_deg=[90.0, 90.0, 90.0],
    )
    wind = WindHistory(
        time_s=[0.0, 1.0, 2.0],
        speed_m_s=[1.0, 1.0, 1.0],
        direction_from_deg=[270.0, 270.0, 270.0],
    )
    receptor = FixedReceptor("east", 1.0, 0.0, 0.0)
    seen_axes = []

    def packet_kernel(packet, age_s, relative_east, relative_north, height_m):
        seen_axes.append(packet.source_direction_to_deg)
        # A deliberately narrow test kernel: the packet is observed only when
        # its causally advected centre reaches the receptor.
        return packet.mass_kg if abs(relative_east) < 1.0e-12 else 0.0

    trace = replay_fixed_receptors_with_source_history(
        source, wind, [receptor], packet_kernel,
        observation_time_s=[0.0, 0.5, 1.5, 2.0],
    )[0]
    assert trace.true_mole_fraction == pytest.approx([0.0, 0.0, 1.0, 0.0])
    assert trace.wind_direction_from_deg == pytest.approx([270.0] * 4)
    assert seen_axes and all(value == pytest.approx(90.0) for value in seen_axes)


def test_source_history_requires_explicit_shutdown_and_wind_coverage():
    with pytest.raises(ValueError, match="explicit zero rate"):
        SourceHistory([0.0, 1.0], [1.0, 1.0]).to_packets()
    source = SourceHistory([0.0, 1.0], [1.0, 0.0])
    wind = WindHistory([0.0, 0.25], [1.0, 1.0], [270.0, 270.0])
    with pytest.raises(ValueError, match="cover every source packet release"):
        replay_fixed_receptors_with_source_history(
            source, wind, [FixedReceptor("east", 1.0, 0.0, 0.0)],
            lambda packet, age, x, y, z: 0.0,
        )
