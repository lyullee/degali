import numpy as np
import pytest

import degali.validation.preslhy_e31 as e31

openpyxl = pytest.importorskip("openpyxl")

from degali.validation.preslhy_e31 import (
    axis_temperature_envelope,
    compare_axis_temperature_envelope,
    compare_blowdown_pressure_stations,
    homogeneous_equilibrium_source_bound,
    ideal_choked_hydrogen_source_bound,
    nozzle_pressure_rise,
    read_e31_high_pressure_run,
    read_e31_pressure_run,
    valve_open_interval,
    vessel_inventory_change,
    vessel_inventory_release_timing,
)
from degali.lh2 import run_lh2_near_field_research


def _workbook(path, *, include_t7=True):
    book = openpyxl.Workbook()
    pressure = book.active
    pressure.title = "sample-Press"
    pressure.append(["dNoz [mm]", None, None, None, None, "Pini [bar]"])
    pressure.append([4.0, None, None, None, None, 200.0])
    pressure.append([])
    pressure.append([None, "Unit", "bar", "bar", "Volt"])
    pressure.append([])
    pressure.append(["Time [s]", "X_Value", "Druck_1", "Druck_2", "Valve-Relay"])
    pressure.append([-0.1, 0, 200.0, 0.5, 0.0])
    pressure.append([0.0, 1, 199.0, 198.0, 10.0])
    pressure.append([0.1, 2, 198.0, 197.0, 10.0])

    temperature = book.create_sheet("sample-Temp")
    temperature.append([])
    temperature.append([])
    temperature.append([])
    temperature.append([None, "Unit", "Kelvin", "Kelvin", "Kelvin", "Kelvin"])
    temperature.append([])
    header = ["Time [s]", "X_Value", "T5", "T6"]
    if include_t7:
        header.append("T7")
    temperature.append(header)
    temperature.append([-0.1, 0, 290.0, 291.0, 292.0] if include_t7 else [-0.1, 0, 290.0, 291.0])
    temperature.append([0.0, 1, 260.0, 270.0, 280.0] if include_t7 else [0.0, 1, 260.0, 270.0])
    temperature.append([0.1, 2, 255.0, 275.0, 285.0] if include_t7 else [0.1, 2, 255.0, 275.0])
    book.save(path)


def test_reads_only_synchronized_pressure_and_axis_temperature(tmp_path):
    path = tmp_path / "e31.xlsx"
    _workbook(path)

    run = read_e31_high_pressure_run(path)

    assert run.path == path
    assert run.nozzle_diameter_mm == pytest.approx(4.0)
    assert run.initial_vessel_pressure_bar == pytest.approx(200.0)
    assert np.allclose(run.pressure.time_s, [-0.1, 0.0, 0.1])
    assert np.allclose(run.pressure.nozzle_pressure_bar, [0.5, 198.0, 197.0])
    assert np.allclose(run.axis_temperature.t6_k, [291.0, 270.0, 275.0])
    assert run.axis_temperature.t5_k[-1] == pytest.approx(255.0)
    assert run.axis_temperature.axial_distances_m == (0.25, 0.75, 1.75)
    assert run.axis_temperature.cold_bath_offset_applied_k == 0.0


def test_reads_pressure_only_without_requiring_the_temp_sheet(tmp_path):
    path = tmp_path / "e31_pressure_only.xlsx"
    _workbook(path)
    book = openpyxl.load_workbook(path)
    del book["sample-Temp"]
    book.save(path)

    run = read_e31_pressure_run(path)

    assert run.nozzle_diameter_mm == pytest.approx(4.0)
    assert np.allclose(run.pressure.vessel_pressure_bar, [200.0, 199.0, 198.0])


def test_large_xml_fast_reader_preserves_the_pressure_history(tmp_path, monkeypatch):
    path = tmp_path / "e31_fast.xlsx"
    _workbook(path)
    monkeypatch.setattr(e31, "_FAST_XML_MINIMUM_BYTES", 0)

    run = e31.read_e31_pressure_run(path)

    assert np.allclose(run.pressure.time_s, [-0.1, 0.0, 0.1])
    assert np.allclose(run.pressure.vessel_pressure_bar, [200.0, 199.0, 198.0])
    assert np.allclose(run.pressure.nozzle_pressure_bar, [0.5, 198.0, 197.0])
    assert np.allclose(run.pressure.valve_relay_v, [0.0, 10.0, 10.0])


def test_rejects_workbook_without_a_required_axis_channel(tmp_path):
    path = tmp_path / "e31_missing_t7.xlsx"
    _workbook(path, include_t7=False)

    with pytest.raises(ValueError, match="T7"):
        read_e31_high_pressure_run(path)


def test_reads_the_optional_welded_release_line_thermocouple(tmp_path):
    path = tmp_path / "e31_release_line.xlsx"
    _workbook(path)
    book = openpyxl.load_workbook(path)
    sheet = book["sample-Temp"]
    sheet.cell(row=6, column=7, value="T-Stück")
    sheet.cell(row=7, column=7, value=82.0)
    sheet.cell(row=8, column=7, value=51.0)
    sheet.cell(row=9, column=7, value=44.0)
    book.save(path)

    run = read_e31_high_pressure_run(path)

    assert run.release_line_temperature is not None
    assert run.release_line_temperature.workbook_heading == "T-Stück"
    assert np.allclose(run.release_line_temperature.temperature_k, [82.0, 51.0, 44.0])


def test_relay_interval_and_temperature_envelope_are_fixed_by_the_record(tmp_path):
    path = tmp_path / "e31_envelope.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)

    interval = valve_open_interval(run)
    envelope = axis_temperature_envelope(run, interval)

    assert interval.start_s == pytest.approx(0.0)
    assert interval.end_s == pytest.approx(0.1)
    assert interval.samples == 2
    assert interval.threshold_v == pytest.approx(5.5)
    assert envelope.samples == 2
    assert envelope.minimum_k == pytest.approx((255.0, 270.0, 280.0))
    assert envelope.percentile_05_k == pytest.approx((255.25, 270.25, 280.25))
    assert envelope.median_k == pytest.approx((257.5, 272.5, 282.5))


def test_relay_interval_rejects_a_non_varying_signal(tmp_path):
    path = tmp_path / "e31_flat_relay.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)
    run.pressure.valve_relay_v[:] = 0.0

    with pytest.raises(ValueError, match="no resolvable"):
        valve_open_interval(run)


def test_nozzle_pressure_rise_keeps_the_relay_delay_explicit(tmp_path):
    path = tmp_path / "e31_nozzle_response.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)
    interval = valve_open_interval(run)

    response = nozzle_pressure_rise(run, interval, required_rise_bar=5.0)

    assert response.baseline_pressure_bar == pytest.approx(0.5)
    assert response.threshold_pressure_bar == pytest.approx(5.5)
    assert response.response_time_s == pytest.approx(0.0)
    assert response.delay_after_relay_s == pytest.approx(0.0)
    with pytest.raises(ValueError, match="positive"):
        nozzle_pressure_rise(run, interval, required_rise_bar=0.0)
    with pytest.raises(ValueError, match="does not cross"):
        nozzle_pressure_rise(run, interval, required_rise_bar=1000.0)


def test_axis_temperature_comparison_keeps_the_time_envelope_visible(tmp_path):
    path = tmp_path / "e31_compare.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)
    envelope = axis_temperature_envelope(run, valve_open_interval(run))

    comparison = compare_axis_temperature_envelope(
        envelope, (255.0, 272.5, 282.5)
    )

    assert comparison.distances_m == (0.25, 0.75, 1.75)
    assert comparison.residual_to_percentile_05_k == pytest.approx(
        (0.25, -2.25, -2.25)
    )
    assert comparison.inside_recorded_minimum_median == (True, True, True)
    with pytest.raises(ValueError, match="three positive"):
        compare_axis_temperature_envelope(envelope, (100.0, 200.0))


def test_blowdown_pressure_comparison_uses_declared_source_time_stations(tmp_path):
    path = tmp_path / "e31_pressure_stations.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)

    comparison = compare_blowdown_pressure_stations(
        run,
        source_start_time_s=0.0,
        model_time_s=[0.0, 0.1],
        model_pressure_pa=[199.0e5, 198.0e5],
        model_tank_temperature_k=[80.0, 79.0],
        stations_after_source_start_s=[0.0, 0.1],
    )

    assert comparison.measured_pressure_bar == pytest.approx((199.0, 198.0))
    assert comparison.pressure_residual_bar == pytest.approx((0.0, 0.0))
    assert comparison.model_tank_temperature_k == pytest.approx((80.0, 79.0))
    with pytest.raises(ValueError, match="outside"):
        compare_blowdown_pressure_stations(
            run,
            source_start_time_s=0.0,
            model_time_s=[0.0, 0.1],
            model_pressure_pa=[199.0e5, 198.0e5],
            model_tank_temperature_k=[80.0, 79.0],
            stations_after_source_start_s=[0.2],
        )


def test_vessel_inventory_change_uses_declared_volume_and_temperature(tmp_path):
    path = tmp_path / "e31_inventory.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)
    interval = valve_open_interval(run)

    inventory = vessel_inventory_change(
        run,
        interval,
        vessel_volume_m3=0.002815,
        storage_temperature_k=80.0,
    )

    assert inventory.samples == 2
    assert inventory.initial_mass_kg > inventory.final_mass_kg > 0.0
    assert inventory.released_mass_kg > 0.0
    assert inventory.mean_release_rate_kg_s == pytest.approx(
        inventory.released_mass_kg / 0.1
    )

    timing = vessel_inventory_release_timing(
        run,
        interval,
        vessel_volume_m3=0.002815,
        storage_temperature_k=80.0,
        released_fraction=0.95,
    )
    assert timing.time_after_open_s == pytest.approx(0.1)
    assert timing.released_mass_kg >= 0.95 * timing.total_released_mass_kg


def test_explicit_80k_temperature_builds_an_ideal_choked_source_upper_bound(tmp_path):
    path = tmp_path / "e31_source_bound.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)

    source = ideal_choked_hydrogen_source_bound(run, storage_temperature_k=80.0)

    assert source.storage_pressure_bar == pytest.approx(200.0)
    assert source.orifice_diameter_mm == pytest.approx(4.0)
    assert source.ideal_mass_flow_kg_s > 0.0
    assert source.throat_pressure_bar < source.storage_pressure_bar
    assert source.throat_temperature_k < source.storage_temperature_k
    assert source.atmospheric_diameter_m > 0.004
    assert source.atmospheric_temperature_k > 0.0
    assert source.atmospheric_density_kg_m3 > 0.0
    assert source.atmospheric_velocity_m_s > source.throat_velocity_m_s


def test_80k_upper_source_reaches_a_phase_safe_hydrogen_gas_handoff(tmp_path):
    path = tmp_path / "e31_phase_safe_source.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)

    handoff = homogeneous_equilibrium_source_bound(run, storage_temperature_k=80.0)

    assert 0.0 < handoff.postflash.postflash_quality < 1.0
    assert handoff.postflash.liquid_mass_flow > 0.0
    assert handoff.formation_distance > 0.0
    assert 0.0 < handoff.source.mass_fraction < 1.0
    assert handoff.hydrogen_mass_residual < 1e-12
    assert handoff.total_mass_residual < 1e-12
    assert handoff.momentum_residual < 1e-12
    assert handoff.energy_residual < 1e-12


def test_80k_phase_safe_upper_bound_enters_the_conserved_nearfield_path(tmp_path):
    pytest.importorskip("CoolProp")
    path = tmp_path / "e31_nearfield_bound.xlsx"
    _workbook(path)
    run = read_e31_high_pressure_run(path)
    handoff = homogeneous_equilibrium_source_bound(run, storage_temperature_k=80.0)

    result = run_lh2_near_field_research(
        handoff.source,
        ambient_temperature=293.0,
        maximum_distance=0.25,
        radial_points=41,
        maximum_step=0.001,
        relative_tolerance=2.0e-6,
        establishment="scalar_peak",
    )

    assert result.solution.S[-1] == pytest.approx(0.25)
    assert result.solution.temperature[-1] > handoff.source.temperature
    assert result.maximum_boundary_residual < 1e-8
    assert result.maximum_species_drift < 2e-4
    assert result.maximum_energy_drift < 2e-4
    assert result.conservative
