import csv
from pathlib import Path

import numpy as np
import pytest

from degali.validation.h2safe import (
    LAB1_VENT_BOUNDARIES,
    LAB2_SUPPLY_BOUNDARIES,
    read_h2safe_ascii_stl_bounds,
    read_h2safe_run,
    read_h2safe_step_product_names,
)


def _write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerows(rows)


def _lab1(root: Path) -> Path:
    lab = root / "Lab1_Small"
    _write_csv(
        lab / "sensors.csv",
        [["id", "X", "Y", "Z"], ["P1", "0", "1", "2"], ["p2", "3", "4", "5"]],
    )
    _write_csv(
        lab / "test1.csv",
        [["time", "p1", "P2"], ["0", "0.0", "#N/A"], ["1", "1.2", "#N/A"]],
    )
    return lab


def test_h2safe_reader_preserves_native_values_and_missing_channels(tmp_path):
    _lab1(tmp_path)
    run = read_h2safe_run(tmp_path, lab_id="Lab1_Small", test_id="test1")

    assert run.event_id == "H2SAFE-Lab1_Small-test1"
    assert run.sample_interval_s == 1.0
    assert run.sensor_ids == ("p1", "p2")
    assert run.release.release_duration_basis == "declared_seconds"
    assert run.release.solver_position_m == (
        -5.186701676, -4.322127, 1.580460097,
    )
    assert run.registry[0].solver_position_m == (0.0, 2.0, 1.0)
    assert run.fully_missing_sensor_ids == ("p2",)
    assert run.native_maximum == 1.2
    assert np.isnan(run.concentration_raw[0, 1])
    assert not run.vol_percent_conversion_allowed
    assert not run.quantitative_hydrogen_validation_allowed


def test_h2safe_reader_rejects_sensor_absent_from_registry(tmp_path):
    lab = _lab1(tmp_path)
    _write_csv(
        lab / "test1.csv",
        [["time", "p1", "p3"], ["0", "0", "1"], ["1", "1", "2"]],
    )
    with pytest.raises(ValueError, match="absent from the registry"):
        read_h2safe_run(tmp_path, lab_id="Lab1_Small", test_id="test1")


def test_h2safe_reader_rejects_clock_duplicates(tmp_path):
    lab = _lab1(tmp_path)
    _write_csv(
        lab / "test1.csv",
        [["time", "p1", "p2"], ["0", "0", "1"], ["0", "1", "2"]],
    )
    with pytest.raises(ValueError, match="strictly increasing"):
        read_h2safe_run(tmp_path, lab_id="Lab1_Small", test_id="test1")


def test_h2safe_ascii_stl_bounds_and_containment(tmp_path):
    path = tmp_path / "lab.STL"
    path.write_text(
        "solid lab\n"
        "facet normal 0 0 1\nouter loop\n"
        "vertex -10 0 -9\nvertex 1 0 -9\nvertex 1 11 1\n"
        "endloop\nendfacet\nendsolid lab\n",
        encoding="ascii",
    )
    bounds = read_h2safe_ascii_stl_bounds(path)
    assert bounds.vertex_record_count == 3
    assert bounds.minimum_m == (-10.0, 0.0, -9.0)
    assert bounds.maximum_m == (1.0, 11.0, 1.0)
    assert bounds.contains((-5.0, 1.0, -4.0))
    assert not bounds.contains((-11.0, 1.0, -4.0))


def test_h2safe_step_product_name_inventory(tmp_path):
    path = tmp_path / "lab.STEP"
    path.write_text(
        "#1 = PRODUCT ( 'east_supply_vent', 'east_supply_vent', '', (#2) );\n"
        "#3 = PRODUCT ( 'room_box', 'room_box', '', (#4) );\n"
        "#5 = PRODUCT ( 'east_supply_vent', 'east_supply_vent', '', (#6) );\n",
        encoding="ascii",
    )
    assert read_h2safe_step_product_names(path) == (
        "east_supply_vent", "room_box",
    )


def test_h2safe_hvac_tables_preserve_reported_flow_inputs():
    inlet = sum(
        item.volumetric_flow_m3_s
        for item in LAB1_VENT_BOUNDARIES if item.role == "inlet"
    )
    outlet = sum(
        item.volumetric_flow_m3_s
        for item in LAB1_VENT_BOUNDARIES if item.role == "outlet"
    )
    assert inlet > 0.0 and outlet > 0.0
    assert np.isclose(inlet, outlet, rtol=0.01)
    assert not np.isclose(inlet, outlet, rtol=1.0e-4)
    assert len(LAB2_SUPPLY_BOUNDARIES) == 8
    assert sum(item.speed_m_s == 0.0 for item in LAB2_SUPPLY_BOUNDARIES) == 3
