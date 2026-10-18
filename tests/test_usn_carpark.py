import math

import numpy as np
import pytest

from degali.validation.usn_carpark import (
    compare_well_mixed_baseline,
    integrated_positive_source_mass_g,
    read_usn_carpark_mat,
    threshold_metrics,
)


h5py = pytest.importorskip("h5py")


def _char(handle, name, value):
    data = np.asarray([[ord(char)] for char in value], dtype=np.uint16)
    dataset = handle.create_dataset(name, data=data)
    dataset.attrs["MATLAB_class"] = np.bytes_("char")
    return dataset.ref


def _numeric(handle, name, value):
    return handle.create_dataset(name, data=np.asarray(value, dtype=float)).ref


def _mat_file(path):
    with h5py.File(path, "w") as handle:
        refs = handle.create_group("#refs#")
        mfm = handle.create_group("MFM")
        tank = handle.create_group("Tank")
        vent = handle.create_group("Vent")
        sensor = handle.create_group("S")
        source_clock = np.arange(0.1, 12.1, 0.1)
        mfm.create_dataset("time", data=source_clock[:, None])
        mfm.create_dataset(
            "mfr", data=np.where((source_clock >= 1.0) & (source_clock < 5.0), 1.0, 0.0)[None, :],
        )
        mfm.create_dataset("t0", data=[[1.0]])
        tank.create_dataset("time", data=source_clock[:, None])
        tank.create_dataset("p", data=np.linspace(100.0, 10.0, source_clock.size)[None, :])
        tank.create_dataset("T", data=np.full((1, source_clock.size), 20.0))
        vent.create_dataset("time", data=source_clock[:, None])
        vent.create_dataset("vfr", data=np.full((1, source_clock.size), 600.0))

        sensor_clock = np.arange(0.0, 12.0, 1.0)
        count = 2
        for field in ("ID", "pos", "time", "conc"):
            sensor.create_dataset(field, shape=(count, 1), dtype=h5py.ref_dtype)
        for index, sensor_id in enumerate(("D139", "D140")):
            sensor["ID"][index, 0] = _char(handle, f"#refs#/id{index}", sensor_id)
            sensor["pos"][index, 0] = _numeric(
                handle, f"#refs#/pos{index}", [[-3000.0], [2170.0], [0.0]],
            )
            sensor["time"][index, 0] = _numeric(
                handle, f"#refs#/time{index}", sensor_clock[None, :],
            )
            concentration = np.asarray([0, 0, 2, 4, 5, 3, 0, 0, 0, 0, 0, 0], dtype=float)
            sensor["conc"][index, 0] = _numeric(
                handle, f"#refs#/conc{index}", concentration[None, :],
            )


def test_usn_carpark_reader_and_threshold_operator(tmp_path):
    path = tmp_path / "Exp19.mat"
    _mat_file(path)
    run = read_usn_carpark_mat(path)

    assert run.experiment_id == "USN-Exp19"
    assert len(run.room_cloud_sensors) == 2
    assert run.room_cloud_sensors[0].position_raw_mm == (-3000.0, 2170.0, 0.0)
    assert not run.quantitative_lh2_field_validation_allowed
    metric = threshold_metrics(run)
    assert metric.first_arrival_after_source_s == 2.0
    assert metric.exceedance_duration_s == 2.0
    assert metric.observed_peak_vol_percent == 5.0
    assert math.isclose(integrated_positive_source_mass_g(run), 4.0, rel_tol=0.05)


def test_usn_carpark_well_mixed_comparison_is_explicitly_no_fit(tmp_path):
    path = tmp_path / "Exp20.mat"
    _mat_file(path)
    run = read_usn_carpark_mat(path)
    comparison = compare_well_mixed_baseline(run)

    assert comparison.experiment_id == "USN-Exp20"
    assert comparison.model_peak_vol_percent > 0.0
    assert comparison.observed_peak_vol_percent == 5.0
    assert comparison.mae_against_spatial_mean_vol_percent >= 0.0


def test_usn_carpark_rejects_clock_drift(tmp_path):
    path = tmp_path / "Exp21.mat"
    _mat_file(path)
    with h5py.File(path, "r+") as handle:
        handle["Vent/time"][3, 0] += 0.01
    with pytest.raises(ValueError, match="clocks differ"):
        read_usn_carpark_mat(path)
