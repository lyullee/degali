import pytest
import numpy as np

openpyxl = pytest.importorskip("openpyxl")

from degali.validation.preslhy_e34 import (
    E34FillPrehistoryPolicy,
    E34PoolHistory,
    E34SourceWindowPolicy,
    evaporation_window,
    read_e34_dynamic_wind_hydrogen_history,
    read_e34_pool_history,
    read_e34_wind_history,
    source_only_evaporation_candidates,
    source_only_fill_prehistory_candidate,
)


def _workbook(path):
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = "Concrete"
    sheet.cell(7, 1, "X_Value")
    sheet.cell(7, 2, "TA000-00")
    sheet.cell(7, 31, "Sync. Time [s]")
    sheet.cell(7, 33, "m(LH2) [g]")
    for i, (time, mass, temp) in enumerate(((0, 1000, 293), (5, 950, 290), (10, 900, 288), (15, 850, 287)), 8):
        sheet.cell(i, 1, i - 8)
        sheet.cell(i, 2, temp)
        sheet.cell(i, 31, time)
        sheet.cell(i, 33, mass)
    book.save(path)


def test_reads_corrected_lh2_mass_and_predeclared_window(tmp_path):
    path = tmp_path / "20200320-Concrete01-Final.xlsx"
    _workbook(path)
    history = read_e34_pool_history(path)
    window = evaporation_window(history, start_s=0, end_s=15, thermocouples=("TA000-00",))
    assert history.substrate == "concrete"
    assert history.liquid_mass_kg.tolist() == pytest.approx([1.0, .95, .9, .85])
    assert window.evaporation_rate_kg_s == pytest.approx(0.01)
    assert window.mass_slope_r_squared == pytest.approx(1.0)
    assert window.temperature_span_k == pytest.approx(6.0)


def test_rejects_a_non_evaporating_window(tmp_path):
    path = tmp_path / "20200320-Concrete01-Final.xlsx"
    _workbook(path)
    history = read_e34_pool_history(path)
    history.liquid_mass_kg[:] = [0.5, 0.6, 0.7, 0.8]
    with pytest.raises(ValueError, match="no net evaporative"):
        evaporation_window(history, start_s=0, end_s=15)


def test_reads_synchronized_dynamic_wind_and_delay_with_quality_flags(tmp_path):
    path = tmp_path / "20200416-Sand03Wind-Final.xlsx"
    book = openpyxl.Workbook()
    book.active.title = "TE-WG"
    sheet = book.create_sheet("cH2-Wind")
    headers = (
        "Orig.Time [s]", "Sync. Time [s]", "CH35", "CH45", "CH55",
        "Young_T", "Young_U", "Young_V", "Young_W", "Trigger",
        "cH2 Time [s]",
    )
    for column, label in enumerate(headers, 1):
        sheet.cell(6, column, label)
        sheet.cell(5, column, 0.1)
    sheet["C5"] = "defective"
    sheet["C5"].fill = openpyxl.styles.PatternFill(
        fill_type="solid", fgColor="FFFFCCFF",
    )
    for row, time in enumerate((0.0, 0.1, 0.2), 7):
        values = (
            time, time - 0.42, 0.4 + time, 0.5 + time, 0.6 + time,
            295.0, 1.0 + time, 2.0 + time, 0.1, 0.0,
            time - 35.42,
        )
        for column, value in enumerate(values, 1):
            sheet.cell(row, column, value)
    book.save(path)

    history = read_e34_dynamic_wind_hydrogen_history(path)

    assert history.original_time_s.tolist() == pytest.approx([0.0, 0.1, 0.2])
    assert history.synchronized_time_s.tolist() == pytest.approx(
        [-0.42, -0.32, -0.22]
    )
    assert history.hydrogen_sample_time_s.tolist() == pytest.approx(
        [-35.42, -35.32, -35.22]
    )
    assert history.declared_sample_transport_delay_s == pytest.approx(35.0)
    assert history.median_sample_interval_s == pytest.approx(0.1)
    assert history.sensor_quality == {
        "CH35": "workbook_flagged_defective",
        "CH45": "usable",
        "CH55": "usable",
    }
    assert history.wind_quality == "usable"
    assert history.wind_v_m_s.tolist() == pytest.approx([2.0, 2.1, 2.2])

    wind_only = read_e34_wind_history(path)
    assert wind_only.synchronized_time_s.tolist() == pytest.approx(
        [-0.42, -0.32, -0.22]
    )
    assert wind_only.wind_quality == "usable"
    assert wind_only.wind_v_m_s.tolist() == pytest.approx([2.0, 2.1, 2.2])
    assert not hasattr(wind_only, "hydrogen_vol_percent")


def test_source_only_policy_freezes_every_post_fill_decay(tmp_path):
    time = np.arange(0.0, 501.0)
    mass = np.zeros_like(time)
    mass[40:61] = np.linspace(0.0, 2.0, 21)
    mass[61:201] = np.linspace(2.0, 0.02, 140)
    mass[240:261] = np.linspace(0.0, 1.5, 21)
    mass[261:451] = np.linspace(1.5, 0.02, 190)
    history = E34PoolHistory(
        path=tmp_path / "synthetic-ConcreteWind-Final.xlsx",
        time_s=time,
        liquid_mass_kg=mass,
        thermocouple_k={},
        substrate="sand",
    )

    candidates = source_only_evaporation_candidates(
        history,
        policy=E34SourceWindowPolicy(
            smoothing_sigma_s=3.0,
            minimum_peak_prominence_kg=0.2,
            minimum_peak_separation_s=100.0,
        ),
    )

    assert len(candidates) == 2
    assert candidates[0].start_s == pytest.approx(60.0, abs=5.0)
    assert candidates[0].end_s < candidates[1].start_s
    assert candidates[1].start_s == pytest.approx(260.0, abs=5.0)
    assert all(candidate.source_only_selection for candidate in candidates)
    assert all(candidate.smoothed_mass_lost_kg > 1.0 for candidate in candidates)
    assert all(candidate.fitted_evaporation_rate_kg_s > 0.0 for candidate in candidates)
    assert all(
        "liquid_mass_loss_is_not_a_direct_measurement_of_atmospheric_vapour_flux"
        in candidate.manual_review_reasons
        for candidate in candidates
    )
    assert all(
        "porous_substrate_mass_artefact_requires_review"
        in candidate.manual_review_reasons
        for candidate in candidates
    )


def test_source_only_policy_rejects_an_irregular_clock(tmp_path):
    history = E34PoolHistory(
        path=tmp_path / "irregular.xlsx",
        time_s=np.asarray([0.0, 1.0, 2.2, 3.2]),
        liquid_mass_kg=np.asarray([0.0, 1.0, 0.5, 0.0]),
        thermocouple_k={},
        substrate="concrete",
    )
    with pytest.raises(ValueError, match="uniform time grid"):
        source_only_evaporation_candidates(history)


def test_source_only_fill_prehistory_uses_declared_fraction_before_peak(tmp_path):
    time = np.arange(0.0, 301.0)
    mass = np.zeros_like(time)
    mass[50:101] = np.linspace(0.0, 2.0, 51)
    mass[101:261] = np.linspace(2.0, 0.02, 160)
    history = E34PoolHistory(
        path=tmp_path / "synthetic-ConcreteWind-Final.xlsx",
        time_s=time,
        liquid_mass_kg=mass,
        thermocouple_k={},
        substrate="concrete",
    )
    source_policy = E34SourceWindowPolicy(
        smoothing_sigma_s=2.0,
        minimum_peak_prominence_kg=0.2,
        minimum_peak_separation_s=100.0,
    )

    candidate = source_only_fill_prehistory_candidate(
        history,
        window_ordinal=1,
        source_window_policy=source_policy,
        prehistory_policy=E34FillPrehistoryPolicy(
            smoothing_sigma_s=2.0,
            relative_fill_onset_fraction=0.10,
            minimum_prehistory_s=10.0,
        ),
    )

    assert candidate.source_only_selection is True
    assert candidate.prehistory_start_s == pytest.approx(55.0, abs=3.0)
    assert candidate.score_start_s == pytest.approx(100.0, abs=3.0)
    assert candidate.prehistory_duration_s >= 40.0
    assert candidate.full_duration_s == pytest.approx(
        candidate.prehistory_duration_s + candidate.score_duration_s
    )
    assert "prehistory_source_rate_requires_a_separate_declared_assumption" in (
        candidate.limitations
    )
