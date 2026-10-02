import pytest

openpyxl = pytest.importorskip("openpyxl")

from degali.validation.preslhy_e34 import evaporation_window, read_e34_pool_history


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
