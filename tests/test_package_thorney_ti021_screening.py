import json
from pathlib import Path

from tools.package_thorney_ti021_screening import package


ROOT = Path(__file__).resolve().parents[1]


def test_screening_package_reports_full_horizon_and_withholds_promotion(tmp_path):
    run = ROOT / "outputs" / "thorney-ti021-coarse-600s-dx2p5m-dt0p0625-checkpoint1s-tol1e3-2026-10-11"
    score = ROOT / "outputs" / "thorney-ti021-coarse-600s-score-2026-10-12"
    result = package(run, tmp_path / "screening", score=score)
    assert result["run"]["duration_s"] == 600.0
    assert result["run"]["cell_size_m"] == 2.5
    assert result["concentration_unit"] == "vol%"
    assert result["hazard"]["threshold_vol_percent"] == 4.0
    assert result["promotion_allowed"] is False
    assert result["independent_error_metrics"]["concentration_unit"] == "vol%"
    assert result["independent_error_metrics"]["dose_mae_vol_percent_s"] is not None
    assert result["independent_error_metrics"]["dose_pearson_r"] is not None
    assert "pressure_pa_gauge" in result["physics"]["state_fields_present"]
    assert result["physics"]["conservation_residuals"]["species_mass_relative"] is not None
    assert result["physics"]["thermal_buoyancy_transition"]["status"] == "not_exercised_in_this_screening_run"
    assert (tmp_path / "screening" / "ti021-600s-screening-package.json").is_file()
