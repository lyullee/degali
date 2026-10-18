import numpy as np
import pytest

from tools.score_thorney_ti021_inner_transient import (
    _resolution_classification,
    _run_completion_metadata,
    trailing_boxcar_average,
)


def test_resolution_classification_does_not_mislabel_dx2p5_as_dx5():
    coarse = _resolution_classification(5.0)
    fine = _resolution_classification(2.5)
    finest = _resolution_classification(1.25)
    assert coarse["status"] == "scored_no_fit_inner_dx5_not_grid_converged"
    assert fine["status"] == "scored_no_fit_inner_dx2p5_pending_finest_grid"
    assert fine["wall_height_cells"] == 2.0
    assert fine["0p4m_and_2p4m_sensors_share_vertical_cell"] is True
    assert finest["wall_height_cells"] == 4.0
    assert finest["0p4m_and_2p4m_sensors_share_vertical_cell"] is False


def test_trailing_boxcar_average_integrates_piecewise_linear_trace():
    time = np.asarray((0.0, 0.5, 1.0, 1.5))
    trace = np.asarray((0.0, 1.0, 2.0, 3.0))

    averaged = trailing_boxcar_average(time, trace, window_s=1.0)

    assert averaged[0] == 0.0
    assert averaged[1] == pytest.approx(0.5)
    assert averaged[2] == pytest.approx(1.0)
    assert averaged[3] == pytest.approx(2.0)


def test_trailing_boxcar_average_preserves_multiple_trace_shape():
    time = np.asarray((0.0, 1.0, 2.0))
    traces = np.asarray(((2.0, 2.0, 2.0), (0.0, 1.0, 0.0)))

    averaged = trailing_boxcar_average(time, traces, window_s=0.5)

    assert averaged.shape == traces.shape
    assert averaged[0] == pytest.approx((2.0, 2.0, 2.0))
    assert averaged[1, 1] == pytest.approx(0.75)


def test_completed_qualification_is_distinct_from_complete_observation_horizon():
    qualification = _run_completion_metadata(
        "accepted_inner_transient_qualification"
    )
    complete = _run_completion_metadata("accepted_completed_inner_transient")
    local_patch = _run_completion_metadata(
        "accepted_local_patch_transient_qualification"
    )

    assert qualification == {
        "scoring_opened_after_completed_run": True,
        "observation_horizon_complete": False,
    }
    assert complete == {
        "scoring_opened_after_completed_run": True,
        "observation_horizon_complete": True,
    }
    assert local_patch == {
        "scoring_opened_after_completed_run": True,
        "observation_horizon_complete": False,
    }
    with pytest.raises(ValueError, match="has not completed"):
        _run_completion_metadata("running")
