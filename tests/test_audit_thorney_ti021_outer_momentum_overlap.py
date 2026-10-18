from tools.audit_thorney_ti021_outer_momentum_overlap import build_audit


def test_outer_momentum_overlap_audit_closes_handoff_but_not_grid_convergence():
    record, arrays = build_audit()

    assert record["status"] == "accepted_outer_momentum_and_overlap_component"
    assert record["all_gates_pass"]
    assert record["configuration"]["observation_value_fields_accessed"] == []
    assert record["overlap_receptors"]["maximum_peak_relative_change"] < 0.10
    assert record["overlap_receptors"]["maximum_dose_relative_change"] < 0.10
    assert record["separate_grid_change_diagnostic"]["status"] == (
        "diagnostic_rejects_2p5_to_5m_as_converged"
    )
    assert record["separate_grid_change_diagnostic"][
        "maximum_peak_relative_change"
    ] > 0.10
    assert record["separate_grid_change_diagnostic"][
        "maximum_dose_relative_change"
    ] > 0.10
    assert arrays["reference_traces"].shape == arrays["outer_traces"].shape

