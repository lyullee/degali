from tools.audit_thorney_ti021_nested_handoff import build_audit


def test_ti021_nested_handoff_audit_passes_all_flux_integrals_without_observations():
    record, arrays = build_audit()

    assert record["status"] == "accepted_conservative_handoff_remap_operator_only"
    assert record["all_gates_pass"]
    assert record["configuration"]["observation_value_fields_accessed"] == []
    assert record["diagnostics"]["maximum_relative_integral_residual"] < 1.0e-12
    assert arrays["source_species"].shape == (5, 32, 160)
    assert arrays["target_species"].shape == (5, 8, 128)
