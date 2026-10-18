from tools.audit_thorney_ti021_outer_inlet_coupling import build_audit


def test_outer_inlet_coupling_closes_binary_species_energy_and_continuity():
    record, arrays = build_audit()

    assert record["status"] == (
        "accepted_outer_binary_boundary_component_pending_momentum_overlap"
    )
    assert record["all_gates_pass"]
    assert record["configuration"]["observation_value_fields_accessed"] == []
    assert abs(record["ledger"]["species_mass_residual_kg"]) < 1.0e-10
    assert abs(record["ledger"]["sensible_enthalpy_residual_j"]) < 1.0e-6
    assert record["ledger"]["continuity_residual_kg_m3_s"] < 1.0e-7
    assert arrays["inlet_species_flux"].shape == (8, 128)
    assert arrays["first_cell_species"].shape == (8, 128)
