from tools.audit_thorney_ti021_checkpoint_equivalence import (
    FINAL_FIELDS,
    HANDOFF_FIELDS,
    RECEPTOR_FIELDS,
    SCHEMA,
)


def test_checkpoint_equivalence_contract_covers_all_published_state_groups():
    assert SCHEMA.endswith(".v1")
    assert "concentration_vol_percent" in RECEPTOR_FIELDS
    assert "maximum_concentration_vol_percent" in FINAL_FIELDS
    assert {
        "velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s",
    }.issubset(FINAL_FIELDS)
    assert "species_mass_flux_kg_m2_s" in HANDOFF_FIELDS
    assert "sensible_enthalpy_flux_w_m2" in HANDOFF_FIELDS
