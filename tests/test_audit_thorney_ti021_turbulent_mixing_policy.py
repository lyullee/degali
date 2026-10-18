import numpy as np
import pytest

from tools.audit_thorney_ti021_turbulent_mixing_policy import (
    SCHEMA,
    build_record,
)


def test_ti021_mixing_policy_is_concentration_blind_and_keeps_full_prior_envelope():
    record = build_record(heights_m=np.asarray((0.5, 1.5, 2.5)))

    assert record["schema"] == SCHEMA
    assert record["status"] == (
        "accepted_concentration_blind_turbulent_mixing_envelope_component"
    )
    assert record["all_gates_pass"]
    assert record["evidence_provenance"][
        "concentration_observation_fields_accessed"
    ] == []
    assert len(record["frozen_policy"]["ensemble"]) == 6
    assert record["frozen_policy"]["required_turbulent_schmidt_numbers"] == [
        0.2, 0.7, 1.3,
    ]
    assert record["frozen_policy"]["required_turbulent_prandtl_numbers"] == [
        0.74, 1.0,
    ]
    assert record["profile"]["eddy_viscosity_m2_s"] == pytest.approx(
        (0.041, 0.123, 0.205)
    )


def test_ti021_mixing_policy_rejects_invalid_vertical_grid():
    with pytest.raises(ValueError, match="strictly increasing"):
        build_record(heights_m=np.asarray((0.5, 0.5)))
    with pytest.raises(ValueError, match="positive"):
        build_record(heights_m=np.asarray((0.0, 0.5)))
