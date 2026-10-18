import runpy

import numpy as np
import pytest


TOOL = runpy.run_path("tools/audit_hytunnel_plume_domain_sensitivity.py")


def _case(scale=1.0):
    return {
        "common_maximum_field": [[[1.0 * scale, 4.0 * scale]]],
        "common_maximum_hydrogen_vol_percent": 4.0 * scale,
        "common_4vol_maximum_volume_m3": (
            TOOL["GRID_SPACING_M"] ** 3 * scale
        ),
        "time_s": [0.0, 0.01, 0.02],
        "receptor_metrics": {
            "sensor": {
                "peak_vol_percent": 2.0 * scale,
                "dose_vol_percent_s": 0.02 * scale,
                "arrival_4vol_s": None,
                "duration_above_4vol_s": 0.0,
            },
        },
        "maximum_relative_mass_residual": 1.0e-12,
        "maximum_relative_energy_residual": 1.0e-12,
        "active_scalar_applicability_passed": True,
        "boundary_clearance": {
            "boundary_clearance_gate_passed": True,
        },
    }


def test_domain_classifier_passes_below_frozen_limits():
    result = TOOL["classify_domain_sensitivity"](_case(1.04), _case(1.0))
    assert result["passed"] is True
    assert result["promotion_allowed"] is False
    assert result["maximum_primary_relative_change"] == pytest.approx(
        1.0 - 1.0 / 1.04
    )


def test_domain_classifier_fails_primary_output_change():
    result = TOOL["classify_domain_sensitivity"](_case(1.20), _case(1.0))
    assert result["passed"] is False
    assert result["status"] == "same_grid_plume_domain_sensitivity_gate_failed"


def test_domain_cases_keep_spacing_and_source_relative_receptors():
    published = TOOL["_case_specification"]("published_box")
    expanded = TOOL["_case_specification"]("expanded_box")
    assert published["cells"] == (56, 20, 20)
    assert expanded["cells"] == (70, 30, 30)
    for domain, cells in (
        (published["domain_m"], published["cells"]),
        (expanded["domain_m"], expanded["cells"]),
    ):
        assert tuple(d / n for d, n in zip(domain, cells)) == pytest.approx(
            (TOOL["GRID_SPACING_M"],) * 3
        )
    published_offsets = [
        (
            position[0] - published["source_x_m"],
            position[1] - published["centre_y_m"],
            position[2] - published["centre_z_m"],
        )
        for _label, position in published["receptors"]
    ]
    expanded_offsets = [
        (
            position[0] - expanded["source_x_m"],
            position[1] - expanded["centre_y_m"],
            position[2] - expanded["centre_z_m"],
        )
        for _label, position in expanded["receptors"]
    ]
    assert np.allclose(expanded_offsets, published_offsets, atol=1.0e-14)


def test_domain_audit_explicitly_selects_cached_direct_projection():
    specification = TOOL["_case_specification"]("expanded_box")
    config = TOOL["_case_config"](specification)
    assert config.variable_density_projection_linear_solver == "cached_direct"
    assert config.flow_boundary_scheme == (
        "characteristic_multidimensional_extrapolated"
    )


def test_single_case_sensitivity_parameters_change_only_declared_inputs():
    specification = TOOL["_case_specification"](
        "published_box", maximum_density_departure=0.05,
    )
    config = TOOL["_case_config"](
        specification, time_steps=4, effective_schmidt_number=1.0,
    )
    baseline = TOOL["_case_specification"]("published_box")
    assert specification["source_x_m"] > baseline["source_x_m"]
    assert config.effective_schmidt_number == 1.0
    assert config.duration_s == pytest.approx(4 * TOOL["ACTIVE_TIME_STEP_S"])


def test_case_config_can_select_local_sgs_gradient_diffusion():
    specification = TOOL["_case_specification"]("published_box")
    config = TOOL["_case_config"](
        specification,
        time_steps=4,
        turbulent_diffusivity_model=(
            "molecular_plus_subgrid_eddy_viscosity"
        ),
    )

    assert config.scalar_turbulent_diffusivity_model == (
        "molecular_plus_subgrid_eddy_viscosity"
    )
    assert config.thermal_turbulent_diffusivity_model == (
        "molecular_plus_subgrid_eddy_viscosity"
    )


def test_case_config_can_use_first_measured_mean_fan_plane():
    specification = TOOL["_case_specification"]("published_box")
    config = TOOL["_case_config"](
        specification,
        time_steps=1,
        fan_inlet_profile="measured_first_station",
    )
    _kwargs, evidence = TOOL["_fan_inlet_profile"](
        specification, "measured_first_station",
    )

    assert config.active_inlet_profile_y_m == pytest.approx(
        (-0.05, 0.20, 0.45, 0.70, 0.95)
    )
    assert config.active_inlet_profile_u_m_s == pytest.approx(
        (3.66247154471545, 2.90246341463415, 3.32214634146341,
         4.24052032520325, 3.73039024390244)
    )
    assert evidence["measured_x_m"] == pytest.approx(0.25)
    assert evidence["workbook_sha256"] == (
        "bdb623e84b3cdc9d67884703762051d40c42e02bbe294df3de97ea062c45cc8a"
    )
    assert evidence["turbulence_statistics_available"] is False


def test_boundary_clearance_rejects_open_face_contact():
    field = np.zeros((8, 8, 12), dtype=float)
    field[4, 4, 5] = 5.0
    clear = TOOL["_boundary_clearance"](field)
    assert clear["boundary_clearance_gate_passed"] is True
    field[4, 0, 5] = 0.11
    contacted = TOOL["_boundary_clearance"](field)
    assert contacted["open_face_concentration_gate_passed"] is False
    assert contacted["boundary_clearance_gate_passed"] is False


def test_boundary_clearance_allows_intentional_x_advective_outflow():
    field = np.zeros((10, 10, 12), dtype=float)
    field[5, 5, -1] = 5.0
    result = TOOL["_boundary_clearance"](field)
    assert result["open_face_concentration_gate_passed"] is True
    assert result["four_vol_boundary_clearance_gate_passed"] is True
    assert result["boundary_clearance_gate_passed"] is True
    assert result["minimum_4vol_open_boundary_clearance_cells_by_face"][
        "x_max"
    ] == 0


def test_runtime_assessment_matches_expanded_case(monkeypatch):
    published = _case()
    published.update({"name": "published_box", "runtime_s": 20.0})
    expanded = _case()
    expanded.update({"name": "expanded_box", "runtime_s": 299.0})
    for case in (published, expanded):
        case.update({
            "common_final_field": [[[0.0, 0.0]]],
            "common_maximum_centre_slice": [[1.0, 4.0]],
        })
    monkeypatch.setitem(
        TOOL["run_audit"].__globals__,
        "_run_case",
        lambda name, *args, **kwargs: {
            "published_box": published, "expanded_box": expanded,
        }[name].copy(),
    )
    result = TOOL["run_audit"](workers=1)
    assert result["runtime_evidence"]["expanded_box_target_passed"] is True
    assert "meets_5min" in result["runtime_evidence"]["assessment"]


def test_receptor_tail_stationarity_accepts_stable_endpoint():
    case = {
        "receptor_traces_vol_percent": {
            "centre": [0.0] * 8 + [4.000, 4.001],
            "side": [0.0] * 8 + [0.1000, 0.1005],
        },
    }
    result = TOOL["_receptor_tail_stationarity"](case)
    assert result["passed"] is True
    assert result["maximum_relative_peak_to_peak"] == pytest.approx(
        0.0005 / 0.10025
    )


def test_receptor_tail_stationarity_rejects_drifting_endpoint():
    case = {
        "receptor_traces_vol_percent": {
            "centre": [0.0] * 8 + [4.0, 4.2],
        },
    }
    result = TOOL["_receptor_tail_stationarity"](case)
    assert result["passed"] is False
