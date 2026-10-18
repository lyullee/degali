"""Audit the public HyTunnel-CS Pro-Science spatial H2 dataset.

The Zenodo record supplies concentration matrices, ventilation velocity
matrices, metadata and CAD geometry.  This audit verifies the published file
digests, inventories the matrices and exposes material metadata gaps.  It
does not promote the room-temperature GH2 tunnel data to LH2 field
validation and it does not invent a transient clock or sensor height.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile

from degali.addons.high_pressure_source_handoff_3d import (
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle
from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    HyTunnelConcentrationPoint,
    HyTunnelFlowPoint,
    read_hytunnel_concentration_workbook,
    read_hytunnel_flow_workbook,
)


SCHEMA = "degali.hytunnel-ps-spatial-audit.v23"
EXPECTED_MD5 = {
    "HTE2440PS000MIXED_300321.pdf": "d4b59189dd4017c695cc652f4f572128",
    "HTE2440PS001MIXED_300321.xlsx": "da176772755b1a2449a76df9e12b9cd3",
    "HTE2440PS002MIXED_300321.xlsx": "20f996be57fc3d8cf7d091d1f982bd75",
    "HTE2440PS003MIXED_300321.xlsx": "953b6a033553288831b1ee02b881dbec",
    "HTE2440PS004MIXED_300321.xlsx": "093d3d8e6ef030615704ee9096d3d84f",
    "HTE2440PS005FLMT00300321.xlsx": "c4b9195a2bb89213eda9aa2ac9e9ec07",
    "HTE2440PS006MIXED_300321.zip": "ade3911c0df50ee70486153455f4bd12",
}
CONCENTRATION_FILES = tuple(
    f"HTE2440PS00{number}MIXED_300321.xlsx" for number in range(1, 5)
)
FLOW_FILE = "HTE2440PS005FLMT00300321.xlsx"
CAD_FILE = "HTE2440PS006MIXED_300321.zip"
DECLARED_CONFIGURATION_COUNT = 40
ICHS_EXPERIMENT_PAPER = "https://hysafe.info/uploads/papers/2021/41.pdf"
HYTUNNEL_FINAL_REPORT = (
    "https://hytunnel.net/wordpress/wp-content/uploads/2022/09/"
    "HyTunnel-CS_D2.3-Final-report-on-dispersion-in-tunnels.pdf"
)


def _matched_source_adapter_audit() -> dict[str, object]:
    """Reproduce the declared DEGALI source route for the matched CFD case."""
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=286.15,
        storage_pressure=11.515e6,
        mass_flow=0.005,
        orifice_diameter=0.001,
        ambient_pressure=98_500.0,
        allow_supercritical_gas=True,
    )
    handoff = low_mach_handoff_from_energy_nozzle(
        nozzle,
        source_id="hytunnel-1mm-5g-3p5ms-coflow",
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_axial_velocity_m_s=3.5,
        maximum_mach_number=0.30,
    )
    fields = handoff.resolved_source_fields(
        nozzle_position_m=(0.0, 0.45, 0.45),
        grid_spacing_m=0.0064,
        duration_s=0.707,
    )
    assert nozzle.diameter is not None
    assert nozzle.temperature is not None
    assert nozzle.density is not None
    return {
        "status": (
            "implemented_profile_preflight_nominal_failed_0p06_candidate_only"
            if handoff.numerically_admissible and fields.spatially_resolved
            else "withheld"
        ),
        "input_basis": "matched conditions transcribed from the HyTunnel public report",
        "storage_temperature_k": 286.15,
        "storage_pressure_pa_abs": 11.515e6,
        "ambient_temperature_k": 286.15,
        "ambient_pressure_pa_abs": 98_500.0,
        "orifice_diameter_m": 0.001,
        "hydrogen_mass_flow_kg_s": 0.005,
        "ambient_coflow_velocity_m_s": 3.5,
        "degali_notional_nozzle": {
            "model": handoff.source.source_model,
            "diameter_m": nozzle.diameter,
            "velocity_m_s": nozzle.velocity,
            "temperature_k": nozzle.temperature,
            "density_kg_m3": nozzle.density,
            "momentum_flux_n": nozzle.momentum_flux,
            "discharge_coefficient": nozzle.discharge_coefficient,
            "relative_mass_residual": nozzle.relative_mass_residual,
            "relative_momentum_residual": nozzle.relative_momentum_residual,
            "relative_energy_residual": nozzle.relative_energy_residual,
        },
        "official_cfd_notional_nozzle": {
            "model": "Molkov et al. (2012), as declared by HyTunnel D2.3",
            "diameter_m": 0.00725,
            "complete_mass_momentum_energy_plane_available_in_public_report": False,
        },
        "notional_diameter_ratio_degali_to_official": nozzle.diameter / 0.00725,
        "source_model_match_to_official_cfd": False,
        "source_model_difference_policy": (
            "report both routes; do not tune the DEGALI nozzle diameter to the "
            "published concentration contours"
        ),
        "low_mach_handoff": {
            "maximum_mach_number": handoff.maximum_mach_number,
            "mach_number": handoff.mach_number,
            "maximum_relative_density_departure": (
                handoff.maximum_relative_density_departure
            ),
            "relative_density_departure": handoff.relative_density_departure,
            "distance_from_orifice_m": handoff.handoff_distance_m,
            "entrained_air_mass_flow_kg_s": handoff.entrained_air_mass_flow_kg_s,
            "total_mass_flow_kg_s": handoff.total_mass_flow_kg_s,
            "hydrogen_mass_fraction": handoff.hydrogen_mass_fraction,
            "mixture_velocity_m_s": handoff.mixture_velocity_m_s,
            "mixture_temperature_k": handoff.mixture_temperature_k,
            "mixture_density_kg_m3": handoff.mixture_density_kg_m3,
            "mixture_sound_speed_m_s": handoff.mixture_sound_speed_m_s,
            "equivalent_diameter_m": handoff.equivalent_diameter_m,
            "moment_matched_radial_sigma_m": (
                handoff.moment_matched_radial_sigma_m
            ),
            "plane_momentum_flux_n": handoff.plane_momentum_flux_n,
            "source_increment_momentum_flux_n": (
                handoff.source_increment_momentum_flux_n
            ),
            "plane_sensible_enthalpy_flux_w": (
                handoff.plane_sensible_enthalpy_flux_w
            ),
            "source_increment_total_energy_flux_w": (
                handoff.source_increment_total_energy_flux_w
            ),
            "equivalent_hydrogen_temperature_k": (
                handoff.equivalent_hydrogen_temperature_k
            ),
            "maximum_relative_conservation_residual": (
                handoff.maximum_relative_conservation_residual
            ),
        },
        "resolved_source_mapping": {
            "position_m": list(fields.source_position_m),
            "sigma_xyz_m": list(fields.source_sigma_xyz_m),
            "source_specific_momentum_vector_m_s": list(
                fields.source_momentum_velocity_m_s
            ),
            "equivalent_hydrogen_temperature_k": fields.source_temperature_k,
            "official_grid_cell_size_m": 0.0064,
            "cells_across_equivalent_diameter": (
                fields.cells_across_equivalent_diameter
            ),
            "minimum_cells_across_equivalent_diameter": 4.0,
            "spatially_resolved": fields.spatially_resolved,
            "top_hat_handoff_numerically_admissible": (
                handoff.numerically_admissible
            ),
            "mapped_gaussian_profile_density_admissibility": (
                "nominal_failed_0p06_candidate_only_not_promoted"
            ),
            "preflight_artifact": (
                "outputs/hytunnel-resolved-source-preflight-2026-10-10-v6.json"
            ),
            "source_mass_over_official_steady_state_kg": (
                fields.source_schedule.released_mass_kg
            ),
        },
        "limitations": [
            "the entrainment distance is a source-momentum closure, not a measured station",
            "the receiving low-Mach solver does not yet close a resolved total-energy equation",
            "the official CFD source diameter cannot be reconstructed as a complete flux plane from diameter alone",
            "top-hat mean-plane admissibility does not guarantee Gaussian mapped-profile admissibility",
            "the raw second-order-extrapolated angle-aware candidate passes the declared 8-cell and 10-cell linear-acoustic 3-D matrices but becomes unstable in the finite-vortex audit; a separate TVD-limited acoustic/vortical candidate remains stable and preserves vorticity but fails the frozen velocity and density-contamination gates and its 10-cell acoustic matrix, so no single pressure-open candidate is qualified; the official 0.707 s same-grid domain audit passes but does not waive this component failure",
        ],
    }


def _condition_record(key: tuple[float, float, float, str, str]) -> dict[str, object]:
    nozzle, mass_flow, wind, direction, gas = key
    return {
        "nozzle_diameter_mm": nozzle,
        "nominal_mass_flow_g_s": mass_flow,
        "nominal_ventilation_speed_m_s": wind,
        "flow_direction": direction,
        "gas_label": gas,
    }


def _expected_test_matrix() -> set[tuple[float, float, float, str, str]]:
    expected: set[tuple[float, float, float, str, str]] = set()
    speeds = (1.5, 3.5, 5.0)
    for nozzle in (1.0, 4.0):
        expected.add((nozzle, 1.0, 0.0, "quiescent", "H2"))
        for direction in ("co-flow", "counter-flow", "cross-flow"):
            for speed in speeds:
                expected.add((nozzle, 1.0, speed, direction, "H2"))
        expected.add((nozzle, 5.0, 0.0, "quiescent", "H2"))
        for direction in ("co-flow", "counter-flow"):
            for speed in speeds:
                expected.add((nozzle, 5.0, speed, direction, "H2"))
    for nozzle, mass_flow in ((1.0, 1.5), (4.0, 2.5)):
        for speed in speeds:
            expected.add((nozzle, mass_flow, speed, "cross-flow", "H2"))
    assert len(expected) == DECLARED_CONFIGURATION_COUNT
    return expected


def _test_matrix_reconciliation(configurations) -> dict[str, object]:
    expected = _expected_test_matrix()
    observed = {
        (
            float(configuration.nozzle_diameter_mm),
            float(configuration.nominal_mass_flow_g_s),
            float(configuration.nominal_ventilation_speed_m_s),
            configuration.flow_direction,
            configuration.gas_label,
        )
        for configuration in configurations
        if configuration.nozzle_diameter_mm is not None
        and configuration.nominal_mass_flow_g_s is not None
    }
    matched = expected & observed
    missing = expected - observed
    additional_hydrogen = {
        key for key in observed - expected if key[-1] == "H2"
    }
    non_hydrogen = {key for key in observed if key[-1] != "H2"}
    order = lambda key: (key[0], key[1], key[3], key[2], key[4])
    return {
        "declared_hydrogen_case_count": len(expected),
        "matched_declared_hydrogen_case_count": len(matched),
        "missing_declared_hydrogen_cases": [
            _condition_record(key) for key in sorted(missing, key=order)
        ],
        "additional_hydrogen_baselines": [
            _condition_record(key)
            for key in sorted(additional_hydrogen, key=order)
        ],
        "non_hydrogen_cases": [
            _condition_record(key) for key in sorted(non_hydrogen, key=order)
        ],
        "reconciled": bool(
            len(matched) == 39
            and len(missing) == 1
            and len(additional_hydrogen) == 2
            and len(non_hydrogen) == 1
        ),
    }


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _range(values: Iterable[float]) -> list[float] | None:
    items = [float(value) for value in values]
    if not items:
        return None
    if not all(math.isfinite(value) for value in items):
        raise ValueError("HyTunnel audit contains a non-finite value")
    return [min(items), max(items)]


def _point_summary(
    points: Iterable[HyTunnelConcentrationPoint | HyTunnelFlowPoint],
    value_name: str,
) -> dict[str, object]:
    items = list(points)
    return {
        "point_count": len(items),
        "statistic_counts": dict(sorted(Counter(point.statistic for point in items).items())),
        "x_range_m": _range(point.x_m for point in items),
        "y_range_m": _range(point.y_m for point in items),
        f"{value_name}_range": _range(getattr(point, value_name) for point in items),
    }


def _cad_inventory(path: Path) -> dict[str, object]:
    with ZipFile(path) as archive:
        names = archive.namelist()
        unsafe = [
            name for name in names
            if Path(name).is_absolute() or ".." in Path(name).parts
        ]
        suffixes = Counter(Path(name).suffix.lower() for name in names if name)
    return {
        "archive_entry_count": len(names),
        "unsafe_path_count": len(unsafe),
        "inventor_part_count": suffixes.get(".ipt", 0),
        "inventor_assembly_count": suffixes.get(".iam", 0),
        "suffix_counts": dict(sorted((key or "<none>", value) for key, value in suffixes.items())),
        "geometry_extracted": False,
    }


def audit_hytunnel_ps(root: str | Path) -> dict[str, Any]:
    directory = Path(root)
    if not directory.is_dir():
        return {
            "schema": SCHEMA,
            "status": "not_available",
            "dataset_doi": DATASET_DOI,
            "root": str(directory.resolve()),
            "promotion_allowed": False,
            "gate_codes": ["dataset_not_available", "not_lh2_field_validation"],
        }

    files: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    for name, expected_md5 in EXPECTED_MD5.items():
        path = directory / name
        try:
            if not path.is_file():
                raise FileNotFoundError(path)
            md5 = _digest(path, "md5")
            if md5 != expected_md5:
                raise ValueError(
                    f"published MD5 mismatch: expected {expected_md5}, got {md5}"
                )
            files.append({
                "file": name,
                "bytes": path.stat().st_size,
                "md5": md5,
                "published_md5_verified": True,
                "sha256": _digest(path, "sha256"),
            })
        except (OSError, ValueError) as error:
            errors.append({"file": name, "error": str(error)})

    configurations = []
    concentration_points: list[HyTunnelConcentrationPoint] = []
    workbook_summaries: list[dict[str, object]] = []
    for name in CONCENTRATION_FILES:
        try:
            dataset = read_hytunnel_concentration_workbook(directory / name)
            configurations.extend(dataset.configurations)
            concentration_points.extend(dataset.points)
            workbook_summaries.append({
                "file": name,
                "configuration_count": len(dataset.configurations),
                **_point_summary(dataset.points, "concentration_vol_percent"),
            })
        except (OSError, RuntimeError, TypeError, ValueError) as error:
            errors.append({"file": name, "error": str(error)})

    flow_points: list[HyTunnelFlowPoint] = []
    flow_configuration_count = 0
    try:
        flow = read_hytunnel_flow_workbook(directory / FLOW_FILE)
        flow_configuration_count = len(flow.configuration_ids)
        flow_points.extend(flow.points)
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        errors.append({"file": FLOW_FILE, "error": str(error)})

    cad = None
    try:
        cad = _cad_inventory(directory / CAD_FILE)
    except (OSError, BadZipFile) as error:
        errors.append({"file": CAD_FILE, "error": str(error)})

    actual_count = len(configurations)
    hydrogen_count = sum(configuration.is_direct_hydrogen for configuration in configurations)
    count_discrepancy = actual_count != DECLARED_CONFIGURATION_COUNT
    matrix_reconciliation = _test_matrix_reconciliation(configurations)
    all_hashes_verified = len(files) == len(EXPECTED_MD5) and not any(
        "MD5 mismatch" in item["error"] for item in errors
    )
    parsed = bool(concentration_points and flow_points and cad is not None)
    status = (
        "component_evidence_conditional"
        if all_hashes_verified and parsed and not errors
        else "withheld"
    )
    source_adapter = _matched_source_adapter_audit()
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "dataset_doi": DATASET_DOI,
        "root": str(directory.resolve()),
        "license": "CC BY 4.0",
        "files": files,
        "all_published_md5_verified": all_hashes_verified,
        "errors": errors,
        "concentration": {
            "declared_configuration_count": DECLARED_CONFIGURATION_COUNT,
            "workbook_sheet_count": actual_count,
            "hydrogen_labeled_sheet_count": hydrogen_count,
            "helium_labeled_sheet_count": actual_count - hydrogen_count,
            "metadata_configuration_count_discrepancy": count_discrepancy,
            "test_matrix_reconciliation": matrix_reconciliation,
            "nozzle_diameter_counts_mm": dict(sorted(Counter(
                str(configuration.nozzle_diameter_mm) for configuration in configurations
            ).items())),
            "nominal_mass_flow_counts_g_s": dict(sorted(Counter(
                str(configuration.nominal_mass_flow_g_s) for configuration in configurations
            ).items())),
            "nominal_ventilation_speed_counts_m_s": dict(sorted(Counter(
                str(configuration.nominal_ventilation_speed_m_s)
                for configuration in configurations
            ).items())),
            "flow_direction_counts": dict(sorted(Counter(
                configuration.flow_direction for configuration in configurations
            ).items())),
            "workbooks": workbook_summaries,
            **_point_summary(concentration_points, "concentration_vol_percent"),
        },
        "ventilation_flow": {
            "configuration_count": flow_configuration_count,
            **_point_summary(flow_points, "velocity_x_m_s"),
        },
        "cad": cad,
        "scope": "gh2_ventilated_jet_spatial_component_only",
        "measurement_protocol": {
            "source": ICHS_EXPERIMENT_PAPER,
            "concentration_sensor_count": 9,
            "sampling_frequency_hz": 2.0,
            "five_g_s_steady_average_window_s": 10.0,
            "one_g_s_steady_average_window_s": 25.0,
            "maximum_axial_sensor_sleigh_distance_m": 3.5,
            "measurement_plane": "horizontal x-y plane through the jet axis",
            "raw_transient_series_in_public_workbooks": False,
        },
        "matched_reference_case": {
            "case_id": "hytunnel-1mm-5g-3p5ms-coflow",
            "experiment_workbook": "HTE2440PS002MIXED_300321.xlsx",
            "experiment_sheet": "1 mm 5g H2 W = 3.5 co-flow",
            "nozzle_diameter_mm": 1.0,
            "nominal_hydrogen_mass_flow_g_s": 5.0,
            "nominal_ventilation_speed_m_s": 3.5,
            "flow_direction": "co-flow",
            "official_cfd_source": HYTUNNEL_FINAL_REPORT,
            "official_cfd_domain_m": [2.5, 0.9, 0.9],
            "official_cfd_source_axis_m": [0.0, 0.45, 0.45],
            "official_cfd_cell_size_m": 0.0064,
            "official_cfd_cell_count": 7_644_000,
            "official_cfd_notional_nozzle_diameter_mm": 7.25,
            "official_cfd_turbulence_model": "RNG k-epsilon",
            "official_cfd_boundary_configuration": {
                "x_zero": "3.5 m/s prescribed air velocity inlet",
                "remaining_five_faces": "open non-reflecting",
                "solid_tunnel_walls_in_computational_box": False,
            },
            "official_cfd_steady_state_time_s": 0.707,
            "published_experiment_4volpct_length_m": 1.65,
            "published_cfd_4volpct_length_m": 2.25,
            "published_experiment_10volpct_length_m": 0.60,
            "published_cfd_10volpct_length_m": 0.50,
            "degali_source_adapter": source_adapter,
            "degali_boundary_configuration": {
                "status": (
                    "official_duration_domain_passed_but_unified_vortical_"
                    "boundary_and_runtime_gates_failed"
                ),
                "flow_boundary_scheme": (
                    "characteristic_multidimensional_extrapolated"
                ),
                "inlet_boundary": "prescribed co-flow velocity with adjacent-plane non-equilibrium stress",
                "outlet_boundary": "fixed ambient pressure with adjacent-plane non-equilibrium stress",
                "nonstreamwise_wall_boundary": "all_pressure_open",
                "nonstreamwise_boundary_interpretation": (
                    "second-order-extrapolated angle-aware isothermal characteristic incoming-wave suppression with extrapolated non-equilibrium stress"
                ),
                "hydrogen_scalar_boundary": (
                    "outflow-only ambient-zero advection with zero normal diffusive gradient"
                ),
                "external_domain_sensitivity_required": False,
                "short_horizon_plume_domain_sensitivity_passed": True,
                "short_horizon_plume_duration_s": (
                    0.12034285714285714
                ),
                "short_horizon_maximum_primary_relative_change": (
                    0.038646077041392043
                ),
                "short_horizon_common_field_normalised_l1": (
                    0.010921300571742875
                ),
                "short_horizon_expanded_runtime_s": 220.395580100012,
                "short_horizon_runtime_target_s": 300.0,
                "short_horizon_runtime_target_passed": True,
                "short_horizon_projection_linear_solver": "cached_direct",
                "official_0p707s_domain_sensitivity_required": False,
                "official_0p707s_domain_sensitivity_passed": True,
                "official_0p707s_duration_s": 0.7072457142857143,
                "official_0p707s_time_steps": 3056,
                "official_0p707s_maximum_primary_relative_change": (
                    0.009900990099009877
                ),
                "official_0p707s_common_field_normalised_l1": (
                    0.012446618980418217
                ),
                "official_0p707s_conservation_gate_passed": True,
                "official_0p707s_active_scalar_applicability_gate_passed": True,
                "official_0p707s_expanded_boundary_clearance_passed": True,
                "official_0p707s_published_boundary_clearance_passed": False,
                "official_0p707s_expanded_lateral_face_maximum_vol_percent": (
                    0.09474100565342727
                ),
                "official_0p707s_expanded_4vol_lateral_clearance_cells": 11,
                "official_0p707s_published_runtime_s": 372.10203780001,
                "official_0p707s_expanded_runtime_s": 1402.1545282999286,
                "official_0p707s_runtime_target_s": 300.0,
                "official_0p707s_runtime_target_passed": False,
                "oblique_angle_wavelength_sweep_required": False,
                "oblique_angle_wavelength_sweep_gate_passed": True,
                "corner_interaction_validation_required": False,
                "canonical_xy_corner_reflection_gate_passed": True,
                "canonical_xy_corner_reflection_ratio": (
                    0.02594083208986574
                ),
                "three_dimensional_corner_orientation_validation_required": False,
                "extrapolated_wavelength8_3d_corner_gate_passed": True,
                "extrapolated_wavelength10_3d_corner_gate_passed": True,
                "multiwavelength_3d_corner_validation_required": False,
                "extrapolated_wavelength8_scaled_xy_control_reflection_ratio": 0.024176960300331944,
                "extrapolated_wavelength8_xz_edge_reflection_ratio": 0.02438657311899492,
                "extrapolated_wavelength8_yz_edge_reflection_ratio": 0.024198275880320055,
                "extrapolated_wavelength8_xyz_corner_reflection_ratio": 0.0184597483357786,
                "extrapolated_wavelength10_scaled_xy_control_reflection_ratio": 0.031921774015279666,
                "extrapolated_wavelength10_xz_edge_reflection_ratio": 0.03178734824695349,
                "extrapolated_wavelength10_yz_edge_reflection_ratio": 0.03190833330649736,
                "extrapolated_wavelength10_xyz_corner_reflection_ratio": 0.03151622268551212,
                "predecessor_boundary_scheme": "characteristic_multidimensional",
                "predecessor_wavelength8_3d_corner_gate_passed": True,
                "predecessor_wavelength10_3d_corner_gate_passed": False,
                "predecessor_wavelength8_xyz_corner_reflection_ratio": 0.0364427770370331,
                "predecessor_wavelength10_xyz_corner_reflection_ratio": 0.05998672094052834,
                "vortical_outflow_boundary_validation_required": True,
                "vortical_limited_boundary_scheme": (
                    "characteristic_multidimensional_vortical_limited"
                ),
                "vortical_limited_candidate_stable": True,
                "vortical_limited_candidate_gate_passed": False,
                "vortical_limited_maximum_normalised_velocity_l2": (
                    0.1982484210333665
                ),
                "vortical_limited_maximum_normalised_vorticity_l2": (
                    0.009429916268370354
                ),
                "vortical_limited_maximum_normalised_density_l2": (
                    2.4322976513928833
                ),
                "vortical_limited_residual_enstrophy_ratio": (
                    8.944661359127176e-05
                ),
                "vortical_limited_final_mean_density_error": (
                    7.5693809974408666e-06
                ),
                "vortical_limited_normalised_error_limit": 0.05,
                "vortical_limited_wavelength8_3d_corner_gate_passed": True,
                "vortical_limited_wavelength10_3d_corner_gate_passed": False,
                "vortical_limited_wavelength8_worst_reflection_ratio": (
                    0.04214847153253591
                ),
                "vortical_limited_wavelength10_worst_reflection_ratio": (
                    0.05564303573851233
                ),
                "unified_pressure_open_boundary_gate_passed": False,
                "streamwise_linear_acoustic_reflection_ratio": (
                    0.01638783115949678
                ),
                "streamwise_linear_acoustic_reflection_gate_passed": True,
                "lateral_normal_linear_acoustic_reflection_ratio": (
                    0.016387831157066126
                ),
                "lateral_normal_linear_acoustic_reflection_gate_passed": True,
                "oblique_30deg_linear_acoustic_reflection_ratio": (
                    0.0219289505449804
                ),
                "oblique_30deg_linear_acoustic_reflection_gate_passed": True,
                "reflection_audit_artifact": (
                    "outputs/pressure-open-boundary-reflection-audit-"
                    "2026-10-10-v2.json"
                ),
                "oblique_reflection_audit_artifact": (
                    "outputs/pressure-open-multidimensional-sweep-"
                    "2026-10-10.json"
                ),
                "angle_wavelength_worst_reflection_ratio": (
                    0.03906875314425576
                ),
                "corner_reflection_audit_artifact": (
                    "outputs/pressure-open-corner-reflection-audit-"
                    "2026-10-10.json"
                ),
                "three_dimensional_corner_reflection_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10.json"
                ),
                "three_dimensional_corner_wavelength10_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10-wavelength10.json"
                ),
                "three_dimensional_corner_extrapolated_wavelength8_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10-extrapolated-wavelength8.json"
                ),
                "three_dimensional_corner_extrapolated_wavelength10_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10-extrapolated-wavelength10.json"
                ),
                "plume_domain_sensitivity_artifact": (
                    "outputs/hytunnel-plume-domain-sensitivity-"
                    "2026-10-11-official-0p707-v6.json"
                ),
                "vortical_outflow_audit_artifact": (
                    "outputs/pressure-open-vortical-outflow-audit-"
                    "2026-10-10.json"
                ),
                "vortical_limited_wavelength8_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10-mode-aware-wavelength8-v2.json"
                ),
                "vortical_limited_wavelength10_audit_artifact": (
                    "outputs/pressure-open-3d-corner-reflection-audit-"
                    "2026-10-10-mode-aware-wavelength10-v2.json"
                ),
            },
            "degali_score_status": (
                "withheld_resolved_observation_operator_and_accuracy_gates_failed"
            ),
            "degali_resolved_overlap_score": {
                "artifact": (
                    "outputs/hytunnel-resolved-transient-score-2026-10-11-v2.json"
                ),
                "endpoint_replay_artifact": (
                    "outputs/hytunnel-resolved-published-box-0p707-final-"
                    "2026-10-11-v1.json"
                ),
                "status": "resolved_transient_validation_withheld",
                "published_average_point_count": 153,
                "supported_downstream_point_count": 54,
                "excluded_upstream_of_handoff_count": 81,
                "excluded_beyond_common_domain_count": 18,
                "model_operator": (
                    "final_snapshot_horizontal_source_axis_slice"
                ),
                "measurement_operator": "published_steady_window_average",
                "operator_match": False,
                "qualified_steady_endpoint_proxy": True,
                "maximum_receptor_tail_relative_peak_to_peak": (
                    6.865646387495322e-05
                ),
                "maximum_slice_replay_normalised_l1": 0.0,
                "mae_vol_pct_point": 2.17922030234363,
                "rmse_vol_pct_point": 2.49050143281193,
                "mean_bias_vol_pct_point": 2.16970774137294,
                "pearson_r": 0.733438352676071,
                "fac2_floor_0p10": 0.5555555555555556,
                "four_vol_sensitivity": 1.0,
                "four_vol_specificity": 0.5217391304347826,
                "interpretation": (
                    "stationary endpoint proxy versus published steady average; "
                    "accuracy gates fail and the unsupported 99 points are not "
                    "extrapolated"
                ),
            },
            "degali_source_mixing_sensitivity": {
                "artifact": (
                    "outputs/hytunnel-source-mixing-sensitivity-"
                    "2026-10-11-v1.json"
                ),
                "status": "diagnostic_source_mixing_sensitivity_not_calibration",
                "common_supported_point_count": 45,
                "baseline_handoff0p06_sc0p7_mae_vol_pct_point": (
                    2.1130664701948967
                ),
                "handoff0p06_sc0p5_mae_vol_pct_point": 1.6697072118487113,
                "handoff0p06_sc1p0_mae_vol_pct_point": 2.5290350563213404,
                "handoff0p05_sc0p7_mae_vol_pct_point": 2.185452923053046,
                "handoff0p06_sc0p5_fac2": 0.6888888888888889,
                "handoff0p06_sc0p5_four_vol_specificity": 0.675,
                "calibration_or_best_variant_selection_allowed": False,
                "interpretation": (
                    "the error is more sensitive to the unvalidated effective "
                    "scalar-mixing closure than to tightening the admissible "
                    "handoff from 6% to 5%; no variant is selected from observations"
                ),
            },
            "degali_local_sgs_closure_candidate": {
                "case_artifact": (
                    "outputs/hytunnel-resolved-sgs-diffusivity-handoff0p05-"
                    "official-0p707-2026-10-11-v1.json"
                ),
                "score_artifact": (
                    "outputs/hytunnel-resolved-sgs-diffusivity-handoff0p05-"
                    "score-2026-10-11-v1.json"
                ),
                "comparison_artifact": (
                    "outputs/hytunnel-scalar-closure-comparison-"
                    "2026-10-11-v1.json"
                ),
                "status": "local_sgs_candidate_rejected_by_frozen_accuracy_gates",
                "common_supported_point_count": 45,
                "handoff_density_departure_limit": 0.05,
                "effective_schmidt_number": 0.7,
                "scalar_diffusivity_model": (
                    "molecular_plus_subgrid_eddy_viscosity"
                ),
                "maximum_relative_density_departure": 0.0925393129270758,
                "maximum_relative_mass_residual": 2.37920391415447e-14,
                "maximum_relative_energy_residual": 8.94343577421222e-15,
                "runtime_s": 366.1471571000293,
                "mae_vol_pct_point": 3.4248066947066316,
                "rmse_vol_pct_point": 3.9756890344428157,
                "mean_bias_vol_pct_point": 3.3804055972032074,
                "pearson_r": 0.3964719367969042,
                "fac2_floor_0p10": 0.37777777777777777,
                "four_vol_sensitivity": 1.0,
                "four_vol_specificity": 0.25,
                "mae_increase_vs_same_handoff_legacy_vol_pct_point": (
                    1.239353771653584
                ),
                "calibration_or_best_variant_selection_allowed": False,
                "interpretation": (
                    "local WALE alone under-mixes the coarse steady carrier; "
                    "the legacy uniform diffusion is compensating numerical "
                    "mixing and is not promoted as validated physics"
                ),
            },
        },
        "promotion_allowed": False,
        "quantitative_lh2_field_validation_allowed": False,
        "usable_comparisons": [
            "published average concentration at matched x-y coordinates",
            "published maximum and minimum concentration envelopes",
            "measured ventilation x-velocity matrices",
            "CAD-controlled obstacle geometry after independent conversion verification",
        ],
        "gate_codes": [
            "component_evidence_only",
            "not_lh2_field_validation",
            "no_transient_common_clock",
            "physical_measurement_elevation_not_independently_declared",
            "sensor_calibration_certificate_not_in_record",
            "one_declared_hydrogen_crossflow_case_is_helium",
            "two_additional_stagnant_hydrogen_baselines",
            "cad_geometry_not_yet_extracted",
            "matched_high_pressure_source_adapter_external_validation_required",
            "nominal_resolved_source_profile_preflight_failed",
            "matched_source_model_differs_from_official_cfd",
            "matched_pressure_open_unified_boundary_candidate_failed",
            "matched_official_duration_runtime_target_failed",
            "matched_resolved_concentration_score_failed",
            "matched_resolved_endpoint_proxy_only",
            "matched_scalar_mixing_external_basis_required",
            "local_sgs_scalar_closure_accuracy_failed",
        ],
    }


def write_audit(record: dict[str, Any], output: str | Path) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    record = audit_hytunnel_ps(args.root)
    if args.output is None:
        print(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        write_audit(record, args.output)
        print(
            f"wrote {args.output} status={record['status']} "
            f"sheets={record.get('concentration', {}).get('workbook_sheet_count', 0)}"
        )
    return 0 if record["status"] in {
        "component_evidence_conditional", "not_available"
    } else 1


if __name__ == "__main__":
    raise SystemExit(main())
