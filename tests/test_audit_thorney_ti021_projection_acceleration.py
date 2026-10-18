from tools.audit_thorney_ti021_projection_acceleration import compare_records


def _record(solver: str, runtime: float, iterations: int):
    return {
        "configuration": {
            "cell_size_m": 2.0,
            "domain_m": [130.0, 140.0, 16.0],
            "cells": [65, 70, 8],
            "cell_count": 36_400,
            "active_steps": 1,
            "active_time_step_s": 0.002,
            "active_duration_s": 0.002,
            "wind_speed_m_s": 3.9,
            "streamwise_wind_m_s": 3.88,
            "crosswind_m_s": 0.414,
            "carrier_flow_initialization": "projected_inlet",
            "turbulence_closure": "none_preflight_only",
        },
        "runtime": {"total_s": runtime},
        "all_numerical_gates_pass": True,
        "numerics": {
            "projection_solver_used": solver,
            "projection_maximum_linear_iterations_per_step": iterations,
            "maximum_final_species_mole_fraction": 0.3200003,
            "maximum_relative_density_departure": 1.0158,
            "maximum_continuity_residual_kg_m3_s": 1.0e-11,
            "maximum_relative_species_mass_residual": 1.0e-16,
            "projection_preconditioner": (
                "symmetric_smoothed_aggregation_v_cycle"
                if solver == "amg_cg" else "jacobi_diagonal"
            ),
            "projection_preconditioner_level_count": (
                4 if solver == "amg_cg" else None
            ),
            "projection_preconditioner_operator_complexity": (
                1.5 if solver == "amg_cg" else None
            ),
            "pyamg_version": "5.3.0" if solver == "amg_cg" else None,
        },
    }


def test_projection_acceleration_requires_equivalence_and_material_speedup():
    result = compare_records(
        _record("cg", 4.5, 4_300),
        _record("amg_cg", 1.5, 110),
    )

    assert result["status"] == "accepted_projection_acceleration_component"
    assert result["all_gates_pass"]
    assert result["comparison"]["runtime_speedup"] == 3.0
    assert result["comparison"]["projection_iteration_reduction_factor"] > 30
