"""Run one full-600 s TI021 wake-residence candidate against the frozen baseline."""

from __future__ import annotations

import argparse
import inspect
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from tools.run_ti021_measured_wind_candidate import _group_metrics
from tools.score_ti021_dense_layer_hybrid import run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-step-s", type=float, default=None)
    parser.add_argument("--capture-rate-s-inv", type=float, required=True)
    parser.add_argument("--return-rate-s-inv", type=float, required=True)
    parser.add_argument(
        "--reservoir-return-target", choices=("main", "reattached"),
        default=None,
    )
    parser.add_argument(
        "--reservoir-reattached-fraction", type=float, default=None,
        help=(
            "fraction of local-reservoir return routed to the conserved "
            "reattached carrier when --reservoir-return-target=reattached"
        ),
    )
    parser.add_argument(
        "--reservoir-capture-mode", choices=("cavity_rate", "wall_impact"),
        default=None,
    )
    parser.add_argument(
        "--reservoir-vertical-capture-mode",
        choices=("column", "below_wall_core"), default=None,
    )
    parser.add_argument(
        "--reservoir-spatial-mode", choices=("well_mixed", "local"),
        default=None,
    )
    parser.add_argument(
        "--reservoir-content-mode", choices=("conserved", "species_only"),
        default=None,
    )
    parser.add_argument(
        "--reservoir-transport-fraction", type=float, default=None,
    )
    parser.add_argument(
        "--reservoir-transport-height-m", type=float, default=None,
    )
    parser.add_argument(
        "--reservoir-transport-mask", choices=("cavity", "wake"), default=None,
    )
    parser.add_argument("--reservoir-wake-length-m", type=float, default=None)
    parser.add_argument(
        "--reservoir-wake-half-width-m", type=float, default=None,
    )
    parser.add_argument(
        "--reservoir-crosswind-convergence-rate-s-inv", type=float,
        default=None,
    )
    parser.add_argument(
        "--reservoir-crosswind-target-offset-m", type=float, default=None,
    )
    parser.add_argument(
        "--local-reservoir-observation-fraction", type=float, default=None,
    )
    parser.add_argument(
        "--local-reservoir-observation-center-m", type=float, default=None,
    )
    parser.add_argument(
        "--local-reservoir-observation-sigma-m", type=float, default=None,
    )
    parser.add_argument("--return-delay-s", type=float, default=None)
    parser.add_argument("--return-radius-m", type=float, default=None)
    parser.add_argument("--return-offset-m", type=float, default=None)
    parser.add_argument("--return-sigma-m", type=float, default=None)
    parser.add_argument("--shelter-strength", type=float, default=None)
    parser.add_argument("--entrainment-coefficient", type=float, default=None)
    parser.add_argument("--entrainment-far-coefficient", type=float, default=None)
    parser.add_argument("--entrainment-relaxation-time-s", type=float, default=None)
    parser.add_argument("--horizontal-diffusivity-m2-s", type=float, default=None)
    parser.add_argument("--crosswind-diffusivity-factor", type=float, default=None)
    parser.add_argument("--wake-meander-strength", type=float, default=None)
    parser.add_argument("--wake-shedding-strouhal", type=float, default=None)
    parser.add_argument("--wake-meander-phase-deg", type=float, default=None)
    parser.add_argument("--wake-length-m", type=float, default=None)
    parser.add_argument("--wake-base-width-m", type=float, default=None)
    parser.add_argument("--wake-spread-slope", type=float, default=None)
    parser.add_argument("--wind-angle-deg", type=float, default=None)
    parser.add_argument("--wind-history-direction-offset-deg", type=float, default=None)
    parser.add_argument("--wind-history-advection-reference-x-m", type=float, default=None)
    parser.add_argument("--wind-history-advection-speed-m-s", type=float, default=None)
    parser.add_argument("--layer-wind-speed-m-s", type=float, default=None)
    parser.add_argument("--friction-velocity-m-s", type=float, default=None)
    parser.add_argument("--gravity-scale", type=float, default=None)
    parser.add_argument("--source-crosswind-scale", type=float, default=None)
    parser.add_argument("--source-longitudinal-scale", type=float, default=None)
    parser.add_argument("--puff-gravity-spreading-coefficient", type=float, default=None)
    parser.add_argument("--puff-form-drag-coefficient", type=float, default=None)
    parser.add_argument(
        "--receptor-interpolation", choices=("nearest", "bilinear"),
        default=None,
    )
    parser.add_argument(
        "--vertical-observation-profile",
        choices=("top_hat", "dense_power", "exponential_tail", "diffusive_tail",
                 "handoff_gaussian", "multilayer"), default=None,
    )
    parser.add_argument("--vertical-layers", type=int, default=None)
    parser.add_argument("--vertical-mixing-diffusivity-m2-s", type=float, default=None)
    parser.add_argument("--multilayer-horizontal-shear-fraction", type=float, default=None)
    parser.add_argument(
        "--multilayer-obstacle-crosswind-shear-rate-s-inv", type=float,
        default=None,
    )
    parser.add_argument(
        "--multilayer-obstacle-crosswind-decay-length-m", type=float,
        default=None,
    )
    parser.add_argument(
        "--multilayer-relative-transport-interval-s", type=float, default=None,
    )
    parser.add_argument("--vertical-tail-diffusivity-m2-s", type=float, default=None)
    parser.add_argument("--initial-vertical-gaussian", action="store_true")
    parser.add_argument("--height-resolved-wall-flux", action="store_true")
    parser.add_argument(
        "--vertical-observation-component-diagnostics",
        action="store_true",
        help=(
            "store additive concentration traces for the ground, elevated, "
            "reattached, slow-wake, and local-reservoir observation branches"
        ),
    )
    parser.add_argument("--cavity-vertical-exponent", type=float, default=None)
    parser.add_argument(
        "--cavity-vertical-radial-decay-power", type=float, default=None
    )
    parser.add_argument("--exterior-entrainment-far-factor", type=float, default=None)
    parser.add_argument(
        "--exterior-entrainment-decay-length-m", type=float, default=None
    )
    parser.add_argument("--wall-overflow-capture-fraction", type=float, default=None)
    parser.add_argument("--wall-overflow-return-rate-s-inv", type=float, default=None)
    parser.add_argument("--wall-overflow-side-gate-length-m", type=float, default=None)
    parser.add_argument("--wall-tangential-bypass-fraction", type=float, default=None)
    parser.add_argument("--wall-tip-bypass-fraction", type=float, default=None)
    parser.add_argument("--wall-tip-bypass-power", type=float, default=None)
    parser.add_argument("--wall-impulse-elevated-transfer-fraction", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-elevated-state-mode",
        choices=("species_marker", "conserved", "hybrid"),
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-fraction", type=float, default=None,
    )
    parser.add_argument("--wall-flux-elevated-transfer-fraction", type=float, default=None)
    parser.add_argument("--wall-impulse-elevated-transport-fraction", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-velocity-m-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-diffusivity-m2-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-outward-velocity-m-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-outward-scale-m", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-outward-decay-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-crosswind-spinup-s", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-transport-interval-s", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-advection-scale", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-advection-duration-s", type=float, default=None,
    )
    parser.add_argument("--wall-impulse-elevated-return-rate-s-inv", type=float, default=None)
    parser.add_argument("--wall-impulse-elevated-center-m", type=float, default=None)
    parser.add_argument("--wall-impulse-elevated-sigma-m", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-elevated-observation-ramp-m", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-settling-velocity-m-s", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-vertical-diffusivity-m2-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-min-center-m", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-rise-velocity-m-s", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-rise-time-s", type=float, default=None,
    )
    parser.add_argument("--wall-impulse-elevated-conserved-center-m", type=float, default=None)
    parser.add_argument("--wall-impulse-elevated-conserved-sigma-m", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-elevated-conserved-settling-velocity-m-s",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-vertical-diffusivity-m2-s",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-min-center-m", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-transport-fraction",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-crosswind-velocity-m-s",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-crosswind-outward-velocity-m-s",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-crosswind-outward-scale-m",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-elevated-conserved-center-mode",
        choices=("fixed", "impact_head"), default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-fraction", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-transport-fraction", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-transport-height-m", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-crosswind-velocity-m-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-crosswind-outward-velocity-m-s",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-crosswind-outward-scale-m",
        type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-return-rate-s-inv", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-center-m", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-sigma-m", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-settling-velocity-m-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-vertical-diffusivity-m2-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-reattached-min-center-m", type=float, default=None,
    )
    parser.add_argument("--wall-impulse-wake-capture-rate-s-inv", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-wake-capture-mode",
        choices=("bulk", "species_weighted"), default=None,
    )
    parser.add_argument("--wall-impulse-wake-return-rate-s-inv", type=float, default=None)
    parser.add_argument("--wall-impulse-wake-transport-fraction", type=float, default=None)
    parser.add_argument("--wall-impulse-wake-crosswind-velocity-m-s", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-wake-crosswind-outward-velocity-m-s", type=float,
        default=None,
    )
    parser.add_argument(
        "--wall-impulse-wake-crosswind-outward-scale-m", type=float,
        default=None,
    )
    parser.add_argument("--wall-impulse-wake-observation-fraction", type=float, default=None)
    parser.add_argument("--wall-impulse-wake-observation-center-m", type=float, default=None)
    parser.add_argument("--wall-impulse-wake-observation-sigma-m", type=float, default=None)
    parser.add_argument(
        "--wall-impulse-wake-settling-velocity-m-s", type=float, default=None,
    )
    parser.add_argument(
        "--wall-impulse-wake-vertical-diffusivity-m2-s", type=float, default=None,
    )
    parser.add_argument("--wall-impulse-wake-min-center-m", type=float, default=None)
    parser.add_argument("--initial-elevated-species-fraction", type=float, default=None)
    parser.add_argument(
        "--cavity-backflow-frame", choices=("cavity", "legacy_x"), default=None,
    )
    parser.add_argument(
        "--source-initial-radius-m", type=float, default=None,
        help="legacy name: full initial 12-sided-prism diameter in metres",
    )
    parser.add_argument(
        "--observation-species-mode",
        choices=("r12_component", "released_mixture"),
        default=None,
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=ROOT / "outputs/ti021-baseline-after-align-default-2026-10-16.json",
    )
    parser.add_argument(
        "--comparison-baseline", type=Path, default=None,
        help="optional frozen promotion baseline; --baseline still supplies run settings",
    )
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    comparison_baseline_path = args.comparison_baseline or args.baseline
    comparison_baseline = json.loads(
        comparison_baseline_path.read_text(encoding="utf-8")
    )
    accepted = inspect.signature(run).parameters
    params = {key: value for key, value in baseline.items() if key in accepted}
    params.update(
        h=baseline["horizontal_cell_size_m"],
        dt=baseline["max_dt_s"],
        source_dt=baseline["source_time_step_s"],
        entrainment=baseline["dense_entrainment_model"],
        cavity_reservoir_capture_rate_s_inv=args.capture_rate_s_inv,
        cavity_reservoir_return_rate_s_inv=args.return_rate_s_inv,
    )
    aliased_overrides = {
        "cavity_reservoir_return_delay_s": args.return_delay_s,
        "cavity_reservoir_return_radius_m": args.return_radius_m,
        "cavity_reservoir_return_offset_m": args.return_offset_m,
        "cavity_reservoir_return_sigma_m": args.return_sigma_m,
        "cavity_shelter_strength": args.shelter_strength,
    }
    direct_overrides = (
        "wall_overflow_capture_fraction",
        "wall_overflow_return_rate_s_inv",
        "wall_overflow_side_gate_length_m",
        "wall_impulse_elevated_transfer_fraction",
        "wall_impulse_elevated_state_mode",
        "wall_impulse_elevated_conserved_fraction",
        "wall_flux_elevated_transfer_fraction",
        "wall_impulse_elevated_transport_fraction",
        "wall_impulse_elevated_crosswind_velocity_m_s",
        "wall_impulse_elevated_crosswind_diffusivity_m2_s",
        "wall_impulse_elevated_crosswind_outward_velocity_m_s",
        "wall_impulse_elevated_crosswind_outward_scale_m",
        "wall_impulse_elevated_crosswind_outward_decay_s",
        "wall_impulse_elevated_crosswind_spinup_s",
        "wall_impulse_elevated_transport_interval_s",
        "wall_impulse_elevated_advection_scale",
        "wall_impulse_elevated_advection_duration_s",
        "wall_impulse_elevated_return_rate_s_inv",
        "wall_impulse_elevated_center_m",
        "wall_impulse_elevated_sigma_m",
        "wall_impulse_elevated_observation_ramp_m",
        "wall_impulse_elevated_settling_velocity_m_s",
        "wall_impulse_elevated_vertical_diffusivity_m2_s",
        "wall_impulse_elevated_min_center_m",
        "wall_impulse_elevated_rise_velocity_m_s",
        "wall_impulse_elevated_rise_time_s",
        "wall_impulse_elevated_conserved_center_m",
        "wall_impulse_elevated_conserved_sigma_m",
        "wall_impulse_elevated_conserved_settling_velocity_m_s",
        "wall_impulse_elevated_conserved_vertical_diffusivity_m2_s",
        "wall_impulse_elevated_conserved_min_center_m",
        "wall_impulse_elevated_conserved_transport_fraction",
        "wall_impulse_elevated_conserved_crosswind_velocity_m_s",
        "wall_impulse_elevated_conserved_crosswind_outward_velocity_m_s",
        "wall_impulse_elevated_conserved_crosswind_outward_scale_m",
        "wall_impulse_elevated_conserved_center_mode",
        "wall_impulse_reattached_fraction",
        "wall_impulse_reattached_transport_fraction",
        "wall_impulse_reattached_transport_height_m",
        "wall_impulse_reattached_crosswind_velocity_m_s",
        "wall_impulse_reattached_crosswind_outward_velocity_m_s",
        "wall_impulse_reattached_crosswind_outward_scale_m",
        "wall_impulse_reattached_return_rate_s_inv",
        "wall_impulse_reattached_center_m",
        "wall_impulse_reattached_sigma_m",
        "wall_impulse_reattached_settling_velocity_m_s",
        "wall_impulse_reattached_vertical_diffusivity_m2_s",
        "wall_impulse_reattached_min_center_m",
        "wall_impulse_wake_capture_rate_s_inv",
        "wall_impulse_wake_capture_mode",
        "wall_impulse_wake_return_rate_s_inv",
        "wall_impulse_wake_transport_fraction",
        "wall_impulse_wake_crosswind_velocity_m_s",
        "wall_impulse_wake_crosswind_outward_velocity_m_s",
        "wall_impulse_wake_crosswind_outward_scale_m",
        "wall_impulse_wake_observation_fraction",
        "wall_impulse_wake_observation_center_m",
        "wall_impulse_wake_observation_sigma_m",
        "wall_impulse_wake_settling_velocity_m_s",
        "wall_impulse_wake_vertical_diffusivity_m2_s",
        "wall_impulse_wake_min_center_m",
        "initial_elevated_species_fraction",
        "cavity_backflow_frame",
        "wake_meander_strength",
        "wake_shedding_strouhal",
        "wake_meander_phase_deg",
        "wake_length_m",
        "wake_base_width_m",
        "wake_spread_slope",
    )
    aliased_overrides.update({
        name: getattr(args, name) for name in direct_overrides
    })
    params.update({
        name: value for name, value in aliased_overrides.items()
        if value is not None
    })
    if args.time_step_s is not None:
        params["dt"] = args.time_step_s
    if args.vertical_observation_component_diagnostics:
        params["vertical_observation_component_diagnostics"] = True
    if args.reservoir_spatial_mode is not None:
        params["cavity_reservoir_spatial_mode"] = args.reservoir_spatial_mode
    if args.reservoir_return_target is not None:
        params["cavity_reservoir_return_target"] = args.reservoir_return_target
    if args.reservoir_reattached_fraction is not None:
        params["cavity_reservoir_reattached_fraction"] = (
            args.reservoir_reattached_fraction
        )
    if args.reservoir_capture_mode is not None:
        params["cavity_reservoir_capture_mode"] = args.reservoir_capture_mode
    if args.reservoir_vertical_capture_mode is not None:
        params["cavity_reservoir_vertical_capture_mode"] = (
            args.reservoir_vertical_capture_mode
        )
    if args.reservoir_content_mode is not None:
        params["cavity_reservoir_content_mode"] = args.reservoir_content_mode
    if args.reservoir_transport_fraction is not None:
        params["cavity_reservoir_transport_fraction"] = (
            args.reservoir_transport_fraction
        )
    if args.reservoir_transport_height_m is not None:
        params["cavity_reservoir_transport_height_m"] = (
            args.reservoir_transport_height_m
        )
    if args.reservoir_transport_mask is not None:
        params["cavity_reservoir_transport_mask"] = args.reservoir_transport_mask
    if args.reservoir_wake_length_m is not None:
        params["cavity_reservoir_wake_length_m"] = args.reservoir_wake_length_m
    if args.reservoir_wake_half_width_m is not None:
        params["cavity_reservoir_wake_half_width_m"] = (
            args.reservoir_wake_half_width_m
        )
    if args.reservoir_crosswind_convergence_rate_s_inv is not None:
        params["cavity_reservoir_crosswind_convergence_rate_s_inv"] = (
            args.reservoir_crosswind_convergence_rate_s_inv
        )
    if args.reservoir_crosswind_target_offset_m is not None:
        params["cavity_reservoir_crosswind_target_offset_m"] = (
            args.reservoir_crosswind_target_offset_m
        )
    if args.local_reservoir_observation_fraction is not None:
        params["local_reservoir_observation_fraction"] = (
            args.local_reservoir_observation_fraction
        )
    if args.local_reservoir_observation_center_m is not None:
        params["local_reservoir_observation_center_m"] = (
            args.local_reservoir_observation_center_m
        )
    if args.local_reservoir_observation_sigma_m is not None:
        params["local_reservoir_observation_sigma_m"] = (
            args.local_reservoir_observation_sigma_m
        )
    if args.source_initial_radius_m is not None:
        params["source_initial_radius_m"] = args.source_initial_radius_m
    if args.observation_species_mode is not None:
        params["observation_species_mode"] = args.observation_species_mode
    if args.entrainment_coefficient is not None:
        params["entrainment_coefficient"] = args.entrainment_coefficient
    if args.entrainment_far_coefficient is not None:
        params["entrainment_far_coefficient"] = args.entrainment_far_coefficient
    if args.entrainment_relaxation_time_s is not None:
        params["entrainment_relaxation_time_s"] = args.entrainment_relaxation_time_s
    if args.horizontal_diffusivity_m2_s is not None:
        params["horizontal_diffusivity_m2_s"] = args.horizontal_diffusivity_m2_s
    if args.crosswind_diffusivity_factor is not None:
        params["crosswind_diffusivity_factor"] = args.crosswind_diffusivity_factor
    if args.wind_angle_deg is not None:
        params["wind_angle_deg"] = args.wind_angle_deg
    if args.wind_history_direction_offset_deg is not None:
        params["wind_history_direction_offset_deg"] = (
            args.wind_history_direction_offset_deg
        )
    if args.wind_history_advection_reference_x_m is not None:
        params["wind_history_advection_reference_x_m"] = (
            args.wind_history_advection_reference_x_m
        )
    if args.wind_history_advection_speed_m_s is not None:
        params["wind_history_advection_speed_m_s"] = (
            args.wind_history_advection_speed_m_s
        )
    if args.layer_wind_speed_m_s is not None:
        params["layer_wind_speed_m_s"] = args.layer_wind_speed_m_s
    if args.friction_velocity_m_s is not None:
        params["friction_velocity_m_s"] = args.friction_velocity_m_s
    if args.gravity_scale is not None:
        params["gravity_scale"] = args.gravity_scale
    if args.source_crosswind_scale is not None:
        params["source_crosswind_scale"] = args.source_crosswind_scale
    if args.source_longitudinal_scale is not None:
        params["source_longitudinal_scale"] = args.source_longitudinal_scale
    if args.puff_gravity_spreading_coefficient is not None:
        params["puff_gravity_spreading_coefficient"] = (
            args.puff_gravity_spreading_coefficient
        )
    if args.puff_form_drag_coefficient is not None:
        params["puff_form_drag_coefficient"] = args.puff_form_drag_coefficient
    if args.receptor_interpolation is not None:
        params["receptor_interpolation"] = args.receptor_interpolation
    if args.vertical_observation_profile is not None:
        params["vertical_observation_profile"] = args.vertical_observation_profile
    if args.vertical_layers is not None:
        params["vertical_layers"] = args.vertical_layers
    if args.vertical_mixing_diffusivity_m2_s is not None:
        params["vertical_mixing_diffusivity_m2_s"] = (
            args.vertical_mixing_diffusivity_m2_s
        )
    if args.multilayer_horizontal_shear_fraction is not None:
        params["multilayer_horizontal_shear_fraction"] = (
            args.multilayer_horizontal_shear_fraction
        )
    if args.multilayer_obstacle_crosswind_shear_rate_s_inv is not None:
        params["multilayer_obstacle_crosswind_shear_rate_s_inv"] = (
            args.multilayer_obstacle_crosswind_shear_rate_s_inv
        )
    if args.multilayer_obstacle_crosswind_decay_length_m is not None:
        params["multilayer_obstacle_crosswind_decay_length_m"] = (
            args.multilayer_obstacle_crosswind_decay_length_m
        )
    if args.multilayer_relative_transport_interval_s is not None:
        params["multilayer_relative_transport_interval_s"] = (
            args.multilayer_relative_transport_interval_s
        )
    if args.vertical_tail_diffusivity_m2_s is not None:
        params["vertical_tail_diffusivity_m2_s"] = (
            args.vertical_tail_diffusivity_m2_s
        )
    if args.initial_vertical_gaussian:
        params["initial_vertical_gaussian"] = True
    if args.height_resolved_wall_flux:
        params["height_resolved_wall_flux"] = True
    if args.cavity_vertical_exponent is not None:
        params["cavity_vertical_exponent"] = args.cavity_vertical_exponent
    if args.cavity_vertical_radial_decay_power is not None:
        params["cavity_vertical_radial_decay_power"] = (
            args.cavity_vertical_radial_decay_power
        )
    if args.exterior_entrainment_far_factor is not None:
        params["exterior_entrainment_far_factor"] = (
            args.exterior_entrainment_far_factor
        )
    if args.exterior_entrainment_decay_length_m is not None:
        params["exterior_entrainment_decay_length_m"] = (
            args.exterior_entrainment_decay_length_m
        )
    if args.wall_tangential_bypass_fraction is not None:
        params["wall_tangential_bypass_fraction"] = (
            args.wall_tangential_bypass_fraction
        )
    if args.wall_tip_bypass_fraction is not None:
        params["wall_tip_bypass_fraction"] = args.wall_tip_bypass_fraction
    if args.wall_tip_bypass_power is not None:
        params["wall_tip_bypass_power"] = args.wall_tip_bypass_power
    args.output.parent.mkdir(parents=True, exist_ok=True)
    config_path = args.output.with_suffix(".configuration.json")
    config_path.write_text(json.dumps(params, indent=2), encoding="utf-8")

    started = time.perf_counter()
    result = run(args.output, **params)
    predicted = np.asarray([
        row["predicted_dose_vol_percent_s"] for row in result["sensor_rows"]
    ])
    old = np.asarray([
        row["predicted_dose_vol_percent_s"]
        for row in comparison_baseline["sensor_rows"]
    ])
    observed = np.asarray([
        row["observed_dose_vol_percent_s"]
        for row in comparison_baseline["sensor_rows"]
    ])
    labels_match = (
        [row["label"] for row in result["sensor_rows"]]
        == [row["label"] for row in comparison_baseline["sensor_rows"]]
    )
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    baseline_groups = _group_metrics(old, observed, dataset.dose_sensors)
    candidate_groups = _group_metrics(predicted, observed, dataset.dose_sensors)
    observation_modes_match = (
        comparison_baseline.get("observation_species_mode", "r12_component")
        == result.get("observation_species_mode", "r12_component")
    )
    comparison = {}
    if observation_modes_match:
        for key in baseline_groups:
            before = baseline_groups[key]["mape_percent"]
            after = candidate_groups[key]["mape_percent"]
            comparison[key] = {
                "baseline_mape_percent": before,
                "candidate_mape_percent": after,
                "delta_percentage_points": after - before,
                "nondegraded": after <= before + 1.0e-12,
            }
    positive = observed > 0.0
    baseline_mape = float(np.mean(np.abs(old[positive] / observed[positive] - 1.0)) * 100.0)
    candidate_mape = float(
        np.mean(np.abs(predicted[positive] / observed[positive] - 1.0)) * 100.0
    )
    diagnostics = result["diagnostics"]
    conservation_ok = bool(
        abs(diagnostics.get("maximum_relative_species_mass_residual", np.inf)) <= 1.0e-8
        and abs(diagnostics.get("mixture_mass_residual_kg", np.inf)) <= 1.0e-6
        and diagnostics.get("momentum_conservation_gate") is True
    )
    result["paired_baseline"] = {
        "source": str(comparison_baseline_path.resolve()),
        "configuration_source": str(args.baseline.resolve()),
        "identical_sensor_order": labels_match,
        "observation_species_modes_match": observation_modes_match,
        "baseline_all_positive_mape_percent": (
            baseline_mape if observation_modes_match else None
        ),
        "candidate_all_positive_mape_percent": candidate_mape,
        "delta_percentage_points": (
            candidate_mape - baseline_mape if observation_modes_match else None
        ),
        "configuration_path": str(config_path.resolve()),
        "group_comparison": comparison,
        "conservation_gate_passed": conservation_ok,
        "promotion_allowed": bool(
            labels_match
            and observation_modes_match
            and conservation_ok
            and candidate_mape < baseline_mape
            and all(item["nondegraded"] for item in comparison.values())
        ),
        "target_10_percent_reached": candidate_mape <= 10.0,
    }
    result["candidate_basis"] = (
        "conservative sub-grid cavity capture, residence delay, and distributed return"
    )
    result["runtime_with_audit_s"] = time.perf_counter() - started
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(result["paired_baseline"], indent=2))


if __name__ == "__main__":
    main()
