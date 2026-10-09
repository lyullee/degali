"""Command-line interface.

The original ships six executables that pass files between them; ``degali``
runs the equivalent stages in one process and prints a summary:

.. code-block:: console

    degali steady    B9.INP --er1 EXAMPLE.ER1 --er2 EXAMPLE.ER2
    degali transient B9T.INP --snapshot 60 --snapshot 120
    degali jet       EX1.INO
    degali jet       EX2.INO --bridge EX2.IN
    degali dose      B9T.INP --at 200 --at 400 --at 800
    degali steady    B9.INP --backend coolprop
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from .run import (
    Receptor,
    run_jet,
    run_jet_to_ground,
    run_steady,
    run_transient,
)


def _write_new_text(path: Path, rendered: str, *, artifact: str) -> None:
    """Create a CLI artifact without replacing a path won by another writer."""
    try:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(rendered)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing {artifact}: {path}"
        ) from error


def _parse_ffi_bound(
    text: str | None,
    *,
    unit: str,
    name: str,
    circular: bool = False,
):
    """Parse one exact value or a lower,nominal,upper CLI bound."""
    if text is None:
        return None
    from .addons.field_contracts import BoundedValue, CircularBoundedValue

    try:
        values = tuple(float(item.strip()) for item in text.split(","))
    except ValueError as error:
        raise ValueError(
            f"{name} must be one value or lower,nominal,upper"
        ) from error
    if len(values) == 1:
        lower = nominal = upper = values[0]
    elif len(values) == 3:
        lower, nominal, upper = values
    else:
        raise ValueError(f"{name} must be one value or lower,nominal,upper")
    source = f"CLI FFI source-state bound: {name}"
    if circular:
        return CircularBoundedValue(nominal, lower, upper, unit=unit, source=source)
    return BoundedValue(nominal, lower, upper, unit=unit, source=source)


def _cmd_ffi_source_state(args) -> int:
    """Run the public FFI source-state envelope audit."""
    from .validation.ffi_source_state import (
        ffi_reference_provenance,
        ffi_source_state_envelope_report,
        FfiSourceState,
        run_ffi_source_state_envelope,
    )
    from .validation.spadeadam import DEFAULT_ROOT, load

    reference_root = Path(args.reference_root or DEFAULT_ROOT)
    trials = load(reference_root)
    try:
        trial = next(item for item in trials if item.test == args.test)
    except StopIteration as error:
        raise ValueError(
            f"FFI reference data contains no test {args.test}"
        ) from error
    state = FfiSourceState.from_trial(trial)
    overrides = {}
    for name, unit in (
        ("rate_kg_s", "kg/s"), ("orifice_m", "m"),
        ("release_height_m", "m"), ("storage_pressure_barg", "barg"),
        ("wind_m_s", "m/s"), ("ambient_temperature_k", "K"),
        ("relative_humidity_pct", "%"), ("ambient_pressure_pa", "Pa"),
        ("wind_reference_height_m", "m"),
    ):
        value = _parse_ffi_bound(
            getattr(args, name), unit=unit, name=name,
        )
        if value is not None:
            overrides[name] = value
    direction = _parse_ffi_bound(
        args.wind_direction_from_deg, unit="deg",
        name="wind_direction_from_deg", circular=True,
    )
    if direction is not None:
        overrides["wind_direction_from_deg"] = direction
    state = replace(state, **overrides)
    readings = None
    if args.radius is not None:
        readings = tuple(
            reading for reading in trial.readings
            if abs(reading.radius - args.radius) < 0.5
        )
    result = run_ffi_source_state_envelope(
        trial, state, corrections=not args.no_corrections,
        readings=readings, max_cases=args.max_cases,
    )
    envelope = ffi_source_state_envelope_report(result)
    payload = {
        "schema": "degali.ffi-source-state-execution.v1",
        "input": {
            "reference_root": str(reference_root.resolve()),
            "reference_provenance": ffi_reference_provenance(reference_root),
            "test": trial.test,
            "radius_m": args.radius,
            "corrections": not args.no_corrections,
            "max_cases": args.max_cases,
        },
        "source_state_envelope": envelope,
    }
    rendered = json.dumps(
        payload, indent=2, ensure_ascii=False, allow_nan=False,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(
                f"refusing to overwrite existing FFI source-state report: {output}"
            )
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="FFI source-state report")
        print(f"FFI source-state report: {output}")
        print(f"source-state status: {result.status}")
    if args.require_complete and result.status != "complete":
        return 2
    return 0


def _cmd_field_screen(args) -> int:
    """Run a strict nominal field case and emit its limited disposition."""
    from .addons.field_case_io import read_field_screening_case_json
    from .addons.field_source_io import read_field_source_schedule_json
    from .addons.field_decision import (
        evaluate_field_operational_screening,
        field_operational_screening_decision_record,
    )
    from .addons.field_report import field_refinement_report, field_screening_report
    from .addons.field_operational_envelope import (
        field_operational_uncertainty_envelope_report,
        run_field_operational_uncertainty_envelope,
    )
    from .addons.field_operational_history import (
        field_operational_joint_measured_history_envelope_report,
        run_field_operational_joint_measured_history_envelope,
        run_field_operational_joint_pressure_driven_history_envelope,
    )
    from .addons.field_historian_io import ImportedPressureDrivenMeasuredHistory
    from .addons.field_operational_phase_transport import (
        field_operational_phase_routing_transport_envelope_report,
        run_field_operational_phase_routing_transport_envelope,
    )
    from .addons.field_operational_sensor_array import (
        field_operational_sensor_array_envelope_report,
        run_field_operational_sensor_array_envelope,
    )
    from .addons.field_operational_source_sensor import (
        field_operational_source_sensor_envelope_report,
        run_field_operational_source_sensor_envelope,
    )
    from .addons.field_workflow import (
        run_field_semi_fv_refinement_study,
        run_field_semi_fv_screening,
    )

    parsed_case = read_field_screening_case_json(args.case)
    request = parsed_case.request
    phase_input = parsed_case.phase_routing_transport
    conditional_review = parsed_case.conditional_review
    atmospheric_source_schedule = None
    atmospheric_source_case = None
    if args.atmospheric_source_case is not None:
        atmospheric_source_case = Path(args.atmospheric_source_case)
        if not (args.uncertainty_envelope or args.joint_source_sensor_envelope):
            raise ValueError(
                "--atmospheric-source-case requires --uncertainty-envelope or "
                "--joint-source-sensor-envelope so every declared source corner "
                "is refined and gated"
            )
        if phase_input is not None or parsed_case.imported_history is not None:
            raise ValueError(
                "--atmospheric-source-case cannot be combined with phase-routing or "
                "measured-history source cases"
            )
        atmospheric_source_schedule = read_field_source_schedule_json(
            atmospheric_source_case
        )
    if args.allow_conditional and conditional_review is None:
        raise ValueError(
            "--allow-conditional requires a strict conditional_review record "
            "in the case JSON"
        )
    if args.sensor_array_envelope and (
        args.uncertainty_envelope
        or args.joint_source_sensor_envelope
        or args.atmospheric_source_case is not None
        or phase_input is not None
        or parsed_case.imported_history is not None
    ):
        raise ValueError(
            "--sensor-array-envelope is a standalone detector-calibration envelope "
            "and cannot be combined with source, historian or phase envelopes"
        )
    if args.joint_source_sensor_envelope and (
        args.uncertainty_envelope
        or args.sensor_array_envelope
        or phase_input is not None
        or parsed_case.imported_history is not None
    ):
        raise ValueError(
            "--joint-source-sensor-envelope cannot be combined with standalone, "
            "historian or phase envelopes"
        )
    if phase_input is not None and not args.uncertainty_envelope:
        raise ValueError(
            "phase-routing field cases require --uncertainty-envelope so every "
            "phase, pool, weather, surface and sensor corner is retained"
        )
    if args.joint_source_sensor_envelope:
        result = run_field_operational_source_sensor_envelope(
            request,
            atmospheric_source_schedule=atmospheric_source_schedule,
            max_cases=args.max_uncertainty_cases,
            table_nodes=161,
            include_refinement=args.refine,
            refinement_factors=tuple(args.refinement_factor),
            relative_tolerance=args.relative_tolerance,
            max_cell_steps=args.max_cell_steps,
            allow_conditional=args.allow_conditional,
        )
        report = field_operational_source_sensor_envelope_report(result)
        decision = result.operational_decision
    elif args.sensor_array_envelope:
        result = run_field_operational_sensor_array_envelope(
            request,
            max_cases=args.max_uncertainty_cases,
            include_refinement=args.refine,
            refinement_factors=tuple(args.refinement_factor),
            relative_tolerance=args.relative_tolerance,
            max_cell_steps=args.max_cell_steps,
            allow_conditional=args.allow_conditional,
        )
        report = field_operational_sensor_array_envelope_report(result)
        decision = result.operational_decision
    elif args.uncertainty_envelope:
        if phase_input is not None:
            result = run_field_operational_phase_routing_transport_envelope(
                request,
                phase_input.config,
                phase_input.pool_launch,
                pool_vertical_sigma_m=phase_input.pool_vertical_sigma_m,
                pool_vertical_sigma_uncertainty=phase_input.pool_vertical_sigma_uncertainty,
                phase_uncertainty=phase_input.phase_uncertainty,
                max_cases=args.max_uncertainty_cases,
                table_nodes=phase_input.table_nodes,
                include_refinement=args.refine,
                refinement_factors=tuple(args.refinement_factor),
                relative_tolerance=args.relative_tolerance,
                max_cell_steps=args.max_cell_steps,
                allow_conditional=args.allow_conditional,
            )
            report = field_operational_phase_routing_transport_envelope_report(result)
            decision = result.operational_decision
        elif parsed_case.imported_history is not None:
            schedule = parsed_case.measured_schedule
            source_request = parsed_case.source_request
            if schedule is None or schedule.quality_assessment is None or source_request is None:
                raise ValueError("measured-history case is missing its approved source-envelope evidence")
            operational_kwargs = {
                "quality_criteria": schedule.quality_assessment.criteria,
                "history_provenance": request.measured_history_provenance,
                "max_cases": args.max_uncertainty_cases,
                "include_refinement": args.refine,
                "refinement_factors": tuple(args.refinement_factor),
                "relative_tolerance": args.relative_tolerance,
                "max_cell_steps": args.max_cell_steps,
                "allow_conditional": args.allow_conditional,
            }
            if isinstance(parsed_case.imported_history, ImportedPressureDrivenMeasuredHistory):
                result = run_field_operational_joint_pressure_driven_history_envelope(
                    source_request, parsed_case.imported_history.history, **operational_kwargs,
                )
            else:
                result = run_field_operational_joint_measured_history_envelope(
                    source_request, parsed_case.imported_history.history, **operational_kwargs,
                )
            report = field_operational_joint_measured_history_envelope_report(result)
            decision = result.operational_decision
        else:
            result = run_field_operational_uncertainty_envelope(
                request,
                atmospheric_source_schedule=atmospheric_source_schedule,
                max_cases=args.max_uncertainty_cases,
                include_refinement=args.refine,
                refinement_factors=tuple(args.refinement_factor),
                relative_tolerance=args.relative_tolerance,
                max_cell_steps=args.max_cell_steps,
                allow_conditional=args.allow_conditional,
            )
            report = field_operational_uncertainty_envelope_report(result)
            decision = result.operational_decision
    else:
        result = (
            run_field_semi_fv_refinement_study(
                request,
                refinement_factors=tuple(args.refinement_factor),
                relative_tolerance=args.relative_tolerance,
                max_cell_steps=args.max_cell_steps,
            )
            if args.refine
            else run_field_semi_fv_screening(request)
        )
        decision = evaluate_field_operational_screening(
            result,
            require_refinement=True,
            allow_conditional=args.allow_conditional,
        )
        report = field_refinement_report(result) if args.refine else field_screening_report(result)
    case = Path(args.case)
    payload = {
        "schema": "degali.field-screening-execution.v1",
        "input": {
            "path": str(case.resolve()),
            "sha256": hashlib.sha256(case.read_bytes()).hexdigest(),
            "refinement_requested": bool(args.refine),
            "uncertainty_envelope_requested": bool(args.uncertainty_envelope),
            "sensor_array_envelope_requested": bool(args.sensor_array_envelope),
            "joint_source_sensor_envelope_requested": bool(
                args.joint_source_sensor_envelope
            ),
            "uncertainty_envelope_kind": (
                None if not (
                    args.uncertainty_envelope
                    or args.sensor_array_envelope
                    or args.joint_source_sensor_envelope
                )
                else (
                    "source_sensor" if args.joint_source_sensor_envelope
                    else "sensor_calibration" if args.sensor_array_envelope
                    else "phase_routing_transport" if phase_input is not None
                    else (
                        (
                            "joint_pressure_driven_history"
                            if isinstance(parsed_case.imported_history, ImportedPressureDrivenMeasuredHistory)
                            else "joint_measured_history"
                        ) if parsed_case.imported_history is not None
                        else (
                            "atmospheric_source" if atmospheric_source_schedule is not None
                            else "nominal_field"
                        )
                    )
                )
            ),
            "atmospheric_source_schedule_case": (
                None if atmospheric_source_case is None else {
                    "path": str(atmospheric_source_case.resolve()),
                    "sha256": hashlib.sha256(
                        atmospheric_source_case.read_bytes()
                    ).hexdigest(),
                    "schedule": atmospheric_source_schedule.as_record(),
                }
            ),
            "measured_history_imported": parsed_case.imported_history is not None,
            "phase_routing_transport_declared": phase_input is not None,
            "validation_evidence": (
                None if request.validation_evidence is None
                else request.validation_evidence.as_record()
            ),
            "conditional_review": (
                None if conditional_review is None else {
                    **conditional_review.as_record(),
                    "applied_to_this_execution": bool(args.allow_conditional),
                }
            ),
            "phase_routing_transport": (
                None if phase_input is None else {
                    "phase_routing": {
                        "post_release_duration_s": phase_input.config.post_release_duration_s,
                        "puff_duration_s": phase_input.config.puff_duration_s,
                        "pool_area_m2": phase_input.config.pool_area_m2,
                        "pool_time_step_s": phase_input.config.pool_time_step_s,
                        "evaporation_coefficient_m2_s": (
                            phase_input.config.evaporation_coefficient_m2_s
                        ),
                        "relative_humidity_pct": phase_input.config.relative_humidity_pct,
                        "roughness_m": phase_input.config.roughness_m,
                        "averaging_time_s": phase_input.config.averaging_time_s,
                        "wind_reference_height_m": (
                            phase_input.config.wind_reference_height_m
                        ),
                        "pool_model": phase_input.config.pool_model,
                        "evidence_id": phase_input.phase_routing_evidence_id,
                        "gas_model_options": dict(phase_input.config.gas_model_options),
                        "phase_model_options": {
                            "maximum_droplet_time_s": (
                                phase_input.config.phase_model_options
                                ["maximum_droplet_time_s"]
                            ),
                        },
                        "droplet_population": {
                            "classes": [
                                {
                                    "diameter_m": item.diameter_m,
                                    "mass_fraction": item.mass_fraction,
                                }
                                for item in phase_input.droplet_population
                            ],
                            "evidence_id": (
                                phase_input.droplet_population_evidence_id
                            ),
                        },
                        "uncertainty": (
                            None if phase_input.phase_uncertainty is None
                            else phase_input.phase_uncertainty.as_record()
                        ),
                    },
                    "pool_launch": {
                        "source_height_m": phase_input.pool_launch.source_height_m,
                        "closure_id": phase_input.pool_launch.closure_id,
                        "evidence_id": phase_input.pool_launch.evidence_id,
                    },
                    "pool_vertical_sigma_m": phase_input.pool_vertical_sigma_m,
                    "pool_vertical_sigma_uncertainty": (
                        None if phase_input.pool_vertical_sigma_uncertainty is None
                        else phase_input.pool_vertical_sigma_uncertainty.as_dict()
                    ),
                    "table_nodes": phase_input.table_nodes,
                }
            ),
            "measured_history_event_id": (
                None if parsed_case.imported_history is None
                else parsed_case.imported_history.provenance.event_id
            ),
            "refinement_factors": list(args.refinement_factor),
            "relative_tolerance": args.relative_tolerance,
            "max_cell_steps": args.max_cell_steps,
            "max_uncertainty_cases": args.max_uncertainty_cases,
            "allow_conditional": bool(args.allow_conditional),
        },
        "field_report": report,
        "operational_screening": field_operational_screening_decision_record(decision),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing field report: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="field report")
        print(f"field report: {output}")
        print(f"operational screening: {decision.status}")
    if args.require_screening and not decision.screening_allowed:
        return 2
    return 0


def _cmd_field_compare(args) -> int:
    """Compare two supplied sensor-prediction CSVs under an explicit basis."""
    from .addons.field_comparison import (
        compare_field_model_sensor_sets,
        field_model_comparison_report,
    )
    from .addons.field_comparison_io import read_field_model_comparison_case_json

    case_path = Path(args.case)
    case = read_field_model_comparison_case_json(case_path)
    comparison = compare_field_model_sensor_sets(
        case.left,
        case.right,
        threshold_mole_fraction=case.threshold_mole_fraction,
        position_tolerance_m=case.position_tolerance_m,
    )
    report = field_model_comparison_report(
        comparison, comparison_evidence=case.comparison_evidence,
    )
    payload = {
        "schema": "degali.field-model-comparison-execution.v1",
        "input": {
            "path": str(case_path.resolve()),
            "sha256": hashlib.sha256(case_path.read_bytes()).hexdigest(),
        },
        "field_model_comparison": report,
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing model comparison: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="model comparison")
        print(f"field model comparison: {output}")
        print(
            "comparison applicability: "
            f"{comparison.applicability.status}; model-selection impact: "
            f"{report['model_selection_impact']['status']}"
        )
    if args.require_comparable and comparison.applicability.status != "accepted":
        return 2
    return 0


def _cmd_field_batch(args) -> int:
    """Run a strict named field batch and retain input-artifact provenance."""
    from .addons.field_batch import export_field_screening_batch, field_batch_decision_summary
    from .addons.field_batch_io import read_field_batch_input_json

    batch_path = Path(args.case)
    parsed = read_field_batch_input_json(batch_path)
    output_directory = Path(args.output_directory)
    exported = export_field_screening_batch(
        parsed.requests,
        output_directory,
        atmospheric_source_schedules=parsed.atmospheric_source_schedules,
        conditional_review_authorizations=parsed.conditional_review_authorizations,
        **parsed.options.as_kwargs(),
    )
    execution_path = output_directory / "batch-execution.json"
    payload = {
        "schema": "degali.field-batch-execution.v1",
        "input": {
            "path": str(batch_path.resolve()),
            "sha256": hashlib.sha256(batch_path.read_bytes()).hexdigest(),
            "batch": parsed.as_record(),
        },
        "manifest_path": str(exported.manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(exported.manifest_path.read_bytes()).hexdigest(),
        "summary": field_batch_decision_summary(exported.cases),
        "cases": [
            {
                "label": item.label,
                "report_file": item.report_path.name,
                "report_sha256": hashlib.sha256(item.report_path.read_bytes()).hexdigest(),
                "applicability_status": item.applicability_status,
                "operational_decision": item.operational_decision.status,
                "screening_allowed": item.operational_decision.screening_allowed,
                "gate_codes": list(item.operational_decision.gate_codes),
            }
            for item in exported.cases
        ],
    }
    if execution_path.exists():
        raise FileExistsError(f"refusing to overwrite existing field batch execution: {execution_path}")
    execution_payload = json.dumps(
        payload, indent=2, ensure_ascii=False, allow_nan=False,
    ) + "\n"
    try:
        with execution_path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(execution_payload)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing to overwrite existing field batch execution: {execution_path}"
        ) from error
    print(f"field batch manifest: {exported.manifest_path}")
    print(f"field batch execution: {execution_path}")
    if args.require_screening and any(
        not item.operational_decision.screening_allowed for item in exported.cases
    ):
        return 2
    return 0


def _cmd_field_validate(args) -> int:
    """Score one declared model set against a fingerprinted observation set."""
    from .addons.field_validation import (
        field_validation_score_record,
        score_field_model_against_validation,
    )
    from .addons.field_validation_io import read_field_validation_case_json

    case_path = Path(args.case)
    case = read_field_validation_case_json(case_path)
    score = score_field_model_against_validation(
        case.model,
        case.dataset,
        threshold_mole_fraction=case.threshold_mole_fraction,
        position_tolerance_m=case.position_tolerance_m,
    )
    payload = {
        "schema": "degali.field-validation-execution.v1",
        "input": {
            "path": str(case_path.resolve()),
            "sha256": hashlib.sha256(case_path.read_bytes()).hexdigest(),
        },
        "field_validation_score": field_validation_score_record(score),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing field validation: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="field validation")
        print(f"field validation: {output}")
        print(f"validation status: {score.status}")
    if args.require_qualified and score.status != "qualified":
        return 2
    return 0


def _cmd_field_audit(args) -> int:
    """Inventory external field-evidence channels without importing observations."""
    from .addons.field_evidence_audit import (
        FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA,
        audit_field_evidence,
        field_evidence_audit_record,
        write_field_evidence_audit_json,
    )

    root = Path(args.root)
    audit = audit_field_evidence(root, max_files=args.max_files)
    payload = {
        "schema": FIELD_EVIDENCE_AUDIT_EXECUTION_SCHEMA,
        "input": {"root": str(root.resolve()), "max_files": args.max_files},
        "field_evidence_audit": field_evidence_audit_record(audit),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        write_field_evidence_audit_json(audit, output)
        print(f"field evidence audit: {output}")
        print(f"audit status: {audit.status}")
    if args.require_complete and audit.status != "candidate_complete":
        return 2
    return 0


def _cmd_field_audit_verify(args) -> int:
    """Verify a saved field-evidence audit without importing observations."""
    from .addons.field_evidence_audit import (
        field_evidence_audit_record,
        read_field_evidence_audit_json,
    )

    artifact = Path(args.artifact)
    audit = read_field_evidence_audit_json(artifact)
    payload = {
        "schema": "degali.field-evidence-audit-verification.v1",
        "input": {
            "path": str(artifact.resolve()),
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
        },
        "field_evidence_audit": field_evidence_audit_record(audit),
        "verification": {
            "candidate_files_verified": True,
            "scan_limit_verified": True,
            "promotion_allowed": False,
        },
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False))
    if args.require_complete and audit.status != "candidate_complete":
        return 2
    return 0


def _cmd_field_evidence_manifest_verify(args) -> int:
    """Verify a manual cross-file evidence manifest and its selected files."""
    from .addons.field_evidence_manifest import (
        field_evidence_manifest_record,
        read_field_evidence_manifest_json,
    )

    artifact = Path(args.artifact)
    try:
        manifest = read_field_evidence_manifest_json(artifact)
    except (OSError, TypeError, ValueError) as exc:
        # Keep the failure machine-readable.  A missing/changed input is a
        # withheld evidence package, not an ambiguous successful verification.
        digest = None
        if artifact.is_file():
            digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        payload = {
            "schema": "degali.field-evidence-manifest-verification.v1",
            "input": {
                "path": str(artifact.resolve()),
                "sha256": digest,
            },
            "field_evidence_manifest": None,
            "evidence_readiness": {
                "status": "withheld",
                "missing_requirements": [str(exc)],
                "promotion_allowed": False,
                "gate_codes": ["manifest_files_withheld", "manifest_not_promoted"],
            },
            "verification": {
                "selected_files_verified": False,
                "promotion_allowed": False,
            },
        }
        print(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False))
        return 2
    payload = {
        "schema": "degali.field-evidence-manifest-verification.v1",
        "input": {
            "path": str(artifact.resolve()),
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
        },
        "field_evidence_manifest": field_evidence_manifest_record(manifest),
        "evidence_readiness": manifest.readiness_record(),
        "verification": {
            "selected_files_verified": True,
            "promotion_allowed": False,
        },
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


def _parse_manifest_selected_paths(values: object) -> dict[str, str]:
    """Parse repeated ``channel=path`` selections for manifest creation."""
    if values is None:
        raise ValueError("at least one --selected-path channel=path is required")
    if not isinstance(values, list):
        raise TypeError("manifest selected paths must be a list")
    selected: dict[str, str] = {}
    for item in values:
        if not isinstance(item, str) or "=" not in item:
            raise ValueError("--selected-path must use channel=path")
        channel, path = item.split("=", 1)
        channel = channel.strip()
        path = path.strip()
        if not channel or not path:
            raise ValueError("--selected-path must use a non-empty channel=path")
        if channel in selected:
            raise ValueError(f"duplicate --selected-path channel: {channel}")
        selected[channel] = path
    return selected


def _cmd_field_evidence_manifest_create(args) -> int:
    """Create a manual five-channel manifest from a verified audit artifact."""
    from .addons.field_evidence_audit import read_field_evidence_audit_json
    from .addons.field_evidence_manifest import (
        FieldEvidenceManifest,
        write_field_evidence_manifest_json,
    )

    audit = read_field_evidence_audit_json(args.audit)
    manifest = FieldEvidenceManifest.from_audit(
        audit,
        manifest_id=args.manifest_id,
        event_id=args.event_id,
        selected_paths=_parse_manifest_selected_paths(args.selected_path),
        dataset_id=args.dataset_id,
        observed_row_count=args.observed_row_count,
        source_boundary_id=args.source_boundary_id,
        weather_id=args.weather_id,
        obstacle_geometry_id=args.obstacle_geometry_id,
        receptor_geometry_id=args.receptor_geometry_id,
        temporal_operator_id=args.temporal_operator_id,
        common_clock_id=args.common_clock_id,
        scope=args.scope,
        sensor_set_id=args.sensor_set_id,
        operator_id=args.operator_id,
        sensor_registry_path=(
            None if args.sensor_registry_path is None
            else str(args.sensor_registry_path)
        ),
        sensor_calibration_status=args.sensor_calibration_status,
    )
    output = write_field_evidence_manifest_json(manifest, args.output)
    print(f"field evidence manifest: {output}")
    print(f"evidence readiness: {manifest.readiness_record()['status']}")
    print("promotion allowed: false")
    return 0


def _cmd_field_verify(args) -> int:
    """Re-verify a deterministic source/model/validation execution artifact."""
    from .addons.field_execution import verify_field_execution_artifact

    artifact = Path(args.artifact)
    verification = verify_field_execution_artifact(artifact)
    rendered = json.dumps(
        verification, indent=2, ensure_ascii=False, allow_nan=False,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing field verification: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="field verification")
        print(f"field execution verification: {output}")
    if args.require_recomputed and not verification["report_recomputed"]:
        return 2
    return 0


def _cmd_field_source(args) -> int:
    """Validate and fingerprint one already-atmospheric source schedule."""
    from .addons.field_source_io import (
        FIELD_SOURCE_SCHEDULE_EXECUTION_SCHEMA,
        field_source_schedule_record,
        read_field_source_schedule_json,
    )

    case_path = Path(args.case)
    schedule = read_field_source_schedule_json(case_path)
    payload = {
        "schema": FIELD_SOURCE_SCHEDULE_EXECUTION_SCHEMA,
        "input": {
            "path": str(case_path.resolve()),
            "sha256": hashlib.sha256(case_path.read_bytes()).hexdigest(),
        },
        "atmospheric_source_schedule": field_source_schedule_record(schedule),
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = Path(args.output)
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing atmospheric source execution: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_new_text(output, rendered, artifact="atmospheric source execution")
        print(f"atmospheric source schedule: {output}")
        print(f"released mass: {schedule.schedule.released_mass_kg:g} kg")
    return 0


def _profile_table(rows: np.ndarray, limit: int = 20) -> str:
    head = f"{'dist (m)':>10} {'mole frac':>11} {'kg/m3':>11} {'T (K)':>8} {'Sz (m)':>8} {'Sy (m)':>8}"
    step = max(1, len(rows) // limit)
    lines = [head, "-" * len(head)]
    for r in rows[::step]:
        lines.append(
            f"{r[0]:10.4g} {r[2]:11.4g} {r[3]:11.4g} {r[5]:8.1f} {r[7]:8.3g} {r[8]:8.3g}"
        )
    return "\n".join(lines)


def _cmd_steady(args) -> int:
    profile, src = run_steady(
        args.deck, er1=args.er1, er2=args.er2, backend=args.backend
    )
    gas = src.case.gas
    fluid = gas.coolprop_name
    props = f"CoolProp ({fluid})" if fluid else src.thermo.backend.name
    print(f"Steady release: {gas.name}, {src.blanket.ess:.4g} kg/s "
          f"[properties: {props}]")
    print(f"  wind-profile exponent alpha = {src.alpha:.5f}")
    print(f"  secondary source: {src.blanket.outl:.3g} m long, "
          f"{src.blanket.outb:.3g} m half-width")
    print(f"  profile: {len(profile.rows)} points, "
          f"{profile.n_dense} in the dense phase, "
          f"transition at {profile.transition:.4g} m")
    for name, level in (("upper", gas.ulc), ("lower", gas.llc)):
        d = profile.distance_to(level)
        shown = "not reached" if np.isnan(d) else f"{d:.4g} m"
        print(f"  distance to the {name} level of concern "
              f"({level * 100:.4g} mol %): {shown}")
    print(f"  mass above the lower level of concern: "
          f"{profile.mass_above_lfl:.5g} kg")
    print(f"  mass between the two levels: {profile.mass_between:.5g} kg")
    print()
    print(_profile_table(profile.rows))
    return 0


def _cmd_transient(args) -> int:
    times = np.array(args.snapshot, dtype=float) if args.snapshot else None
    out = run_transient(
        args.deck, er1=args.er1, er2=args.er2, times=times,
        backend=args.backend,
    )
    gas = out.source.case.gas
    print(f"Transient release: {gas.name}, "
          f"{len(out.field.observers)} observers, {len(out.snapshots)} snapshots")
    head = (f"{'t (s)':>8} {'points':>7} {'extent (m)':>11} "
            f"{'peak mol frac':>14} {'mass > LLC (kg)':>16}")
    print()
    print(head)
    print("-" * len(head))
    for s in out.snapshots:
        print(f"{s.time:8.4g} {len(s):7d} {s.column('dist')[-1]:11.4g} "
              f"{s.column('yc').max():14.5g} {s.mass_above_llc:16.5g}")
    peak = max(out.snapshots, key=lambda s: s.mass_above_llc, default=None)
    if peak is not None:
        print(f"\npeak flammable mass {peak.mass_above_llc:.5g} kg at "
              f"t = {peak.time:.4g} s")
    return 0


def _cmd_dose(args) -> int:
    out = run_transient(
        args.deck, er1=args.er1, er2=args.er2, backend=args.backend
    )
    receptors = [Receptor(x=x) for x in args.at]
    histories = out.dose(receptors)
    print(f"Concentration histories at {len(histories)} receptors")
    head = (f"{'x (m)':>9} {'points':>7} {'peak mol frac':>14} "
            f"{'at t (s)':>9} {'dose (mol frac.s)':>18}")
    print()
    print(head)
    print("-" * len(head))
    for h in histories:
        peak, t = h.peak
        print(f"{h.receptor.x:9.4g} {len(h.rows):7d} {peak:14.5g} "
              f"{t:9.4g} {h.dose():18.5g}")
    return 0


def _cmd_jet(args) -> int:
    if args.bridge is None:
        jet, deck = run_jet(args.deck, backend=args.backend)
        print(f"Jet release: {deck.erate:.4g} kg/s through a "
              f"{deck.diajet:.3g} m orifice at {deck.elejet:.3g} m")
        if not jet.touchdown:
            print("  the plume never reaches the ground at a concentration "
                  "of interest")
            print(f"  final elevation {jet.rows[-1, 1]:.4g} m at "
                  f"{jet.rows[-1, 0]:.4g} m downwind")
            return 0
        print(f"  touchdown at {jet.distance:.5g} m, "
              f"{jet.concentration:.5g} kg/m3, "
              f"half-width {jet.halfwidth:.5g} m")
        print("  pass --bridge with the .IN deck to continue downwind")
        return 0

    profile, jet, src = run_jet_to_ground(
        args.deck, args.bridge, er1=args.er1, er2=args.er2,
        backend=args.backend,
    )
    gas = src.case.gas
    print(f"Jet release: touchdown at {jet.distance:.5g} m")
    print(f"  bridged to a {jet.halfwidth:.4g} m source, "
          f"diluted to {src.case.source.wc[0]:.4g} mass fraction")
    d = profile.distance_to(gas.llc)
    shown = "not reached" if np.isnan(d) else f"{d:.4g} m"
    print(f"  distance to the lower level of concern "
          f"({gas.llc * 100:.4g} mol %): {shown}")
    print()
    print(_profile_table(profile.rows))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="degali",
        description=(
            "DEGALI - Dense Gas Dispersion for Liquid Hydrogen, built on a "
            "verified DEGADIS 2.1 reimplementation"
        ),
    )
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, deck_help):
        sp.add_argument("deck", type=Path, help=deck_help)
        sp.add_argument("--er1", type=Path, help="source-model parameter file")
        sp.add_argument("--er2", type=Path, help="downwind parameter file")
        sp.add_argument(
            "--backend", choices=("legacy", "coolprop"), default="legacy",
            help="property model: 'legacy' reproduces DEGADIS 2.1 exactly, "
                 "'coolprop' uses equations of state and solves the "
                 "thermodynamic inversions accurately (default: legacy)",
        )

    s = sub.add_parser("steady", help="a steady ground-level release")
    common(s, "the .INP input deck")
    s.set_defaults(func=_cmd_steady)

    t = sub.add_parser("transient", help="an unsteady ground-level release")
    common(t, "the .INP input deck")
    t.add_argument(
        "--snapshot", type=float, action="append",
        help="time to report the cloud at, in seconds; repeatable. "
             "Defaults to the window DEGADIS chooses itself.",
    )
    t.set_defaults(func=_cmd_transient)

    d = sub.add_parser("dose", help="concentration history at fixed receptors")
    common(d, "the .INP input deck")
    d.add_argument(
        "--at", type=float, action="append", required=True,
        help="downwind distance of a receptor, in metres; repeatable",
    )
    d.set_defaults(func=_cmd_dose)

    j = sub.add_parser("jet", help="a pressurised release")
    common(j, "the .INO jet deck")
    j.add_argument(
        "--bridge", type=Path,
        help="the .IN jet deck; supplying it continues past touchdown "
             "into the ground-level model",
    )
    j.set_defaults(func=_cmd_jet)

    f = sub.add_parser(
        "field-screen",
        help="run one strict JSON field case with a fail-safe operational disposition",
    )
    f.add_argument(
        "case", type=Path,
        help=(
            "UTF-8 degali.field-screening-input.v1 or "
            "degali.field-phase-routing-screening-input.v1 JSON case"
        ),
    )
    f.add_argument("--output", type=Path, help="new JSON execution record; never overwritten")
    f.add_argument(
        "--no-refinement", dest="refine", action="store_false", default=True,
        help="skip the grid/time study for diagnosis; operational screening remains withheld",
    )
    f.add_argument(
        "--uncertainty-envelope", action="store_true",
        help=(
            "run/refine every declared source, weather, ambient, obstacle and "
            "detector corner (including measured-history corners)"
        ),
    )
    f.add_argument(
        "--sensor-array-envelope", action="store_true",
        help=(
            "run/refine every declared detector calibration corner on one shared "
            "scalar transport field; incompatible with source/history/phase envelopes"
        ),
    )
    f.add_argument(
        "--joint-source-sensor-envelope", action="store_true",
        help=(
            "run/refine the Cartesian source/weather/geometry and supplemental "
            "detector-calibration envelope; incompatible with historian/phase modes"
        ),
    )
    f.add_argument(
        "--atmospheric-source-case", type=Path,
        help=(
            "separate degali.field-atmospheric-schedule-input.v1 JSON; requires "
            "--uncertainty-envelope and is mutually exclusive with historian/phase routing"
        ),
    )
    f.add_argument(
        "--max-uncertainty-cases", type=int, default=64,
        help="hard case limit for --uncertainty-envelope (default: 64)",
    )
    f.add_argument(
        "--refinement-factor", type=int, action="append", default=[1, 2],
        help="positive grid/time refinement factor; repeatable (default: 1, 2)",
    )
    f.add_argument(
        "--relative-tolerance", type=float, default=0.05,
        help="peak/dose refinement tolerance (default: 0.05)",
    )
    f.add_argument(
        "--max-cell-steps", type=int, default=20_000_000,
        help="hard work budget for the grid/time study (default: 20000000)",
    )
    f.add_argument(
        "--allow-conditional", action="store_true",
        help=(
            "apply the case JSON conditional_review authorization to screening only; "
            "never design basis"
        ),
    )
    f.add_argument(
        "--require-screening", action="store_true",
        help="exit 2 when the fail-safe operational screen is withheld",
    )
    f.set_defaults(func=_cmd_field_screen)

    b = sub.add_parser(
        "field-batch",
        help="run named strict field cases and emit a fail-safe batch manifest",
    )
    b.add_argument(
        "case", type=Path,
        help="UTF-8 degali.field-batch-input.v1 JSON case",
    )
    b.add_argument(
        "--output-directory", type=Path, required=True,
        help="new or empty directory for per-case reports and manifest",
    )
    b.add_argument(
        "--require-screening", action="store_true",
        help="exit 2 when any named case remains withheld",
    )
    b.set_defaults(func=_cmd_field_batch)

    c = sub.add_parser(
        "field-compare",
        help="compare two declared fixed-sensor model CSVs without hidden operator changes",
    )
    c.add_argument(
        "case", type=Path,
        help="UTF-8 degali.field-model-comparison-input.v1 JSON case",
    )
    c.add_argument("--output", type=Path, help="new JSON comparison record; never overwritten")
    c.add_argument(
        "--require-comparable", action="store_true",
        help="exit 2 unless source, weather, sensor geometry and temporal operator match",
    )
    c.set_defaults(func=_cmd_field_compare)

    v = sub.add_parser(
        "field-validate",
        help="score a declared field model against fingerprinted observed sensors",
    )
    v.add_argument(
        "case", type=Path,
        help="UTF-8 degali.field-validation-input.v1 JSON case",
    )
    v.add_argument("--output", type=Path, help="new JSON validation record; never overwritten")
    v.add_argument(
        "--require-qualified", action="store_true",
        help="exit 2 unless the matched observed-sensor score is qualified",
    )
    v.set_defaults(func=_cmd_field_validate)

    a = sub.add_parser(
        "field-audit",
        help="scan external CSV/JSON files for field-validation evidence channels",
    )
    a.add_argument(
        "root", type=Path,
        help="read-only directory to scan; no raw observations are imported",
    )
    a.add_argument("--output", type=Path, help="new JSON audit record; never overwritten")
    a.add_argument(
        "--max-files", type=int, default=20000,
        help="hard structured-file scan limit (default: 20000)",
    )
    a.add_argument(
        "--require-complete", action="store_true",
        help="exit 2 unless all five channels appear somewhere in the scan",
    )
    a.set_defaults(func=_cmd_field_audit)

    av = sub.add_parser(
        "field-audit-verify",
        help="verify a saved field-evidence audit artifact and its candidate files",
    )
    av.add_argument(
        "artifact", type=Path,
        help="saved degali.field-evidence-audit-execution.v1 JSON artifact",
    )
    av.add_argument(
        "--require-complete", action="store_true",
        help="exit 2 unless the saved audit status is candidate_complete",
    )
    av.set_defaults(func=_cmd_field_audit_verify)

    am = sub.add_parser(
        "field-evidence-manifest-verify",
        help="verify a hash-pinned cross-file field-evidence manifest",
    )
    am.add_argument(
        "artifact", type=Path,
        help="saved degali.field-evidence-manifest-execution.v1 JSON artifact",
    )
    am.set_defaults(func=_cmd_field_evidence_manifest_verify)

    ac = sub.add_parser(
        "field-evidence-manifest-create",
        help=(
            "create a hash-pinned five-channel evidence manifest from a "
            "verified field-audit artifact"
        ),
    )
    ac.add_argument(
        "audit", type=Path,
        help="saved degali.field-evidence-audit-execution.v1 JSON artifact",
    )
    ac.add_argument(
        "--selected-path", action="append", required=True,
        metavar="CHANNEL=PATH",
        help=(
            "explicit audited candidate selection; repeat exactly once for "
            "each of the five channels"
        ),
    )
    ac.add_argument("--manifest-id", required=True)
    ac.add_argument("--event-id", required=True)
    ac.add_argument("--dataset-id", required=True)
    ac.add_argument("--observed-row-count", type=int, required=True)
    ac.add_argument("--source-boundary-id", required=True)
    ac.add_argument("--weather-id", required=True)
    ac.add_argument("--obstacle-geometry-id", required=True)
    ac.add_argument("--receptor-geometry-id", required=True)
    ac.add_argument("--temporal-operator-id", required=True)
    ac.add_argument("--common-clock-id", required=True)
    ac.add_argument(
        "--sensor-set-id", default="legacy-unspecified",
        help=(
            "sensor registry/calibration set ID covering the receptor file; "
            "omit only for a legacy conditional manifest"
        ),
    )
    ac.add_argument(
        "--operator-id", default="legacy-unspecified",
        help=(
            "accountable operator/curator ID; omit only for a legacy conditional manifest"
        ),
    )
    ac.add_argument(
        "--sensor-registry-path", type=Path,
        help=(
            "optional sensor registry/calibration file under the audit root; "
            "its SHA-256 is pinned when supplied"
        ),
    )
    ac.add_argument(
        "--sensor-calibration-status",
        choices=("certified", "specification_only", "missing", "legacy-unspecified"),
        default="legacy-unspecified",
        help=(
            "event-level sensor calibration status; specification_only or "
            "missing keeps the manifest conditional"
        ),
    )
    ac.add_argument(
        "--scope", choices=("lh2_free_field", "lh2_obstacle_transport"),
        default="lh2_obstacle_transport",
    )
    ac.add_argument(
        "--output", type=Path, required=True,
        help="new JSON manifest execution record; never overwritten",
    )
    ac.set_defaults(func=_cmd_field_evidence_manifest_create)

    ev = sub.add_parser(
        "field-verify",
        help=(
            "re-verify a source, model-comparison, or validation execution "
            "artifact and its referenced files"
        ),
    )
    ev.add_argument(
        "artifact", type=Path,
        help="saved deterministic field execution JSON artifact",
    )
    ev.add_argument(
        "--output", type=Path,
        help="new JSON verification record; never overwritten",
    )
    ev.add_argument(
        "--require-recomputed", action="store_true",
        help="exit 2 unless the deterministic execution report was recomputed",
    )
    ev.set_defaults(func=_cmd_field_verify)

    q = sub.add_parser(
        "field-source",
        help="validate and fingerprint an already-atmospheric source-rate CSV",
    )
    q.add_argument(
        "case", type=Path,
        help="UTF-8 degali.field-atmospheric-schedule-input.v1 JSON case",
    )
    q.add_argument("--output", type=Path, help="new JSON execution record; never overwritten")
    q.set_defaults(func=_cmd_field_source)

    fs = sub.add_parser(
        "ffi-source-state",
        help=(
            "audit explicitly bounded source/weather/ambient corners against "
            "public FFI horizontal-release sensors"
        ),
    )
    fs.add_argument("--test", type=int, default=6, help="FFI outdoor test number (default: 6)")
    fs.add_argument(
        "--reference-root", type=Path,
        help="reference/spadeadam directory (default: package reference root)",
    )
    fs.add_argument(
        "--radius", type=float,
        help="retain readings within 0.5 m of this arc radius",
    )
    fs.add_argument(
        "--max-cases", type=int, default=64,
        help="hard deterministic corner limit (default: 64)",
    )
    fs.add_argument(
        "--no-corrections", action="store_true",
        help="run the as-shipped free-field branch instead of the corrected branch",
    )
    for flag, dest, unit in (
        ("--rate-bounds", "rate_kg_s", "kg/s"),
        ("--orifice-bounds", "orifice_m", "m"),
        ("--release-height-bounds", "release_height_m", "m"),
        ("--storage-pressure-bounds", "storage_pressure_barg", "barg"),
        ("--wind-bounds", "wind_m_s", "m/s"),
        ("--ambient-temperature-bounds", "ambient_temperature_k", "K"),
        ("--humidity-bounds", "relative_humidity_pct", "%"),
        ("--ambient-pressure-bounds", "ambient_pressure_pa", "Pa"),
        ("--wind-reference-height-bounds", "wind_reference_height_m", "m"),
    ):
        fs.add_argument(
            flag, dest=dest,
            help=f"{dest}: one value or lower,nominal,upper ({unit})",
        )
    fs.add_argument(
        "--wind-direction-from-deg", dest="wind_direction_from_deg",
        help="wind-from direction: one value or lower,nominal,upper (deg)",
    )
    fs.add_argument("--output", type=Path, help="new JSON execution record; never overwritten")
    fs.add_argument(
        "--require-complete", action="store_true",
        help="exit 2 unless every selected sensor is available in every corner",
    )
    fs.set_defaults(func=_cmd_ffi_source_state)

    args = p.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        print(f"degali: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
