"""Evidence-safe intake for Thorney Island Phase-II Trial 021.

The generic SMEDIS reader cannot safely interpret this workbook: the
``Concentration_Sensor`` table contains integrated dose rather than a mean
concentration, and the ``Arc Positions`` table contains both maximum dose and
maximum concentration.  This module therefore gives TI021 a narrow, explicit
schema and preserves the published site coordinates, timing summaries and
known conflicts with the final HSE report.

TI021 is an R12 dense-gas puff behind a semi-circular impermeable fence.  It is
valuable obstacle-transport evidence, but it is not hydrogen validation and it
does not contain complete concentration time histories.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence


THORNEY_TI021_EVIDENCE_SCHEMA = "degali.thorney-ti021-evidence.v12"


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _optional(value: object) -> float | None:
    if value in (None, ""):
        return None
    result = float(value)
    if not math.isfinite(result) or result in {-999.0, -99.0, -9.99}:
        return None
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Ti021ArcObservation:
    """One published radial-arc reduction."""

    distance_m: float
    height_m: float
    maximum_dose_percent_s: float
    maximum_concentration_percent: float
    peak_time_s: float
    lateral_sigma_m: float
    arrival_time_s: float
    departure_time_s: float | None

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if value is None:
                continue
            number = _finite(value, name)
            if number < 0.0:
                raise ValueError(f"{name} cannot be negative")


@dataclass(frozen=True)
class Ti021DoseSensor:
    """One dose/timing sensor in release-centred site-grid coordinates."""

    site_x_m: float
    site_y_m: float
    local_x_m: float
    local_y_m: float
    height_m: float
    dose_percent_s: float
    arrival_time_s: float | None
    departure_time_s: float | None

    def __post_init__(self) -> None:
        for name in (
            "site_x_m", "site_y_m", "local_x_m", "local_y_m", "height_m",
            "dose_percent_s",
        ):
            _finite(getattr(self, name), name)
        if self.height_m < 0.0 or self.dose_percent_s < 0.0:
            raise ValueError("sensor height and dose cannot be negative")
        for name in ("arrival_time_s", "departure_time_s"):
            value = getattr(self, name)
            if value is not None and _finite(value, name) < 0.0:
                raise ValueError(f"{name} cannot be negative")


@dataclass(frozen=True)
class ThorneyTi021Dataset:
    """Typed TI021 workbook content plus immutable file provenance."""

    workbook_path: Path
    workbook_sha256: str
    release_origin_site_m: tuple[float, float, float]
    substance: str
    release_type: str
    initial_r12_mole_fraction: float
    effective_molecular_weight_kg_mol: float
    initial_temperature_k: float
    initial_pressure_pa: float
    workbook_mass_proxy_kg: float
    workbook_mass_proxy_duration_s: float
    workbook_wind_speed_m_s: float
    workbook_wind_reference_height_m: float
    workbook_wind_from_north_deg: float
    ideal_wind_from_north_deg: float
    site_y_axis_from_north_deg: float
    friction_velocity_m_s: float
    roughness_length_m: float
    stability_class: str
    workbook_wall_radius_m: float
    workbook_wall_center_site_m: tuple[float, float]
    workbook_wall_height_m: float
    workbook_porosity_field: float
    arc_observations: tuple[Ti021ArcObservation, ...]
    dose_sensors: tuple[Ti021DoseSensor, ...]

    @property
    def workbook_release_rate_kg_s(self) -> float:
        """Raw workbook field labelled ``release rate (kg/s)``.

        The field is retained as evidence only.  It is not silently converted
        into an atmospheric source boundary.
        """

        return self.workbook_mass_proxy_kg

    @property
    def workbook_release_duration_s(self) -> float:
        """Raw workbook release-duration field, retained without reinterpretation."""

        return self.workbook_mass_proxy_duration_s

    @property
    def workbook_release_mass_proxy_kg(self) -> float:
        """Rate-times-duration mass proxy, not a validated atmospheric mass."""

        return (
            self.workbook_release_rate_kg_s
            * self.workbook_release_duration_s
        )

    @property
    def positive_dose_sensor_count(self) -> int:
        return sum(sensor.dose_percent_s > 0.0 for sensor in self.dose_sensors)

    @property
    def validation_scope(self) -> str:
        return "r12_curved_impermeable_fence_dose_and_timing_component"

    @property
    def quantitative_hydrogen_prediction_allowed(self) -> bool:
        return False


def _row_label(row: Sequence[object]) -> str:
    return str(row[0]).strip().lower() if row else ""


def _lookup(rows: Sequence[Sequence[object]], fragment: str) -> object:
    fragment = fragment.lower()
    for row in rows:
        if fragment in _row_label(row):
            return row[1] if len(row) > 1 else ""
    raise ValueError(f"TI021 workbook lacks field containing {fragment!r}")


def _lookup_exact(rows: Sequence[Sequence[object]], label: str) -> object:
    label = label.strip().lower()
    for row in rows:
        if _row_label(row) == label:
            return row[1] if len(row) > 1 else ""
    raise ValueError(f"TI021 workbook lacks field {label!r}")


def _section_index(rows: Sequence[Sequence[object]], fragment: str) -> int:
    fragment = fragment.lower()
    for index, row in enumerate(rows):
        if fragment in _row_label(row):
            return index
    raise ValueError(f"TI021 workbook lacks section containing {fragment!r}")


def _split_concatenated_xy(value: object) -> tuple[float, float]:
    """Recover the known six-digit ``400200`` centre cell without guessing.

    The SMEDIS spreadsheet collapsed the comma between the two three-digit
    Thorney site coordinates.  Other layouts are rejected rather than split
    heuristically.
    """

    number = _finite(value, "wall centre coordinates")
    integer = int(round(number))
    if not math.isclose(number, integer, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError("TI021 wall centre coordinate field is not integral")
    text = f"{integer:06d}"
    if len(text) != 6:
        raise ValueError("TI021 wall centre coordinate field is not six digits")
    return float(text[:3]), float(text[3:])


def _parse_rows(
    rows: Sequence[Sequence[object]], *, workbook_path: Path,
    workbook_sha256: str,
) -> ThorneyTi021Dataset:
    dataset = str(_lookup(rows, "dataset_reference")).strip()
    trial = str(_lookup(rows, "trial_identification")).strip()
    if dataset.lower() != "thorney island" or trial.lower() != "trial 021":
        raise ValueError("workbook is not Thorney Island Trial 021")

    origin_x = _finite(_lookup(rows, "release point x"), "release point x")
    # In this workbook the following row is labelled only ``y``.
    origin_y = _finite(_lookup_exact(rows, "y"), "release point y")
    origin_z = _finite(_lookup_exact(rows, "z"), "release point z")

    arc_start = _section_index(rows, "#arc positions") + 2
    sensor_header = _section_index(rows, "# concentration_sensor")
    arcs: list[Ti021ArcObservation] = []
    for row in rows[arc_start:sensor_header]:
        if not row or row[0] in ("", None):
            continue
        values = list(row) + [""] * (8 - len(row))
        arcs.append(Ti021ArcObservation(
            distance_m=_finite(values[0], "arc distance"),
            height_m=_finite(values[1], "arc height"),
            maximum_dose_percent_s=_finite(values[2], "arc maximum dose"),
            maximum_concentration_percent=_finite(
                values[3], "arc maximum concentration"
            ),
            peak_time_s=_finite(values[4], "arc peak time"),
            lateral_sigma_m=_finite(values[5], "arc lateral sigma"),
            arrival_time_s=_finite(values[6], "arc arrival time"),
            departure_time_s=_optional(values[7]),
        ))

    sensor_start = sensor_header + 2
    sensors: list[Ti021DoseSensor] = []
    for row in rows[sensor_start:]:
        if not row or row[0] in ("", None):
            continue
        values = list(row) + [""] * (6 - len(row))
        site_x = _finite(values[0], "sensor site x")
        site_y = _finite(values[1], "sensor site y")
        sensors.append(Ti021DoseSensor(
            site_x_m=site_x,
            site_y_m=site_y,
            local_x_m=site_x - origin_x,
            local_y_m=site_y - origin_y,
            height_m=_finite(values[2], "sensor height"),
            dose_percent_s=_finite(values[3], "sensor dose"),
            arrival_time_s=_optional(values[4]),
            departure_time_s=_optional(values[5]),
        ))

    if len(arcs) != 6 or len(sensors) != 145:
        raise ValueError(
            "TI021 workbook layout changed: expected 6 arc rows and 145 sensors"
        )

    return ThorneyTi021Dataset(
        workbook_path=workbook_path,
        workbook_sha256=workbook_sha256,
        release_origin_site_m=(origin_x, origin_y, origin_z),
        substance=str(_lookup(rows, "substance")).strip(),
        release_type=str(_lookup(rows, "release type")).strip(),
        initial_r12_mole_fraction=_finite(
            _lookup(rows, "initial concentration"), "initial concentration"
        ),
        effective_molecular_weight_kg_mol=1.0e-3 * _finite(
            _lookup(rows, "effective mol wt"), "effective molecular weight"
        ),
        initial_temperature_k=273.15 + _finite(
            _lookup(rows, "exit temperature"), "initial temperature"
        ),
        initial_pressure_pa=100.0 * _finite(
            _lookup(rows, "exit pressure"), "initial pressure"
        ),
        workbook_mass_proxy_kg=_finite(
            _lookup(rows, "release rate"), "workbook mass proxy"
        ),
        workbook_mass_proxy_duration_s=_finite(
            _lookup(rows, "release duration"), "workbook mass proxy duration"
        ),
        workbook_wind_speed_m_s=_finite(
            _lookup(rows, "site average windspeed"), "workbook wind speed"
        ),
        workbook_wind_reference_height_m=_finite(
            _lookup(rows, "reference height for wind"), "wind reference height"
        ),
        workbook_wind_from_north_deg=_finite(
            _lookup(rows, "site av wind direction"), "site wind direction"
        ),
        ideal_wind_from_north_deg=_finite(
            _lookup(rows, "ideal wind direction"), "ideal wind direction"
        ),
        site_y_axis_from_north_deg=_finite(
            _lookup(rows, "y-axis orientation"), "site y-axis orientation"
        ),
        friction_velocity_m_s=_finite(
            _lookup(rows, "friction velocity"), "friction velocity"
        ),
        roughness_length_m=_finite(
            _lookup(rows, "surface roughness"), "roughness length"
        ),
        stability_class=str(_lookup(rows, "stability class (from data")).strip(),
        workbook_wall_radius_m=_finite(_lookup(rows, "radius"), "wall radius"),
        workbook_wall_center_site_m=_split_concatenated_xy(
            _lookup(rows, "centre coordinates")
        ),
        workbook_wall_height_m=_finite(_lookup(rows, "height (m)"), "wall height"),
        workbook_porosity_field=_finite(_lookup(rows, "porosity"), "porosity field"),
        arc_observations=tuple(arcs),
        dose_sensors=tuple(sensors),
    )


def read_thorney_ti021_workbook(path: str | Path) -> ThorneyTi021Dataset:
    """Read the user-held SMEDIS TI021 workbook without unit reinterpretation."""

    import xlrd

    source = Path(path).resolve()
    workbook = xlrd.open_workbook(str(source))
    if workbook.nsheets != 1:
        raise ValueError("TI021 workbook must contain exactly one sheet")
    sheet = workbook.sheet_by_index(0)
    rows = tuple(
        tuple(sheet.cell_value(row, column) for column in range(sheet.ncols))
        for row in range(sheet.nrows)
    )
    return _parse_rows(
        rows, workbook_path=source, workbook_sha256=_sha256(source)
    )


def thorney_ti021_evidence_record(
    dataset: ThorneyTi021Dataset, *, report_path: str | Path,
    curved_wall_geometry_audit_path: str | Path | None = None,
    initial_cloud_audit_path: str | Path | None = None,
    r12_property_audit_path: str | Path | None = None,
    full_scale_preflight_path: str | Path | None = None,
    turbulent_mixing_policy_path: str | Path | None = None,
    full_transient_plan_path: str | Path | None = None,
    nested_handoff_audit_path: str | Path | None = None,
    outer_inlet_coupling_audit_path: str | Path | None = None,
    outer_momentum_overlap_audit_path: str | Path | None = None,
    handoff_streaming_preflight_path: str | Path | None = None,
    coarse_inner_transient_path: str | Path | None = None,
    coarse_inner_score_path: str | Path | None = None,
    completed_inner_transient_path: str | Path | None = None,
    inner_transient_score_path: str | Path | None = None,
    corrector_tolerance_audit_path: str | Path | None = None,
    strict_peak_window_audit_path: str | Path | None = None,
    full_event_grid_comparison_path: str | Path | None = None,
) -> dict[str, object]:
    """Build the fail-closed evidence decision against the final HSE report.

    The report facts below are deliberately named and versioned.  They were
    transcribed from sections 8.3, 18.3 and the Trial-021 summary rather than
    inferred from the concentration observations.
    """

    report = Path(report_path).resolve()
    geometry_audit = None
    geometry_audit_provenance = None
    if curved_wall_geometry_audit_path is not None:
        geometry_path = Path(curved_wall_geometry_audit_path).resolve()
        geometry_audit = json.loads(geometry_path.read_text(encoding="utf-8"))
        if geometry_audit.get("schema") != (
            "degali.thorney-ti021-curved-wall-geometry-audit.v1"
        ):
            raise ValueError("unsupported TI021 curved-wall geometry audit schema")
        geometry_audit_provenance = {
            "path": str(geometry_path),
            "sha256": _sha256(geometry_path),
            "status": geometry_audit.get("status"),
        }
    geometry_verified = bool(
        geometry_audit is not None
        and geometry_audit.get("status") == "accepted_geometry_component"
        and geometry_audit.get("all_gates_pass") is True
    )
    initial_cloud_audit = None
    initial_cloud_audit_provenance = None
    if initial_cloud_audit_path is not None:
        initial_path = Path(initial_cloud_audit_path).resolve()
        initial_cloud_audit = json.loads(initial_path.read_text(encoding="utf-8"))
        if initial_cloud_audit.get("schema") != (
            "degali.thorney-ti021-initial-cloud-audit.v1"
        ):
            raise ValueError("unsupported TI021 initial-cloud audit schema")
        initial_cloud_audit_provenance = {
            "path": str(initial_path),
            "sha256": _sha256(initial_path),
            "status": initial_cloud_audit.get("status"),
        }
    initial_cloud_verified = bool(
        initial_cloud_audit is not None
        and initial_cloud_audit.get("status")
        == "accepted_initial_state_component"
        and all(initial_cloud_audit.get("gates", {}).values())
    )
    property_audit = None
    property_audit_provenance = None
    if r12_property_audit_path is not None:
        property_path = Path(r12_property_audit_path).resolve()
        property_audit = json.loads(property_path.read_text(encoding="utf-8"))
        if property_audit.get("schema") != (
            "degali.thorney-ti021-r12-property-audit.v1"
        ):
            raise ValueError("unsupported TI021 R12 property audit schema")
        property_audit_provenance = {
            "path": str(property_path),
            "sha256": _sha256(property_path),
            "status": property_audit.get("status"),
        }
    molecular_properties_verified = bool(
        property_audit is not None
        and property_audit.get("status")
        == "accepted_molecular_property_component"
        and all(property_audit.get("gates", {}).values())
    )
    full_scale_preflight = None
    full_scale_preflight_provenance = None
    if full_scale_preflight_path is not None:
        preflight_path = Path(full_scale_preflight_path).resolve()
        full_scale_preflight = json.loads(
            preflight_path.read_text(encoding="utf-8")
        )
        if full_scale_preflight.get("schema") != (
            "degali.thorney-ti021-full-scale-preflight.v1"
        ):
            raise ValueError("unsupported TI021 full-scale preflight schema")
        full_scale_preflight_provenance = {
            "path": str(preflight_path),
            "sha256": _sha256(preflight_path),
            "status": full_scale_preflight.get("status"),
        }
    full_scale_preflight_verified = bool(
        full_scale_preflight is not None
        and full_scale_preflight.get("status")
        == "accepted_full_scale_single_step_preflight"
        and full_scale_preflight.get("all_gates_pass") is True
    )
    turbulent_mixing_policy = None
    turbulent_mixing_policy_provenance = None
    if turbulent_mixing_policy_path is not None:
        mixing_path = Path(turbulent_mixing_policy_path).resolve()
        turbulent_mixing_policy = json.loads(
            mixing_path.read_text(encoding="utf-8")
        )
        if turbulent_mixing_policy.get("schema") != (
            "degali.thorney-ti021-turbulent-mixing-policy.v1"
        ):
            raise ValueError("unsupported TI021 turbulent-mixing policy schema")
        turbulent_mixing_policy_provenance = {
            "path": str(mixing_path),
            "sha256": _sha256(mixing_path),
            "status": turbulent_mixing_policy.get("status"),
        }
    turbulent_mixing_verified = bool(
        turbulent_mixing_policy is not None
        and turbulent_mixing_policy.get("status")
        == "accepted_concentration_blind_turbulent_mixing_envelope_component"
        and turbulent_mixing_policy.get("all_gates_pass") is True
        and len(
            turbulent_mixing_policy.get("frozen_policy", {}).get(
                "ensemble", []
            )
        ) == 6
    )
    full_transient_plan = None
    full_transient_plan_provenance = None
    if full_transient_plan_path is not None:
        plan_path = Path(full_transient_plan_path).resolve()
        full_transient_plan = json.loads(plan_path.read_text(encoding="utf-8"))
        if full_transient_plan.get("schema") not in {
            "degali.thorney-ti021-full-transient-plan.v1",
            "degali.thorney-ti021-full-transient-plan.v2",
        }:
            raise ValueError("unsupported TI021 full-transient plan schema")
        full_transient_plan_provenance = {
            "path": str(plan_path),
            "sha256": _sha256(plan_path),
            "status": full_transient_plan.get("status"),
        }
    full_transient_plan_verified = bool(
        full_transient_plan is not None
        and full_transient_plan.get("status") in {
            "accepted_nested_execution_plan_pending_conservative_handoff",
            "accepted_coordinate_consistent_nested_execution_plan",
        }
        and full_transient_plan.get("all_gates_pass") is True
        and (
            full_transient_plan.get("coverage", {}).get(
                "nested_union", {}
            ).get("combined_unique_sensor_count") == 145
            or full_transient_plan.get("coverage", {}).get(
                "outer_sensor_domain", {}
            ).get("covered_sensor_count") == 145
        )
    )
    nested_handoff_audit = None
    nested_handoff_audit_provenance = None
    if nested_handoff_audit_path is not None:
        handoff_path = Path(nested_handoff_audit_path).resolve()
        nested_handoff_audit = json.loads(
            handoff_path.read_text(encoding="utf-8")
        )
        if nested_handoff_audit.get("schema") != (
            "degali.thorney-ti021-nested-handoff-audit.v1"
        ):
            raise ValueError("unsupported TI021 nested-handoff audit schema")
        nested_handoff_audit_provenance = {
            "path": str(handoff_path),
            "sha256": _sha256(handoff_path),
            "status": nested_handoff_audit.get("status"),
        }
    handoff_remap_verified = bool(
        nested_handoff_audit is not None
        and nested_handoff_audit.get("status")
        == "accepted_conservative_handoff_remap_operator_only"
        and nested_handoff_audit.get("all_gates_pass") is True
    )
    outer_inlet_audit = None
    outer_inlet_audit_provenance = None
    if outer_inlet_coupling_audit_path is not None:
        outer_inlet_path = Path(outer_inlet_coupling_audit_path).resolve()
        outer_inlet_audit = json.loads(
            outer_inlet_path.read_text(encoding="utf-8")
        )
        if outer_inlet_audit.get("schema") != (
            "degali.thorney-ti021-outer-inlet-coupling-audit.v1"
        ):
            raise ValueError("unsupported TI021 outer-inlet coupling audit schema")
        outer_inlet_audit_provenance = {
            "path": str(outer_inlet_path),
            "sha256": _sha256(outer_inlet_path),
            "status": outer_inlet_audit.get("status"),
        }
    outer_binary_boundary_verified = bool(
        outer_inlet_audit is not None
        and outer_inlet_audit.get("status")
        == "accepted_outer_binary_boundary_component_pending_momentum_overlap"
        and outer_inlet_audit.get("all_gates_pass") is True
    )
    outer_momentum_audit = None
    outer_momentum_audit_provenance = None
    if outer_momentum_overlap_audit_path is not None:
        outer_momentum_path = Path(
            outer_momentum_overlap_audit_path
        ).resolve()
        outer_momentum_audit = json.loads(
            outer_momentum_path.read_text(encoding="utf-8")
        )
        if outer_momentum_audit.get("schema") != (
            "degali.thorney-ti021-outer-momentum-overlap-audit.v1"
        ):
            raise ValueError(
                "unsupported TI021 outer-momentum overlap audit schema"
            )
        outer_momentum_audit_provenance = {
            "path": str(outer_momentum_path),
            "sha256": _sha256(outer_momentum_path),
            "status": outer_momentum_audit.get("status"),
        }
    outer_momentum_overlap_verified = bool(
        outer_momentum_audit is not None
        and outer_momentum_audit.get("status")
        == "accepted_outer_momentum_and_overlap_component"
        and outer_momentum_audit.get("all_gates_pass") is True
        and outer_momentum_audit.get("configuration", {}).get(
            "observation_value_fields_accessed"
        ) == []
    )
    conservative_nested_handoff_verified = bool(
        handoff_remap_verified
        and outer_binary_boundary_verified
        and outer_momentum_overlap_verified
    )
    streaming_preflight = None
    streaming_preflight_provenance = None
    if handoff_streaming_preflight_path is not None:
        streaming_path = Path(handoff_streaming_preflight_path).resolve()
        streaming_preflight = json.loads(
            streaming_path.read_text(encoding="utf-8")
        )
        if streaming_preflight.get("schema") != (
            "degali.thorney-ti021-handoff-streaming-preflight.v1"
        ):
            raise ValueError(
                "unsupported TI021 handoff-streaming preflight schema"
            )
        streaming_preflight_provenance = {
            "path": str(streaming_path),
            "sha256": _sha256(streaming_path),
            "status": streaming_preflight.get("status"),
        }
    streaming_handoff_verified = bool(
        streaming_preflight is not None
        and streaming_preflight.get("status")
        == "accepted_coordinate_consistent_streaming_handoff_preflight"
        and streaming_preflight.get("all_gates_pass") is True
        and streaming_preflight.get("configuration", {}).get(
            "observation_value_fields_accessed"
        ) == []
        and streaming_preflight.get("gates", {}).get(
            "durable_hdf5_reopened_with_expected_frames"
        ) is True
        and streaming_preflight.get("gates", {}).get(
            "all_five_fluxes_conserved_by_area_remap"
        ) is True
    )
    completed_inner_transient = None
    completed_inner_transient_provenance = None
    if completed_inner_transient_path is not None:
        inner_path = Path(completed_inner_transient_path).resolve()
        completed_inner_transient = json.loads(
            inner_path.read_text(encoding="utf-8")
        )
        if completed_inner_transient.get("schema") != (
            "degali.thorney-ti021-inner-transient.v1"
        ):
            raise ValueError("unsupported TI021 inner-transient schema")
        completed_inner_transient_provenance = {
            "path": str(inner_path),
            "sha256": _sha256(inner_path),
            "status": completed_inner_transient.get("status"),
        }
    completed_inner_transient_verified = bool(
        completed_inner_transient is not None
        and completed_inner_transient.get("status")
        == "accepted_completed_inner_transient"
        and completed_inner_transient.get("all_gates_pass") is True
        and completed_inner_transient.get("configuration", {}).get(
            "duration_s"
        ) == 600.0
        and completed_inner_transient.get("configuration", {}).get(
            "inner_receptor_count"
        ) == 60
        and completed_inner_transient.get("configuration", {}).get(
            "initial_binary_cloud_momentum"
        ) == "quiescent"
        and completed_inner_transient.get("configuration", {}).get(
            "observation_values_used_to_configure_run"
        ) == []
    )
    coarse_inner_transient = None
    coarse_inner_transient_provenance = None
    if coarse_inner_transient_path is not None:
        coarse_path = Path(coarse_inner_transient_path).resolve()
        coarse_inner_transient = json.loads(coarse_path.read_text(encoding="utf-8"))
        if coarse_inner_transient.get("schema") != (
            "degali.thorney-ti021-inner-transient.v1"
        ):
            raise ValueError("unsupported TI021 coarse inner-transient schema")
        coarse_inner_transient_provenance = {
            "path": str(coarse_path),
            "sha256": _sha256(coarse_path),
            "status": coarse_inner_transient.get("status"),
        }
    explicit_coarse_inner_transient_verified = bool(
        coarse_inner_transient is not None
        and coarse_inner_transient.get("status")
        == "accepted_completed_inner_transient"
        and coarse_inner_transient.get("all_gates_pass") is True
        and coarse_inner_transient.get("configuration", {}).get("duration_s")
        == 600.0
        and coarse_inner_transient.get("configuration", {}).get(
            "cell_size_m", 5.0
        ) == 5.0
        and coarse_inner_transient.get("configuration", {}).get(
            "inner_receptor_count"
        ) == 60
        and coarse_inner_transient.get("configuration", {}).get(
            "initial_binary_cloud_momentum"
        ) == "quiescent"
        and coarse_inner_transient.get("configuration", {}).get(
            "observation_values_used_to_configure_run"
        ) == []
    )
    completed_inner_cell_size_m = (
        None if completed_inner_transient is None else float(
            completed_inner_transient.get("configuration", {}).get(
                "cell_size_m", 5.0
            )
        )
    )
    completed_dx5_inner_transient_verified = bool(
        explicit_coarse_inner_transient_verified or (
            completed_inner_transient_verified
            and math.isclose(completed_inner_cell_size_m, 5.0)
        )
    )
    completed_dx2p5_inner_transient_verified = bool(
        completed_inner_transient_verified
        and math.isclose(completed_inner_cell_size_m, 2.5)
    )
    any_completed_inner_transient_verified = bool(
        completed_inner_transient_verified
        or explicit_coarse_inner_transient_verified
    )
    inner_transient_score = None
    inner_transient_score_provenance = None
    if inner_transient_score_path is not None:
        score_path = Path(inner_transient_score_path).resolve()
        inner_transient_score = json.loads(
            score_path.read_text(encoding="utf-8")
        )
        if inner_transient_score.get("schema") != (
            "degali.thorney-ti021-inner-transient-score.v1"
        ):
            raise ValueError("unsupported TI021 inner-transient score schema")
        inner_transient_score_provenance = {
            "path": str(score_path),
            "sha256": _sha256(score_path),
            "status": inner_transient_score.get("status"),
            "decision": inner_transient_score.get("decision", {}).get(
                "result"
            ),
        }
    coarse_inner_score_available = bool(
        inner_transient_score is not None
        and inner_transient_score.get("status")
        == "scored_no_fit_inner_dx5_not_grid_converged"
        and inner_transient_score.get("protocol", {}).get(
            "configuration_accessed_observation_values"
        ) == []
        and inner_transient_score.get("resolution_limit", {}).get(
            "promotion_allowed"
        ) is False
        and inner_transient_score.get("decision", {}).get("result")
        == "withheld_as_quantitative_validation"
    )
    coarse_inner_score = None
    coarse_inner_score_provenance = None
    if coarse_inner_score_path is not None:
        coarse_score_path = Path(coarse_inner_score_path).resolve()
        coarse_inner_score = json.loads(
            coarse_score_path.read_text(encoding="utf-8")
        )
        if coarse_inner_score.get("schema") != (
            "degali.thorney-ti021-inner-transient-score.v1"
        ):
            raise ValueError("unsupported TI021 coarse inner-score schema")
        coarse_inner_score_provenance = {
            "path": str(coarse_score_path),
            "sha256": _sha256(coarse_score_path),
            "status": coarse_inner_score.get("status"),
            "decision": coarse_inner_score.get("decision", {}).get("result"),
        }
    explicit_coarse_inner_score_available = bool(
        coarse_inner_score is not None
        and coarse_inner_score.get("status")
        == "scored_no_fit_inner_dx5_not_grid_converged"
        and coarse_inner_score.get("protocol", {}).get(
            "configuration_accessed_observation_values"
        ) == []
        and coarse_inner_score.get("resolution_limit", {}).get(
            "promotion_allowed"
        ) is False
        and coarse_inner_score.get("decision", {}).get("result")
        == "withheld_as_quantitative_validation"
    )
    coarse_inner_score_available = bool(
        coarse_inner_score_available or explicit_coarse_inner_score_available
    )
    refined_inner_score_available = bool(
        inner_transient_score is not None
        and inner_transient_score.get("status")
        == "scored_no_fit_inner_dx2p5_pending_finest_grid"
        and inner_transient_score.get("protocol", {}).get(
            "configuration_accessed_observation_values"
        ) == []
        and inner_transient_score.get("resolution_limit", {}).get(
            "cell_size_m"
        ) == 2.5
        and inner_transient_score.get("resolution_limit", {}).get(
            "promotion_allowed"
        ) is False
        and inner_transient_score.get("decision", {}).get("result")
        == "withheld_as_quantitative_validation"
    )
    inner_score_available = bool(
        coarse_inner_score_available or refined_inner_score_available
    )
    corrector_tolerance_audit = None
    corrector_tolerance_audit_provenance = None
    if corrector_tolerance_audit_path is not None:
        tolerance_path = Path(corrector_tolerance_audit_path).resolve()
        corrector_tolerance_audit = json.loads(
            tolerance_path.read_text(encoding="utf-8")
        )
        if corrector_tolerance_audit.get("schema") != (
            "degali.thorney-ti021-corrector-tolerance-audit.v1"
        ):
            raise ValueError("unsupported TI021 corrector-tolerance audit schema")
        corrector_tolerance_audit_provenance = {
            "path": str(tolerance_path),
            "sha256": _sha256(tolerance_path),
            "status": corrector_tolerance_audit.get("status"),
        }
    dx2p5_corrector_tolerance_verified = bool(
        corrector_tolerance_audit is not None
        and corrector_tolerance_audit.get("status")
        == "accepted_1e4_corrector_tolerance_for_dx2p5_production_with_postcheck"
        and corrector_tolerance_audit.get("all_gates_pass") is True
        and corrector_tolerance_audit.get("configuration", {}).get(
            "cell_size_m"
        ) == 2.5
        and corrector_tolerance_audit.get("configuration", {}).get(
            "reference_tolerance"
        ) == 1.0e-6
        and corrector_tolerance_audit.get("configuration", {}).get(
            "candidate_tolerance"
        ) == 1.0e-4
        and corrector_tolerance_audit.get("configuration", {}).get(
            "observation_value_fields_accessed"
        ) == []
        and corrector_tolerance_audit.get("decision", {}).get(
            "mandatory_postcheck"
        )
    )
    strict_peak_window_audit = None
    strict_peak_window_audit_provenance = None
    if strict_peak_window_audit_path is not None:
        strict_path = Path(strict_peak_window_audit_path).resolve()
        strict_peak_window_audit = json.loads(
            strict_path.read_text(encoding="utf-8")
        )
        if strict_peak_window_audit.get("schema") != (
            "degali.thorney-ti021-strict-peak-window-audit.v1"
        ):
            raise ValueError("unsupported TI021 strict peak-window audit schema")
        strict_peak_window_audit_provenance = {
            "path": str(strict_path),
            "sha256": _sha256(strict_path),
            "status": strict_peak_window_audit.get("status"),
        }
    strict_peak_window_verified = bool(
        strict_peak_window_audit is not None
        and strict_peak_window_audit.get("status")
        == "accepted_1e4_corrector_tolerance_at_model_peak"
        and strict_peak_window_audit.get("configuration", {}).get(
            "cell_size_m"
        ) == 2.5
        and strict_peak_window_audit.get("configuration", {}).get(
            "candidate_tolerance"
        ) == 1.0e-4
        and strict_peak_window_audit.get("configuration", {}).get(
            "strict_tolerance"
        ) == 1.0e-6
        and strict_peak_window_audit.get("selection", {}).get(
            "observation_value_fields_accessed"
        ) == []
        and all(strict_peak_window_audit.get("gates", {}).values())
    )
    full_event_grid_comparison = None
    full_event_grid_comparison_provenance = None
    if full_event_grid_comparison_path is not None:
        comparison_path = Path(full_event_grid_comparison_path).resolve()
        full_event_grid_comparison = json.loads(
            comparison_path.read_text(encoding="utf-8")
        )
        if full_event_grid_comparison.get("schema") != (
            "degali.thorney-ti021-inner-grid-comparison.v1"
        ):
            raise ValueError("unsupported TI021 inner-grid comparison schema")
        full_event_grid_comparison_provenance = {
            "path": str(comparison_path),
            "sha256": _sha256(comparison_path),
            "status": full_event_grid_comparison.get("status"),
        }
    full_event_two_grid_diagnostic_available = bool(
        full_event_grid_comparison is not None
        and full_event_grid_comparison.get("status")
        == "diagnostic_two_grid_comparison_pending_1p25m_window"
        and full_event_grid_comparison.get("grid_convergence_promoted") is False
        and full_event_grid_comparison.get("configuration", {}).get(
            "coarse_cell_size_m"
        ) == 5.0
        and full_event_grid_comparison.get("configuration", {}).get(
            "fine_cell_size_m"
        ) == 2.5
        and full_event_grid_comparison.get("configuration", {}).get(
            "duration_s"
        ) == 600.0
        and full_event_grid_comparison.get("provenance", {}).get(
            "observation_value_fields_accessed"
        ) == []
    )
    report_facts = {
        "initial_cloud_volume_m3": 2050.0,
        "initial_relative_density_ratio": 2.02,
        "container_geometry": {
            "shape": "twelve-sided flexible container",
            "across_m": 14.0,
            "nominal_full_height_m": 14.0,
            "nominal_2000_m3_fill_height_m": 13.0,
        },
        "wall_geometry": {
            "shape": "semi-circular impermeable tarpaulin fence",
            "centre_site_m": [400.0, 200.0],
            "radius_m": 50.0,
            "height_m": 5.0,
            "arc_start_math_deg": 0.0,
            "arc_end_math_deg": 180.0,
            "thickness_m": None,
        },
        "event_mean_wind_speed_m_s": 3.9,
        "event_mean_heading_relative_to_array_deg": -6.1,
        "sensor_averaging_time_s": 0.6,
        "nominal_sensor_lower_resolution_percent": 0.1,
        "arrival_operator": {
            "definition": "first concentration reaching 5 percent of the local peak",
            "departure_definition": "last concentration exceeding 5 percent of the local peak",
            "time_reading_resolution_s": 5.0,
            "concentration_time_series_available": False,
            "local_peak_concentration_available_for_sensor_rows": False,
            "source": "CPP Project 88-0470 / OSTI 6934392, Thorney Island data reduction",
        },
    }
    geometry_conflict = not math.isclose(
        dataset.workbook_wall_height_m,
        float(report_facts["wall_geometry"]["height_m"]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    wind_speed_conflict = not math.isclose(
        dataset.workbook_wind_speed_m_s,
        float(report_facts["event_mean_wind_speed_m_s"]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    workbook_heading_relative_to_array = (
        (dataset.workbook_wind_from_north_deg + 180.0) % 360.0
        - dataset.site_y_axis_from_north_deg
    )
    workbook_heading_relative_to_array = (
        (workbook_heading_relative_to_array + 180.0) % 360.0 - 180.0
    )
    heading_conflict = not math.isclose(
        workbook_heading_relative_to_array,
        float(report_facts["event_mean_heading_relative_to_array_deg"]),
        rel_tol=0.0,
        abs_tol=0.1,
    )
    conflicts = []
    if geometry_conflict:
        conflicts.append("wall_height_report_5m_vs_workbook_6m")
    if wind_speed_conflict:
        conflicts.append("event_wind_report_3p9mps_vs_workbook_3p4mps")
    if heading_conflict:
        conflicts.append(
            "wind_heading_report_minus6p1deg_vs_workbook_minus15p8deg"
        )
    conflicts.extend((
        "wall_thickness_not_reported",
        "workbook_porosity_field_is_1_for_a_reported_impermeable_fence",
        "full_sensor_concentration_time_series_not_in_workbook",
        "sensor_arrival_peak_relative_operator_not_recoverable_from_workbook",
    ))

    return {
        "schema": THORNEY_TI021_EVIDENCE_SCHEMA,
        "dataset_id": "thorney-island-phase2-trial-021",
        "status": (
            "conditional_pending_outer_transient_and_grid_convergence"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and conservative_nested_handoff_verified
                and streaming_handoff_verified
                and any_completed_inner_transient_verified
                and inner_score_available
            ) else (
            "conditional_pending_full_transient_convergence_and_no_fit_scoring"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and conservative_nested_handoff_verified
                and streaming_handoff_verified
            ) else (
            "conditional_pending_streaming_handoff_and_full_transient"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and conservative_nested_handoff_verified
            ) else (
            "conditional_pending_outer_momentum_overlap_and_full_transient"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and handoff_remap_verified
                and outer_binary_boundary_verified
            ) else (
            "conditional_pending_outer_boundary_coupling_and_full_transient"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and handoff_remap_verified
            ) else (
            "conditional_pending_conservative_nested_handoff_and_full_transient"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
            ) else (
            "conditional_pending_full_transient_convergence_and_no_fit_run"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
            ) else (
            "conditional_pending_turbulent_closure_and_no_fit_run"
            if (
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
            ) else (
                "conditional_pending_full_scale_no_fit_preflight"
                if (
                    geometry_verified and initial_cloud_verified
                    and molecular_properties_verified
                ) else (
                    "conditional_pending_property_evidence_and_full_scale_no_fit_run"
                    if geometry_verified and initial_cloud_verified else (
                        "withheld_pending_r12_full_mixture_and_initial_cloud_solver"
                        if geometry_verified
                        else "withheld_pending_curved_wall_geometry_verification"
                    )
                )
            ))))))))
        ),
        "provenance": {
            "workbook_path": str(dataset.workbook_path),
            "workbook_sha256": dataset.workbook_sha256,
            "report_path": str(report),
            "report_sha256": _sha256(report),
            "curved_wall_geometry_audit": geometry_audit_provenance,
            "initial_cloud_audit": initial_cloud_audit_provenance,
            "r12_property_audit": property_audit_provenance,
            "full_scale_preflight": full_scale_preflight_provenance,
            "turbulent_mixing_policy": turbulent_mixing_policy_provenance,
            "full_transient_plan": full_transient_plan_provenance,
            "nested_handoff_audit": nested_handoff_audit_provenance,
            "outer_inlet_coupling_audit": outer_inlet_audit_provenance,
            "outer_momentum_overlap_audit": (
                outer_momentum_audit_provenance
            ),
            "handoff_streaming_preflight": (
                streaming_preflight_provenance
            ),
            "coarse_inner_transient": coarse_inner_transient_provenance,
            "coarse_inner_score": coarse_inner_score_provenance,
            "completed_inner_transient": (
                completed_inner_transient_provenance
            ),
            "inner_transient_score": inner_transient_score_provenance,
            "corrector_tolerance_audit": (
                corrector_tolerance_audit_provenance
            ),
            "strict_peak_window_audit": strict_peak_window_audit_provenance,
            "full_event_grid_comparison": (
                full_event_grid_comparison_provenance
            ),
        },
        "scope": {
            "admissible": dataset.validation_scope,
            "hydrogen_validation": False,
            "dynamic_waveform_validation": False,
            "obstacle_component_validation_after_freeze": True,
        },
        "source": {
            "interpretation": "instantaneous finite-volume initial cloud",
            "continuous_point_source_interpretation_allowed": False,
            "release_origin_site_m": list(dataset.release_origin_site_m),
            "substance": dataset.substance,
            "initial_r12_mole_fraction": dataset.initial_r12_mole_fraction,
            "effective_molecular_weight_kg_mol": (
                dataset.effective_molecular_weight_kg_mol
            ),
            "initial_temperature_k": dataset.initial_temperature_k,
            "initial_pressure_pa": dataset.initial_pressure_pa,
            "workbook_mass_proxy_kg": dataset.workbook_mass_proxy_kg,
            "workbook_mass_proxy_duration_s": (
                dataset.workbook_mass_proxy_duration_s
            ),
            "workbook_release_rate_kg_s": dataset.workbook_release_rate_kg_s,
            "workbook_release_duration_s": dataset.workbook_release_duration_s,
            "workbook_release_mass_proxy_kg": dataset.workbook_release_mass_proxy_kg,
            "workbook_release_mass_proxy_is_validated_source": False,
            "report_initial_cloud_volume_m3": (
                report_facts["initial_cloud_volume_m3"]
            ),
            "report_initial_relative_density_ratio": (
                report_facts["initial_relative_density_ratio"]
            ),
            "container_geometry": report_facts["container_geometry"],
        },
        "meteorology": {
            "workbook_wind_speed_m_s": dataset.workbook_wind_speed_m_s,
            "report_event_mean_wind_speed_m_s": (
                report_facts["event_mean_wind_speed_m_s"]
            ),
            "reference_height_m": dataset.workbook_wind_reference_height_m,
            "friction_velocity_m_s": dataset.friction_velocity_m_s,
            "roughness_length_m": dataset.roughness_length_m,
            "stability_class": dataset.stability_class,
            "workbook_wind_from_north_deg": dataset.workbook_wind_from_north_deg,
            "ideal_wind_from_north_deg": dataset.ideal_wind_from_north_deg,
            "site_y_axis_from_north_deg": dataset.site_y_axis_from_north_deg,
            "workbook_heading_relative_to_array_deg": (
                workbook_heading_relative_to_array
            ),
            "report_event_mean_heading_relative_to_array_deg": (
                report_facts["event_mean_heading_relative_to_array_deg"]
            ),
        },
        "obstacle": {
            "report_geometry": report_facts["wall_geometry"],
            "workbook_radius_m": dataset.workbook_wall_radius_m,
            "workbook_center_site_m": list(dataset.workbook_wall_center_site_m),
            "workbook_height_m": dataset.workbook_wall_height_m,
            "workbook_porosity_field": dataset.workbook_porosity_field,
        },
        "observations": {
            "dose_sensor_count": len(dataset.dose_sensors),
            "positive_dose_sensor_count": dataset.positive_dose_sensor_count,
            "arc_observation_count": len(dataset.arc_observations),
            "reported_sensor_averaging_time_s": (
                report_facts["sensor_averaging_time_s"]
            ),
            "nominal_sensor_lower_resolution_percent": (
                report_facts["nominal_sensor_lower_resolution_percent"]
            ),
            "arrival_operator": report_facts["arrival_operator"],
            "dose_sensors": [asdict(sensor) for sensor in dataset.dose_sensors],
            "arc_observations": [
                asdict(observation) for observation in dataset.arc_observations
            ],
        },
        "conflicts_and_missing_evidence": conflicts,
        "frozen_source_of_truth_policy": {
            "basis": (
                "use the final HSE report for physical obstacle design and "
                "event-summary meteorology; retain SMEDIS values as declared "
                "sensitivity alternatives; no concentration value selects the policy"
            ),
            "primary_wall_height_m": 5.0,
            "wall_height_sensitivity_m": 6.0,
            "primary_wind_speed_m_s": 3.9,
            "wind_speed_sensitivity_m_s": 3.4,
            "primary_heading_relative_to_array_deg": -6.1,
            "heading_sensitivity_relative_to_array_deg": (
                workbook_heading_relative_to_array
            ),
            "primary_wall_solidity": "impermeable",
            "workbook_porosity_field_treatment": (
                "not interpreted because the final report directly states impermeable"
            ),
            "wall_thickness_policy": (
                "one-cell numerical surface regularisation that shrinks with grid"
            ),
        },
        "promotion_gates": {
            "source_semantics_frozen": True,
            "site_coordinate_transform_recovered": True,
            "curved_wall_planform_recovered": True,
            "wall_height_conflict_resolved": True,
            "wall_thickness_resolved": geometry_verified,
            "meteorological_operator_resolved": True,
            "curved_wall_solver_representation_verified": geometry_verified,
            "full_binary_r12_air_mixture_solver_available": (
                initial_cloud_verified
            ),
            "finite_volume_initial_cloud_state_supported": (
                initial_cloud_verified
            ),
            "r12_molecular_property_inputs_frozen": (
                molecular_properties_verified
            ),
            "full_scale_single_step_preflight_passed": (
                full_scale_preflight_verified
            ),
            "turbulent_scalar_mixing_closure_frozen": (
                turbulent_mixing_verified
            ),
            "sensor_complete_execution_plan_frozen": (
                full_transient_plan_verified
            ),
            "conservative_handoff_remap_operator_verified": (
                handoff_remap_verified
            ),
            "outer_binary_species_energy_continuity_boundary_verified": (
                outer_binary_boundary_verified
            ),
            "outer_momentum_and_overlap_verified": (
                outer_momentum_overlap_verified
            ),
            "conservative_nested_handoff_verified": (
                conservative_nested_handoff_verified
            ),
            "coordinate_consistent_streaming_handoff_verified": (
                streaming_handoff_verified
            ),
            "primary_dx5_inner_600s_transient_completed": (
                completed_dx5_inner_transient_verified
            ),
            "refined_dx2p5_inner_600s_transient_completed": (
                completed_dx2p5_inner_transient_verified
            ),
            "coarse_dx5_inner_no_fit_score_available": (
                coarse_inner_score_available
            ),
            "refined_dx2p5_inner_no_fit_score_available": (
                refined_inner_score_available
            ),
            "coarse_dx5_inner_quantitative_validation_promoted": False,
            "dx2p5_corrector_tolerance_verified_with_peak_postcheck_pending": (
                dx2p5_corrector_tolerance_verified
            ),
            "dx2p5_strict_peak_window_postcheck_passed": (
                strict_peak_window_verified
            ),
            "full_event_5m_to_2p5m_grid_diagnostic_available": (
                full_event_two_grid_diagnostic_available
            ),
            "full_event_grid_convergence_promoted": False,
            "independent_evaluation_observations_available": True,
            "eligible_for_full_scale_no_fit_preflight": bool(
                geometry_verified and initial_cloud_verified
            ),
            "eligible_for_inner_arc_no_fit_run": bool(
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
            ),
            "eligible_for_full_145_sensor_no_fit_run": bool(
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and conservative_nested_handoff_verified
                and streaming_handoff_verified
            ),
            "eligible_for_no_fit_model_run": bool(
                geometry_verified and initial_cloud_verified
                and molecular_properties_verified
                and full_scale_preflight_verified
                and turbulent_mixing_verified
                and full_transient_plan_verified
                and conservative_nested_handoff_verified
                and streaming_handoff_verified
            ),
        },
        "decision": {
            "obstacle_dense_gas_component": (
                "conditional" if geometry_verified and initial_cloud_verified
                else "withheld"
            ),
            "hydrogen_accuracy": "not_applicable",
            "reason": (
                (
                    "curved-wall geometry and a conservative full-mixture finite-"
                    "volume initial-cloud path now pass their component gates; "
                    + (
                        (
                            "molecular R12 properties and the report-scale "
                            "single-step coupled preflight are also accepted; "
                            + (
                                "the concentration-blind mixing envelope is "
                                "frozen; " + (
                                    "a sensor-complete nested execution plan "
                                    "shows that " + (
                                        "the conservative plane-remap operator "
                                        "passes; " + (
                                            (
                                                "the time-interpolated outer binary "
                                                "species/energy/continuity boundary, "
                                                "three-component momentum coupling "
                                                "and same-grid overlap comparison all "
                                                "pass; " + (
                                                    "the actual-inner coordinate-"
                                                    "consistent HDF5 streaming handoff "
                                                    "also passes; " + (
                                                        "the 600 s dx=5 m inner run and "
                                                        "no-fit score are complete but "
                                                        "quantitative promotion is "
                                                        "withheld; the 2.5/1.25 m grid "
                                                        "checks and outer 145-sensor "
                                                        "transient remain"
                                                        if (
                                                            completed_inner_transient_verified
                                                            and coarse_inner_score_available
                                                        ) else
                                                        "the 600 s three-grid and 145-"
                                                        "sensor no-fit runs remain"
                                                    )
                                                    if streaming_handoff_verified else
                                                    "the actual-inner coordinate-"
                                                    "consistent streaming handoff "
                                                    "preflight remains before the 600 s "
                                                    "three-grid and 145-sensor no-fit runs"
                                                )
                                                if outer_momentum_overlap_verified else
                                                "the time-interpolated outer binary "
                                                "species/energy/continuity boundary "
                                                "also passes, but momentum coupling "
                                                "and overlap validation remain"
                                            )
                                            if outer_binary_boundary_verified else
                                            "outer transient boundary coupling "
                                            "and overlap validation remain before "
                                            "the full 145-sensor no-fit run"
                                        )
                                        if handoff_remap_verified else
                                        "conservative inner-to-outer handoff must "
                                        "be verified before the full 145-sensor "
                                        "no-fit run"
                                    )
                                    if full_transient_plan_verified else
                                    "a sensor-complete execution plan, the six-"
                                    "member full-transient grid/time/domain audit "
                                    "and no-fit sensor score are next"
                                )
                                if turbulent_mixing_verified else
                                "turbulent mixing, full-transient grid/time/"
                                "domain convergence and the no-fit sensor score "
                                "remain required"
                            )
                            if full_scale_preflight_verified else
                            "molecular R12 properties are also frozen, but "
                            "the report-scale coupled preflight, turbulent "
                            "mixing, full-transient grid/time/domain "
                            "convergence and the no-fit sensor score remain "
                            "required"
                        )
                        if molecular_properties_verified else
                        "R12 transport-property evidence, full-scale grid/time/"
                        "domain preflight and the no-fit sensor score remain required"
                    )
                    if geometry_verified and initial_cloud_verified else
                    "TI021 requires both the curved-wall geometry audit and the "
                    "conservative full-mixture finite-volume initial-cloud audit"
                )
            ),
        },
    }


__all__ = [
    "THORNEY_TI021_EVIDENCE_SCHEMA", "Ti021ArcObservation",
    "Ti021DoseSensor", "ThorneyTi021Dataset", "read_thorney_ti021_workbook",
    "thorney_ti021_evidence_record",
]
