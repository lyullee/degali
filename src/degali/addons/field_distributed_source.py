"""Explicit global atmospheric-vapour sources for the local field plane."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

from .field_contracts import BoundedValue
from .semi_fv_obstacle import SourceRateSchedule


_FIELD_DISTRIBUTED_SOURCE_KINDS = frozenset({
    "declared_atmospheric_vapour",
    "post_flash_atmospheric_vapour",
    "pool_vapour",
    "droplet_evaporation",
    "in_flight_droplet_evaporation",
})


@dataclass(frozen=True)
class FieldDistributedVapourSource:
    """One already-atmospheric H2 source with a global location and schedule.

    ``source_kind`` documents the upstream physical ledger but does not select
    a flash, pool or droplet model. In particular, callers may use
    ``in_flight_droplet_evaporation`` only when a time/space-resolved droplet
    ledger has been supplied. Optional position, vertical-width, and common-
    clock lower/upper rate bounds are deterministic corners; they never imply
    a probability density or a lateral-wake closure.
    """

    label: str
    position_m: tuple[float, float, float]
    vertical_sigma_m: float
    schedule: SourceRateSchedule
    evidence_id: str
    source_kind: str = "declared_atmospheric_vapour"
    warnings: tuple[str, ...] = ()
    position_uncertainty_m: tuple[BoundedValue, BoundedValue, BoundedValue] | None = None
    vertical_sigma_uncertainty: BoundedValue | None = None
    lower_schedule: SourceRateSchedule | None = None
    upper_schedule: SourceRateSchedule | None = None

    def __post_init__(self) -> None:
        for name, value in {
            "label": self.label,
            "evidence_id": self.evidence_id,
            "source_kind": self.source_kind,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"field distributed vapour source requires a declared {name}")
        if not isinstance(self.position_m, tuple):
            raise TypeError("field distributed vapour source position_m must be a tuple")
        if len(self.position_m) != 3 or any(isinstance(value, bool) for value in self.position_m) or not all(math.isfinite(float(value)) for value in self.position_m):
            raise ValueError("field distributed vapour source position_m must contain three finite values")
        if self.position_m[2] < 0.0:
            raise ValueError("field distributed vapour source height cannot be negative")
        if isinstance(self.vertical_sigma_m, bool) or not math.isfinite(float(self.vertical_sigma_m)) or self.vertical_sigma_m <= 0.0:
            raise ValueError("field distributed vapour source vertical_sigma_m must be positive and finite")
        if self.position_uncertainty_m is not None:
            if not isinstance(self.position_uncertainty_m, tuple):
                raise TypeError("position_uncertainty_m must be a tuple")
            if len(self.position_uncertainty_m) != 3 or not all(
                isinstance(value, BoundedValue) for value in self.position_uncertainty_m
            ):
                raise TypeError(
                    "position_uncertainty_m must contain three BoundedValue values"
                )
            for index, (nominal, bound) in enumerate(
                zip(self.position_m, self.position_uncertainty_m)
            ):
                if bound.unit != "m":
                    raise ValueError(
                        f"position_uncertainty_m[{index}] must use unit 'm'"
                    )
                if (
                    not isinstance(bound.source, str)
                    or not bound.source.strip()
                    or bound.source.strip().lower() == "unspecified"
                ):
                    raise ValueError(
                        f"position_uncertainty_m[{index}] requires an explicit source"
                    )
                if not math.isclose(
                    bound.nominal, float(nominal), rel_tol=1.0e-9, abs_tol=1.0e-12
                ):
                    raise ValueError(
                        f"position_uncertainty_m[{index}].nominal must match position_m"
                    )
                if index == 2 and bound.lower < 0.0:
                    raise ValueError(
                        "position_uncertainty_m height bounds cannot be negative"
                    )
        if self.vertical_sigma_uncertainty is not None:
            bound = self.vertical_sigma_uncertainty
            if not isinstance(bound, BoundedValue):
                raise TypeError("vertical_sigma_uncertainty must be a BoundedValue or None")
            if bound.unit != "m":
                raise ValueError("vertical_sigma_uncertainty must use unit 'm'")
            if (
                not isinstance(bound.source, str)
                or not bound.source.strip()
                or bound.source.strip().lower() == "unspecified"
            ):
                raise ValueError("vertical_sigma_uncertainty requires an explicit source")
            if bound.lower <= 0.0:
                raise ValueError("vertical_sigma_uncertainty bounds must be positive")
            if not math.isclose(
                bound.nominal, float(self.vertical_sigma_m),
                rel_tol=1.0e-9, abs_tol=1.0e-12,
            ):
                raise ValueError(
                    "vertical_sigma_uncertainty.nominal must match vertical_sigma_m"
                )
        if not isinstance(self.schedule, SourceRateSchedule):
            raise TypeError("field distributed vapour source schedule must be a SourceRateSchedule")
        if (self.lower_schedule is None) != (self.upper_schedule is None):
            raise ValueError("distributed source lower_schedule and upper_schedule must be supplied together")
        if self.lower_schedule is not None and self.upper_schedule is not None:
            for name, corner in (
                ("lower_schedule", self.lower_schedule),
                ("upper_schedule", self.upper_schedule),
            ):
                if not isinstance(corner, SourceRateSchedule):
                    raise TypeError(f"distributed source {name} must be a SourceRateSchedule")
                if (
                    corner.source_id != self.schedule.source_id
                    or corner.time_s != self.schedule.time_s
                    or corner.rate_operator != self.schedule.rate_operator
                ):
                    raise ValueError(
                        f"distributed source {name} must use the nominal source ID, time axis, and rate operator"
                    )
            for lower, nominal, upper in zip(
                self.lower_schedule.rate_kg_s,
                self.schedule.rate_kg_s,
                self.upper_schedule.rate_kg_s,
            ):
                if not lower <= nominal <= upper:
                    raise ValueError(
                        "distributed source rate bounds must contain nominal rates"
                    )
            self.lower_schedule.require_zero_endpoint(
                "distributed source lower_schedule"
            )
            self.upper_schedule.require_zero_endpoint(
                "distributed source upper_schedule"
            )
        if self.source_kind not in _FIELD_DISTRIBUTED_SOURCE_KINDS:
            raise ValueError(
                "field distributed vapour source_kind must be one of "
                + ", ".join(sorted(_FIELD_DISTRIBUTED_SOURCE_KINDS))
            )
        if self.schedule.source_id.strip().lower() == "declared":
            raise ValueError(
                "field distributed vapour source schedule requires explicit source_id provenance"
            )
        if self.schedule.released_mass_kg <= 0.0:
            raise ValueError("field distributed vapour source schedule must release positive mass")
        self.schedule.require_zero_endpoint(
            "field distributed vapour source schedule"
        )
        if not isinstance(self.warnings, tuple):
            raise TypeError("field distributed vapour source warnings must be a tuple")
        if any(not isinstance(item, str) or not item.strip() for item in self.warnings):
            raise ValueError("field distributed vapour source warnings must contain non-empty strings")

    def uncertainty_fields(self) -> dict[str, BoundedValue]:
        """Return explicit geometry/width bounds with stable envelope keys."""
        values: dict[str, BoundedValue] = {}
        if self.position_uncertainty_m is not None:
            values.update({
                name: bound
                for name, bound in zip(
                    (
                        f"distributed_source.{self.label}.position_x_m",
                        f"distributed_source.{self.label}.position_y_m",
                        f"distributed_source.{self.label}.position_z_m",
                    ),
                    self.position_uncertainty_m,
                )
            })
        if self.vertical_sigma_uncertainty is not None:
            values[
                f"distributed_source.{self.label}.vertical_sigma_m"
            ] = self.vertical_sigma_uncertainty
        return values

    @property
    def has_schedule_uncertainty(self) -> bool:
        """Whether lower/nominal/upper source-rate histories were declared."""
        return self.lower_schedule is not None and self.upper_schedule is not None

    def schedule_corners(self) -> tuple[tuple[str, SourceRateSchedule], ...]:
        """Return deterministic source-rate corners on the same clock."""
        if self.lower_schedule is None or self.upper_schedule is None:
            return (("nominal", self.schedule),)
        return (
            ("lower", self.lower_schedule),
            ("nominal", self.schedule),
            ("upper", self.upper_schedule),
        )

    def at_schedule_corner(self, label: str) -> "FieldDistributedVapourSource":
        """Fix one declared source-rate corner before a transport run."""
        choices = dict(self.schedule_corners())
        if label not in choices:
            raise ValueError(
                f"unknown distributed source schedule corner {label!r}; "
                f"expected one of {tuple(choices)}"
            )
        return replace(
            self,
            schedule=choices[label],
            lower_schedule=None,
            upper_schedule=None,
        )

    def at_corner(self, values: dict[str, float]) -> "FieldDistributedVapourSource":
        """Fix every declared geometry/width bound at one deterministic corner."""
        fields = self.uncertainty_fields()
        if set(values) != set(fields):
            missing = sorted(set(fields) - set(values))
            extra = sorted(set(values) - set(fields))
            raise ValueError(
                "distributed source corner keys must match uncertainty fields; "
                f"missing={missing}, extra={extra}"
            )
        position = self.position_m
        if self.position_uncertainty_m is not None:
            position = (
                float(values[f"distributed_source.{self.label}.position_x_m"]),
                float(values[f"distributed_source.{self.label}.position_y_m"]),
                float(values[f"distributed_source.{self.label}.position_z_m"]),
            )
        sigma = self.vertical_sigma_m
        if self.vertical_sigma_uncertainty is not None:
            sigma = float(values[f"distributed_source.{self.label}.vertical_sigma_m"])
        return replace(
            self,
            position_m=position,
            vertical_sigma_m=sigma,
            position_uncertainty_m=None,
            vertical_sigma_uncertainty=None,
        )


__all__ = ["FieldDistributedVapourSource"]
