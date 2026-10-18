"""Conservative resolved-inner to outer-domain face-flux handoff.

The handoff transfers *flux densities* on a downstream y-z plane.  Exact
rectangular overlap areas map the resolved inner plane into a wider and/or
coarser outer plane.  Values outside the inner footprint remain zero anomalies
and are later combined with the independently declared ambient outer inflow.
No clipping, renormalisation or concentration-derived coefficient is applied.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
    BinaryIdealGasMixtureState,
)


CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA = (
    "degali.conservative-nested-handoff-3d.v1"
)


def _edges(values: np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if (
        result.ndim != 1 or result.size < 2
        or np.any(~np.isfinite(result))
        or np.any(np.diff(result) <= 0.0)
    ):
        raise ValueError(f"{name} must be finite strictly increasing cell edges")
    return result


def _overlap_matrix(target_edges: np.ndarray, source_edges: np.ndarray) -> np.ndarray:
    left = np.maximum(target_edges[:-1, None], source_edges[None, :-1])
    right = np.minimum(target_edges[1:, None], source_edges[None, 1:])
    return np.maximum(right - left, 0.0)


def conservative_overlap_remap_flux_density_2d(
    source_flux_density: np.ndarray, *,
    source_z_edges_m: np.ndarray,
    source_y_edges_m: np.ndarray,
    target_z_edges_m: np.ndarray,
    target_y_edges_m: np.ndarray,
) -> np.ndarray:
    """Area-conservatively remap one or more y-z face flux-density fields.

    The final two source dimensions must be ``(z, y)``.  Any leading dimensions
    (for example time or ensemble member) are preserved.
    """

    source_z = _edges(source_z_edges_m, "source_z_edges_m")
    source_y = _edges(source_y_edges_m, "source_y_edges_m")
    target_z = _edges(target_z_edges_m, "target_z_edges_m")
    target_y = _edges(target_y_edges_m, "target_y_edges_m")
    if (
        target_z[0] > source_z[0] or target_z[-1] < source_z[-1]
        or target_y[0] > source_y[0] or target_y[-1] < source_y[-1]
    ):
        raise ValueError("target handoff plane must completely cover source plane")
    flux = np.asarray(source_flux_density, dtype=float)
    if flux.shape[-2:] != (source_z.size - 1, source_y.size - 1):
        raise ValueError("source flux final dimensions must match source z-y cells")
    if np.any(~np.isfinite(flux)):
        raise ValueError("source flux density must be finite")
    overlap_z = _overlap_matrix(target_z, source_z)
    overlap_y = _overlap_matrix(target_y, source_y)
    integrated = np.einsum(
        "ai,...ij,bj->...ab", overlap_z, flux, overlap_y,
        optimize=True,
    )
    target_area = (
        np.diff(target_z)[:, None] * np.diff(target_y)[None, :]
    )
    return integrated / target_area


def _integral(field: np.ndarray, z_edges: np.ndarray, y_edges: np.ndarray) -> np.ndarray:
    area = np.diff(z_edges)[:, None] * np.diff(y_edges)[None, :]
    return np.sum(field * area, axis=(-2, -1))


def _relative_integral_residual(
    source: np.ndarray, target: np.ndarray,
    source_z: np.ndarray, source_y: np.ndarray,
    target_z: np.ndarray, target_y: np.ndarray,
) -> float:
    source_integral = _integral(source, source_z, source_y)
    target_integral = _integral(target, target_z, target_y)
    source_area = np.diff(source_z)[:, None] * np.diff(source_y)[None, :]
    scale = np.sum(np.abs(source) * source_area, axis=(-2, -1))
    residual = np.abs(target_integral - source_integral) / np.maximum(scale, 1.0e-30)
    return float(np.max(residual))


@dataclass(frozen=True)
class ConservativeNestedHandoff3DResult:
    schema: str
    time_s: np.ndarray
    target_z_edges_m: np.ndarray
    target_y_edges_m: np.ndarray
    total_mass_flux_anomaly_kg_m2_s: np.ndarray
    species_mass_flux_kg_m2_s: np.ndarray
    sensible_enthalpy_flux_w_m2: np.ndarray
    lateral_momentum_flux_anomaly_pa: np.ndarray
    vertical_momentum_flux_anomaly_pa: np.ndarray
    diagnostics: dict[str, object]
    advective_species_mass_flux_kg_m2_s: np.ndarray | None = None
    advective_sensible_enthalpy_flux_w_m2: np.ndarray | None = None


@dataclass(frozen=True)
class NestedOuterInletFrame3D:
    time_s: float
    total_mass_flux_kg_m2_s: np.ndarray
    species_mass_flux_kg_m2_s: np.ndarray
    sensible_enthalpy_flux_w_m2: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_v_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    thermodynamic_state: BinaryIdealGasMixtureState
    diagnostics: dict[str, object]


def _interpolate_history(
    time_s: np.ndarray, field: np.ndarray, query_time_s: float,
) -> np.ndarray:
    if query_time_s < time_s[0] or query_time_s > time_s[-1]:
        raise ValueError("outer inlet query time lies outside handoff history")
    upper = int(np.searchsorted(time_s, query_time_s, side="right"))
    if upper == 0:
        return field[0].copy()
    if upper >= time_s.size:
        return field[-1].copy()
    lower = upper - 1
    fraction = (
        (query_time_s - time_s[lower])
        / (time_s[upper] - time_s[lower])
    )
    return (1.0 - fraction) * field[lower] + fraction * field[upper]


def nested_outer_inlet_frame_at(
    handoff: ConservativeNestedHandoff3DResult, *,
    time_s: float,
    ambient_total_mass_flux_kg_m2_s: float | np.ndarray,
    closure: BinaryIdealGasMixtureClosure,
) -> NestedOuterInletFrame3D:
    """Combine ambient inflow and net fluxes with a primitive inlet state.

    Version-1 handoffs lack advective/diffusive decomposition; those retain the
    historical net-flux-ratio state reconstruction and are marked in diagnostics.
    """

    if not isinstance(handoff, ConservativeNestedHandoff3DResult):
        raise TypeError("handoff must be ConservativeNestedHandoff3DResult")
    if not isinstance(closure, BinaryIdealGasMixtureClosure):
        raise TypeError("closure must be BinaryIdealGasMixtureClosure")
    query_time = float(time_s)
    if not np.isfinite(query_time):
        raise ValueError("time_s must be finite")
    shape = handoff.species_mass_flux_kg_m2_s.shape[1:]
    ambient = np.broadcast_to(
        np.asarray(ambient_total_mass_flux_kg_m2_s, dtype=float), shape,
    ).copy()
    if np.any(~np.isfinite(ambient)) or np.any(ambient <= 0.0):
        raise ValueError("ambient total mass flux must be positive and finite")
    mass_anomaly = _interpolate_history(
        handoff.time_s, handoff.total_mass_flux_anomaly_kg_m2_s, query_time,
    )
    species_flux = _interpolate_history(
        handoff.time_s, handoff.species_mass_flux_kg_m2_s, query_time,
    )
    enthalpy_flux = _interpolate_history(
        handoff.time_s, handoff.sensible_enthalpy_flux_w_m2, query_time,
    )
    lateral_momentum = _interpolate_history(
        handoff.time_s, handoff.lateral_momentum_flux_anomaly_pa, query_time,
    )
    vertical_momentum = _interpolate_history(
        handoff.time_s, handoff.vertical_momentum_flux_anomaly_pa, query_time,
    )
    total_mass_flux = ambient + mass_anomaly
    if np.any(total_mass_flux <= 0.0):
        raise ValueError(
            "ambient plus handoff anomaly must remain a unidirectional inlet"
        )
    # Independent conservative remapping/interpolation of the total and H2
    # fluxes can leave round-off-sized violations (observed at ~1e-21 kg/m2/s)
    # when the physical H2 flux is zero. Repair only within a local floating
    # point error bound; materially invalid inputs still fail closed.
    admissibility_tolerance = 64.0 * np.finfo(float).eps * np.maximum(
        np.maximum(np.abs(total_mass_flux), np.abs(species_flux)), 1.0e-30,
    )
    lower_roundoff = (species_flux < 0.0) & (
        species_flux >= -admissibility_tolerance
    )
    upper_roundoff = (species_flux > total_mass_flux) & (
        species_flux <= total_mass_flux + admissibility_tolerance
    )
    roundoff_correction = np.zeros_like(species_flux)
    roundoff_correction[lower_roundoff] = -species_flux[lower_roundoff]
    roundoff_correction[upper_roundoff] = (
        total_mass_flux[upper_roundoff] - species_flux[upper_roundoff]
    )
    if np.any(lower_roundoff) or np.any(upper_roundoff):
        species_flux = species_flux + roundoff_correction
    if np.any(species_flux < 0.0) or np.any(species_flux > total_mass_flux):
        raise ValueError(
            "outer inlet species flux must lie between zero and total mass flux"
        )
    # Net species/energy fluxes include molecular diffusion/conduction. Keep
    # those exact fluxes for conservation, but do not mistake their ratios to
    # total mass flux for the primitive thermodynamic state.
    advective_species_flux = handoff.advective_species_mass_flux_kg_m2_s
    advective_enthalpy_flux = handoff.advective_sensible_enthalpy_flux_w_m2
    legacy_state_from_net_flux = (
        advective_species_flux is None or advective_enthalpy_flux is None
    )
    if legacy_state_from_net_flux:
        mass_fraction = species_flux / total_mass_flux
        specific_enthalpy = enthalpy_flux / total_mass_flux
    else:
        mass_fraction = _interpolate_history(
            handoff.time_s, advective_species_flux, query_time,
        ) / total_mass_flux
        specific_enthalpy = _interpolate_history(
            handoff.time_s, advective_enthalpy_flux, query_time,
        ) / total_mass_flux
        if np.any(mass_fraction < 0.0) or np.any(mass_fraction > 1.0):
            raise ValueError("advective inlet species fraction must lie in [0, 1]")
    state = closure.state_from_mass_fraction_and_specific_sensible_enthalpy(
        mass_fraction, specific_enthalpy,
    )
    velocity_u = total_mass_flux / state.mixture_density_kg_m3
    velocity_v = lateral_momentum / total_mass_flux
    velocity_w = vertical_momentum / total_mass_flux
    return NestedOuterInletFrame3D(
        time_s=query_time,
        total_mass_flux_kg_m2_s=total_mass_flux,
        species_mass_flux_kg_m2_s=species_flux,
        sensible_enthalpy_flux_w_m2=enthalpy_flux,
        velocity_u_m_s=velocity_u,
        velocity_v_m_s=velocity_v,
        velocity_w_m_s=velocity_w,
        thermodynamic_state=state,
        diagnostics={
            "temporal_interpolation": "linear",
            "ambient_plus_anomaly": True,
            "minimum_total_mass_flux_kg_m2_s": float(np.min(total_mass_flux)),
            "maximum_total_mass_flux_kg_m2_s": float(np.max(total_mass_flux)),
            "minimum_species_mass_fraction": float(np.min(mass_fraction)),
            "maximum_species_mass_fraction": float(np.max(mass_fraction)),
            "minimum_temperature_k": float(np.min(state.temperature_k)),
            "maximum_temperature_k": float(np.max(state.temperature_k)),
            "unidirectional_inlet": True,
            "thermodynamic_state_source": (
                "legacy_net_flux_ratio" if legacy_state_from_net_flux
                else "separate_advective_species_and_enthalpy_fluxes"
            ),
            # Stable public gate name used by audit consumers.  Keep the more
            # explicit alias below for serialized reports produced by newer
            # callers.
            "clipping_or_posthoc_renormalisation": False,
            "physical_flux_clipping_or_renormalisation": False,
            "roundoff_admissibility_correction_cell_count": int(
                np.count_nonzero(roundoff_correction)
            ),
            "roundoff_admissibility_correction_l1_kg_m2_s": float(
                np.sum(np.abs(roundoff_correction))
            ),
        },
    )


def remap_conservative_nested_handoff_3d(
    *,
    time_s: np.ndarray,
    source_z_edges_m: np.ndarray,
    source_y_edges_m: np.ndarray,
    target_z_edges_m: np.ndarray,
    target_y_edges_m: np.ndarray,
    total_mass_flux_anomaly_kg_m2_s: np.ndarray,
    species_mass_flux_kg_m2_s: np.ndarray,
    sensible_enthalpy_flux_w_m2: np.ndarray,
    lateral_momentum_flux_anomaly_pa: np.ndarray,
    vertical_momentum_flux_anomaly_pa: np.ndarray,
    source_id: str,
    advective_species_mass_flux_kg_m2_s: np.ndarray | None = None,
    advective_sensible_enthalpy_flux_w_m2: np.ndarray | None = None,
) -> ConservativeNestedHandoff3DResult:
    """Remap a time history of conserved downstream-plane flux densities."""

    times = np.asarray(time_s, dtype=float)
    if (
        times.ndim != 1 or times.size < 1 or np.any(~np.isfinite(times))
        or np.any(np.diff(times) <= 0.0)
    ):
        raise ValueError("time_s must be finite and strictly increasing")
    if not isinstance(source_id, str) or not source_id.strip():
        raise ValueError("source_id must be a non-empty string")
    source_z = _edges(source_z_edges_m, "source_z_edges_m")
    source_y = _edges(source_y_edges_m, "source_y_edges_m")
    target_z = _edges(target_z_edges_m, "target_z_edges_m")
    target_y = _edges(target_y_edges_m, "target_y_edges_m")
    fields = {
        "total_mass_flux_anomaly": np.asarray(
            total_mass_flux_anomaly_kg_m2_s, dtype=float,
        ),
        "species_mass_flux": np.asarray(species_mass_flux_kg_m2_s, dtype=float),
        "sensible_enthalpy_flux": np.asarray(
            sensible_enthalpy_flux_w_m2, dtype=float,
        ),
        "lateral_momentum_flux_anomaly": np.asarray(
            lateral_momentum_flux_anomaly_pa, dtype=float,
        ),
        "vertical_momentum_flux_anomaly": np.asarray(
            vertical_momentum_flux_anomaly_pa, dtype=float,
        ),
    }
    expected_shape = (times.size, source_z.size - 1, source_y.size - 1)
    for name, field in fields.items():
        if field.shape != expected_shape or np.any(~np.isfinite(field)):
            raise ValueError(f"{name} must be finite with shape {expected_shape}")
    mapped = {
        name: conservative_overlap_remap_flux_density_2d(
            field,
            source_z_edges_m=source_z,
            source_y_edges_m=source_y,
            target_z_edges_m=target_z,
            target_y_edges_m=target_y,
        )
        for name, field in fields.items()
    }
    expected_shape = (times.size, source_z.size - 1, source_y.size - 1)
    advective_fields: dict[str, np.ndarray | None] = {}
    for name, value in (
        ("advective_species_mass_flux", advective_species_mass_flux_kg_m2_s),
        (
            "advective_sensible_enthalpy_flux",
            advective_sensible_enthalpy_flux_w_m2,
        ),
    ):
        if value is None:
            advective_fields[name] = None
            continue
        field = np.asarray(value, dtype=float)
        if field.shape != expected_shape or np.any(~np.isfinite(field)):
            raise ValueError(f"{name} must be finite with shape {expected_shape}")
        advective_fields[name] = conservative_overlap_remap_flux_density_2d(
            field,
            source_z_edges_m=source_z,
            source_y_edges_m=source_y,
            target_z_edges_m=target_z,
            target_y_edges_m=target_y,
        )
    residuals = {
        name: _relative_integral_residual(
            fields[name], mapped[name],
            source_z, source_y, target_z, target_y,
        )
        for name in fields
    }
    for name, field in advective_fields.items():
        if field is not None:
            residuals[name] = _relative_integral_residual(
                np.asarray(
                    advective_species_mass_flux_kg_m2_s
                    if name == "advective_species_mass_flux"
                    else advective_sensible_enthalpy_flux_w_m2,
                    dtype=float,
                ),
                field,
                source_z, source_y, target_z, target_y,
            )
    maximum_residual = max(residuals.values())
    return ConservativeNestedHandoff3DResult(
        schema=CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA,
        time_s=times.copy(),
        target_z_edges_m=target_z.copy(),
        target_y_edges_m=target_y.copy(),
        total_mass_flux_anomaly_kg_m2_s=mapped["total_mass_flux_anomaly"],
        species_mass_flux_kg_m2_s=mapped["species_mass_flux"],
        sensible_enthalpy_flux_w_m2=mapped["sensible_enthalpy_flux"],
        lateral_momentum_flux_anomaly_pa=mapped[
            "lateral_momentum_flux_anomaly"
        ],
        vertical_momentum_flux_anomaly_pa=mapped[
            "vertical_momentum_flux_anomaly"
        ],
        diagnostics={
            "source_id": source_id,
            "mapping": "exact rectangular face-overlap area weighting",
            "representation": (
                "plume flux anomalies embedded into an independently declared "
                "ambient outer inlet"
            ),
            "source_shape": list(expected_shape),
            "target_shape": [
                times.size, target_z.size - 1, target_y.size - 1,
            ],
            "relative_integral_residuals": residuals,
            "maximum_relative_integral_residual": maximum_residual,
            "source_species_backflow_cell_count": int(np.count_nonzero(
                fields["species_mass_flux"] < 0.0
            )),
            "target_species_backflow_cell_count": int(np.count_nonzero(
                mapped["species_mass_flux"] < 0.0
            )),
            "clipping_or_posthoc_renormalisation": False,
            "outer_ambient_state_supplied_by_this_operator": False,
            "limitations": [
                "numerical remap component only; outer transient boundary coupling remains required",
                "source and target planes must share the same physical z-y frame",
                "total mass flux is transferred as an anomaly relative to declared ambient inflow",
            ],
        },
        advective_species_mass_flux_kg_m2_s=advective_fields[
            "advective_species_mass_flux"
        ],
        advective_sensible_enthalpy_flux_w_m2=advective_fields[
            "advective_sensible_enthalpy_flux"
        ],
    )


def remap_conservative_nested_handoff_result_3d(
    handoff: ConservativeNestedHandoff3DResult,
    *,
    target_z_edges_m: np.ndarray,
    target_y_edges_m: np.ndarray,
    source_id: str | None = None,
) -> ConservativeNestedHandoff3DResult:
    """Area-remap an already recorded handoff history to an outer plane."""

    if not isinstance(handoff, ConservativeNestedHandoff3DResult):
        raise TypeError("handoff must be ConservativeNestedHandoff3DResult")
    inherited_source = handoff.diagnostics.get("source_id", "recorded handoff")
    return remap_conservative_nested_handoff_3d(
        time_s=handoff.time_s,
        source_z_edges_m=handoff.target_z_edges_m,
        source_y_edges_m=handoff.target_y_edges_m,
        target_z_edges_m=target_z_edges_m,
        target_y_edges_m=target_y_edges_m,
        total_mass_flux_anomaly_kg_m2_s=(
            handoff.total_mass_flux_anomaly_kg_m2_s
        ),
        species_mass_flux_kg_m2_s=handoff.species_mass_flux_kg_m2_s,
        sensible_enthalpy_flux_w_m2=handoff.sensible_enthalpy_flux_w_m2,
        lateral_momentum_flux_anomaly_pa=(
            handoff.lateral_momentum_flux_anomaly_pa
        ),
        vertical_momentum_flux_anomaly_pa=(
            handoff.vertical_momentum_flux_anomaly_pa
        ),
        source_id=(
            str(inherited_source) if source_id is None else source_id
        ),
        advective_species_mass_flux_kg_m2_s=(
            handoff.advective_species_mass_flux_kg_m2_s
        ),
        advective_sensible_enthalpy_flux_w_m2=(
            handoff.advective_sensible_enthalpy_flux_w_m2
        ),
    )


__all__ = [
    "CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA",
    "ConservativeNestedHandoff3DResult",
    "NestedOuterInletFrame3D",
    "conservative_overlap_remap_flux_density_2d",
    "nested_outer_inlet_frame_at",
    "remap_conservative_nested_handoff_result_3d",
    "remap_conservative_nested_handoff_3d",
]
