"""Append-only HDF5 persistence for conservative handoff plane histories."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .binary_mixture_handoff_plane_3d import (
    BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
    BinaryMixtureHandoffPlane3DFrame,
)
from .conservative_nested_handoff_3d import (
    CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA,
    ConservativeNestedHandoff3DResult,
    NestedOuterInletFrame3D,
    conservative_overlap_remap_flux_density_2d,
    nested_outer_inlet_frame_at,
)
from .binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure


LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA = "degali.conservative-handoff-hdf5-3d.v1"
CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA = "degali.conservative-handoff-hdf5-3d.v2"
LEGACY_BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA = "degali.binary-mixture-handoff-plane-3d.v1"

_FLUX_FIELDS = (
    "total_mass_flux_anomaly_kg_m2_s",
    "species_mass_flux_kg_m2_s",
    "sensible_enthalpy_flux_w_m2",
    "lateral_momentum_flux_anomaly_pa",
    "vertical_momentum_flux_anomaly_pa",
)
_STATE_FLUX_FIELDS = (
    "advective_species_mass_flux_kg_m2_s",
    "advective_sensible_enthalpy_flux_w_m2",
)
_FIELDS = _FLUX_FIELDS + _STATE_FLUX_FIELDS
_LAYER_PACKET_SCHEMA = 'degali.marked-layer-donor-packet.v1'
_LAYER_FIELDS = ('mixture_mass_flux_kg_m2_s', 'species_mass_flux_kg_m2_s',
                 'momentum_x_flux_pa', 'momentum_y_flux_pa', 'momentum_z_flux_pa')


def _persisted_fields(handle):
    fields = list(_FIELDS)
    if 'dense_layer_packet_schema' in handle.attrs:
        if handle.attrs['dense_layer_packet_schema'] != _LAYER_PACKET_SCHEMA:
            raise ValueError('unsupported dense-layer packet extension')
        fields += ['dense_layer_packet/'+name for name in _LAYER_FIELDS]
    return fields


def _edges(value: np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if (
        result.ndim != 1 or result.size < 2
        or np.any(~np.isfinite(result)) or np.any(np.diff(result) <= 0.0)
    ):
        raise ValueError(f"{name} must be finite strictly increasing edges")
    return result


class ConservativeHandoffHdf5Writer3D:
    """Durably append one frame at a time without retaining history in RAM."""

    def __init__(
        self,
        path: str | Path,
        *,
        source_z_edges_m: np.ndarray,
        source_y_edges_m: np.ndarray,
        source_id: str,
        metadata: dict[str, Any] | None = None,
        resume: bool = False,
    ) -> None:
        import h5py

        self.path = Path(path).resolve()
        self.source_z_edges_m = _edges(source_z_edges_m, "source_z_edges_m")
        self.source_y_edges_m = _edges(source_y_edges_m, "source_y_edges_m")
        self.plane_shape = (
            self.source_z_edges_m.size - 1,
            self.source_y_edges_m.size - 1,
        )
        if not isinstance(source_id, str) or not source_id.strip():
            raise ValueError("source_id must be a non-empty string")
        self.source_id = source_id
        self.metadata = {} if metadata is None else dict(metadata)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if resume:
            if not self.path.is_file():
                raise FileNotFoundError(f"handoff checkpoint does not exist: {self.path}")
            self._file = h5py.File(self.path, "r+")
            self._validate_existing()
        else:
            if self.path.exists():
                raise FileExistsError(
                    f"handoff checkpoint already exists: {self.path}"
                )
            self._file = h5py.File(self.path, "x")
            self._file.attrs["schema"] = CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA
            self._file.attrs["frame_schema"] = BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA
            self._file.attrs["source_id"] = source_id
            self._file.attrs["metadata_json"] = json.dumps(
                self.metadata, sort_keys=True, ensure_ascii=False,
            )
            self._file.create_dataset("source_z_edges_m", data=self.source_z_edges_m)
            self._file.create_dataset("source_y_edges_m", data=self.source_y_edges_m)
            self._file.create_dataset(
                "time_s", shape=(0,), maxshape=(None,), chunks=(256,), dtype="f8",
            )
            chunks = (1, *self.plane_shape)
            for name in _FIELDS:
                self._file.create_dataset(
                    name,
                    shape=(0, *self.plane_shape),
                    maxshape=(None, *self.plane_shape),
                    chunks=chunks,
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                    dtype="f8",
                )
            self._file.flush()

    def _validate_existing(self) -> None:
        if self._file.attrs.get("schema") != CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA:
            raise ValueError("unsupported conservative handoff HDF5 schema")
        if self._file.attrs.get("frame_schema") != BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA:
            raise ValueError("handoff frame schema does not match")
        if self._file.attrs.get("source_id") != self.source_id:
            raise ValueError("handoff source_id does not match checkpoint")
        stored_metadata = json.loads(self._file.attrs.get("metadata_json", "{}"))
        if stored_metadata != self.metadata:
            raise ValueError("handoff metadata does not match checkpoint")
        if not np.array_equal(
            self._file["source_z_edges_m"][:], self.source_z_edges_m,
        ) or not np.array_equal(
            self._file["source_y_edges_m"][:], self.source_y_edges_m,
        ):
            raise ValueError("handoff plane edges do not match checkpoint")
        count = self._file["time_s"].shape[0]
        if any(self._file[name].shape != (count, *self.plane_shape) for name in _persisted_fields(self._file)):
            raise ValueError("handoff checkpoint datasets have inconsistent lengths")

    @property
    def frame_count(self) -> int:
        return int(self._file["time_s"].shape[0])

    @property
    def last_time_s(self) -> float | None:
        if self.frame_count == 0:
            return None
        return float(self._file["time_s"][-1])

    def truncate_after(self, committed_time_s: float) -> int:
        """Discard frames newer than an atomically committed solver restart.

        The committed time must already exist in the history.  Requiring an
        exact frame (within floating-point roundoff) prevents recovery from
        silently inventing or interpolating a handoff boundary.
        """

        if self._file.id.valid == 0:
            raise RuntimeError("handoff checkpoint writer is closed")
        committed = float(committed_time_s)
        if not np.isfinite(committed):
            raise ValueError("committed handoff time must be finite")
        times = np.asarray(self._file["time_s"][:], dtype=float)
        tolerance = 1.0e-10 * max(1.0, abs(committed))
        matches = np.flatnonzero(np.abs(times - committed) <= tolerance)
        if matches.size != 1:
            raise ValueError(
                "committed solver restart time is not a unique handoff frame"
            )
        retained = int(matches[0]) + 1
        removed = int(times.size - retained)
        if removed:
            self._file["time_s"].resize((retained,))
            for name in _persisted_fields(self._file):
                self._file[name].resize((retained, *self.plane_shape))
            self._file.flush()
        return removed

    def __call__(self, frame: BinaryMixtureHandoffPlane3DFrame) -> None:
        if not isinstance(frame, BinaryMixtureHandoffPlane3DFrame):
            raise TypeError("frame must be BinaryMixtureHandoffPlane3DFrame")
        if self._file.id.valid == 0:
            raise RuntimeError("handoff checkpoint writer is closed")
        if self.last_time_s is not None and frame.time_s <= self.last_time_s:
            raise ValueError("handoff frame times must be strictly increasing")
        fields = {
            name: np.asarray(getattr(frame, name), dtype=float)
            for name in _FIELDS
        }
        has_packet = frame.dense_layer_packet is not None
        stored_packet = 'dense_layer_packet_schema' in self._file.attrs
        if self.frame_count and has_packet != stored_packet:
            raise ValueError('dense-layer packet presence must remain constant across frames')
        if has_packet:
            if set(frame.dense_layer_packet) != set(_LAYER_FIELDS):
                raise ValueError('dense-layer packet must contain named mass and xyz momentum fields')
            fields.update({'dense_layer_packet/'+name: np.asarray(frame.dense_layer_packet[name],float)
                           for name in _LAYER_FIELDS})
        for name, field in fields.items():
            if field.shape != self.plane_shape or np.any(~np.isfinite(field)):
                raise ValueError(
                    f"{name} must be finite with shape {self.plane_shape}"
                )
        if has_packet and not stored_packet:
            self._file.attrs['dense_layer_packet_schema'] = _LAYER_PACKET_SCHEMA
            group = self._file.create_group('dense_layer_packet')
            for name in _LAYER_FIELDS:
                group.create_dataset(name, shape=(0,*self.plane_shape),
                    maxshape=(None,*self.plane_shape), chunks=(1,*self.plane_shape),
                    compression='gzip',compression_opts=1,shuffle=True,dtype='f8')
        index = self.frame_count
        self._file["time_s"].resize((index + 1,))
        self._file["time_s"][index] = frame.time_s
        for name, field in fields.items():
            dataset = self._file[name]
            dataset.resize((index + 1, *self.plane_shape))
            dataset[index] = field
        self._file.flush()

    def close(self) -> None:
        if self._file.id.valid:
            self._file.flush()
            self._file.close()

    def __enter__(self) -> "ConservativeHandoffHdf5Writer3D":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


class ConservativeHandoffHdf5Reader3D:
    """Lazily interpolate and remap a durable handoff history.

    Only the one or two source frames bracketing each requested outer time are
    read. This avoids materialising the full transient history in RAM. Version
    2 also carries advective species/enthalpy components separately from the
    net conservative fluxes; version-1 files remain readable with legacy state
    reconstruction explicitly identified in diagnostics.
    """

    def __init__(
        self, path: str | Path, *,
        target_z_edges_m: np.ndarray,
        target_y_edges_m: np.ndarray,
    ) -> None:
        import h5py

        self.path = Path(path).resolve()
        self._file = h5py.File(self.path, "r")
        try:
            schema = self._file.attrs.get("schema")
            if schema not in {
                CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
                LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
            }:
                raise ValueError("unsupported conservative handoff HDF5 schema")
            legacy = schema == LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA
            expected_frame_schema = (
                LEGACY_BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA if legacy
                else BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA
            )
            if self._file.attrs.get("frame_schema") != expected_frame_schema:
                raise ValueError("handoff frame schema does not match history schema")
            self.time_s = np.asarray(self._file["time_s"][:], dtype=float)
            if (
                self.time_s.size < 1 or np.any(~np.isfinite(self.time_s))
                or np.any(np.diff(self.time_s) <= 0.0)
            ):
                raise ValueError("handoff checkpoint must contain increasing frames")
            self.source_z_edges_m = _edges(
                self._file["source_z_edges_m"][:], "source_z_edges_m",
            )
            self.source_y_edges_m = _edges(
                self._file["source_y_edges_m"][:], "source_y_edges_m",
            )
            self.target_z_edges_m = _edges(target_z_edges_m, "target_z_edges_m")
            self.target_y_edges_m = _edges(target_y_edges_m, "target_y_edges_m")
            expected = (
                self.time_s.size,
                self.source_z_edges_m.size - 1,
                self.source_y_edges_m.size - 1,
            )
            required_fields = _FLUX_FIELDS if legacy else _FIELDS
            if any(self._file[name].shape != expected for name in required_fields):
                raise ValueError("handoff checkpoint datasets have inconsistent shapes")
            conservative_overlap_remap_flux_density_2d(
                np.zeros(expected[1:], dtype=float),
                source_z_edges_m=self.source_z_edges_m,
                source_y_edges_m=self.source_y_edges_m,
                target_z_edges_m=self.target_z_edges_m,
                target_y_edges_m=self.target_y_edges_m,
            )
            self.source_id = str(self._file.attrs["source_id"])
            self.metadata = json.loads(self._file.attrs.get("metadata_json", "{}"))
            self.legacy_state_from_net_flux = legacy
            self.maximum_relative_integral_residual = 0.0
            self.query_count = 0
            self.source_frame_read_count = 0
        except BaseException:
            self._file.close()
            raise

    @property
    def frame_count(self) -> int:
        return int(self.time_s.size)

    def _bracket(self, query_time_s: float) -> tuple[int, int, float]:
        query = float(query_time_s)
        if not np.isfinite(query):
            raise ValueError("time_s must be finite")
        tolerance = 1.0e-12 * max(1.0, abs(query))
        if query < self.time_s[0] - tolerance or query > self.time_s[-1] + tolerance:
            raise ValueError("outer inlet query time lies outside handoff history")
        nearest = int(np.argmin(np.abs(self.time_s - query)))
        if abs(float(self.time_s[nearest]) - query) <= tolerance:
            return nearest, nearest, 0.0
        upper = int(np.searchsorted(self.time_s, query, side="right"))
        lower = upper - 1
        fraction = float(
            (query - self.time_s[lower])
            / (self.time_s[upper] - self.time_s[lower])
        )
        return lower, upper, fraction

    def _source_field_at(
        self, name: str, lower: int, upper: int, fraction: float,
    ) -> np.ndarray:
        left = np.asarray(self._file[name][lower], dtype=float)
        self.source_frame_read_count += 1
        if upper == lower:
            return left
        right = np.asarray(self._file[name][upper], dtype=float)
        self.source_frame_read_count += 1
        return (1.0 - fraction) * left + fraction * right

    def frame_at(
        self, *, time_s: float,
        ambient_total_mass_flux_kg_m2_s: float | np.ndarray,
        closure: BinaryIdealGasMixtureClosure,
    ) -> NestedOuterInletFrame3D:
        """Return one remapped inlet while retaining only bracketing frames."""

        if self._file.id.valid == 0:
            raise RuntimeError("handoff checkpoint reader is closed")
        lower, upper, fraction = self._bracket(time_s)
        mapped = {}
        residuals = {}
        source_area = (
            np.diff(self.source_z_edges_m)[:, None]
            * np.diff(self.source_y_edges_m)[None, :]
        )
        target_area = (
            np.diff(self.target_z_edges_m)[:, None]
            * np.diff(self.target_y_edges_m)[None, :]
        )
        for name in _FLUX_FIELDS:
            source = self._source_field_at(name, lower, upper, fraction)
            target = conservative_overlap_remap_flux_density_2d(
                source,
                source_z_edges_m=self.source_z_edges_m,
                source_y_edges_m=self.source_y_edges_m,
                target_z_edges_m=self.target_z_edges_m,
                target_y_edges_m=self.target_y_edges_m,
            )
            mapped[name] = target
            scale = max(float(np.sum(np.abs(source) * source_area)), 1.0e-30)
            residuals[name] = float(
                abs(np.sum(target * target_area) - np.sum(source * source_area))
                / scale
            )
        for name in _STATE_FLUX_FIELDS:
            if self.legacy_state_from_net_flux:
                mapped[name] = None
                continue
            source = self._source_field_at(name, lower, upper, fraction)
            mapped[name] = conservative_overlap_remap_flux_density_2d(
                source,
                source_z_edges_m=self.source_z_edges_m,
                source_y_edges_m=self.source_y_edges_m,
                target_z_edges_m=self.target_z_edges_m,
                target_y_edges_m=self.target_y_edges_m,
            )
            scale = max(float(np.sum(np.abs(source) * source_area)), 1.0e-30)
            residuals[name] = float(
                abs(np.sum(mapped[name] * target_area) - np.sum(source * source_area))
                / scale
            )
        maximum_residual = max(residuals.values())
        self.maximum_relative_integral_residual = max(
            self.maximum_relative_integral_residual, maximum_residual,
        )
        self.query_count += 1
        handoff = ConservativeNestedHandoff3DResult(
            schema=CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA,
            time_s=np.asarray((float(time_s),)),
            target_z_edges_m=self.target_z_edges_m.copy(),
            target_y_edges_m=self.target_y_edges_m.copy(),
            total_mass_flux_anomaly_kg_m2_s=mapped[
                "total_mass_flux_anomaly_kg_m2_s"
            ][None, :, :],
            species_mass_flux_kg_m2_s=mapped[
                "species_mass_flux_kg_m2_s"
            ][None, :, :],
            sensible_enthalpy_flux_w_m2=mapped[
                "sensible_enthalpy_flux_w_m2"
            ][None, :, :],
            lateral_momentum_flux_anomaly_pa=mapped[
                "lateral_momentum_flux_anomaly_pa"
            ][None, :, :],
            vertical_momentum_flux_anomaly_pa=mapped[
                "vertical_momentum_flux_anomaly_pa"
            ][None, :, :],
            diagnostics={
                "source_id": self.source_id,
                "source_file": str(self.path),
                "mapping": "lazy exact rectangular face-overlap area weighting",
                "source_frame_indices": [lower, upper],
                "temporal_interpolation_fraction": fraction,
                "relative_integral_residuals": residuals,
                "maximum_relative_integral_residual": maximum_residual,
                "clipping_or_posthoc_renormalisation": False,
            },
            advective_species_mass_flux_kg_m2_s=(
                None if mapped["advective_species_mass_flux_kg_m2_s"] is None
                else mapped["advective_species_mass_flux_kg_m2_s"][None, :, :]
            ),
            advective_sensible_enthalpy_flux_w_m2=(
                None if mapped["advective_sensible_enthalpy_flux_w_m2"] is None
                else mapped["advective_sensible_enthalpy_flux_w_m2"][None, :, :]
            ),
        )
        frame = nested_outer_inlet_frame_at(
            handoff, time_s=float(time_s),
            ambient_total_mass_flux_kg_m2_s=ambient_total_mass_flux_kg_m2_s,
            closure=closure,
        )
        frame.diagnostics.update({
            "durable_handoff_lazy_read": True,
            "source_frame_indices": [lower, upper],
            "temporal_interpolation_fraction": fraction,
            "maximum_relative_integral_remap_residual": maximum_residual,
        })
        return frame

    def close(self) -> None:
        if self._file.id.valid:
            self._file.close()

    def __enter__(self) -> "ConservativeHandoffHdf5Reader3D":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def load_conservative_handoff_hdf5_3d(
    path: str | Path,
) -> ConservativeNestedHandoff3DResult:
    """Load a completed or checkpointed history into the remap contract."""

    import h5py

    source = Path(path).resolve()
    with h5py.File(source, "r") as handle:
        hdf5_schema = handle.attrs.get("schema")
        if hdf5_schema not in {
            CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
            LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
        }:
            raise ValueError("unsupported conservative handoff HDF5 schema")
        legacy = hdf5_schema == LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA
        times = np.asarray(handle["time_s"][:], dtype=float)
        if times.size < 1 or np.any(np.diff(times) <= 0.0):
            raise ValueError("handoff checkpoint must contain increasing frames")
        z_edges = np.asarray(handle["source_z_edges_m"][:], dtype=float)
        y_edges = np.asarray(handle["source_y_edges_m"][:], dtype=float)
        required_fields = _FLUX_FIELDS if legacy else _FIELDS
        fields = {
            name: np.asarray(handle[name][:], dtype=float)
            for name in required_fields
        }
        source_id = str(handle.attrs["source_id"])
        metadata = json.loads(handle.attrs.get("metadata_json", "{}"))
    return ConservativeNestedHandoff3DResult(
        schema=CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA,
        time_s=times,
        target_z_edges_m=z_edges,
        target_y_edges_m=y_edges,
        total_mass_flux_anomaly_kg_m2_s=fields[
            "total_mass_flux_anomaly_kg_m2_s"
        ],
        species_mass_flux_kg_m2_s=fields["species_mass_flux_kg_m2_s"],
        sensible_enthalpy_flux_w_m2=fields["sensible_enthalpy_flux_w_m2"],
        lateral_momentum_flux_anomaly_pa=fields[
            "lateral_momentum_flux_anomaly_pa"
        ],
        vertical_momentum_flux_anomaly_pa=fields[
            "vertical_momentum_flux_anomaly_pa"
        ],
        diagnostics={
            "source_id": source_id,
            "source_file": str(source),
            "frame_count": int(times.size),
            "metadata": metadata,
            "mapping": "identity load from durable source-plane history",
            "clipping_or_posthoc_renormalisation": False,
        },
        advective_species_mass_flux_kg_m2_s=fields.get(
            "advective_species_mass_flux_kg_m2_s"
        ),
        advective_sensible_enthalpy_flux_w_m2=fields.get(
            "advective_sensible_enthalpy_flux_w_m2"
        ),
    )


__all__ = [
    "CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA",
    "ConservativeHandoffHdf5Reader3D",
    "ConservativeHandoffHdf5Writer3D",
    "load_conservative_handoff_hdf5_3d",
]
