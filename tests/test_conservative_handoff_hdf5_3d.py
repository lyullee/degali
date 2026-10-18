import numpy as np
import pytest

from degali.addons.binary_mixture_handoff_plane_3d import (
    BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
    BinaryMixtureHandoffPlane3DFrame,
)
from degali.addons.conservative_handoff_hdf5_3d import (
    CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
    LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA,
    ConservativeHandoffHdf5Reader3D,
    ConservativeHandoffHdf5Writer3D,
    load_conservative_handoff_hdf5_3d,
)
from degali.addons.conservative_nested_handoff_3d import (
    nested_outer_inlet_frame_at,
    remap_conservative_nested_handoff_result_3d,
)
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.hdf5_dense_layer_interface import Hdf5DenseLayerPacketSource


def _frame(time_s: float, scale: float = 1.0):
    shape = (2, 3)
    base = scale * np.arange(1.0, 7.0).reshape(shape)
    return BinaryMixtureHandoffPlane3DFrame(
        schema=BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
        time_s=time_s,
        x_face_index=4,
        total_mass_flux_anomaly_kg_m2_s=0.1 * base,
        advective_species_mass_flux_kg_m2_s=0.019 * base,
        species_mass_flux_kg_m2_s=0.02 * base,
        advective_sensible_enthalpy_flux_w_m2=-99.0 * base,
        sensible_enthalpy_flux_w_m2=-100.0 * base,
        lateral_momentum_flux_anomaly_pa=0.03 * base,
        vertical_momentum_flux_anomaly_pa=-0.01 * base,
        diagnostics={},
    )


def test_hdf5_writer_flushes_resumes_and_loads_conservative_history(tmp_path):
    path = tmp_path / "handoff.h5"
    z_edges = np.asarray((0.0, 1.0, 2.0))
    y_edges = np.asarray((-1.5, -0.5, 0.5, 1.5))
    metadata = {"case": "checkpoint-test", "x_face_m": 4.0}
    with ConservativeHandoffHdf5Writer3D(
        path,
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        source_id="test source",
        metadata=metadata,
    ) as writer:
        writer(_frame(0.0, 0.0))
        writer(_frame(1.0, 1.0))
        assert writer.frame_count == 2
        assert writer.last_time_s == 1.0

    with ConservativeHandoffHdf5Writer3D(
        path,
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        source_id="test source",
        metadata=metadata,
        resume=True,
    ) as writer:
        assert writer.frame_count == 2
        writer(_frame(2.0, 0.5))
        assert writer.frame_count == 3

    with ConservativeHandoffHdf5Writer3D(
        path,
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        source_id="test source",
        metadata=metadata,
        resume=True,
    ) as writer:
        assert writer.truncate_after(1.0) == 1
        assert writer.last_time_s == 1.0
        writer(_frame(2.0, 0.5))

    loaded = load_conservative_handoff_hdf5_3d(path)
    assert loaded.time_s.tolist() == [0.0, 1.0, 2.0]
    assert loaded.species_mass_flux_kg_m2_s[1] == pytest.approx(
        _frame(1.0).species_mass_flux_kg_m2_s
    )
    assert loaded.advective_species_mass_flux_kg_m2_s[1] == pytest.approx(
        _frame(1.0).advective_species_mass_flux_kg_m2_s
    )
    assert loaded.diagnostics["frame_count"] == 3
    assert loaded.diagnostics["metadata"] == metadata

    remapped = remap_conservative_nested_handoff_result_3d(
        loaded,
        target_z_edges_m=np.asarray((0.0, 2.0)),
        target_y_edges_m=np.asarray((-2.0, 0.0, 2.0)),
    )
    assert remapped.diagnostics["maximum_relative_integral_residual"] < 1.0e-14
    assert remapped.species_mass_flux_kg_m2_s.shape == (3, 1, 2)
    assert remapped.advective_species_mass_flux_kg_m2_s.shape == (3, 1, 2)


def test_reader_keeps_legacy_v1_histories_readable_as_legacy_state(tmp_path):
    import h5py
    import json

    path = tmp_path / "legacy.h5"
    z_edges = np.asarray((0.0, 1.0, 2.0))
    y_edges = np.asarray((-1.5, -0.5, 0.5, 1.5))
    shape = (2, 3)
    with h5py.File(path, "w") as handle:
        handle.attrs["schema"] = LEGACY_CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA
        handle.attrs["frame_schema"] = "degali.binary-mixture-handoff-plane-3d.v1"
        handle.attrs["source_id"] = "legacy source"
        handle.attrs["metadata_json"] = json.dumps({"legacy": True})
        handle.create_dataset("time_s", data=np.asarray((0.0, 1.0)))
        handle.create_dataset("source_z_edges_m", data=z_edges)
        handle.create_dataset("source_y_edges_m", data=y_edges)
        handle.create_dataset(
            "total_mass_flux_anomaly_kg_m2_s", data=np.zeros((2, *shape)),
        )
        handle.create_dataset(
            "species_mass_flux_kg_m2_s", data=np.full((2, *shape), 0.02),
        )
        handle.create_dataset(
            "sensible_enthalpy_flux_w_m2", data=np.zeros((2, *shape)),
        )
        handle.create_dataset(
            "lateral_momentum_flux_anomaly_pa", data=np.zeros((2, *shape)),
        )
        handle.create_dataset(
            "vertical_momentum_flux_anomaly_pa", data=np.zeros((2, *shape)),
        )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0,
    )
    with ConservativeHandoffHdf5Reader3D(
        path, target_z_edges_m=z_edges, target_y_edges_m=y_edges,
    ) as reader:
        frame = reader.frame_at(
            time_s=0.5,
            ambient_total_mass_flux_kg_m2_s=np.ones(shape),
            closure=closure,
        )
    assert frame.diagnostics["thermodynamic_state_source"] == "legacy_net_flux_ratio"
    assert frame.species_mass_flux_kg_m2_s == pytest.approx(
        np.full(shape, 0.02)
    )


def test_hdf5_writer_rejects_overwrite_and_nonincreasing_resume_time(tmp_path):
    path = tmp_path / "handoff.h5"
    z_edges = np.asarray((0.0, 1.0, 2.0))
    y_edges = np.asarray((-1.5, -0.5, 0.5, 1.5))
    with ConservativeHandoffHdf5Writer3D(
        path,
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        source_id="test source",
    ) as writer:
        writer(_frame(0.0, 0.0))
    with pytest.raises(FileExistsError):
        ConservativeHandoffHdf5Writer3D(
            path,
            source_z_edges_m=z_edges,
            source_y_edges_m=y_edges,
            source_id="test source",
        )
    with ConservativeHandoffHdf5Writer3D(
        path,
        source_z_edges_m=z_edges,
        source_y_edges_m=y_edges,
        source_id="test source",
        resume=True,
    ) as writer:
        with pytest.raises(ValueError, match="strictly increasing"):
            writer(_frame(0.0, 1.0))
        with pytest.raises(ValueError, match="not a unique handoff frame"):
            writer.truncate_after(0.5)


def test_hdf5_schema_constant_is_stable():
    assert CONSERVATIVE_HANDOFF_HDF5_3D_SCHEMA.endswith(".v2")


def test_hdf5_dense_layer_interface_collapses_height_without_ambient_double_count(tmp_path):
    from dataclasses import replace
    path = tmp_path / "handoff-interface.h5"
    z_edges = np.asarray((0.0, 1.0, 2.0))
    y_edges = np.asarray((-1.5, -0.5, 0.5, 1.5))
    with ConservativeHandoffHdf5Writer3D(
        path, source_z_edges_m=z_edges, source_y_edges_m=y_edges,
        source_id="interface-test",
    ) as writer:
        for t,scale in ((0.,1.),(1.,.5)):
            frame = _frame(t,scale)
            mixture = np.maximum(frame.total_mass_flux_anomaly_kg_m2_s,frame.species_mass_flux_kg_m2_s)
            writer(replace(frame,dense_layer_packet={
                'mixture_mass_flux_kg_m2_s':mixture,
                'species_mass_flux_kg_m2_s':frame.species_mass_flux_kg_m2_s,
                'momentum_x_flux_pa':mixture*3.,
                'momentum_y_flux_pa':mixture*.5,
                'momentum_z_flux_pa':mixture*100.,
            }))
    with Hdf5DenseLayerPacketSource(path, target_y_edges_m=y_edges) as source:
        packet = source(0.5, 0.25, 3, 1.0)
        assert packet.shape == (4, 3)
        assert np.all(packet[0] >= packet[1])
        assert np.all(packet[1] >= 0.0)
        assert np.all(np.isfinite(packet[2:4]))
        assert packet[2] == pytest.approx(3.*packet[0])
        assert packet[3] == pytest.approx(.5*packet[0])
        assert source.discarded_vertical_impulse_kg_m_s == pytest.approx(100.*source.integrated_total_mass_kg)
        assert source.calls == 1


def test_legacy_handoff_cannot_invent_normal_momentum_for_layer(tmp_path):
    path=tmp_path/'old.h5'
    with ConservativeHandoffHdf5Writer3D(path,source_z_edges_m=np.array((0.,1.,2.)),
            source_y_edges_m=np.array((-1.5,-.5,.5,1.5)),source_id='old') as writer:
        writer(_frame(0.,1.));writer(_frame(1.,.5))
    with pytest.raises(ValueError,match='lacks named donor x/y momentum'):
        Hdf5DenseLayerPacketSource(path,target_y_edges_m=np.array((-1.5,-.5,.5,1.5)))


def test_lazy_reader_matches_eager_interpolation_and_remap(tmp_path):
    path = tmp_path / "handoff.h5"
    z_edges = np.asarray((0.0, 1.0, 2.0))
    y_edges = np.asarray((-1.5, -0.5, 0.5, 1.5))
    with ConservativeHandoffHdf5Writer3D(
        path, source_z_edges_m=z_edges, source_y_edges_m=y_edges,
        source_id="lazy-test",
    ) as writer:
        writer(_frame(0.0, 0.0))
        writer(_frame(1.0, 1.0))
        writer(_frame(2.0, 0.5))
    target_z = np.asarray((0.0, 2.0))
    target_y = np.asarray((-2.0, 0.0, 2.0))
    eager = remap_conservative_nested_handoff_result_3d(
        load_conservative_handoff_hdf5_3d(path),
        target_z_edges_m=target_z, target_y_edges_m=target_y,
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0,
    )
    ambient = np.full((1, 2), closure.ambient_density_kg_m3)
    expected = nested_outer_inlet_frame_at(
        eager, time_s=0.5, ambient_total_mass_flux_kg_m2_s=ambient,
        closure=closure,
    )
    with ConservativeHandoffHdf5Reader3D(
        path, target_z_edges_m=target_z, target_y_edges_m=target_y,
    ) as reader:
        actual = reader.frame_at(
            time_s=0.5, ambient_total_mass_flux_kg_m2_s=ambient,
            closure=closure,
        )
        assert reader.frame_count == 3
        assert reader.source_frame_read_count == 14
        assert reader.maximum_relative_integral_residual < 1.0e-14
    for name in (
        "total_mass_flux_kg_m2_s", "species_mass_flux_kg_m2_s",
        "sensible_enthalpy_flux_w_m2", "velocity_u_m_s",
        "velocity_v_m_s", "velocity_w_m_s",
    ):
        assert getattr(actual, name) == pytest.approx(getattr(expected, name))
    assert actual.diagnostics["durable_handoff_lazy_read"]
