import numpy as np
import pytest
from types import SimpleNamespace

from degali.addons.binary_mixture_handoff_plane_3d import (
    BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
    BinaryMixtureHandoffPlane3DFrame,
)
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.conservative_handoff_hdf5_3d import ConservativeHandoffHdf5Writer3D
from degali.addons.nested_outer_binary_flow_3d import initialise_nested_outer_binary_flow_3d
import tools.run_thorney_ti021_outer_transient_checkpointed as outer_runner
from tools.run_thorney_ti021_outer_transient_checkpointed import (
    _configuration_sha256,
    _read_restart,
    _write_restart,
    sample_trilinear_many,
)


def test_outer_trilinear_sampling_clamps_physical_boundaries():
    x = np.asarray((0.5, 1.5, 2.5))
    y = np.asarray((0.5, 1.5, 2.5))
    z = np.asarray((0.5, 1.5, 2.5))
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    field = xx + 2.0 * yy + 3.0 * zz
    positions = np.asarray(((1.0, 1.0, 1.0), (0.0, 0.0, 0.0), (3.0, 3.0, 3.0)))
    sampled = sample_trilinear_many(field, x, y, z, positions)
    assert sampled == pytest.approx((6.0, 3.0, 15.0))


def test_outer_restart_round_trip_reconstructs_thermodynamics(tmp_path):
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0,
    )
    state = initialise_nested_outer_binary_flow_3d(
        cells=(3, 2, 2), closure=closure, ambient_u_m_s=1.0,
    )
    state = state.__class__(**{**state.__dict__, "time_s": 10.0})
    maximum = np.ones((2, 2, 3))
    config_hash = _configuration_sha256({"case": "round-trip", "dx": 1.0})
    path = tmp_path / "restart.npz"
    _write_restart(path, state, maximum, configuration_sha256=config_hash)
    restored, restored_maximum = _read_restart(
        path, closure=closure, configuration_sha256=config_hash,
    )
    assert restored.time_s == pytest.approx(10.0)
    assert restored.u_face_m_s == pytest.approx(state.u_face_m_s)
    assert restored.thermodynamic_state.mixture_density_kg_m3 == pytest.approx(
        state.thermodynamic_state.mixture_density_kg_m3,
    )
    assert restored_maximum == pytest.approx(maximum)
    with pytest.raises(ValueError, match="configuration SHA-256"):
        _read_restart(path, closure=closure, configuration_sha256="0" * 64)


def test_outer_checkpointed_runner_executes_lazy_two_step_path(tmp_path, monkeypatch):
    inner = tmp_path / "inner"
    inner.mkdir()
    workbook = tmp_path / "workbook.xls"
    workbook.write_bytes(b"synthetic coordinate-only workbook placeholder")
    property_path = tmp_path / "properties.json"
    property_path.write_text(
        '{"solver_policy":{"molecular_species_diffusivity_m2_s":9.8e-6},'
        '"pure_r12_ambient_properties":{"thermal_diffusivity_m2_s":3.1e-6}}',
        encoding="utf-8",
    )
    (inner / "ti021-inner-transient.json").write_text(
        __import__("json").dumps({
            "all_gates_pass": True,
            "configuration": {"duration_s": 0.25},
            "provenance": {
                "workbook_path": str(workbook),
                "property_audit_path": str(property_path),
            },
        }), encoding="utf-8",
    )
    (inner / "ti021-inner-checkpoint.json").write_text(
        '{"status":"finalized"}', encoding="utf-8",
    )
    z_edges = np.arange(0.0, 45.0, 5.0)
    y_edges = np.arange(-100.0, 105.0, 5.0)
    shape = (z_edges.size - 1, y_edges.size - 1)
    zeros = np.zeros(shape)
    with ConservativeHandoffHdf5Writer3D(
        inner / "ti021-inner-handoff.h5",
        source_z_edges_m=z_edges, source_y_edges_m=y_edges,
        source_id="synthetic zero-anomaly handoff",
    ) as writer:
        for time_s in (0.0, 0.125, 0.25):
            writer(BinaryMixtureHandoffPlane3DFrame(
                schema=BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
                time_s=time_s, x_face_index=1,
                total_mass_flux_anomaly_kg_m2_s=zeros,
                advective_species_mass_flux_kg_m2_s=zeros,
                species_mass_flux_kg_m2_s=zeros,
                advective_sensible_enthalpy_flux_w_m2=zeros,
                sensible_enthalpy_flux_w_m2=zeros,
                lateral_momentum_flux_anomaly_pa=zeros,
                vertical_momentum_flux_anomaly_pa=zeros,
                diagnostics={},
            ))
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0,
    )
    monkeypatch.setattr(outer_runner, "read_thorney_ti021_workbook", lambda path: SimpleNamespace())
    labels = np.asarray([f"TI021-S{i + 1:03d}" for i in range(94)])
    positions = np.tile(np.asarray((10.0, 320.0, 2.5)), (94, 1))
    metadata = [{"label": label} for label in labels]
    monkeypatch.setattr(
        outer_runner, "outer_receptors",
        lambda dataset, **kwargs: (labels, positions, metadata),
    )
    monkeypatch.setattr(
        outer_runner, "build_configuration",
        lambda *args, **kwargs: (None, closure, {}),
    )
    record, paths = outer_runner.run_outer_transient(
        inner_run_directory=inner,
        output_directory=tmp_path / "outer",
        duration_s=0.25, checkpoint_interval_s=0.125,
        cell_size_m=5.0, time_step_s=0.125,
        progress_interval_s=0.125,
    )
    assert record["all_gates_pass"]
    assert record["numerics"]["lazy_handoff_query_count"] == 2
    assert record["checkpointing"]["segment_count"] == 2
    assert paths["record_json"].is_file()
