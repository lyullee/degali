import numpy as np
import pytest

from degali.addons.conservative_nested_handoff_3d import (
    CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA,
    conservative_overlap_remap_flux_density_2d,
    nested_outer_inlet_frame_at,
    remap_conservative_nested_handoff_3d,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)


def _integral(field, z_edges, y_edges):
    area = np.diff(z_edges)[:, None] * np.diff(y_edges)[None, :]
    return np.sum(field * area, axis=(-2, -1))


def test_overlap_remap_conserves_signed_non_aligned_face_flux():
    source_z = np.asarray((0.0, 1.0, 2.0, 4.0))
    source_y = np.asarray((-2.0, -0.5, 0.5, 2.0))
    target_z = np.asarray((0.0, 2.5, 4.0))
    target_y = np.asarray((-4.0, -2.0, -1.0, 1.0, 2.0, 4.0))
    source = np.asarray((
        (1.0, -0.5, 2.0),
        (0.2, 3.0, -0.4),
        (1.5, 0.7, 0.1),
    ))

    target = conservative_overlap_remap_flux_density_2d(
        source,
        source_z_edges_m=source_z, source_y_edges_m=source_y,
        target_z_edges_m=target_z, target_y_edges_m=target_y,
    )

    assert _integral(target, target_z, target_y) == pytest.approx(
        _integral(source, source_z, source_y), abs=1.0e-14,
    )
    assert np.all(target[:, 0] == 0.0)
    assert np.all(target[:, -1] == 0.0)


def test_nested_handoff_preserves_every_flux_at_every_time_without_clipping():
    source_z = np.linspace(0.0, 4.0, 5)
    source_y = np.linspace(-2.0, 2.0, 9)
    target_z = np.linspace(0.0, 4.0, 3)
    target_y = np.linspace(-4.0, 4.0, 5)
    time = np.asarray((0.0, 1.0, 2.0))
    zc = 0.5 * (source_z[:-1] + source_z[1:])
    yc = 0.5 * (source_y[:-1] + source_y[1:])
    shape = np.exp(-((zc[:, None] - 1.0) / 0.8) ** 2 - (yc[None, :] / 0.7) ** 2)
    species = np.stack((shape, 2.0 * shape, -0.1 * shape))

    result = remap_conservative_nested_handoff_3d(
        time_s=time,
        source_z_edges_m=source_z, source_y_edges_m=source_y,
        target_z_edges_m=target_z, target_y_edges_m=target_y,
        total_mass_flux_anomaly_kg_m2_s=0.8 * species,
        advective_species_mass_flux_kg_m2_s=species.copy(),
        species_mass_flux_kg_m2_s=species,
        advective_sensible_enthalpy_flux_w_m2=-0.9e6 * species,
        sensible_enthalpy_flux_w_m2=-1.0e6 * species,
        lateral_momentum_flux_anomaly_pa=0.2 * species,
        vertical_momentum_flux_anomaly_pa=-0.05 * species,
        source_id="synthetic nested handoff",
    )

    assert result.schema == CONSERVATIVE_NESTED_HANDOFF_3D_SCHEMA
    assert result.diagnostics["maximum_relative_integral_residual"] < 1.0e-13
    assert not result.diagnostics["clipping_or_posthoc_renormalisation"]
    assert result.diagnostics["source_species_backflow_cell_count"] > 0
    assert result.diagnostics["target_species_backflow_cell_count"] > 0
    assert _integral(
        result.species_mass_flux_kg_m2_s, target_z, target_y,
    ) == pytest.approx(_integral(species, source_z, source_y), abs=1.0e-13)


def test_nested_handoff_rejects_target_that_does_not_cover_source():
    with pytest.raises(ValueError, match="completely cover"):
        conservative_overlap_remap_flux_density_2d(
            np.ones((2, 2)),
            source_z_edges_m=np.asarray((0.0, 1.0, 2.0)),
            source_y_edges_m=np.asarray((-1.0, 0.0, 1.0)),
            target_z_edges_m=np.asarray((0.0, 2.0)),
            target_y_edges_m=np.asarray((-0.5, 0.5)),
        )


def test_outer_inlet_adds_ambient_and_linearly_interpolates_conserved_fluxes():
    edges = np.asarray((0.0, 1.0, 2.0))
    time = np.asarray((0.0, 2.0))
    pulse = np.ones((2, 2, 2))
    pulse[0] *= 0.02
    pulse[1] *= 0.06
    handoff = remap_conservative_nested_handoff_3d(
        time_s=time,
        source_z_edges_m=edges, source_y_edges_m=edges,
        target_z_edges_m=edges, target_y_edges_m=edges,
        total_mass_flux_anomaly_kg_m2_s=2.0 * pulse,
        advective_species_mass_flux_kg_m2_s=pulse,
        species_mass_flux_kg_m2_s=pulse,
        advective_sensible_enthalpy_flux_w_m2=-9_000.0 * pulse,
        sensible_enthalpy_flux_w_m2=-10_000.0 * pulse,
        lateral_momentum_flux_anomaly_pa=0.5 * pulse,
        vertical_momentum_flux_anomaly_pa=-0.25 * pulse,
        source_id="positive inlet pulse",
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )

    frame = nested_outer_inlet_frame_at(
        handoff, time_s=1.0, ambient_total_mass_flux_kg_m2_s=1.0,
        closure=closure,
    )

    assert frame.total_mass_flux_kg_m2_s == pytest.approx(
        np.full((2, 2), 1.08)
    )
    assert frame.species_mass_flux_kg_m2_s == pytest.approx(
        np.full((2, 2), 0.04)
    )
    assert frame.thermodynamic_state.species_mass_fraction == pytest.approx(
        np.full((2, 2), 0.04 / 1.08)
    )
    assert frame.diagnostics["temporal_interpolation"] == "linear"


def test_outer_eos_state_uses_advective_components_but_keeps_net_fluxes():
    edges = np.asarray((0.0, 1.0, 2.0))
    shape = (2, 2, 2)
    handoff = remap_conservative_nested_handoff_3d(
        time_s=np.asarray((0.0, 1.0)),
        source_z_edges_m=edges, source_y_edges_m=edges,
        target_z_edges_m=edges, target_y_edges_m=edges,
        total_mass_flux_anomaly_kg_m2_s=np.zeros(shape),
        advective_species_mass_flux_kg_m2_s=np.full(shape, 0.01),
        species_mass_flux_kg_m2_s=np.full(shape, 0.09),
        advective_sensible_enthalpy_flux_w_m2=np.full(shape, 100.0),
        sensible_enthalpy_flux_w_m2=np.full(shape, 900.0),
        lateral_momentum_flux_anomaly_pa=np.zeros(shape),
        vertical_momentum_flux_anomaly_pa=np.zeros(shape),
        source_id="advective-state versus net-flux test",
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )

    frame = nested_outer_inlet_frame_at(
        handoff, time_s=0.5, ambient_total_mass_flux_kg_m2_s=1.0,
        closure=closure,
    )

    assert frame.species_mass_flux_kg_m2_s == pytest.approx(0.09)
    assert frame.thermodynamic_state.species_mass_fraction == pytest.approx(0.01)
    assert frame.diagnostics["thermodynamic_state_source"] == (
        "separate_advective_species_and_enthalpy_fluxes"
    )


def test_outer_inlet_repairs_only_roundoff_sized_species_flux_bound_violation():
    edges = np.asarray((0.0, 1.0, 2.0))
    shape = (2, 2, 2)
    species = np.zeros(shape)
    species[0, 0, 0] = -1.0e-21
    handoff = remap_conservative_nested_handoff_3d(
        time_s=np.asarray((0.0, 1.0)),
        source_z_edges_m=edges, source_y_edges_m=edges,
        target_z_edges_m=edges, target_y_edges_m=edges,
        total_mass_flux_anomaly_kg_m2_s=np.zeros(shape),
        advective_species_mass_flux_kg_m2_s=np.zeros(shape),
        species_mass_flux_kg_m2_s=species,
        advective_sensible_enthalpy_flux_w_m2=np.zeros(shape),
        sensible_enthalpy_flux_w_m2=np.zeros(shape),
        lateral_momentum_flux_anomaly_pa=np.zeros(shape),
        vertical_momentum_flux_anomaly_pa=np.zeros(shape),
        source_id="roundoff-bound pulse",
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12", species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=611.0,
    )

    frame = nested_outer_inlet_frame_at(
        handoff, time_s=0.5, ambient_total_mass_flux_kg_m2_s=1.0,
        closure=closure,
    )

    assert np.min(frame.species_mass_flux_kg_m2_s) == 0.0
    assert frame.diagnostics["roundoff_admissibility_correction_cell_count"] == 1
    assert frame.diagnostics["roundoff_admissibility_correction_l1_kg_m2_s"] == pytest.approx(1.0e-21)
    assert not frame.diagnostics["physical_flux_clipping_or_renormalisation"]
