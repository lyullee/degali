"""Natural flux boundary, work bookkeeping, geometry and short-step tests."""

import math
import numpy as np
import pytest
from scipy.special import roots_legendre

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_transverse_mixing import mixing
from degali.addons.reservoir_thermal import (
    ReservoirThermalMoments, ReservoirShortSegment, encode_section, decode_section,
)


@pytest.fixture
def reservoir(mixing):
    return ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")


def test_encode_decode_keeps_log_mass_density_distinct_from_mass_fraction(mixing):
    par = encode_section(mixing.state, mixing.section.thermal_width_ratio)
    state, beta = decode_section(par)
    assert par[1] == pytest.approx(math.log(state[0]*state[1]))
    assert state == pytest.approx(mixing.state, rel=1e-12)
    assert beta == pytest.approx(mixing.section.thermal_width_ratio)


def test_reduced_work_is_added_once_without_changing_other_sources(reservoir):
    m = reservoir.mixing
    a = reservoir.source_ledger()
    assert a["sources"][:4] == pytest.approx(a["original"][:4], rel=1e-14)
    assert a["sources"][4] == pytest.approx(.5*m.wind**2*a["sources"][0]+a["buoyancy_work"], abs=1e-12)
    expected = m.partition.integrate(lambda q: m.area*reservoir.axial_force_density(q, m.local(q))*m.local(q)["u"], 16, square=True)
    assert a["buoyancy_work"] == pytest.approx(expected, abs=1e-12)


def test_reservoir_energy_is_mass_inflow_energy_minus_stress_work(reservoir):
    m = reservoir.mixing
    family = m.tangent_family(reservoir.source_ledger()["sources"])
    q = np.linspace(m.partition.q0, m.partition.qmax, 31)
    d = reservoir.edge_fields(q, family)
    assert d["reservoir_enthalpy"]+d["kinetic_flux"] == pytest.approx(.5*m.wind**2*d["mass_flux"], abs=1e-10)
    # No temperature overwrite or disappearance of the unresolved edge defect.
    assert np.max(abs(d["enthalpy_gradient_defect"])) > 1.


def test_five_sources_and_weak_heat_moments_are_consistent(reservoir):
    out = reservoir.evaluate()
    assert out["family"].jacobian@out["rates"] == pytest.approx(out["ledger"]["sources"], abs=1e-7)
    assert out["weak_budget_scaled_error"] < 1e-8
    b = out["scalar_budgets"]
    assert b["reservoir_energy_outward"] == pytest.approx(-out["ledger"]["original"][4], abs=1e-8)
    assert b["mass_outward"] == pytest.approx(-out["ledger"]["sources"][0], abs=1e-10)
    assert b["species_outward"] == pytest.approx(0., abs=1e-10)
    assert out["edge_gradient_defects"]["heat"] > 1e-3


def test_direct_four_face_reservoir_flux_and_second_moment(reservoir):
    m, order = reservoir.mixing, 128
    result = reservoir.evaluate()
    family, gamma = result["family"], result["gamma"]
    at = np.array([1., gamma])
    nodes, weights = roots_legendre(order)
    ly, ln = math.sqrt(2*m.partition.q0)*np.array([m.sy, m.sn])
    outward = np.zeros(2)
    for axis in (0, 1):
        for sign in (-1., 1.):
            y = np.full_like(nodes, sign*ly) if axis == 0 else ly*nodes
            n = ln*nodes if axis == 0 else np.full_like(nodes, sign*ln)
            q = .5*((y/m.sy)**2+(n/m.sn)**2)
            h = reservoir.edge_fields(q, family)["reservoir_enthalpy"]@at
            normal = sign*(y if axis == 0 else n)*h/m.area
            dw = weights*(ln if axis == 0 else ly)
            outward += [np.sum(dw*normal), np.sum(dw*normal*(y*y+n*n))]
    expected = np.array([result["scalar_budgets"][key] for key in ("reservoir_heat_outward", "reservoir_weighted_heat_outward")])
    assert outward == pytest.approx(expected, rel=1e-5, abs=1e-7)


def test_sources_and_boundary_refine_without_changing_phase_table(reservoir):
    original = reservoir.mixing.section.phase_inverse.h_values.copy()
    coarse, fine = reservoir.evaluate(order=8), reservoir.evaluate(order=16)
    assert coarse["gamma"] == pytest.approx(fine["gamma"], rel=1e-7, abs=1e-10)
    assert coarse["rates"] == pytest.approx(fine["rates"], rel=1e-7, abs=1e-8)
    assert np.array_equal(original, reservoir.mixing.section.phase_inverse.h_values)


@pytest.mark.parametrize("ratio,work", [(0., "reduced_buoyancy_work"), (math.nan, "reduced_buoyancy_work"), (1., "none")])
def test_unsupported_physics_choices_are_rejected(mixing, ratio, work):
    with pytest.raises(ValueError):
        ReservoirThermalMoments(mixing, thermal_species_ratio=ratio, mechanical_work=work)


def test_missing_choices_and_ground_cannot_silently_enable_driver(mixing):
    with pytest.raises(TypeError):
        ReservoirThermalMoments(mixing)
    mixing.section.ground_interaction = "geometry"
    with pytest.raises(NotImplementedError):
        ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")


def test_driver_rejects_failed_local_gate(mixing, monkeypatch):
    driver = ReservoirShortSegment(mixing.section.jetplume, mixing.section.thermodynamics,
                                  thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
    monkeypatch.setattr(ReservoirThermalMoments, "evaluate", lambda *a, **k: {"valid": False})
    with pytest.raises(ValueError, match="left the valid"):
        driver.evaluate(encode_section(mixing.state, mixing.section.thermal_width_ratio))


def test_midpoint_joint_state_and_source_integration_has_second_order_convergence():
    driver = object.__new__(ReservoirShortSegment)
    def evaluate(parameters):
        rate = np.zeros(8)
        rate[0] = parameters[0]
        return dict(rates=rate, ledger={"sources": np.r_[parameters[0], np.zeros(4)]},
            moment_rate=2*parameters[0], weak_budget_scaled_error=0., edge_gradient_defects={"heat": .3},
            minimum_chi_species=1., minimum_chi_momentum=2.)
    driver.evaluate = evaluate
    initial = np.r_[1., np.zeros(7)]
    results = [driver.integrate(initial, .1, steps) for steps in (2, 4, 8)]
    errors = [abs(r["parameters"][0]-math.exp(.1)) for r in results]
    assert 3.8 < errors[0]/errors[1] < 4.1
    assert 3.8 < errors[1]/errors[2] < 4.1
    for r in results:
        assert r["cumulative_sources"][0] == pytest.approx(r["parameters"][0]-1., abs=1e-14)
        assert r["cumulative_sources"][5] == pytest.approx(2*(r["parameters"][0]-1.), abs=1e-14)
        assert r["maximum_edge_heat_defect"] == .3
    with pytest.raises(ValueError):
        driver.integrate(initial, -.1, 2)
    with pytest.raises(ValueError):
        driver.integrate(initial, .1, 2.5)


@pytest.mark.parametrize("par", [np.zeros(7), np.full(8, math.nan)])
def test_invalid_parameter_vectors_are_rejected(par):
    with pytest.raises(ValueError):
        decode_section(par)
