from degali.addons.hybrid_amr import make_obstacle_patch_bounds
from degali.addons.hybrid_amr_plan import HybridAMRRunPlan
from degali.addons.hybrid_amr_executor import HybridAMRExecutor
from degali.addons.hybrid_amr_executor import run_solver_with_patch_flux_collector
from degali.addons.hybrid_amr_executor import prepare_patch_flux_collector
from degali.addons.hybrid_amr_executor import run_solver_with_prepared_patch_flux
from degali.addons.hybrid_conservative_exchange import (
    ConservativeControlVolumeState,
    ConservativeInterfacePacket,
    conserved_totals,
    exchange_conservative_pair,
)
from degali.addons.fast_obstacle_transport_3d import FastObstacleTransport3DConfig, solve_fast_obstacle_transport_3d
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid

def test_shared_clock_preserves_total_species_mass():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
                                       downstream_m=4, lateral_m=3, vertical_m=2,
                                       domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 2, 1, .25)
    result = HybridAMRExecutor(plan).run(
        coarse_mass_kg=10, patch_mass_kg=1,
        request=lambda _t, _c, _p: (0.2, 0.05))
    assert result.patch_steps == 8
    assert abs(result.coarse_mass_kg + result.patch_mass_kg - 11) < 1e-12
    assert abs(result.interface_residual_kg) <= 1e-12

def test_flux_schedule_drives_dynamic_exchange_on_patch_clock():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
                                       downstream_m=4, lateral_m=3, vertical_m=2,
                                       domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    result = HybridAMRExecutor(plan).run_with_flux_rates(
        coarse_mass_kg=2, patch_mass_kg=0,
        time_s=[0, .25, .5, .75, 1.0],
        coarse_to_patch_kg_s=[1, 1, 1, 1, 1],
        patch_to_coarse_kg_s=[0, 0, 0, 0, 0])
    assert result.patch_mass_kg == 1
    assert abs(result.coarse_mass_kg + result.patch_mass_kg - 2) < 1e-12

def test_solver_wrapper_requires_explicit_flux_contract():
    collector = lambda *_args: None
    assert run_solver_with_patch_flux_collector(
        lambda config, face_flux_callback: (config, face_flux_callback),
        "cfg", collector)[0] == "cfg"
    with pytest.raises(TypeError, match="does not implement"):
        run_solver_with_patch_flux_collector(lambda config: config, "cfg", collector)

def test_prepare_collector_rejects_misaligned_patch():
    class Config:
        domain_m = (4., 2.4, 2.)
        cells = (16, 10, 8)
    patch = make_obstacle_patch_bounds((1, 2, .8, 1.6, 0, 1), upstream_m=.01,
        downstream_m=.01, lateral_m=.01, vertical_m=.01, domain_m=Config.domain_m)
    with pytest.raises(ValueError, match="solver face"):
        prepare_patch_flux_collector(Config(), patch)

def test_prepare_collector_accepts_aligned_patch():
    class Config:
        domain_m = (4., 2.4, 2.)
        cells = (16, 10, 8)
    patch = make_obstacle_patch_bounds((1.5, 2., .72, 1.68, 0., 1.),
        upstream_m=0, downstream_m=0, lateral_m=0, vertical_m=0,
        domain_m=Config.domain_m)
    collector = prepare_patch_flux_collector(Config(), patch)
    assert collector.patch == patch

def test_prepared_flux_runner_returns_result_and_collector():
    class Config:
        domain_m = (4., 2.4, 2.)
        cells = (16, 10, 8)
    patch = make_obstacle_patch_bounds((1.5, 2., .72, 1.68, 0., 1.),
        upstream_m=0, downstream_m=0, lateral_m=0, vertical_m=0,
        domain_m=Config.domain_m)
    result, collector = run_solver_with_prepared_patch_flux(
        lambda config, face_flux_callback: ("ok", face_flux_callback), Config(), patch)
    assert result[0] == "ok"
    assert collector.patch == patch

def test_bidirectional_runner_passes_previous_interface_mass_to_callbacks():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
                                       downstream_m=4, lateral_m=3, vertical_m=2,
                                       domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    seen = []
    def step(state, t, dt, incoming):
        seen.append(incoming)
        return state, 0.0
    _, _, history = HybridAMRExecutor(plan).run_bidirectional(
        coarse_state={}, patch_state={}, coarse_step=step, patch_step=step,
        coarse_mass_kg=1.0, patch_mass_kg=0.0)
    assert len(history) == 4
    assert seen == [0.0] * 8

def test_bidirectional_runner_limits_transfer_to_initial_donor_mass():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
        downstream_m=4, lateral_m=3, vertical_m=2, domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    def coarse(state, t, dt, incoming): return state, 10.0
    def local(state, t, dt, incoming): return state, 0.0
    _, _, history = HybridAMRExecutor(plan).run_bidirectional(
        coarse_state={}, patch_state={}, coarse_step=coarse, patch_step=local,
        coarse_mass_kg=1.0, patch_mass_kg=0.0)
    assert history[0][1] == 1.0

def test_bidirectional_runner_delivers_equal_opposite_mass_next_substep():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
        downstream_m=4, lateral_m=3, vertical_m=2, domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 2, 1, .25)
    coarse_seen, patch_seen = [], []
    def coarse(state, t, dt, incoming):
        coarse_seen.append(incoming)
        return state, 0.5 if t == 0 else 0.0
    def local(state, t, dt, incoming):
        patch_seen.append(incoming)
        return state, 0.0
    HybridAMRExecutor(plan).run_bidirectional(
        coarse_state={}, patch_state={}, coarse_step=coarse, patch_step=local,
        coarse_mass_kg=1.0, patch_mass_kg=0.0)
    assert patch_seen[1] == 0.5

def test_transactional_runner_keeps_callback_state_contract_bounded():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
        downstream_m=4, lateral_m=3, vertical_m=2, domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    seen = []
    def coarse(state, t, dt, incoming, available):
        seen.append(available)
        return state, min(.1, available)
    def local(state, t, dt, incoming, available):
        return state, 0.0
    _, _, history = HybridAMRExecutor(plan).run_bidirectional_transactional(
        coarse_state={}, patch_state={}, coarse_step=coarse, patch_step=local,
        coarse_mass_kg=.15, patch_mass_kg=0.)
    assert history[-1][3] + history[-1][4] == pytest.approx(.15)
    assert seen[0] == .15

def test_transactional_runner_rejects_outgoing_above_available_mass():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
        downstream_m=4, lateral_m=3, vertical_m=2, domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    def bad(state, t, dt, incoming, available): return state, available + 1
    with pytest.raises(ValueError, match="above available"):
        HybridAMRExecutor(plan).run_bidirectional_transactional(
            coarse_state={}, patch_state={}, coarse_step=bad, patch_step=bad,
            coarse_mass_kg=.1, patch_mass_kg=0.)


def test_conservative_field_exchange_preserves_species_energy_and_momentum():
    coarse = ConservativeControlVolumeState(2.0, 0.5, 10.0, (1.0, 0.0, 0.0))
    patch = ConservativeControlVolumeState(1.0, 0.2, -2.0, (0.0, 2.0, 0.0))
    initial = conserved_totals(coarse, patch)
    coarse, patch = exchange_conservative_pair(
        coarse, patch,
        coarse_to_patch=ConservativeInterfacePacket(.25, .1, 1.5, (.2, .0, .0)),
        patch_to_coarse=ConservativeInterfacePacket(.1, .02, -.5, (.0, .1, .0)),
    )
    final = conserved_totals(coarse, patch)
    assert final.total_mass_kg == pytest.approx(initial.total_mass_kg)
    assert final.species_mass_kg == pytest.approx(initial.species_mass_kg)
    assert final.sensible_enthalpy_j == pytest.approx(initial.sensible_enthalpy_j)
    assert final.momentum_kg_m_s == pytest.approx(initial.momentum_kg_m_s)


def test_conservative_executor_runs_shared_clock_without_posthoc_clipping():
    patch = make_obstacle_patch_bounds((10, 20, 5, 15, 0, 8), upstream_m=4,
                                       downstream_m=4, lateral_m=3, vertical_m=2,
                                       domain_m=(30, 25, 20))
    plan = HybridAMRRunPlan((30, 25, 20), 5, 1, patch, 1, 1, .25)
    executor = HybridAMRExecutor(plan)
    initial_coarse = ConservativeControlVolumeState(1.0, .4, 3.0, (1., 0., 0.))
    initial_patch = ConservativeControlVolumeState(.5, .1, 2.0, (0., 1., 0.))

    def coarse_step(state, _t, _dt, _incoming):
        return state, ConservativeInterfacePacket(.05, .02, .1, (.01, 0., 0.))

    def patch_step(state, _t, _dt, _incoming):
        return state, ConservativeInterfacePacket(.01, .002, -.02, (0., .01, 0.))

    coarse, local, history = executor.run_bidirectional_conservative(
        coarse_state=initial_coarse, patch_state=initial_patch,
        coarse_step=coarse_step, patch_step=patch_step,
    )
    assert len(history) == 4
    assert conserved_totals(coarse, local).total_mass_kg == pytest.approx(1.5)
    assert conserved_totals(coarse, local).species_mass_kg == pytest.approx(.5)


def test_conservative_exchange_rejects_species_above_donor():
    coarse = ConservativeControlVolumeState(.1, .01)
    patch = ConservativeControlVolumeState(.1, .01)
    with pytest.raises(ValueError, match="exceeds donor species"):
        exchange_conservative_pair(
            coarse, patch,
            coarse_to_patch=ConservativeInterfacePacket(.02, .02),
            patch_to_coarse=ConservativeInterfacePacket.zero(),
        )

def test_prepared_wrapper_runs_real_fast_obstacle_solver_flux_callback():
    config = FastObstacleTransport3DConfig(
        domain_m=(4., 2.4, 2.), cells=(16, 10, 8), duration_s=.001,
        time_step_s=.001, wind_speed_m_s=2., diffusivity_m2_s=.03,
        source_schedule=SourceRateSchedule((0., .001), (0., 0.)),
        source_position_m=(.5, 1.2, .4), source_sigma_m=.15,
        obstacles=(AxisAlignedCuboid(1.5, 2., .8, 1.6, 0., 1., "block"),))
    patch = make_obstacle_patch_bounds((1.5, 2., .72, 1.68, 0., 1.),
        upstream_m=0, downstream_m=0, lateral_m=0, vertical_m=0,
        domain_m=config.domain_m)
    result, collector = run_solver_with_prepared_patch_flux(
        solve_fast_obstacle_transport_3d, config, patch)
    assert result.time_s[-1] == .001
    assert len(collector.rates_kg_s) == 1
import pytest

