from tools.run_hybrid_callback_smoke import run


def test_hybrid_callback_smoke_exchanges_every_shared_clock_step(tmp_path):
    record = run(tmp_path / "smoke.json")
    assert record["status"] == "conditional_callback_ledger_smoke"
    assert record["exchanged_steps"] == record["shared_clock_steps"]
    assert record["fast_diagnostics"]["maximum_mass_residual_kg"] == 0.0
    assert record["resolved_diagnostics"]["maximum_relative_mass_residual"] == 0.0
