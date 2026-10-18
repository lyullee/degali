from types import SimpleNamespace

from degali import cli
from degali.core import transient


def test_observer_failure_is_retained_with_stage_and_exception(monkeypatch):
    run = transient.TransientRun.__new__(transient.TransientRun)
    run.kin = SimpleNamespace(position=lambda _time, _t0: 0.0)
    run.v = object()
    run.h = {"tend": 1.0}
    run.nobs = 1
    run._crossings = lambda _t0: (0.0, 1.0, True, True)

    def fail_collection(*_args):
        raise ValueError("synthetic observer failure")

    run._collect = fail_collection
    monkeypatch.setattr(transient, "release_times", lambda *_args: [0.25])

    result = run.run()

    assert len(result.observers) == 1
    assert result.failures == result.observers
    failure = result.failures[0]
    assert failure.rows.shape == (0, 10)
    assert failure.failure_stage == "collection"
    assert failure.failure_type == "ValueError"
    assert failure.failure_message == "synthetic observer failure"


def test_empty_observer_without_exception_is_not_a_failure(monkeypatch):
    run = transient.TransientRun.__new__(transient.TransientRun)
    run.kin = SimpleNamespace(position=lambda _time, _t0: 0.0)
    run.v = object()
    run.h = {"tend": 1.0}
    run.nobs = 1
    run._crossings = lambda _t0: (0.0, 1.0, True, True)
    run._collect = lambda *_args: (None, None, None)
    monkeypatch.setattr(transient, "release_times", lambda *_args: [0.25])

    result = run.run()

    assert len(result.observers) == 1
    assert result.failures == []


def test_cli_warns_that_failed_observer_is_not_a_zero_prediction(capsys):
    failure = transient.ObserverResult(
        index=3,
        t0=0.0,
        tup=0.0,
        tdown=1.0,
        xup=0.0,
        xdown=1.0,
        erate=0.0,
        b=0.0,
        sz0=0.0,
        cc0=0.0,
        failure_stage="downwind_profile",
        failure_type="OverflowError",
        failure_message="synthetic overflow",
    )

    cli._warn_observer_failures(transient.TransientResult([failure]))

    error = capsys.readouterr().err
    assert "not zero-concentration predictions" in error
    assert "observer 3: downwind_profile: OverflowError" in error
