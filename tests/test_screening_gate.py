import pytest

from degali.validation.screening_gate import evaluate_screening


class _Result:
    def __init__(self, scope, warnings=()):
        self.screening_scope = scope
        self.warnings = list(warnings)


class _Geometry:
    def __init__(self, clear):
        self.is_clear = clear


def test_qualified_clear_result_is_allowed_only_for_screening():
    decision = evaluate_screening(_Result("qualified"), obstacle_screen=_Geometry(True))
    assert decision.screening_allowed
    assert not decision.design_basis_allowed
    assert not decision.approval_allowed


def test_conditional_result_requires_explicit_opt_in():
    result = _Result("conditional", ["wind range is conditional"])
    blocked = evaluate_screening(result)
    assert not blocked.screening_allowed
    allowed = evaluate_screening(result, allow_conditional=True)
    assert allowed.screening_allowed
    with pytest.raises(ValueError, match="rejected"):
        blocked.require_screening()


def test_out_of_scope_or_obstacle_contact_is_rejected():
    assert not evaluate_screening(_Result("out_of_scope")).screening_allowed
    decision = evaluate_screening(_Result("qualified"), obstacle_screen=_Geometry(False))
    assert not decision.screening_allowed
    assert "obstacle" in " ".join(decision.reasons)


def test_unknown_result_scope_is_rejected_before_use():
    with pytest.raises(ValueError, match="screening_scope"):
        evaluate_screening(_Result("unknown"))
