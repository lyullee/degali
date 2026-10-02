import pytest

from degali.validation.integration import (
    ValidationEvidence,
    assess_validation_portfolio,
    validation_portfolio,
    validation_portfolio_dict,
)


def test_current_portfolio_keeps_reconstruction_field_screen_and_research_distinct():
    decisions = {item.branch: item for item in assess_validation_portfolio(validation_portfolio())}
    assert decisions["legacy_degadis_reconstruction"].classification == "reconstruction_only"
    assert decisions["fast_single_velocity_lh2"].classification == "qualified_field_screen"
    assert decisions["finite_tke_transport"].classification == "research_only"
    assert decisions["mixed_phase_droplet_transport"].classification == "research_only"
    assert all(not item.default_promotion_allowed for item in decisions.values())
    assert all(not item.field_accuracy_claim_allowed for item in decisions.values())


def test_only_independent_uncalibrated_observed_experiment_can_qualify_a_field_screen():
    common = dict(branch="branch", component="component", campaign="campaign",
                  evidence_kind="experimental", outcome="passed", limitation="none")
    blocked = ValidationEvidence("blocked", **common, uncalibrated_prediction=False,
                                 independent_experiment=True, decisive_state_observed=True)
    accepted = ValidationEvidence("accepted", **common, uncalibrated_prediction=True,
                                  independent_experiment=True, decisive_state_observed=True)
    assert assess_validation_portfolio([blocked])[0].classification == "research_only"
    assert assess_validation_portfolio([accepted])[0].classification == "qualified_field_screen"


def test_portfolio_rejects_duplicate_identifiers_and_invalid_evidence_claims():
    row = validation_portfolio()[0]
    with pytest.raises(ValueError, match="unique"):
        assess_validation_portfolio([row, row])
    with pytest.raises(ValueError, match="only experimental"):
        ValidationEvidence("bad", "branch", "component", "campaign", "simulation", "passed",
                           True, True, True, "none")


def test_json_ready_portfolio_cannot_claim_a_composite_score_or_default_promotion():
    report = validation_portfolio_dict()
    assert not report["third_party_observations_distributed"]
    assert not report["composite_accuracy_score_calculated"]
    assert not report["automatic_default_promotion_allowed"]
    assert len(report["evidence"]) == 9
