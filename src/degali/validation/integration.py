"""Claim-safe integration of independent DEGALI validation evidence.

This module deliberately does not pool MG/VG/FAC2 values across campaigns.
They describe different source boundaries, sensors, averaging windows and
applicability domains.  It instead records what each evidence item can and
cannot establish, then makes the promotion boundary executable:

* a legacy oracle establishes reconstruction only;
* a numerical LES or manufactured test can bound/check an implementation but
  cannot establish experimental accuracy;
* an experimental screen must be an uncalibrated prediction on the declared
  model branch before it can support a *qualified* field claim; and
* a new physical branch remains research-only when its decisive state (for
  example particle inventory or epsilon) is unmeasured.

The fixed portfolio contains provenance and status, not external observations.
It is safe to distribute and can be regenerated without controlled data.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


_EVIDENCE_KINDS = {"oracle", "experimental", "simulation", "manufactured"}
_OUTCOMES = {"passed", "qualified", "provisional", "rejected", "unmeasured"}


@dataclass(frozen=True)
class ValidationEvidence:
    """One non-interchangeable validation item and its claim boundary."""

    identifier: str
    branch: str
    component: str
    campaign: str
    evidence_kind: str
    outcome: str
    uncalibrated_prediction: bool
    independent_experiment: bool
    decisive_state_observed: bool
    limitation: str

    def __post_init__(self):
        if not self.identifier or not self.branch or not self.component or not self.campaign:
            raise ValueError("identifier, branch, component and campaign are required")
        if self.evidence_kind not in _EVIDENCE_KINDS:
            raise ValueError(f"unknown evidence kind: {self.evidence_kind!r}")
        if self.outcome not in _OUTCOMES:
            raise ValueError(f"unknown validation outcome: {self.outcome!r}")
        if self.evidence_kind != "experimental" and self.independent_experiment:
            raise ValueError("only experimental evidence can be marked independent")
        if self.outcome == "unmeasured" and self.decisive_state_observed:
            raise ValueError("an unmeasured result cannot claim its decisive state was observed")


@dataclass(frozen=True)
class BranchDecision:
    """Integrated claim class for one branch; never an accuracy score."""

    branch: str
    classification: str
    experimental_items: tuple[str, ...]
    blockers: tuple[str, ...]
    default_promotion_allowed: bool = False
    field_accuracy_claim_allowed: bool = False


def assess_validation_portfolio(
    evidence: Iterable[ValidationEvidence],
) -> tuple[BranchDecision, ...]:
    """Classify branches without mixing their errors or inventing weights.

    ``qualified_field_screen`` means a declared branch has at least one
    independent, uncalibrated experimental result with its decisive state
    observed. It is intentionally weaker than a universal-accuracy claim and
    is never sufficient for automatic default promotion. All other branch
    classes make their missing evidence explicit.
    """
    items = tuple(evidence)
    if not items:
        raise ValueError("at least one validation evidence item is required")
    identifiers = [item.identifier for item in items]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("validation evidence identifiers must be unique")

    decisions = []
    for branch in dict.fromkeys(item.branch for item in items):
        rows = tuple(item for item in items if item.branch == branch)
        experimental = tuple(item for item in rows if item.evidence_kind == "experimental")
        qualifying = tuple(
            item for item in experimental
            if item.outcome in {"passed", "qualified"}
            and item.uncalibrated_prediction
            and item.independent_experiment
            and item.decisive_state_observed
        )
        blockers = []
        if not experimental:
            blockers.append("no experimental evidence on this branch")
        if experimental and not qualifying:
            blockers.append("no independent uncalibrated experimental screen observes the decisive state")
        for row in rows:
            if row.outcome in {"rejected", "unmeasured"}:
                blockers.append(f"{row.identifier}: {row.limitation}")
        if qualifying:
            classification = "qualified_field_screen"
        elif all(row.evidence_kind == "oracle" for row in rows):
            classification = "reconstruction_only"
        else:
            classification = "research_only"
        decisions.append(BranchDecision(
            branch=branch,
            classification=classification,
            experimental_items=tuple(row.identifier for row in experimental),
            blockers=tuple(dict.fromkeys(blockers)),
        ))
    return tuple(decisions)


def validation_portfolio() -> tuple[ValidationEvidence, ...]:
    """Return the frozen 2026-09-17 LH2 validation integration register.

    The register names only evidence already documented in the repository.
    It does not expose a third-party observation, value or figure.
    """
    return (
        ValidationEvidence(
            "degadis_fortran_oracle", "legacy_degadis_reconstruction",
            "original DEGADIS state reproduction", "controlled source-built oracle",
            "oracle", "passed", True, False, True,
            "reproduction is not a liquid-hydrogen field-accuracy test",
        ),
        ValidationEvidence(
            "hecht_panda_mean_profiles", "cryogenic_gas_axisymmetric",
            "near-field mean species and thermal profiles", "Sandia Raman free jets",
            "experimental", "provisional", True, True, True,
            "published aggregate fits lack the environmental metadata needed for a final predictive claim",
        ),
        ValidationEvidence(
            "preslhy_e35_mean_arcs", "fast_single_velocity_lh2",
            "outdoor concentration/trajectory screen", "PRESLHY E3.5 horizontal releases",
            "experimental", "qualified", True, True, True,
            "sensor ceiling, wind meander and uncertain cold phase inventory limit the scope",
        ),
        ValidationEvidence(
            "spadeadam_reported_arcs", "fast_single_velocity_lh2",
            "independent downstream concentration/trajectory screen", "DNV/FFI Spadeadam releases",
            "experimental", "qualified", True, True, True,
            "reported arc maxima and mast winds do not resolve transient plume meander",
        ),
        ValidationEvidence(
            "preslhy_e31_blowdown", "cryocompressed_source",
            "vessel/nozzle pressure history", "PRESLHY E3.1 80 K release",
            "experimental", "qualified", True, True, True,
            "unmeasured wall, withdrawal and two-phase pipe states prevent a universal source closure",
        ),
        ValidationEvidence(
            "two_normal_rms_les", "finite_tke_transport",
            "TKE positive-semidefinite lower-bound screen", "published cryogenic-hydrogen LES",
            "simulation", "provisional", True, False, False,
            "no matched experimental velocity RMS, epsilon or integral-length profile",
        ),
        ValidationEvidence(
            "sandia_particle_piv", "finite_tke_transport",
            "mean particle-PIV velocity capability", "Sandia PIV/Raman presentation",
            "experimental", "unmeasured", False, True, False,
            "particle tracer response and velocity fluctuation arrays are not public gas-TKE data",
        ),
        ValidationEvidence(
            "finite_rate_droplet_bound", "mixed_phase_droplet_transport",
            "droplet evaporation/relaxation limiting bounds", "PRESLHY and Spadeadam source states",
            "manufactured", "unmeasured", False, False, False,
            "no airborne particle-size, phase-fraction or gas-particle-slip observation",
        ),
        ValidationEvidence(
            "li2026_equilibrium_limit", "equilibrium_air_condensation",
            "equilibrium air-condensation limiting comparison", "Li et al. integral model",
            "simulation", "rejected", False, False, False,
            "sub-triple nitrogen liquid assumption and missing N2/O2 mixed-phase closure",
        ),
    )


def validation_portfolio_dict() -> dict:
    """JSON-ready portfolio and non-promoting decisions for audit tools."""
    evidence = validation_portfolio()
    decisions = assess_validation_portfolio(evidence)
    return {
        "portfolio_date": "2026-09-17",
        "third_party_observations_distributed": False,
        "composite_accuracy_score_calculated": False,
        "automatic_default_promotion_allowed": False,
        "evidence": [asdict(item) for item in evidence],
        "decisions": [asdict(item) for item in decisions],
    }


__all__ = [
    "ValidationEvidence", "BranchDecision", "assess_validation_portfolio",
    "validation_portfolio", "validation_portfolio_dict",
]
