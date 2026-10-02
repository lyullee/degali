"""Operational gate between an LH2 result and an engineering screen.

The gate is intentionally stricter than the numerical solver.  It lets a
workflow use a qualified or explicitly accepted conditional result for
screening while refusing out-of-scope results and any trajectory that enters
declared solid geometry.  No result is ever promoted to a design basis or an
approval decision by this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ScreeningDecision:
    """Auditable disposition of one calculated scenario."""

    screening_allowed: bool
    design_basis_allowed: bool
    scope: str
    geometry_clear: bool | None
    reasons: tuple[str, ...]

    @property
    def approval_allowed(self) -> bool:
        """Always false: this package is not an approval/permitting tool."""
        return False

    def require_screening(self) -> None:
        """Raise when a caller attempts to consume a disallowed result."""
        if not self.screening_allowed:
            raise ValueError("LH2 screening gate rejected the result: " + "; ".join(self.reasons))


def evaluate_screening(
    result: Any,
    *,
    obstacle_screen: Any | None = None,
    allow_conditional: bool = False,
) -> ScreeningDecision:
    """Evaluate a steady or time-resolved result for engineering screening.

    ``result`` must expose ``screening_scope`` and may expose ``warnings``.
    ``obstacle_screen`` is the return value of
    :func:`degali.addons.site_geometry.screen_trajectory`; any encounter
    rejects downstream free-plume screening.  Conditional physics can be
    accepted only explicitly, and even then this function never grants a
    design-basis or permitting status.
    """
    scope = getattr(result, "screening_scope", None)
    if scope not in {"qualified", "conditional", "out_of_scope"}:
        raise ValueError("result must expose a valid screening_scope")
    warnings = tuple(str(item) for item in getattr(result, "warnings", ()) or ())
    reasons: list[str] = []
    allowed = scope == "qualified" or (scope == "conditional" and allow_conditional)
    if scope == "out_of_scope":
        reasons.append("calculation is outside the validated physical envelope")
    elif scope == "conditional" and not allow_conditional:
        reasons.append("conditional result requires allow_conditional=True")
    if warnings:
        reasons.append("; ".join(warnings))

    geometry_clear: bool | None = None
    if obstacle_screen is not None:
        geometry_clear = bool(getattr(obstacle_screen, "is_clear", False))
        if not geometry_clear:
            allowed = False
            reasons.append("declared obstacle/wall contact requires an obstacle-resolved model")

    if not reasons and not allowed:
        reasons.append("screening criteria were not satisfied")
    return ScreeningDecision(
        screening_allowed=bool(allowed),
        design_basis_allowed=False,
        scope=scope,
        geometry_clear=geometry_clear,
        reasons=tuple(dict.fromkeys(reasons)),
    )


__all__ = ["ScreeningDecision", "evaluate_screening"]
