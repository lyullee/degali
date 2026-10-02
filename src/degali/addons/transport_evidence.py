"""Evidence boundaries for unresolved turbulence and condensed-particle transport.

These records distinguish a conservation-verified numerical kernel from a
validated physical closure.  They intentionally contain no tunable number and
cannot promote a research calculation into a quantitative LH2 prediction.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class TransportEvidenceBoundary:
    """Immutable status of a transport mechanism awaiting identifying data."""

    identifier: str
    numerical_boundary_verified: bool
    physical_closure_validated: bool
    default_prediction_enabled: bool
    quantitative_lh2_prediction_allowed: bool
    required_measurements: tuple[str, ...]
    admitted_use: str

    def audit_fields(self) -> dict[str, bool | str]:
        """Return serialisable result fields without weakening this boundary."""
        return {
            "transport_evidence_boundary": self.identifier,
            "numerical_boundary_verified": self.numerical_boundary_verified,
            "physical_closure_validated": self.physical_closure_validated,
            "default_prediction_enabled": self.default_prediction_enabled,
            "quantitative_lh2_prediction_allowed": (
                self.quantitative_lh2_prediction_allowed
            ),
        }


TURBULENCE_BOUNDARY = TransportEvidenceBoundary(
    identifier="turbulence_conservation_and_realizability_boundary",
    numerical_boundary_verified=True,
    physical_closure_validated=False,
    default_prediction_enabled=False,
    quantitative_lh2_prediction_allowed=False,
    required_measurements=(
        "gas-phase velocity-fluctuation profile",
        "Reynolds-stress information or an explicitly bounded subset",
        "dissipation rate or an independently measured length/time scale",
        "matched thermal/species observations at a separate validation location",
    ),
    admitted_use=(
        "explicit-input conservation, positive-semidefinite realizability, "
        "and predeclared sensitivity rejection only"
    ),
)


PARTICLE_SLIP_BOUNDARY = TransportEvidenceBoundary(
    identifier="particle_slip_conservation_and_kinematic_limit_boundary",
    numerical_boundary_verified=True,
    physical_closure_validated=False,
    default_prediction_enabled=False,
    quantitative_lh2_prediction_allowed=False,
    required_measurements=(
        "condensed N2/O2 mass fraction or phase inventory",
        "mass-weighted airborne particle-size distribution",
        "co-located gas and particle velocity measurements",
        "gas properties and measurement locations/time windows for the drag law",
    ),
    admitted_use=(
        "explicit-input drag conservation plus no-slip and stationary-condensate "
        "kinematic limits only"
    ),
)
