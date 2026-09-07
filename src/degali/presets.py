"""Two models in one package, chosen by name.

``degali`` is a port of DEGADIS 2.1, and it is also a place where physics
DEGADIS lacks has been added.  Those are different models and should not be
confused, so each is a named configuration rather than a set of flags a caller
has to remember to set.

.. code-block:: python

    from degali.presets import DEGADIS_21, LIQUID_HYDROGEN

    DEGADIS_21.describe()       # what it is and what it was checked against
    LIQUID_HYDROGEN.applies(rate=0.28, wind=2.5, exit_velocity=105.0)

:data:`DEGADIS_21`
    The original, bit for bit.  Every parity test runs against this and
    nothing in it has been changed: the closure is DEGADIS's own, the
    thermodynamics are the 1989 correlations, and the numerics reproduce
    ``RKGST`` including its approximations.

:data:`LIQUID_HYDROGEN`
    The same code with the additions hydrogen needs: equations of state rather
    than the 1989 correlations, a buoyant closure so a cloud that becomes
    lighter than air can leave the ground, and the orifice-plane source for a
    flashing release.

Neither is a default that silently applies.  A caller names one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .evidence import NEAR_FIELD_CONCENTRATION, NEAR_FIELD_CORRECTED, RANGE


@dataclass(frozen=True)
class Applicability:
    """Where a configuration has been checked, and against what."""

    summary: str
    limits: dict = field(default_factory=dict)
    evidence: str = ""

    def check(self, **values) -> list[str]:
        """Warnings for values outside the checked range."""
        out = []
        for name, value in values.items():
            span = self.limits.get(name)
            if span is None or value is None:
                continue
            lo, hi = span
            if not (lo <= value <= hi):
                out.append(
                    f"{name} {value:g} is outside the checked range "
                    f"{lo:g} to {hi:g}"
                )
        return out


@dataclass(frozen=True)
class Preset:
    """A named model configuration."""

    name: str
    backend: str
    legacy_numerics: bool
    closure: str  #: ``"degadis"``, ``"buoyant"`` or ``"unified"``
    source: str  #: ``"deck"`` or ``"flashing"``
    applicability: Applicability
    #: Corrections applied to the jet path: published relations omitted by the
    #: shipped model plus conservation fixes at an inconsistent source plane.
    jet_corrections: tuple = ()

    def describe(self) -> str:
        parts = [
            f"{self.name}",
            f"  thermodynamics : {self.backend}",
            f"  numerics       : {'1989, bit-faithful' if self.legacy_numerics else 'modern'}",
            f"  ground closure : {self.closure}",
            f"  source term    : {self.source}",
            f"  checked        : {self.applicability.summary}",
        ]
        if self.jet_corrections:
            parts.append("  jet path       :")
            parts += [f"    - {c}" for c in self.jet_corrections]
        if self.applicability.evidence:
            parts.append(f"  evidence       : {self.applicability.evidence}")
        return "\n".join(parts)

    def applies(self, **values) -> list[str]:
        return self.applicability.check(**values)


#: The original model, unchanged.
DEGADIS_21 = Preset(
    name="DEGADIS 2.1",
    backend="legacy",
    legacy_numerics=True,
    closure="degadis",
    source="deck",
    applicability=Applicability(
        summary=(
            "reproduces the 1989 Fortran to 1e-12 on all five EPA test cases; "
            "evaluated against LNG pool spills and ammonia jets"
        ),
        limits={},
        evidence=(
            "Burro, lowest instrumented height: MG 0.81, 95 % CI [0.63, 1.02], "
            "FAC2 0.56, n=61. Over every height the vertical structure fails "
            "(MG 7.3 at 3 m, 5500 at 8 m), and the lateral spread runs one to "
            "three times wide."
        ),
    ),
)

#: The hydrogen configuration.
#:
#: The applicability limits are not decoration.  A release whose exit velocity
#: is comparable with the wind is steered by the wind rather than by its own
#: momentum, and a steady jet model does not describe it: on the PRESLHY 1 barg
#: trials, nominally identical releases gave arc maxima of 83 % and 4 %.
LIQUID_HYDROGEN = Preset(
    name="liquid hydrogen",
    backend="coolprop",
    legacy_numerics=False,
    closure="buoyant",
    source="flashing",
    jet_corrections=(
        "start the plume at the expanded plane, not the orifice "
        "(EPA-450/4-90-018 §4.3: the source area should be the cross-section "
        "of the fully expanded jet)",
        "scale the shear entrainment as sqrt(rho_a/rho_jet) "
        "(Ricou & Spalding 1961; Panda & Hecht for cryogenic hydrogen)",
        "use the measured plume entrainment coefficient, 0.0875, not the jet "
        "one (Papanicolaou & List 1988: an LH2 release is buoyant within a metre)",
        "preserve the flash-entrained air by rebuilding the thermodynamic "
        "table from the mixed equivalent-source state (mass and enthalpy "
        "conservation at the ODE boundary)",
        "conserve total momentum while initially stationary air is entrained "
        "(pressure thrust belongs to the preceding notional-nozzle zone)",
    ),
    applicability=Applicability(
        summary=(
            "momentum-dominated flashing releases, 0.79 to 6 m downwind; "
            "buoyancy regime checked against the NASA spills"
        ),
        limits={
            **RANGE["jet"],
            "velocity_ratio": (10.0, 1.0e6),  # exit speed over wind speed
        },
        evidence=(
            "PRESLHY E3.5 raw data, uncensored 0-100 vol % sensors, "
            "recomputed from the reduction rather than quoted. Arc maxima "
            "against the model evaluated at the sensor locations over 69 "
            "historical arcs and 62 arcs downstream of the corrected source "
            "plane from nine trials give " + NEAR_FIELD_CONCENTRATION.cite() +
            " as shipped and " + NEAR_FIELD_CORRECTED.cite() +
            " with the corrections. Corrected VG and FAC2 pass Hanna's "
            "thresholds; MG is close to unity. Measured separately "
            "from Gaussian fits using the same standard-deviation convention "
            "as JETPLU, the corrections take the vertical spread from 0.901 "
            "of the measured value to 1.033 over 23 "
            "momentum-driven fits. Mean signed centre-height error is 0.036 m. At "
            "30--100 m, Spadeadam test 4 remains low while test 6 detaches, "
            "matching the observed regime split; forcing ground contact is "
            "therefore not part of the preset. Pool buoyancy regime: 4 of 4 "
            "on the NASA spills."
        ),
    ),
)

PRESETS = {p.name: p for p in (DEGADIS_21, LIQUID_HYDROGEN)}
