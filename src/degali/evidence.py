"""Every validated number, in one place.

A review of this package found the same quantity reported with four different
values across the README, the handover, a test assertion and the code that
prints it to a user: the neutral-buoyancy concentration was 85-90 mol % in one
document, 99.9 in another, ``> 0.99`` in a test and 100.0 in the program's own
output.  The near-field concentration statistic appeared as ``MG 0.74, VG 1.22,
FAC2 0.92`` over 53 arcs in the README and in a hard-coded f-string, and as
``MG 0.738, VG 1.41, FAC2 0.83`` over 69 arcs everywhere else. Those historical
configurations remain reproducible, while the current conserved-source claim
is held separately below.

That is a documentation defect rather than a modelling one, but it is the kind
that survives review, because there is no single place a reader can check.  So
there is one now.  Anything that reports a validated number -- the one-call
API, the presets, the test suite, the generated tables in ``docs/`` -- reads it
from here, and :func:`degali.evidence.check_documents` fails the suite if a
retired figure reappears in prose.

Each entry carries the grade from ``docs/claim-grading.md``:

``A``
    deterministic; reproducible exactly, and sample size does not enter.
``B``
    statistical, with an interval resampled over *trials* rather than points.
``C``
    directional only.  Consistent-with, never asserted.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Claim:
    """One validated number, with everything needed to cite it."""

    #: What was measured, in words.
    quantity: str
    #: The headline value, formatted for display.
    value: str
    #: Sample size, or ``None`` for a deterministic result.
    n: int | None
    #: ``"A"``, ``"B"`` or ``"C"``.
    grade: str
    #: The dataset or document it comes from.
    source: str
    #: Supporting statistics, keyed by name.
    detail: dict = field(default_factory=dict)

    def cite(self) -> str:
        """A one-line citation, for printing next to a computed answer."""
        bits = [self.value]
        if self.detail:
            bits += [f"{k} {v}" for k, v in self.detail.items()]
        if self.n is not None:
            bits.append(f"n={self.n}")
        return ", ".join(bits)


# ==========================================================================
# the claims
# ==========================================================================

#: Near-field concentration on momentum-driven flashing jets.
#:
#: Arc maxima against the model centreline.  The interval is bootstrapped over
#: trials, not points: readings within one trial share a release, a wind and a
#: source estimate.  Point-level resampling would give [0.68, 0.80] and
#: overstate the evidence by roughly the square root of the sensor count.
NEAR_FIELD_CONCENTRATION = Claim(
    quantity="near-field concentration, momentum-driven releases",
    value="MG 1.121",
    n=69,
    grade="B",
    source="PRESLHY E3.5 (DOI 10.35097/1481), 9 trials, uncensored, arc "
           "maximum against the model evaluated at each sensor -- the "
           "convention Burro already used. Against the centreline the same "
           "arcs give MG 0.729, VG 1.42, FAC2 0.83, which "
           "reproduces the published figures and validated the rebuild.",
    detail={"95 % CI": "[0.736, 2.264]", "VG": "16.85", "FAC2": "0.80"},
)

#: Corrected source, restricted to arcs downstream of its establishment plane.
NEAR_FIELD_CORRECTED = Claim(
    quantity="near-field concentration, corrections applied",
    value="MG 1.047",
    n=62,
    grade="B",
    source="as above; the five adopted jet corrections switched on together; "
           "seven 0.35/0.53 m arcs upstream of the conserved equivalent-"
           "source plane are outside the established-plume model",
    detail={"95 % CI": "[0.759, 1.402]", "VG": "1.42", "FAC2": "0.84"},
)

#: Lowest flammable height on the NASA Langley pool spills.
#:
#: Four trials, and three of the four "measured" values are inversions of
#: thermocouple data through a mixing model rather than concentration
#: measurements.  Grade C: the *regime* is right 4 of 4, which is strong for
#: what it is; the magnitude is not established.
LIFTOFF_HEIGHT = Claim(
    quantity="lowest flammable height at the 33.8 m tower row",
    value="RMS 2.9 m",
    n=4,
    grade="C",
    source="Witcofski & Chirivella 1980, NASA Langley, Table 4",
    detail={"mean error": "-1.0 m", "regime": "4 of 4"},
)

#: Buoyancy regime across the mixing line.
#:
#: A calculation on the adiabatic mixing line, with no fitted parameter.  The
#: cloud is denser than ambient only above 99.97 mole per cent; the lower
#: flammable limit, stoichiometric and the *upper* flammable limit are all
#: buoyant, which is why a dense-gas closure alone cannot carry an LH2 release.
#:
#: The figure is a property of the top of the mixing table, where air is
#: treated as a non-condensing ideal gas (``CoolPropBackend.cp_air`` clamps at
#: 100 K) and one linear segment spans 0.99894 to 1.0.  Quote the flammable
#: limits, which are far from that segment, rather than the crossover itself.
NEUTRAL_BUOYANCY = Claim(
    quantity="denser than ambient only above",
    value="99.97 mol %",
    n=None,
    grade="A",
    source="adiabatic mixing line, CoolProp equations of state",
    detail={"rho/rho_a at LFL": "0.98", "at stoich": "0.87", "at UFL": "0.69"},
)

#: Reproduction of the original Fortran.
REFERENCE_PARITY = Claim(
    quantity="agreement with DEGADIS 2.1",
    value="1e-12 relative",
    n=None,
    grade="A",
    source="the original Fortran, built from source and probed at full precision",
    detail={"programs": "6 of 6", "EPA test cases": "5 of 5"},
)

#: The five adopted jet corrections, measured against the quantity each was meant to
#: affect rather than against concentration.  Judging a change on concentration
#: alone is how four of the five comparison errors in this project happened.
JET_CORRECTIONS = Claim(
    quantity="the liquid-hydrogen jet corrections, on their target sub-models",
    value="sigma_z ratio 0.901 -> 1.033",
    n=23,
    grade="B",
    source="PRESLHY E3.5, well-constrained vertical fits on momentum-driven "
           "horizontal trials, mean of per-fit ratios; both configurations "
           "validated against parameter dumps to four significant figures",
    detail={
        "remaining mean width bias": "+3.3 %",
        "Gaussian convention": "standard deviation in model and measurement",
        "mean signed centre-height error": "0.036 m",
        "median signed centre-height error": "0.002 m",
        "concentration MG": "1.156 -> 1.047 on 62 common sensor arcs",
    },
)

#: Burro, at the lowest instrumented height -- what the literature reports.
#: **Unsupported.** The only claim in this package with nothing computing it.
#:
#: The complete REDIPHEM reduction now builds runnable cases and the model
#: runs on all eight trials, but the published statistic does not come out:
#: against arc maxima at the lowest instrumented height, comparing the model's
#: ground-level centreline and merging arcs within fifteen per cent, MG is
#: 0.209 over 40 arcs. Running EPA's own B9 deck through the same path gives
#: the same behaviour, so the model is not being run wrongly -- the published
#: figure must pair something other than the centreline against the arc
#: maximum, and what that is has not been established.
BURRO_LOWEST_HEIGHT = Claim(
    quantity="Burro LNG spills, lowest instrumented height",
    value="MG 0.620",
    n=61,
    grade="B",
    source="REDIPHEM complete reduction; DEGADIS's own vertical profile "
           "evaluated at the 1 m sensor height, arcs unmerged, readings at "
           "or below 0.1 vol % dropped",
    detail={"VG": "3.32", "FAC2": "0.54",
            "published": "MG 0.811, VG 2.53, FAC2 0.56, same n"},
)

CLAIMS = {
    "near_field_concentration": NEAR_FIELD_CONCENTRATION,
    "near_field_corrected": NEAR_FIELD_CORRECTED,
    "liftoff_height": LIFTOFF_HEIGHT,
    "neutral_buoyancy": NEUTRAL_BUOYANCY,
    "reference_parity": REFERENCE_PARITY,
    "jet_corrections": JET_CORRECTIONS,
    "burro_lowest_height": BURRO_LOWEST_HEIGHT,
}


# ==========================================================================
# where the model has been checked
# ==========================================================================

#: The conditions each family of releases was checked over.
#:
#: These were previously one dictionary, which conflated two campaigns: the
#: wind range 1.5 to 6.3 m/s is the NASA pool spills, and the distance 35 m is
#: their tower row, but both were being applied to jets whose corrected
#: evidence spans 0.79 to 6 m at winds of 0.6 to 4.2 m/s.  A jet asked about 30 m was
#: therefore reported as inside the checked range when nothing had checked it.
RANGE = {
    "jet": {
        "rate": (0.084, 0.285),  # kg/s, the PRESLHY E3.5 flow meter
        "wind": (0.6, 4.2),  # m/s
        # Seven sensors at 0.35/0.53 m lie upstream of the conserved
        # equivalent-source plane.  The corrected jet is established by the
        # next instrumented arc, 0.79 m; do not imply validation before it.
        "distance": (0.79, 6.0),  # m
        "orifice": (0.006, 0.0254),  # m
        "height": (0.5, 1.5),  # m, release elevation
    },
    "pool": {
        "rate": (9.2, 10.3),  # kg/s, the NASA Langley spills
        "wind": (1.55, 6.30),  # m/s
        "distance": (0.0, 33.8),  # m, the furthest tower row
        "diameter": (9.1, 9.1),  # m, one pond
    },
}

#: Where a "range" is really one point, say so rather than let the warning
#: read as though the caller narrowly missed a boundary.
#:
#: The four NASA spills vary the wind by a factor of four and the release size
#: by twelve per cent: they are one spill in four winds, not a spill-rate
#: series. Reporting 9.0 kg/s as "outside 9.2 to 10.3" without that context
#: makes a warning out of a rounding difference, and warnings that fire on
#: nothing are how a reader learns to skip them.
NOT_A_RANGE = {
    ("pool", "rate"): "the four NASA spills are one release size, not a series",
    ("pool", "diameter"): "one 9.1 m pond; no other pool size has been checked",
    ("jet", "height"): "two elevations were tested, 0.5 and 1.5 m",
}

#: Below this ratio of exit velocity to wind speed the plume is steered by the
#: wind rather than by its own momentum.  On the PRESLHY 1 barg trials,
#: nominally identical releases gave arc maxima of 83 % and 4 %; a steady jet
#: model has nothing to say about that.
MOMENTUM_RATIO = 10.0


def check_range(kind: str, **values) -> list[str]:
    """Warnings for values outside the range ``kind`` was checked over.

    ``values`` of ``None`` are skipped, so a caller can pass everything it has
    without testing each one first.
    """
    limits = RANGE[kind]
    # A value of None means "not supplied", so callers can pass everything
    # they have without testing each one. Only a supplied value under a name
    # that is not a limit is a mistake.
    supplied = {k: v for k, v in values.items() if v is not None}
    unknown = set(supplied) - set(limits)
    if unknown:
        # A misspelled keyword used to be silently ignored, which turns a
        # scope check into no check at all without anything saying so. The
        # caller passes everything it has and lets `None` be skipped, so a
        # name that is not a limit is a mistake rather than a shorthand.
        raise TypeError(
            f"no {kind} limit for {sorted(unknown)}; "
            f"known limits are {sorted(limits)}"
        )
    out = []
    for name, value in supplied.items():
        span = limits[name]
        lo, hi = span
        if not (lo <= value <= hi):
            message = (
                f"{name} {value:g} is outside the {kind} range "
                f"{lo:g} to {hi:g}"
            )
            note = NOT_A_RANGE.get((kind, name))
            if note:
                message += f" -- {note}"
            out.append(message)
    return out


# ==========================================================================
# drift guard
# ==========================================================================

#: Figures that were true once and are not now.  A test greps the prose for
#: these; each maps to what should be written instead.
#:
#: The near-field statistic was superseded when the far-field arcs were
#: recovered by fitting the plume position out, which took the sample from 53
#: arcs on 7 trials to 69 on 9.  The buoyancy figures were superseded when the
#: crossover was computed from the mixing table rather than estimated.
RETIRED = {
    "MG 0.74, VG 1.22": "MG 0.738, VG 1.41 -- 69 arcs, not 53",
    "FAC2 0.92": "FAC2 0.83",
    "53 arc maxima": "69 arc maxima",
    "n=53": "n=69",
    "85-90 mol": "99.97 mol % -- computed, not estimated",
    "85\u201390 mol": "99.97 mol % -- computed, not estimated",
    "85 to 90 mol": "99.97 mol %",
    "RMS of 2.8 m": "RMS 2.9 m -- regenerate from the suite",
    "mean of \u22120.1 m": "mean error -1.0 m",
    "168 pass": "count the suite instead of quoting it",
    "MG 0.855": "MG 1.132 at-sensor; 0.722 on the centreline",
    "0.672, 1.064": "[0.726, 2.303] at-sensor",
    "1.228": "MG 1.113 corrected, same convention, same arcs",
    "1.07 m of rise": "0.877 m at 6 m on trial 10; state the band",
    "1.07 m | **0.19 m**": "withdrawn at both ends; the pair is 0.877 -> 0.658 m",
    "103 tests": "count the suite instead of quoting it",
    "167 tests": "count the suite instead of quoting it",
    "192 pass": "count the suite instead of quoting it",
    "alfa1 = 0.0833": "alfa1 = 0.0875 -- Papanicolaou & List's measured value",
    "Ri_p = 0.557": "Ri_p = 0.716 -- Papanicolaou & List's measured value",
}
