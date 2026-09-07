"""The SMEDIS trial spreadsheets.

The EU Scientific Model Evaluation of Dense Gas Dispersion (SMEDIS) exercise
reduced a large set of field and wind-tunnel trials into one spreadsheet
each, in a fixed layout.  They are worth reading in preference to the raw
archives for three reasons.

**The channels are already identified.**  REDIPHEM stores concentrations as
numbered channels whose meaning is per series and defined in a file some
series do not ship.  FLADIS is one of those, and guessing its numbering
selected channels reading 302 and 23.5 at 20 m downwind — not concentrations
at all — and produced a clean-looking result that meant nothing.  Here every
sensor is a row with its position and its value.

**The wind direction is there, with its standard deviation.**  Every
comparison in this package against a fixed sensor array has had to work around
not knowing how much the plume wandered; these files record it.

**The arc reductions are there too.**  Distance, height, arc-maximum
concentration and the crosswind spread, already computed, so a lateral
comparison does not depend on reconstructing them.

Layout
------
Each sheet is a sequence of ``# Section`` headers followed by rows.  The
sections used here:

============================  =================================================
``# Dataset_Reference``       the series name, in the second column
``# Trial_Identification``    the trial name
``# Release_Conditions``      substance, rate, geometry, one per row
``# Ambient_Conditions``      wind, roughness, Monin-Obukhov, stability
``# Wind_Direction_Sensors``  ``x, y, z, mean dir, std dir``
``# Concentration_Sensors``   ``x, y, z, mean C (%), std C (%)``
``# Arc_Positions``           ``distance, height, C(arcmax, %), sy``
============================  =================================================

Column headers are themselves prefixed with ``#``, so a ``#`` line is a
section header only when its second column is empty.

Values of ``-999``, ``-99`` and blanks mean "not recorded" and become ``nan``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: The sentinels SMEDIS uses for a missing value.
MISSING = (-999.0, -99.0, -9.99)


def _number(value) -> float:
    if value is None or value == "":
        return float("nan")
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return float("nan") if any(abs(out - m) < 1e-9 for m in MISSING) else out


@dataclass
class Sensor:
    """One concentration sensor."""

    x: float
    y: float
    z: float
    mean: float  #: per cent
    std: float  #: per cent


@dataclass
class Arc:
    """One sampling arc, as SMEDIS reduced it."""

    distance: float
    height: float
    concentration: float  #: arc maximum, per cent
    sigma_y: float  #: crosswind spread, m


@dataclass
class SmedisTrial:
    """One trial."""

    path: Path
    dataset: str = ""
    trial: str = ""
    substance: str = ""
    conditions: dict = field(default_factory=dict)
    sensors: list[Sensor] = field(default_factory=list)
    arcs: list[Arc] = field(default_factory=list)
    #: Where the release sits in the file's coordinate system, m.  Most
    #: trials put it at the origin, but Thorney Island uses site grid
    #: coordinates with the release at (400, 200), so sensor positions there
    #: are not downwind distances until this is subtracted.
    origin: tuple = (0.0, 0.0)
    #: Circular standard deviation of the recorded wind directions, degrees.
    #: This is the meander that makes a fixed sensor array hard to compare
    #: against a steady plume, and it is usually the largest uncertainty in
    #: such a comparison.
    wind_direction_spread: float = float("nan")

    @property
    def wind(self) -> float:
        return self.conditions.get("wind speed", float("nan"))

    @property
    def roughness(self) -> float:
        return self.conditions.get("surface roughness", float("nan"))

    @property
    def suspect_units(self) -> bool:
        """Whether the concentrations can be read as volume per cent.

        Thorney Island peaks at 2060 and Prairie Grass at 235 in these files;
        both are recorded under a ``mean_C(%)`` header. Something other than
        a volume percentage is being reported and the comparison would be
        meaningless, so it is flagged rather than silently scaled.
        """
        return any(s.mean > 100.0 for s in self.sensors)

    def heights(self) -> list[float]:
        return sorted({s.z for s in self.sensors})

    def arc_distances(self) -> list[float]:
        return sorted({a.distance for a in self.arcs})


_HEADER = re.compile(r"^#\s*(.+?)\s*$")


def read_trial(path: str | Path) -> SmedisTrial:
    """Read one SMEDIS spreadsheet."""
    import xlrd

    path = Path(path)
    sheet = xlrd.open_workbook(str(path)).sheet_by_index(0)
    rows = [
        [sheet.cell_value(r, c) for c in range(sheet.ncols)]
        for r in range(sheet.nrows)
    ]

    out = SmedisTrial(path=path)
    section = ""
    directions: list[float] = []
    for row in rows:
        first = str(row[0]).strip() if row else ""
        head = _HEADER.match(first)
        if head:
            rest = "".join(str(v).strip() for v in row[1:])
            name = head.group(1).lower()
            if name.startswith("dataset_reference"):
                out.dataset = str(row[1]).strip()
            elif name.startswith("trial_identification"):
                out.trial = str(row[1]).strip()
            if not rest:
                section = name
            # otherwise this is a column-header line, also prefixed with '#'
            continue
        if not first:
            continue

        label = first.lower().strip()
        if "release_condition" in section or "release condition" in section:
            if "substance" in label or "chemical released" in label:
                out.substance = str(row[1]).strip()
            else:
                value = _number(row[1] if len(row) > 1 else None)
                if not math.isnan(value):
                    out.conditions[label] = value
        elif "ambient" in section or "condition" in section:
            value = _number(row[1] if len(row) > 1 else None)
            if not math.isnan(value):
                out.conditions[label] = value
        elif "wind_direction" in section or "wind direction" in section:
            # x, y, z, mean direction, standard deviation
            mean = _number(row[3] if len(row) > 3 else None)
            if not math.isnan(mean):
                directions.append(mean)
        elif "concentration_sensor" in section:
            values = [_number(v) for v in row[:5]]
            if len(values) >= 4 and not any(math.isnan(v) for v in values[:4]):
                out.sensors.append(
                    Sensor(*values[:4], std=values[4] if len(values) > 4 else float("nan"))
                )
        elif "source" in section or "release_point" in section:
            pass
        elif "arc_position" in section:
            values = [_number(v) for v in row[:4]]
            if not any(math.isnan(v) for v in values[:3]):
                out.arcs.append(
                    Arc(*values[:3], sigma_y=values[3] if len(values) > 3 else float("nan"))
                )

    # shift sensor positions so that x is downwind distance from the release
    ox = out.conditions.get("release point x", 0.0)
    oy = out.conditions.get("y", 0.0)
    if not math.isnan(ox) and (ox or oy):
        out.origin = (ox, oy)
        out.sensors = [
            Sensor(s.x - ox, s.y - oy, s.z, s.mean, s.std) for s in out.sensors
        ]

    if directions:
        radians = np.radians(directions)
        resultant = math.hypot(
            float(np.cos(radians).mean()), float(np.sin(radians).mean())
        )
        out.wind_direction_spread = math.degrees(
            math.sqrt(-2.0 * math.log(max(resultant, 1e-12)))
        )
    return out


def load(directory: str | Path) -> list[SmedisTrial]:
    """Read every spreadsheet under a directory, skipping unreadable ones."""
    out = []
    for path in sorted(Path(directory).rglob("*.xls")):
        try:
            trial = read_trial(path)
        except Exception:
            continue
        if trial.sensors or trial.arcs:
            out.append(trial)
    return out
