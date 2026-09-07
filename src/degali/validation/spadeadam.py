"""The Spadeadam liquid hydrogen tests, at 30 to 100 m.

FFI-RAPPORT 20/03101 (Aaneby, Gjesdal and Voie, 2021), tests run at DNV
Spadeadam in December 2019. Tests 1--7 are outdoor releases. Two of them --
tests 4 and 6 -- are horizontal from 0.50 m through a 25.4 mm orifice, which
is the geometry the PRESLHY trials use, at three times the flow and five to
seventeen times the distance. Tests 8--15 are a separate closed-room and
ventilation-mast experiment and must not be interpreted as outdoor releases
from the supply nozzle.

Why this module exists
----------------------
Every liquid hydrogen number in this package rests on an array reaching 6 m.
A facility assessment needs tens to hundreds of metres. That is one to two
orders of extrapolation, and nothing tested it.

These tests do. They also carry a constraint the near-field array cannot give:
**three heights on each arc**, so a plume that has lifted out of the array
shows up as a near-zero reading rather than being invisible.

Reading the measurements
------------------------
Two properties of the arrays matter more than they look.

**The arc maxima are censored from below.** Five sensors 20 m apart at 50 m,
two 35 m apart at 100 m, and the plume is documented as missing them -- test 6
reports 2 % at 50 m between arcs that read 21 % and nothing, which its own
authors read as the plume passing between the sensors. A reported value is a
lower bound on what was there. That is enough for the comparison here, which is
whether the model puts *anything* at ground level.

**The vertical structure is the signal.** At 30 m test 6 reads 13.5, 15.4 and
13.3 vol % at 0.1, 1.0 and 1.8 m. A plume that has lifted clear would fall off
steeply with height; one sitting on the ground is flat. Flat is what the data
shows, and it is a statement about the trajectory that needs no absolute
calibration.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DEFAULT_ROOT = Path(__file__).resolve().parents[3] / "reference" / "spadeadam"

#: Outdoor releases in the open-field part of the campaign.
OUTDOOR = tuple(range(1, 8))

#: Closed-room releases exhausted through the ventilation mast. The tabulated
#: nozzle orientation describes the inlet inside the test container, not the
#: atmospheric source at the mast outlet.
CLOSED_ROOM = tuple(range(8, 16))
VENTILATION_MAST = CLOSED_ROOM

#: Outdoor tests released horizontally from 0.50 m -- the PRESLHY geometry.
HORIZONTAL = (4, 6)

#: Outdoor tests released downwards and impinging on the ground.
OUTDOOR_DOWNWARD = (1, 2, 3, 5, 7)

#: Release elevation, m. Not in the sensor tables; from the report text and
#: from Hansen and Hansen (2023), who model the same two tests.
RELEASE_HEIGHT = 0.50

#: Ambient conditions during the December campaign. The tables record no
#: temperature or humidity; 4 C and 90 % are what Hansen and Hansen assume for
#: the site in winter, following the DNV modelling team.
AMBIENT_TEMPERATURE = 277.15
RELATIVE_HUMIDITY = 90.0

#: Ventilation-mast outlet geometry relative to the common TCS-floor origin,
#: from DNV GL 902696's instrument layout. Bearings are clockwise from north.
MAST_DIAMETER = 0.450
MAST_EAST = 5.205
MAST_NORTH = 0.0
MAST_HEIGHT = 11.5

#: DNV reports that field drift below about this level was not removed.
FIELD_DETECTION_LIMIT = 0.5  # vol %


@dataclass
class Reading:
    """One sensor, one test."""

    sensor: str
    radius: float  #: m from the release
    height: float  #: m above ground
    bearing: float  #: degrees
    mean: float  #: vol %
    peak: float  #: vol %
    over_range: bool


@dataclass(frozen=True)
class MastSource:
    """Steady outlet estimate reported for one ventilation-mast test.

    DNV GL inferred these values from the tanker mass loss and three stack
    temperatures, assuming that all released hydrogen left through the
    450 mm mast and that the stack contained 100 vol % hydrogen. They are an
    atmospheric source estimate, unlike the supply-nozzle fields on
    :class:`Trial`.
    """

    start: float  #: averaging-window start, s
    end: float  #: averaging-window end, s
    temperature: float  #: K at the mast outlet
    density: float  #: kg/m3 at the mast outlet
    volume_flow: float  #: m3/s
    velocity: float  #: m/s through the 450 mm mast
    mass_flow: float  #: kg/s
    ambient_pressure: float  #: Pa
    hydrogen_fraction: float  #: assumed volume fraction
    report_pdf_page: int
    note: str = ""


@dataclass
class Trial:
    """One release."""

    test: int
    orifice: float  #: m
    orientation: str
    rate: float  #: kg/s
    wind_high: float  #: m/s
    wind_low: float  #: m/s
    wind_direction: float  #: meteorological direction wind blows from, deg
    line_pressure: float  #: barg at P04
    mast_source: MastSource | None = None
    readings: list[Reading] = field(default_factory=list)

    @property
    def horizontal(self) -> bool:
        return self.test in HORIZONTAL

    @property
    def outdoor(self) -> bool:
        return self.test in OUTDOOR

    @property
    def closed_room(self) -> bool:
        return self.test in CLOSED_ROOM

    @property
    def downward_outdoor(self) -> bool:
        return self.test in OUTDOOR_DOWNWARD

    def arc(self, radius: float, *, tolerance: float = 0.5) -> dict:
        """Peak reading at each height on one arc, maximised over bearing.

        Maximising over bearing is what makes the value comparable with a
        centreline: it is the closest the array came to the plume. It stays a
        lower bound, because the plume can pass between sensors.
        """
        out: dict[float, float] = {}
        for r in self.readings:
            if r.over_range or abs(r.radius - radius) > tolerance:
                continue
            out[r.height] = max(out.get(r.height, 0.0), r.peak)
        return dict(sorted(out.items()))

    def radii(self) -> list[float]:
        return sorted({r.radius for r in self.readings if not r.over_range})


def load(root=None) -> list[Trial]:
    """Read the extracted tables."""
    root = Path(root or DEFAULT_ROOT)
    trials = {}
    with open(root / "conditions.csv") as fh:
        for row in csv.DictReader(fh):
            n = int(row["test"])
            trials[n] = Trial(
                test=n,
                orifice=float(row["orifice_mm"]) / 1000.0,
                orientation=row["orientation"],
                rate=float(row["flow_kgs"]),
                wind_high=float(row["wind_high_ms"]),
                wind_low=float(row["wind_low_ms"]),
                wind_direction=float(row["wind_dir_deg"]),
                line_pressure=float(row["P04_barg"] or 0.0),
            )
    mast_table = root / "mast_conditions.csv"
    if mast_table.exists():
        with open(mast_table, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                n = int(row["test"])
                if n not in trials:
                    continue
                trials[n].mast_source = MastSource(
                    start=float(row["averaging_start_s"]),
                    end=float(row["averaging_end_s"]),
                    temperature=float(row["outlet_temperature_C"]) + 273.15,
                    density=float(row["outlet_density_kgm3"]),
                    volume_flow=float(row["outlet_volumetric_flow_m3s"]),
                    velocity=float(row["outlet_velocity_ms"]),
                    mass_flow=float(row["mass_flow_kgs"]),
                    ambient_pressure=float(row["ambient_pressure_mbar"]) * 100.0,
                    hydrogen_fraction=float(row["h2_vol_fraction_assumed"]),
                    report_pdf_page=int(row["report_pdf_page"]),
                    note=row["notes"],
                )
    with open(root / "sensors.csv") as fh:
        for row in csv.DictReader(fh):
            n = int(row["test"])
            if n not in trials:
                continue
            trials[n].readings.append(Reading(
                sensor=row["sensor"], radius=float(row["radius_m"]),
                height=float(row["height_m"]), bearing=float(row["bearing_deg"]),
                mean=float(row["mean_pct"]), peak=float(row["peak_pct"]),
                over_range=bool(int(row["over_range"])),
            ))
    return [trials[k] for k in sorted(trials)]


# ==========================================================================
# comparison
# ==========================================================================


@dataclass
class ArcComparison:
    """One arc: what the array saw, and what the model puts there."""

    test: int
    radius: float
    #: Measured peak at each height, vol %.
    observed: dict
    #: Modelled concentration at the same heights, vol %.
    predicted: dict
    #: Where the model puts the plume centre, m.
    centre: float
    sigma_z: float

    @property
    def top_sensor(self) -> float:
        return max(self.observed)

    @property
    def lifted(self) -> bool:
        """Whether the model has put the plume clear of the array."""
        return self.centre > self.top_sensor

    @property
    def observed_max(self) -> float:
        return max(self.observed.values())

    @property
    def predicted_max(self) -> float:
        return max(self.predicted.values())

    @property
    def observed_flatness(self) -> float:
        """Ratio of the smallest to the largest reading across the arc.

        Near one means the profile is flat over 0.1 to 1.8 m, which is what a
        grounded, well-mixed plume looks like. Near zero means it falls off
        with height. This is the trajectory signal that needs no absolute
        calibration.
        """
        values = [v for v in self.observed.values() if v > 0.0]
        if len(values) < 2:
            return float("nan")
        return min(values) / max(values)

    @property
    def predicted_flatness(self) -> float:
        values = [v for v in self.predicted.values() if v > 0.0]
        if len(values) < 2:
            return float("nan")
        return min(values) / max(values)


@dataclass
class MastComparison:
    """One external field sensor relative to the ventilation-mast plume."""

    test: int
    sensor: str
    x: float  #: downwind of mast outlet, m
    y: float  #: crosswind of mast outlet, m
    z: float  #: height above ground, m
    observed: float  #: measured peak, vol %
    predicted: float  #: model at the sensor, vol %
    centre: float  #: modelled plume-centre height at x, m

    @property
    def censored(self) -> bool:
        return self.observed < FIELD_DETECTION_LIMIT

    @property
    def false_flammable(self) -> bool:
        return self.censored and self.predicted >= 4.0


def mast_coordinates(reading: Reading, wind_from: float) -> tuple[float, float]:
    """Return a sensor's downwind/crosswind coordinates from the mast outlet.

    ``wind_from`` follows the meteorological convention. Sensor bearings are
    clockwise from north, as in the DNV tables.
    """
    bearing = math.radians(reading.bearing)
    east = reading.radius * math.sin(bearing) - MAST_EAST
    north = reading.radius * math.cos(bearing) - MAST_NORTH
    direction = math.radians(wind_from)
    down_east, down_north = -math.sin(direction), -math.cos(direction)
    right_east, right_north = down_north, -down_east
    return (
        east * down_east + north * down_north,
        east * right_east + north * right_north,
    )


def _external_mast_sensor(reading: Reading) -> bool:
    """Field sensors belonging to the mast comparison, fixed in the prereg."""
    try:
        number = int(reading.sensor.rsplit("_", 1)[1])
    except (IndexError, ValueError):
        return False
    return 1 <= number <= 24 or 28 <= number <= 30


def compare_mast(trial: Trial, *, wind=None, **kw) -> list[MastComparison]:
    """Compare one closed-room test against its vertical mast-outlet source.

    Tests 8--14 only. This is a negative-constraint comparison: most field
    readings sit below the documented 0.5 vol % drift level, so callers should
    count :attr:`MastComparison.false_flammable` rather than compute logarithmic
    performance statistics.
    """
    if not trial.closed_room:
        raise ValueError(f"test {trial.test} is not a ventilation-mast test")
    if trial.mast_source is None:
        raise ValueError(
            f"test {trial.test} has no intact ventilation-mast source"
        )

    from .nearfield import STEP, Trajectory, hydrogen_gas_jet

    source = trial.mast_source
    jp, y0 = hydrogen_gas_jet(
        rate=source.mass_flow, diameter=MAST_DIAMETER,
        velocity=source.velocity,
        wind=trial.wind_high if wind is None else wind,
        height=MAST_HEIGHT, source_temperature=source.temperature,
        source_density=source.density,
        ambient_temperature=AMBIENT_TEMPERATURE,
        ambient_pressure=source.ambient_pressure,
        relative_humidity=RELATIVE_HUMIDITY,
        wind_reference_height=10.0,
        **kw,
    )
    rows = jp.run(y0, distmx=STEP, smax=250.0).rows
    trajectory = Trajectory(jp.th.table, rows)

    out = []
    for reading in trial.readings:
        if reading.over_range or not _external_mast_sensor(reading):
            continue
        x, y = mast_coordinates(reading, trial.wind_direction)
        if x <= 0.0:
            continue
        state = trajectory.at(x)
        if state is None:
            continue
        out.append(MastComparison(
            test=trial.test, sensor=reading.sensor,
            x=x, y=y, z=reading.height, observed=reading.peak,
            predicted=trajectory.concentration_at(x, y, reading.height),
            centre=state.z,
        ))
    return out


def compare(trial: Trial, *, corrections: bool, wind=None, **kw) -> list:
    """Run the model over one trial and pair it with the arcs.

    ``wind`` defaults to the low mast, which is the one Hansen and Hansen use
    as representative. The tables give a high and a low mast without stating
    their reference heights, so the choice is recorded rather than derived.
    """
    if trial.test not in HORIZONTAL:
        if trial.closed_room:
            raise ValueError(
                f"test {trial.test} is a closed-room/ventilation-mast test; "
                "the mast-outlet source term is not available"
            )
        raise ValueError(
            f"test {trial.test} is an outdoor downward release; "
            "use `compare_downward`"
        )

    from .nearfield import REACH, STEP, Trajectory, hydrogen_jet

    jp, y0 = hydrogen_jet(
        rate=trial.rate, diameter=trial.orifice,
        wind=wind if wind is not None else trial.wind_low,
        height=RELEASE_HEIGHT,
        ambient_temperature=AMBIENT_TEMPERATURE,
        relative_humidity=RELATIVE_HUMIDITY,
        storage_pressure_barg=trial.line_pressure,
        corrections=corrections, **kw,
    )
    rows = jp.run(y0, distmx=STEP, smax=max(REACH, 250.0)).rows
    traj = Trajectory(jp.th.table, rows)

    out = []
    for radius in trial.radii():
        observed = trial.arc(radius)
        if not observed:
            continue
        state = traj.at(radius)
        if state is None:
            continue
        out.append(ArcComparison(
            test=trial.test, radius=radius, observed=observed,
            predicted={
                z: traj.concentration_at(radius, 0.0, z) for z in observed
            },
            centre=state.z, sigma_z=state.sz,
        ))
    return out


def report(arcs: list) -> str:
    lines = [
        f"{'test':>5}{'x (m)':>7}{'centre':>8}{'sig_z':>7}"
        f"{'obs max':>9}{'pred max':>10}{'obs/pred':>10}"
        f"{'obs flat':>10}{'pred flat':>11}"
    ]
    for a in arcs:
        ratio = a.observed_max / a.predicted_max if a.predicted_max > 0 else float("inf")
        lines.append(
            f"{a.test:>5}{a.radius:>7.1f}{a.centre:>8.2f}{a.sigma_z:>7.2f}"
            f"{a.observed_max:>9.2f}{a.predicted_max:>10.3f}{ratio:>10.1f}"
            f"{a.observed_flatness:>10.2f}{a.predicted_flatness:>11.2f}"
        )
    return "\n".join(lines)


# ==========================================================================
# the downward releases
# ==========================================================================

#: Source footprint handed to the ground-level plume, m.
#:
#: **It does not matter, and that is the point.** A downward jet at these
#: conditions reaches the ground inside its own development length --
#: `initial_conditions_directed` refuses it, correctly, with "the jet reaches
#: the ground within its development length (0.42 m from a 0.32 m release)" --
#: so there is no jet to integrate and the fluid arrives essentially
#: undiluted. What is left unknown is the radius of the impingement footprint,
#: and sweeping it over eighty-fold, 0.05 m to 4 m, moves the concentration at
#: 30 m by seven per cent and the plume centre by half a metre.
#:
#: So the unknown was measured before a value was chosen for it, rather than
#: after.
FOOTPRINT = 0.5

#: Release height for the outdoor downward tests, m. From the narrative of
#: Hansen and Hansen for test 5; the tables do not carry it.
DOWNWARD_HEIGHT = 0.32


def compare_downward(trial: Trial, radius: float = 30.0, *,
                     footprint: float = FOOTPRINT, wind=None) -> tuple:
    """Model a downward release through the ground-level buoyant plume.

    Returns ``(observed_max, predicted_centreline, centre_height, regime)``
    in vol %, or ``None`` if the flammable cloud does not reach ``radius``.

    This goes through `addons.LiftoffPlume`, the URAHFREP model, and not
    through `JetPlume`. That makes the two halves of the Spadeadam comparison
    **independent**: the horizontal tests exercise the jet path and these
    exercise the ground-level buoyant path, so a fault that shows in both is
    not an artefact of either.
    """
    if trial.test not in OUTDOOR_DOWNWARD:
        if trial.horizontal:
            raise ValueError(f"test {trial.test} is horizontal; use `compare`")
        if trial.closed_room:
            raise ValueError(
                f"test {trial.test} is a closed-room/ventilation-mast test; "
                "the mast-outlet source term is not available"
            )
        raise ValueError(f"test {trial.test} is not an outdoor downward release")

    import numpy as np

    from ..lh2 import assess

    arc = trial.arc(radius)
    if not arc:
        return None
    result = assess(
        rate=trial.rate, pool_diameter=footprint,
        wind=wind if wind is not None else trial.wind_low,
        ambient_temperature=AMBIENT_TEMPERATURE,
        relative_humidity=RELATIVE_HUMIDITY,
        max_distance=5.0 * radius, at_distance=radius,
    )
    traj = result.trajectory
    if len(traj) == 0 or traj[-1, 0] < radius:
        return None
    return (
        max(arc.values()),
        float(np.interp(radius, traj[:, 0], traj[:, 2])) * 100.0,
        float(np.interp(radius, traj[:, 0], traj[:, 1])),
        result.regime,
    )
