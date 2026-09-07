"""The PRESLHY E3.5 liquid hydrogen trials, from the raw dataset.

Twenty-four LH₂ releases at the HSL facility in September 2019, published as
DOI 10.35097/1481.  Each trial is one workbook with five sheets; two of them
matter here.

``Xensor``
    Thermal-conductivity sensors supplied by NREL, ranging **0 to 100 vol %**.
    They are the reason to use the raw data rather than the report: the
    Dräger devices in the far field top out at 4 %, which is hydrogen's lower
    flammable limit, so every reading that matters is censored.  These are
    not.  Sampled at about 3 Hz.
``Flowmeter``
    Coriolis mass flow at 1 Hz, so the source is measured through the release
    rather than taken from a summary table.

Sensor positions come from Table A3 of the report (deliverable D3.6), keyed by
the same serial numbers the workbook columns carry.

Averaging
---------
Both a mean and a peak are computed over the same window, and both are
returned.  That is deliberate.  Comparing a model centreline against a
*time-averaged* point measurement understates it, because the plume meanders
and is not always over the sensor; comparing against the peak overstates the
averaging the model actually does.  The right answer lies between, and
choosing after seeing which one agrees is how a result gets manufactured.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Column header pattern: ``X2019_09_10_02EC42Output`` carries serial 02EC42.
_SERIAL = re.compile(r"_(\d{2}[A-Z]{2}\d+)Output")

#: Table A3 row: serial, mount/sensor id, then x, y, z in metres.  Two
#: sensors have seven-character serials (``10DC100`` and ``10DC101``); using
#: a fixed six-character field silently dropped both ground-level readings.
_TABLE_A3 = re.compile(
    r"(\d{2}[A-Z]{2}\d+)\s+(M\d+S\d+)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
)

#: Converted near-field thermocouple headings in the Flexlogger sheet.
_NEAR_FIELD_TC = re.compile(
    r"^(Centreline|Above1|Above2|Below1|Below2|Left|Right)_"
    r"(\d{3})_TC(\d+)C$"
)

_TC_OFFSETS = {
    "Centreline": (0.0, 0.50),
    "Above1": (0.0, 0.75),
    "Above2": (0.0, 1.00),
    "Below1": (0.0, 0.25),
    "Below2": (0.0, 0.00),
    "Left": (-1.0, 0.50),
    "Right": (1.0, 0.50),
}

# NIST ITS-90 Type-T forward polynomial, -270 to 0 degC.  The coefficients
# return millivolts and are published in Monograph 175.  PRESLHY D3.6 uses
# the -253 degC value to correct a batch gain error in TC1/TC2/TC3.
_TYPE_T_NEGATIVE_COEFFICIENTS = (
    0.0,
    0.387481063640e-1,
    0.441944343470e-4,
    0.118443231050e-6,
    0.200329735540e-7,
    0.901380195590e-9,
    0.226511565930e-10,
    0.360711542050e-12,
    0.384939398830e-14,
    0.282135219250e-16,
    0.142515947790e-18,
    0.487686622860e-21,
    0.107955392700e-23,
    0.139450270620e-26,
    0.797951539270e-30,
)


def type_t_emf_mv(temperature_c: float) -> float:
    """NIST ITS-90 Type-T thermocouple EMF below 0 degC, in mV."""
    temperature_c = float(temperature_c)
    if not -270.0 <= temperature_c <= 0.0:
        raise ValueError("negative-range Type-T polynomial requires -270..0 degC")
    emf = 0.0
    for coefficient in reversed(_TYPE_T_NEGATIVE_COEFFICIENTS):
        emf = emf * temperature_c + coefficient
    return emf


def type_t_temperature_c(emf_mv: float) -> float:
    """Invert the NIST negative-range Type-T relation."""
    from scipy.optimize import brentq

    emf_mv = float(emf_mv)
    lower = type_t_emf_mv(-270.0)
    if not lower <= emf_mv <= 0.0:
        raise ValueError(
            f"negative-range Type-T EMF must lie between {lower:.6f} and 0 mV"
        )
    return float(brentq(
        lambda temperature: type_t_emf_mv(temperature) - emf_mv,
        -270.0, 0.0, xtol=1.0e-12,
    ))


def compensate_pipe_temperature_c(
    voltage_uv: float,
    tc1_min_voltage_uv: float,
    *,
    reference_temperature_c: float = -253.0,
) -> float:
    """Apply PRESLHY D3.6 equation 1 and return temperature in degC.

    ``voltage_uv`` is the cold-junction-corrected TC1/TC2/TC3 signal and
    ``tc1_min_voltage_uv`` is the minimum TC1 value from the same trial.  The
    report fixes that flat inlet sensor at liquid-hydrogen boiling
    temperature, -253 degC, and applies the resulting batch gain correction
    to all three pipe thermocouples.
    """
    voltage_uv = float(voltage_uv)
    tc1_min_voltage_uv = float(tc1_min_voltage_uv)
    if voltage_uv > 0.0 or tc1_min_voltage_uv >= 0.0:
        raise ValueError("PRESLHY pipe thermocouple voltages must be negative")
    reference_uv = 1000.0 * type_t_emf_mv(reference_temperature_c)
    error_fraction = reference_uv / tc1_min_voltage_uv
    corrected_uv = voltage_uv * (
        1.0
        + (error_fraction - 1.0)
        * abs(voltage_uv) / abs(tc1_min_voltage_uv)
    )
    return type_t_temperature_c(corrected_uv / 1000.0)


@dataclass(frozen=True)
class PipeSourceObservation:
    """One steady-window reconstruction of the PRESLHY pipe source."""

    trial: int
    path: Path
    window: tuple[int, int]
    samples: int
    tc1_min_voltage_uv: float
    nozzle_temperature_min_k: float
    nozzle_temperature_median_k: float
    nozzle_temperature_mean_k: float
    nozzle_temperature_p05_k: float
    nozzle_temperature_p95_k: float
    pt1_median_barg: float
    pt1_mean_barg: float
    pt2_median_barg: float
    pt2_mean_barg: float
    flex_mass_flow_median_g_s: float
    flex_mass_flow_mean_g_s: float
    coriolis_mass_flow_median_g_s: float
    coriolis_mass_flow_mean_g_s: float
    drive_gain_median_pct: float
    calculated_density_median_kg_m3: float
    calculated_density_mean_kg_m3: float


def _finite_column(rows, column: int) -> np.ndarray:
    values = []
    for row in rows:
        try:
            value = float(row[column])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(value):
            values.append(value)
    if not values:
        raise ValueError(f"source column {column + 1} has no finite samples")
    return np.asarray(values, dtype=float)


def read_pipe_source(
    workbook: str | Path,
    window: tuple[int, int],
) -> PipeSourceObservation:
    """Read PT1/PT2/TC3 and Coriolis diagnostics over a frozen window.

    ``window`` uses the same zero-based, end-exclusive convention returned by
    :func:`read_trial`.  Both source sheets are sampled at 1 Hz, so sample
    index zero is Excel row 2.  Pipe TC3 is reconstructed from the recorded
    microvolt signal using the trial's TC1 minimum and D3.6 equation 1; the
    workbook's converted TC columns are intentionally not used.
    """
    import openpyxl

    path = Path(workbook)
    match = re.search(r"trial_(\d+)_", path.name)
    if match is None:
        raise ValueError(f"cannot identify PRESLHY trial from {path.name}")
    first, last = map(int, window)
    if first < 0 or last <= first:
        raise ValueError("source window must be non-empty and non-negative")

    book = openpyxl.load_workbook(path, read_only=True, data_only=True)

    def sheet_rows(name: str):
        sheet = book[name]
        rows = sheet.iter_rows(
            min_row=first + 2,
            max_row=last + 1,
            values_only=True,
        )
        return list(rows)

    def headings(name: str) -> dict[str, int]:
        sheet = book[name]
        row = next(sheet.iter_rows(min_row=1, max_row=1, values_only=True))
        return {str(value): index for index, value in enumerate(row) if value}

    flex_header = headings("Flexlogger")
    flow_header = headings("Flowmeter")
    flex_rows = sheet_rows("Flexlogger")
    flow_rows = sheet_rows("Flowmeter")
    required_flex = {
        "flow": "MFM1_Mass_Flow_Rate",
        "pt1": "PT1_Pipe_Pressure",
        "pt2": "PT2_Nozzle_Pressure",
        "tc1": "TC1_MFM_In",
        "tc3": "TC3_Release_Nozzle_Flow",
    }
    required_flow = {
        "mass": "FlMassFlowRategsecR0247",
        "drive": "FlDriveGainR0291",
        "density": "FlCalcDensitykgm3",
    }
    missing = [
        heading
        for heading in (*required_flex.values(), *required_flow.values())
        if heading not in (flex_header if heading in required_flex.values() else flow_header)
    ]
    if missing:
        raise KeyError(f"missing PRESLHY source headings: {missing}")

    flex = {
        key: _finite_column(flex_rows, flex_header[heading])
        for key, heading in required_flex.items()
    }
    flow = {
        key: _finite_column(flow_rows, flow_header[heading])
        for key, heading in required_flow.items()
    }
    tc1_min = float(np.min(flex["tc1"]))
    nozzle_k = np.asarray([
        compensate_pipe_temperature_c(value, tc1_min) + 273.15
        for value in flex["tc3"]
    ])
    book.close()

    return PipeSourceObservation(
        trial=int(match.group(1)),
        path=path,
        window=(first, last),
        samples=len(flex_rows),
        tc1_min_voltage_uv=tc1_min,
        nozzle_temperature_min_k=float(np.min(nozzle_k)),
        nozzle_temperature_median_k=float(np.median(nozzle_k)),
        nozzle_temperature_mean_k=float(np.mean(nozzle_k)),
        nozzle_temperature_p05_k=float(np.percentile(nozzle_k, 5.0)),
        nozzle_temperature_p95_k=float(np.percentile(nozzle_k, 95.0)),
        pt1_median_barg=float(np.median(flex["pt1"])),
        pt1_mean_barg=float(np.mean(flex["pt1"])),
        pt2_median_barg=float(np.median(flex["pt2"])),
        pt2_mean_barg=float(np.mean(flex["pt2"])),
        flex_mass_flow_median_g_s=float(np.median(flex["flow"])),
        flex_mass_flow_mean_g_s=float(np.mean(flex["flow"])),
        coriolis_mass_flow_median_g_s=float(np.median(flow["mass"])),
        coriolis_mass_flow_mean_g_s=float(np.mean(flow["mass"])),
        drive_gain_median_pct=float(np.median(flow["drive"])),
        calculated_density_median_kg_m3=float(np.median(flow["density"])),
        calculated_density_mean_kg_m3=float(np.mean(flow["density"])),
    )


def _table_a3_positions(text: str) -> dict[str, tuple[float, float, float]]:
    """Parse Table A3 after normalising PDF minus-sign variants.

    The published PDF extracts its negative lateral coordinates with a
    Unicode hyphen, sometimes separated from the number by a space.  Python's
    ordinary ``-?`` pattern does not match it, so the whole negative-y side of
    the array used to disappear from the reduction without an error.
    """
    for dash in ("\N{HYPHEN}", "\N{NON-BREAKING HYPHEN}",
                 "\N{EN DASH}", "\N{MINUS SIGN}"):
        text = text.replace(dash, "-")
    text = re.sub(r"-\s+(?=\d)", "-", text)
    return {
        serial: (float(x), float(y), float(z))
        for serial, _mount, x, y, z in _TABLE_A3.findall(text)
    }


def _gaussian(coordinate, peak, centre, sigma):
    """Gaussian profile whose width is the statistical standard deviation.

    The former fits used ``exp(-((x-centre)/width)**2)`` but stored ``width``
    under the name ``sigma``.  That parameter is the e-folding width
    ``sqrt(2) * sigma`` and is not directly comparable with JETPLU, whose
    concentration profile is ``exp(-0.5*(x/sigma)**2)``.  Keeping the
    standard form here makes the fitted and modelled widths use one contract.
    """
    return peak * np.exp(-0.5 * (((coordinate - centre) / sigma) ** 2))


def sensor_positions(report: str | Path, pages=range(42, 46)) -> dict:
    """Read Table A3 from the D3.6 report: serial -> ``(x, y, z)`` in metres."""
    import pypdf

    reader = pypdf.PdfReader(str(report))
    text = "".join((reader.pages[i].extract_text() or "") for i in pages)
    return _table_a3_positions(text)


@dataclass
class Reading:
    """One sensor in one trial, reduced over the release window."""

    serial: str
    x: float
    y: float
    z: float  #: height above the ground, m
    z_axis: float  #: height above the release axis, m, as Table A3 gives it
    mean: float  #: vol % over the window
    peak: float  #: vol % , the maximum in the window
    samples: int


@dataclass
class TemperatureReading:
    """One converted near-field thermocouple over a frozen release window."""

    channel: str
    serial: str
    x: float
    y: float
    z: float  #: height above ground, m
    z_axis: float  #: array coordinate before the release-height shift, m
    minimum_c: float
    percentile_05_c: float
    median_c: float
    samples: int


def read_nearfield_temperatures(
    path: str | Path,
    *,
    release_height: float = 0.5,
    threshold: float = 0.5,
) -> tuple[tuple[int, int], list[TemperatureReading]]:
    """Reduce the collocated type-T array over the mass-flow window.

    The returned half-open row interval is the exact interval used by
    :func:`read_trial`.  Temperatures come only from the already converted
    ``...TC<serial>C`` channels; raw thermovoltages, pipework sensors and
    far-field stand thermocouples are deliberately excluded.
    """
    import openpyxl

    path = Path(path)
    window = read_trial(
        path, {}, threshold=threshold, release_height=release_height
    ).window
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = list(book["Flexlogger"].iter_rows(values_only=True))
    header = rows[0]
    lo, hi = window
    selected = rows[1:][lo:hi]
    readings = []
    for column, heading in enumerate(header):
        match = _NEAR_FIELD_TC.match(str(heading))
        if match is None:
            continue
        placement, distance, serial_number = match.groups()
        values = _column(selected, column)
        finite = values[np.isfinite(values)]
        if not finite.size:
            continue
        y, z_axis = _TC_OFFSETS[placement]
        readings.append(TemperatureReading(
            channel=str(heading),
            serial=f"TC{serial_number}",
            x=float(distance) / 100.0,
            y=y,
            z=z_axis + (float(release_height) - 0.5),
            z_axis=z_axis,
            minimum_c=float(np.min(finite)),
            percentile_05_c=float(np.percentile(finite, 5.0)),
            median_c=float(np.median(finite)),
            samples=int(finite.size),
        ))
    return window, readings


@dataclass
class Timing:
    """When a sensor saw the cloud, relative to the start of the release."""

    serial: str
    x: float
    y: float
    z: float
    arrival: float  #: s after the flow reached half its peak
    peak_time: float  #: s, when the maximum was recorded
    departure: float  #: s, last sample above the threshold
    peak: float  #: vol %

    @property
    def duration(self) -> float:
        return self.departure - self.arrival


@dataclass
class Trial:
    """One release, with its measured source."""

    number: int
    path: Path
    flow_peak: float  #: g/s
    flow_mean: float  #: g/s, over the window used
    window: tuple[int, int]  #: first and last flow sample in the window
    readings: list[Reading] = field(default_factory=list)
    #: Arrival, peak and departure times, one per sensor that saw the cloud.
    #: The transient path of the model -- the observers, the time sort, the
    #: receptor histories -- has never been compared against measurement in
    #: this package, and these are what such a comparison needs.
    timings: list[Timing] = field(default_factory=list)
    #: Seconds per sample of the concentration record.
    sample_interval: float = 0.3
    #: Sample index in the *gas* record at which the flow reached half peak.
    release_start: int = 0

    def at_height(self, z: float, tolerance: float = 0.01) -> list[Reading]:
        return [r for r in self.readings if abs(r.z - z) <= tolerance]

    def centreline(self) -> list[Reading]:
        return [r for r in self.readings if abs(r.y) < 0.01]


def _column(rows, index) -> np.ndarray:
    out = np.full(len(rows), np.nan)
    for i, row in enumerate(rows):
        value = row[index] if index < len(row) else None
        if value not in (None, ""):
            try:
                out[i] = float(value)
            except (TypeError, ValueError):
                pass
    return out


def read_trial(
    path: str | Path,
    positions: dict,
    *,
    threshold: float = 0.5,
    floor: float = 0.05,
    release_height: float = 0.5,
) -> Trial:
    """Reduce one workbook.

    Parameters
    ----------
    threshold
        Fraction of the peak mass flow that defines the release window.  The
        flow meter ramps and tails off, so a fixed window would average the
        cloud with the empty record either side -- the same mistake that moved
        a REDIPHEM statistic by an order of magnitude.
    floor
        Readings whose peak stays below this are treated as not having seen
        the cloud and are dropped.
    release_height
        Height of the nozzle above the ground, m.

        Table A3's ``z`` is measured from the *release axis*, not from the
        ground: the report's Figure 4 is captioned "near-field array
        configured for releases at 1.5 m height", and the array was rebuilt
        for each height. The data says the same thing -- at 0.35 m downwind
        the 0.5 m sensor reads 80.8 vol % for a 0.5 m release and 87.9 for a
        1.5 m one, and a jet cannot fall a metre in 35 centimetres of travel.

        Reading it as a height above ground puts the 1.5 m releases' sensors
        a metre below their own axis, where the Gaussian factor is
        ``exp(-(1.0/0.15)**2)`` and the prediction is ten to the minus
        nineteen. That is what an earlier version of this work reported as the
        model failing catastrophically on those ten trials.
    """
    import openpyxl

    path = Path(path)
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)

    flow_rows = list(book["Flowmeter"].iter_rows(values_only=True))
    header = flow_rows[0]
    index = next(
        i for i, h in enumerate(header) if h and "MassFlow" in str(h)
    )
    flow = _column(flow_rows[1:], index)
    if not np.any(np.isfinite(flow)) or np.nanmax(flow) <= 0.0:
        raise ValueError(f"{path.name}: no mass flow recorded")
    # A single spike is enough to stretch a first-to-last window across the
    # whole record: on trial 2 the meter oscillates between -18.8 and +30.5
    # g/s and only six per cent of such a window is actually flowing, which
    # gives a mean mass flow of -0.5 g/s. Smooth first, then take the longest
    # *run* above the threshold rather than the outermost crossings.
    smooth = np.convolve(
        np.nan_to_num(flow), np.ones(5) / 5.0, mode="same"
    )
    peak = float(np.nanmax(smooth))
    if peak <= 0.0:
        raise ValueError(f"{path.name}: mass flow never rises")
    live = smooth > threshold * peak
    best, first, last, run = 0, 0, 0, 0
    for i, on in enumerate(live):
        run = run + 1 if on else 0
        if run > best:
            best, first, last = run, i - run + 1, i + 1
    if best < 5:
        raise ValueError(f"{path.name}: no sustained flow (longest run {best})")

    gas_rows = list(book["Xensor"].iter_rows(values_only=True))
    gas_header = gas_rows[0]
    columns = {
        _SERIAL.search(str(h)).group(1): i
        for i, h in enumerate(gas_header)
        if h and _SERIAL.search(str(h))
    }
    # the two sheets run at different rates; map the flow window onto the gas
    # record by fraction of the record rather than by clock, because the
    # timestamps are strings of differing precision
    n = len(gas_rows) - 1
    lo = int(round(first / len(flow) * n))
    hi = max(int(round(last / len(flow) * n)), lo + 1)

    readings = []
    for serial, column in columns.items():
        if serial not in positions:
            continue
        series = _column(gas_rows[1 + lo : 1 + hi], column)
        if not np.any(np.isfinite(series)):
            continue
        peak = float(np.nanmax(series))
        if peak < floor or peak > 100.0:
            # over-range: the sensors read to 100 vol % and a handful exceed it
            continue
        x, y, z_axis = positions[serial]
        # Table A3's z is above the release axis; the array moves with it
        z = z_axis + (release_height - 0.5)
        readings.append(
            Reading(
                serial=serial, x=x, y=y, z=z, z_axis=z_axis,
                mean=float(np.nanmean(series)), peak=peak,
                samples=int(np.sum(np.isfinite(series))),
            )
        )

    # arrival, peak and departure, relative to the start of the release
    interval = 0.3
    timings = []
    for serial, column in columns.items():
        if serial not in positions:
            continue
        whole = _column(gas_rows[1:], column)
        if not np.any(np.isfinite(whole)):
            continue
        above = np.flatnonzero(whole > max(floor, 1.0))
        if above.size == 0:
            continue
        x, y, z_axis = positions[serial]
        timings.append(
            Timing(
                serial=serial, x=x, y=y, z=z_axis + (release_height - 0.5),
                arrival=(int(above[0]) - lo) * interval,
                peak_time=(int(np.nanargmax(whole)) - lo) * interval,
                departure=(int(above[-1]) - lo) * interval,
                peak=float(np.nanmax(whole)),
            )
        )

    number = int(re.search(r"trial_(\d+)_", path.name).group(1))
    return Trial(
        number=number, path=path,
        flow_peak=float(np.nanmax(flow)),
        flow_mean=float(np.nanmean(flow[first:last])),
        window=(first, last), readings=readings, timings=timings,
        sample_interval=interval, release_start=lo,
    )


#: Far-field stand geometry, from Figure A4 of the report: two arcs, at 10 and
#: 14 m, with the on-axis stand and pairs either side.
FAR_FIELD_STANDS = {
    8: (10.0, 0.0), 7: (10.0, 2.0), 9: (10.0, -2.0),
    6: (10.0, 4.0), 10: (10.0, -4.0),
    3: (14.0, 0.0), 2: (14.0, 2.5), 4: (14.0, -2.5),
    1: (14.0, 5.0), 5: (14.0, -5.0),
}


@dataclass
class ArcFit:
    """A Gaussian fitted across one far-field arc.

    The far-field stands are fixed while the plume follows the wind, and the
    wind direction is recorded only as a compass point every five minutes --
    coarser than the 10 to 12 degrees that separate adjacent stands. So the
    plume cannot be *placed* from the wind record.

    It does not need to be. Five stands across an arc constrain a Gaussian,
    and fitting one recovers the centreline concentration and the lateral
    spread directly from the measurements, with the plume's position as a
    nuisance parameter. The position is fitted out; what is compared is the
    magnitude and the width, which is what the model predicts anyway.
    """

    trial: int
    radius: float  #: arc radius, m
    height: float  #: sensor height, m
    centreline: float  #: fitted peak, vol %
    offset: float  #: fitted plume centre, m across wind
    sigma_y: float  #: fitted lateral standard deviation, m
    r_squared: float
    points: int

    @property
    def well_constrained(self) -> bool:
        """Whether the plume was actually over the array.

        A fit whose centre lands on the edge of the stand pattern is
        extrapolating, and its peak is unbounded above.
        """
        return (
            self.r_squared > 0.85
            and abs(self.offset) < 4.0
            and self.points >= 4
            and self.centreline < 100.0
        )


def fit_arcs(
    far_field_rows,
    *,
    floor: float = 0.02,
    ceiling: float = 3.9,
) -> list[ArcFit]:
    """Fit a Gaussian across each far-field arc.

    Parameters
    ----------
    far_field_rows
        Dicts with ``trial``, ``stand``, ``height_m``, ``h2_max`` and
        ``saturated``, as the reduced far-field table provides.
    ceiling
        Readings at or above this are treated as censored and dropped.  The
        Dräger devices range to 4 vol %, which is hydrogen's lower flammable
        limit, and they compress well before it.
    """
    from collections import defaultdict

    from scipy.optimize import curve_fit

    grouped = defaultdict(dict)
    for row in far_field_rows:
        if str(row.get("saturated", "0")) == "1":
            continue
        try:
            value = float(row["h2_max"])
            stand = int(row["stand"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (floor < value < ceiling) or stand not in FAR_FIELD_STANDS:
            continue
        radius, across = FAR_FIELD_STANDS[stand]
        grouped[(int(row["trial"]), radius, float(row["height_m"]))][across] = value

    out = []
    for (trial, radius, height), points in sorted(grouped.items()):
        if len(points) < 4:
            continue
        across = np.array(sorted(points))
        values = np.array([points[a] for a in across])
        if values.max() < 0.2:
            continue
        try:
            fitted, _ = curve_fit(
                _gaussian, across, values,
                p0=[
                    values.max(), across[int(np.argmax(values))],
                    3.0 / math.sqrt(2.0),
                ],
                # These are the former e-folding-width limits converted to
                # standard deviations, preserving the accepted fit set.
                bounds=(
                    [0.0, -8.0, 0.3 / math.sqrt(2.0)],
                    [100.0, 8.0, 40.0 / math.sqrt(2.0)],
                ),
                maxfev=20000,
            )
        except Exception:
            continue
        residual = np.sum((values - _gaussian(across, *fitted)) ** 2)
        spread = np.sum((values - values.mean()) ** 2)
        out.append(
            ArcFit(
                trial=trial, radius=radius, height=height,
                centreline=float(fitted[0]), offset=float(fitted[1]),
                sigma_y=float(fitted[2]),
                r_squared=float(1.0 - residual / max(spread, 1e-12)),
                points=len(across),
            )
        )
    return out


@dataclass
class VerticalFit:
    """A Gaussian fitted down one near-field arc.

    The array carries four or five heights at each distance, which is enough
    to locate the plume's centre and measure its vertical spread directly.
    That matters because those two are what an integral model computes, and
    reading them separately separates faults that a concentration comparison
    mixes together: near the source this model reads a third high and at 10 to
    14 m several times low, which looks like one confusing signal and is
    actually a trajectory error and a spread error pulling opposite ways.
    """

    trial: int
    x: float  #: downwind distance, m
    centre: float  #: fitted plume centre height above the ground, m
    sigma_z: float  #: fitted vertical standard deviation, m
    peak: float  #: fitted centreline concentration, vol %
    r_squared: float
    points: int

    @property
    def well_constrained(self) -> bool:
        # Preserve the former 3 m e-folding-width rejection limit after
        # expressing the stored width as a standard deviation.
        return (
            self.r_squared > 0.85
            and self.points >= 4
            and self.sigma_z < 3.0 / math.sqrt(2.0)
        )


def fit_vertical(trial: "Trial", *, floor: float = 1.0) -> list[VerticalFit]:
    """Fit a Gaussian in height at each distance on the plume axis."""
    from collections import defaultdict

    from scipy.optimize import curve_fit

    grouped = defaultdict(dict)
    for reading in trial.readings:
        if abs(reading.y) < 0.01 and reading.peak > 0.05:
            grouped[round(reading.x, 2)][round(reading.z, 2)] = reading.peak

    out = []
    for x, levels in sorted(grouped.items()):
        if len(levels) < 4:
            continue
        heights = np.array(sorted(levels))
        values = np.array([levels[h] for h in heights])
        if values.max() < floor:
            continue
        try:
            fitted, _ = curve_fit(
                _gaussian, heights, values,
                p0=[
                    values.max(), heights[int(np.argmax(values))],
                    0.4 / math.sqrt(2.0),
                ],
                # Preserve the former e-folding-width limits after the
                # width-definition correction.
                bounds=(
                    [0.0, -1.0, 0.03 / math.sqrt(2.0)],
                    [110.0, 6.0, 10.0 / math.sqrt(2.0)],
                ),
                maxfev=20000,
            )
        except Exception:
            continue
        residual = float(np.sum((values - _gaussian(heights, *fitted)) ** 2))
        spread = float(np.sum((values - values.mean()) ** 2))
        out.append(
            VerticalFit(
                trial=trial.number, x=x, peak=float(fitted[0]),
                centre=float(fitted[1]), sigma_z=float(fitted[2]),
                r_squared=1.0 - residual / max(spread, 1e-12),
                points=len(heights),
            )
        )
    return out


def load(
    directory: str | Path, report: str | Path, heights: dict | None = None
) -> list[Trial]:
    """Read every trial in a directory, skipping those with no flow record.

    ``heights`` maps trial number to release height in metres; trials not
    listed are taken as 0.5 m.
    """
    positions = sensor_positions(report)
    heights = heights or {}
    out = []
    for path in sorted(Path(directory).glob("trial_*alldata.xlsx")):
        number = int(re.search(r"trial_(\d+)_", path.name).group(1))
        try:
            out.append(
                read_trial(
                    path, positions,
                    release_height=heights.get(number, 0.5),
                )
            )
        except (ValueError, KeyError, StopIteration):
            continue
    return sorted(out, key=lambda t: t.number)
