"""Reading the REDIPHEM field-trial database.

REDIPHEM (Risø National Laboratory) collects the heavy-gas dispersion field
trials of the 1980s and 90s into one format: Burro and Coyote (LNG pool
spills), Desert Tortoise and FLADIS (pressurised ammonia jets), Eagle
(nitrogen tetroxide), Thorney Island (Freon puffs), Lathen and the wind-tunnel
series.  Each trial is a directory holding

===============  ==========================================================
``SPECS.DAT``    release and meteorological conditions, one per line
``SETUP.DAT``    sensor positions and channel types
``DATA.DBF``     the measurements
``CHANDEF.DAT``  what each channel type is, one level up
===============  ==========================================================

``DATA.DBF`` is not a dBase file despite the extension.  It is a raw
little-endian ``float32`` stream: a channel count, then that many channel
identifiers, then a matrix of one row per sample whose first column is the
time of day in seconds.  ``-1234`` marks a missing reading.

The point of this module is that the trial conditions map almost one to one
onto a DEGADIS input deck, so :func:`to_case` can build one without anybody
transcribing numbers by hand -- which is where model-evaluation studies
usually go wrong.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Channel type numbers are *per series*: Burro uses 21 and 22 for
#: concentration, Eagle 21 and 24.  There is no global list, and guessing one
#: is worse than refusing: FLADIS ships no ``CHANDEF.DAT``, and a plausible
#: guess there selects channels reading 302 and 23.5 at 20 m, which are not
#: concentrations at all and produce a clean-looking but meaningless result.
CONCENTRATION_CHANNELS: tuple[int, ...] = ()

#: Sentinel for a missing reading.
MISSING = -1234.0

#: Pasquill-Gifford class letter to the 1-based index DEGADIS uses.  Trials
#: recorded as a boundary case ("D-E") take the more stable of the two, which
#: is the conservative reading for a dispersion calculation.
_STABILITY = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6}


def _number(text: str) -> float | None:
    """First number in a REDIPHEM field.

    Values often carry a qualifier -- ``-164 est``, ``0 approx note``,
    ``15 spurious`` -- and a few use a comma as the decimal separator.
    """
    if not text:
        return None
    cleaned = text.replace(",", ".")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", cleaned)
    return float(m.group()) if m else None


@dataclass
class Trial:
    """One REDIPHEM trial."""

    series: str
    name: str
    path: Path
    specs: dict[str, str]
    #: Qualifiers stripped from numeric fields, e.g. ``{"exit temperature":
    #: "est"}``. Worth surfacing: an "estimated" release temperature is a
    #: different kind of input from a measured one.
    qualifiers: dict[str, str] = field(default_factory=dict)

    # -- typed access to the specification ---------------------------------

    def value(self, key: str) -> float | None:
        for k, v in self.specs.items():
            if k.startswith(key):
                return _number(v)
        return None

    def text(self, key: str) -> str:
        for k, v in self.specs.items():
            if k.startswith(key):
                return v
        return ""

    @property
    def substance(self) -> str:
        return self.text("substance").split()[0] if self.text("substance") else ""

    @property
    def release_type(self) -> str:
        return (self.text("release type") or "").split()[0].lower()

    @property
    def stability(self) -> int | None:
        letter = (self.text("stability class") or "").strip().upper()
        if not letter:
            return None
        # "D-E" and similar: take the more stable end
        letters = [c for c in letter if c in _STABILITY]
        return max(_STABILITY[c] for c in letters) if letters else None

    @property
    def wind_height(self) -> float | None:
        """Anemometer height, m.

        Usually its own field.  Thorney Island instead puts it in the *name*
        of the wind-speed field -- ``site average windspeed at 10m`` -- so
        that spelling is read too.
        """
        explicit = self.value("reference height for wind")
        if explicit is None:
            explicit = self.value("reference height")
        if explicit is not None:
            return explicit
        for key in self.specs:
            if "windspeed" in key.lower():
                m = re.search(r"at\s*(\d+(?:\.\d+)?)\s*m", key, re.IGNORECASE)
                if m:
                    return float(m.group(1))
        return None

    @property
    def duration(self) -> float | None:
        return self.value("release duration")

    @property
    def rate(self) -> float | None:
        return self.value("release rate")

    def __repr__(self) -> str:
        return f"<Trial {self.series}/{self.name} {self.substance} {self.release_type}>"

    # -- measurements -------------------------------------------------------

    def channel_types(self) -> set[int]:
        """Channel types that carry a concentration, from ``CHANDEF.DAT``.

        The file sits one level above the trial and lists, in blocks of five
        lines, a type number followed by the quantity, its unit and the
        instrument.  Series that do not ship one fall back on
        :data:`CONCENTRATION_CHANNELS`.
        """
        path = self.path.parent / "CHANDEF.DAT"
        if not path.exists():
            return set()
        lines = path.read_text(encoding="latin1").replace("\r", "").split("\n")
        found: set[int] = set()
        for i in range(0, len(lines) - 1, 5):
            head = lines[i].split()
            if not head or not head[0].lstrip("-").isdigit():
                continue
            if "concentration" in lines[i + 1].strip().lower():
                found.add(int(head[0]))
        return found

    def sensors(self) -> dict[int, tuple[float, float, float, int]]:
        """``SETUP.DAT``: channel id -> ``(x, y, z, channel type)``, metres."""
        out: dict[int, tuple[float, float, float, int]] = {}
        text = (self.path / "SETUP.DAT").read_text(encoding="latin1")
        # some files carry the DOS end-of-file byte on the last data line
        for line in text.replace("\r", "").replace("\x1a", "").split("\n"):
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                out[int(parts[0])] = (
                    float(parts[1]), float(parts[2]), float(parts[3]), int(parts[4])
                )
        return out

    def data(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``DATA.DBF``: ``(channel ids, time, values)``.

        ``time`` is seconds since midnight, as recorded.  Missing readings
        come back as ``nan``.
        """
        a = np.fromfile(self.path / "DATA.DBF", dtype="<f4")
        n = int(a[0]) - 1
        ids = a[1 : 1 + n].astype(int)
        block = a[1 + n :]
        rows = block.size // (n + 1)
        block = block[: rows * (n + 1)].reshape(rows, n + 1)
        values = np.where(block[:, 1:] <= MISSING + 1.0, np.nan, block[:, 1:])
        return ids, block[:, 0].astype(float), values

    @property
    def start_time(self) -> float | None:
        """Release start, seconds since midnight, on the data's own clock.

        The hour, minute and second are three separate lines, and only the
        first carries a label -- the other two are indented continuations. So
        they have to be read by position, not by name: taking the labelled
        line alone gives the hour and silently loses 37 minutes.
        """
        keys = list(self.specs)
        for i, key in enumerate(keys):
            if not key.startswith("release starttime"):
                continue
            parts = [_number(self.specs[k]) for k in keys[i : i + 3]]
            if parts[0] is None:
                return None
            h, m, sec = (p or 0.0 for p in parts)
            return h * 3600.0 + m * 60.0 + sec
        return None

    #: Arc maxima carried by a reduction, keyed by height. Set by
    #: :func:`from_reduced`; ``None`` for a trial read from the archive.
    _reduced_arcs: dict | None = None

    def arc_maxima(
        self, *, height: float | None = None, tolerance: float = 0.05,
        averaging: float = 0.0,
    ) -> dict[float, float]:
        """Peak concentration at each downwind distance, mole per cent.

        Parameters
        ----------
        height
            Restrict to sensors at this elevation.  Trials instrument several
            heights on the same mast, and comparing a ground-level model
            prediction against whichever height happened to read highest is a
            common way to manufacture agreement.
        averaging
            Report the largest *running* mean of this length rather than the
            instantaneous peak.  Field data is sampled at 1 s while a model
            reports a concentration averaged over its own averaging time, so
            the two are not comparable untouched.

            A running mean, not a fixed block average: the cloud passes in a
            few tens of seconds within a record several minutes long, so
            fixed blocks straddle the arrival and dilute the peak with the
            empty record either side of it.
        """
        if self._reduced_arcs is not None:
            # a reduction carries the maxima already computed, at the
            # averaging time recorded in the file. Asking for a different one
            # is a question the reduction cannot answer, so it says so rather
            # than silently returning the wrong average.
            if height is None:
                raise ValueError(
                    "a reduced trial carries maxima per height; give one"
                )
            for z, arcs in self._reduced_arcs.items():
                if abs(z - height) <= tolerance:
                    return dict(sorted(arcs.items()))
            return {}

        ids, time, values = self.data()
        setup = self.sensors()
        kinds = self.channel_types()
        dt = float(np.median(np.diff(time))) if time.size > 1 else 1.0
        window = max(int(round(averaging / dt)), 1) if averaging > dt else 1

        out: dict[float, list[float]] = {}
        for column, channel in enumerate(ids):
            info = setup.get(int(channel))
            if info is None or info[3] not in kinds:
                continue
            x, _y, z, _kind = info
            if height is not None and abs(z - height) > tolerance:
                continue
            series = values[:, column]
            if np.all(np.isnan(series)):
                continue
            out.setdefault(x, []).append(_running_peak(series, window))
        return {x: max(v) for x, v in sorted(out.items()) if x > 0.0}

    def arc_maxima_by_height(
        self, *, averaging: float = 0.0
    ) -> dict[tuple[float, float], float]:
        """Peak concentration keyed by ``(distance, height)``.

        Comparing only one height and correcting the model down to it with the
        model's own vertical profile is circular: the profile is doing the
        work and never gets tested. Using every instrumented height breaks
        that, and multiplies the sample size by the number of levels on the
        mast.
        """
        out: dict[tuple[float, float], float] = {}
        for h in self.sensor_heights():
            for x, c in self.arc_maxima(height=h, averaging=averaging).items():
                out[(x, h)] = c
        return out

    def mast_profiles(
        self, *, averaging: float = 0.0, minimum_levels: int = 2
    ) -> dict[tuple[float, float], dict[float, float]]:
        """Peak concentration by elevation, grouped by mast.

        Keyed by ``(x, y)``, so every value in one entry comes from the *same*
        position on the ground.  This is the only way to test a vertical
        profile: taking the arc maximum at each height separately mixes masts,
        and can put the highest reading at 8 m and the lowest at 3 m simply
        because different sensors were nearest the plume centre.
        """
        if self._reduced_arcs is not None:
            raise NotImplementedError(
                "a reduction carries arc maxima only; this needs the archive"
            )

        ids, time, values = self.data()
        setup = self.sensors()
        kinds = self.channel_types()
        dt = float(np.median(np.diff(time))) if time.size > 1 else 1.0
        window = max(int(round(averaging / dt)), 1) if averaging > dt else 1

        out: dict[tuple[float, float], dict[float, float]] = {}
        for column, channel in enumerate(ids):
            info = setup.get(int(channel))
            if info is None or info[3] not in kinds:
                continue
            x, y, z, _kind = info
            if x <= 0.0:
                continue
            series = values[:, column]
            if np.all(np.isnan(series)):
                continue
            peak = _running_peak(series, window)
            if not np.isfinite(peak):
                continue
            out.setdefault((x, y), {})[z] = max(
                peak, out.setdefault((x, y), {}).get(z, 0.0)
            )
        return {k: v for k, v in sorted(out.items()) if len(v) >= minimum_levels}

    def crosswind_profiles(
        self, *, height: float = 1.0, tolerance: float = 0.05,
        averaging: float = 0.0, minimum_sensors: int = 3,
    ) -> dict[float, dict[float, float]]:
        """Peak concentration across each arc, keyed by distance then ``y``.

        The lateral counterpart of :meth:`mast_profiles`.  All the sensors in
        one entry sit at the same downwind distance and the same height, so
        the spread across them is the plume's crosswind structure rather than
        an artefact of comparing different positions.
        """
        if self._reduced_arcs is not None:
            raise NotImplementedError(
                "a reduction carries arc maxima only; this needs the archive"
            )

        ids, time, values = self.data()
        setup = self.sensors()
        kinds = self.channel_types()
        dt = float(np.median(np.diff(time))) if time.size > 1 else 1.0
        window = max(int(round(averaging / dt)), 1) if averaging > dt else 1

        out: dict[float, dict[float, float]] = {}
        for column, channel in enumerate(ids):
            info = setup.get(int(channel))
            if info is None or info[3] not in kinds:
                continue
            x, y, z, _kind = info
            if x <= 0.0 or abs(z - height) > tolerance:
                continue
            series = values[:, column]
            if np.all(np.isnan(series)):
                continue
            peak = _running_peak(series, window)
            if np.isfinite(peak):
                arc = out.setdefault(x, {})
                arc[y] = max(peak, arc.get(y, 0.0))
        return {x: v for x, v in sorted(out.items())
                if len(v) >= minimum_sensors}

    def time_series(
        self, *, height: float | None = None, tolerance: float = 0.05
    ) -> dict[tuple[float, float], tuple[np.ndarray, np.ndarray]]:
        """Concentration histories, keyed by ``(distance, height)``.

        Returns seconds *since the release started* against mole per cent, so
        arrival and departure times can be compared as well as magnitudes.
        A trial that records no release time returns times as given.
        """
        if self._reduced_arcs is not None:
            raise NotImplementedError(
                "a reduction carries arc maxima only; this needs the archive"
            )

        ids, time, values = self.data()
        setup = self.sensors()
        kinds = self.channel_types()
        origin = self.start_time
        if origin is None or not (time[0] - 600.0 <= origin <= time[-1]):
            # the recorded start does not lie on the data's clock; fall back
            # to the start of the record rather than inventing an offset
            origin = float(time[0])
        out: dict[tuple[float, float], tuple[np.ndarray, np.ndarray]] = {}
        for column, channel in enumerate(ids):
            info = setup.get(int(channel))
            if info is None or info[3] not in kinds:
                continue
            x, _y, z, _kind = info
            if x <= 0.0:
                continue
            if height is not None and abs(z - height) > tolerance:
                continue
            series = values[:, column]
            if np.all(np.isnan(series)):
                continue
            key = (x, z)
            if key in out and np.nanmax(out[key][1]) >= np.nanmax(series):
                continue  # keep the sensor that saw most of the cloud
            out[key] = (time - origin, series)
        return out

    def sensor_heights(self) -> list[float]:
        setup = self.sensors()
        kinds = self.channel_types()
        return sorted(
            {z for (_x, _y, z, kind) in setup.values() if kind in kinds}
        )


def _running_peak(series: np.ndarray, window: int) -> float:
    """Largest running mean of ``window`` samples, ignoring gaps."""
    if window <= 1:
        return float(np.nanmax(series))
    filled = np.nan_to_num(series, nan=0.0)
    valid = (~np.isnan(series)).astype(float)
    kernel = np.ones(window)
    total = np.convolve(filled, kernel, mode="valid")
    count = np.convolve(valid, kernel, mode="valid")
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(count > 0, total / np.maximum(count, 1.0), np.nan)
    return float(np.nanmax(means)) if np.any(np.isfinite(means)) else float("nan")


def read_specs(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Parse a ``SPECS.DAT``, separating values from their qualifiers."""
    specs: dict[str, str] = {}
    qualifiers: dict[str, str] = {}
    text = path.read_text(encoding="latin1").replace("\r", "")
    for line in text.split("\n"):
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key, value = key.strip(), value.strip()
        specs[key] = value
        number = _number(value)
        if number is not None:
            trailing = value[value.find(str(int(number)) if number == int(number)
                                        else "") :]
            words = [w for w in re.sub(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", " ", value).split()
                     if w]
            if words:
                qualifiers[key] = " ".join(words)
    return specs, qualifiers


#: Environment variable naming the unpacked REDIPHEM database.  The data is
#: not redistributable with this package, so the tests that use it skip when
#: it is absent.
ENV_VAR = "REDIPHEM_ROOT"


# ==========================================================================
# reduction
# ==========================================================================

#: Specification fields a prediction needs, in the spelling ``Trial.value``
#: matches on. A reduction that omits any of these carries the observations
#: without the conditions, and a statistic needs both.
#:
#: This list is not decorative. A reduced REDIPHEM table was written to make
#: ``Burro MG 0.811`` reproducible without the archive, and it could not,
#: because it carried 840 arc maxima and none of these. Enumerating them here
#: means the next reduction is checked against the requirement rather than
#: against someone's memory of it.
CASE_FIELDS = (
    "pool diameter",
    "release rate",
    "release duration",
    "site average windspeed",
    "reference height for wind",
    "surface roughness",
    "Monin-Obukov length",
    "ambient temperature",
    "ambient pressure",
    "relative humidity",
    "exit temperature",
    "nozzle diameter",
    "initial concentration",
)


def reduce(root=None, *, heights=(0.1, 1.0, 3.0, 8.0), averaging=18.4) -> dict:
    """Reduce the archive to what the reported statistics need.

    Carries, per trial, both halves of a comparison: the arc maxima at each
    instrumented height, **and** the specification fields a prediction is
    built from. A reduction with only the first is a table of measurements.

    The result is JSON-serialisable and a few tens of kilobytes, so it can
    travel with the repository while the archive -- which is not ours to
    redistribute -- does not.
    """
    trials = []
    for trial in load(root):
        arcs = {}
        for z in heights:
            maxima = trial.arc_maxima(height=z, averaging=averaging)
            if maxima:
                arcs[f"z={z}"] = {str(k): v for k, v in sorted(maxima.items())}
        trials.append({
            "series": trial.series,
            "name": trial.name,
            "substance": trial.substance,
            "release_type": trial.release_type,
            "stability": trial.stability,
            "conditions": {
                key: trial.value(key) for key in CASE_FIELDS
            },
            "qualifiers": dict(trial.qualifiers),
            "arcs": arcs,
        })
    return {
        "source": "REDIPHEM archive (JRC); not redistributable, this is a reduction",
        "note": (
            "arc maxima and the specification fields a prediction needs. "
            "An earlier reduction carried only the maxima and could not "
            "reproduce a single statistic."
        ),
        "averaging_time_s": averaging,
        "trials": trials,
    }


def cases_from_reduced(path) -> list:
    """Read a complete reduction as ready-to-run cases.

    The complete reduction carries the built `Case` rather than the trial
    specification it was built from, which is the right thing to store: the
    specification needs `to_case` to interpret it, and `to_case` makes
    judgements -- which source route a release type takes, what to assume for
    a missing pool diameter -- that should be made once and recorded, not
    remade by every reader.

    Trials `to_case` refused carry a `case_error` instead, and are returned
    with `usable` false and the reason attached. Those are the trials the
    original statistics skipped, so skipping them here reproduces the
    published population rather than silently choosing a different one.
    """
    import json

    import numpy as np

    from ..core.thermo import AmbientConditions, GasProperties
    from ..io.inp import Case, SourceTable
    from .trialcase import TrialCase

    data = json.loads(Path(path).read_text())
    out = []
    for record in data["trials"]:
        trial = Trial(
            series=record["series"], name=record["name"], path=Path(path),
            specs={"release type": record.get("release_type") or ""},
        )
        trial._reduced_arcs = {
            float(k.removeprefix("z=")): {float(d): v for d, v in arcs.items()}
            for k, arcs in (record.get("arcs") or {}).items()
        }
        if record.get("case_error"):
            out.append(TrialCase(
                trial=trial, case=None, substance=None,
                problems=[record["case_error"]],
            ))
            continue

        c, g, a, st = (
            record["case"], record["gas"], record["ambient"],
            record["source_term"],
        )
        source = SourceTable(
            time=np.asarray(st["time"], float),
            rate=np.asarray(st["rate"], float),
            radius=np.asarray(st["radius"], float),
            wc=np.asarray(st["wc"], float),
            temp=np.asarray(st["temp"], float),
            fracv=np.asarray(st["fracv"], float),
            enthalpy=np.asarray(st["enthalpy"], float),
            rho=np.asarray(st["rho"], float),
        )
        # The dispersion coefficients are derived, not stored: they follow
        # from roughness, stability class and averaging time, and recomputing
        # them keeps the reduction from carrying a stale copy. The
        # Monin-Obukhov length *is* stored, because the archive measures it
        # and it does not follow from the class.
        from ..core.atmosphere import stability_defaults

        import dataclasses as _dc

        defaults = _dc.replace(
            stability_defaults(
                c["zr"], "ABCDEF"[int(c["istab"]) - 1], c["avtime"]
            ),
            rml=c["rml"],
        )

        case = Case(
            titles=[f"{record['series']} {record['name']}"],
            stability=defaults,
            u0=c["u0"], z0=c["z0"], zr=c["zr"], istab=int(c["istab"]),
            oodist=c["oodist"], avtime=c["avtime"], indvel=2, rml=c["rml"],
            # The reduction stores four gas fields. The other five --
            # `ulc`, `llc`, `zzc`, `cpk`, `cpp` -- are not in it, and their
            # dataclass defaults are not usable: with `llc` left at its
            # default the flammable-mass integral evaluates
            # `log(cc / clow)` at 28 and `SERIES` overflows, which is how
            # this was found. They are levels of concern and a heat-capacity
            # correlation, properties of the substance rather than of the
            # trial, so they come from the substance table.
            gas=_gas_properties(g),
            ambient=AmbientConditions(
                tamb=a["tamb"], pamb=a["pamb"], humid=a["humid"],
                tsurf=a["tsurf"], ihtfl=int(a["ihtfl"]),
                iwtfl=int(a["iwtfl"]),
            ),
            relhum=c["relhum"], yclow=c["yclow"], gmass0=c["gmass0"],
            source=source, steady_state=bool(c["steady_state"]),
            instantaneous=bool(c["instantaneous"]),
        )
        out.append(TrialCase(trial=trial, case=case, substance=None))
    return out


def _gas_properties(g: dict):
    """Fill a reduction's four gas fields out to the nine a case needs."""
    from ..core.thermo import GasProperties
    from .trialcase import SUBSTANCES

    known = next(
        (s for s in SUBSTANCES.values() if s.label == g["name"]), None
    ) or next(
        (s for s in SUBSTANCES.values() if abs(s.mw - g["mw"]) < 0.05), None
    )
    if known is None:
        raise ValueError(
            f"no levels of concern for {g['name']!r} (mw {g['mw']}); "
            "add it to trialcase.SUBSTANCES rather than letting the "
            "dataclass defaults through"
        )
    return GasProperties(
        name=g["name"], mw=g["mw"], temp=g["temp"], rho=g["rho"],
        ulc=known.ulc, llc=known.llc, cpk=known.cpk, cpp=known.cpp,
    )


def from_reduced(path) -> list:
    """Read a reduction back as objects `trialcase` can use.

    The specification is rebuilt as a ``specs`` dict with the same keys the
    archive uses, so every existing reader path works unchanged: the reduction
    is a different *transport*, not a different data model.
    """
    import json

    data = json.loads(Path(path).read_text())
    out = []
    for record in data["trials"]:
        specs = {}
        for key, value in (record.get("conditions") or {}).items():
            if value is not None:
                specs[key] = str(value)
        if record.get("substance"):
            specs["substance"] = record["substance"]
        if record.get("release_type"):
            specs["release type"] = record["release_type"]
        if record.get("stability") is not None:
            specs["stability class"] = "ABCDEF"[int(record["stability"]) - 1]
        trial = Trial(
            series=record["series"], name=record["name"], path=Path(path),
            specs=specs, qualifiers=dict(record.get("qualifiers") or {}),
        )
        trial._reduced_arcs = {
            float(k.removeprefix("z=")): {float(d): v for d, v in arcs.items()}
            for k, arcs in record["arcs"].items()
        }
        out.append(trial)
    return out

def default_root() -> Path | None:
    """The database location from the environment, if it is set and exists."""
    value = os.environ.get(ENV_VAR)
    if not value:
        return None
    path = Path(value)
    return path if path.exists() else None


def load(root: str | Path | None = None) -> list[Trial]:
    """Find every trial under ``root``.

    ``root`` is the directory containing ``REDIPHEM/DATA``, or any ancestor.
    Defaults to ``$REDIPHEM_ROOT``.
    """
    if root is None:
        root = default_root()
        if root is None:
            raise FileNotFoundError(
                f"set {ENV_VAR} to the unpacked REDIPHEM database, or pass a path"
            )
    root = Path(root)
    trials: list[Trial] = []
    for spec_path in sorted(root.rglob("SPECS.DAT")):
        directory = spec_path.parent
        if not (directory / "DATA.DBF").exists():
            continue
        specs, qualifiers = read_specs(spec_path)
        trials.append(
            Trial(
                series=directory.parent.name,
                name=directory.name,
                path=directory,
                specs=specs,
                qualifiers=qualifiers,
            )
        )
    return trials
