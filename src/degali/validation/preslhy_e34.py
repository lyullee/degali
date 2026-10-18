"""PRESLHY E3.4 pool mass and substrate-temperature observation operator.

E3.4 workbooks are optional local validation inputs (DOI 10.35097/1319).  No
measurement is shipped with DEGALI.  The reader extracts the experimenter's
corrected ``m(LH2) [g]`` channel rather than reconstructing pool inventory
from the total facility scale.  This matters for sand/gravel records where
solid displacement or material loss is visible in the total-weight signal.

The code intentionally requires a caller-declared evaporation interval.  The
record includes filling, sloshing and, on some substrates, mass artefacts;
selecting a favourable interval automatically from the mass slope would turn
the validation operator into a fitted source term.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re

import numpy as np


@dataclass(frozen=True)
class E34PoolHistory:
    """Corrected pool mass plus co-recorded thermocouple histories."""

    path: Path
    time_s: np.ndarray
    liquid_mass_kg: np.ndarray
    thermocouple_k: dict[str, np.ndarray]
    substrate: str


@dataclass(frozen=True)
class E34EvaporationWindow:
    """Mass-loss observation over a fixed, caller-declared stable interval."""

    start_s: float
    end_s: float
    samples: int
    initial_mass_kg: float
    final_mass_kg: float
    evaporation_rate_kg_s: float
    mass_slope_r_squared: float
    temperature_span_k: float
    notes: tuple[str, ...]

    @property
    def mass_lost_kg(self) -> float:
        return self.initial_mass_kg - self.final_mass_kg


@dataclass(frozen=True)
class E34DynamicWindHydrogenHistory:
    """Synchronized 10 Hz hydrogen and three-component sonic record.

    ``synchronized_time_s`` is the acquisition clock aligned with the first
    thermocouple/scale sheet. ``hydrogen_sample_time_s`` additionally removes
    the campaign-declared approximately 35 s sample-tube transport delay and
    is therefore the physical observation time at the three sampling points.
    Quality flags are read from the workbook's own defective/pink cells rather
    than inferred from the signal values.
    """

    path: Path
    original_time_s: np.ndarray
    synchronized_time_s: np.ndarray
    hydrogen_sample_time_s: np.ndarray
    hydrogen_vol_percent: dict[str, np.ndarray]
    sonic_temperature_k: np.ndarray
    wind_u_m_s: np.ndarray
    wind_v_m_s: np.ndarray
    wind_w_m_s: np.ndarray
    trigger_v: np.ndarray
    sensor_quality: dict[str, str]
    wind_quality: str
    declared_sample_transport_delay_s: float
    median_sample_interval_s: float


@dataclass(frozen=True)
class E34WindHistory:
    """Wind-only view of the synchronized dynamic worksheet.

    This reader exists so source/prehistory preparation can use the sonic
    record without loading any H2 concentration values into memory.
    """

    path: Path
    synchronized_time_s: np.ndarray
    wind_u_m_s: np.ndarray
    wind_v_m_s: np.ndarray
    wind_w_m_s: np.ndarray
    wind_quality: str
    median_sample_interval_s: float


@dataclass(frozen=True)
class E34SourceWindowPolicy:
    """Predeclared mass-only rule for identifying every post-fill decay.

    The policy deliberately has no hydrogen, wind or model-error inputs.  A
    caller can therefore freeze all detected source windows before any
    concentration calculation is run.  Absolute defaults are intentionally
    conservative for the approximately 1.5--5.7 kg public E3.4 records.
    """

    smoothing_sigma_s: float = 10.0
    minimum_peak_prominence_kg: float = 0.25
    relative_peak_prominence: float = 0.08
    minimum_peak_separation_s: float = 150.0
    active_mass_floor_kg: float = 0.10
    relative_end_mass_fraction: float = 0.10
    minimum_mass_loss_kg: float = 0.10
    minimum_duration_s: float = 30.0


@dataclass(frozen=True)
class E34SourceWindowCandidate:
    """One mass-only post-fill interval retained by the frozen policy."""

    window_id: str
    start_s: float
    end_s: float
    duration_s: float
    samples: int
    smoothed_initial_mass_kg: float
    smoothed_final_mass_kg: float
    smoothed_mass_lost_kg: float
    fitted_evaporation_rate_kg_s: float
    raw_mass_slope_r_squared: float
    end_reason: str
    source_only_selection: bool
    manual_review_reasons: tuple[str, ...]


@dataclass(frozen=True)
class E34FillPrehistoryPolicy:
    """Mass-only rule for starting a diagnostic fill prehistory."""

    smoothing_sigma_s: float = 10.0
    relative_fill_onset_fraction: float = 0.10
    minimum_prehistory_s: float = 10.0


@dataclass(frozen=True)
class E34FillPrehistoryCandidate:
    """One mass-only fill onset preceding a frozen post-fill score window."""

    window_id: str
    prehistory_start_s: float
    score_start_s: float
    score_end_s: float
    prehistory_duration_s: float
    score_duration_s: float
    full_duration_s: float
    preceding_baseline_time_s: float
    smoothed_baseline_mass_kg: float
    smoothed_peak_mass_kg: float
    onset_threshold_mass_kg: float
    relative_fill_onset_fraction: float
    source_only_selection: bool
    limitations: tuple[str, ...]


_DYNAMIC_REQUIRED_COLUMNS = (
    "Orig.Time [s]", "Sync. Time [s]", "CH35", "CH45", "CH55",
    "Young_T", "Young_U", "Young_V", "Young_W", "Trigger",
    "cH2 Time [s]",
)

_WIND_REQUIRED_COLUMNS = (
    "Sync. Time [s]", "Young_U", "Young_V", "Young_W",
)


def _campaign_flagged(cell: object) -> bool:
    value = getattr(cell, "value", None)
    if isinstance(value, str) and value.strip().lower() == "defective":
        return True
    fill = getattr(cell, "fill", None)
    if getattr(fill, "fill_type", None) != "solid":
        return False
    colour = getattr(getattr(fill, "fgColor", None), "rgb", None)
    return isinstance(colour, str) and colour.upper() in {
        "FFFFCCFF", "FFFFC7CE", "FFFF0000",
    }


def read_e34_dynamic_wind_hydrogen_history(
    workbook: str | Path,
) -> E34DynamicWindHydrogenHistory:
    """Read the campaign-synchronized H2 and sonic time-series worksheet."""

    import openpyxl

    path = Path(workbook)
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        if len(book.worksheets) < 2:
            raise ValueError(f"{path.name}: dynamic cH2/wind worksheet missing")
        sheet = book.worksheets[1]
        header = next(sheet.iter_rows(min_row=6, max_row=6, values_only=True))
        columns = {
            str(value).strip(): index
            for index, value in enumerate(header) if value is not None
        }
        missing = [name for name in _DYNAMIC_REQUIRED_COLUMNS if name not in columns]
        if missing:
            raise ValueError(
                f"{path.name}: dynamic columns missing: {', '.join(missing)}"
            )
        quality_row = next(sheet.iter_rows(min_row=5, max_row=5))
        sensor_quality = {
            label: (
                "workbook_flagged_defective"
                if _campaign_flagged(quality_row[columns[label]]) else "usable"
            )
            for label in ("CH35", "CH45", "CH55")
        }
        wind_quality = (
            "workbook_flagged_defective"
            if any(
                _campaign_flagged(quality_row[columns[label]])
                for label in ("Young_U", "Young_V", "Young_W")
            )
            else "usable"
        )
        rows = list(sheet.iter_rows(min_row=7, values_only=True))
    finally:
        book.close()

    arrays = {
        name: _finite(
            [row[index] if index < len(row) else None for row in rows]
        )
        for name, index in columns.items() if name in _DYNAMIC_REQUIRED_COLUMNS
    }
    valid = np.logical_and.reduce([
        np.isfinite(arrays[name]) for name in _DYNAMIC_REQUIRED_COLUMNS
    ])
    if np.count_nonzero(valid) < 3:
        raise ValueError(f"{path.name}: fewer than three complete dynamic samples")
    arrays = {name: values[valid] for name, values in arrays.items()}
    original = arrays["Orig.Time [s]"]
    synchronized = arrays["Sync. Time [s]"]
    hydrogen_time = arrays["cH2 Time [s]"]
    if np.any(np.diff(original) <= 0.0):
        raise ValueError(f"{path.name}: dynamic time is not strictly increasing")
    if not np.allclose(
        np.diff(synchronized), np.diff(original), rtol=0.0, atol=1.0e-8,
    ):
        raise ValueError(f"{path.name}: synchronized clock changes sample spacing")
    delay = synchronized - hydrogen_time
    median_delay = float(np.median(delay))
    if not np.allclose(delay, median_delay, rtol=0.0, atol=1.0e-8):
        raise ValueError(f"{path.name}: H2 transport delay is not constant")
    return E34DynamicWindHydrogenHistory(
        path=path,
        original_time_s=original,
        synchronized_time_s=synchronized,
        hydrogen_sample_time_s=hydrogen_time,
        hydrogen_vol_percent={
            label: arrays[label] for label in ("CH35", "CH45", "CH55")
        },
        sonic_temperature_k=arrays["Young_T"],
        wind_u_m_s=arrays["Young_U"],
        wind_v_m_s=arrays["Young_V"],
        wind_w_m_s=arrays["Young_W"],
        trigger_v=arrays["Trigger"],
        sensor_quality=sensor_quality,
        wind_quality=wind_quality,
        declared_sample_transport_delay_s=median_delay,
        median_sample_interval_s=float(np.median(np.diff(original))),
    )


def read_e34_wind_history(workbook: str | Path) -> E34WindHistory:
    """Read only synchronized Young U/V/W; no H2 value column is loaded."""

    import openpyxl

    path = Path(workbook)
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        if len(book.worksheets) < 2:
            raise ValueError(f"{path.name}: dynamic wind worksheet missing")
        sheet = book.worksheets[1]
        header = next(sheet.iter_rows(min_row=6, max_row=6, values_only=True))
        columns = {
            str(value).strip(): index
            for index, value in enumerate(header) if value is not None
        }
        missing = [name for name in _WIND_REQUIRED_COLUMNS if name not in columns]
        if missing:
            raise ValueError(
                f"{path.name}: wind columns missing: {', '.join(missing)}"
            )
        quality_row = next(sheet.iter_rows(min_row=5, max_row=5))
        wind_quality = (
            "workbook_flagged_defective"
            if any(
                _campaign_flagged(quality_row[columns[label]])
                for label in ("Young_U", "Young_V", "Young_W")
            )
            else "usable"
        )
        rows = list(sheet.iter_rows(min_row=7, values_only=True))
    finally:
        book.close()

    arrays = {
        name: _finite([
            row[index] if index < len(row) else None for row in rows
        ])
        for name, index in columns.items() if name in _WIND_REQUIRED_COLUMNS
    }
    valid = np.logical_and.reduce([
        np.isfinite(arrays[name]) for name in _WIND_REQUIRED_COLUMNS
    ])
    if np.count_nonzero(valid) < 3:
        raise ValueError(f"{path.name}: fewer than three complete wind samples")
    arrays = {name: values[valid] for name, values in arrays.items()}
    synchronized = arrays["Sync. Time [s]"]
    if np.any(np.diff(synchronized) <= 0.0):
        raise ValueError(f"{path.name}: synchronized wind time is not increasing")
    return E34WindHistory(
        path=path,
        synchronized_time_s=synchronized,
        wind_u_m_s=arrays["Young_U"],
        wind_v_m_s=arrays["Young_V"],
        wind_w_m_s=arrays["Young_W"],
        wind_quality=wind_quality,
        median_sample_interval_s=float(np.median(np.diff(synchronized))),
    )


def _finite(values) -> np.ndarray:
    out = []
    for value in values:
        try:
            value = float(value)
        except (TypeError, ValueError):
            value = float("nan")
        out.append(value if math.isfinite(value) else float("nan"))
    return np.asarray(out, dtype=float)


def read_e34_pool_history(workbook: str | Path) -> E34PoolHistory:
    """Read one raw E3.4 workbook without altering its time origin.

    The first sheet contains the thermocouples and the scale-derived pool
    record.  The required headings live on row seven in the public files;
    locating them by name is safer than hard-coding their columns.  Mass is
    converted from the supplied corrected grams to kilograms.
    """

    import openpyxl

    path = Path(workbook)
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = book.worksheets[0]
        header = next(sheet.iter_rows(min_row=7, max_row=7, values_only=True))
        try:
            time_column = next(
                i for i, value in enumerate(header)
                if str(value).strip() == "Sync. Time [s]"
            )
        except StopIteration:
            # Concrete01 is the lone public file without a synchronised
            # export column; its thermocouple sample index is exactly seconds.
            time_column = 0
        try:
            mass_column = next(
                i for i, value in enumerate(header)
                if str(value).strip() == "m(LH2) [g]"
            )
        except StopIteration as error:
            raise ValueError(f"{path.name}: corrected m(LH2) channel missing") from error

        tc_columns = {
            str(value).strip(): i
            for i, value in enumerate(header)
            if value is not None and i < 29 and str(value).strip() not in {"", "X_Value"}
        }
        rows = list(sheet.iter_rows(min_row=8, values_only=True))
    finally:
        book.close()

    time_s = _finite([row[time_column] if time_column < len(row) else None for row in rows])
    mass_kg = 1.0e-3 * _finite([row[mass_column] if mass_column < len(row) else None for row in rows])
    valid = np.isfinite(time_s) & np.isfinite(mass_kg)
    if np.count_nonzero(valid) < 3:
        raise ValueError(f"{path.name}: fewer than three synchronized pool-mass samples")
    time_s, mass_kg = time_s[valid], mass_kg[valid]
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError(f"{path.name}: pool time is not strictly increasing")
    thermocouple_k = {
        name: _finite([row[column] if column < len(row) else None for row in rows])[valid]
        for name, column in tc_columns.items()
    }
    name_parts = path.name.split("-", 2)
    raw_substrate = name_parts[1] if len(name_parts) >= 2 else "unknown"
    match = re.match(r"[A-Za-z]+", raw_substrate)
    substrate = match.group(0).lower() if match else raw_substrate.lower()
    if substrate.endswith("wind"):
        substrate = substrate[:-4]
    return E34PoolHistory(
        path=path, time_s=time_s, liquid_mass_kg=mass_kg,
        thermocouple_k=thermocouple_k,
        substrate=substrate,
    )


def evaporation_window(
    history: E34PoolHistory,
    *,
    start_s: float,
    end_s: float,
    thermocouples: tuple[str, ...] = (),
) -> E34EvaporationWindow:
    """Reduce a predeclared post-fill mass-loss interval.

    A negative fitted slope is required.  ``temperature_span_k`` records the
    maximum span of selected thermal channels; it is a screening diagnostic,
    not an automatic acceptance threshold because each substrate and depth
    has a physically different cooling transient.
    """

    if end_s <= start_s:
        raise ValueError("end_s must be later than start_s")
    selected = (history.time_s >= start_s) & (history.time_s <= end_s)
    if np.count_nonzero(selected) < 3:
        raise ValueError("evaporation window needs at least three samples")
    time = history.time_s[selected]
    mass = history.liquid_mass_kg[selected]
    slope, intercept = np.polyfit(time, mass, 1)
    predicted = slope * time + intercept
    variance = float(np.sum((mass - mass.mean()) ** 2))
    r2 = 1.0 - float(np.sum((mass - predicted) ** 2)) / max(variance, 1.0e-30)
    if slope >= 0.0:
        raise ValueError("selected interval has no net evaporative mass loss")
    requested = thermocouples or tuple(history.thermocouple_k)
    missing = [name for name in requested if name not in history.thermocouple_k]
    if missing:
        raise KeyError(f"thermocouples not present: {missing}")
    spans = [
        float(np.nanmax(history.thermocouple_k[name][selected]) - np.nanmin(history.thermocouple_k[name][selected]))
        for name in requested
        if np.any(np.isfinite(history.thermocouple_k[name][selected]))
    ]
    notes = [
        "rate is the least-squares slope of the supplied corrected m(LH2) channel",
        "the time interval is caller-declared; no interval was selected to improve agreement",
    ]
    if history.substrate in {"sand", "gravel"}:
        notes.append("sand/gravel mass artefacts remain a manual review item")
    return E34EvaporationWindow(
        start_s=float(time[0]), end_s=float(time[-1]), samples=int(time.size),
        initial_mass_kg=float(mass[0]), final_mass_kg=float(mass[-1]),
        evaporation_rate_kg_s=float(-slope), mass_slope_r_squared=float(r2),
        temperature_span_k=max(spans, default=float("nan")), notes=tuple(notes),
    )


def source_only_evaporation_candidates(
    history: E34PoolHistory,
    *,
    policy: E34SourceWindowPolicy = E34SourceWindowPolicy(),
) -> tuple[E34SourceWindowCandidate, ...]:
    """Return every source-mass post-fill decay selected by a fixed rule.

    The corrected pool-mass channel is Gaussian-smoothed solely to locate
    prominent fill maxima.  Each candidate begins at a maximum and ends at
    the lowest smoothed mass before the next detected maximum.  The last
    candidate ends when it first reaches a declared residual-mass threshold,
    or at the record boundary when it never does.  No wind, hydrogen sensor,
    temperature or model output is read by this function.

    These are frozen *candidate* source intervals.  They do not by themselves
    prove that liquid-mass loss equals atmospheric hydrogen vapour flux.
    """

    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks

    if history.time_s.size < 3:
        raise ValueError("source-window selection needs at least three samples")
    if policy.smoothing_sigma_s <= 0.0:
        raise ValueError("smoothing_sigma_s must be positive")
    if policy.minimum_peak_prominence_kg <= 0.0:
        raise ValueError("minimum_peak_prominence_kg must be positive")
    if not 0.0 <= policy.relative_peak_prominence <= 1.0:
        raise ValueError("relative_peak_prominence must be between zero and one")
    if policy.minimum_peak_separation_s <= 0.0:
        raise ValueError("minimum_peak_separation_s must be positive")
    if policy.active_mass_floor_kg < 0.0:
        raise ValueError("active_mass_floor_kg must be non-negative")
    if not 0.0 < policy.relative_end_mass_fraction < 1.0:
        raise ValueError("relative_end_mass_fraction must be between zero and one")
    if policy.minimum_mass_loss_kg <= 0.0:
        raise ValueError("minimum_mass_loss_kg must be positive")
    if policy.minimum_duration_s <= 0.0:
        raise ValueError("minimum_duration_s must be positive")

    spacing = np.diff(history.time_s)
    median_spacing = float(np.median(spacing))
    if median_spacing <= 0.0:
        raise ValueError("pool time spacing must be positive")
    # A highly irregular clock makes a sample-count Gaussian ambiguous.  The
    # public workbooks are uniform; fail closed rather than silently resample.
    if not np.allclose(spacing, median_spacing, rtol=1.0e-4, atol=1.0e-8):
        raise ValueError("source-window selection requires a uniform time grid")

    mass = np.asarray(history.liquid_mass_kg, dtype=float)
    sigma_samples = policy.smoothing_sigma_s / median_spacing
    smoothed = gaussian_filter1d(mass, sigma_samples, mode="nearest")
    lower = float(np.quantile(smoothed, 0.05))
    upper = float(np.quantile(smoothed, 0.99))
    excursion = max(upper - lower, 0.0)
    prominence = max(
        policy.minimum_peak_prominence_kg,
        policy.relative_peak_prominence * excursion,
    )
    separation_samples = max(
        1, int(round(policy.minimum_peak_separation_s / median_spacing)),
    )
    peaks, _ = find_peaks(
        smoothed, prominence=prominence, distance=separation_samples,
    )
    if peaks.size == 0:
        raise ValueError("no prominent pool-fill maxima found by source-only policy")

    candidates: list[E34SourceWindowCandidate] = []
    for ordinal, peak_index in enumerate(peaks, start=1):
        next_peak = int(peaks[ordinal]) if ordinal < peaks.size else None
        end_reason: str
        if next_peak is not None:
            # The minimum before the next fill is deterministic even when the
            # pool did not fully empty between successive top-ups.
            search = smoothed[peak_index:next_peak]
            end_index = int(peak_index + np.argmin(search))
            end_reason = "minimum_before_next_detected_fill"
        else:
            residual_threshold = lower + max(
                policy.active_mass_floor_kg,
                policy.relative_end_mass_fraction
                * max(float(smoothed[peak_index]) - lower, 0.0),
            )
            crossings = np.flatnonzero(
                smoothed[peak_index + 1:] <= residual_threshold
            )
            if crossings.size:
                end_index = int(peak_index + 1 + crossings[0])
                end_reason = "declared_residual_mass_threshold"
            else:
                end_index = int(smoothed.size - 1)
                end_reason = "record_boundary_before_residual_threshold"

        duration = float(history.time_s[end_index] - history.time_s[peak_index])
        smoothed_loss = float(smoothed[peak_index] - smoothed[end_index])
        if duration < policy.minimum_duration_s:
            continue
        if smoothed_loss < policy.minimum_mass_loss_kg:
            continue
        reduced = evaporation_window(
            history,
            start_s=float(history.time_s[peak_index]),
            end_s=float(history.time_s[end_index]),
        )
        review: list[str] = [
            "liquid_mass_loss_is_not_a_direct_measurement_of_atmospheric_vapour_flux"
        ]
        if history.substrate in {"sand", "gravel"}:
            review.append("porous_substrate_mass_artefact_requires_review")
        if end_reason == "record_boundary_before_residual_threshold":
            review.append("decay_is_right_censored_by_record_boundary")
        if reduced.mass_slope_r_squared < 0.90:
            review.append("raw_mass_decay_is_not_well_represented_by_one_linear_rate")
        candidates.append(E34SourceWindowCandidate(
            window_id=f"{history.path.stem}:source-window-{ordinal:02d}",
            start_s=float(history.time_s[peak_index]),
            end_s=float(history.time_s[end_index]),
            duration_s=duration,
            samples=int(end_index - peak_index + 1),
            smoothed_initial_mass_kg=float(smoothed[peak_index]),
            smoothed_final_mass_kg=float(smoothed[end_index]),
            smoothed_mass_lost_kg=smoothed_loss,
            fitted_evaporation_rate_kg_s=reduced.evaporation_rate_kg_s,
            raw_mass_slope_r_squared=reduced.mass_slope_r_squared,
            end_reason=end_reason,
            source_only_selection=True,
            manual_review_reasons=tuple(review),
        ))
    if not candidates:
        raise ValueError("source-only policy found no qualifying decay windows")
    return tuple(candidates)


def source_only_fill_prehistory_candidate(
    history: E34PoolHistory,
    *,
    window_ordinal: int,
    source_window_policy: E34SourceWindowPolicy = E34SourceWindowPolicy(),
    prehistory_policy: E34FillPrehistoryPolicy = E34FillPrehistoryPolicy(),
) -> E34FillPrehistoryCandidate:
    """Locate a fill onset before one frozen post-fill window using mass only.

    The local baseline is the lowest smoothed corrected mass after the prior
    frozen window (or the record start for window one) and before the selected
    fill maximum. The onset is the first upward crossing of a declared fraction
    of the baseline-to-peak excursion. This identifies timing only. It does not
    infer atmospheric vapour flow while the vessel is filling.
    """

    from scipy.ndimage import gaussian_filter1d

    if isinstance(window_ordinal, bool) or not isinstance(window_ordinal, int):
        raise TypeError("window_ordinal must be an integer")
    if window_ordinal < 1:
        raise ValueError("window_ordinal must be at least one")
    if prehistory_policy.smoothing_sigma_s <= 0.0:
        raise ValueError("prehistory smoothing_sigma_s must be positive")
    if not 0.0 < prehistory_policy.relative_fill_onset_fraction < 1.0:
        raise ValueError("relative_fill_onset_fraction must lie between zero and one")
    if prehistory_policy.minimum_prehistory_s <= 0.0:
        raise ValueError("minimum_prehistory_s must be positive")

    windows = source_only_evaporation_candidates(
        history, policy=source_window_policy,
    )
    if window_ordinal > len(windows):
        raise ValueError("window_ordinal is outside the frozen source-window set")
    selected = windows[window_ordinal - 1]
    spacing = float(np.median(np.diff(history.time_s)))
    smoothed = gaussian_filter1d(
        np.asarray(history.liquid_mass_kg, dtype=float),
        prehistory_policy.smoothing_sigma_s / spacing,
        mode="nearest",
    )
    peak_index = int(np.argmin(np.abs(history.time_s - selected.start_s)))
    if not math.isclose(
        float(history.time_s[peak_index]), selected.start_s, abs_tol=spacing / 2.0,
    ):
        raise ValueError("selected source-window start is not on the mass clock")
    if window_ordinal == 1:
        lower_index = 0
    else:
        previous_end = windows[window_ordinal - 2].end_s
        lower_index = int(np.searchsorted(history.time_s, previous_end, side="left"))
    if lower_index >= peak_index:
        raise ValueError("no pre-fill mass history precedes the selected peak")
    local = smoothed[lower_index:peak_index + 1]
    baseline_index = int(lower_index + np.argmin(local))
    baseline_mass = float(smoothed[baseline_index])
    peak_mass = float(smoothed[peak_index])
    excursion = peak_mass - baseline_mass
    if excursion <= 0.0:
        raise ValueError("selected window has no positive preceding fill excursion")
    threshold = baseline_mass + (
        prehistory_policy.relative_fill_onset_fraction * excursion
    )
    crossings = np.flatnonzero(
        smoothed[baseline_index:peak_index + 1] >= threshold
    )
    if not crossings.size:
        raise ValueError("declared fill-onset threshold is not crossed")
    onset_index = int(baseline_index + crossings[0])
    onset_time = float(history.time_s[onset_index])
    prehistory_duration = selected.start_s - onset_time
    if prehistory_duration < prehistory_policy.minimum_prehistory_s:
        raise ValueError("mass-only fill prehistory is shorter than declared minimum")
    limitations = [
        "fill_onset_timing_uses_only_smoothed_corrected_pool_mass",
        "positive_liquid_fill_slope_is_not_atmospheric_vapour_flow",
        "prehistory_source_rate_requires_a_separate_declared_assumption",
        "candidate_is_diagnostic_and_cannot_replace_the_failed_primary_score",
    ]
    if baseline_mass < 0.0:
        limitations.append("corrected_mass_baseline_is_negative")
    if history.substrate in {"sand", "gravel"}:
        limitations.append("porous_substrate_mass_artefact_requires_review")
    return E34FillPrehistoryCandidate(
        window_id=selected.window_id,
        prehistory_start_s=onset_time,
        score_start_s=selected.start_s,
        score_end_s=selected.end_s,
        prehistory_duration_s=prehistory_duration,
        score_duration_s=selected.duration_s,
        full_duration_s=selected.end_s - onset_time,
        preceding_baseline_time_s=float(history.time_s[baseline_index]),
        smoothed_baseline_mass_kg=baseline_mass,
        smoothed_peak_mass_kg=peak_mass,
        onset_threshold_mass_kg=threshold,
        relative_fill_onset_fraction=(
            prehistory_policy.relative_fill_onset_fraction
        ),
        source_only_selection=True,
        limitations=tuple(limitations),
    )


__all__ = [
    "E34PoolHistory", "E34EvaporationWindow", "E34DynamicWindHydrogenHistory",
    "E34WindHistory", "E34SourceWindowPolicy", "E34SourceWindowCandidate",
    "E34FillPrehistoryPolicy", "E34FillPrehistoryCandidate",
    "read_e34_pool_history", "read_e34_dynamic_wind_hydrogen_history",
    "read_e34_wind_history", "evaporation_window",
    "source_only_evaporation_candidates", "source_only_fill_prehistory_candidate",
]
