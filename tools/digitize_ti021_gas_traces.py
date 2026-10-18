"""Digitize representative TI021 0.6 s gas-sensor traces and compare DEGALI.

The scanned report plots concentration in vol.% against ``seconds*10**-2``;
one plotted abscissa unit is therefore 100 s.  This script uses only the
printed axes to recover the trace.  It does not scale a trace to the workbook
dose.  The independently tabulated dose is used afterwards as a raster-
digitization check.

The selected pages span the near-wall array and progressively farther/higher
receptors.  Their workbook rows are matched by printed site x/y/z coordinates,
not by assuming that the report's G-number equals workbook row order.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class TracePage:
    pdf_page: int
    gas_id: str
    site_x_m: float
    site_y_m: float
    height_m: float
    ordinate_max_vol_percent: float


TRACE_PAGES = (
    TracePage(17, "G01", 386.0, 214.0, 0.4, 100.0),
    TracePage(20, "G04", 400.0, 220.0, 0.4, 100.0),
    TracePage(21, "G05", 400.0, 220.0, 2.4, 50.0),
    TracePage(55, "G39", 400.0, 275.0, 2.4, 5.0),
    TracePage(59, "G43", 400.0, 275.0, 4.4, 2.0),
    TracePage(67, "G51", 300.0, 300.0, 6.4, 1.0),
    TracePage(89, "G73", 400.0, 350.0, 4.4, 1.0),
    TracePage(100, "G84", 300.0, 500.0, 2.4, 1.0),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _find_sensor_index(dataset: object, page: TracePage) -> int:
    matches = [
        index
        for index, sensor in enumerate(dataset.dose_sensors)
        if (
            sensor.site_x_m == page.site_x_m
            and sensor.site_y_m == page.site_y_m
            and sensor.height_m == page.height_m
        )
    ]
    if len(matches) != 1:
        raise ValueError(f"{page.gas_id}: expected one workbook coordinate match")
    return matches[0]


def _detect_axes(binary: np.ndarray) -> tuple[int, int, int, int]:
    """Locate the printed -2..15 / positive-ordinate plot rectangle."""

    x_left = 850 + int(np.argmax(binary[250:1700, 850:1030].sum(axis=0)))
    x_right = 2750 + int(np.argmax(binary[250:1700, 2750:3000].sum(axis=0)))
    y_top = 220 + int(np.argmax(binary[220:500, x_left:x_right + 1].sum(axis=1)))
    y_zero = 1250 + int(
        np.argmax(binary[1250:1750, x_left:x_right + 1].sum(axis=1))
    )
    if not (850 <= x_left < 1030 and 2750 <= x_right < 3000):
        raise ValueError("gas plot abscissa axes were not found")
    if not (220 <= y_top < 500 and 1250 <= y_zero < 1750):
        raise ValueError("gas plot ordinate axes were not found")
    return x_left, x_right, y_top, y_zero


def _digitize_trace(
    image_path: Path, ordinate_max: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[int, int, int, int]]:
    image = Image.open(image_path).convert("L").transpose(Image.Transpose.ROTATE_270)
    pixels = np.asarray(image)
    if pixels.shape != (2481, 3509):
        raise ValueError(f"expected rotated 300 dpi page 2481x3509; got {pixels.shape}")
    binary = pixels < 130
    x_left, x_right, y_top, y_zero = _detect_axes(binary)

    # The printed plot spans -2..15, with units of 100 s.  Retain 0..600 s.
    x_start = int(round(x_left + 2.0 / 17.0 * (x_right - x_left)))
    x_stop = int(round(x_left + 8.0 / 17.0 * (x_right - x_left)))
    x_pixels = np.arange(x_start, x_stop + 1, dtype=int)
    times = -200.0 + (x_pixels - x_left) / (x_right - x_left) * 1700.0

    # Top-axis ticks extend into the graph.  The upper 8% is kept out of the
    # search; none of these selected traces approaches the printed full scale.
    search_top = y_top + int(round(0.08 * (y_zero - y_top)))
    path_y = np.full(x_pixels.size, y_zero, dtype=float)
    for column, x_pixel in enumerate(x_pixels):
        rows = np.flatnonzero(binary[search_top:y_zero - 4, x_pixel]) + search_top
        if not rows.size:
            continue
        # Join scan gaps up to three pixels, select the ink group nearest the
        # zero line, then take its upper edge.  This retains real fast rises
        # while rejecting isolated paper/scan specks above the trace.
        groups = np.split(rows, np.flatnonzero(np.diff(rows) > 3) + 1)
        selected = max(groups, key=lambda group: (float(np.mean(group)), group.size))
        path_y[column] = float(np.min(selected))

    values = ordinate_max * (y_zero - path_y) / (y_zero - y_top)
    values = np.clip(values, 0.0, ordinate_max)
    return times, values, x_pixels, path_y, (x_left, x_right, y_top, y_zero)


def _save_overlay(
    image_path: Path,
    output_path: Path,
    x_pixels: np.ndarray,
    y_pixels: np.ndarray,
) -> None:
    image = Image.open(image_path).convert("RGB").transpose(Image.Transpose.ROTATE_270)
    draw = ImageDraw.Draw(image)
    draw.line(
        [(int(x), int(round(y))) for x, y in zip(x_pixels, y_pixels)],
        fill=(220, 0, 0),
        width=4,
    )
    image.save(output_path)


def _trace_metrics(time_s: np.ndarray, values: np.ndarray, threshold: float) -> dict[str, float | None]:
    dose = float(np.trapezoid(values, time_s))
    peak_index = int(np.argmax(values))
    hits = np.flatnonzero(values >= threshold)
    return {
        "dose_vol_percent_s": dose,
        "peak_vol_percent": float(values[peak_index]),
        "peak_time_s": float(time_s[peak_index]),
        "threshold_vol_percent": float(threshold),
        "threshold_arrival_s": None if not hits.size else float(time_s[hits[0]]),
        "threshold_departure_s": None if not hits.size else float(time_s[hits[-1]]),
        "threshold_duration_s": None
        if not hits.size
        else float(time_s[hits[-1]] - time_s[hits[0]]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--image-dir", type=Path, default=ROOT / "tmp/pdfs/ti021")
    parser.add_argument("--workbook", type=Path, default=ROOT / "ti21_input_ascii.xls")
    parser.add_argument(
        "--baseline",
        type=Path,
        default=ROOT / "outputs/ti021-baseline-after-align-default-2026-10-16.npz",
    )
    parser.add_argument(
        "--pdf", type=Path, default=ROOT / "reference/admlc-trial-021.pdf"
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=ROOT / "outputs/ti021-gas-traces-digitized-2026-10-17",
    )
    args = parser.parse_args()

    dataset = read_thorney_ti021_workbook(args.workbook)
    baseline = np.load(args.baseline)
    model_time = np.asarray(baseline["time_s"], dtype=float)
    model_concentration = np.asarray(
        baseline["concentration_vol_percent"], dtype=float
    )
    common_time = np.arange(0.0, 601.0, 1.0)

    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    observed_traces = []
    model_traces = []
    rows = []
    for page in TRACE_PAGES:
        image_path = args.image_dir / f"gas-{page.pdf_page:03d}.png"
        raster_t, raster_c, path_x, path_y, axes = _digitize_trace(
            image_path, page.ordinate_max_vol_percent
        )
        observed = np.interp(common_time, raster_t, raster_c)
        sensor_index = _find_sensor_index(dataset, page)
        sensor = dataset.dose_sensors[sensor_index]
        model = np.interp(common_time, model_time, model_concentration[sensor_index])
        # Three raster pixels or 0.5% of full-scale, whichever is larger.
        threshold = max(
            0.005 * page.ordinate_max_vol_percent,
            3.0 * page.ordinate_max_vol_percent / (axes[3] - axes[2]),
        )
        observed_metric = _trace_metrics(common_time, observed, threshold)
        model_metric = _trace_metrics(common_time, model, threshold)
        tabulated_dose = float(sensor.dose_percent_s)
        digitized_dose = float(observed_metric["dose_vol_percent_s"])
        model_dose = float(model_metric["dose_vol_percent_s"])
        overlay_path = prefix.with_name(prefix.name + f"-{page.gas_id}-overlay.png")
        _save_overlay(image_path, overlay_path, path_x, path_y)
        rows.append({
            **asdict(page),
            "workbook_sensor_index_one_based": sensor_index + 1,
            "workbook_arrival_s": sensor.arrival_time_s,
            "workbook_departure_s": sensor.departure_time_s,
            "workbook_dose_vol_percent_s": tabulated_dose,
            "digitized": observed_metric,
            "model": model_metric,
            "digitized_dose_error_vs_workbook_percent": 100.0 * (
                digitized_dose / tabulated_dose - 1.0
            ),
            "model_dose_error_vs_workbook_percent": 100.0 * (
                model_dose / tabulated_dose - 1.0
            ),
            "model_to_digitized_peak_ratio": float(
                model_metric["peak_vol_percent"] / observed_metric["peak_vol_percent"]
            ),
            "detected_axes_pixels": {
                "x_left_minus_200s": axes[0],
                "x_right_1500s": axes[1],
                "y_top": axes[2],
                "y_zero": axes[3],
            },
            "source_image": str(image_path.resolve()),
            "source_image_sha256": _sha256(image_path),
            "overlay": str(overlay_path),
        })
        observed_traces.append(observed)
        model_traces.append(model)

    digitized_abs = np.abs([
        row["digitized_dose_error_vs_workbook_percent"] for row in rows
    ])
    model_abs = np.abs([row["model_dose_error_vs_workbook_percent"] for row in rows])
    report = {
        "schema": "degali.ti021-digitized-gas-traces.v1",
        "status": "diagnostic evidence; no gas-trace amplitude calibration applied",
        "source_pdf": str(args.pdf.resolve()),
        "source_pdf_sha256": _sha256(args.pdf),
        "workbook": str(args.workbook.resolve()),
        "workbook_sha256": dataset.workbook_sha256,
        "baseline": str(args.baseline.resolve()),
        "baseline_sha256": _sha256(args.baseline),
        "digitization_contract": {
            "render": "pdftoppm -r 300 -png; rotate 90 degrees clockwise",
            "time_axis": "printed -2..15 means -200..1500 s",
            "retained_time_s": [0.0, 600.0],
            "ordinate": "printed sensor-specific concentration in vol.%",
            "trace_scaling_to_workbook_dose": False,
            "reported_sensor_averaging_time_s": 0.6,
        },
        "summary": {
            "selected_sensor_count": len(rows),
            "digitized_vs_workbook_dose_median_absolute_percent": float(
                np.median(digitized_abs)
            ),
            "digitized_vs_workbook_dose_maximum_absolute_percent": float(
                np.max(digitized_abs)
            ),
            "baseline_vs_workbook_selected_mape_percent": float(np.mean(model_abs)),
        },
        "sensor_rows": rows,
    }
    np.savez_compressed(
        prefix.with_suffix(".npz"),
        time_s=common_time,
        observed_concentration_vol_percent=np.asarray(observed_traces),
        model_concentration_vol_percent=np.asarray(model_traces),
        workbook_sensor_index_one_based=np.asarray([
            row["workbook_sensor_index_one_based"] for row in rows
        ]),
    )
    prefix.with_suffix(".json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
