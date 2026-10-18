"""Digitize the independent TI021 wind plots without using gas receptors.

The HSE scan contains a 0.6 s wind-speed plot on PDF page 178 and a wind-
heading plot on page 181.  The plots also contain a three-minute running mean.
This tool follows that smooth trace with a bounded-slope dynamic programme,
then preserves only its time variation: speed is re-centred to the frozen
3.4 m/s workbook mean and heading to the report-priority -6.1 degree mean.

The re-centring is deliberate.  The remote plotted vane and the report summary
do not have the same absolute mean/reference, so treating the raster ordinate
as an exact absolute heading would silently choose one conflicting source.
No gas concentration, dose, arrival, or sensor location enters this tool.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _trace_smooth_curve(
    image_path: Path,
    *,
    axis: tuple[int, int, int, int],
    ordinate: tuple[float, float],
    search_range: tuple[float, float],
    time_range_s: tuple[float, float] = (0.0, 600.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return time/value and path pixels from one 300 dpi rotated scan.

    The declared graph abscissa is -300..1500 s.  Only 0..600 s is retained.
    A four-pixel maximum inter-column jump excludes the jagged raw trace while
    following the separately plotted three-minute running mean.
    """

    original = Image.open(image_path).convert("L")
    rotated = original.transpose(Image.Transpose.ROTATE_270)
    pixels = np.asarray(rotated, dtype=float)
    if pixels.shape != (2481, 3509):
        raise ValueError(
            f"expected a 300 dpi 2481x3509 rotated page, got {pixels.shape}"
        )

    x0, x1, y0, y1 = axis
    value_top, value_bottom = ordinate
    x_start = int(round(x0 + (float(time_range_s[0]) + 300.0) / 1800.0 * (x1 - x0)))
    x_stop = int(round(x0 + (float(time_range_s[1]) + 300.0) / 1800.0 * (x1 - x0)))
    x_pixels = np.arange(x_start, x_stop + 1, dtype=int)

    value_low, value_high = search_range
    y_start = int(round(
        y0 + (value_top - value_high) / (value_top - value_bottom) * (y1 - y0)
    ))
    y_stop = int(round(
        y0 + (value_top - value_low) / (value_top - value_bottom) * (y1 - y0)
    ))
    y_pixels = np.arange(y_start, y_stop + 1, dtype=int)
    raster = pixels[np.ix_(y_pixels, x_pixels)]

    # Black ink is cheap; white background is expensive.  A small slope cost
    # chooses the continuous running-mean trace when raw spikes intersect it.
    observation_cost = np.clip((raster - 40.0) / 180.0, 0.0, 1.0)
    height, width = observation_cost.shape
    cost = observation_cost[:, 0].copy()
    back = np.zeros((height, width), dtype=np.int16)
    maximum_jump = 4
    slope_cost = 0.09
    indices = np.arange(height)
    for column in range(1, width):
        next_cost = np.empty(height, dtype=float)
        for row in range(height):
            lower = max(0, row - maximum_jump)
            upper = min(height, row + maximum_jump + 1)
            candidates = (
                cost[lower:upper]
                + slope_cost * np.abs(indices[lower:upper] - row)
            )
            selected = int(np.argmin(candidates))
            back[row, column] = lower + selected
            next_cost[row] = observation_cost[row, column] + candidates[selected]
        cost = next_cost

    path = np.empty(width, dtype=int)
    path[-1] = int(np.argmin(cost))
    for column in range(width - 1, 0, -1):
        path[column - 1] = back[path[column], column]
    path_y = y_pixels[path]

    time_s = -300.0 + (x_pixels - x0) / (x1 - x0) * 1800.0
    values = value_top - (path_y - y0) / (y1 - y0) * (
        value_top - value_bottom
    )
    return time_s, values, x_pixels, path_y


def _trace_raw_curve(
    image_path: Path,
    *,
    axis: tuple[int, int, int, int],
    ordinate: tuple[float, float],
    search_range: tuple[float, float],
    time_range_s: tuple[float, float] = (0.0, 600.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Recover the raster-resolvable 0.6 s trace, independent of gas data.

    At 300 dpi the graph supplies only about 1.15 horizontal pixels per
    second, so the original 0.6 s samples cannot be recovered one-for-one.
    We therefore report a one-column median of the jagged trace.  Pixels near
    the separately traced three-minute running mean are removed whenever an
    alternative ink path exists.  Vertical connecting strokes are represented
    by their midpoint rather than selecting an artificial extreme.
    """

    smooth_t, _, x_pixels, smooth_y = _trace_smooth_curve(
        image_path, axis=axis, ordinate=ordinate, search_range=search_range,
        time_range_s=time_range_s,
    )
    original = Image.open(image_path).convert("L")
    rotated = original.transpose(Image.Transpose.ROTATE_270)
    pixels = np.asarray(rotated, dtype=float)
    x0, x1, y0, y1 = axis
    value_top, value_bottom = ordinate
    value_low, value_high = search_range
    y_start = int(round(
        y0 + (value_top - value_high) / (value_top - value_bottom) * (y1 - y0)
    ))
    y_stop = int(round(
        y0 + (value_top - value_low) / (value_top - value_bottom) * (y1 - y0)
    ))
    y_pixels = np.arange(y_start, y_stop + 1, dtype=int)
    raster = pixels[np.ix_(y_pixels, x_pixels)]
    raw_y = np.full(x_pixels.size, np.nan, dtype=float)
    for column in range(x_pixels.size):
        candidates = y_pixels[raster[:, column] < 105.0]
        if candidates.size == 0:
            continue
        separated = candidates[np.abs(candidates - smooth_y[column]) >= 4]
        selected = separated if separated.size else candidates
        raw_y[column] = float(np.median(selected))
    valid = np.isfinite(raw_y)
    if np.count_nonzero(valid) < 0.8 * raw_y.size:
        raise RuntimeError("raw meteorology trace is too sparse for reconstruction")
    raw_y = np.interp(np.arange(raw_y.size), np.flatnonzero(valid), raw_y[valid])
    values = value_top - (raw_y - y0) / (y1 - y0) * (
        value_top - value_bottom
    )
    return smooth_t, values, x_pixels, np.rint(raw_y).astype(int)


def _save_overlay(
    image_path: Path,
    output_path: Path,
    x_pixels: np.ndarray,
    y_pixels: np.ndarray,
) -> None:
    image = Image.open(image_path).convert("RGB").transpose(Image.Transpose.ROTATE_270)
    draw = ImageDraw.Draw(image)
    points = [(int(x), int(y)) for x, y in zip(x_pixels, y_pixels)]
    draw.line(points, fill=(220, 0, 0), width=4)
    image.save(output_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, default=ROOT / "reference/admlc-trial-021.pdf")
    parser.add_argument("--speed-image", type=Path, required=True)
    parser.add_argument("--heading-image", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=ROOT / "outputs/ti021-meteorology-digitized-2026-10-17",
    )
    parser.add_argument("--target-speed-m-s", type=float, default=3.4)
    parser.add_argument("--target-report-heading-deg", type=float, default=-6.1)
    parser.add_argument("--heading-fluctuation-scale", type=float, default=1.0)
    parser.add_argument(
        "--trace-mode", choices=("running_mean", "raw_column_median"),
        default="running_mean",
    )
    parser.add_argument("--time-start-s", type=float, default=0.0)
    parser.add_argument("--time-stop-s", type=float, default=600.0)
    parser.add_argument(
        "--output-time-shift-s", type=float, default=0.0,
        help="add this shift to the stored time coordinate after digitization",
    )
    parser.add_argument(
        "--pad-stored-time-stop-s", type=float, default=None,
        help="optionally extend the final sample as a constant boundary value",
    )
    args = parser.parse_args()
    if not -300.0 <= args.time_start_s < args.time_stop_s <= 1500.0:
        raise ValueError("trace time range must satisfy -300 <= start < stop <= 1500")
    time_range = (args.time_start_s, args.time_stop_s)

    tracer = (
        _trace_raw_curve if args.trace_mode == "raw_column_median"
        else _trace_smooth_curve
    )
    speed_t, speed_raw, speed_x, speed_y = tracer(
        args.speed_image,
        axis=(863, 2935, 345, 1689),
        ordinate=(5.0, 0.0),
        search_range=(2.0, 4.7),
        time_range_s=time_range,
    )
    heading_t, heading_raw, heading_x, heading_y = tracer(
        args.heading_image,
        axis=(839, 2911, 347, 1689),
        ordinate=(50.0, -50.0),
        search_range=(-35.0, 5.0),
        time_range_s=time_range,
    )

    digitized_time_s = np.arange(args.time_start_s, args.time_stop_s + 0.5, 1.0)
    speed_raw = np.interp(digitized_time_s, speed_t, speed_raw)
    heading_raw = np.interp(digitized_time_s, heading_t, heading_raw)
    time_s = digitized_time_s + float(args.output_time_shift_s)
    speed = speed_raw * float(args.target_speed_m_s) / float(np.mean(speed_raw))
    report_heading = (
        float(args.heading_fluctuation_scale)
        * (heading_raw - float(np.mean(heading_raw)))
        + float(args.target_report_heading_deg)
    )
    solver_heading = -report_heading
    wind_xy = np.column_stack((
        speed * np.cos(np.deg2rad(solver_heading)),
        speed * np.sin(np.deg2rad(solver_heading)),
    ))
    if (args.pad_stored_time_stop_s is not None
            and float(args.pad_stored_time_stop_s) > time_s[-1]):
        time_s = np.r_[time_s, float(args.pad_stored_time_stop_s)]
        speed = np.r_[speed, speed[-1]]
        report_heading = np.r_[report_heading, report_heading[-1]]
        solver_heading = np.r_[solver_heading, solver_heading[-1]]
        wind_xy = np.vstack((wind_xy, wind_xy[-1]))
        speed_raw = np.r_[speed_raw, speed_raw[-1]]
        heading_raw = np.r_[heading_raw, heading_raw[-1]]

    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    npz_path = prefix.with_suffix(".npz")
    json_path = prefix.with_suffix(".json")
    speed_overlay = prefix.with_name(prefix.name + "-speed-overlay.png")
    heading_overlay = prefix.with_name(prefix.name + "-heading-overlay.png")
    np.savez_compressed(
        npz_path,
        time_s=time_s,
        wind_speed_m_s=speed,
        report_heading_relative_deg=report_heading,
        solver_heading_deg=solver_heading,
        wind_xy_m_s=wind_xy,
        raw_digitized_speed_m_s=speed_raw,
        raw_digitized_heading_deg=heading_raw,
    )
    _save_overlay(args.speed_image, speed_overlay, speed_x, speed_y)
    _save_overlay(args.heading_image, heading_overlay, heading_x, heading_y)

    report = {
        "schema": "degali.ti021-digitized-meteorology.v1",
        "status": "independent-input sensitivity; raster-derived",
        "source_pdf": str(args.pdf.resolve()),
        "source_pdf_sha256": _sha256(args.pdf),
        "source_pages": {
            "wind_speed_pdf_page": 178,
            "wind_heading_pdf_page": 181,
            "report_summary_pdf_page": 3,
            "history_pdf_page": 4,
            "sensor_and_turbulence_notes_pdf_page": 9,
        },
        "render_contract": "pdftoppm -r 300 -png, then rotate 90 degrees clockwise",
        "speed_image_sha256": _sha256(args.speed_image),
        "heading_image_sha256": _sha256(args.heading_image),
        "digitization": {
            "method": (
                "one-column median of jagged 0.6 s trace after running-mean exclusion"
                if args.trace_mode == "raw_column_median"
                else "bounded-slope minimum-cost trace of plotted three-minute running mean"
            ),
            "trace_mode": args.trace_mode,
            "time_resolution_s": 1.0,
            "declared_plot_averaging_time_s": 0.6,
            "estimated_speed_digitization_uncertainty_m_s": 0.05,
            "estimated_heading_digitization_uncertainty_deg": 1.0,
            "gas_receptor_values_used": False,
            "digitized_time_range_s": [
                float(args.time_start_s), float(args.time_stop_s)
            ],
            "stored_time_shift_s": float(args.output_time_shift_s),
            "constant_end_padding_to_s": (
                None if args.pad_stored_time_stop_s is None
                else float(args.pad_stored_time_stop_s)
            ),
            "heading_fluctuation_scale": float(args.heading_fluctuation_scale),
        },
        "independent_report_facts": {
            "main_data_mean_wind_speed_m_s": 3.9,
            "report_priority_mean_heading_relative_to_array_deg": -6.1,
            "stability_at_release": "D/E",
            "turbulence_intensity_percent": {
                "longitudinal": 15.7,
                "crosswind": 10.1,
                "vertical": 8.0,
            },
            "gas_sensor_error_band": "about +/-10% of reading; per-sensor estimate stated as separately issued",
        },
        "conflict_policy": {
            "speed_target_mean_m_s": float(args.target_speed_m_s),
            "heading_target_mean_deg": float(args.target_report_heading_deg),
            "reason": (
                "retain raster time variation while preserving the frozen workbook speed "
                "and report-priority heading; do not silently select the remote-vane absolute mean"
            ),
        },
        "raw_digitized_summary": {
            "speed_mean_m_s": float(np.mean(speed_raw)),
            "speed_min_m_s": float(np.min(speed_raw)),
            "speed_max_m_s": float(np.max(speed_raw)),
            "heading_mean_deg": float(np.mean(heading_raw)),
            "heading_min_deg": float(np.min(heading_raw)),
            "heading_max_deg": float(np.max(heading_raw)),
        },
        "normalized_summary": {
            "speed_mean_m_s": float(np.mean(speed)),
            "heading_mean_deg": float(np.mean(report_heading)),
            "solver_heading_min_deg": float(np.min(solver_heading)),
            "solver_heading_max_deg": float(np.max(solver_heading)),
        },
        "npz": str(npz_path),
        "speed_overlay": str(speed_overlay),
        "heading_overlay": str(heading_overlay),
        "promotion_allowed": False,
    }
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
