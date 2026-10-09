"""Decision-oriented LH2 release cases for the Applied Energy manuscript.

The cases are reference calculations, not site-specific operating limits.  A
screening-line result answers only whether a modelled 4 vol % contour reaches
a declared line under the frozen inputs.  It does not create a legal setback,
an operating permit, or a probabilistic risk estimate.

Run from the repository root with ``PYTHONPATH=src``.  The output directory
must be empty.  No third-party raw data are copied into the derived outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics as statlib
import time
from pathlib import Path

import numpy as np

from degali.lh2 import assess
from degali.validation import nearfield


ROOT = Path(__file__).resolve().parents[1]
BOUNDARIES_M = (20.0, 30.0, 50.0)
LFL_VOL_PCT = 4.0
SLABX_MANIFEST = (
    ROOT / "outputs/applied-energy-evidence-2026-10-03/"
    "test4-external/test4_external_manifest.json"
)

JET_CASES = (
    {
        "case_id": "FFI_T4_HIGH_WIND_TRANSFER",
        "description": "horizontal transfer release; higher-wind FFI Test 4",
        "rate_kg_s": 0.828,
        "orifice_m": 0.0254,
        "height_m": 0.5,
        "winds_m_s": (5.0, 5.85, 6.7),
        "wind_labels": ("reported_low", "midpoint", "reported_high"),
        "ambient_temperature_k": 276.45,
        "relative_humidity_pct": 75.0,
        "storage_pressure_barg": 2.12,
        "roughness_m": 0.03,
        "averaging_s": 275.0,
        "wind_reference_height_m": 10.0,
        "receptor_height_m": 1.5,
        "evidence": "FFI/DNV Test 4; conditional far-field screen",
    },
    {
        "case_id": "FFI_T6_LOW_WIND_TRANSFER",
        "description": "horizontal transfer release; low-wind FFI Test 6",
        "rate_kg_s": 0.833,
        "orifice_m": 0.0254,
        "height_m": 0.5,
        "winds_m_s": (2.3, 2.5, 2.7),
        "wind_labels": ("reported_low", "midpoint", "reported_high"),
        "ambient_temperature_k": 277.15,
        "relative_humidity_pct": 90.0,
        "storage_pressure_barg": 2.53,
        "roughness_m": 0.001,
        "averaging_s": 60.0,
        "wind_reference_height_m": 10.0,
        "receptor_height_m": 1.5,
        "evidence": "FFI/DNV Test 6; unresolved 30 m low-wind residual",
    },
)

# The four published NASA pool-spill states already used by the package's
# independent regime test.  The labels are deliberately neutral because the
# transferred source records do not support reconstructing a site operation.
POOL_CASES = (
    ("NASA_POOL_1", 9.23, 1.55, 297.15, 49.0, 18.3, "aloft"),
    ("NASA_POOL_2", 10.29, 3.35, 288.15, 43.0, 6.4, "low"),
    ("NASA_POOL_3", 9.95, 6.30, 285.15, 43.0, 0.3, "grounded"),
    ("NASA_POOL_4", 9.48, 2.20, 288.15, 29.0, 3.4, "low"),
)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_jet(case: dict, wind: float, corrected: bool):
    jet, initial = nearfield.hydrogen_jet(
        rate=case["rate_kg_s"],
        diameter=case["orifice_m"],
        wind=wind,
        height=case["height_m"],
        ambient_temperature=case["ambient_temperature_k"],
        relative_humidity=case["relative_humidity_pct"],
        storage_pressure_barg=case["storage_pressure_barg"],
        roughness=case["roughness_m"],
        stability="D",
        averaging=case["averaging_s"],
        wind_reference_height=case["wind_reference_height_m"],
        corrections=corrected,
        ground_effect=False,
    )
    run = jet.run(initial, distmx=nearfield.STEP, smax=250.0)
    trajectory = nearfield.Trajectory(jet.th.table, run.rows)
    if not trajectory.ok:
        raise RuntimeError(f"{case['case_id']} integration failed at wind={wind}")
    return trajectory


def concentration_at(trajectory, x_m: float, height_m: float) -> float:
    if trajectory.at(x_m) is None:
        return float("nan")
    return float(trajectory.concentration_at(x_m, 0.0, height_m))


def lfl_reach(trajectory, *, height_m: float, limit_m: float = 250.0):
    """Furthest fixed-height exit from the 4 vol % contour."""
    xs = np.geomspace(0.05, limit_m, 1600)
    points = [
        (float(x), concentration_at(trajectory, float(x), height_m))
        for x in xs
        if trajectory.at(float(x)) is not None
    ]
    above = [i for i, (_, value) in enumerate(points) if value >= LFL_VOL_PCT]
    if not above:
        return float("nan"), "never_reaches_lfl"
    last = above[-1]
    if last == len(points) - 1:
        return float("nan"), "right_censored"
    left_x, left_c = points[last]
    right_x, right_c = points[last + 1]
    crossing = left_x + (LFL_VOL_PCT - left_c) * (
        right_x - left_x
    ) / (right_c - left_c)
    return float(crossing), "resolved"


def interpolation_at(trajectory: np.ndarray, distance_m: float) -> float:
    if trajectory.ndim != 2 or trajectory.shape[0] < 2:
        return float("nan")
    rows = trajectory[np.argsort(trajectory[:, 0])]
    if distance_m < rows[0, 0] or distance_m > rows[-1, 0]:
        return float("nan")
    return float(np.interp(distance_m, rows[:, 0], rows[:, 2]) * 100.0)


def external_benchmark_rows() -> list[dict]:
    """Literature/model lanes and the exact reason each can or cannot rank."""
    return [
        {
            "benchmark": "FFI/DNV Test 4 SLABx",
            "model": "SLABx-LH2",
            "same_campaign_as_local_case": "yes",
            "quantitative_values_available": "yes",
            "ranking_allowed": "no",
            "reason": "same sensors and bulk inputs, but native source mapping and temporal operator differ",
            "doi_or_url": "local accepted 640/64 run; see test4_external_manifest.json",
        },
        {
            "benchmark": "FFI/DNV Test 6 external model",
            "model": "PHAST/EFFECTS/CFD",
            "same_campaign_as_local_case": "not located",
            "quantitative_values_available": "no",
            "ranking_allowed": "no",
            "reason": "no licensed same-input run or published matched receptor table is present",
            "doi_or_url": "pending external run",
        },
        {
            "benchmark": "HSL Buxton Tests 6 and 7",
            "model": "FLACS",
            "same_campaign_as_local_case": "no",
            "quantitative_values_available": "figures/time-series only",
            "ranking_allowed": "no",
            "reason": "different 2010 HSL campaign: 0.071 kg/s and 26.3 mm; not FFI 2019 Test 6 at 0.833 kg/s",
            "doi_or_url": "https://doi.org/10.1016/j.ijhydene.2012.05.145",
        },
        {
            "benchmark": "HSL Buxton Tests 5, 6 and 7",
            "model": "ADREA-HF",
            "same_campaign_as_local_case": "no",
            "quantitative_values_available": "published figures",
            "ranking_allowed": "no",
            "reason": "different campaign and source; useful for mechanism discussion, not paired score",
            "doi_or_url": "https://doi.org/10.1016/j.ijhydene.2014.07.042",
        },
        {
            "benchmark": "NASA WSTF large pool releases",
            "model": "FLACS",
            "same_campaign_as_local_case": "related pool evidence",
            "quantitative_values_available": "hazard-distance trends; no imported paired table",
            "ranking_allowed": "no",
            "reason": "supports trend and model-form discussion; source/observation operators are not paired here",
            "doi_or_url": "https://doi.org/10.1016/j.ijhydene.2020.06.131",
        },
    ]


def external_metric_rows(path: Path = SLABX_MANIFEST) -> list[dict]:
    """Read only derived external scores; never copy the third-party table."""
    if not path.exists():
        return [{
            "case_id": "FFI_T4_HIGH_WIND_TRANSFER",
            "model": "SLABx-LH2",
            "n_sensors": "",
            "mae_vol_pct_point": "",
            "rmse_vol_pct_point": "",
            "lfl_true_positive": "",
            "lfl_false_negative": "",
            "lfl_false_positive": "",
            "ranking_allowed": "no",
            "status": "derived external manifest absent",
        }]
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for model in ("SLABx", "DEGALI"):
        score = payload["scores"][model]
        rows.append({
            "case_id": "FFI_T4_HIGH_WIND_TRANSFER",
            "model": model,
            "n_sensors": score["n_sensors"],
            "mae_vol_pct_point": score["mae_vol_pct_point"],
            "rmse_vol_pct_point": score["rmse_vol_pct_point"],
            "lfl_true_positive": score["true_positive"],
            "lfl_false_negative": score["false_negative"],
            "lfl_false_positive": score["false_positive"],
            "ranking_allowed": "no",
            "status": "conditional: source mapping and temporal operator differ",
        })
    return rows


def make_svg(path: Path, jet_rows: list[dict], pool_rows: list[dict]) -> None:
    width, height = 1080, 650
    left, top, plot_w, plot_h = 220, 70, 780, 470
    rows = []
    for item in jet_rows:
        rows.append((
            f"{item['case_id'].replace('FFI_', '')} {item['wind_label']} {item['model']}",
            item["lfl_distance_m"],
            "#0069aa" if item["model"] == "corrected" else "#777777",
        ))
    for item in pool_rows:
        rows.append((item["case_id"] + " DEGALI", item["lfl_distance_m"], "#2c8b57"))
    rows = [row for row in rows if math.isfinite(float(row[1]))]
    xmax = max(60.0, math.ceil(max(float(row[1]) for row in rows) / 10.0) * 10.0)
    row_h = plot_h / max(len(rows), 1)

    def sx(value: float) -> float:
        return left + min(max(value, 0.0), xmax) / xmax * plot_w

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.label{font-size:12px}.small{font-size:11px}</style>',
        '<text x="70" y="35" class="title">Reference LH2 cases: modelled 4 vol % reach</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#222"/>',
    ]
    for boundary, colour in ((20, "#e6a700"), (30, "#d55e00"), (50, "#aa3344")):
        parts.append(
            f'<line x1="{sx(boundary):.1f}" y1="{top}" x2="{sx(boundary):.1f}" y2="{top + plot_h}" stroke="{colour}" stroke-dasharray="5 4"/>'
        )
        parts.append(
            f'<text x="{sx(boundary) + 3:.1f}" y="{top - 8}" class="small">{boundary} m line</text>'
        )
    for index, (label, value, colour) in enumerate(rows):
        y = top + (index + 0.5) * row_h
        parts.append(f'<text x="10" y="{y + 4:.1f}" class="label">{label}</text>')
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{sx(float(value)):.1f}" y2="{y:.1f}" stroke="{colour}" stroke-width="5"/>')
        parts.append(f'<circle cx="{sx(float(value)):.1f}" cy="{y:.1f}" r="5" fill="{colour}"/>')
        parts.append(f'<text x="{sx(float(value)) + 8:.1f}" y="{y + 4:.1f}" class="label">{float(value):.1f} m</text>')
    for tick in np.linspace(0, xmax, 7):
        parts.append(f'<text x="{sx(float(tick)) - 8:.1f}" y="{top + plot_h + 24}" class="label">{tick:.0f}</text>')
    parts += [
        f'<text x="{left + plot_w / 2 - 90}" y="{height - 55}" class="label">Downwind/fixed-height LFL reach (m)</text>',
        '<text x="70" y="625" class="small">Jet rows: 1.5 m receptor height. Pool rows: plume-axis reach. These are conditional screening diagnostics, not setbacks.</text>',
        '</svg>',
    ]
    path.write_text("\n".join(parts), encoding="utf-8")


def build_markdown(jet_rows: list[dict], decision_rows: list[dict], pool_rows: list[dict], external_rows: list[dict]) -> str:
    changed = [row for row in decision_rows if row["historical_vs_corrected_change"] == "yes"]
    test6 = [row for row in jet_rows if row["case_id"] == "FFI_T6_LOW_WIND_TRANSFER" and row["model"] == "corrected"]
    pool_errors = [abs(row["lowest_flammable_height_error_m"]) for row in pool_rows]
    lines = [
        "# Decision-oriented LH2 scenario evidence (2026-10-03)",
        "",
        "## Scope",
        "",
        "These are reference calculations that connect a dispersion result to a declared",
        "screening line. They are not real-site operating data, statutory separation",
        "distances, probabilistic risk estimates or permission to transfer LH2.",
        "",
        "## Main result",
        "",
        f"Across the horizontal-jet wind cases and the 20/30/50 m lines, **{len(changed)}**",
        "screening classifications change between the historical reconstruction and the",
        "corrected path. This is the model-choice consequence relevant to the manuscript.",
        "",
        "For FFI Test 6, the corrected 1.5 m LFL reach over the reported wind interval is",
        f"**{min(row['lfl_distance_m'] for row in test6):.1f}–{max(row['lfl_distance_m'] for row in test6):.1f} m**.",
        "The unresolved 30 m concentration residual means that this interval remains a",
        "conditional diagnostic, not an operational exclusion zone.",
        "",
        "## Reference cases",
        "",
        "- FFI Test 4: horizontal, higher-wind transfer analogue.",
        "- FFI Test 6: horizontal, low-wind transfer analogue and present failure boundary.",
        "- Four NASA pool-spill states: storage/spill archetype and buoyancy-regime check.",
        "",
        "The NASA pool states reproduce the published regime in all four cases; the mean",
        f"absolute error in lowest flammable height is **{statlib.mean(pool_errors):.2f} m**.",
        "Pool-axis LFL distance is reported separately from the jet's fixed 1.5 m receptor",
        "operator and must not be pooled into one validation score.",
        "",
        "## External-model gate",
        "",
        "The local SLABx Test 4 calculation is retained as a conditional fixed-sensor",
        "comparison. No same-input PHAST, EFFECTS or CFD table is available for FFI Test 6.",
        "The published HSL Tests 6/7 CFD studies concern a different 2010 campaign",
        "(0.071 kg/s), not the 2019 FFI Test 6 (0.833 kg/s), and are therefore used only",
        "for mechanism discussion. Renumbering coincidence is not treated as a match.",
        "",
        "## Files",
        "",
        "- `jet_scenario_summary.csv`: model, wind, runtime and fixed-height LFL reach.",
        "- `screening_line_decisions.csv`: concentration and category at 20/30/50 m.",
        "- `pool_reference_scenarios.csv`: pool-axis distance and regime evidence.",
        "- `external_benchmark_register.csv`: comparison eligibility and provenance.",
        "- `external_model_metrics.csv`: derived SLABx/DEGALI Test 4 sensor scores.",
        "- `decision_lfl_distances.svg`: manuscript-ready diagnostic figure.",
        "- `manifest.json`: frozen inputs, limitations and file hashes.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "$env:PYTHONPATH = 'src'",
        ".venv\\Scripts\\python.exe tools\\audit_decision_scenarios.py `",
        "  --output-dir tmp\\decision-scenarios-reproduction",
        "```",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"output directory not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    jet_rows: list[dict] = []
    decision_native: list[dict] = []
    for case in JET_CASES:
        for wind_label, wind in zip(case["wind_labels"], case["winds_m_s"]):
            by_model = {}
            for model, corrected in (("historical_reconstruction", False), ("corrected", True)):
                start = time.perf_counter()
                trajectory = build_jet(case, wind, corrected)
                distance, status = lfl_reach(
                    trajectory, height_m=case["receptor_height_m"]
                )
                runtime = time.perf_counter() - start
                concentrations = {
                    boundary: concentration_at(
                        trajectory, boundary, case["receptor_height_m"]
                    ) for boundary in BOUNDARIES_M
                }
                by_model[model] = concentrations
                jet_rows.append({
                    "case_id": case["case_id"],
                    "description": case["description"],
                    "model": model,
                    "wind_label": wind_label,
                    "wind_m_s": wind,
                    "rate_kg_s": case["rate_kg_s"],
                    "orifice_m": case["orifice_m"],
                    "release_height_m": case["height_m"],
                    "receptor_height_m": case["receptor_height_m"],
                    "lfl_distance_m": distance,
                    "lfl_status": status,
                    "runtime_s": runtime,
                    "evidence_scope": case["evidence"],
                })
            for boundary in BOUNDARIES_M:
                old = by_model["historical_reconstruction"][boundary]
                new = by_model["corrected"][boundary]
                old_exceeds = bool(math.isfinite(old) and old >= LFL_VOL_PCT)
                new_exceeds = bool(math.isfinite(new) and new >= LFL_VOL_PCT)
                decision_native.append({
                    "case_id": case["case_id"],
                    "wind_label": wind_label,
                    "wind_m_s": wind,
                    "screening_line_m": boundary,
                    "receptor_height_m": case["receptor_height_m"],
                    "historical_vol_pct": old,
                    "corrected_vol_pct": new,
                    "historical_line_exceeded": "yes" if old_exceeds else "no",
                    "corrected_line_exceeded": "yes" if new_exceeds else "no",
                    "historical_vs_corrected_change": "yes" if old_exceeds != new_exceeds else "no",
                    "interpretation": "conditional model-choice diagnostic; not an operating rule",
                })

    pool_rows: list[dict] = []
    for case_id, rate, wind, ambient, rh, measured_height, measured_regime in POOL_CASES:
        start = time.perf_counter()
        result = assess(
            rate=rate,
            wind=wind,
            pool_diameter=9.1,
            ambient_temperature=ambient,
            relative_humidity=rh,
            max_distance=100.0,
            at_distance=33.8,
        )
        pool_rows.append({
            "case_id": case_id,
            "model": "DEGALI_LiftoffPlume",
            "rate_kg_s": rate,
            "wind_m_s": wind,
            "pool_diameter_m": 9.1,
            "ambient_temperature_k": ambient,
            "relative_humidity_pct": rh,
            "lfl_distance_m": result.distance_to_lfl,
            "stoichiometric_distance_m": result.distance_to_stoichiometric,
            "measured_lowest_flammable_height_m": measured_height,
            "predicted_lowest_flammable_height_m": result.lowest_flammable_height,
            "lowest_flammable_height_error_m": result.lowest_flammable_height - measured_height,
            "measured_regime": measured_regime,
            "predicted_regime": result.regime,
            "regime_match": "yes" if result.regime == measured_regime else "no",
            "screening_scope": result.screening_scope,
            "runtime_s": time.perf_counter() - start,
            "warnings": " | ".join(result.warnings),
        })

    external_rows = external_benchmark_rows()
    external_metrics = external_metric_rows()
    write_csv(out / "jet_scenario_summary.csv", jet_rows)
    write_csv(out / "screening_line_decisions.csv", decision_native)
    write_csv(out / "pool_reference_scenarios.csv", pool_rows)
    write_csv(out / "external_benchmark_register.csv", external_rows)
    write_csv(out / "external_model_metrics.csv", external_metrics)
    make_svg(out / "decision_lfl_distances.svg", jet_rows, pool_rows)
    (out / "README.md").write_text(
        build_markdown(jet_rows, decision_native, pool_rows, external_rows),
        encoding="utf-8",
    )

    source_paths = [
        ROOT / "src/degali/lh2.py",
        ROOT / "src/degali/validation/nearfield.py",
        Path(__file__),
    ]
    if SLABX_MANIFEST.exists():
        source_paths.append(SLABX_MANIFEST)
    manifest = {
        "purpose": "decision-oriented Applied Energy reference cases",
        "qualification": "conditional screening diagnostics; not setbacks, permits, QRA or operating limits",
        "threshold_vol_pct": LFL_VOL_PCT,
        "screening_lines_m": BOUNDARIES_M,
        "jet_operator": "centreline crosswind coordinate, fixed receptor height 1.5 m, steady model",
        "pool_operator": "plume-axis LFL reach; lowest flammable height evaluated at 33.8 m",
        "external_model_rule": "rank only with identical source, weather, geometry, receptors and time operator",
        "runtime_environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "input_sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(out),
        "jet_runs": len(jet_rows),
        "decision_rows": len(decision_native),
        "changed_decisions": sum(
            row["historical_vs_corrected_change"] == "yes"
            for row in decision_native
        ),
        "pool_runs": len(pool_rows),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
