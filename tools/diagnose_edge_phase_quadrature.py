"""Post-screen independent phase-cell audit; original coefficients unchanged."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from audit_transverse_mixing import ROOT, actual_source_model
from audit_reservoir_thermal_segments import serial
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.reservoir_thermal import ReservoirThermalMoments
from degali.addons.edge_enrichment import FixedTransportEdgeProjection
from degali.addons.edge_phase_quadrature import PhaseSplitSquareMoments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite final or partial phase diagnosis")
    ref = ROOT/"reference/preslhy"
    upstream_path = ref/"edge_enrichment_screen_2026-09-05.json"
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    if not upstream["completed"]:
        raise RuntimeError("enrichment screen is not complete")
    for relative, digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"screen source/input changed: {relative}")
    paths = [ROOT/p for p in upstream["sha256"]]
    paths += [upstream_path, Path(__file__).resolve(), ROOT/"src/degali/addons/edge_phase_quadrature.py",
              ROOT/"tests/test_edge_phase_quadrature.py", ROOT/"docs/edge-enrichment-postscreen-numerics.md"]
    def hashes():
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes, rows, caches = hashes(), [], {}
    def payload(completed):
        return dict(phase="postscreen_phase_split_square_moment_diagnosis", completed=completed,
            rows=rows, sha256=initial_hashes, coefficients_refitted=False,
            original_screen_all_passed=upstream["all_passed"], promoted=False,
            limitations=["Selected lower independent edge-defect candidate per trial for numerical diagnosis only.",
                "No original failure is overwritten; no new shape transport or field scoring."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    boundary = json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced = json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured = json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured_rows = {r["trial"]: r for r in measured}
    checkpoint()
    for old in upstream["rows"]:
        start = time.perf_counter()
        n = old["trial"]
        eligible = [f for f in old["fits"] if "enriched_edge_defects" in f]
        selected = min(eligible, key=lambda f: max(f["enriched_edge_defects"].values()))
        row = dict(trial=n, degree=selected["degree"], parameters=selected["parameters"],
                   selected_edge_defects=selected["enriched_edge_defects"], checks=[])
        rows.append(row)
        print(f"trial {n}, degree {selected['degree']}: phase-cell quadrature diagnosis", flush=True)
        try:
            trial = next(t for t in reduced if t["trial"] == n)
            jp, th, _ = actual_source_model(trial, measured_rows)
            entry = boundary["interfaces"][str(n)]
            section = BuoyancyConstrainedEnthalpySection(jp, th, thermal_width_ratio=entry["thermal_width_ratio"],
                                                        quadrature_points=entry["quadrature_points"])
            section._quadrature_cache = caches.setdefault(section.k.delta, {})
            mixing = ConservativeTransverseMixing(section, np.array(entry["state"]))
            reservoir = ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work="reduced_buoyancy_work")
            baseline = reservoir.evaluate(order=16, probes=1025)
            projection = FixedTransportEdgeProjection(reservoir, baseline, degree=selected["degree"])
            if max(abs(projection.target-selected["target"])/projection.scales) > 1e-12:
                raise RuntimeError("original projection target does not replay")
            quad = PhaseSplitSquareMoments(projection)
            zero = quad.moments(np.zeros(projection.count), order=16, angular_order=128)
            row.update(target=projection.target, scales=projection.scales, zero_moments=zero,
                zero_scaled_error=float(max(abs(zero-projection.target)/projection.scales)))
            checkpoint()
            for radial, angular in ((8, 32), (16, 64), (16, 128)):
                actual = quad.moments(np.array(selected["parameters"]), order=radial, angular_order=angular)
                residual = (actual-projection.target)/projection.scales
                check = dict(radial_order=radial, angular_order=angular, moments=actual,
                             scaled_residual=residual, maximum_scaled_error=float(max(abs(residual))))
                row["checks"].append(check)
                checkpoint()
                print(f"trial {n}: {radial}/{angular} moments error={check['maximum_scaled_error']:.3e}", flush=True)
            values = [c["moments"] for c in row["checks"]]
            row["adjacent_changes"] = [float(max(abs(values[i+1]-values[i])/projection.scales)) for i in (0, 1)]
            tensor = np.array(selected["actual_moments"][-1])
            row["tensor_384_to_phase_split_change"] = float(max(abs(tensor-values[-1])/projection.scales))
            row["independent_moment_passed"] = bool(row["checks"][-1]["maximum_scaled_error"] <= 1e-8)
            row["diagnostic_quadrature_passed"] = bool(row["zero_scaled_error"] <= 1e-8 and row["adjacent_changes"][-1] <= 1e-5)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["elapsed_seconds"] = time.perf_counter()-start
        checkpoint()
        print(f"trial {n}: zero={row.get('zero_scaled_error')}, refined={row.get('adjacent_changes')}, "
              f"failure={row.get('failure')}, elapsed={row['elapsed_seconds']:.1f}s", flush=True)
    if hashes() != initial_hashes:
        raise RuntimeError("code/input changed during the phase diagnosis")
    args.output.write_text(json.dumps(payload(True), default=serial, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
