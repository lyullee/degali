"""Independent 16/32 and direct-formula checks after retained pilot failure."""

import argparse,hashlib,json
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from degali.addons.enriched_transport import PrescribedRadialMixing


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite independent mixing verification")
    ref=ROOT/"reference/preslhy"
    pilot_path=ref/"enriched_transport_pilot_2026-09-05.json"
    pilot=json.loads(pilot_path.read_text(encoding="utf-8"))
    if not pilot["completed"]:
        raise ValueError("pilot must be complete")
    for name,digest in pilot["sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f"pilot source/input changed: {name}")
    paths=[ROOT/n for n in pilot["sha256"]]+[pilot_path,Path(__file__).resolve(),ROOT/"docs/enriched-mixing-input-verification-note.md"]
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial=hashes()
    upstream=json.loads((ref/"edge_conservative_refit_combined_2026-09-05.json").read_text(encoding="utf-8"))
    boundary=json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced=json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured=json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured={r["trial"]:r for r in measured}
    rows,caches=[],{}
    for old in upstream["rows"]:
        p,_=construct_projection(old,boundary,reduced,measured,caches)
        supplied=[PrescribedRadialMixing.from_weak_baseline(p,order=k) for k in (16,32)]
        knots=supplied[0].knots
        q=(knots[:-1,None]+np.diff(knots)[:,None]*np.array([.08,.31,.67,.94])).ravel()
        data=p.mixing.mixing_family(q,p.baseline["family"],order=16)
        direct=data["chi"]@np.array([1.,p.baseline["gamma"]])
        v16,v32=[v.evaluate(q)[:,0] for v in supplied]
        scale=np.maximum(abs(direct),1e-12)
        errors=[float(max(abs(v-direct)/scale)) for v in (v16,v32)]
        change=float(max(abs(v16-v32)/scale))
        row=dict(trial=old["trial"],samples=len(q),direct_errors=errors,relative_16_32_change=change,
            passed=bool(max(errors+[change])<=1e-5))
        rows.append(row)
        print(row,flush=True)
    if hashes()!=initial:
        raise RuntimeError("code/input changed during mixing verification")
    args.output.write_text(json.dumps(dict(phase="independent_prescribed_mixing_input_verification",
        completed=True,rows=rows,all_passed=all(r["passed"] for r in rows),sha256=initial,
        changed_physics=False,original_failures_preserved=True),default=serial,indent=2,allow_nan=False)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
