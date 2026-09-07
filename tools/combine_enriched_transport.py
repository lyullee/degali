"""Retain original gates, separately combine independently verified fine inputs."""

import argparse,hashlib,json
from pathlib import Path
from audit_transverse_mixing import ROOT


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    parser.add_argument("mixing_verification",type=Path)
    parser.add_argument("groups",type=Path,nargs="+")
    args=parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite combined transport evidence")
    records=[json.loads(p.read_text(encoding="utf-8")) for p in [args.mixing_verification,*args.groups]]
    sha={}
    for record in records:
        if not record["completed"]:
            raise ValueError("only completed records can be combined")
        for key,digest in record["sha256"].items():
            if hashlib.sha256((ROOT/key).read_bytes()).hexdigest()!=digest or (key in sha and sha[key]!=digest):
                raise ValueError(f"source version mismatch: {key}")
            sha[key]=digest
    mix=records[0]
    if mix["phase"]!="independent_prescribed_mixing_input_verification":
        raise ValueError("wrong independent input verification")
    rows=[]
    for group in records[1:]:
        if group["phase"]!="gauge_free_conditional_enriched_modal_transport":
            raise ValueError("wrong transport group")
        if [r["trial"] for r in group["rows"]]!=group["selected_trials"]:
            raise ValueError("incomplete group coverage")
        for old in group["rows"]:
            row=old.copy()
            verification=next(r for r in mix["rows"] if r["trial"]==row["trial"])
            row["original_numerics_passed"]=row["numerics_passed"]
            row["fine_input_verification"]=verification
            checked=False
            if "failure" not in row and len(row.get("evaluations",[]))==2:
                fine=row["evaluations"][-1]
                checked=bool(verification["passed"] and fine["linear_scaled_error"]<=1e-8
                    and fine["weak_heat_scaled_error"]<=1e-8
                    and abs(fine["mass_boundary_error"])/max(abs(fine["source"][0]),1.)<=1e-8
                    and max(row[k] for k in ("rate_refinement","matrix_refinement","fd_scaled_error","fd_step_change"))<=1e-5)
            row["independently_verified_numerics"]=checked
            row["combined_constitutive_passed"]=bool(checked and row["evaluations"][-1]["constitutive_passed"])
            rows.append(row)
    if sorted(r["trial"] for r in rows)!=[10,11,12,22,23,24,25]:
        raise ValueError("all seven trials are required exactly once")
    rows.sort(key=lambda r:r["trial"])
    sources=[args.mixing_verification,*args.groups,Path(__file__).resolve()]
    for p in sources:
        sha[str(p.resolve().relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    payload=dict(phase="combined_conditional_enriched_modal_transport",completed=True,rows=rows,
        sha256=sha,source_records=[str(p.resolve().relative_to(ROOT)) for p in sources[:-1]],
        original_numerics_passed_trials=[r["trial"] for r in rows if r["original_numerics_passed"]],
        independently_verified_numerics_trials=[r["trial"] for r in rows if r["independently_verified_numerics"]],
        constitutive_passed_trials=[r["trial"] for r in rows if r["combined_constitutive_passed"]],
        original_failures_preserved=True,downstream_integrated=False,field_scored=False,promoted=False)
    args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print({k:v for k,v in payload.items() if k.endswith("_trials")})


if __name__=="__main__":
    main()
