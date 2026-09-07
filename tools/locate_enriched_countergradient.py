"""Locate higher-order interior witnesses for rejected positive-viscosity cases."""

import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_edge_conservative_refit import construct_projection
from degali.addons.enriched_transport import EnrichedModalTransport,PrescribedRadialMixing


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("combined",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite counter-gradient witnesses")
    data=json.loads(args.combined.read_text(encoding="utf-8"))
    for k,v in data["sha256"].items():
        if hashlib.sha256((ROOT/k).read_bytes()).hexdigest()!=v:
            raise ValueError(f"changed input/source: {k}")
    paths=[ROOT/k for k in data["sha256"]]+[args.combined.resolve(),Path(__file__).resolve()]
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial=hashes()
    ref=ROOT/"reference/preslhy"
    boundary=json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced=json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured=json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured={r["trial"]:r for r in measured}
    upstream=json.loads((ref/"edge_conservative_refit_combined_2026-09-05.json").read_text(encoding="utf-8"))
    rows,caches=[],{}
    for row in data["rows"]:
        fine=row["evaluations"][-1]
        if fine["minimum_chi_momentum"]>=0.:
            continue
        old=next(r for r in upstream["rows"] if r["trial"]==row["trial"])
        p,_=construct_projection(old,boundary,reduced,measured,caches)
        supplied=PrescribedRadialMixing.from_weak_baseline(p,order=32)
        model=EnrichedModalTransport(p,np.array(old["parameters"]),scalar_mixing=supplied,
            thermal_species_ratio=1.,mechanical_work="reduced_buoyancy_work_immediate_shear_heat")
        at=np.r_[1.,fine["rates"]]
        angles=np.linspace(0.,math.pi/4,9)
        partitions=model.quadrature.partitions(model.parameters,angles)
        best=None
        for angle,k in zip(angles,partitions):
            extra=supplied.knots[(supplied.knots>0.)&(supplied.knots<k[-1])]
            k=np.sort(np.r_[k,extra])
            k=k[np.r_[True,np.diff(k)>k[-1]*2e-13]]
            d,fm,fp,edge,_=model._ray(angle,k,16)
            rp=(fp-d["u"][:,None]*fm)@at
            chi=-rp/(model.area*d["rho"]*d["uq"])
            index=int(np.argmin(chi))
            witness=dict(trial=row["trial"],angle_rad=float(angle),q=float(d["q"][index]),
                a=float(d["a"][index]),b=float(d["b"][index]),density=float(d["rho"][index]),
                temperature_K=float(d["temperature"][index]),hydrogen_mass_fraction=float(d["y"][index]),
                velocity_m_s=float(d["u"][index]),chi_momentum_per_s=float(chi[index]),
                reduced_production_W_m3=float(-2*d["q"][index]*d["uq"][index]*rp[index]/model.area))
            if best is None or witness["chi_momentum_per_s"]<best["chi_momentum_per_s"]:
                best=witness
        rows.append(best)
        print(best,flush=True)
    if hashes()!=initial:
        raise ValueError("code/input changed during witness replay")
    args.output.write_text(json.dumps(dict(phase="interior_countergradient_witnesses",completed=True,
        rows=rows,ray_order=16,angles=9,sha256=initial,
        interpretation="Violates the selected positive-eddy-viscosity/immediate-heating model, not every possible turbulence model.",
        field_scored=False,promoted=False),indent=2,allow_nan=False)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
