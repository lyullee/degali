"""Conditional new-shape modal rates; separate numerical and physical gates."""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from degali.addons.enriched_transport import (PrescribedRadialMixing,EnrichedModalTransport,
    shifted_advective_moments)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    parser.add_argument("--trials",type=int,nargs="+")
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+".partial.json")
    if args.output.exists() or partial.exists():
        parser.error("refusing to overwrite modal-transport evidence")
    ref=ROOT/"reference/preslhy"
    previous=ref/"edge_conservative_refit_combined_2026-09-05.json"
    upstream=json.loads(previous.read_text(encoding="utf-8"))
    if not upstream["all_passed"]:
        raise ValueError("all frozen enriched sections must pass before transport reconstruction")
    for name,digest in upstream["sha256"].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f"frozen input/source changed: {name}")
    paths=[ROOT/n for n in upstream["sha256"]]
    paths += [previous,Path(__file__).resolve(),ROOT/"src/degali/addons/enriched_transport.py",
              ROOT/"tests/test_enriched_transport.py",ROOT/"docs/prereg-enriched-shape-transport.md"]
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    selected=args.trials or [10,11,12,22,23,24,25]
    if len(set(selected))!=len(selected) or not set(selected).issubset({r["trial"] for r in upstream["rows"]}):
        parser.error("select unique frozen trials")
    boundary=json.loads((ref/"buoyancy_constrained_enthalpy_width_interface_2026-09-05.json").read_text(encoding="utf-8"))
    reduced=json.loads((ref/"e35_reduced.json").read_text(encoding="utf-8"))["trials"]
    measured=json.loads((ref/"measured_pipe_source_2026-09-05.json").read_text(encoding="utf-8"))["trials"]
    measured={r["trial"]:r for r in measured}
    rows,caches=[],{}
    def payload(completed):
        return dict(phase="gauge_free_conditional_enriched_modal_transport",completed=completed,
            selected_trials=selected,rows=rows,sha256=initial_hashes,
            all_numerics_passed=len(rows)==len(selected) and all(r.get("numerics_passed",False) for r in rows),
            all_constitutive_passed=len(rows)==len(selected) and all(r.get("constitutive_passed",False) for r in rows),
            mixing_input="prescribed_original_positive_phase_resolved_chi_C",thermal_species_ratio=1.,
            downstream_integrated=False,field_scored=False,promoted=False,
            limitations=["Local rates conditional on supplied scalar mixing and reduced immediate-shear-heating assumptions.",
                "Radial conservative mass/momentum flux ansatz does not include arbitrary divergence-free transverse flow.",
                "Weak boundary enforcement is not pointwise constitutive agreement; inferred negative shear diffusivity rejects adoption."])
    def checkpoint():
        partial.write_text(json.dumps(payload(False),default=serial,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    checkpoint()
    for n in selected:
        started=time.perf_counter()
        row=dict(trial=n,numerics_passed=False,constitutive_passed=False)
        rows.append(row)
        print(f"trial {n}: reconstruct enriched mass/stress and coupled modal rates",flush=True)
        try:
            old=next(r for r in upstream["rows"] if r["trial"]==n)
            p,row["source_reconstruction"]=construct_projection(old,boundary,reduced,measured,caches)
            parameters=np.array(old["parameters"])
            outputs,models,mixings=[],[],[]
            for mixing_order,order,angular_order in ((8,4,8),(16,8,16)):
                supplied=PrescribedRadialMixing.from_weak_baseline(p,order=mixing_order)
                model=EnrichedModalTransport(p,parameters,scalar_mixing=supplied,
                    thermal_species_ratio=1.,mechanical_work="reduced_buoyancy_work_immediate_shear_heat")
                out=model.assemble(order=order,angular_order=angular_order)
                outputs.append(out)
                models.append(model)
                mixings.append(supplied)
                row["evaluations"]=outputs
                checkpoint()
                print(f"trial {n}: {order}/{angular_order}, linear={out['linear_scaled_error']:.3e}, "
                      f"chi_P min={out['minimum_chi_momentum']:.4g}, edge={out['edge_defects']}",flush=True)
            coarse,fine=outputs
            rates=fine["rates"]
            rate_change=float(max(abs(rates-coarse["rates"])/np.maximum(abs(rates),1.)))
            matrix_change=float(np.max(abs(fine["matrix"]-coarse["matrix"])/np.maximum(abs(fine["matrix"]),1.)))
            # Interior quarter points in every original phase panel, not the
            # interpolation nodes, verify the supplied mixing representation.
            knots=mixings[-1].knots
            q=(knots[:-1,None]+np.diff(knots)[:,None]*np.array([.17,.43,.79])).ravel()
            v8,v16=[s.evaluate(q)[:,0] for s in mixings]
            mixing_change=float(max(abs(v8-v16)/np.maximum(abs(v16),1e-12)))
            scale=np.maximum(abs(fine["fluxes"]),1.)
            t=np.linspace(0.,1.,33)
            aa,bb=np.meshgrid(t,t,indexing="ij")
            prep=p.prepare(aa.ravel(),bb.ravel())
            psi=prep["psi"]
            h=2e-5/max(1.,max(abs(rates)))
            for par,der in ((parameters[:p.basis.size],rates[7:7+p.basis.size]),(parameters[p.basis.size:],rates[7+p.basis.size:])):
                current,speed=psi@par,abs(psi@der)
                mask=speed>1e-12
                h=min(h,float(np.min(.15*(.1-abs(current[mask]))/speed[mask])))
            if not h>0.:
                raise ValueError("no two-sided finite-change step inside the frozen shape trust region")
            derivatives=[]
            row.update(rate_refinement=rate_change,matrix_refinement=matrix_change,mixing_refinement=mixing_change,fd_step=h)
            checkpoint()
            print(f"trial {n}: independent actual-flux finite changes, step={h:.3e} m",flush=True)
            for ds in (h,h/2):
                plus=shifted_advective_moments(p,parameters,rates,ds)
                minus=shifted_advective_moments(p,parameters,rates,-ds)
                derivatives.append((plus-minus)/(2*ds))
            expected=fine["moment_jacobian"]@rates
            fd_error=float(max(abs(derivatives[-1]-expected)/scale))
            fd_change=float(max(abs(derivatives[-1]-derivatives[0])/scale))
            numerical=bool(fine["linear_scaled_error"]<=1e-8 and fine["weak_heat_scaled_error"]<=1e-8
                and max(rate_change,matrix_change,mixing_change,fd_error,fd_change)<=1e-5
                and abs(fine["mass_boundary_error"])/max(abs(fine["source"][0]),1.)<=1e-8)
            row.update(fd_derivatives=derivatives,expected_flux_rate=expected,fd_scaled_error=fd_error,
                fd_step_change=fd_change,numerics_passed=numerical,
                constitutive_passed=bool(numerical and fine["constitutive_passed"]))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row["failure"]=str(exc)
        row["elapsed_seconds"]=time.perf_counter()-started
        checkpoint()
        print(f"trial {n}: numerics={row['numerics_passed']}, constitutive={row['constitutive_passed']}, "
              f"refinements={[row.get(k) for k in ['rate_refinement','matrix_refinement','mixing_refinement','fd_scaled_error','fd_step_change']]}, "
              f"failure={row.get('failure')}, seconds={row['elapsed_seconds']:.1f}",flush=True)
    if hashes()!=initial_hashes:
        raise RuntimeError("code/input changed during the transport audit")
    args.output.write_text(json.dumps(payload(True),default=serial,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
