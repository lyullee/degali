"""Retain and resolve flagged near-knot tangent samples by finer differences."""
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

from stage2_matched_source import ROOT,REF,read,digest,write_new
from run_preslhy_ambient_profile_audit import replay
from audit_actual_enthalpy_tangents import caloric_tangent,profile_q
from degali.addons.axisymmetric_jet import _air_phase_property_table
from stage1 import verify


def main():
    output=REF/'actual_enthalpy_tangent_refinement_2026-09-06.json'
    if output.exists(): raise FileExistsError(output)
    paths=[REF/'actual_enthalpy_tangents_2026-09-06.json',
        REF/'phase_ambient_consistency_field_complete_2026-09-05.json',REF/'e35_reduced.json',
        ROOT/'tools/audit_actual_enthalpy_tangents.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    prior,field=read(paths[0]),read(paths[1]);trials={t['trial']:t for t in read(paths[2])['trials']}
    for path,expected in prior['sha256'].items():
        if digest(ROOT/path)!=expected: raise ValueError('prior tangent inputs changed')
    verify()
    flagged=[r for r in prior['rows'] if r['status']!='resolved']
    result=[]
    for old in flagged:
        n=old['trial'];trajectory,check=replay(field,trials[n])
        state=np.array(field['interfaces'][str(n)]['downstream']['states'][old['index']])
        rho,c,y,t,h=(float(v[0]) for v in profile_q(trajectory.model,state,old['q']))
        grid=_air_phase_property_table()['temperature']
        distance=float(min(abs(grid-t)))
        q=old['q'];yq=-c*trajectory.model.rhoa/rho**2
        values=[]
        for scale in (1.,.5,.1,.05,.02,.01):
            ht,hy=caloric_tangent(trajectory.model.thermodynamics,t,y,scale)
            dq=1e-4*scale
            _,_,_,tp,hp=profile_q(trajectory.model,state,np.array([q-dq,q+dq]))
            tq,hq=(float((a[1]-a[0])/(2*dq)) for a in (tp,hp))
            norm=max(abs(hq),abs(ht*tq),abs(hy*yq),1.)
            values.append(dict(scale=scale,temperature_step_K=.002*scale,q_step=dq,
                h_T_J_kg_K=ht,h_Y_J_kg=hy,T_q=tq,h_q=hq,chain_scaled_error=abs(hq-ht*tq-hy*yq)/norm,
                composition_gradient_fraction=hy*yq/hq,thermal_gradient_fraction=ht*tq/hq,
                central_T_difference_crosses_nearest_table_knot=.002*scale>distance,
                q_stencil_T_values_K=tp.tolist()))
        a,b=values[-2:]
        changes={k:abs(a[k]-b[k])/max(abs(b[k]),1.) for k in ('h_T_J_kg_K','h_Y_J_kg','T_q','h_q')}
        passed=max(changes.values())<=1e-3 and b['chain_scaled_error']<=1e-3 and b['h_T_J_kg_K']>0
        result.append(dict(trial=n,x_m=old['x_m'],q=q,old_status=old['status'],T_K=t,
            nearest_phase_temperature_knot_K=float(grid[np.argmin(abs(grid-t))]),distance_to_knot_K=distance,
            original_replay=check,evaluations=values,last_two_scaled_changes=changes,resolved=passed,
            note='Original coarse failure retained. A near table knot is not a demonstrated physical singularity.'))
    verify()
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('refinement input changed')
    write_new(output,dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        original_rows_preserved=True,rows=result,all_flagged_resolved=all(r['resolved'] for r in result),
        no_new_physical_parameters=True,new_dispersion_field=False,candidate_promoted=False))
    print(result,flush=True)


if __name__=='__main__': main()
