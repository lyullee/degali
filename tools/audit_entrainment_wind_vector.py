"""Direction-only ambient momentum audit on unchanged actual field ledgers."""
from datetime import datetime,timezone
from pathlib import Path
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from stage2_matched_source import ROOT,REF,read,digest,write_new
from run_preslhy_ambient_profile_audit import replay
from degali.addons.ambient_vector_exchange import entrained_ambient_rates
from stage1 import verify


def main():
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',REF/'e35_reduced.json',
        REF/'lateral_temperature_symmetry_2026-09-06.json',REF/'stage1_single_image_2026-09-06.json',
        ROOT/'docs/prereg-entrainment-wind-vector.md',ROOT/'src/degali/addons/ambient_vector_exchange.py',
        ROOT/'tests/test_ambient_vector_exchange.py',ROOT/'tools/run_preslhy_ambient_profile_audit.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    directory=REF/'entrainment_wind_vector_2026-09-06';directory.mkdir(exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes))
    verify();field,reduced,evidence,frozen=(read(p) for p in paths[:4])
    trials={r['trial']:r for r in reduced['trials']};summary={};checks={}
    for n in (10,23):
        trajectory,checks[str(n)]=replay(field,trials[n]);model=trajectory.model
        down=field['interfaces'][str(n)]['downstream']
        states,arc,sources,flux=(np.array(down[k]) for k in ('states','arc_length','sources','fluxes'))
        m=sources[:,0];u=np.array([model._wind(s) for s in states])
        if np.any(m<0): raise ValueError('frozen field contains ambient detrainment')
        base=entrained_ambient_rates(m,np.column_stack([u,0*u,0*u]),specific_enthalpy=0.)
        if max(abs(base[:,2]-m*u))>1e-12: raise ValueError('coflow entrainment replay mismatch')
        cases=[];wind=evidence['wind'][str(n)]
        if wind['fast']['status']=='valid_report_no_fault':
            cases.append(dict(label='farfield_frozen_window_vector_mean',direction_from_deg=wind['fast']['mean_direction_from_deg'],
                direction_source=wind['fast']['range'],measurement_height_m=wind['fast']['height_m']))
        cases.extend(dict(label=f'local_row_{r["excel_row"]}',direction_from_deg=r['direction_from_deg'],
            direction_source=r['source_range'],time=r['time'],measurement_height_m=1.5)
            for r in wind['local_bracketing_records'] if not r['missing_direction'])
        record=[]
        for case in cases:
            a=math.radians(case['direction_from_deg']+180.-75.)
            velocity=np.column_stack([u*math.cos(a),u*math.sin(a),u*0])
            new=entrained_ambient_rates(m,velocity,specific_enthalpy=0.)
            change=new-base
            energy_relative=float(max(abs(change[:,5])/np.maximum(abs(base[:,5]),1.)))
            if energy_relative>1e-12: raise ValueError('direction-only kinetic energy changed')
            cumulative=cumulative_trapezoid(change,arc,axis=0,initial=0.)
            ip=np.unique(np.r_[np.arange(0,len(arc),2),len(arc)-1])
            coarse=cumulative_trapezoid(change[ip],arc[ip],axis=0,initial=0.)
            scale=np.maximum(np.hypot(flux[:,2],flux[:,3]),1.)
            refine=float(np.max(abs(coarse[:,2:5]-cumulative[ip,2:5])/scale[ip,None]))
            transverse=cumulative[:,3]/scale
            stream=cumulative[:,2]/scale
            sensors=[]
            for x in sorted({r['x'] for r in frozen['temperature_rows'] if r['trial']==n}):
                k=int(np.argmin(abs(states[:,5]-x)))
                sensors.append(dict(requested_x_m=x,stored_x_m=float(states[k,5]),
                    cumulative_delta_Pxyz_kg_m_s2=cumulative[k,2:5],
                    transverse_over_old_momentum=float(transverse[k]),streamwise_change_over_old_momentum=float(stream[k])))
            case_summary=dict(**case,maximum_abs_transverse_over_old_momentum=float(max(abs(transverse))),
                maximum_abs_streamwise_change_over_old_momentum=float(max(abs(stream))),
                maximum_scaled_coarse_integral_difference=refine,maximum_relative_energy_change=energy_relative,
                sensor_sections=sensors)
            record.append(dict(**case_summary,arc_length_m=arc,x_m=states[:,5],mass_rate_kg_s_m=m,
                frozen_local_wind_speed_m_s=u,cumulative_delta_flux=cumulative))
            summary[f'{n}_{case["label"]}']=case_summary
        write_new(directory/f'trial_{n}.json',dict(trial=n,cases=record,original_entrainment_geometry_and_drag_unchanged=True))
        print('Wind vector forcing',n,[{k:v for k,v in r.items() if k.startswith('maximum') or k=='label'} for r in record],flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('input changed')
    verify()
    write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        summary=summary,replay_checks=checks,new_field=False,three_dimensional_solver_completed=False,
        source_formation_recomputed=False,candidate_promoted=False,direction_only_diagnostic=True))


if __name__=='__main__': main()
