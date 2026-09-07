"""Conditional measured-binary transport scales on unchanged seven plume fields."""
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

from stage2_matched_source import ROOT,REF,TRIALS,read,digest,write_new
from run_preslhy_ambient_profile_audit import replay
from degali.addons.binary_molecular_transport import (
    TEMPERATURE_K,DP_M2_S_ATM,DP_PRINTED_90_INTERVAL,
    constant_reference_spreading,equimolar_reference_diffusivity,
    binary_gas_species_flux)
from stage1 import verify


def main():
    directory=REF/'molecular_transport_scale_2026-09-06'
    paths=[REF/'phase_ambient_consistency_field_complete_2026-09-05.json',
        REF/'e35_reduced.json',REF/'stage1_single_image_2026-09-06.json',
        REF/'actual_enthalpy_tangents_2026-09-06.json',
        REF/'van_heijningen_1967_diffusion.pdf',ROOT/'docs/prereg-molecular-transport-scale.md',
        ROOT/'src/degali/addons/binary_molecular_transport.py',
        ROOT/'tests/test_binary_molecular_transport.py',ROOT/'tools/run_preslhy_ambient_profile_audit.py',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    directory.mkdir(exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,
        no_default_change=True,source_url='https://pure.tue.nl/ws/files/112230102/101206954.pdf',
        printed_table_page=11,pdf_table_page_one_based=20,source_warning_page=15))
    verify()
    field,reduced,frozen,tangents=(read(p) for p in paths[:4])
    trials={r['trial']:r for r in reduced['trials']}
    summary={};soret=[];checks={}
    for n in TRIALS:
        trajectory,checks[str(n)]=replay(field,trials[n]);model=trajectory.model
        down=field['interfaces'][str(n)]['downstream']
        states=np.array(down['states']);arc=np.array(down['arc_length'])
        if len(states)!=len(arc) or np.any(np.diff(arc)<=0): raise ValueError('invalid saved path')
        widths=np.array([model.section_widths(s) for s in states])
        winds=np.array([model._wind(s) for s in states])
        # Largest printed D*p plus its printed interval, not a general upper bound.
        dref=float((DP_M2_S_ATM[-1]+DP_PRINTED_90_INTERVAL[-1])*101325/model.thermodynamics.ambient_pressure)
        out=[]
        for q in (0.,.5,1.,2.,4.):
            u=winds*np.cos(states[:,3])+states[:,4]*np.exp(-model.velocity_shape_exponent*q)
            if np.any(u<=0): raise ValueError('reverse reference path velocity')
            elapsed=cumulative_trapezoid(1/u,arc,initial=0)
            coarse_idx=np.unique(np.r_[np.arange(0,len(arc),2),len(arc)-1])
            coarse=cumulative_trapezoid(1/u[coarse_idx],arc[coarse_idx],initial=0)
            denom=np.maximum(elapsed[coarse_idx],1e-12)
            refinement=float(max(abs(coarse-elapsed[coarse_idx])/denom))
            for k,(state,sigma,time) in enumerate(zip(states,widths,elapsed)):
                scale=constant_reference_spreading(time,sigma,dref)
                out.append(dict(index=k,x_m=float(state[5]),q=q,elapsed_reference_s=float(time),
                    axial_velocity_m_s=float(u[k]),sigma_yz_m=sigma,
                    reference_diffusivity_m2_s=dref,
                    diffusion_length_m=float(scale['diffusion_length'][0]),
                    relative_width_increment_yz=scale['relative_width_increment'],
                    diffusivity_for_one_percent_width_yz_m2_s=None if k==0 else scale['diffusivity_for_one_percent_width'],
                    reference_path_time_coarse_relative_difference=refinement))
        sensor_rows=[]
        for old in frozen['temperature_rows']:
            if old['trial']==n:
                k=int(np.argmin(abs(states[:,5]-old['x'])))
                sensor_rows.append(dict(channel=old['channel'],requested_x_m=old['x'],
                    nearest_stored_x_m=float(states[k,5]),reference_scales=[r for r in out if r['index']==k]))
        summary[str(n)]=dict(stored_nodes=len(states),reference_nodes=len(out),
            maximum_reference_time_s=max(r['elapsed_reference_s'] for r in out),
            maximum_reference_width_increment=max(float(max(r['relative_width_increment_yz'])) for r in out),
            maximum_reference_time_refinement=max(r['reference_path_time_coarse_relative_difference'] for r in out),
            reference_D_m2_s=dref,
            fixed_q_reference_not_actual_streamline_age=True)
        write_new(directory/f'trial_{n}.json',dict(rows=out,same_temperature_sensor_sections=sensor_rows,summary=summary[str(n)]))
        for row in tangents['rows']:
            if row['trial']!=n: continue
            record=dict(trial=n,x_m=row['x_m'],q=row['q'],original_phase=row.get('phase_classification'))
            if row.get('status')!='resolved' or row.get('phase_classification')!='all_gas':
                record.update(status='unsupported_phase_or_tangent');soret.append(record);continue
            t,y=row['T_K'],row['Y_H2'];gy,gt=row['Y_q'],row['T_q']/t
            try: d=float(equimolar_reference_diffusivity(t,model.thermodynamics.ambient_pressure))
            except ValueError:
                record.update(status='unsupported_reference_temperature');soret.append(record);continue
            f=binary_gas_species_flux(density=row['density_kg_m3'],diffusivity=d,mass_fraction=y,
                mass_fraction_gradient=gy,log_temperature_gradient=gt,thermal_diffusion_factor=1.)
            # q-gradient ratios cancel the common spatial metric, including the
            # limiting radial-curvature ratio at q=0 (both physical gradients 0).
            record.update(status='conditional_binary_unit_response',T_K=t,Y_H2=y,
                reference_D_m2_s=d,signed_Soret_over_Fick_per_alpha=float(f['soret']/f['fick']),
                reference_diffusion_is_equimolar_N2_H2_not_actual_humid_gas=True,
                no_physical_alpha_selected=True)
            soret.append(record)
        print('Molecular scale',n,summary[str(n)],flush=True)
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in paths}: raise ValueError('input changed')
    verify()
    evaluated=[r for r in soret if r['status']=='conditional_binary_unit_response']
    write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),
        sha256=hashes,summary=summary,replay_checks=checks,soret_unit_rows=soret,
        soret_summary=dict(total=len(soret),evaluated=len(evaluated),unsupported=len(soret)-len(evaluated),
            maximum_abs_Soret_over_Fick_per_alpha=max(abs(r['signed_Soret_over_Fick_per_alpha']) for r in evaluated)),
        all_reference_width_increments_below_one_percent=all(r['maximum_reference_width_increment']<.01 for r in summary.values()),
        new_field=False,candidate_promoted=False,physical_Soret_coefficient_selected=False,
        strict_molecular_upper_bound=False))


if __name__=='__main__': main()
