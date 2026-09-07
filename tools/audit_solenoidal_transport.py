"""Conditional conservative redistribution on two fixed trial10 fields."""

import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK
from degali.addons.enriched_transport import PrescribedRadialMixing
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.solenoidal_transport import SolenoidalTransport,redistributed_fluxes


def independent_vectors(model,out):
    p,m=model.projection,model.projection.mixing
    angles=np.linspace(0.,math.pi/4,65)
    defects=np.zeros(3)
    minimum_prod,minimum_diff,maximum_mass,nonradial=np.inf,np.inf,-np.inf,0.
    signatures=[]
    for angle,cuts in zip(angles,model.quadrature.partitions(model.parameters,angles)):
        extra=model.mixing.knots[(model.mixing.knots>0.)&(model.mixing.knots<cuts[-1])]
        cuts=np.sort(np.r_[cuts,extra]); cuts=cuts[np.r_[True,np.diff(cuts)>cuts[-1]*2e-13]]
        d,fm,fp,fe,_=model._ray(angle,cuts,16)
        vectors=redistributed_fluxes(model,d,fm,fp,out['rates'],out['mass_modes'],out['momentum_modes'])
        t=math.tan(angle)
        ed=model.local(np.array([1.]),np.array([t]))
        edge=redistributed_fluxes(model,ed,fe[None,0,:],fe[None,1,:],out['rates'],out['mass_modes'],out['momentum_modes'])
        em,ep=edge['mass'][0,0],edge['momentum'][0,0]
        u,rho,y,h=[ed[k][0] for k in ('u','rho','y','specific_h')]
        q=np.array([p.q0*(1+t*t)])
        chi,oldfm=model.mixing.evaluate(q)[0]
        old=m.local(q)
        target=(.5*m.wind**2+.5*u*u)*em-u*ep
        residual=np.array([y*em-model.area*rho*chi*ed['grad_y'][0,0]/(2*p.q0),
            h*em-model.area*rho*chi*ed['grad_h'][0,0]/(2*p.q0)-target,ep-m.wind*math.cos(m.state[3])*em])
        scales=np.array([max(abs(old['y'][0]*oldfm),1e-12),max(abs(old['h'][0]/old['rho'][0]*oldfm),1.),max(abs(old['u'][0]*oldfm),1.)])
        defects=np.maximum(defects,abs(residual)/scales)
        for v in (vectors,edge):
            minimum_prod=min(minimum_prod,float(min(v['production'])))
            minimum_diff=min(minimum_diff,float(min(v['radial_momentum_diffusivity'])))
            nonradial=max(nonradial,float(max(v['nonradial_stress_fraction'])))
        maximum_mass=max(maximum_mass,float(em))
        signatures.append(np.column_stack([vectors['mass'],vectors['momentum']]))
    return dict(edge_defects=dict(zip(('hydrogen','heat','momentum'),defects)),
        minimum_shear_production=minimum_prod,minimum_radial_momentum_diffusivity=minimum_diff,
        maximum_nonradial_stress_fraction=nonradial,maximum_outward_mass=maximum_mass,
        signature=np.vstack(signatures),angles=65,radial_order=16)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite solenoidal evidence')
    ref=ROOT/'reference/preslhy'
    origin=ref/'enriched_boundary_exit_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed independent failed endpoint required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen upstream changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/solenoidal_transport.py',ROOT/'tests/test_solenoidal_transport.py',
        ROOT/'docs/prereg-solenoidal-transport.md']
    def hashes():
        return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial_hashes=hashes()
    def read(name):
        return json.loads((ref/name).read_text(encoding='utf-8'))
    seed=next(r for r in read('edge_conservative_refit_combined_2026-09-05.json')['rows'] if r['trial']==10)
    initial=next(r for r in read('coupled_initialization_combined_2026-09-06.json')['rows'] if r['trial']==10)
    measured={r['trial']:r for r in read('measured_pipe_source_2026-09-05.json')['trials']}
    p,replay=construct_projection(seed,read('buoyancy_constrained_enthalpy_width_interface_2026-09-05.json'),
        read('e35_reduced.json')['trials'],measured,{})
    par=np.array(initial['parameters'])
    supplied=PrescribedRadialMixing.from_weak_baseline(p,order=32)
    driver=EnrichedShortSegment(p,par,scalar_mixing=supplied,thermal_species_ratio=1.,mechanical_work=WORK,
        mixing_update='equilibrium_velocity_width')
    failed=np.array(next(r for r in upstream['rows'] if r['divisions']==8)['endpoint'])
    data=dict(phase='conditional_solenoidal_flux_pilot',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    for label,encoded in (('initial',driver.initial),('rejected_8step_endpoint',failed)):
        row=dict(trial=10,state=label,levels=[],necessary_screen_passed=False)
        data['rows'].append(row); checkpoint()
        try:
            current,parameters,mix,factor=driver.context(encoded)
            row['mixing_amplitude_ratio']=factor
            for k in (4,6):
                level=dict(flux_modes=k,evaluations=[],passed=False)
                row['levels'].append(level); checkpoint()
                for order,angular in ((8,48),(16,96)):
                    model=SolenoidalTransport(current,parameters,scalar_mixing=mix,thermal_species_ratio=1.,
                        mechanical_work=WORK,flux_modes=k,order=order,angular_order=angular)
                    try:
                        out=model.assemble()
                    except ValueError as exc:
                        level['linear_failure']=getattr(model,'linear_diagnostics',{})
                        raise ValueError(f'K{k} order{order}/{angular}: {exc}') from exc
                    level['evaluations'].append(out); checkpoint()
                    print('solenoidal',label,'K',k,'order',order,angular,'edges',out['edge_defects'],
                        'production',out['minimum_shear_production'],'nonradial',out['maximum_nonradial_stress_fraction'],
                        'numerics',out['local_numerics_passed'],flush=True)
                coarse,fine=level['evaluations']
                level['independent']=independent_vectors(model.model,fine)
                rays=level['independent']
                level['rate_refinement']=float(max(abs(fine['solution']-coarse['solution'])/np.maximum(abs(fine['solution']),1.)))
                level['matrix_refinement']=float(np.max(abs(fine['matrix']-coarse['matrix'])/np.maximum(abs(fine['matrix']),1.)))
                level['passed']=bool(all(o['local_numerics_passed'] for o in (coarse,fine))
                    and max(level['rate_refinement'],level['matrix_refinement'])<=1e-5
                    and max(*fine['edge_defects'].values(),*rays['edge_defects'].values())<=.05
                    and min(fine['minimum_shear_production'],rays['minimum_shear_production'])>=0.
                    and min(fine['minimum_radial_momentum_diffusivity'],rays['minimum_radial_momentum_diffusivity'])>=0.
                    and max(fine['maximum_outward_mass'],rays['maximum_outward_mass'])<0.
                    and fine['curvature_half_width']<.1)
                checkpoint()
                if not level['passed']:
                    break
            if len(row['levels'])==2 and all(l['passed'] for l in row['levels']):
                a,b=[l['evaluations'][-1]['rates'] for l in row['levels']]
                row['flow_mode_rate_change']=float(max(abs(a-b)/np.maximum(abs(b),1.)))
                a,b=[l['independent']['signature'] for l in row['levels']]
                row['flow_mode_vector_change']=float(np.max(abs(a-b)/np.maximum(abs(b),1.)))
                row['necessary_screen_passed']=bool(max(row['flow_mode_rate_change'],row['flow_mode_vector_change'])<=1e-5)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row['failure']=str(exc)
        row['completed']=True; checkpoint()
        print('solenoidal',label,'completed',row['necessary_screen_passed'],row.get('failure'),flush=True)
    if hashes()!=initial_hashes:
        raise RuntimeError('solenoidal dependencies changed during screen')
    data['completed']=True
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
