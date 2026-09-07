"""Value-based smooth enthalpy component: cancellation and geometry isolation."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import mpmath
from audit_transverse_mixing import ROOT
from audit_reservoir_thermal_segments import serial
from audit_edge_conservative_refit import construct_projection
from audit_coupled_shape_initialization import WORK
from degali.addons.enriched_transport import PrescribedRadialMixing,EnrichedModalTransport
from degali.addons.enriched_segments import EnrichedShortSegment
from degali.addons.stable_enthalpy_difference import stable_enthalpy_difference
from degali.addons.phase_radial_quadrature import gauss_rule


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('output',type=Path); args=parser.parse_args()
    partial=args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite smooth component evidence')
    ref=ROOT/'reference/preslhy'; origin=ref/'polished_energy_difference_2026-09-06.json'
    upstream=json.loads(origin.read_text(encoding='utf-8'))
    if not upstream['completed']:
        raise ValueError('completed root-polishing diagnosis required')
    for name,digest in upstream['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'frozen dependency changed: {name}')
    paths=[ROOT/n for n in upstream['sha256']]+[origin,Path(__file__).resolve(),
        ROOT/'src/degali/addons/stable_enthalpy_difference.py',ROOT/'tests/test_stable_enthalpy_difference.py',
        ROOT/'docs/prereg-smooth-enthalpy-conditioning.md']
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
    par=np.array(initial['parameters']); mix=PrescribedRadialMixing.from_weak_baseline(p,order=32)
    driver=EnrichedShortSegment(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK,
        mixing_update='equilibrium_velocity_width')
    failed=np.array(next(r for r in read('enriched_boundary_exit_2026-09-06.json')['rows'] if r['divisions']==8)['endpoint'])
    p,par,mix,_=driver.context(failed)
    model=EnrichedModalTransport(p,par,scalar_mixing=mix,thermal_species_ratio=1.,mechanical_work=WORK)
    witness=read('phase_consistent_flux_witness_2026-09-06.json')
    rates=np.array(witness['evaluations'][-1]['rates']); step=witness['finite_difference']['step_m']
    energy_expected=witness['finite_difference']['expected'][4]; scale=max(abs(energy_expected),1.)
    data=dict(phase='smooth_enthalpy_conditioning',completed=False,rows=[],sha256=initial_hashes,
        started_utc=datetime.now(timezone.utc).isoformat(),source_reconstruction=replay,
        original_step=step,total_energy_expected=energy_expected,mpmath_version=mpmath.__version__,
        full_energy_verified=False,adopted=False,field_scored=False)
    def checkpoint():
        partial.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    checkpoint()
    try:
        expected=[]
        for order in (64,96):
            x,w=gauss_rule(order); a,b=np.meshgrid(x,x,indexing='ij')
            local=model.local(a.ravel(),b.ravel())
            expected.append(float((8*p.q0*np.outer(w,w)).ravel()@local['bh']@rates))
        data['enthalpy_expected']=expected
        data['analytic_quadrature_refinement_energy_scaled']=abs(expected[1]-expected[0])/scale
        center_hr,center_hc=p.section._phase_partials(np.asarray(p.mixing.state[0]),np.asarray(p.mixing.state[0]*p.mixing.state[1]))
        data['center_enthalpy_expected']=float(center_hr*p.mixing.state[0]*rates[0]+center_hc*p.mixing.state[0]*p.mixing.state[1]*rates[1])
        data['wind_expected']=float(p.mixing.wind_partials[:7]@rates[:7]); checkpoint()
        for policy in ('binary_legacy','exact_constraint'):
            for factor in (.5,1.):
                for order,digits in ((32,50),(64,50),(64,70)):
                    out=stable_enthalpy_difference(p,par,rates,step*factor,order=order,precision=digits,wind_policy=policy)
                    error=out['derivative']-expected[-1]
                    out.update(factor=factor,absolute_error=error,energy_scaled_error=abs(error)/scale,
                        enthalpy_scaled_error=abs(error)/max(abs(expected[-1]),1.))
                    pair=out['scalar_pair']
                    out['center_difference_error']=pair['h']['difference']/(2*step*factor)-data['center_enthalpy_expected']
                    out['wind_difference_error']=pair['wind']['difference']/(2*step*factor)-data['wind_expected']
                    data['rows'].append(out); checkpoint()
                    print('smooth H',policy,factor,order,digits,out['energy_scaled_error'],flush=True)
        data['component_checks']={}
        for policy in ('binary_legacy','exact_constraint'):
            rows=[r for r in data['rows'] if r['scalar_pair']['wind_policy']==policy]
            values=[r['derivative'] for r in rows]
            data['component_checks'][policy]=dict(maximum_energy_scaled_error=max(r['energy_scaled_error'] for r in rows),
                overall_spread_energy_scaled=(max(values)-min(values))/scale,
                component_passed=bool(max(r['energy_scaled_error'] for r in rows)<=1e-5 and (max(values)-min(values))/scale<=1e-5))
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        data['failure']=str(exc)
    if hashes()!=initial_hashes:
        raise RuntimeError('smooth component dependencies changed')
    data['completed']=True
    args.output.write_text(json.dumps(data,default=serial,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('smooth component completed',data.get('component_checks'),data.get('failure'),flush=True)


if __name__=='__main__':
    main()
