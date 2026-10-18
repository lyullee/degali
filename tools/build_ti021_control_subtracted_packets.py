"""Subtract the SAME obstacle carrier rather than a nominal inlet flux.

The resulting marked branch is still a one-way positive-anomaly closure;
record discarded negative mixture anomaly and vertical impulse explicitly.
This is an interface sensitivity, not a full signed binary-mixture coupling.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools.audit_thorney_ti021_handoff_streaming_preflight import build_configuration
from degali.addons.resolved_obstacle_flow_3d import solve_resolved_obstacle_flow_3d
from degali.addons.binary_mixture_handoff_plane_3d import BinaryMixtureHandoffPlane3DRequest
from degali.addons.conservative_handoff_hdf5_3d import ConservativeHandoffHdf5Writer3D


def audit_existing(donor, directory):
    """Recover the audit from complete durable histories without rerunning."""
    control_path=directory/'control.h5';output=directory/'control-subtracted.h5'
    with h5py.File(donor) as source,h5py.File(control_path) as control,h5py.File(output) as target:
        times=source['time_s'][:]
        for handle in (control,target):
            if not np.array_equal(times,handle['time_s'][:]) or times[-1]!=600.:
                raise ValueError('recovery requires all three complete identical 600s histories')
        if np.max(np.abs(control['species_mass_flux_kg_m2_s'][:]))!=0.:
            raise ValueError('air-only control contains nonzero species flux')
        if not np.array_equal(source['species_mass_flux_kg_m2_s'][:],
                              target['dense_layer_packet/species_mass_flux_kg_m2_s'][:]):
            raise ValueError('recovery found changed species packets')
        area=np.diff(source['source_z_edges_m'][:])[:,None]*np.diff(source['source_y_edges_m'][:])[None,:]
        integrate=lambda values:float(np.trapezoid(np.sum(values*area,axis=(1,2)),times))
        anomaly=source['total_mass_flux_anomaly_kg_m2_s'][:]-control['total_mass_flux_anomaly_kg_m2_s'][:]
        report={'schema':'degali.ti021-control-subtracted-packets.v1',
            'source':str(donor.resolve()),'output':str(output.resolve()),
            'donor_sha256':hashlib.sha256(donor.read_bytes()).hexdigest(),
            'control_sha256':hashlib.sha256(control_path.read_bytes()).hexdigest(),
            'original_marked_mixture_kg':integrate(source['dense_layer_packet/mixture_mass_flux_kg_m2_s'][:]),
            'control_subtracted_marked_mixture_kg':integrate(target['dense_layer_packet/mixture_mass_flux_kg_m2_s'][:]),
            'species_kg':integrate(target['dense_layer_packet/species_mass_flux_kg_m2_s'][:]),
            'discarded_negative_mixture_anomaly_kg':integrate(np.maximum(-anomaly,0.)),
            'complete_durable_histories_verified':True,
            'control_solver_diagnostics_available':False,
            'audit_recovered_after_reporting_attribute_error':True,
            'full_signed_mixture_handoff':False,'promotion_allowed':False}
    (directory/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(report,indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--donor',type=Path,required=True)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--audit-existing',action='store_true')
    args=parser.parse_args()
    record_path=args.donor.parent/'ti021-inner-transient.json'
    donor_config=json.loads(record_path.read_text(encoding='utf-8'))['configuration']
    frozen={'duration_s':600.,'time_step_s':.125,'cell_size_m':5.,
        'domain_m':[140.,140.,40.],'wind_speed_m_s':3.4,'wind_heading_relative_deg':6.1,
        'circular_arc_representation':'zero_thickness_surface','obstacle_height_m':6.,
        'stratification_damping_alpha':5.,'stratification_damping_exponent':1.,
        'stratification_shear_regularization_s_inv':.041,'finite_source_release_duration_s':None}
    for key,value in frozen.items():
        if donor_config.get(key)!=value:
            raise ValueError(f'donor {key} differs from the frozen control configuration')
    args.directory.mkdir(parents=True,exist_ok=True)
    if args.audit_existing:
        audit_existing(args.donor,args.directory)
        return
    control_path=args.directory/'control.h5'
    output=args.directory/'control-subtracted.h5'
    if control_path.exists() or output.exists():
        raise FileExistsError('choose a fresh directory; control histories are not overwritten')
    properties=json.loads((ROOT/'outputs/thorney-ti021-r12-property-audit-2026-10-11-v1.json').read_text())
    config,closure,coordinates=build_configuration(properties,cell_size_m=5.,time_step_s=.125,
        active_steps=4800,inner_domain_m=(140.,140.,40.),handoff_site_downwind_y_m=70.,
        wind_speed_m_s=3.4,wind_heading_relative_deg=6.1,
        circular_arc_representation='zero_thickness_surface',obstacle_height_m=6.)
    config=replace(config,initial_binary_cloud=replace(config.initial_binary_cloud,species_mole_fraction=0.),
        stratification_damping_alpha=5.,stratification_damping_exponent=1.,
        stratification_shear_regularization_s_inv=.041)
    with h5py.File(args.donor,'r') as source:
        times=source['time_s'][:];ze=source['source_z_edges_m'][:];ye=source['source_y_edges_m'][:]
        if not np.allclose(times,np.arange(4801)*.125):
            raise ValueError('control configuration requires the matching frozen 600s/0.125s donor')
    started=time.perf_counter()
    with ConservativeHandoffHdf5Writer3D(control_path,source_z_edges_m=ze,source_y_edges_m=ye,
            source_id='TI021 same-obstacle air-only control',metadata={'species_mole_fraction':0.,'duration_s':600.}) as writer:
        def observe(frame):
            writer(frame)
            if abs(frame.time_s/60-round(frame.time_s/60))<1e-10:
                print(json.dumps({'control_time_s':frame.time_s,'elapsed_s':time.perf_counter()-started}),flush=True)
        result=solve_resolved_obstacle_flow_3d(config,
            handoff_request=BinaryMixtureHandoffPlane3DRequest(130.,.125),handoff_observer=observe)
    control_runtime=time.perf_counter()-started
    (args.directory/'control-solver.json').write_text(json.dumps({
        'runtime_s':control_runtime,'maximum_domain_species_kg':float(np.max(result.domain_mass_kg)),
        'diagnostics':result.diagnostics},indent=2,allow_nan=False),encoding='utf-8')
    with h5py.File(args.donor,'r') as source,h5py.File(control_path,'r') as control,h5py.File(output,'x') as target:
        for key,value in source.attrs.items():target.attrs[key]=value
        for name in source:
            if name!='dense_layer_packet':source.copy(name,target)
        group=target.create_group('dense_layer_packet')
        species=source['dense_layer_packet/species_mass_flux_kg_m2_s'][:]
        old_mass=source['dense_layer_packet/mixture_mass_flux_kg_m2_s'][:]
        anomaly=source['total_mass_flux_anomaly_kg_m2_s'][:]-control['total_mass_flux_anomaly_kg_m2_s'][:]
        new_mass=np.maximum(anomaly,np.maximum(species,0.))
        for name in source['dense_layer_packet']:
            if name=='mixture_mass_flux_kg_m2_s':values=new_mass
            elif name=='species_mass_flux_kg_m2_s':values=species
            else:
                velocity=np.divide(source['dense_layer_packet/'+name][:],old_mass,
                    out=np.zeros_like(old_mass),where=old_mass>1e-30)
                values=new_mass*velocity
            group.create_dataset(name,data=values,compression='gzip',compression_opts=1,shuffle=True)
        target.attrs['marked_branch_reference']='matched air-only obstacle control'
        target.attrs['control_sha256']=hashlib.sha256(control_path.read_bytes()).hexdigest()
        area=np.diff(ze)[:,None]*np.diff(ye)[None,:]
        integrate=lambda values:float(np.trapezoid(np.sum(values*area,axis=(1,2)),times))
        report={'schema':'degali.ti021-control-subtracted-packets.v1',
            'source':str(args.donor.resolve()),'output':str(output.resolve()),
            'donor_sha256':hashlib.sha256(args.donor.read_bytes()).hexdigest(),
            'control_runtime_s':control_runtime,
            'original_marked_mixture_kg':integrate(old_mass),
            'control_subtracted_marked_mixture_kg':integrate(new_mass),
            'species_kg':integrate(species),
            'discarded_negative_mixture_anomaly_kg':integrate(np.maximum(-anomaly,0.)),
            'control_species_mass_max_kg':float(np.max(result.domain_mass_kg)),
            'control_diagnostics':result.diagnostics,
            'full_signed_mixture_handoff':False,'promotion_allowed':False}
    (args.directory/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({key:value for key,value in report.items() if key!='control_diagnostics'},indent=2))


if __name__=='__main__':main()
