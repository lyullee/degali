"""Audit a completed binary 3D iteration without scoring a prefix as 600 s."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT = Path(__file__).resolve().parents[1]


def audit(directory: Path, output: Path) -> dict:
    record_path = directory / 'ti021-inner-transient.json'
    record = json.loads(record_path.read_text(encoding='utf-8'))
    dataset = read_thorney_ti021_workbook(ROOT/'ti21_input_ascii.xls')
    config = record['configuration']
    with np.load(directory/'ti021-inner-receptor-traces.npz', allow_pickle=False) as traces:
        times = traces['time_s']
        labels = traces['labels']
        positions = traces['solver_positions_m']
        concentrations = traces['concentration_vol_percent']
    full = bool(np.isclose(times[0],0.) and np.isclose(times[-1],600.))
    dose = np.trapezoid(concentrations, times, axis=1)
    sensor_indices = np.array([int(str(label).rsplit('S',1)[1])-1 for label in labels])
    observed = np.array([dataset.dose_sensors[i].dose_percent_s for i in sensor_indices])
    positive = observed > 0.
    with h5py.File(directory/'ti021-inner-handoff.h5','r') as handoff:
        ht = handoff['time_s'][:]
        flux = handoff['species_mass_flux_kg_m2_s'][:]
        ze, ye = handoff['source_z_edges_m'][:], handoff['source_y_edges_m'][:]
    dz, dy = np.diff(ze), np.diff(ye)
    area = dz[:,None]*dy[None,:]
    positive_flux = np.maximum(flux,0.)
    signed_mass = float(np.trapezoid(np.sum(flux*area,axis=(1,2)),ht))
    forward_mass = float(np.trapezoid(np.sum(positive_flux*area,axis=(1,2)),ht))
    above = np.clip((ze[1:]-np.maximum(ze[:-1],config['obstacle_height_m']))/dz,0.,1.)
    above_mass = float(np.trapezoid(np.sum(positive_flux*area*above[None,:,None],axis=(1,2)),ht))
    with np.load(directory/'ti021-inner-final-field.npz',allow_pickle=False) as final:
        fluid = ~final['solid_mask']
        u,w = final['velocity_u_m_s'], final['velocity_w_m_s']
        x,y,z = final['x_m'],final['y_m'],final['z_m']
        xx,yy = np.meshgrid(x,y)
        cx,cy = 60., config['domain_m'][1]/2
        near_wall = ((xx-cx)**2+(yy-cy)**2)**.5
        near_wall = (near_wall>=35.) & (near_wall<=65.) & (xx>=cx)
        near = fluid & near_wall[None,:,:]
        lee = fluid & (x[None,None,:] > cx+50.) & (z[:,None,None] < config['obstacle_height_m'])
        flow = {'near_wall_w_min_m_s': float(w[near].min()),
                'near_wall_w_max_m_s': float(w[near].max()),
                'near_wall_upward_cell_count': int(np.sum((w>1e-4)&near)),
                'near_wall_downward_cell_count': int(np.sum((w< -1e-4)&near)),
                'lee_reverse_u_cell_count': int(np.sum((u<0.)&lee)),
                'lee_minimum_u_m_s': float(u[lee].min()) if lee.any() else None}
    groups = {}
    for height in np.unique(positions[:,2]):
        mask = positive & np.isclose(positions[:,2],height)
        groups[f'{height:g}m'] = {'count': int(mask.sum()),
            'dose_mape_percent': float(np.mean(np.abs(dose[mask]/observed[mask]-1))*100)
                                 if full and mask.any() else None}
    baseline_path = ROOT/'outputs/ti021-baseline-after-align-default-2026-10-16.json'
    baseline = json.loads(baseline_path.read_text(encoding='utf-8'))
    baseline_by_label = {r['label']:r for r in baseline['sensor_rows']}
    baseline_dose = np.array([baseline_by_label[str(label)]['predicted_dose_vol_percent_s'] for label in labels])
    report = {
        'schema': 'degali.ti021-3d-iteration-audit.v1',
        'source_record': str(record_path.resolve()),
        'source_record_sha256': hashlib.sha256(record_path.read_bytes()).hexdigest(),
        'duration_s': float(times[-1]), 'complete_600s': full,
        'comparison_scope': 'local matched sensor subset; full-domain score unavailable',
        'sensor_count': len(labels), 'positive_sensor_count': int(positive.sum()),
        'dose_mape_percent': float(np.mean(np.abs(dose[positive]/observed[positive]-1))*100) if full else None,
        'normalized_mean_bias': float(np.sum(dose[positive]-observed[positive])/np.sum(observed[positive])) if full else None,
        'baseline_same_sensor_subset_mape_percent': float(np.mean(np.abs(baseline_dose[positive]/observed[positive]-1))*100),
        'baseline_all_84_sensor_mape_percent': baseline['dose_mape_percent'],
        'height_groups': groups,
        'runtime_s': record['runtime']['wall_clock_s'], 'numerics': record['numerics'],
        'downstream_plane_signed_species_mass_kg': signed_mass,
        'downstream_plane_forward_species_mass_kg': forward_mass,
        'plane_forward_mass_above_wall_kg': above_mass,
        'plane_forward_mass_above_wall_fraction': above_mass/forward_mass if forward_mass>0 else None,
        'plane_upper_layer_note': 'downstream plane flux; does not uniquely identify wall-crossing route',
        'final_flow': flow,
        'all_internal_execution_gates_pass': record['all_gates_pass'],
        'full_goal_achieved': False, 'promotion_allowed': False,
        'unresolved': ['full 145-sensor exterior continuation and matched path/height/distance nondegradation',
                       '3D-to-layer adapter maps lateral/vertical momentum into x/y; donor normal momentum unavailable',
                       'grid/time refinement near 0.4m and 2.4m sensors',
                       '30s end-to-end runtime target on the complete coupled path'],
        'sensor_rows': [{'label': str(label),'height_m': float(positions[i,2]),
            'predicted_dose_vol_percent_s': float(dose[i]),
            'observed_600s_dose_vol_percent_s': float(observed[i]),
            'absolute_percentage_error': float(abs(dose[i]/observed[i]-1)*100) if full and positive[i] else None}
            for i,label in enumerate(labels)]}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in {'sensor_rows','numerics'}},indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    audit(args.directory,args.output)
