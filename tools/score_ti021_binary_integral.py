"""Independent free-field dense-puff ablation on the same 145 TI021 sensors.

No obstacle claim: isolate restoring density/gravity/initial resting momentum
using existing integral coefficients, with no observed-data fitting.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import time
import numpy as np
from degali.addons.binary_integral_puff import initial_prism_puff, binary_puff_traces
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.finite_volume_initial_cloud_3d import RegularPolygonPrismInitialCloud
from degali.addons.finite_puff import GaussianPuffConfig
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook

ROOT=Path(__file__).resolve().parents[1]


def run(output,*,dt=.25,gravity=True,entrainment='constant'):
    started=time.perf_counter()
    workbook=ROOT/'ti21_input_ascii.xls'
    dataset=read_thorney_ti021_workbook(workbook)
    closure=BinaryIdealGasMixtureClosure('R12',.120913,606.0748985077907,
        ambient_temperature_k=291.15,ambient_pressure_pa=100500.)
    cloud=RegularPolygonPrismInitialCloud((60.,320.),14.,0.,13.,.32,291.15)
    angle=np.deg2rad(6.1)
    config=GaussianPuffConfig(duration_s=600.,time_step_s=dt,
        wind_velocity_m_s=(float(3.9*np.cos(angle)),float(3.9*np.sin(angle)),0.),
        width_cutoff_sigma=math.sqrt(3),gravity_spreading_coefficient=.75 if gravity else 0.,
        dense_entrainment_model=entrainment)
    puff=initial_prism_puff(cloud,closure,config)
    positions=np.array([(60.+s.local_y_m,320.-s.local_x_m,s.height_m) for s in dataset.dose_sensors])
    traces=binary_puff_traces(puff,positions)
    t=np.array([s.elapsed_s for s in puff.states])
    prediction_runtime=time.perf_counter()-started
    observed=np.array([s.dose_percent_s for s in dataset.dose_sensors])
    dose=np.trapezoid(traces,t,axis=1)
    positive=observed>0
    report={
        'schema':'degali.ti021-binary-integral-ablation.v1','gas':'R12',
        'status':'free_field_ablation_not_obstacle_validation',
        'obstacle_physics':False,'source_quiescent':True,'density_feedback':True,
        'gravity_spreading':gravity,'duration_s':600.,'time_step_s':dt,
        'dense_entrainment_model':entrainment,
        'species_mass_kg':puff.states[0].hydrogen_mass_kg,
        'maximum_relative_species_mass_residual':puff.maximum_relative_hydrogen_residual,
        'maximum_relative_total_mass_residual':puff.maximum_relative_mass_residual,
        'sensor_count':145,'positive_observed_sensor_count':int(positive.sum()),
        'dose_mape_percent':float(np.abs(dose[positive]/observed[positive]-1).mean()*100),
        'normalized_mean_bias':float((dose[positive]-observed[positive]).sum()/observed[positive].sum()),
        'prediction_runtime_s':prediction_runtime,
        'observation_unit':'vol%; dose vol% s','experimental_time_series_available':False,
        'configuration_used_observed_dose_peak_arrival_values':False,
        'workbook_sha256':hashlib.sha256(workbook.read_bytes()).hexdigest(),
        'coefficients':{name:getattr(config,name) for name in (
            'shear_entrainment_coefficient','ambient_turbulence_coefficient',
            'gravity_spreading_coefficient','form_drag_coefficient','width_cutoff_sigma')},
        'source_mapping':'exact prism inventory and second moments; Gaussian tails, not an exact top-hat',
        'sensor_rows':[{'label':f'TI021-S{i+1:03d}','observed_dose_vol_percent_s':float(observed[i]),
                       'predicted_dose_vol_percent_s':float(dose[i])} for i in range(len(dose))],
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),time_s=t,concentration_vol_percent=traces,
        observed_dose_vol_percent_s=observed,predicted_dose_vol_percent_s=dose,
        centre_m=np.array([s.centre_position_m for s in puff.states]),
        half_width_m=np.array([(s.longitudinal_half_length_m,s.lateral_half_width_m,s.vertical_half_width_m) for s in puff.states]),
        density_kg_m3=np.array([s.bulk_density_kg_m3 for s in puff.states]))
    report['runtime_including_output_s']=time.perf_counter()-started
    output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='sensor_rows'},indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--time-step-s',type=float,default=.25)
    p.add_argument('--no-gravity-spreading',action='store_true')
    p.add_argument('--entrainment',choices=['constant','parker_bulk_suppression'],default='constant')
    a=p.parse_args()
    run(a.output,dt=a.time_step_s,gravity=not a.no_gravity_spreading,entrainment=a.entrainment)
