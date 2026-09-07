"""Frozen manufactured operator/energy verification, never a trial score."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from CoolProp.CoolProp import PropsSI
from degali.core.jetplume import JetCoefficients
from degali.addons.axisymmetric_jet import AxisymmetricJetSource, ConservedGaussianJet, phase_ambient_from_rh
from degali.addons.buoyancy_profile import BuoyancyConstrainedEnthalpySection
from degali.addons.transverse_mixing import ConservativeTransverseMixing
from degali.addons.reservoir_thermal import ReservoirThermalMoments
from degali.addons.edge_enrichment import FixedTransportEdgeProjection
from degali.addons.enriched_transport import EnrichedModalTransport, PrescribedRadialMixing
from degali.addons.enriched_segments import encode_enriched
from degali.addons.exact_transverse_geometry import exact_geometry_field_view
from degali.addons.finite_tke_transport import FiniteTkeModalTransport
from degali.addons.finite_tke_energy_difference import actual_tke_moment_difference
from degali.addons.stable_enthalpy_difference import stable_enthalpy_difference
from degali.addons.paired_exact_kinetic_difference import paired_exact_kinetic_difference

ROOT = Path(__file__).resolve().parents[1]


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def manufactured_model():
    mw, humidity, rho = phase_ambient_from_rh(288.65, 101325., 53.6666667)
    phase = ConservedGaussianJet(
        AxisymmetricJetSource(diameter=.01, velocity=100., density=1., temperature=100.),
        ambient_temperature=288.65, ambient_pressure=101325., ambient_density=rho,
        fuel_molecular_weight=PropsSI('M', 'Hydrogen'), ambient_molecular_weight=mw,
        fuel_heat_capacity=14300., ambient_heat_capacity=1006.,
        ambient_absolute_humidity=humidity, consistent_phase_ambient=True,
        equilibrium_air_condensation=True, temperature_dependent_phase_enthalpy=True,
        conservative_establishment='entrained_mass', radial_points=41)
    jp = SimpleNamespace(th=None, k=JetCoefficients(sc=1.16**2), deltay=0., deltaz=0.,
        betay=1., betaz=1., gammaz=0., _split=lambda area, *_: (math.sqrt(area), math.sqrt(area)),
        _wind=lambda *_: 2., _wind_profile=lambda *_: (2., 0.))
    section = BuoyancyConstrainedEnthalpySection(jp, phase, thermal_width_ratio=1.04, quadrature_points=256)
    mixing = ConservativeTransverseMixing(section, np.array([1.1312, .1537, .02, .01, 25., 1., 1.5]))
    reservoir = ReservoirThermalMoments(mixing, thermal_species_ratio=1., mechanical_work='reduced_buoyancy_work')
    p = FixedTransportEdgeProjection(reservoir, reservoir.evaluate(), degree=4)
    supplied = PrescribedRadialMixing.from_weak_baseline(p, order=8)
    for name, value in dict(ustar=.2, zr=.1, rml=0., spread_floor=False).items():
        setattr(jp, name, value)
    p, par, geometry = exact_geometry_field_view(p, encode_enriched(p, np.zeros(p.count)))
    mean = EnrichedModalTransport(p, par, scalar_mixing=supplied, thermal_species_ratio=1.,
                                 mechanical_work='reduced_buoyancy_work_immediate_shear_heat')
    model = FiniteTkeModalTransport(mean, np.r_[math.log(2.), np.zeros(mean.size)], ambient_tke=.07,
        tke_diffusivity=.3, dissipation_time=.2, circulation_amplitudes=[1e-4, -1e-4, 1e-4, 0.])
    return model, geometry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    partial = args.output.with_name(args.output.stem+'.partial.json')
    if args.output.exists() or partial.exists():
        parser.error('refusing to overwrite frozen manufactured evidence')
    paths = list((ROOT/'src').rglob('*.py')) + [Path(__file__).resolve()]
    paths += [ROOT/name for name in ('tests/test_finite_tke_transport.py', 'tests/test_finite_tke_energy_difference.py',
        'docs/prereg-finite-tke-modal-operator.md', 'docs/prereg-finite-tke-manufactured-audit.md')]
    def hashes():
        return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    initial = hashes()
    data = dict(phase='manufactured_finite_tke_operator', completed=False, numerical_passed=False,
        physical_closure_passed=False, field_scored=False, adopted=False,
        manufactured_inputs_not_trial_estimates=True, sha256=initial,
        started_utc=datetime.now(timezone.utc).isoformat(), operators=[], differences=[])
    def checkpoint():
        partial.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    def progress(item):
        if item['completed_angles'] % 256 == 0 or item['completed_angles'] == item['total_angles']:
            data['progress'] = item
            checkpoint()
            print(data['stage'], item, flush=True)
    checkpoint()
    try:
        model, geometry = manufactured_model()
        data['geometry'] = geometry
        data['explicit_inputs'] = dict(tke_parameters=model.parameters, ambient_tke=model.ambient_tke,
            tke_diffusivity=model.diffusivity, dissipation_time=model.dissipation_time,
            circulation_amplitudes=model.amplitudes)
        for order, angular in ((4, 8), (8, 16)):
            data['stage'] = f'operator_{order}_{angular}'
            data['operators'].append(model.assemble(order=order, angular_order=angular, callback=progress))
            checkpoint()
        coarse, fine = data['operators']
        data['refinements'] = {key: float(np.max(abs(fine[key]-coarse[key])/np.maximum(abs(fine[key]), 1.)))
                               for key in ('matrix', 'right', 'rates')}
        rates = fine['rates']
        data['cross_grid_errors'] = [float(np.max(abs(op['matrix'] @ rates-op['right'])/np.maximum(abs(op['right']), 1.)))
                                      for op in data['operators']]
        data['independent_diagnostics'] = model.independent_diagnostics(rates, angles=65, order=16)
        ledger = fine['ledgers']
        terms = ledger['terms']
        boundary_scale = max(1., *(abs(terms[key]) for key in ('heat_boundary', 'mean_kinetic_boundary', 'tke_boundary')))
        data['boundary_identity_scaled_error'] = abs(ledger['natural_boundary_identity'])/boundary_scale
        data['mass_boundary_scaled_error'] = abs(ledger['mass_boundary_error'])/max(1., abs(fine['source'][0]))
        data['operator_passed'] = bool(max(data['refinements'].values()) <= 1e-5
            and max(data['cross_grid_errors']) <= 1e-5 and max(ledger['scaled_errors'].values()) <= 1e-8
            and max(data['boundary_identity_scaled_error'], data['mass_boundary_scaled_error']) <= 1e-8)
        checkpoint()
        if not data['operator_passed']:
            raise ValueError('manufactured coupled operator or conservation gate failed')
        expected = float(fine['advective_jacobian'][4] @ rates)
        scale = max(abs(expected), 1.)
        h = 2e-6/max(1., float(max(abs(rates))))
        data.update(expected_total_energy_derivative=expected, original_step=h)
        b = model.base
        for factor in (.5, 1.):
            enthalpy = stable_enthalpy_difference(b.projection, b.parameters, rates[:model.mean_count], h*factor,
                order=64, precision=70, wind_policy='exact_constraint')
            tke = actual_tke_moment_difference(model, rates, h*factor, order=64)
            for order, angular in ((4, 8), (8, 16)):
                data['stage'] = f'actual_energy_{factor}_{order}_{angular}'
                checkpoint()
                kinetic = paired_exact_kinetic_difference(b.projection, b.parameters, rates[:model.mean_count], h*factor,
                    order=order, angular_order=angular, callback=progress)
                derivative = math.fsum([enthalpy['derivative'], kinetic['derivative'], tke['derivative']])
                row = dict(factor=factor, enthalpy=enthalpy, kinetic=kinetic, tke=tke,
                    total_derivative=derivative, energy_scaled_error=abs(derivative-expected)/scale)
                data['differences'].append(row)
                checkpoint()
                print('actual H+meanKE+Q', factor, order, angular, row['energy_scaled_error'], flush=True)
        values = [r['total_derivative'] for r in data['differences']]
        data['all_difference_refinement'] = (max(values)-min(values))/scale
        data['energy_difference_passed'] = bool(max(data['all_difference_refinement'],
            *(r['energy_scaled_error'] for r in data['differences'])) <= 1e-5)
        data['numerical_passed'] = data['operator_passed'] and data['energy_difference_passed']
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        data['failure'] = str(exc)
    if hashes() != initial:
        raise RuntimeError('frozen finite-TKE dependencies changed during audit')
    data.update(completed=True, stage='completed', finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('manufactured finite-TKE completed', data['numerical_passed'], data.get('failure'), flush=True)


if __name__ == '__main__':
    main()
