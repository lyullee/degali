"""Record the actual low-T condensate composition needed by mixed-phase work."""
from datetime import datetime,timezone

import numpy as np

from stage2_matched_source import ROOT,REF,read,digest,write_new
from audit_mixed_air_liquid import independent_split
from run_preslhy_ambient_profile_audit import thermodynamics
from degali.addons.axisymmetric_jet import _air_phase_property_table
from stage1 import verify


def main():
    output=REF/'cold_condensate_inventory_2026-09-06.json'
    if output.exists(): raise FileExistsError(output)
    verify()
    paths=[REF/'mixed_air_liquid_property_screen_2026-09-06.json',REF/'e35_reduced.json',
        ROOT/'reference/air_mixture/Cambridge_MPT_Part6.pdf',ROOT/'reference/air_mixture/KIT_CryoPHAEQTS_2018.pdf',
        ROOT/'tools/audit_mixed_air_liquid.py',ROOT/'tools/run_preslhy_ambient_profile_audit.py',
        ROOT/'tools/audit_cold_condensate_inventory.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    previous=read(paths[0])
    for name,sha in previous['sha256'].items():
        if digest(ROOT/name)!=sha: raise ValueError(f'previous sample input changed:{name}')
    trials={r['trial']:r for r in read(paths[1])['trials']}
    table=_air_phase_property_table()
    cache,rows={},[]
    for old in previous['rows']:
        if old['status']!='below_liquid_domain': continue
        n,t,y=old['trial'],old['original_exact_T_K'],old['fraction']
        if n not in cache: cache[n]=thermodynamics(trials[n],consistent=True)
        th=cache[n]
        dry=(1-y)/(1+th.ambient_absolute_humidity)
        mass=dry*np.array([th._dry_nitrogen_mass_fraction,th._dry_oxygen_mass_fraction,th.ambient_absolute_humidity])
        mw=np.array([.0280134,.0319988,th.water_molecular_weight])
        k=np.array([np.interp(t,table['temperature'],table[s+'_saturation']) for s in ('nitrogen','oxygen','water')])/th.ambient_pressure
        gas,cond=independent_split(y/th.fuel_molecular_weight,mass/mw,k)
        if sum(cond[:2])<=0: raise ValueError('expected cold condensed inventory absent')
        rows.append(dict(trial=n,x_m=old['x_m'],kind=old['kind'],Gaussian_shape=old.get('Gaussian_shape'),
            exact_legacy_T_K=t,legacy_lookup_T_K=old['original_interpolated_T_K'],Y_H2=y,
            legacy_condensed_N2_kg_kg=float(cond[0]*mw[0]),legacy_condensed_O2_kg_kg=float(cond[1]*mw[1]),
            combined_condensed_N2_mole_fraction=float(cond[0]/sum(cond[:2])),
            gas_nitrogen_moles_kg=float(gas[0]),gas_oxygen_moles_kg=float(gas[1]),
            interpretation='Composition of the existing separate condensates combined on paper; NOT a new mixed-phase equilibrium composition.'))
    assert len(rows)==9
    verify()
    assert hashes=={str(p.relative_to(ROOT)):digest(p) for p in paths}
    write_new(output,dict(completed=True,created_utc=datetime.now(timezone.utc).isoformat(),input_sha256=hashes,
        rows=rows,summary=dict(n=len(rows),trials=sorted(cache),
            temperature_range_K=[min(r['exact_legacy_T_K'] for r in rows),max(r['exact_legacy_T_K'] for r in rows)],
            combined_condensed_N2_mole_fraction_range=[min(r['combined_condensed_N2_mole_fraction'] for r in rows),
                max(r['combined_condensed_N2_mole_fraction'] for r in rows)]),
        primary_author_source=dict(path='reference/air_mixture/Cambridge_MPT_Part6.pdf',
            title='Transformations in Solidified Gases',author='Charles S. Barrett',
            manuscript_received='1968-02-12',printed_pages='313-318',PDF_pages='45-50',
            relevant_figure='Fig.5, printed315/PDF47',
            note='Author account of experiments, NOT a retrieved copy of JCP48,2670. Figure is not digitized as a new quantitative closure; later1972 low-T erratum must be considered for that region.'),
        missing=['A quantitatively supported common mixed-solid/liquid chemical-potential and caloric representation for59--63K.',
            'Original1935 SLE table and its method/uncertainty; the retrieved2018 KIT plot cites it but is not original measurements.',
            'Kinetic justification that separately condensed air particles mix/equilibrate on the actual plume residence time.'],
        new_dispersion_field=False,candidate_promoted=False))
    print('Saved9 cold inventories:',output,flush=True)


if __name__=='__main__': main()
