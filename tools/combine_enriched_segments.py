"""Combine the complete seven-case guarded-segment screen, including failures."""

import argparse
import hashlib
import json
from pathlib import Path
from audit_transverse_mixing import ROOT


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('groups',type=Path,nargs='+')
    args=parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite combined segment evidence')
    rows,sha=[],{}
    for path in args.groups:
        d=json.loads(path.read_text(encoding='utf-8'))
        if not d['completed'] or d['phase']!='guarded_enriched_short_segments':
            raise ValueError('completed compatible segment groups required')
        if d['selected_trials']!=[r['trial'] for r in d['rows']]:
            raise ValueError('incomplete selected group')
        for name,digest in d['sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest or (name in sha and sha[name]!=digest):
                raise ValueError(f'inconsistent frozen dependency: {name}')
            sha[name]=digest
        rows.extend(d['rows'])
        sha[str(path.resolve().relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    if sorted(r['trial'] for r in rows)!=[10,11,12,22,23,24,25]:
        raise ValueError('all7 trials required exactly once, including failures')
    rows.sort(key=lambda r:r['trial'])
    sha[str(Path(__file__).resolve().relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(phase='combined_guarded_enriched_short_segments',completed=True,rows=rows,sha256=sha,
        passed_trials=[r['trial'] for r in rows if r['primary_passed']],selected_trials=[r['trial'] for r in rows],
        all_passed=all(r['primary_passed'] for r in rows),full_downstream_integrated=False,field_scored=False,promoted=False)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('primary passed',out['passed_trials'],'hashes',len(sha),flush=True)


if __name__=='__main__':
    main()
