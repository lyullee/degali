"""Save an auditable non-slow regression snapshot without modifying Stage1."""
import argparse
from datetime import datetime,timezone
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

from stage2_matched_source import ROOT,digest,write_new
from stage1 import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    directory=args.directory.resolve()
    if directory.exists(): raise FileExistsError(directory)
    verify()
    paths=sorted(set(p for name in ('src','tests','tools') for p in (ROOT/name).rglob('*.py')))
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    command=[sys.executable,'-X','utf8','-m','pytest','-q','-m','not slow',
             '--junitxml='+str(directory/'pytest.xml')]
    directory.mkdir(parents=True,exist_ok=False)
    write_new(directory/'inputs.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),
        command=command,sha256=hashes,selection='all collected tests except marker slow',
        skipped_tests_are_not_passes=True,not_a_linux_oracle_completion=True))
    env=os.environ.copy()
    env['PYTHONPATH']=str(ROOT/'src')+os.pathsep+str(ROOT/'tools')
    env['OPENBLAS_NUM_THREADS']='1'
    env['PYTHONIOENCODING']='utf-8'
    print('Regression snapshot started:',directory,flush=True)
    with (directory/'pytest.log').open('x',encoding='utf-8') as log:
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        for line in process.stdout:
            log.write(line);log.flush()
            print(line,end='',flush=True)
        code=process.wait()
    mismatches=[name for name,value in hashes.items() if not (ROOT/name).exists() or digest(ROOT/name)!=value]
    added=sorted(str(p.relative_to(ROOT)) for name in ('src','tests','tools')
                 for p in (ROOT/name).rglob('*.py') if str(p.relative_to(ROOT)) not in hashes)
    report=ET.parse(directory/'pytest.xml').getroot()
    suites=[report] if report.tag=='testsuite' else list(report.iter('testsuite'))
    counts={key:sum(int(s.get(key,'0')) for s in suites) for key in ('tests','failures','errors','skipped')}
    counts['passed']=counts['tests']-counts['failures']-counts['errors']-counts['skipped']
    verify()
    write_new(directory/'complete.json',dict(completed=True,finished_utc=datetime.now(timezone.utc).isoformat(),
        exit_code=code,counts=counts,initial_code_files=len(hashes),changed_snapshot_files=mismatches,
        added_files_not_in_snapshot=added,original_snapshot_passed=code==0 and not mismatches,
        all_current_code_covered=code==0 and not mismatches and not added,
        slow_tests_excluded=True,linux_oracle_completed=False,candidate_promoted=False,
        output_sha256={p.name:digest(p) for p in (directory/'inputs.json',directory/'pytest.xml',directory/'pytest.log')}))
    print('Regression snapshot completed:',counts,'exit',code,flush=True)


if __name__=='__main__': main()
