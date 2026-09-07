#!/usr/bin/env python3
"""Report historical campaign timing separately from current-code A/B timing."""
import hashlib
import json
import os
import statistics
import subprocess
import h5py
import numpy as np

here=os.path.dirname(os.path.abspath(__file__))
def toml(path):
    return json.loads(subprocess.check_output(['python3.11','-c',
        'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',path]).decode())
verification=json.load(open(os.path.join(here,'verification.json')))
results=[]
for c in verification['checks']:
    if c['reference']['synthetic_sif']:
        continue
    old_path=c['reference']['input']
    with open(old_path,'rb') as f:
        assert hashlib.sha256(f.read()).hexdigest()==c['reference']['input_sha256']
    with h5py.File(old_path) as f:
        terminal=float(np.asarray(f.attrs['final_evaluation_seconds']).reshape(-1)[0])
        old=float(f['evaluation_seconds'][:].sum()+terminal)
        assert np.array_equal(f['trial_accepted'][:].astype(bool),
            np.array([True,True,True] if '001' in c['case'] else [True,False,True,True,True]))
    new=c['optimized']['seconds']
    results.append(dict(case=c['case'],historical_evaluation_total_seconds=old,
        optimized_full_solve_seconds=new,historical_ratio=old/new,
        current_reference_full_solve_seconds=c['reference']['seconds'],
        current_reference_ratio=c['reference']['seconds']/new,
        historical_terminal_seconds=terminal,
        first_pair_includes_compilation=c['case']=='state001_corrected_siffalse'))
lut=toml(os.path.join(here,'../shared_luts/timings.toml'))
shared=next(r for r in lut['records'] if r['mode']=='shared')
new_evaluation=statistics.median(r['seconds'] for r in shared['samples'])
aerosol=next(r for r in results if r['case']=='state035_corrected_siffalse')
print(json.dumps(dict(cases=results,complete_evaluation=dict(
    historical_seconds=aerosol['historical_terminal_seconds'],
    optimized_warmed_median_seconds=new_evaluation,
    historical_ratio=aerosol['historical_terminal_seconds']/new_evaluation),
    caveat='Archived campaign versus isolated replay; original branch was not rerun under matched worker load. Historical totals sum model evaluations including terminal evaluation, excluding startup and campaign I/O.'),indent=2))
