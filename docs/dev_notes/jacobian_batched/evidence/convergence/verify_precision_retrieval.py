#!/usr/bin/env python3
"""Verify paired full-inversion inputs, outcomes, costs, and state differences."""
import hashlib
import json
import os
import subprocess
import sys
import h5py
import numpy as np

root=sys.argv[1]
metadata=json.loads(subprocess.check_output(['python3.11','-c',
    'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
    os.path.join(root,'precision-retrieval.toml')]).decode())
assert len(metadata['records'])==3
def load(name):
    with h5py.File(os.path.join(root,name)) as f:
        return {k:f[k][:] for k in ('state','y','K','posterior','states','costs',
                                   'observation','variance','xa','Sa')}
baseline=load('state035_corrected_siffalse-optimized.jld2')
records=[]
initial=None
for r in metadata['records']:
    data=load(r['file'])
    assert all(np.all(np.isfinite(v)) for v in data.values())
    for k in ('observation','variance','xa','Sa'):
        assert np.array_equal(data[k],baseline[k])
    if initial is None:
        initial=data['states'][0]
    assert np.array_equal(data['states'][0],initial)
    if r['configuration']=='Float32':
        for k in ('state','y','K','posterior','states','costs'):
            assert np.array_equal(data[k],baseline[k])
    dx=data['state']-baseline['state']
    assert np.array_equal(dx,np.array(r['state_shift']))
    delta=data['state']-data['xa']
    cost=np.sum((data['y']-data['observation'])**2/data['variance']) + delta.dot(np.linalg.solve(data['Sa'],delta))
    assert np.isclose(cost,r['cost'],rtol=1e-12,atol=1e-9)
    assert np.isclose(np.max(np.abs(dx)/np.sqrt(np.diag(data['Sa']))),r['max_state_prior_sigma'],rtol=1e-12,atol=1e-14)
    assert r['trials']==len(data['costs'])==len(r['accepted'])
    with open(os.path.join(root,r['file']),'rb') as f:
        assert hashlib.sha256(f.read()).hexdigest()==r['sha256']
    records.append({k:r[k] for k in ('configuration','converged','outcome','trials','accepted',
        'xco2_common64_ppm','xco2_shift_common64_ppm','posterior_xco2_sigma_ppm',
        'psurf_shift_hpa','active_layer_co2_shift_ppm','max_state_prior_sigma','chi_squared','cost')})
for r in records:
    assert np.isclose(r['xco2_common64_ppm']-records[0]['xco2_common64_ppm'],
                      r['xco2_shift_common64_ppm'],rtol=0,atol=1e-12)
print(json.dumps(dict(identical_observation_noise_prior_initial_state=True,
    Float32_reproduces_archived_optimized_exactly=True,
    all_converged=all(r['converged'] for r in records),
    records=records,linearized_projections=metadata['linearized_projections']),indent=2))
