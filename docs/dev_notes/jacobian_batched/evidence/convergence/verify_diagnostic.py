#!/usr/bin/env python3
"""Compare both solvers at both terminal states and with fixed Fourier order."""
import os,sys,json,subprocess,hashlib
import h5py
import numpy as np
root=sys.argv[1]
case='state035_corrected_siffalse'

def load(name):
    with h5py.File(os.path.join(root,name+'.jld2')) as f:
        result = {key:np.array(f[key]) for key in ('state','y','K')}
        assert all(np.all(np.isfinite(v)) for v in result.values())
        return result
with h5py.File(os.path.join(root,case+'-reference.jld2')) as f:
    noise=np.sqrt(np.array(f['variance']))

def compare(a,b):
    dy=b['y']-a['y']
    # HDF5 presents Julia's K in (state,measurement) order.
    prediction=a['K'].T.dot(b['state']-a['state'])
    return dict(max_noise_sigma=float(np.max(np.abs(dy)/noise)),
        prediction_max_noise_sigma=float(np.max(np.abs(prediction)/noise)),
        remainder_max_noise_sigma=float(np.max(np.abs(dy-prediction)/noise)))
result={'same_state':{},'between_states':{},'repeat_exact':{}}
for state in ('reference','optimized'):
    a,b=[load('diagnostic-'+state+'-'+mode) for mode in ('reference','optimized')]
    assert np.array_equal(a['state'],b['state'])
    result['same_state'][state]=compare(a,b)
    reference_norm=np.linalg.norm(a['K'],axis=1)
    difference_norm=np.linalg.norm(a['K']-b['K'],axis=1)
    assert np.all(difference_norm[reference_norm==0]==0)
    result['same_state'][state]['max_jacobian_column_relative_l2']=float(
        np.max(difference_norm[reference_norm>0]/reference_norm[reference_norm>0]))
    original=load(case+'-'+state)
    repeated=load('diagnostic-'+state+'-'+state)
    exact=all(np.array_equal(original[k],repeated[k]) for k in ('state','y','K'))
    assert exact
    result['repeat_exact'][state]=exact
for mode in ('reference','optimized','allmoments'):
    a,b=[load('diagnostic-'+state+'-'+mode) for state in ('reference','optimized')]
    result['between_states'][mode]=compare(a,b)
result['fourier']=json.loads(subprocess.check_output(['python3.11','-c',
    'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
    os.path.join(root,'diagnostic.toml')]).decode())
result['doubling']=json.loads(subprocess.check_output(['python3.11','-c',
    'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
    os.path.join(root,'diagnostic-doubling.toml')]).decode())
result['output_sha256']={}
for name in sorted(os.listdir(root)):
    if name.startswith('diagnostic-') and name.endswith('.jld2'):
        with open(os.path.join(root,name),'rb') as f:
            result['output_sha256'][name]=hashlib.sha256(f.read()).hexdigest()
print(json.dumps(result,indent=2))
