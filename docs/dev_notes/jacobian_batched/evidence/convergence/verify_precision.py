#!/usr/bin/env python3
"""Compare the same two states in Float32/64 in the affected O2 band."""
import sys,os,json,hashlib,subprocess
import h5py
import numpy as np
root=sys.argv[1]
prefix=sys.argv[2] if len(sys.argv)>2 else 'precision64'
def load(name):
    with h5py.File(os.path.join(root,name+'.jld2')) as f:
        r={k:np.array(f[k]) for k in ('state','y','K')}
        assert all(np.all(np.isfinite(v)) for v in r.values())
        return r
high=[load(prefix+'-'+s) for s in ('reference','optimized')]
low=[load('diagnostic-'+s+'-optimized') for s in ('reference','optimized')]
n=len(high[0]['y'])
with h5py.File(os.path.join(root,'state035_corrected_siffalse-reference.jld2')) as f:
    noise=np.sqrt(f['variance'][:n])
result={'configuration':prefix}
for label,pair in (('Float32',low),('Float64',high)):
    a,b=pair
    dy=b['y'][:n]-a['y'][:n]
    pred=a['K'][:,:n].T.dot(b['state']-a['state'])
    result[label]=dict(max_noise_sigma=float(np.max(np.abs(dy)/noise)),
        prediction_max_noise_sigma=float(np.max(np.abs(pred)/noise)),
        remainder_max_noise_sigma=float(np.max(np.abs(dy-pred)/noise)))
result['fixed_state_precision_difference']={}
for s,lo,hi in zip(('reference','optimized'),low,high):
    assert np.array_equal(lo['state'],hi['state'])
    result['fixed_state_precision_difference'][s]=float(np.max(np.abs(lo['y'][:n]-hi['y'])/noise))
result['n_measurements']=n
result['output_sha256']={}
for s in ('reference','optimized'):
    name=prefix+'-'+s+'.jld2'
    with open(os.path.join(root,name),'rb') as f:
        result['output_sha256'][name]=hashlib.sha256(f.read()).hexdigest()
if prefix=='precision64-matched':
    def toml(name):
        return json.loads(subprocess.check_output(['python3.11','-c',
            'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
            os.path.join(root,name)]).decode())
    original={r['state']:r['ndoubl'] for r in toml('diagnostic-doubling.toml')['records'] if r['band']==1}
    matched=toml('precision-matched-doubling.toml')['records']
    assert all(r['ndoubl']==original[r['state']] for r in matched)
    result['matched_doubling_counts']=True
    result['numerical_settings']=matched
print(json.dumps(result,indent=2))
