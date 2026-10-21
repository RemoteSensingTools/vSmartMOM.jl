#!/usr/bin/env python3
"""Independent HDF5 verification of paired full OE results (no Julia types)."""
import sys, json, subprocess, os, hashlib
import h5py
import numpy as np

root = sys.argv[1]
records = json.loads(subprocess.check_output(['python3.11','-c',
    'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
    os.path.join(root,'results.toml')]).decode())['records']

def arrays(path):
    with h5py.File(path) as f:
        result = {k:np.array(f[k]) for k in ('state','y','K','posterior','averaging_kernel',
                    'states','costs','observation','variance','xa','Sa')}
        # JLD2 stores Julia Bool as an HDF5 bitfield, unsupported by old
        # h5py's dtype converter. Read the one-byte values without conversion.
        d = f['accepted']
        assert d.id.get_type().get_size() == 1
        accepted = np.empty(d.shape, dtype='u1')
        d.id.read(h5py.h5s.ALL,h5py.h5s.ALL,accepted,mtype=d.id.get_type())
        assert np.all((accepted == 0) | (accepted == 1))
        result['accepted'] = accepted
        return result

# Comparison thresholds set before observing paired results. They concern
# implementation equivalence, not acceptance of the underlying scientific fit.
limits = dict(xco2_ppm=0.01,state_prior_sigma=0.01,measurement_noise_sigma=0.01,
              relative_cost=1e-3,scaled_posterior_relative_l2=1e-3)
checks=[]
for case in sorted(set(r['case'] for r in records)):
    rr={r['mode']:r for r in records if r['case']==case}
    if set(rr)!=set(('reference','optimized')):
        raise RuntimeError('Incomplete pair '+case)
    a,b=[arrays(os.path.join(root,case+'-'+mode+'.jld2'))
         for mode in ('reference','optimized')]
    assert all(np.all(np.isfinite(v)) for data in (a,b) for v in data.values())
    for key in ('observation','variance','xa','Sa'):
        assert np.array_equal(a[key],b[key]),key
    for mode in rr:
        with open(rr[mode]['input'],'rb') as f:
            assert hashlib.sha256(f.read()).hexdigest()==rr[mode]['input_sha256']
    sigma=np.sqrt(np.diag(a['Sa']))
    scale=np.outer(sigma,sigma)
    pa,pb=a['posterior']/scale,b['posterior']/scale
    metrics=dict(xco2_ppm=abs(rr['reference']['xco2_ppm']-rr['optimized']['xco2_ppm']),
        state_prior_sigma=float(np.max(np.abs(a['state']-b['state'])/sigma)),
        measurement_noise_sigma=float(np.max(np.abs(a['y']-b['y'])/np.sqrt(a['variance']))),
        relative_cost=abs(rr['reference']['cost']-rr['optimized']['cost'])/max(1,abs(rr['reference']['cost'])),
        scaled_posterior_relative_l2=float(np.linalg.norm(pa-pb)/np.linalg.norm(pa)))
    same_decisions=(rr['reference']['outcome']==rr['optimized']['outcome'] and
                    np.array_equal(a['accepted'],b['accepted']))
    passed=same_decisions and all(metrics[k]<=limits[k] for k in limits)
    archive_comparison = None
    if not rr['reference']['synthetic_sif']:
        with h5py.File(rr['reference']['input']) as f:
            archived_xco2=float(np.array(f['XCO2']).reshape(-1)[0])
            archive_comparison=dict(
                xco2_change_ppm=rr['reference']['xco2_ppm']-archived_xco2,
                archived_outcome=int(np.array(f.attrs['outcome']).reshape(-1)[0]),
                archived_trials=int(np.array(f['trial_accepted']).size),
                archived_rejected=int(np.count_nonzero(np.array(f['trial_accepted'])==0)))
    checks.append(dict(case=case,metrics=metrics,same_decisions=same_decisions,passed=passed,
        archive_comparison=archive_comparison,
        reference=rr['reference'],optimized=rr['optimized']))
output_hashes={}
for name in sorted(os.listdir(root)):
    if name.endswith('.jld2'):
        with open(os.path.join(root,name),'rb') as f:
            output_hashes[name]=hashlib.sha256(f.read()).hexdigest()
result=dict(limits=limits,checks=checks,output_sha256=output_hashes,
            passed=len(checks)==6 and all(c['passed'] for c in checks))
print(json.dumps(result,indent=2))
sys.exit(0 if result['passed'] else 1)
