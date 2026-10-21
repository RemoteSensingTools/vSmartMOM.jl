#!/usr/bin/env python3
"""Compare Float64 on native Float64 versus exact promoted Float32 grid nodes."""
import hashlib
import json
import os
import subprocess
import sys
import h5py
import numpy as np

root = sys.argv[1]
def toml(name):
    return json.loads(subprocess.check_output(['python3.11', '-c',
        'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
        os.path.join(root, name)]).decode())
metadata = toml('precision-grid.toml')
frozen = toml('precision-frozen.toml')
assert len(metadata['records']) == 2
assert all(r['m_used'] == 3 and r['ndoubl'] == frozen['records'][0]['ndoubl']
           for r in metadata['records'])
def load(name):
    with h5py.File(os.path.join(root, name + '.jld2')) as f:
        data = {k: np.array(f[k]) for k in ('y', 'K', 'state')}
    assert all(np.all(np.isfinite(v)) for v in data.values())
    return data
with h5py.File(os.path.join(root, 'state035_corrected_siffalse-reference.jld2')) as f:
    noise = np.sqrt(f['variance'][:934])
def compare(a,b):
    difference = (b['y']-a['y']) / noise
    return dict(max_noise_sigma=float(np.max(np.abs(difference))),
                rms_noise_sigma=float(np.sqrt(np.mean(difference**2))))
result = dict(metadata=metadata, fixed_state={})
grid = {}
for s in ('reference','optimized'):
    grid[s] = g = load('precision64-grid32-' + s)
    old = load('frozen-{}-prep32-rt64'.format(s))
    native = load('frozen-{}-prep64-rt64'.format(s))
    low = load('frozen-{}-prep32-rt32'.format(s))
    result['fixed_state'][s] = dict(
        grid_change_Float64=compare(native,g),
        preparation_RT64_with_grid_matched=compare(old,g),
        Float32_vs_Float64_with_grid_matched=compare(low,g))
a,b = grid['reference'],grid['optimized']
pred = a['K'].T.dot(b['state']-a['state'])
dy = b['y']-a['y']
result['between_states'] = dict(max_noise_sigma=float(np.max(np.abs(dy)/noise)),
    prediction_max_noise_sigma=float(np.max(np.abs(pred)/noise)),
    remainder_max_noise_sigma=float(np.max(np.abs(dy-pred)/noise)))
for r in metadata['records']:
    with open(os.path.join(root,r['file']),'rb') as f:
        assert hashlib.sha256(f.read()).hexdigest() == r['sha256']
print(json.dumps(result,indent=2))
