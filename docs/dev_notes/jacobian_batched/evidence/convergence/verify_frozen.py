#!/usr/bin/env python3
"""Verify frozen-boundary controls and quantify O2 precision contributions."""
import hashlib
import json
import os
import subprocess
import sys

import h5py
import numpy as np

root = sys.argv[1]
metadata = json.loads(subprocess.check_output([
    'python3.11', '-c',
    'import tomllib,json,sys; print(json.dumps(tomllib.load(open(sys.argv[1],"rb"))))',
    os.path.join(root, 'precision-frozen.toml')]).decode())
assert len(metadata['records']) == 8
assert sum(r['identity_checked'] for r in metadata['records']) == 4
assert all(r['m_used'] == 3 for r in metadata['records'])
assert all(r['ndoubl'] == metadata['records'][0]['ndoubl'] for r in metadata['records'])

def load(name):
    with h5py.File(os.path.join(root, name + '.jld2')) as f:
        result = {k: np.array(f[k]) for k in ('state', 'y', 'K')}
    assert all(np.all(np.isfinite(v)) for v in result.values())
    return result

data = {(s, p, r): load('frozen-{}-prep{}-rt{}'.format(s, p, r))
        for s in ('reference', 'optimized') for p in (32, 64) for r in (32, 64)}
n = len(data['reference', 32, 32]['y'])
with h5py.File(os.path.join(root, 'state035_corrected_siffalse-reference.jld2')) as f:
    noise = np.sqrt(f['variance'][:n])

def measure(a, b):
    difference = (b['y'] - a['y']) / noise
    norm = np.linalg.norm(a['K'], axis=1)
    active = norm > 0
    return dict(max_noise_sigma=float(np.max(np.abs(difference))),
                rms_noise_sigma=float(np.sqrt(np.mean(difference**2))),
                worst_column_relative_l2=float(np.max(
                    np.linalg.norm(b['K']-a['K'], axis=1)[active] / norm[active])))

result = dict(n_measurements=n, fourier_and_doubling_counts_match=True,
              native_identity_checks=4, fixed_state={}, between_states={}, sha256={})
for s in ('reference', 'optimized'):
    # The clones must also reproduce the independently saved earlier probes.
    for ft in (32, 64):
        historical = load(('diagnostic-{}-optimized' if ft == 32 else
                           'precision64-matched-{}').format(s))
        assert np.array_equal(data[s, ft, ft]['y'], historical['y'][:n])
        assert np.array_equal(data[s, ft, ft]['K'], historical['K'][:, :n])
    result['fixed_state'][s] = {
        'RT_only_preparation32': measure(data[s, 32, 32], data[s, 32, 64]),
        'RT_only_preparation64': measure(data[s, 64, 32], data[s, 64, 64]),
        'preparation_and_instrument_RT64': measure(data[s, 32, 64], data[s, 64, 64]),
        'full_workflow': measure(data[s, 32, 32], data[s, 64, 64])}
for p in (32, 64):
    for r in (32, 64):
        a, b = (data[s, p, r] for s in ('reference', 'optimized'))
        dy = b['y'] - a['y']
        prediction = a['K'].T.dot(b['state'] - a['state'])
        result['between_states']['prep{}_rt{}'.format(p, r)] = dict(
            max_noise_sigma=float(np.max(np.abs(dy) / noise)),
            prediction_max_noise_sigma=float(np.max(np.abs(prediction) / noise)),
            remainder_max_noise_sigma=float(np.max(np.abs(dy-prediction) / noise)))
for record in metadata['records']:
    path = os.path.join(root, record['file'])
    with open(path, 'rb') as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    assert digest == record['sha256']
    result['sha256'][record['file']] = digest
print(json.dumps(result, indent=2))
