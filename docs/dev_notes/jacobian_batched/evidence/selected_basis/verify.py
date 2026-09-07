"""Independently verify replay arrays and timings using NumPy/HDF5.

Usage: python3 verify.py BASELINE_DIR SELECTED_DIR > verification.json
On the benchmark host Python 3.6 supplies h5py/NumPy; Python 3.11 parses TOML.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


def timing_file(path):
    program = ('import json,sys,tomllib\n'
               'with open(sys.argv[1],"rb") as f: print(json.dumps(tomllib.load(f)))')
    return json.loads(subprocess.check_output(['python3.11', '-c', program, str(path)]))


def load_output(path):
    with h5py.File(str(path), 'r') as data:
        # Julia's column-major K is exposed as (parameter, measurement) by HDF5.
        return {name: np.array(data[name]) for name in ('y', 'K', 'state')}


def metrics(result, reference, noise):
    column_errors = np.linalg.norm(result['K'] - reference['K'], axis=1) / np.maximum(
        np.linalg.norm(reference['K'], axis=1), np.finfo(float).eps)
    return dict(max_noise_units=float(np.max(np.abs((result['y'] - reference['y']) / noise))),
                max_column_relative_l2=float(np.max(column_errors)),
                measurements_bitwise_equal=(result['y'].dtype == reference['y'].dtype and
                                            result['y'].tobytes() == reference['y'].tobytes()),
                jacobian_bitwise_equal=(result['K'].dtype == reference['K'].dtype and
                                       result['K'].tobytes() == reference['K'].tobytes()),
                column_relative_l2=column_errors.tolist())


def main():
    baseline, selected = map(Path, sys.argv[1:3])
    directories = {'baseline': baseline, 'selected': selected}
    times = {name: timing_file(path / 'timings.toml') for name, path in directories.items()}
    assert times['baseline']['saved_state'] == times['selected']['saved_state']
    with h5py.File(times['baseline']['saved_state'], 'r') as data:
        noise = np.array(data['noise_standard_deviation'])
        archived = dict(y=np.array(data['final_forward_model']),
                        K=np.array(data['final_jacobian']))
        state = np.array(data['final_state'])
    outputs = {name + '_' + mode: load_output(path / (mode + '-output.jld2'))
               for name, path in directories.items() for mode in ('matrix', 'source')}
    for output in outputs.values():
        assert np.array_equal(output['state'], state)
        assert output['K'].shape == (30, 2742)
        assert np.all(np.isfinite(output['K'])) and np.all(np.isfinite(output['y']))
    summary = {'timings': {}, 'comparisons': {}, 'sha256': {}}
    for name, records in times.items():
        for record in records['records']:
            samples = record['samples']
            summary['timings'][name + '_' + record['mode']] = {
                'seconds': float(np.median([s['seconds'] for s in samples])),
                'rt_seconds': float(np.median([s['timing']['rt_linearized_seconds'] for s in samples])),
                'host_bytes': float(np.median([s['host_bytes'] for s in samples])),
            }
    for result, reference in [('selected_matrix', 'baseline_matrix'),
                              ('selected_source', 'baseline_source'),
                              ('selected_source', 'selected_matrix')]:
        values = metrics(outputs[result], outputs[reference], noise)
        assert values['max_noise_units'] < 1e-3, (result, reference, values)
        assert values['max_column_relative_l2'] < 3e-4, (result, reference, values)
        summary['comparisons'][result + '_vs_' + reference] = values
    summary['comparisons']['selected_source_vs_archive'] = metrics(
        outputs['selected_source'], archived, noise)
    for name, path in directories.items():
        for filename in ('matrix-output.jld2', 'source-output.jld2', 'timings.toml'):
            summary['sha256'][name + '/' + filename] = hashlib.sha256((path / filename).read_bytes()).hexdigest()
    summary['verification'] = 'NumPy/HDF5 comparison of saved arrays; warmed three-sample medians from TOML.'
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
