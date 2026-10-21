"""Independent verification: python3 verify.py REPLAY_OUTPUT > verification.json."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np

helper_path = Path(__file__).resolve().parent.parent / 'selected_basis' / 'verify.py'
spec = importlib.util.spec_from_file_location('basis_verification', str(helper_path))
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def main():
    directory = Path(sys.argv[1])
    timings = helper.timing_file(directory / 'timings.toml')
    isolation = helper.timing_file(directory / 'isolation.toml')
    outputs = {name: helper.load_output(directory / (name + '-output.jld2'))
               for name in ('copied', 'shared', 'copied-b', 'shared-b', 'shared-a')}
    with h5py.File(timings['saved_state'], 'r') as data:
        archived_state = np.array(data['final_state'])
    assert np.array_equal(outputs['copied']['state'], archived_state)
    comparisons = {}
    for left, right in [('copied', 'shared'), ('copied-b', 'shared-b'), ('shared', 'shared-a')]:
        a, b = outputs[left], outputs[right]
        assert a['K'].shape == b['K'].shape == (30, 2742)
        assert np.array_equal(a['state'], b['state'])
        for field in ('y', 'K'):
            assert np.all(np.isfinite(a[field])) and np.all(np.isfinite(b[field]))
            assert a[field].dtype == b[field].dtype
            assert a[field].tobytes() == b[field].tobytes(), (left, right, field)
        comparisons[left + '_vs_' + right] = dict(measurements_bitwise_equal=True,
                                                  jacobian_bitwise_equal=True)
    assert not np.array_equal(outputs['shared']['state'], outputs['shared-b']['state'])
    assert not np.array_equal(outputs['shared']['y'], outputs['shared-b']['y'])
    assert all(isolation[name] for name in ('template_unchanged', 'table_coefficients_unchanged',
        'exact_copy_policy_parity', 'exact_perturbed_state_parity', 'exact_repeat_state_parity'))
    summary = dict(copy_probes=timings['copy_probes'], comparisons=comparisons,
                   isolation=isolation, timings={}, sha256={})
    for record in timings['records']:
        samples = record['samples']
        summary['timings'][record['mode']] = dict(
            seconds=float(np.median([s['seconds'] for s in samples])),
            host_bytes=float(np.median([s['host_bytes'] for s in samples])),
            rt_seconds=float(np.median([s['timing']['rt_linearized_seconds'] for s in samples])))
    for filename in [name + '-output.jld2' for name in outputs] + ['timings.toml', 'isolation.toml']:
        summary['sha256'][filename] = hashlib.sha256((directory / filename).read_bytes()).hexdigest()
    summary['verification'] = 'Independent NumPy/HDF5 bitwise checks and three-sample medians; Julia validates template and coefficient hashes.'
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
