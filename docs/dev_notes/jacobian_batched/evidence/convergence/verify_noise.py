#!/usr/bin/env python3
"""Trace noise sigma to the archived covariance and inspect signed precision errors."""
import hashlib
import json
import os
import sys
import h5py
import numpy as np

root, study = sys.argv[1:3]
archive = os.path.join(study, 'bottom_layer_XCO2_retrievals',
    'retrievals_acos_mapped_tapered_vertical_correlation_nosif', 'corrected',
    'retrieval_state035_perturbation10.nc')
with h5py.File(archive) as f:
    variance = f['Se_diagonal'][:]
    noise_path = f.attrs['source_noise_covariance'].decode()
    assert int(f['Se_diagonal'].attrs['frozen_during_retrieval'][0]) == 1
with h5py.File(noise_path) as f:
    assert np.array_equal(variance, f['Se_diagonal_corrected'][:])
    sigma = f['noise_std_corrected'][:]
    assert np.array_equal(sigma**2, variance)
    wavelength = f['wavelength'][:]
    signal = f['signal_photon_corrected'][:]
    band = f['band_index'][:].astype(int)-1
    maximum_signal = f['max_ms'][:][band]
    reconstructed_photon = maximum_signal/100 * np.sqrt(
        np.abs(100*signal/maximum_signal)*f['c_photon'][:]**2 + f['c_background'][:]**2)
    energy = (6.62607015e-34*299792458.0)/(wavelength*1e-9)
    reconstructed = reconstructed_photon*energy
    assert np.allclose(reconstructed,sigma,rtol=5e-15,atol=0)
    snr = f['snr_corrected'][:]
result = dict(definition='sigma_i = sqrt(Se_diagonal_i)',
    units='mW m-2 sr-1 nm-1', noise_frozen_during_retrieval=True,
    diagonal_covariance=True, covariance_matches_noise_source=True,
    photon_background_noise_reconstruction_max_relative=float(np.max(np.abs(reconstructed-sigma)/sigma)),
    o2_sigma_min=float(sigma[:934].min()),o2_sigma_median=float(np.median(sigma[:934])),
    o2_sigma_max=float(sigma[:934].max()),o2_snr_median=float(np.median(snr[:934])),
    signed_residuals=[], input_sha256={})
for state in ('reference','optimized'):
    def read(name):
        path=os.path.join(root,name+'.jld2')
        with h5py.File(path) as f:
            return f['y'][:]
    low=read('frozen-{}-prep32-rt32'.format(state))
    for label,name in (
        ('RT_only','frozen-{}-prep32-rt64'),
        ('full_matched_grid','precision64-grid32-{}'),
        ('full_native_grids','frozen-{}-prep64-rt64')):
        d=(read(name.format(state))-low)/sigma[:934]
        result['signed_residuals'].append(dict(state=state,comparison=label,
            mean_noise_sigma=float(d.mean()),rms_noise_sigma=float(np.sqrt(np.mean(d*d))),
            fraction_negative=float(np.mean(d<0)),lag1_correlation=float(np.corrcoef(d[:-1],d[1:])[0,1])))
for path in (archive,noise_path):
    with open(path,'rb') as f:
        result['input_sha256'][path]=hashlib.sha256(f.read()).hexdigest()
print(json.dumps(result,indent=2))
