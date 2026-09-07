#!/usr/bin/env python3
"""Plot saved pre/post-convolution residuals with common radiance normalization."""
import os
import sys
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root, output = sys.argv[1:3]
with h5py.File(os.path.join(root, 'convolution-stages.jld2')) as f:
    data = {k: f[k][:] for k in (
        'raw_wavelength_nm', 'dense_wavelength_nm', 'detector_wavelength_nm',
        'raw_rt_ppm', 'raw_matched_ppm', 'dense_rt_ppm', 'dense_matched_ppm',
        'detector_rt_noise', 'detector_matched_noise')}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex=True)
for column, key, title in ((0, 'rt', 'RT precision only; optical inputs fixed'),
                           (1, 'matched', 'Full precision change; spectral grid matched')):
    ax = axes[0, column]
    ax.plot(data['raw_wavelength_nm'], data['raw_' + key + '_ppm'],
            color='#aaaaaa', lw=0.6, label='Before convolution')
    ax.plot(data['dense_wavelength_nm'], data['dense_' + key + '_ppm'],
            color='#1766a3', lw=0.9, label='After 0.04 nm Gaussian convolution')
    ax.set_title(title)
    ax.set_ylabel('Radiance difference\n(ppm of common raw peak)')
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    ax = axes[1, column]
    ax.plot(data['detector_wavelength_nm'], data['detector_' + key + '_noise'],
            color='#1766a3', lw=0.9, label='Convolved detector measurements')
    ax.axhline(0.01, color='#ae4938', lw=0.7, ls='--', label='±0.01 noise σ comparison level')
    ax.axhline(-0.01, color='#ae4938', lw=0.7, ls='--')
    ax.set_ylabel('Difference / detector noise σ')
    ax.set_xlabel('Wavelength (nm)')
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    for row in range(2):
        axes[row, column].set_xlim(758, 772)
        axes[row, column].grid(alpha=0.15)
        axes[row, column].axhline(0, color='black', lw=0.4)
fig.suptitle('O₂ precision residuals: Float64 − Float32\nCase 035, reference terminal state; synthetic OCO instrument', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.92))
for extension in ('png', 'pdf'):
    fig.savefig(output + '.' + extension, dpi=180)
