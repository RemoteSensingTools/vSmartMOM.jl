#!/usr/bin/env python3
"""Plot per-Fourier-order TOA Stokes contributions at every angular stream."""

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm


def main():
    root = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path,
                        default=root / "untruncated_fourier_stream_contributions.jld2")
    parser.add_argument("--output", type=Path,
                        default=root / "untruncated_fourier_stream_contributions.png")
    args = parser.parse_args()

    table_path = args.input.with_suffix(".dat")
    table = np.loadtxt(table_path, comments="#")
    m = np.unique(table[:, 0]).astype(int)
    mu = table[table[:, 0] == m[0], 1]
    contrib = table[:, 2:5].reshape(len(m), len(mu), 3)
    elapsed = float(table[0, 5])
    order = np.argsort(mu)
    mu = mu[order]
    contrib = contrib[:, order, :]

    fig, axes = plt.subplots(3, 2, figsize=(14, 11), constrained_layout=True,
                             gridspec_kw={"width_ratios": [4.5, 1.5]})
    labels = ("I cosine contribution", "Q cosine contribution", "U sine amplitude")
    for s, label in enumerate(labels):
        z = contrib[:, :, s].T
        vmax = np.nanmax(np.abs(z))
        nonzero = np.abs(z[np.nonzero(z)])
        linthresh = max(vmax * 1e-7,
                        np.nanpercentile(nonzero, 2) if nonzero.size else 1e-30)
        norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax)
        mesh = axes[s, 0].pcolormesh(m, mu, z, shading="nearest",
                                     cmap="RdBu_r", norm=norm)
        fig.colorbar(mesh, ax=axes[s, 0], label=r"$w_m J_m^-(\mu)$")
        axes[s, 0].set_ylabel(r"Quadrature $\mu$")
        axes[s, 0].set_title(label)

        peak = np.nanmax(np.abs(contrib[:, :, s]), axis=1)
        rms = np.sqrt(np.nanmean(contrib[:, :, s] ** 2, axis=1))
        axes[s, 1].semilogy(m, peak, label="max |stream|", lw=1.2)
        axes[s, 1].semilogy(m, rms, label="stream RMS", lw=1.0)
        axes[s, 1].set_title("Angular magnitude")
        axes[s, 1].legend(fontsize=8)
        axes[s, 1].grid(alpha=.25)

    for ax in axes[-1, :]:
        ax.set_xlabel("Fourier order m")
    fig.suptitle(f"Untruncated O₂ A-band Fourier-stream diagnostic "
                 f"({elapsed / 60:.1f} min)")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print(args.output)


if __name__ == "__main__":
    main()
