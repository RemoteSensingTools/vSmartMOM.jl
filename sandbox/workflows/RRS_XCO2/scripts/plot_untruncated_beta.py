#!/usr/bin/env python3
"""Plot the full untruncated Greek beta_l series."""

from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_beta")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "truth_map_aerosols"
DATA = ROOT / "untruncated_beta_757nm.dat"
OUTPUT = ROOT / "untruncated_beta_757nm.png"
NAMES = ("sulfate", "organic_carbon", "utls_sulfate", "mixture")

dtype = [("l", int), ("species", "U20"), ("weight", float), ("beta", float)]
rows = np.loadtxt(DATA, dtype=dtype, comments="#")
fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True,
                         constrained_layout=True)
for name in NAMES:
    r = rows[rows["species"] == name]
    style = dict(color="black", lw=2.4) if name == "mixture" else dict(lw=1.4)
    label = name.replace("_", " ")
    axes[0].plot(r["l"], r["beta"], label=label, **style)
    axes[1].semilogy(r["l"], np.abs(r["beta"]), label=label, **style)
axes[0].set_ylabel(r"$β_l$")
axes[0].set_title("Signed coefficient")
axes[1].set_ylabel(r"$|β_l|$")
axes[1].set_xlabel("Legendre degree l")
axes[1].set_title("Absolute magnitude")
for ax in axes:
    ax.grid(alpha=.3)
    ax.legend()
    ax.set_xlim(0, rows["l"].max())
fig.suptitle("Full untruncated aerosol Greek coefficient β at 757.001655 nm")
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
