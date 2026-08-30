from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
COLORS = {"original": "#c44e52", "fixed": "#2a6fbb"}


def load(kind, altitude=False):
    suffix = "_alt1p67" if altitude else ""
    return np.genfromtxt(
        ROOT / f"so4_{kind}{suffix}_dry_spectra.csv", delimiter=",", names=True
    )


def ratio(data):
    return data["up_aod0"] / data["up_aod003"]


datasets = {
    "TOA": {kind: load(kind) for kind in ("original", "fixed")},
    "1.67 km above BOA": {
        kind: load(kind, altitude=True) for kind in ("original", "fixed")
    },
}

fig, axes = plt.subplots(
    2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True
)
for ax, (location, cases) in zip(axes, datasets.items()):
    plotted = []
    for kind, data in cases.items():
        values = ratio(data)
        plotted.append(values[np.isfinite(values)])
        ax.plot(
            data["wavenumber_cm1"], values,
            color=COLORS[kind], lw=1.15, label=kind,
        )
        finite = values[np.isfinite(values)]
        print(location, kind, "min median mean max", *map(float, (
            np.min(finite), np.median(finite), np.mean(finite), np.max(finite)
        )))

    limits = np.concatenate(plotted)
    padding = 0.04 * (np.max(limits) - np.min(limits))
    ax.set_ylim(np.min(limits) - padding, np.max(limits) + padding)
    ax.axhline(1.0, color="0.25", lw=0.8, ls=":")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    ax.set_ylabel("upwelling ratio\n(AOD 0 / AOD 0.03)")
    ax.set_title(location)
axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig.suptitle("SO$_4$ upwelling — dry atmosphere (q = 0)")
output = ROOT / "so4_dry_aod0_over_aod003_upwelling_TOA_and_alt1p67km.png"
fig.savefig(output, dpi=180)
print(output)
