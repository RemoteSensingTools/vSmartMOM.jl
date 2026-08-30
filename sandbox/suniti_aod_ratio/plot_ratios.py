from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def load(name):
    return np.genfromtxt(ROOT / f"{name}_spectra.csv", delimiter=",", names=True)


original = load("original")
fixed = load("fixed")

wn = original["wavenumber_cm1"]
assert np.array_equal(wn, fixed["wavenumber_cm1"])

ratios = {
    "original": (
        original["toa_aod003"] / original["toa_aod0"],
        original["boa_aod003"] / original["boa_aod0"],
    ),
    "fixed": (
        fixed["toa_aod003"] / fixed["toa_aod0"],
        fixed["boa_aod003"] / fixed["boa_aod0"],
    ),
}

fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True)
colors = {"original": "#c44e52", "fixed": "#2a6fbb"}
for label, (toa, boa) in ratios.items():
    axes[0].plot(wn, toa, lw=1.0, color=colors[label], label=label)
    axes[1].plot(wn, boa, lw=1.0, color=colors[label], label=label)

for ax in axes:
    ax.axhline(1.0, color="0.25", lw=0.8, ls="--")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    ax.set_ylabel("radiance ratio\n(AOD 0.03 / AOD 0)")

axes[0].set_title("TOA upwelling")
axes[0].set_yscale("log")
axes[1].set_title("BOA downwelling at VZA = 79.1° (off the direct-solar ordinate)")
axes[1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig.savefig(ROOT / "aod003_over_aod0_ratios.png", dpi=180)

with np.errstate(divide="ignore", invalid="ignore"):
    for endpoint, index in (("TOA", 0), ("BOA", 1)):
        old = ratios["original"][index]
        new = ratios["fixed"][index]
        delta = new - old
        finite = np.isfinite(delta)
        print(
            endpoint,
            "original ratio range", (np.min(old[np.isfinite(old)]), np.max(old[np.isfinite(old)])),
            "fixed ratio range", (np.min(new[np.isfinite(new)]), np.max(new[np.isfinite(new)])),
            "fixed-original max abs", np.max(np.abs(delta[finite])),
            "fixed-original RMS", np.sqrt(np.mean(delta[finite] ** 2)),
        )

print(ROOT / "aod003_over_aod0_ratios.png")
