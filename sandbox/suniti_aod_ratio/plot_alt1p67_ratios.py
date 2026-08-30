from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def load(name):
    return np.genfromtxt(
        ROOT / f"{name}_alt1p67_spectra.csv", delimiter=",", names=True
    )


original = load("original")
fixed = load("fixed")
wn = original["wavenumber_cm1"]
assert np.array_equal(wn, fixed["wavenumber_cm1"])

ratios = {
    "original": (
        original["up_aod003"] / original["up_aod0"],
        original["down_aod003"] / original["down_aod0"],
    ),
    "fixed": (
        fixed["up_aod003"] / fixed["up_aod0"],
        fixed["down_aod003"] / fixed["down_aod0"],
    ),
}

fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
colors = {"original": "#c44e52", "fixed": "#2a6fbb"}
for label, (up, _) in ratios.items():
    ax.plot(wn, up, lw=1.15, color=colors[label], label=label)

ax.axhline(1.0, color="0.25", lw=0.8, ls="--")
ax.grid(alpha=0.22)
ax.legend(frameon=False)
ax.set_ylabel("upwelling radiance ratio (AOD 0.03 / AOD 0)")
ax.set_title("Upwelling at 1.67 km above BOA")
ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
fig.savefig(ROOT / "aod003_over_aod0_ratios_alt1p67km.png", dpi=180)

for endpoint, index in (("up", 0), ("down", 1)):
    old = ratios["original"][index]
    new = ratios["fixed"][index]
    delta = new - old
    finite = np.isfinite(delta)
    print(
        endpoint,
        "original quantiles", np.quantile(old[np.isfinite(old)], [0, 0.5, 1]),
        "fixed quantiles", np.quantile(new[np.isfinite(new)], [0, 0.5, 1]),
        "fixed-original max abs", np.max(np.abs(delta[finite])),
        "fixed-original RMS", np.sqrt(np.mean(delta[finite] ** 2)),
    )

print(ROOT / "aod003_over_aod0_ratios_alt1p67km.png")
