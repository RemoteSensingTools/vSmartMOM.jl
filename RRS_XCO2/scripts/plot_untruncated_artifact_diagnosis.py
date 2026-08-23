#!/usr/bin/env python3
"""Diagnose the untruncated-vs-truncated angular discrepancy.

Four panels:
  (a) Greek beta_l decay -> how much physics the "untruncated" run actually adds.
  (b) TOA I vs signed VZA: converged dBGE ladder, the two untruncated runs, and
      the single-scattering angular template.
  (c) Fourier decomposition: the 50- vs 111-stream disagreement lives in m=0.
  (d) Float32 elemental-layer thinness dtau/mu_min vs nstreams, against the
      observed peak error of each run.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "truth_map_aerosols/"

# ---------------------------------------------------------------- load
def load_bym(f):
    d = np.loadtxt(D + f, skiprows=1)
    ms = np.unique(d[:, 0].astype(int)); vs = np.unique(d[:, 1])
    I = np.zeros((len(ms), len(vs))); cI = np.zeros_like(I)
    for r in d:
        j = np.searchsorted(vs, r[1])
        I[int(r[0]), j] = r[2]; cI[int(r[0]), j] = r[5]
    return ms, vs, I, cI

m, sv, I50, cI50 = load_bym("untruncated_stokes_by_m_nstreams50_mmax99.dat")
_, _,  I111, cI111 = load_bym("untruncated_stokes_by_m_nstreams111_mmax99.dat")

t = np.genfromtxt(D + "truncation_iqU.dat", dtype=None, encoding=None,
                  names=["vza", "vaz", "case", "I", "Q", "U"])
def signed(case):
    vz = np.unique(t["vza"])
    g = lambda az: np.array([t["I"][(t["vza"] == v) & (t["case"] == case)
                                    & (t["vaz"] == az)][0] for v in vz])
    a, b = g(180.0), g(0.0)
    return np.concatenate([a[:0:-1], b])

b = np.genfromtxt(D + "untruncated_beta_757nm.dat", dtype=None, encoding=None,
                  names=["l", "sp", "w", "b"])

# single-scattering angular template
sza = 30.0; mu0 = np.cos(np.deg2rad(sza))
mu = np.cos(np.deg2rad(np.abs(sv)))
cosT = -mu0*mu + np.sqrt(1-mu0**2)*np.sqrt(1-mu**2)*np.where(sv < 0, -1.0, 1.0)
Th = np.rad2deg(np.arccos(np.clip(cosT, -1, 1)))
p = np.genfromtxt(D + "aerosol_phase_functions_757nm.dat", dtype=None,
                  encoding=None, names=["a", "sp", "w", "f11", "f12", "q"])
mix = p["sp"] == "mixture"; ray = p["sp"] == "rayleigh"
ss = (0.28*0.99*np.interp(Th, p["a"][mix], p["f11"][mix])
      + 0.0245*np.interp(Th, p["a"][ray], p["f11"][ray])) / (4*(mu+mu0))

fig, ax = plt.subplots(2, 2, figsize=(15.5, 10.5))
fig.suptitle("Why the 'untruncated' TOA radiance disagrees with the converged "
             "$\\delta$BGE ladder\nO$_2$ A band, 757.001655 nm, SZA=30$^\\circ$; "
             "negative signed VZA = $\\Delta\\varphi$=180$^\\circ$", fontsize=13)

# ---- (a) beta decay
a0 = ax[0, 0]
for sp, c in [("mixture", "k"), ("organic_carbon", "tab:orange"),
              ("sulfate", "tab:blue"), ("utls_sulfate", "tab:green")]:
    s = b["sp"] == sp
    a0.semilogy(b["l"][s], np.abs(b["b"][s]), color=c, lw=1.4, label=sp)
for L, c in [(32, "tab:red"), (99, "tab:purple"), (221, "tab:brown")]:
    a0.axvline(L, color=c, ls="--", lw=1.2)
    a0.text(L, 3e-1, f"  l={L}", color=c, rotation=90, va="top", fontsize=9)
a0.axhspan(1e-16, 1e-12, color="0.85", zorder=0)
a0.text(120, 2e-14, "Float64 Mie noise floor", fontsize=9, color="0.35")
a0.set(xlabel="Legendre degree $l$", ylabel="$|\\beta_l|$", ylim=(1e-16, 5),
       xlim=(0, 221), title="(a) The 'untruncated' series adds no physics past $l\\approx50$")
a0.legend(fontsize=8); a0.grid(alpha=.3)

# ---- (b) radiance vs VZA
a1 = ax[0, 1]
for cap, c in zip(["8", "12", "16", "24", "32"],
                  plt.cm.viridis(np.linspace(0, .85, 5))):
    a1.plot(sv, signed(cap), color=c, lw=1.3, label=f"$\\delta$BGE l_trunc={cap}")
a1.plot(sv, cI50[-1], "k-", lw=2.2, label="untruncated, 50 streams")
a1.plot(sv, signed("untruncated"), "r-", lw=2.2,
        label="untruncated, 111 streams")
a1.axvline(-30, color="0.5", ls=":", lw=1.5)
a1.text(-29.5, 6.35, "$\\Theta=180^\\circ$\n(exact backscatter)", fontsize=9, color="0.35")
a1b = a1.twinx()
a1b.plot(sv, ss, color="tab:cyan", ls="--", lw=1.6)
a1b.set_ylabel("single-scattering template (arb.)", color="tab:cyan", fontsize=9)
a1b.tick_params(axis="y", colors="tab:cyan")
a1.set(xlabel="signed VZA [deg]", ylabel="TOA $I$",
       title="(b) The $\\delta$BGE ladder is converged and tracks the SS geometry")
a1.legend(fontsize=8, loc="upper left"); a1.grid(alpha=.3)

# ---- (c) m-decomposition
a2 = ax[1, 0]
a2.plot(sv, I50[0]-I111[0], "k-", lw=2.2,
        label="$m=0$:  50str $-$ 111str")
a2.plot(sv, (cI50[-1]-I50[0])-(cI111[-1]-I111[0]), "tab:orange", lw=2.2,
        label="$\\sum_{m>0}$:  50str $-$ 111str")
a2.plot(sv, cI50[-1]-signed("32"), "b--", lw=1.5,
        label="total: 50str $-$ converged $\\delta$BGE(32)")
a2.plot(sv, signed("untruncated")-signed("32"), "r--", lw=1.5,
        label="total: 111str $-$ converged $\\delta$BGE(32)")
a2.axhline(0, color="0.5", lw=1)
a2.set(xlabel="signed VZA [deg]", ylabel="$\\Delta I$ (absolute)",
       title="(c) The stream-count disagreement lives entirely in $m=0$")
a2.legend(fontsize=8, loc="lower left"); a2.grid(alpha=.3)
d0 = np.max(np.abs(I50[0]-I111[0]))
dg = np.max(np.abs((cI50[-1]-I50[0])-(cI111[-1]-I111[0])))
a2.text(.03, .94, f"max |$\\Delta$|:  $m$=0 $\\to$ {d0:.3f}   "
                  f"($\\sum_{{m>0}}$ $\\to$ {dg:.3f}, {d0/dg:.0f}$\\times$ smaller)",
        transform=a2.transAxes, fontsize=10, va="top",
        bbox=dict(fc="w", ec="0.6"))

# ---- (d) measured Float32 error vs elemental-layer thickness
a3 = ax[1, 1]
sw = np.genfromtxt(D + "doubling_sweep.dat", dtype=None, encoding=None,
                   names=["n", "thr", "floor", "flt", "mumin", "dmax",
                          "ratio", "I", "t"])
ref = sw["I"][(sw["flt"] == "Float64") & (sw["n"] == 111)][0]
f32 = sw[sw["flt"] == "Float32"]
f32 = f32[np.argsort(f32["dmax"])]
err = 100*np.abs(f32["I"] - ref)/ref
a3.loglog(f32["dmax"], err, "o-", color="tab:red", lw=2, ms=7,
          label="Float32 (measured)")
f64 = sw[sw["flt"] == "Float64"]
a3.loglog(f64["dmax"], np.maximum(100*np.abs(f64["I"]-ref)/ref, 1e-7),
          "s", color="tab:blue", ms=8, label="Float64 (measured)")
for r, e in zip(f32, err):
    if r["thr"] == 1e-3 and r["floor"] == 0:
        a3.annotate(f"{int(r['n'])} str\ndefault", (r["dmax"], e),
                    textcoords="offset points", xytext=(2, 12), fontsize=9,
                    color="tab:red", ha="center")
fl = f32[f32["floor"] > 0]
a3.axvline(fl["dmax"][0], color="tab:green", ls="--", lw=1.6)
a3.text(fl["dmax"][0]*0.75, 3e-5, "intended F32 floor\n$1024\\,\\epsilon$ = 1.2e-4",
        color="tab:green", fontsize=9, rotation=90, va="bottom", ha="right")
a3.axhline(0.1, color="0.5", ls=":", lw=1.2)
a3.text(2e-8, 0.12, "0.1% (δBGE ladder spread)", fontsize=9, color="0.35")
a3.set(xlabel="$d\\tau_{max}$ = elemental-layer optical depth",
       ylabel="|error| vs converged Float64  [%]",
       title="(d) Root cause: Float32 fails when the elemental layer is too THIN.\n"
             "Error collapses onto $d\\tau_{max}$, not onto nstreams.")
a3.legend(fontsize=9, loc="upper right"); a3.grid(alpha=.3, which="both")

fig.tight_layout(rect=[0, 0, 1, .945])
out = D + "untruncated_artifact_diagnosis.png"
fig.savefig(out, dpi=130)
print(out)
