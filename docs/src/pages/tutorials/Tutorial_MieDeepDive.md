```@meta
EditURL = "Tutorial_MieDeepDive.jl"
```

# Mie Scattering Deep Dive

This tutorial explains Mie theory in vSmartMOM, the NAI2 vs PCW methods, parameter selection,
numerical precision, phase reconstruction, delta-M truncation, and how aerosol optics feed into RT.

---

### 1. Introduction: What Mie Theory Computes

Mie theory solves Maxwell's equations for electromagnetic scattering by a homogeneous sphere.
For radiative transfer (RT), we need:

- **Extinction cross-section** ``C_{\mathrm{ext}}``: total attenuation (scattering + absorption)
- **Scattering cross-section** ``C_{\mathrm{sca}}``: energy scattered
- **Phase function** ``p(\cos\Theta)``: angular distribution of scattered light

The phase matrix (for polarized RT) is expanded in generalized spherical functions; the expansion
coefficients are the "Greek" coefficients ``(\alpha, \beta, \gamma, \delta, \epsilon, \zeta)``.
RT kernels work directly in this Fourier/Legendre space—angle-space values are reconstructed
only when needed (e.g., for plotting or diagnostics).

```julia
using vSmartMOM.Scattering
using Distributions
```

---

### 2. NAI2 vs PCW Methods

**NAI2 (Numerical Angular Integration)** — Siewert's method: directly integrates the phase
matrix over scattering angles to obtain Greek coefficients. No precomputation; suitable when
you compute a few aerosols or vary parameters often.

**PCW (Precomputed Wigner)** — Domke's method: precomputes Wigner d-functions for all
``(l,n,m)`` indices. The Greek coefficients are then evaluated via fast summations over
Mie indices. **Use PCW when you have many quadrature angles** (e.g., high-resolution RT
with many streams)—the one-time Wigner cost is amortized.

Both methods produce identical Greek coefficients; choose based on workflow:

| Scenario                    | Preferred method |
|-----------------------------|------------------|
| Single aerosol, few runs    | NAI2             |
| Many angles, repeated runs | PCW              |
| Parameter sweeps           | NAI2 (no Wigner) |

Example: NAI2 (no Wigner tables needed)

```julia
rₘ, σ, nᵣ, nᵢ = 0.3, 2.0, 1.3, 0.001
aero = Aerosol(LogNormal(log(rₘ), log(σ)), nᵣ, nᵢ)
λ = 0.55
model_NAI2 = make_mie_model(NAI2(), aero, λ, Stokes_IQUV(), δBGE(20, 2.0), 30.0, 500)
```

aerosol_optics = compute_aerosol_optical_properties(model_NAI2)

PCW requires precomputed Wigner tables (expensive one-time cost):
N_max = 120
wigner_A, wigner_B = compute_wigner_values(N_max)
model_PCW = make_mie_model(PCW(), aero, λ, Stokes_IQUV(), δBGE(20, 2.0), 30.0, 500, wigner_A, wigner_B)

---

### 3. Parameter Selection Guide

**`r_max`** — Upper bound for the size-distribution integral. Must cover the bulk of the
distribution; truncating too low biases cross-sections. Use ``\sim 5\sigma \cdot r_m``
for lognormal; check the info message about the fraction cut.

**`nquad_radius`** — Gauss–Legendre points over ``[0, r_{\max}]``. More points → better
accuracy for broad size distributions. Typical: 500–2000; narrow distributions need fewer.

**`l_max` / `l_trunc`** — Maximum Legendre order in the Greek expansion. Must exceed
the effective size parameter; ``l_{\max} \gtrsim 2\pi r_{\max}/\lambda``. Too low truncates
the phase function; too high adds cost with little gain.

**Trade-off**: Higher `nquad_radius` and `l_max` improve accuracy but increase runtime.
Start with moderate values (e.g., 500, 20) and increase if needed.

```julia
r_max = 30.0
nquad_radius = 500
l_max = 20
Δ_angle = 2.0
truncation_type = δBGE(l_max, Δ_angle)
```

---

### 4. Float32 vs Float64 Behavior

Mie computations use **downward recurrence** for the logarithmic derivative ``D_n`` and
**upward recurrence** for Bessel functions. Float32 can suffer from:

- Loss of precision in the recurrence for large size parameters
- Instability in ``a_n``, ``b_n`` for ``n \gg 1``

**Recommendation**: Use **Float64** for Mie. Float32 is acceptable only for small size
parameters (e.g., ``x = 2\pi r/\lambda \lesssim 10``) and when you have verified stability.
The RT model may use Float32 elsewhere; convert aerosol optics to Float64 for Mie, then
cast back if needed.

---

### 5. Phase Function Reconstruction

Greek coefficients encode the phase matrix in Legendre space. To get angle-dependent values:

```math
f_{11}(\mu) = P \cdot \beta,\quad f_{12} = P^2 \cdot (\mathrm{fac} \odot \gamma),\quad \ldots
```

Use `reconstruct_phase(greek_coefs, μ)` where `μ = cos(Θ)`. The scalar phase function is
`f₁₁`; it is normalized so ``\frac{1}{4\pi}\int p\,d\Omega = 1``.

**Visualization**: Plot `f₁₁` vs `μ` (or scattering angle Θ); use log scale for `f₁₁` to
see the forward peak. The ratio `f₁₂/f₁₁` gives the degree of linear polarization.

Minimal reconstruction example (no heavy compute):
μ, _ = gausslegendre(100)
smat = reconstruct_phase(aerosol_optics.greek_coefs, μ)
# plot(μ, smat.f₁₁, yscale=:log10); plot(μ, smat.f₁₂ ./ smat.f₁₁)

---

### 6. Delta-M Truncation

Strong forward peaks (large particles) require many Legendre terms. **Delta-M** (δ-BGE)
approximates the forward peak and renormalizes, reducing the effective expansion length.

**Truncation factor** ``f^t``: fraction of scattered energy moved to the "delta" term:
``f^t = 1 - c_0``, where ``c_0`` is the retained scattering fraction. The RT solver uses
``f^t`` to adjust the single-scattering albedo and phase function so that the truncated
expansion matches the original phase function outside the exclusion cone ``\Delta_\theta``.

Larger ``\Delta_\mathrm{angle}`` → more of the peak is fitted → smaller ``f^t`` and
fewer required streams. Typical: 1–3°.

---

### 7. Connection to RT

`AerosolOptics` provides exactly what the RT model needs:

- **`greek_coefs`** — Fourier/Legendre coefficients for the phase matrix
- **`ω̃`** — single-scattering albedo (``C_{\mathrm{sca}}/C_{\mathrm{ext}}``)
- **`k`** — extinction cross-section (per particle or per mass, depending on convention)
- **`fᵗ`** — delta-M truncation factor

These feed into the doubling-adding kernels. Aerosol layers are combined with molecular
(Rayleigh) scattering and absorption; the RT solver never needs angle-space phase
functions—only the Greek expansion.

```julia
println("Mie deep-dive tutorial complete. Key types: Aerosol, MieModel, AerosolOptics, GreekCoefs.")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

