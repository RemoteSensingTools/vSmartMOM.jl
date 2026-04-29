```@meta
EditURL = "Tutorial_Scattering.jl"
```

# Scattering: Mie Phase Function Tutorial

### Introduction
This tutorial shows the standard Mie workflow in `vSmartMOM.Scattering`:

1. Define an aerosol size distribution and refractive index
2. Build a `MieModel` (NAI2 in this example)
3. Compute aerosol optical properties (`GreekCoefs`, `ω̃`, `k`, `fᵗ`)
4. Reconstruct phase matrix elements on a custom angular grid
5. Run optional AD and cross-section convenience routines

---

### Load packages

```julia
using vSmartMOM.Scattering
using Distributions
using FastGaussQuadrature
using Parameters
using CairoMakie
```

---

### Define aerosol properties

```julia
rₘ = 0.3                 # median radius [μm]
σ = 2.0                  # geometric stddev [-]
nᵣ = 1.3                 # real refractive index
nᵢ = 0.001               # imaginary refractive index (use positive convention)
r_max = 30.0             # upper integration radius [μm]
nquad_radius = 1500      # size quadrature points

size_distribution = LogNormal(log(rₘ), log(σ))
aero = Aerosol(size_distribution, nᵣ, nᵢ)
```

---

### Create model settings

```julia
λ = 0.55                               # wavelength [μm]
polarization_type = Stokes_IQUV()
l_max = 20
Δ_angle = 2.0
truncation_type = δBGE(l_max, Δ_angle)

model_NAI2 = make_mie_model(
    NAI2(),
    aero,
    λ,
    polarization_type,
    truncation_type,
    r_max,
    nquad_radius,
)
```

---

### Compute aerosol optical properties

```julia
aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2)

println("ω̃ = ", aerosol_optics_NAI2.ω̃)
println("k  = ", aerosol_optics_NAI2.k)
println("fᵗ = ", aerosol_optics_NAI2.fᵗ)
```

---

### Inspect Greek coefficients

```julia
(; α, β, γ, δ, ϵ, ζ) = aerosol_optics_NAI2.greek_coefs

fig = Figure(size=(700, 700))
coef_names = ["α", "β", "γ", "δ", "ϵ", "ζ"]
coef_data  = [α, β, γ, δ, ϵ, ζ]
for i in 1:6
    row, col = fldmod1(i, 2)
    ax = Axis(fig[row, col], title=coef_names[i], xlabel="l")
    lines!(ax, coef_data[i])
    xlims!(ax, 0, 100)
end
fig
```

These coefficients are Fourier-space quantities used by RT kernels.
Reconstruct angle-dependent phase-matrix elements only when needed.

---

### Reconstruct phase matrix elements

```julia
μ_quad, _ = gausslegendre(500)
scattering_matrix = reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ_quad)
(; f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄) = scattering_matrix
```

Plot the phase function and degree of linear polarization as a function of scattering angle:

```julia
Θ = rad2deg.(acos.(μ_quad))

fig = Figure(size=(700, 600))
ax1 = Axis(fig[1,1], ylabel="f₁₁", yscale=log10)
lines!(ax1, Θ, f₁₁)

ax2 = Axis(fig[2,1], ylabel="f₁₂ / f₁₁", xlabel="Scattering angle (°)")
lines!(ax2, Θ, f₁₂ ./ f₁₁)
fig
```

The forward peak in f₁₁ is characteristic of Mie scattering by particles larger than the wavelength.
The f₁₂/f₁₁ ratio gives the degree of linear polarization for unpolarized incident light.

The docs build also writes a standalone Plotly version of this phase-function
preview. It keeps the log-scaled f₁₁ panel and the polarization-ratio panel
interactive on the rendered documentation page.

```@raw html
<iframe title="Mie phase function preview" src="../../assets/plots/scattering_phase_preview.html" loading="lazy" style="width: 100%; height: 560px; border: 1px solid var(--vp-c-divider); border-radius: 8px;"></iframe>
```

Polar view of the phase function and polarization ratio:

```julia
fig = Figure(size=(800, 400))
θ_rad = acos.(μ_quad)
θ_full = vcat(θ_rad, .-reverse(θ_rad))
f₁₁_full = vcat(f₁₁, reverse(f₁₁))
ratio = f₁₂ ./ f₁₁
ratio_full = vcat(abs.(ratio), reverse(abs.(ratio)))

ax1 = PolarAxis(fig[1,1], title="log₁₀(f₁₁)")
lines!(ax1, θ_full, log10.(f₁₁_full))

ax2 = PolarAxis(fig[1,2], title="|f₁₂/f₁₁|")
lines!(ax2, θ_full, ratio_full)
fig
```

---

### AD mode

`autodiff=true` returns one `AerosolOptics` object.
The Jacobian is stored in `derivs` with columns corresponding to:
`[rₘ, σ, nᵣ, nᵢ]`.

```julia
aerosol_optics_ad = compute_aerosol_optical_properties(model_NAI2; autodiff=true)
println("AD derivs size: ", size(aerosol_optics_ad.derivs))
```

---

### Convenience APIs for cross-sections and scalar phase function

```julia
μ_pf, w_μ_pf, p11, C_ext, C_sca, g = phase_function(aero, λ, r_max, nquad_radius)
println("phase_function C_ext = ", C_ext, ", C_sca = ", C_sca, ", g = ", g)

XS_ext, XS_sca, Cext_eff, Csca_eff = compute_aerosol_XS(aero, λ, r_max, nquad_radius)
println("compute_aerosol_XS XS_ext = ", XS_ext, ", XS_sca = ", XS_sca)
println("compute_aerosol_XS Cext_eff = ", Cext_eff, ", Csca_eff = ", Csca_eff)
```

---

### Optional PCW workflow

If you want to run the PCW decomposition directly, precompute (or load) Wigner tables:

```julia
N_max = 120
wigner_A, wigner_B = compute_wigner_values(N_max)   # can be expensive
model_PCW = make_mie_model(
    PCW(),
    aero,
    λ,
    polarization_type,
    truncation_type,
    r_max,
    nquad_radius,
    wigner_A,
    wigner_B,
)
aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW)
```

---

### Optional animation (heavy)

Set `ENV["VSMARTMOM_RUN_HEAVY_DOCS"]="true"` before execution to enable.

```julia
if get(ENV, "VSMARTMOM_RUN_HEAVY_DOCS", "false") == "true"
    fig = Figure(size=(600, 700))
    ax1 = Axis(fig[1,1], yscale=log10, ylabel="f₁₁", xlabel="cos(Θ)")
    ax2 = Axis(fig[2,1], ylabel="f₁₂/f₁₁", xlabel="cos(Θ)")
    record(fig, joinpath(@__DIR__, "scattering_radius_sweep.gif"), 0.03:0.10:4.3; framerate=5) do r
        empty!(ax1); empty!(ax2)
        local size_distribution_i = LogNormal(log(r), log(σ))
        local aero_i = Aerosol(size_distribution_i, nᵣ, nᵢ)
        local model_i = make_mie_model(NAI2(), aero_i, λ, polarization_type, truncation_type, r_max, nquad_radius)
        local optics_i = compute_aerosol_optical_properties(model_i)
        local smat_i = reconstruct_phase(optics_i.greek_coefs, μ_quad)
        lines!(ax1, μ_quad, smat_i.f₁₁)
        ylims!(ax1, 1e-3, 1e3)
        ax1.title = "f₁₁, r = $(round(r, digits=2)) μm"
        lines!(ax2, μ_quad, smat_i.f₁₂ ./ smat_i.f₁₁)
        ylims!(ax2, -1.1, 1.1)
    end
end
```

When run locally with `VSMARTMOM_RUN_HEAVY_DOCS=true`, this writes
`scattering_radius_sweep.gif` next to the tutorial.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

