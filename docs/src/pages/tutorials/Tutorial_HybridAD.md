```@meta
EditURL = "Tutorial_HybridAD.jl"
```

# Hybrid AD: Combining ForwardDiff with Analytic RT

vSmartMOM uses a **hybrid approach** to compute Jacobians efficiently:
the outer derivatives (optical properties w.r.t. state vector parameters)
can use ForwardDiff.jl, while the inner RT solver uses fully analytic
derivatives for performance. This tutorial explains the architecture
and shows how to use both paths on CPU and GPU.

## 1) The AD Boundary

The key architectural idea is a clean separation between two zones:

```
State Vector x = [albedo, τ_aer, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p, VMR_gas, ...]
     │
     │  ← AD zone: ForwardDiff.Dual allowed here
     ▼
Optical Properties: τ(λ,z), ω(λ,z), Z(λ,z) per layer
     │
     │  ← extract values; store ∂(τ,ω,Z)/∂x as arrays
     ▼
RT Kernels: pure Float64/Float32 [analytic ∂R/∂(τ,ω,Z)]
     │
     │  ← chain rule: ∂R/∂x = ∂R/∂(τ,ω,Z) · ∂(τ,ω,Z)/∂x
     ▼
Output Jacobians: ∂R/∂x
```

The RT kernels (elemental, doubling, interaction) never see Dual numbers.
They compute analytic derivatives with respect to 3 core optical
properties (τ, ϖ, Z) using the linearized adding-doubling method.
This is fast, GPU-compatible, and avoids the overhead of tracking
Dual numbers through deep batched matrix operations.

## 2) Linearized RT on CPU

Let's run the linearized model and inspect the Jacobian layout.

```julia
using vSmartMOM
using vSmartMOM.CoreRT

params = parameters_from_yaml(
    joinpath(dirname(dirname(pathof(vSmartMOM))),
             "test", "test_parameters", "JacobianTestFast.yaml"))
params.architecture = vSmartMOM.Architectures.CPU()
```

Build the linearized model. `LinMode()` tells the constructor to
compute derivative containers alongside the forward model.

```julia
model, lin_model = model_from_parameters(LinMode(), params)
```

Run linearized RT: returns (R, T, dR, dT)

```julia
NAer  = length(params.scattering_params.rt_aerosols)
NGas  = size(lin_model.τ̇_abs[1], 1)
NSurf = 1
R_cpu, T_cpu, dR_cpu, dT_cpu = rt_run(model, lin_model, NAer, NGas, NSurf)

println("R  shape: ", size(R_cpu), "  (nVZA × nStokes × nSpec)")
println("dR shape: ", size(dR_cpu), "  (nParams × nVZA × nStokes × nSpec)")
```

## 3) Understanding ParameterLayout

The Jacobian dimension of `dR` is ordered by `ParameterLayout`:
each aerosol gets 7 sub-parameters, followed by gas VMRs, then surface.

```julia
layout = CoreRT.ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                                 n_gases=NGas, n_surface=NSurf)

println("\nParameter layout:")
println("  Total parameters: ", CoreRT.n_total(layout))
for ia in 1:NAer
    rng = CoreRT.aerosol_range(layout, ia)
    println("  Aerosol $ia: indices $rng  [τ_ref, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p]")
end
if NGas > 0
    println("  Gas VMRs:  indices ", CoreRT.gas_range(layout))
end
println("  Surface:   indices ", CoreRT.surface_range(layout))
```

Extract specific Jacobians (first VZA, Stokes-I):

```julia
idx_albedo = CoreRT.surface_index(layout, 1)
idx_tau    = CoreRT.aerosol_range(layout, 1)[1]  # τ_ref
idx_nr     = CoreRT.aerosol_range(layout, 1)[2]  # nᵣ

println("\ndR/d(albedo) first 3 spec points: ",
        round.(dR_cpu[idx_albedo, 1, 1, 1:min(3,end)], digits=6))
println("dR/d(τ_ref)  first 3 spec points: ",
        round.(dR_cpu[idx_tau, 1, 1, 1:min(3,end)], digits=6))
println("dR/d(nᵣ)     first 3 spec points: ",
        round.(dR_cpu[idx_nr, 1, 1, 1:min(3,end)], digits=6))
```

## 4) Linearized RT on GPU

The same code runs on GPU by switching the architecture. All arrays
are automatically moved to `CuArray` and the RT kernels use
`KernelAbstractions.jl` for device-portable dispatch.

```julia
if vSmartMOM.Architectures.has_cuda()
    params_gpu = parameters_from_yaml(
        joinpath(dirname(dirname(pathof(vSmartMOM))),
                 "test", "test_parameters", "JacobianTestFast.yaml"))
    params_gpu.architecture = vSmartMOM.Architectures.GPU()
    model_gpu, lin_model_gpu = model_from_parameters(LinMode(), params_gpu)
    R_gpu, T_gpu, dR_gpu, dT_gpu = rt_run(model_gpu, lin_model_gpu,
                                           NAer, NGas, NSurf)

    R_gpu_a  = Array(R_gpu)
    dR_gpu_a = Array(dR_gpu)

    max_R_diff  = maximum(abs.(R_gpu_a .- R_cpu))
    max_dR_diff = maximum(abs.(dR_gpu_a .- dR_cpu))
    println("\nGPU vs CPU comparison:")
    println("  max |R_gpu - R_cpu|  = ", round(max_R_diff, sigdigits=3))
    println("  max |dR_gpu - dR_cpu| = ", round(max_dR_diff, sigdigits=3))
else
    println("\nCUDA not available — skipping GPU section.")
end
```

## 5) ForwardDiff for Mie Optical Properties

The Mie code supports ForwardDiff through
`compute_aerosol_optical_properties(model; autodiff=true)`.
This wraps the NAI2 computation in `ForwardDiff.jacobian` with respect
to the 4 Mie-sensitive parameters [rₘ, σ, nᵣ, nᵢ].

```julia
using vSmartMOM.Scattering

rt_aer = params.scattering_params.rt_aerosols[1]
truncation_type = Scattering.δBGE{Float64}(params.l_trunc, params.Δ_angle)
mie_model = make_mie_model(params.scattering_params.decomp_type,
                            rt_aer.aerosol, params.scattering_params.λ_ref,
                            params.polarization_type,
                            truncation_type,
                            params.scattering_params.r_max,
                            params.scattering_params.nquad_radius)
```

Analytic derivatives (the default production path)

```julia
aer_optics_analytic, lin_aer_optics = compute_aerosol_optical_properties(
    LinMode(), mie_model, Float64)

println("\nAnalytic Mie derivatives:")
println("  ω̃  = ", round(aer_optics_analytic.ω̃, digits=6))
println("  dω̃/d[nᵣ,nᵢ,rₘ,σ] = ", round.(lin_aer_optics.ω̃̇, digits=6))
```

ForwardDiff derivatives (for validation / custom parameters)

```julia
aer_optics_ad = compute_aerosol_optical_properties(mie_model; autodiff=true)

println("\nForwardDiff Mie derivatives stored in AerosolOptics.derivs:")
println("  derivs shape: ", size(aer_optics_ad.derivs))
println("  (rows = flattened [α;β;γ;δ;ϵ;ζ;ω̃;k], cols = [rₘ,σ,nᵣ,nᵢ])")
```

## 6) The δ-M Truncation Chain Rule

The `RawAerosolJacobian` struct defines the AD boundary for aerosols.
It holds un-truncated derivatives:

```julia
struct RawAerosolJacobian{T}
    τ̇_aer::T   # ∂τ_aer/∂x
    ω̃̇::T       # ∂ω̃/∂x
    ḟᵗ::T       # ∂fᵗ/∂x   ← truncation factor derivative
    Ż⁺⁺::T      # ∂Z⁺⁺/∂x
    Ż⁻⁺::T      # ∂Z⁻⁺/∂x
end
```

The function `delta_m_truncation_lin` then applies the analytic δ-M
chain rule to transform these raw derivatives into the
`CoreScatteringOpticalPropertiesLin` that the RT kernel expects:

```math
τ̇_mod = (1 - fᵗω̃) · τ̇_aer - τ_aer · (fᵗ · ω̃̇ + ω̃ · ḟᵗ)
```

```math
ϖ̇_mod = [ω̃̇ · (1 - fᵗ) - ḟᵗ · ω̃(1 - ω̃)] / (1 - fᵗω̃)²
```

This separation means you can compute `RawAerosolJacobian` via
**any method** (analytic Mie, ForwardDiff, or finite differences)
and the downstream RT chain rule handles the rest identically.

The current production code in `createAero` (compEffectiveLayerProperties_lin.jl)
uses the analytic Mie path. The `delta_m_truncation_lin` function
in `delta_m_truncation.jl` provides a cleaner, composable interface
for new code.

## Summary

The hybrid AD architecture gives you flexibility:

| Component | Method | Why |
|:----------|:-------|:----|
| RT kernels | Analytic (linearized adding-doubling) | Performance, GPU, batched matmul |
| Mie optics | Analytic (keep for speed) or ForwardDiff (for validation) | Mie series is AD-hostile |
| Surface albedo | Analytic (trivial) | Simple scalar derivative |
| Surface pressure | ForwardDiff (planned) | Affects many paths simultaneously |
| Gas VMR | Analytic (trivial: ∂τ_abs/∂VMR = cross_section) | Simple scaling |
| Aerosol profiles | Analytic | Already validated in atmo_prof_lin.jl |

The `ParameterLayout` struct generalizes the index arithmetic, and
the `RawAerosolJacobian` / `delta_m_truncation_lin` pair defines
a clean boundary where any upstream derivative method can plug in.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

