```@meta
EditURL = "Tutorial_QuickStart.jl"
```

# Quick Start: vSmartMOM Radiative Transfer

vSmartMOM is a polarized radiative transfer solver using the adding-doubling method.
This 5-minute tutorial runs a forward RT simulation and computes Jacobians.

## 1) Load the package

```julia
using vSmartMOM
using CairoMakie
```

## 2) Load parameters from YAML

Parameters define the atmosphere, surface, geometry, and spectral bands.
See the IO Schema documentation for the full YAML format.

```julia
yaml_path = joinpath(pkgdir(vSmartMOM),
                     "test", "test_parameters", "ParamsEMIT_fast.yaml")
params = read_parameters(yaml_path)
params.architecture = vSmartMOM.Architectures.CPU()
```

## 3) Build the model and run forward RT

```julia
model = model_from_parameters(params)
R, T = rt_run(model)
```

## 4) Interpret results

- **R** = reflectance at top-of-atmosphere (TOA) `[nVZA × nStokes × nSpec]`
- **T** = transmittance at bottom-of-atmosphere (BOA)

Dimensions: view zenith angles × Stokes components (I,Q,U,V) × spectral points.

```julia
println("R shape: ", size(R))
println("T shape: ", size(T))
println("R(nadir, I, λ₁) = ", R[1, 1, 1])
```

Plot the reflectance spectrum at nadir (Stokes-I):

```julia
fig = Figure(size=(700, 400))
ax = Axis(fig[1,1], xlabel="Spectral index", ylabel="TOA Reflectance (Stokes I)")
for ivza in 1:min(size(R, 1), 4)
    lines!(ax, R[ivza, 1, :], label="VZA #$(ivza)")
end
axislegend(ax, position=:rt)
fig
```

## 5) Linearized RT for Jacobians

Use `model_from_parameters(LinMode(), params)` and `rt_run` to get analytic
derivatives of R and T with respect to atmospheric and surface parameters.
No finite differences required — Jacobians come from a single RT pass.

```julia
using vSmartMOM.CoreRT

model_lin, lin_model = model_from_parameters(LinMode(), params)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R_lin, T_lin, dR, dT = rt_run(model_lin, lin_model, NAer, NGas, NSurf)
```

`dR` has shape `[nVZA × nStokes × nSpec × nParams]` — derivatives w.r.t. aerosols,
gas VMRs, and surface albedo.

```julia
println("dR shape (Jacobian): ", size(dR))
```

Visualize the Jacobian as a heatmap (spectral × parameter):

```julia
fig = Figure(size=(700, 450))
ax = Axis(fig[1,1],
    xlabel = "Parameter index",
    ylabel = "Spectral index",
    title  = "dR/dx at nadir (Stokes I)")
hm = heatmap!(ax, dR[1, 1, :, :]')
Colorbar(fig[1,2], hm, label="dR/dx")
fig
```

## 6) GPU support

For NVIDIA GPUs with CUDA.jl, set `params.architecture = vSmartMOM.Architectures.GPU()`
before building the model. The same `model_from_parameters` and `rt_run` calls apply.
See the GPU tutorial for more.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

