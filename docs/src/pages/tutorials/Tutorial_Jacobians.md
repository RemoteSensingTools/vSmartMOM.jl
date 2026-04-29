```@meta
EditURL = "Tutorial_Jacobians.jl"
```

# Jacobians: Linearized Radiative Transfer

This tutorial shows how to compute analytic Jacobians (derivatives of
top-of-atmosphere radiance with respect to atmospheric and surface parameters)
using the linearized RT mode in vSmartMOM.

The linearized mode propagates derivatives through the adding-doubling
solver alongside the forward radiance, producing exact Jacobians in a
single RT pass — no finite differences needed.

This is a local walkthrough. It builds HITRAN absorption and aerosol Mie
optical properties, so it is intentionally heavier than `docs/test_examples.jl`.

## 1) Setup: load a test case

```julia
using vSmartMOM
using vSmartMOM.CoreRT
using CairoMakie

params = read_parameters(
    joinpath(pkgdir(vSmartMOM),
             "test", "test_parameters", "JacobianTestFast.yaml"))
params.architecture = vSmartMOM.Architectures.CPU()
```

## 2) Build the linearized model

`model_from_parameters_lin` tells the constructor to compute both the forward
model and the derivative containers (`RTModelLin`). It is equivalent to
`model_from_parameters(LinMode(), params)`.

```julia
model, lin_model = model_from_parameters_lin(params)
```

`model` is the same forward model as `model_from_parameters(params)`.
`lin_model` holds pre-computed ∂τ/∂x arrays for gases, aerosols, and
aerosol optics.

```julia
println("τ̇_abs bands: ", length(lin_model.τ̇_abs),
        "  shape[1]: ", size(lin_model.τ̇_abs[1]))
println("τ̇_aer bands: ", length(lin_model.τ̇_aer),
        "  shape[1]: ", size(lin_model.τ̇_aer[1]))
```

## 3) Run linearized RT

The linearized `rt_run` returns `(R, T, dR, dT)`:
- `R`  — reflected Stokes field  `[nVZA × nStokes × nSpec]`
- `T`  — transmitted field
- `dR` — Jacobian of R w.r.t. all parameters  `[nVZA × nStokes × nSpec × nParams]`
- `dT` — Jacobian of T

```julia
NAer  = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
NGas  = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R, T, dR, dT = rt_run_lin(model, lin_model, NAer, NGas, NSurf)

println("R  shape: ", size(R))
println("dR shape: ", size(dR))
```

Plot the forward reflectance spectrum:

```julia
fig = Figure(size=(700, 400))
ax = Axis(fig[1,1], xlabel="Spectral index", ylabel="TOA Reflectance (Stokes I)")
lines!(ax, R[1, 1, :], label="Nadir")
axislegend(ax, position=:rt)
fig
```

## 4) Interpret the Jacobian layout

The derivative dimension of `dR` is ordered by `ParameterLayout`:

| Index range                | Parameter                        |
|:---------------------------|:---------------------------------|
| `1 : NAer*7`               | Aerosol properties (7 per type)  |
| `NAer*7+1 : NAer*7+NGas`   | Gas VMRs (per variable molecule) |
| `NAer*7+NGas+1 : end`      | Surface albedo                   |

For this test case:

```julia
layout = CoreRT.ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                                n_gases=NGas, n_surface=NSurf)
Nparams = CoreRT.n_total(layout)
println("Total Jacobian parameters: ", Nparams,
        " (NAer×7=", NAer*7, ", NGas=", NGas, ", NSurf=", NSurf, ")")
```

## 5) Inspect individual derivatives

Surface albedo Jacobian (surface parameter, first VZA, Stokes-I):

```julia
surface_idx = CoreRT.surface_index(layout)
dR_albedo = dR[:, 1, :, surface_idx]
println("dR/d(albedo) at nadir, first 5 spectral points: ",
        round.(dR_albedo[1, 1:min(5,end)], digits=6))
```

Gas VMR Jacobians:

```julia
if NGas > 0
    igas_start = first(CoreRT.gas_range(layout))
    dR_gas1 = dR[:, 1, :, igas_start]
    println("dR/d(gas₁ VMR) at nadir, first 5 points: ",
            round.(dR_gas1[1, 1:min(5,end)], digits=6))
end
```

Plot the aerosol optical depth Jacobian (`τ_ref`, the first aerosol
parameter) rather than the full parameter matrix. This is the most useful
first diagnostic: it shows how an AOD perturbation changes the reflected
Stokes-I spectrum.

```julia
fig = Figure(size=(700, 500))
ax = Axis(fig[1,1],
    xlabel = "Spectral index",
    ylabel = "dR/dτ_ref (Stokes I)",
    title  = "AOD Jacobian")

if NAer > 0
    aod_idx = first(CoreRT.aerosol_range(layout, 1))
    lines!(ax, dR[1, 1, :, aod_idx], label="VZA 0°")
    if size(dR, 1) > 1
        lines!(ax, dR[2, 1, :, aod_idx], label="VZA 30°")
    end
end
axislegend(ax, position=:rt)
fig
```

The rendered docs use Plotly for the visible figure. The top panel shows the
forward reflectance; the bottom panel shows the `τ_ref` Jacobian for the same
viewing angles.

```@raw html
<iframe title="AOD Jacobian response" src="../../assets/plots/jacobian_aod_response.html" loading="lazy" style="width: 100%; height: 560px; border: 1px solid var(--vp-c-divider); border-radius: 8px;"></iframe>
```

## 6) Verify against finite differences (optional)

For validation, one can perturb a parameter by ε and compare:

```julia
ε = 1e-4
params_pert = deepcopy(params)
# e.g. perturb surface albedo
# ... modify params_pert ...
model_pert = model_from_parameters(params_pert)
R_pert, = rt_run(model_pert)
K_fd = (R_pert .- R) ./ ε
# Compare K_fd with dR[:, :, :, end]
```

The test suite in `test/test_forward_lin.jl` performs this comparison
systematically for all parameter types.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

