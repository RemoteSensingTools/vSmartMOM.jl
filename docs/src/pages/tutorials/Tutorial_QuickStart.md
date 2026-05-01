```@meta
EditURL = "Tutorial_QuickStart.jl"
```

# Quick Start: vSmartMOM Radiative Transfer

vSmartMOM is a polarized radiative transfer solver using the adding-doubling method.
This 5-minute tutorial runs a small forward RT simulation on CPU.

## 1) Load the package

```julia
using vSmartMOM
```

## 2) Load parameters from YAML

Parameters define the atmosphere, surface, geometry, and spectral bands.
See the IO Schema documentation for the full YAML format.

```julia
yaml_path = joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml")
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

The documentation build renders the same two arrays as an interactive Plotly
view. Hovering the bars shows the numeric Stokes-I values returned above.

```@raw html
<iframe title="Quickstart RT output" src="../../assets/plots/quickstart_rt_response.html" loading="lazy" style="width: 100%; height: 420px; border: 1px solid var(--vp-c-divider); border-radius: 8px;"></iframe>
```

## 5) Next steps

The shipped quickstart scene intentionally avoids absorption and aerosol Mie
setup. For richer spectra and plots, use the absorption, scattering, surface,
and Jacobian tutorials. For analytic derivatives, start with the Compute
Jacobians page.

## 6) GPU support

For NVIDIA GPUs with CUDA.jl, set `params.architecture = vSmartMOM.Architectures.GPU()`
before building the model. The same `model_from_parameters` and `rt_run` calls apply.
See the GPU tutorial for more.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

