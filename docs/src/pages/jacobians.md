# Compute Jacobians

**For:** retrieval and inversion developers who need analytic sensitivities.

**Next:** [Linearization (Concepts)](concepts/06_linearization.md), [Library](api_reference.md), [Tutorial: Jacobians](tutorials/Tutorial_Jacobians.md).

vSmartMOM has a linearized radiative-transfer path for Jacobians. It propagates derivative optical properties through the same adding-doubling solver used for the forward model.

## Build A Linearized Model

Use `model_from_parameters_lin(params)` when the call site should make the Jacobian path explicit. It is a convenience alias for `model_from_parameters(LinMode(), params)`.

```julia
using vSmartMOM
using vSmartMOM.CoreRT

params = read_parameters(joinpath(pkgdir(vSmartMOM),
                                  "test", "test_parameters",
                                  "JacobianTestFast.yaml"))
params.architecture = vSmartMOM.Architectures.CPU()

model, lin_model = model_from_parameters_lin(params)
```

`model` is the forward model used by `rt_run(model)`. `lin_model` carries derivative optical properties, including gas absorption derivatives, aerosol optical-property derivatives, and layer optical-depth derivatives.

## Run Linearized RT

The linearized solver needs the number of aerosol, gas, and surface retrieval parameters so it can size the final Jacobian dimension.

```julia
NAer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R, T, dR, dT = rt_run_lin(model, lin_model, NAer, NGas, NSurf)
```

The outputs are:

- `R`: reflected Stokes field, shaped `nVZA × nStokes × nSpec`
- `T`: transmitted Stokes field, shaped like `R`
- `dR`: derivative of `R`, shaped `nVZA × nStokes × nSpec × nParams`
- `dT`: derivative of `T`, shaped like `dR`

## Slice The Parameter Dimension

Use [`ParameterLayout`](@ref) rather than hard-coded offsets.

```julia
layout = CoreRT.ParameterLayout(aerosol_params = 7,
                                n_aerosols = NAer,
                                n_gases = NGas,
                                n_surface = NSurf)

surface_idx = CoreRT.surface_index(layout)
gas_block = CoreRT.gas_range(layout)
aerosol1_block = NAer > 0 ? CoreRT.aerosol_range(layout, 1) : 1:0

dR_surface = dR[:, 1, :, surface_idx]
```

The current ordering is aerosol parameters first, then gas VMR parameters, then surface parameters. Each aerosol contributes seven derivative slots: `τ_ref`, `nᵣ`, `nᵢ`, mean radius, geometric width, vertical-location parameter, and vertical-width parameter.

For a longer walkthrough with plots and finite-difference checks, use [Tutorial: Jacobians](tutorials/Tutorial_Jacobians.md). The hybrid AD notes are in [Tutorial: HybridAD](tutorials/Tutorial_HybridAD.md).

## Standalone Exact Single-Scatter Seam

`StandaloneSS` also exposes a small Jacobian seam for exact first-order
single-scatter diagnostics. This is separate from the full MOM linearized
solver above. It is useful when validating path 1 atmospheric single scatter
or path 2 direct-beam surface reflection against first-order references.

The seam variables are:

- `τ_layer[iz, ispec]`
- `ϖ_eff[iz, ispec]`
- scalar phase `P_eff[igeom, iz, ispec]`
- vector phase `P_eff[igeom, istokes, iz, ispec]`
- scalar surface BRDF `surface_brdf[igeom, ispec]`
- vector surface BRDF `surface_brdf[igeom, istokes, ispec]`

Use `SSMeasurementSelector` to flatten only the retrieval rows you want.
By default it keeps all Stokes components; pass `stokes_indices = 1` for an
I-only retrieval vector.

```julia
using ForwardDiff
using vSmartMOM
using vSmartMOM.StandaloneSS

FT = Float64
geometry = SSGeometry(μ₀ = FT(0.79),
                      μv = FT[0.41, 0.73],
                      Δϕ = FT[0.2, 1.1])
absorber = AbsorptionSSContributor(τ = FT[0.06 0.02; 0.03 0.05])
n_water = Complex{FT}(FT(1.34), FT(1e-8))

surface_from_wind(wind) = CoxMunkSSSurface(
    wind_speed = wind, n_water = n_water,
    include_whitecaps = false, shadowing = true)

config_from_wind(wind) = ExactSSConfig(
    geometry = geometry,
    surface = surface_from_wind(wind),
    contributors = (absorber,),
    I0 = FT[1.0, 0.8],
    polarization_type = vSmartMOM.Scattering.Stokes_IQ{FT}())

selector = SSMeasurementSelector(paths = :path2, stokes_indices = 1:2)

f2 = run_exact_ss_with_jacobians(config_from_wind(FT(4.0));
                                 paths = :path2,
                                 selector = selector)

dρ_dwind = surface_brdf_wind_jacobian(config_from_wind(FT(4.0)))
J_wind = chain_rule_combine_surface_brdf(
    f2.jacobians.path2.surface_brdf,
    dρ_dwind,
    selector)
```

`J_wind` has shape `(nSelectedMeasurement, 1)`. The runnable version with a
finite-difference check is `examples/standalone_ss_vector_jacobian.jl`.
