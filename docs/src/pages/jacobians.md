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

result = rt_run_lin(model, lin_model, NAer, NGas, NSurf)
R, T, dR, dT = result
```

`result` is an [`ObserverRTResultLin`](@ref). Its named endpoint fields are:

- `result.toa`: reflected Stokes field, shaped `nVZA × nStokes × nSpec`
- `result.boa`: transmitted Stokes field, shaped like `result.toa`
- `result.toa_jacobian`: derivative of `result.toa`, shaped
  `nVZA × nStokes × nSpec × nParams`
- `result.boa_jacobian`: derivative of `result.boa`, shaped like the TOA
  Jacobian

An endpoint field and its Jacobian are `nothing` when that endpoint was not
selected by `obs_alt`. Historical four-value destructuring remains available,
as shown above.

## Interior-Height Jacobians

Strict-interior observer heights are supported by the analytic elastic
multiple-scatter solver. Configure them through `obs_alt` before building the
forward and linearized models:

```julia
height_params = deepcopy(params)
height_params.obs_alt = [0.0, 5.0]  # endpoints plus one interior interface

height_model, height_lin_model = model_from_parameters_lin(height_params)
height_result = rt_run_lin(
    height_model, height_lin_model, NAer, NGas, NSurf)

level = only(height_result.levels)
level.height_km
level.boundary_index

level.upwelling
level.downwelling
level.unscattered_downwelling

level.upwelling_jacobian
level.downwelling_jacobian
level.unscattered_downwelling_jacobian

total_down = total_downwelling(level)
total_down_jacobian = total_downwelling_jacobian(level)
```

The radiance fields have shape `nVZA × nStokes × nSpec`; each Jacobian
adds `nParams` as its last axis. The parameter order is stored directly in
`height_result.layout`.

The unscattered solar contribution remains separate from the diffuse MOM
field. For fixed solar irradiance and geometry its atmospheric tangent is

```math
\frac{\partial L_{\mathrm{direct}}}{\partial x_j}
= -\frac{L_{\mathrm{direct}}}{\mu_0}
  \frac{\partial \tau_{\mathrm{above}}}{\partial x_j}.
```

It is zero for non-solar output ordinates and for surface-parameter columns.
The BOA endpoint `result.boa` and `result.boa_jacobian` retain the historical
total-radiance convention and therefore already include the direct carrier;
the separate fields apply to strict-interior records.
Raman/inelastic linearization remains unsupported. The separate
single-scatter driver `rt_run_ss` continues to reject strict-interior observer
requests. Interior linearized runs currently accept `SolarBeam` or `NoSource`;
thermal and surface-emission sources, plus `CanopySurface`, are rejected
explicitly until their multisensor tangent/source-slot propagation is
implemented.

## Slice The Parameter Dimension

Use [`ParameterLayout`](@ref) rather than hard-coded offsets. The last axis of
every full-RT Jacobian has this stable block ordering:

| Block | Accessor | Columns within block | Derivative units |
|:--|:--|:--|:--|
| Surface pressure | `psurf_index(layout)` | one scalar `p_surf` | per hPa |
| Aerosols | `aerosol_range(layout, iaer)` | seven columns per aerosol | parameter-dependent |
| Trace gases | `gas_profile_range(layout, igas, Nz)` | `Nz` columns per gas, TOA to BOA | per unit VMR |
| Surface | `surface_range(layout)` | surface-model parameters | parameter-dependent |
| Canopy | `canopy_range(layout)` | canopy-model parameters, when present | parameter-dependent |

Thus the state vector is

```text
x = [p_surf,
     aerosol₁[1:7], ..., aerosol_NAer[1:7],
     VMR₁(z₁), ..., VMR₁(z_Nz),
     ...,
     VMR_NGas(z₁), ..., VMR_NGas(z_Nz),
     surface parameters,
     canopy parameters]
```

Here `z₁` is the top atmospheric layer and `z_Nz` is the bottom layer. Gas
VMR means dimensionless mole fraction; a derivative “per ppm” is therefore
`1e-6 .* dR_dVMR`. Pressure, temperature, dry-air column, and absorption
cross section are held fixed for a gas-VMR derivative. Surface-pressure
derivatives are a separate state column.

```julia
layout = result.layout
Nz = length(model.profile.p_full)

surface_idx = CoreRT.surface_index(layout)
gas_block = CoreRT.gas_range(layout)
gas1_profile = CoreRT.gas_profile_range(layout, 1, Nz)
gas1_bottom = CoreRT.gas_layer_index(layout, 1, Nz, Nz)
aerosol1_block = NAer > 0 ? CoreRT.aerosol_range(layout, 1) : 1:0

dR_surface = dR[:, 1, :, surface_idx]
dR_dgas1_profile = dR[:, :, :, gas1_profile]
dT_dgas1_bottom = dT[:, :, :, gas1_bottom]
```

`NGas = size(lin_model.τ̇_abs[1], 1)` is the number of flattened gas-layer
columns, not merely the number of molecular species. Consequently,
`NGas % Nz == 0`, and the number of gas slots is `NGas ÷ Nz`. Do not pass the
number of species to `rt_run_lin`; pass the flattened dimension exactly as
shown in the setup example.

Each aerosol contributes seven derivative slots in this order:

1. reference optical depth `τ_ref`;
2. real refractive index `nᵣ`;
3. imaginary refractive index `nᵢ`;
4. median radius `rₘ`;
5. geometric width `σ_g`;
6. vertical-profile location (`p₀` or `z₀`);
7. vertical-profile width (`σ_p` or `σ₀`).

The array dimensions should always satisfy:

```julia
@assert size(dR)[1:3] == size(R)
@assert size(dT)[1:3] == size(T)
@assert size(dR, 4) == CoreRT.n_total(layout)
@assert size(dT, 4) == CoreRT.n_total(layout)
@assert NGas % Nz == 0
```

### Parsing into named retrieval blocks

The following pattern avoids embedding any column offsets in downstream
retrieval or diagnostic code:

```julia
function parse_jacobian(J, model, layout, NAer, NGas)
    Nz = length(model.profile.p_full)
    @assert NGas % Nz == 0
    n_gas_species = NGas ÷ Nz

    aerosols = [@view J[:, :, :, CoreRT.aerosol_range(layout, ia)]
                for ia in 1:NAer]
    gases = [@view J[:, :, :, CoreRT.gas_profile_range(layout, ig, Nz)]
             for ig in 1:n_gas_species]

    return (
        psurf = @view(J[:, :, :, CoreRT.psurf_index(layout)]),
        aerosols = aerosols,
        gases = gases, # each entry: nVZA × nStokes × nSpec × Nz
        surface = @view(J[:, :, :, CoreRT.surface_range(layout)]),
    )
end

dR_named = parse_jacobian(dR, model, layout, NAer, NGas)
dT_named = parse_jacobian(dT, model, layout, NAer, NGas)
```

When applying a physical profile perturbation, contract over the final axis:

```julia
delta_vmr = fill(1e-6, Nz) # +1 ppm in every layer
delta_R = dropdims(sum(dR_named.gases[1] .* reshape(delta_vmr, 1, 1, 1, Nz),
                       dims=4), dims=4)
```

This contraction recovers a column-wide perturbation while retaining the
ability to apply arbitrary nonuniform vertical profiles.

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
