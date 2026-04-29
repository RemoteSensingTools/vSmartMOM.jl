# Add a Raman Mode

**For:** method developers extending inelastic-scattering support.

**Next:** [Inelastic Scattering](../Inelastic/Overview.md), [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [Library](../api_reference.md).

Raman support is selected by dispatch on `AbstractRamanType`. Adding a new mode should mean adding methods for that mode, not adding new `isa` branches inside the main RT loop.

This guide covers Raman/Cabannes extension points only. Solar-induced fluorescence remains a product/data-policy decision and should not get a public extension guide until its fixtures and supported workflows are settled.

## Current Dispatch Layers

| Layer | File | What to extend |
| --- | --- | --- |
| Mode state | `src/Inelastic/types.jl` | New `AbstractRamanType` subtype and trait methods |
| Molecular properties | `src/Inelastic/raman_atmo_prop.jl` or `raman_stellar_prop.jl` | `getRamanSSProp!` methods that fill mode state |
| Rayleigh/Cabannes source choice | `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl` | `_rayleigh_greek_source` if the elastic phase source differs |
| Layer storage | `src/CoreRT/tools/rt_helper_functions.jl` | `make_added_layer` and `make_composite_layer` for required tensors |
| Elemental source | `src/CoreRT/CoreKernel/elemental_inelastic*.jl` | source-term construction for redistribution |
| Core loop | `src/CoreRT/CoreKernel/rt_kernel*.jl` | `rt_kernel!` dispatch for the new mode |
| Doubling/interaction | `src/CoreRT/CoreKernel/doubling_inelastic.jl`, `interaction_inelastic.jl` | propagation of inelastic matrices and source vectors |

## Define The Mode

Start with a concrete mode type in `src/Inelastic/types.jl`.

```julia
Base.@kwdef mutable struct MyRaman{FT<:AbstractFloat} <: AbstractRamanType
    n2::InelasticScattering.MolecularConstants{Float64}
    o2::InelasticScattering.MolecularConstants{Float64}
    greek_raman::GreekCoefs
    weights::Vector{FT}
    indices::Vector{Int}
    F₀::Matrix{FT}
    n_Raman::Int
end
```

Molecular constants should normally stay `Float64`, even for `MyRaman{Float32}`. Raman prefactors can be near the `Float32` subnormal range; the current code deliberately widens molecular state while keeping RT arrays parameterized by `FT`.

## Implement Traits

The default trait methods assume an inelastic mode:

```julia
has_inelastic(::AbstractRamanType) = true
uses_cabannes_phase(rs::AbstractRamanType) = has_inelastic(rs)
needs_interaction_workspace(rs::AbstractRamanType) = has_inelastic(rs)
needs_rayleigh_expansion(rs::AbstractRamanType) = has_inelastic(rs)
```

Override only when the new mode does something different. For example, `noRS` and `noRS_plus` return `false` for `has_inelastic`, which keeps the elastic path on Rayleigh coefficients and avoids inelastic workspaces.

If the mode has redistribution weights that require normalization, add a specific `normalize_raman_weights!(rs::MyRaman, model, iBand)` method.

## Fill Mode State

Add a `getRamanSSProp!` method that populates all arrays consumed by the RT kernels. Existing `RRS`, `VS_0to1`, `VS_1to0`, and `_plus` methods are the reference.

The method usually computes:

- Cabannes/Rayleigh fractions;
- redistribution weights from incident to scattered spectral indices;
- `Z⁺⁺` and `Z⁻⁺` Raman phase matrices;
- `n_Raman`, reference indices, and any band limits needed by plus-mode paths.

Keep the precomputation separate from `rt_run`; `rt_run` should receive a configured Raman mode and then let dispatch carry it through the solver.

## Allocate The Right Layer Type

If the mode needs inelastic matrices, extend:

```julia
make_added_layer(rs::MyRaman, FT, arr_type, dims, nSpec)
make_composite_layer(rs::MyRaman, FT, arr_type, dims, nSpec)
```

The existing Raman modes use `AddedLayerRS` and `CompositeLayerRS`, which carry elastic matrices plus inelastic redistribution tensors and source arrays. Do not allocate these for modes that can use the elastic `AddedLayer`; memory scales quickly with `nSpec * n_Raman`.

## Add Kernel Dispatch

Add the smallest set of `rt_kernel!`, `rt_kernel_ss!`, and multisensor methods required by the mode. The method should follow the existing sequence:

1. build inelastic elemental terms;
2. build elastic elemental terms if the mode still has an elastic component;
3. run the matching inelastic doubling routine;
4. copy into or interact with the composite layer.

If the mode only supports forward RT, make that explicit. The current linearized core has a mature elastic path; inelastic linearization is not yet a general public extension surface.

## Add Tests

At minimum:

- trait tests for `has_inelastic`, Cabannes selection, and workspace needs;
- constructor/type-stability tests for `Float32` and `Float64`;
- a tiny `getRamanSSProp!` smoke test that checks dimensions and finite values;
- a CPU forward smoke test against `noRS` for the same scene;
- a focused regression test if a reference spectrum exists.

Use `test/test_type_stability.jl`, `test/test_forward_raman.jl`, `test/test_forward_ss.jl`, and `test/test_forward_raman_phase1b.jl` as current patterns. GPU tests should remain conditional on CUDA availability.
