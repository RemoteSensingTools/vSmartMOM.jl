#=

This file contains linearized RT interaction (Adding Method) functions —
the tangent-linear partners of `interaction.jl`.

The interaction step combines two adjacent layers — the "composite" layer
(accumulated from above) and the newly "added" layer below — into a new
composite layer. The forward kernel implements Sanghavi et al. (2014,
JQSRT 133:412–433), Eqs. (23)–(28); these `_lin` companions implement the
tangent-linear partners from Sanghavi 2014 App. C, Eqs. (C.11)–(C.16).

Multiple dispatch on `ScatteringInterface` types selects the appropriate
case (see also `interaction.jl` for the forward methods):

- `ScatteringInterface_00`: Neither layer scatters (pure absorption cascade).
- `ScatteringInterface_01`: Only the added (lower) layer scatters.
- `ScatteringInterface_10`: Only the composite (upper) layer scatters.
- `ScatteringInterface_11`: Both layers scatter (general case; two batched
  matrix inversions; full Sanghavi 2014 Eqs. 23–28).

**Forward adding equations (`ScatteringInterface_11`):**

```math
\\mathbf{G} = (\\mathbf{E} - \\mathbf{R}^{-+}_\\text{add} \\, \\mathbf{R}^{+-}_\\text{comp})^{-1}
```
```math
\\mathbf{T}^{++}_\\text{new} = \\mathbf{T}^{++}_\\text{add} \\, \\mathbf{G} \\, \\mathbf{T}^{++}_\\text{comp}
```
```math
\\mathbf{R}^{-+}_\\text{new} = \\mathbf{R}^{-+}_\\text{comp} +
  \\mathbf{T}^{--}_\\text{comp} \\, \\mathbf{G} \\, \\mathbf{R}^{-+}_\\text{add} \\, \\mathbf{T}^{++}_\\text{comp}
```

**Linearized versions (Sanghavi 2014 Eqs. C.11–C.16):** propagate the full
``N_\\text{params}``-axis of derivatives through the same matrix algebra
using the product rule, with the closed-form
``\\dot{\\mathbf{G}} = \\mathbf{G}\\,\\dot{(\\mathbf{R}\\mathbf{R})}\\,\\mathbf{G}``
identity. No AD through the batched matrix inversion; the chain rule on
``(\\mathbf{E} - \\mathbf{R}\\mathbf{R})^{-1}`` is closed-form, which keeps
the Jacobian numerically stable and roughly forward-cost.

The ``\\mathbf{ap\\_*}`` fields on `AddedLayerLin` carry the chain-rule
expansion to the **physical state vector** ``\\mathbf{x}`` (computed by
[`lin_added_layer_all_params!`](@ref) before this kernel runs); the
``\\dot{\\mathbf{*}}`` fields carry the three-core derivatives w.r.t.
``(\\tau, \\varpi, \\mathbf{Z})``.

See [`Concepts/06 — Linearization`](../../docs/src/pages/concepts/06_linearization.md)
for the three-tier Jacobian diagram, the parameter-strategy table, and the
ParameterLayout column ordering.

=#

# No scattering in either the added layer or the composite layer
"""
    interaction_helper!(::ScatteringInterface_00, ...)

**Non-scattering interaction.** Both the composite and added layers are purely absorbing.
Transmission matrices multiply directly (no reflection), and the source functions are
attenuated by the transmission of the adjacent layer. Derivatives propagate linearly.
"""
function interaction_helper!(::ScatteringInterface_00, SFI,
                                #computed_layer_properties,
                                #computed_layer_properties_lin,
                                composite_layer::CompositeLayer{FT},
                                composite_layer_lin::CompositeLayerLin{FT},
                                added_layer::AddedLayer{FT},
                                added_layer_lin::AddedLayerLin{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[end]

    # If SFI, interact source function in no scattering
    if SFI
        for iparam=1:Nparams
            composite_layer_lin.J̇₀⁺[:,:,:,iparam] .= added_layer_lin.ap_J̇₀⁺[:,:,:,iparam] .+
                added_layer.t⁺⁺ ⊠ composite_layer_lin.J̇₀⁺[:,:,:,iparam] .+
                added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.J₀⁺
            composite_layer_lin.J̇₀⁻[:,:,:,iparam] .= composite_layer_lin.J̇₀⁻[:,:,:,iparam] .+
                composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_J̇₀⁻[:,:,:,iparam] .+
                composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠ added_layer.j₀⁻
        end
        composite_layer.J₀⁺ .= added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end
    # Batched multiplication between added and composite
    for iparam=1:Nparams
        composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] = added_layer_lin.ap_ṫ⁻⁻[:,:,:,iparam] ⊠ composite_layer.T⁻⁻ .+
                                added_layer.t⁻⁻ ⊠ composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam]
        composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam] = added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.T⁺⁺ .+
                                added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam]
    end
    composite_layer.T⁻⁻[:] = added_layer.t⁻⁻ ⊠ composite_layer.T⁻⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
end

"""
    interaction_helper!(::ScatteringInterface_01, ...)

**First-scattering interaction.** The composite layer (above) is non-scattering, but the
added layer (below) scatters. The composite layer's reflection matrices are initialized
from the added layer, while transmission is the product of both layers' transmissions.
"""
function interaction_helper!(::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT},
                                composite_layer_lin::CompositeLayerLin{FT},
                                added_layer::AddedLayer{FT},
                                added_layer_lin::AddedLayerLin{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[end]
    if SFI
        #J₀⁺, J₀⁻ = similar(composite_layer.J₀⁺), similar(composite_layer.J₀⁺)
        #J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻)
        #J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        for iparam=1:Nparams
            composite_layer_lin.J̇₀⁻[:,:,:,iparam] .= composite_layer_lin.J̇₀⁻[:,:,:,iparam] .+
                composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠
                (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻) .+
                composite_layer.T⁻⁻ ⊠
                (added_layer_lin.ap_ṙ⁻⁺[:,:,:,iparam] ⊠ composite_layer.J₀⁺ .+
                added_layer.r⁻⁺ ⊠ composite_layer_lin.J̇₀⁺[:,:,:,iparam] .+
                added_layer_lin.ap_J̇₀⁻[:,:,:,iparam])
            composite_layer_lin.J̇₀⁺[:,:,:,iparam] .= added_layer_lin.ap_J̇₀⁺[:,:,:,iparam] .+
                added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.J₀⁺ .+
                added_layer.t⁺⁺ ⊠ composite_layer_lin.J̇₀⁺[:,:,:,iparam]
        end
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+
            composite_layer.T⁻⁻ ⊠
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻)
        composite_layer.J₀⁺ .= added_layer.j₀⁺ .+
            added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
    end

    # Batched multiplication between added and composite
    for iparam = 1:Nparams
        composite_layer_lin.Ṙ⁻⁺[:,:,:,iparam] = composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṙ⁻⁺[:,:,:,iparam] ⊠ composite_layer.T⁺⁺ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam]
        composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] = added_layer_lin.ap_ṙ⁺⁻[:,:,:,iparam]
        composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam] = added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.T⁺⁺ .+
                                    added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam]
        composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] = composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠ added_layer.t⁻⁻ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[:,:,:,iparam]
    end
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
end

"""
    interaction_helper!(::ScatteringInterface_10, ...)

**Composite-only scattering interaction.** The composite layer (above) scatters but the
added layer (below) does not. Transmission and source functions are attenuated by the
non-scattering added layer, while the composite reflection is preserved.
"""
function interaction_helper!(::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT},
                                composite_layer_lin::CompositeLayerLin{FT},
                                added_layer::AddedLayer{FT},
                                added_layer_lin::AddedLayerLin{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[end]
    if SFI
        for iparam=1:Nparams
            composite_layer_lin.J̇₀⁺[:,:,:,iparam] .= added_layer_lin.ap_J̇₀⁺[:,:,:,iparam] .+
                added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠
                (composite_layer.J₀⁺ .+ composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻) .+
                added_layer.t⁺⁺ ⊠
                (composite_layer_lin.J̇₀⁺[:,:,:,iparam] .+
                composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] ⊠ added_layer.j₀⁻ .+
                composite_layer.R⁺⁻ ⊠ added_layer_lin.ap_J̇₀⁻[:,:,:,iparam])
            composite_layer_lin.J̇₀⁻[:,:,:,iparam] .= composite_layer_lin.J̇₀⁻[:,:,:,iparam] .+
                composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠ added_layer.j₀⁻ .+
                composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_J̇₀⁻[:,:,:,iparam]
        end
        composite_layer.J₀⁺ .= added_layer.j₀⁺ .+
            added_layer.t⁺⁺ ⊠
            (composite_layer.J₀⁺ .+ composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻)
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+
            composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    # Batched multiplication between added and composite
    for iparam=1:Nparams
        composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam] = added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.T⁺⁺ .+
                                        added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam]
        composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] = composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] ⊠ added_layer.t⁻⁻ .+
                                        composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[:,:,:,iparam]
        composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] = added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam] ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻ .+
                                            added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] ⊠ added_layer.t⁻⁻ .+
                                            added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[:,:,:,iparam]
    end
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

"""
    interaction_helper!(::ScatteringInterface_11, ...)

**Full scattering interaction (general case).** Both the composite and added layers scatter.
This is the most general case requiring the full Adding Method with geometric series
inversion for inter-layer multiple reflections:

```math
\\mathbf{G} = (\\mathbf{I} - \\mathbf{R}^{-+}_\\text{add} \\, \\mathbf{R}^{+-}_\\text{comp})^{-1}
```

The linearized version propagates `Nparams` derivatives through all matrix operations,
accounting for both layers' contributions to each physical parameter. The source function
interaction includes multiple-reflection terms between layers.

!!! note "Bug 17 fix"
    Corrected variable references: uses `J̇₀⁻, J̇₀⁺` from `composite_layer_lin`
    (not `ap_J̇₀⁻, ap_J̇₀⁺` from `added_layer_lin`) in the base terms of the
    source function interaction, avoiding double-counting of surface derivatives.
"""
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT},
                                composite_layer_lin::CompositeLayerLin{FT},
                                added_layer::AddedLayer{FT},
                                added_layer_lin::AddedLayerLin{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺) = added_layer #these are aliases to the respective struct elements
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) = composite_layer #these are aliases to the respective struct elements
    (; ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺, ap_J̇₀⁺, ap_J̇₀⁻) = added_layer_lin #these are aliases to the respective struct elements
    (; Ṙ⁻⁺, Ṙ⁺⁻, Ṫ⁺⁺, Ṫ⁻⁻, J̇₀⁺, J̇₀⁻) = composite_layer_lin #these are aliases to the respective struct elements

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[end]
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)
    tmp_inv_lin = similar(Ṫ⁺⁺)
    T01_inv_lin = similar(Ṫ⁺⁺)
    tmpṘ⁻⁺ = similar(Ṫ⁺⁺)
    tmpṪ⁻⁻ = similar(Ṫ⁺⁺)
    tmpap_J̇₀⁻ = similar(ap_J̇₀⁻)
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻)
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    # Hoist param-independent product
    r_times_T = r⁻⁺ ⊠ T⁺⁺
    @inbounds for iparam=1:Nparams
        tmp_inv_lin[:,:,:,iparam] .= tmp_inv ⊠ (ap_ṙ⁻⁺[:,:,:,iparam] ⊠ R⁺⁻ .+ r⁻⁺ ⊠ Ṙ⁺⁻[:,:,:,iparam]) ⊠ tmp_inv
        T01_inv_lin[:,:,:,iparam] .= Ṫ⁻⁻[:,:,:,iparam] ⊠ tmp_inv .+ T⁻⁻ ⊠ tmp_inv_lin[:,:,:,iparam]
        # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
        tmpṘ⁻⁺[:,:,:,iparam] .= Ṙ⁻⁺[:,:,:,iparam] .+
                        T01_inv_lin[:,:,:,iparam] ⊠ r_times_T .+
                        T01_inv ⊠ (ap_ṙ⁻⁺[:,:,:,iparam] ⊠ T⁺⁺ .+
                        r⁻⁺ ⊠ Ṫ⁺⁺[:,:,:,iparam])

        # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
        tmpṪ⁻⁻[:,:,:,iparam] .= T01_inv_lin[:,:,:,iparam] ⊠ t⁻⁻ .+ T01_inv ⊠ ap_ṫ⁻⁻[:,:,:,iparam]
    end

    if SFI
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        r_J0p_plus_j0m = r⁻⁺ ⊠ J₀⁺ .+ added_layer.j₀⁻
        tmpJ₀⁻ = J₀⁻ .+ T01_inv ⊠ r_J0p_plus_j0m
        @inbounds for iparam=1:Nparams
            #@show size(tmpap_J̇₀⁻), size(ap_J̇₀⁻)
            #@show size(T01_inv_lin), size(r⁻⁺)
            #@show size(J₀⁺), size(added_layer.j₀⁻)
            tmpap_J̇₀⁻[:,:,:,iparam] .= J̇₀⁻[:,:,:,iparam] .+
                T01_inv_lin[:,:,:,iparam] ⊠ r_J0p_plus_j0m .+
                T01_inv ⊠ (ap_ṙ⁻⁺[:,:,:,iparam] ⊠ J₀⁺ .+ r⁻⁺ ⊠ J̇₀⁺[:,:,:,iparam] .+ ap_J̇₀⁻[:,:,:,iparam])
        end
    end

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    tmpR⁻⁺ = R⁻⁺ .+ T01_inv ⊠ r_times_T

    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ = T01_inv ⊠ t⁻⁻

    # Repeating for mirror-reflected directions
    T21_inv_lin = similar(Ṫ⁺⁺)
    tmpṘ⁺⁻ = similar(Ṫ⁺⁺)
    tmpṪ⁺⁺ = similar(Ṫ⁺⁺)
    tmpap_J̇₀⁺ = similar(ap_J̇₀⁺)
    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺)
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    # Hoist param-independent product
    R_times_t = R⁺⁻ ⊠ t⁻⁻
    @inbounds for iparam=1:Nparams
        tmp_inv_lin[:,:,:,iparam] .= tmp_inv ⊠ (R⁺⁻ ⊠ ap_ṙ⁻⁺[:,:,:,iparam] .+ Ṙ⁺⁻[:,:,:,iparam] ⊠ r⁻⁺) ⊠ tmp_inv
        T21_inv_lin[:,:,:,iparam] .= ap_ṫ⁺⁺[:,:,:,iparam] ⊠ tmp_inv .+ t⁺⁺ ⊠ tmp_inv_lin[:,:,:,iparam]

        # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
        tmpṪ⁺⁺[:,:,:,iparam] .= T21_inv_lin[:,:,:,iparam] ⊠ T⁺⁺ .+ T21_inv ⊠ Ṫ⁺⁺[:,:,:,iparam]

        # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
        tmpṘ⁺⁻[:,:,:,iparam] .= ap_ṙ⁺⁻[:,:,:,iparam] .+ T21_inv_lin[:,:,:,iparam] ⊠ R_times_t .+
                                    T21_inv ⊠ (Ṙ⁺⁻[:,:,:,iparam] ⊠ t⁻⁻ .+ R⁺⁻ ⊠ ap_ṫ⁻⁻[:,:,:,iparam])
    end
    if SFI
        J0p_plus_R_j0m = J₀⁺ .+ R⁺⁻ ⊠ added_layer.j₀⁻
        @inbounds for iparam=1:Nparams
            tmpap_J̇₀⁺[:,:,:,iparam] .= added_layer_lin.ap_J̇₀⁺[:,:,:,iparam] .+
                T21_inv_lin[:,:,:,iparam] ⊠ J0p_plus_R_j0m .+
                T21_inv ⊠ (J̇₀⁺[:,:,:,iparam] .+
                    Ṙ⁺⁻[:,:,:,iparam] ⊠ added_layer.j₀⁻ .+
                    R⁺⁻ ⊠ added_layer_lin.ap_J̇₀⁻[:,:,:,iparam])
        end
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        tmpJ₀⁺ = added_layer.j₀⁺ .+ T21_inv ⊠
            J0p_plus_R_j0m
    end

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺ = T21_inv  ⊠ T⁺⁺

    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻ = r⁺⁻ .+ T21_inv ⊠ R_times_t

    if SFI
        composite_layer.J₀⁺[:] = tmpJ₀⁺
        composite_layer.J₀⁻[:] = tmpJ₀⁻

        @inbounds for iparam=1:Nparams
            #@show size(tmpap_J̇₀⁺), size(composite_layer_lin.J̇₀⁺)
            #@show size(tmpap_J̇₀⁻), size(composite_layer_lin.J̇₀⁻)
            composite_layer_lin.J̇₀⁺[:,:,:,iparam] .= tmpap_J̇₀⁺[:,:,:,iparam]
            composite_layer_lin.J̇₀⁻[:,:,:,iparam] .= tmpap_J̇₀⁻[:,:,:,iparam]
        end
    end
    composite_layer.R⁺⁻[:] = tmpR⁺⁻
    composite_layer.T⁻⁻[:] = tmpT⁻⁻
    composite_layer.R⁻⁺[:] = tmpR⁻⁺
    composite_layer.T⁺⁺[:] = tmpT⁺⁺

    composite_layer_lin.Ṙ⁺⁻[:] = tmpṘ⁺⁻
    composite_layer_lin.Ṫ⁻⁻[:] = tmpṪ⁻⁻
    composite_layer_lin.Ṙ⁻⁺[:] = tmpṘ⁻⁺
    composite_layer_lin.Ṫ⁺⁺[:] = tmpṪ⁺⁺
end

"""
    interaction!(scattering_interface, SFI, composite_layer, composite_layer_lin,
                 added_layer, added_layer_lin, I_static)

Linearized Adding Method: combine composite and added layers, propagating N parameter derivatives.

Dispatches to `interaction_helper!` based on `ScatteringInterface` (00, 01, 10, or 11).
The forward matrices (R, T, J₀) and their derivatives (Ṙ, Ṫ, J̇₀) are updated in-place
in `composite_layer` and `composite_layer_lin`. Uses the Adding equations of
de Haan, Bosma & Hovenier (1987).
"""
function interaction!(scattering_interface::AbstractScatteringInterface,
        SFI,
        composite_layer::CompositeLayer{FT},
        composite_layer_lin::CompositeLayerLin{FT},
        added_layer::AddedLayer{FT},
        added_layer_lin::AddedLayerLin{FT},
        I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    #@show A1[1,1,1], A2[1,1,1]
    interaction_helper!(scattering_interface, SFI,
        composite_layer, composite_layer_lin,
        added_layer, added_layer_lin, I_static)
    #A1 = Array(composite_layer.J₀⁻)
    #A2 = Array(composite_layer.J₀⁺)
    #@show A1[1,1,1], A2[1,1,1]
    synchronize_if_gpu()
end
