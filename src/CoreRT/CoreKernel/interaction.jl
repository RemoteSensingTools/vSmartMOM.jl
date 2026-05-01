#=
 
This file contains RT interaction-related functions
 
=#

"""
    interaction_helper!(::ScatteringInterface_00, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where neither the
composite layer above nor the newly added layer below scatters.

Both R-fields are zero by construction; the geometric-series factor
`(E − r⁻⁺ R⁺⁻)⁻¹` collapses to the identity. The full adding equations
reduce to:

    J₀⁺ ← j₀⁺ + t⁺⁺ · J₀⁺                  (downward source cascade)
    J₀⁻ ← J₀⁻ + T⁻⁻ · j₀⁻                  (upward source cascade)
    T⁻⁻ ← t⁻⁻ · T⁻⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺                        (Beer-law transmission only)

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case and `docs/src/pages/concepts/04_mom_solver.md` for the
prose walkthrough.
"""
function interaction_helper!(::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer     
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) = composite_layer 

    # Source Function
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ J₀⁺
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ j₀⁻

    # Batched multiplication between added and composite
    T⁻⁻  .= t⁻⁻ ⊠ T⁻⁻
    T⁺⁺  .= t⁺⁺ ⊠ T⁺⁺
end

"""
    interaction_helper!(::ScatteringInterface_01, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where only the
newly added (lower) layer scatters; the composite layer above is
non-scattering so `R⁺⁻ = 0` and the geometric-series factor
`(E − r⁻⁺ R⁺⁻)⁻¹` collapses to the identity.

The composite reflectance now picks up the added layer's `r⁻⁺` carried
through the composite transmission both ways:

    J₀⁻ ← J₀⁻ + T⁻⁻ · (r⁻⁺ J₀⁺ + j₀⁻)
    J₀⁺ ← j₀⁺ + t⁺⁺ · J₀⁺
    R⁻⁺ ← T⁻⁻ · r⁻⁺ · T⁺⁺                   (Eq. 24 with R₀₁ = 0)
    R⁺⁻ ← r⁺⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺
    T⁻⁻ ← T⁻⁻ · t⁻⁻

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case.
"""
function interaction_helper!(::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer     
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) = composite_layer 

    # Source Function
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ J₀⁺         

    # Batched multiplication between added and composite
    R⁻⁺ .= T⁻⁻ ⊠ r⁻⁺ ⊠ T⁺⁺
    R⁺⁻ .= r⁺⁻
    T⁺⁺ .= t⁺⁺ ⊠ T⁺⁺
    T⁻⁻ .= T⁻⁻ ⊠ t⁻⁻    
end

"""
    interaction_helper!(::ScatteringInterface_10, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where only the
composite (upper) layer scatters; the newly added layer below is
non-scattering so `r⁻⁺ = 0` and the geometric-series factor
`(E − R⁺⁻ r⁻⁺)⁻¹` collapses to the identity.

`R⁻⁺` is unchanged because nothing below it reflects; `R⁺⁻` is the existing
composite reflectance carried through the new layer's Beer-law transmission:

    J₀⁺ ← j₀⁺ + t⁺⁺ · (J₀⁺ + R⁺⁻ j₀⁻)
    J₀⁻ ← J₀⁻ + T⁻⁻ · j₀⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺
    T⁻⁻ ← T⁻⁻ · t⁻⁻
    R⁺⁻ ← t⁺⁺ · R⁺⁻ · t⁻⁻                  (Eq. 26 with R₂₁ = 0)

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case.
"""
function interaction_helper!(::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer     
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) = composite_layer 

    # Source Function
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ j₀⁻    
   
    # Batched multiplication between added and composite
    T⁺⁺ .= t⁺⁺ ⊠ T⁺⁺
    T⁻⁻ .= T⁻⁻ ⊠ t⁻⁻
    R⁺⁻ .= t⁺⁺ ⊠ R⁺⁻ ⊠ t⁻⁻
end

"""
    interaction_helper!(::ScatteringInterface_11, SFI, composite_layer, added_layer, I_static)

Full matrix-operator-method **adding** step: both the composite layer above
and the newly added layer below scatter. This is the algebraic heart of MOM
and the reason the method handles arbitrarily inhomogeneous atmospheres.

Implements Sanghavi et al. (2014, JQSRT 133:412–433), Eqs. (23)–(28),
restated for the layered atmosphere of Sanghavi & Frankenberg (2023, JQSRT
311:108791) Eq. (12). With the convention that the added (lower) layer is
labelled `21` and the composite (upper) layer is labelled `01`, two batched
matrix inversions form the geometric-series factors that capture infinitely
many reflections between the two layers:

    T01_inv = T⁻⁻ · (E − r⁻⁺ R⁺⁻)⁻¹     # Eq. (24): "T₀₁(I − R₂₁R₀₁)⁻¹"
    T21_inv = t⁺⁺ · (E − R⁺⁻ r⁻⁺)⁻¹     # Eq. (23): "T₂₁(I − R₀₁R₂₁)⁻¹"

Then the composite operators are updated **in place** following the
adding equations (lowercase = added layer; uppercase = composite):

    J₀⁻ ← J₀⁻ + T01_inv · (r⁻⁺ J₀⁺ + j₀⁻)        Eq. (28): J₀₂⁻
    R⁻⁺ ← R⁻⁺ + T01_inv · r⁻⁺ · T⁺⁺              Eq. (24): R₂₀
    T⁻⁻ ← T01_inv · t⁻⁻                          Eq. (25): T₀₂

    J₀⁺ ← j₀⁺ + T21_inv · (J₀⁺ + R⁺⁻ j₀⁻)        Eq. (27): J₂₀⁺
    T⁺⁺ ← T21_inv · T⁺⁺                          Eq. (23): T₂₀
    R⁺⁻ ← r⁺⁻ + T21_inv · R⁺⁻ · t⁻⁻              Eq. (26): R₀₂

The `(E − r⁻⁺ R⁺⁻)⁻¹` factor is the matrix geometric series capturing photon
paths that bounce between the two layers any number of times before
escaping. `batch_inv!` dispatches by array type to threaded BLAS (CPU),
CUBLAS strided (CUDA), or the portable KA LU kernel (Metal) — see
`docs/src/pages/concepts/07_architecture.md`.

This kernel runs once per atmospheric layer (after the TOA layer is copied
directly into the composite). The other three `ScatteringInterface_*`
methods drop terms that are zero by construction; this one is the full case.

# Concepts page
See [The MOM Solver § Adding / Interaction](../../docs/src/pages/concepts/04_mom_solver.md#adding--interaction)
for the prose walkthrough and the dispatch table that selects between the
four `ScatteringInterface_*` cases.
"""
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, temp1, temp2, temp1_ptr, temp2_ptr) = added_layer     #these are aliases to the respective struct elements  
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) = composite_layer #these are aliases to the respective struct elements 
    
    # X₂₁ refers to added layer, X₁₀ to composite layer!

    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    #tmp_inv = similar(t⁺⁺)
    temp2 .= I_static .- r⁻⁺ ⊠ R⁺⁻
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1 bla" batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr) 
    # Temporary arrays:
    
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ temp1;
    
    # J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
    J₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 
 
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    R⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺
    
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    T⁻⁻ .= T01_inv ⊠ t⁻⁻ 

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    #handle = CUBLAS.handle()
    #CUBLAS.math_mode!(handle, CUDA.FAST_MATH)
    #@show typeof(I_static .- R⁺⁻ ⊠ r⁻⁺)
    temp2 .= I_static .- R⁺⁻ ⊠ r⁻⁺
    @timeit "interaction inv2" batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr) 
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ temp1

    # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
    J₀⁺ .= j₀⁺ .+ T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    T⁺⁺ .= T21_inv  ⊠ T⁺⁺ 
    
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    R⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻  
end

"""
    interaction!(scattering_interface, SFI, composite_layer, added_layer, I_static)

Combine the accumulated [`CompositeLayer`](@ref) (from above) with a newly
doubled [`AddedLayer`](@ref) (below) using the adding equations.

Dispatches to the appropriate [`interaction_helper!`](@ref) based on the
`scattering_interface` type (`ScatteringInterface_00`, `_01`, `_10`, or
`_11`), then issues a GPU synchronisation barrier.
"""
function interaction!(scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
end