#=
 
This file contains RT interaction-related functions
 
=#

"""
    interaction_helper!(::ScatteringInterface_00, SFI, composite_layer, added_layer, I_static)

Non-scattering interaction: neither the composite nor the added layer
scatters.  Transmission matrices are simply multiplied and source terms
are propagated without any reflectance coupling.
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

Interaction when only the added (homogeneous) layer scatters.  The composite
layer above is non-scattering, so `R⁻⁺` is built from the added layer's `r⁻⁺`
transported through the composite transmission.
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

Interaction when only the composite layer scatters and the added layer is
non-scattering.  The composite `R⁺⁻` is transported through the added
layer's transmission without requiring a matrix inversion.
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

Full-scattering interaction: both the composite and added layers scatter.
Requires two batched matrix inversions — `(I − r⁻⁺R⁺⁻)⁻¹` and
`(I − R⁺⁻r⁻⁺)⁻¹` — to account for the geometric series of inter-layer
reflections.  Updates all six composite-layer fields (R⁻⁺, R⁺⁻, T⁺⁺,
T⁻⁻, J₀⁺, J₀⁻) in-place.
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