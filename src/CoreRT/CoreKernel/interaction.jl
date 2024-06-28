#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper!(::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 

    # Source Function
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ J₀⁺
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ j₀⁻

    # Batched multiplication between added and composite
    T⁻⁻  .= t⁻⁻ ⊠ T⁻⁻
    T⁺⁺  .= t⁺⁺ ⊠ T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 

    # Source Function
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ J₀⁺         

    # Batched multiplication between added and composite
    R⁻⁺ .= T⁻⁻ ⊠ r⁻⁺ ⊠ T⁺⁺
    R⁺⁻ .= r⁺⁻
    T⁺⁺ .= t⁺⁺ ⊠ T⁺⁺
    T⁻⁻ .= T⁻⁻ ⊠ t⁻⁻    
end

# Scattering in inhomogeneous composite layer.
# no scattering in homogeneous layer which is 
# added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 

    # Source Function
    J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
    J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ j₀⁻    
   
    # Batched multiplication between added and composite
    T⁺⁺ .= t⁺⁺ ⊠ T⁺⁺
    T⁻⁻ .= T⁻⁻ ⊠ t⁻⁻
    R⁺⁻ .= t⁺⁺ ⊠ R⁺⁻ ⊠ t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻,temp1, temp2, temp1_ptr,temp2_ptr = added_layer     #these are aliases to the respective struct elements  
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer #these are aliases to the respective struct elements 
    
    # X₂₁ refers to added layer, X₁₀ to composite layer!

    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    #tmp_inv = similar(t⁺⁺)
    temp2 .= I_static .- r⁻⁺ ⊠ R⁺⁻
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr) 
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

"Compute interaction between composite and added layers"
function interaction!(scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
end