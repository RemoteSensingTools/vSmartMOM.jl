#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper!(::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    if SFI
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ added_layer.J₀⁻
    end

    # Batched multiplication between added and composite
    composite_layer.T⁻⁻[:] = added_layer.t⁻⁻ ⊠ composite_layer.T⁻⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        #J₀⁺, J₀⁻ = similar(composite_layer.J₀⁺), similar(composite_layer.J₀⁺)
        #J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻) 
        #J₀⁺ = added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ 
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻)
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺         
    end

    # Batched multiplication between added and composite
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻    
end

# Scattering in inhomogeneous composite layer.
# no scattering in homogeneous layer which is 
# added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ 
            (composite_layer.J₀⁺ .+ composite_layer.R⁺⁻ ⊠ added_layer.J₀⁻)
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.J₀⁻    
    end

    # Batched multiplication between added and composite
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer #these are aliases to the respective struct elements  
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer #these are aliases to the respective struct elements 
    
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) 
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    
    if SFI
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻ = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.J₀⁻) 
    end 

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    tmpR⁻⁺ = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺
    
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ = T01_inv ⊠ t⁻⁻ 

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) 
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    if SFI
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        tmpJ₀⁺ = added_layer.J₀⁺ .+ T21_inv ⊠ 
            (J₀⁺ .+ R⁺⁻ ⊠ added_layer.J₀⁻)
    end

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺ = T21_inv  ⊠ T⁺⁺ 
    
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻ = r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻  

    if SFI
        composite_layer.J₀⁺[:] = tmpJ₀⁺
        composite_layer.J₀⁻[:] = tmpJ₀⁻
    end
    composite_layer.R⁺⁻[:] = tmpR⁺⁻
    composite_layer.T⁻⁻[:] = tmpT⁻⁻
    composite_layer.R⁻⁺[:] = tmpR⁻⁺
    composite_layer.T⁺⁺[:] = tmpT⁺⁺
end

"Compute interaction between composite and added layers"
function interaction!(scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    
    
    #@show A1[1,1,1], A2[1,1,1]
    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    #A1 = Array(composite_layer.J₀⁻)
    #A2 = Array(composite_layer.J₀⁺)
    #@show A1[1,1,1], A2[1,1,1]
    synchronize_if_gpu()
end