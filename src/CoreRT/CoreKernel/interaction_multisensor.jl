#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper_ms!(::ScatteringInterface_00, SFI,
                                #composite_layer::CompositeLayerMS{FT}, 
                                #added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2},
                                r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻,
                                R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    if SFI
        tmpJ₀⁺, tmpJ₀⁻ = similar(J₀⁺), similar(J₀⁺)
        tmpJ₀⁺ = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺
        tmpJ₀⁻ = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
        J₀⁺ .= tmpJ₀⁺
        J₀⁻ .= tmpJ₀⁻
    end

    # Batched multiplication between added and composite
    T⁻⁻[:] = t⁻⁻ ⊠ T⁻⁻
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper_ms!(::ScatteringInterface_01, SFI,
                                #composite_layer::CompositeLayerMS{FT}, 
                                #added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2},
                                r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻,
                                R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        #J₀⁺, J₀⁻ = similar(composite_layer.J₀⁺), similar(composite_layer.J₀⁺)
        #J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻) 
        #J₀⁺ = added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
        J₀⁺ .= j₀⁺ .+ t⁺⁺ ⊠ J₀⁺ 
        J₀⁻ .= J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
    end

    # Batched multiplication between added and composite
    R⁻⁺[:] = T⁻⁻ ⊠ r⁻⁺ ⊠ T⁺⁺
    R⁺⁻[:] = r⁺⁻
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
    T⁻⁻[:] = T⁻⁻ ⊠ t⁻⁻    
end

# Scattering in inhomogeneous composite layer.
# no scattering in homogeneous layer which is 
# added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper_ms!(::ScatteringInterface_10, SFI,
                                #composite_layer::CompositeLayerMS{FT}, 
                                #added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2},
                                r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻,
                                R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        tmpJ₀⁺, tmpJ₀⁻ = similar(J₀⁺), similar(J₀⁺)
        tmpJ₀⁺ = j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
        tmpJ₀⁻ = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
        J₀⁺ = tmpJ₀⁺
        J₀⁻ = tmpJ₀⁻
    end

    # Batched multiplication between added and composite
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
    T⁻⁻[:] = T⁻⁻ ⊠ t⁻⁻
    R⁺⁻[:] = t⁺⁺ ⊠ R⁺⁻ ⊠ t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper_ms!(::ScatteringInterface_11, SFI,
                                #composite_layer::CompositeLayerMS{FT}, 
                                #added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2},
                                r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺,  j₀⁺, j₀⁻,
                                R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    #@unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    #@unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)
    tmpR⁻⁺ = similar(R⁻⁺)
    tmpT⁻⁻ = similar(T⁻⁻)
    tmpJ₀⁻ = similar(J₀⁻)
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    @show size(T⁻⁻), size(tmp_inv)
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    tmpR⁻⁺[:] = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻[:] = T01_inv ⊠ t⁻⁻ #Suniti

    if SFI
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻[:] = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 
    end 

    # Repeating for mirror-reflected directions
    tmpR⁺⁻ = similar(R⁺⁻)
    tmpT⁺⁺ = similar(T⁺⁺)
    tmpJ₀⁺ = similar(J₀⁺)
    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺[:] = T21_inv  ⊠ T⁺⁺ #Suniti
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻[:] = r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti
    
    if SFI
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        tmpJ₀⁺[:] = j₀⁺ .+ T21_inv ⊠ (J₀⁺ + R⁺⁻ ⊠ j₀⁻)
    end 

    R⁻⁺ = tmpR⁻⁺
    T⁻⁻ = tmpT⁻⁻
    J₀⁻ = tmpJ₀⁻
    R⁺⁻ = tmpR⁺⁻
    T⁺⁺ = tmpT⁺⁺
    J₀⁺ = tmpJ₀⁺
end

"Compute interaction between composite and added layers above the sensor"
function interaction_top!(ims::Int64, 
                        RS_type::Union{noRS, noRS_plus}, 
                        scattering_interface::AbstractScatteringInterface, 
                        SFI,
                        composite_layer::CompositeLayerMS{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @show size(composite_layer.topT⁺⁺)
    #@unpack topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻ = composite_layer
    interaction_helper_ms!(scattering_interface, SFI, #composite_layer, added_layer, 
                        I_static,
                        r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                        composite_layer.topR⁻⁺[ims,:,:,:], composite_layer.topR⁺⁻[ims,:,:,:], 
                        composite_layer.topT⁺⁺[ims,:,:,:], composite_layer.topT⁻⁻[ims,:,:,:], 
                        composite_layer.topJ₀⁺[ims,:,:,:], composite_layer.topJ₀⁻[ims,:,:,:]);
    synchronize_if_gpu()
    #@pack composite_layer = topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻   
end

"Compute interaction between composite and added layers above the sensor"
function interaction_bot!(ims::Int64, 
                        RS_type::Union{noRS, noRS_plus}, 
                        scattering_interface::AbstractScatteringInterface, 
                        SFI,
                        composite_layer::CompositeLayerMS{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    #@unpack botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻ = composite_layer
    interaction_helper_ms!(scattering_interface, SFI, #composite_layer, added_layer, 
                        I_static,
                        r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                        composite_layer.botR⁻⁺[ims,:], composite_layer.botR⁺⁻[ims,:], 
                        composite_layer.botT⁺⁺[ims,:], composite_layer.botT⁻⁻[ims,:], 
                        composite_layer.botJ₀⁺[ims,:], composite_layer.botJ₀⁻[ims,:]);
    synchronize_if_gpu()
    #@pack composite_layer = botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻
    
end