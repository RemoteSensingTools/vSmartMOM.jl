#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper!(::ScatteringInterface_00, SFI,
    composite_layer::CompositeLayer{FT},
    lin_composite_layer::linCompositeLayer{FT}, 
    added_layer::AddedLayer{FT},  
    lin_added_layer::linAddedLayer{FT}, 
    nparams::Int,
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 
    @unpack dxr⁺⁻, dxr⁻⁺, dxt⁻⁻, dxt⁺⁺, dxj₀⁺, dxj₀⁻ = lin_added_layer     
    @unpack dR⁻⁺, dR⁺⁻, dT⁺⁺, dT⁻⁻, dJ₀⁺, dJ₀⁻ = lin_composite_layer 

    for ctr=1:nparams
        # Source Function
        dJ₀⁺[ctr,:,1,:] .= dxj₀⁺[ctr,:,1,:] .+ 
                    t⁺⁺[:,:,:] ⊠ dJ₀⁺[ctr,:,1,:] .+ 
                    dxt⁺⁺[ctr,:,:,:] ⊠ J₀⁺[:,1,:]
        dJ₀⁻[ctr,:,1,:] .= dJ₀⁻[ctr,:,1,:] .+ 
                    T⁻⁻[:,:,:] ⊠ dxj₀⁻[ctr,:,1,:] .+ 
                    dT⁻⁻[ctr,:,:,:] ⊠ j₀⁻[:,1,:]

        # Batched multiplication between added and composite
        dT⁻⁻[ctr,:,:,:]  .= 
                    dxt⁻⁻[ctr,:,:,:] ⊠ T⁻⁻[:,:,:] .+ 
                    t⁻⁻[:,:,:] ⊠ dT⁻⁻[ctr,:,:,:]
        dT⁺⁺[ctr,:,:,:]  .= 
                    dxt⁺⁺[ctr,:,:,:] ⊠ T⁺⁺[:,:,:] .+ 
                    t⁺⁺[:,:,:] ⊠ dT⁺⁺[ctr,:,:,:]
    end
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
    lin_composite_layer::linCompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    lin_added_layer::linAddedLayer{FT}, 
    nparams::Int
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 
    @unpack dxr⁺⁻, dxr⁻⁺, dxt⁻⁻, dxt⁺⁺, dxj₀⁺, dxj₀⁻ = lin_added_layer         
    @unpack dR⁻⁺, dR⁺⁻, dT⁺⁺, dT⁻⁻, dJ₀⁺, dJ₀⁻ = lin_composite_layer 

    for ctr=1:nparams
        # Source Function
        dJ₀⁻[ctr,:,1,:] .= dJ₀⁻[ctr,:,1,:] .+ 
            dT⁻⁻[ctr,:,:,:] ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) .+
            T⁻⁻[ctr,:,:,:] ⊠ 
                (dxr⁻⁺[ctr,:,:,:] ⊠ J₀⁺ .+
                r⁻⁺ ⊠ dJ₀⁺[ctr,:,1,:] .+ 
                dxj₀⁻[ctr,:,1,:])

        dJ₀⁺[ctr,:,1,:] .= dxj₀⁺[ctr,:,1,:] .+ 
            dxt⁺⁺[ctr,:,:,:] ⊠ J₀⁺ .+ 
            t⁺⁺ ⊠ dJ₀⁺[ctr,:,1,:]         

        # Batched multiplication between added and composite
        dR⁻⁺[ctr,:,:,:] .= dT⁻⁻[ctr,:,:,:] ⊠ r⁻⁺ ⊠ T⁺⁺ .+
                    T⁻⁻ ⊠ dxr⁻⁺[ctr,:,:,:] ⊠ T⁺⁺ .+
                    T⁻⁻ ⊠ r⁻⁺ ⊠ dT⁺⁺[ctr,:,:,:]
        dR⁺⁻[ctr,:,:,:] .= dxr⁺⁻[ctr,:,:,:]
        dT⁺⁺[ctr,:,:,:] .= dxt⁺⁺[ctr,:,:,:] ⊠ T⁺⁺ .+
                    t⁺⁺ ⊠ dT⁺⁺[ctr,:,:,:]
        dT⁻⁻[ctr,:,:,:] .= dT⁻⁻[ctr,:,:,:] ⊠ t⁻⁻ .+
                    T⁻⁻ ⊠ dxt⁻⁻[ctr,:,:,:]    
    end
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
    lin_composite_layer::linCompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    lin_added_layer::linAddedLayer{FT},
    nparams::Int,
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer 
    @unpack dxr⁺⁻, dxr⁻⁺, dxt⁻⁻, dxt⁺⁺, dxj₀⁺, dxj₀⁻ = lin_added_layer     
    @unpack dR⁻⁺, dR⁺⁻, dT⁺⁺, dT⁻⁻, dJ₀⁺, dJ₀⁻ = lin_composite_layer 

    for ctr=1:nparams
        # Source Function
        dJ₀⁺[ctr,:,1,:]  .= dxj₀⁺[ctr,:,1,:] .+ 
            dxt⁺⁺[ctr,:,:,:]  ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻) .+
            t⁺⁺ ⊠ (dJ₀⁺[ctr,:,1,:]  .+ dR⁺⁻[ctr,:,:,:]  ⊠ j₀⁻ .+ R⁺⁻ ⊠ dxj₀⁻[ctr,:,1,:])
        dJ₀⁻[ctr,:,1,:]  .= dJ₀⁻[ctr,:,1,:]  .+ 
            dT⁻⁻[ctr,:,:,:] ⊠ j₀⁻ .+ T⁻⁻ ⊠ dxj₀⁻[ctr,:,1,:]    

        # Batched multiplication between added and composite
        dT⁺⁺[ctr,:,:,:] .= dxt⁺⁺[ctr,:,:,:] ⊠ T⁺⁺ .+ t⁺⁺ ⊠ dT⁺⁺[ctr,:,:,:]
        dT⁻⁻[ctr,:,:,:] .= dT⁻⁻[ctr,:,:,:] ⊠ t⁻⁻ .+ T⁻⁻ ⊠ dxt⁻⁻[ctr,:,:,:]
        dR⁺⁻[ctr,:,:,:] .= dxt⁺⁺[ctr,:,:,:] ⊠ R⁺⁻ ⊠ t⁻⁻ .+ 
            t⁺⁺ ⊠ (dR⁺⁻[ctr,:,:,:] ⊠ t⁻⁻ .+ R⁺⁻ ⊠ dxt⁻⁻[ctr,:,:,:])
    end
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
    lin_composite_layer::linCompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    lin_added_layer::linAddedLayer{FT}, 
    nparams::Int,
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer     #these are aliases to the respective struct elements  
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer #these are aliases to the respective struct elements 
    @unpack dxr⁺⁻, dxr⁻⁺, dxt⁻⁻, dxt⁺⁺, dxj₀⁺, dxj₀⁻ = lin_added_layer     
    @unpack dR⁻⁺, dR⁺⁻, dT⁺⁺, dT⁻⁻, dJ₀⁺, dJ₀⁻ = lin_composite_layer 

    # X₂₁ refers to added layer, X₁₀ to composite layer!

    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) 
    # Temporary arrays:

    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;

    dtmp_inv = zeros(nparams,size(t⁺⁺,1), ,size(t⁺⁺,2), ,size(t⁺⁺,3))
    dT01_inv = similar(dtmp_inv)
    for ctr=1:nparams
        dtmp_inv[ctr,:,:,:] = tmp_inv ⊠ 
            (dxr⁻⁺[ctr,:,:,:] ⊠ R⁺⁻ .+ 
            r⁻⁺ ⊠ dR⁺⁻[ctr,:,:,:]) ⊠ tmp_inv
        dT01_inv[ctr,:,:,:] = dT⁻⁻[ctr,:,:,:] ⊠ tmp_inv .+
            T⁻⁻ ⊠ dtmp_inv[ctr,:,:,:]

        # J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        dJ₀⁻[ctr,:,1,:] .= dJ₀⁻[ctr,:,1,:] .+ 
            dT01_inv[ctr,:,:,:] ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) .+
            T01_inv ⊠ (dxr⁻⁺[ctr,:,:,:] ⊠ J₀⁺ .+ 
                        r⁻⁺ ⊠ dJ₀⁺[ctr,:,1,:] .+ 
                        dxj₀⁻[ctr,:,1,:]) 

        # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
        dR⁻⁺[ctr,:,:,:] .= dR⁻⁺[ctr,:,:,:] .+ 
            dT01_inv[ctr,:,:,:] ⊠ r⁻⁺ ⊠ T⁺⁺ .+ 
            T01_inv ⊠ (dxr⁻⁺[ctr,:,:,:] ⊠ T⁺⁺ .+ 
                        r⁻⁺ ⊠ dT⁺⁺[ctr,:,:,:])

        # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
        dT⁻⁻[ctr,:,:,:] .= dT01_inv[ctr,:,:,:] ⊠ t⁻⁻ .+ 
            T01_inv ⊠ dxt⁻⁻[ctr,:,:,:] 
    end

    # J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
    J₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    R⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺

    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    T⁻⁻ .= T01_inv ⊠ t⁻⁻ 

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) 
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    dtmp_inv .= 0.0
    dT21_inv .= similar(dtmp_inv)

    for ctr=1:nparams
        dtmp_inv[ctr,:,:,:] = tmp_inv ⊠ 
            (dR⁺⁻[ctr,:,:,:] ⊠ r⁻⁺ .+ 
            R⁺⁻ ⊠ dxr⁻⁺[ctr,:,:,:]) ⊠ tmp_inv
        dT01_inv[ctr,:,:,:] = dxt⁺⁺[ctr,:,:,:] ⊠ tmp_inv .+
            t⁺⁺ ⊠ dtmp_inv[ctr,:,:,:]

        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        dJ₀⁺[ctr,:,1,:] .= dxj₀⁺[ctr,:,1,:] .+ 
            dT21_inv[ctr,:,:,:] ⊠ (R⁺⁻ ⊠ j₀⁻ .+ J₀⁺) .+
            T21_inv ⊠ (dR⁺⁻[ctr,:,:,:] ⊠ j₀⁻ .+ 
                        R⁺⁻ ⊠ dxj₀⁻[ctr,:,1,:] .+ 
                        dJ₀⁺[ctr,:,1,:]) 

        # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
        dR⁻⁺[ctr,:,:,:] .= dxr⁺⁻[ctr,:,:,:] .+ 
            dT21_inv[ctr,:,:,:] ⊠ R⁺⁻ ⊠ t⁻⁻ .+ 
            T21_inv ⊠ (dR⁺⁻[ctr,:,:,:] ⊠ t⁻⁻ .+ 
                        R⁺⁻ ⊠ dxt⁻⁻[ctr,:,:,:])

        # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
        dT⁺⁺[ctr,:,:,:] .= dT21_inv[ctr,:,:,:] ⊠ T⁺⁺ .+ 
            T21_inv ⊠ dT⁺⁺[ctr,:,:,:] 
    end


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
lin_composite_layer::linCompositeLayer{FT}, 
added_layer::AddedLayer{FT}, 
lin_added_layer::linAddedLayer{FT}, 
nparams::Int,
I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    interaction_helper!(scattering_interface, SFI, 
            composite_layer, 
            lin_composite_layer, 
            added_layer, 
            lin_added_layer,
            nparams,
            I_static)
    synchronize_if_gpu()
end
