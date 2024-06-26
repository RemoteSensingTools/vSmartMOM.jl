#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper!(RS_type,::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    if SFI
        composite_layer.ieJ₀⁺[:] = 0.0 #ieJ₀⁺
        composite_layer.ieJ₀⁻[:] = 0.0 #ieJ₀⁻

        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    # Batched multiplication between added and composite
    composite_layer.T⁻⁻[:] = added_layer.t⁻⁻ ⊠ composite_layer.T⁻⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack i_λ₁λ₀ = RS_type                             
    if SFI   
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                composite_layer.ieJ₀⁻[:,1,n₁,Δn] = 
                    composite_layer.T⁻⁻[:,:,n₁] * 
                    (added_layer.ier⁻⁺[:,:,n₁,Δn] * composite_layer.J₀⁺[:,1,n₀] + 
                    added_layer.ieJ₀⁻[:,1,n₁,Δn]) 
                composite_layer.ieJ₀⁺[:,1,n₁,Δn] = 
                    added_layer.ieJ₀⁺[:,1,n₁,Δn] + 
                    added_layer.iet⁺⁺[:,:,n₁,Δn] * composite_layer.J₀⁺[:,1,n₀]
            end 
        end

        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻) 
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
    end

    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            composite_layer.ieR⁻⁺[:,:,n₁,Δn] = 
                composite_layer.T⁻⁻[:,:,n₁] * added_layer.ier⁻⁺[:,:,n₁,Δn] * 
                composite_layer.T⁺⁺[:,:,n₀]
    
            composite_layer.ieR⁺⁻[:,:,n₁,Δn] = added_layer.ier⁺⁻[:,:,n₁,Δn]
    
            composite_layer.ieT⁺⁺[:,:,n₁,Δn] = 
                added_layer.iet⁺⁺[:,:,n₁,Δn] * composite_layer.T⁺⁺[:,:,n₀]
    
            composite_layer.ieT⁻⁻[:,:,n₁,Δn] = 
                composite_layer.T⁻⁻[:,:,n₁] * added_layer.iet⁻⁻[:,:,n₁,Δn]    
        end
    end
    # Batched multiplication between added and composite
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻    
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                        ::ScatteringInterface_01, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT}, 
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack i_λ₁λ₀_all = RS_type                             
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                # TODO: replace all the n₀ indices with ref_r, ref_t, and ref_J
                composite_layer.ieJ₀⁻[:,1,n₁,n₀] = 
                    composite_layer.T⁻⁻[:,:,n₁] * (added_layer.ier⁻⁺[:,:,n₁,n₀] * 
                    composite_layer.J₀⁺[:,1,n₀] + 
                    added_layer.ieJ₀⁻[:,1,n₁,n₀]) 

                composite_layer.ieJ₀⁺[:,1,n₁,n₀] = 
                    added_layer.ieJ₀⁺[:,1,n₁,n₀] + 
                    added_layer.iet⁺⁺[:,:,n₁,n₀] * composite_layer.J₀⁺[:,1,n₀] 
            end
        end

        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ 
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.j₀⁻) 
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
    end

    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        # Batched multiplication between added and composite
        if (n₁>0)
            composite_layer.ieR⁻⁺[:,:,n₁,n₀] = 
                        composite_layer.T⁻⁻[:,:,n₁] * added_layer.ier⁻⁺[:,:,n₁,n₀] * 
                        composite_layer.T⁺⁺[:,:,n₀]

            composite_layer.ieR⁺⁻[:,:,n₁,n₀] = added_layer.ier⁺⁻[:,:,n₁,n₀]

            composite_layer.ieT⁺⁺[:,:,n₁,n₀] = 
                        added_layer.iet⁺⁺[:,:,n₁,n₀] * composite_layer.T⁺⁺[:,:,n₀]

            composite_layer.ieT⁻⁻[:,:,n₁,n₀] = 
                        composite_layer.T⁻⁻[:,:,n₁] * added_layer.iet⁻⁻[:,:,n₁,n₀] 
        end   
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
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type 
    if SFI
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                
                composite_layer.ieJ₀⁺[:,1,n₁,Δn] = added_layer.t⁺⁺[:,:,n₁] * 
                        (composite_layer.ieJ₀⁺[:,1,n₁,Δn] + 
                        composite_layer.ieR⁺⁻[:,:,n₁,Δn] * added_layer.j₀⁻[:,1,n₀])
                composite_layer.ieJ₀⁻[:,1,n₁,Δn] = composite_layer.ieJ₀⁻[:,1,n₁,Δn] + 
                        composite_layer.ieT⁻⁻[:,:,n₁,Δn] * added_layer.j₀⁻[:,1,n₀]
            end
        end
        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ (composite_layer.J₀⁺ .+ 
                composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻)
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            composite_layer.ieT⁺⁺[:,:,n₁,Δn] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieT⁺⁺[:,:,n₁,Δn]
            composite_layer.ieT⁻⁻[:,:,n₁,Δn] = 
                    composite_layer.ieT⁻⁻[:,:,n₁,Δn] * added_layer.t⁻⁻[:,:,n₀]
            composite_layer.ieR⁺⁻[:,:,n₁,Δn] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieR⁺⁻[:,:,n₁,Δn] * 
                    added_layer.t⁻⁻[:,:,n₀]
        end
    end
    # Batched multiplication between added and composite
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
    ::ScatteringInterface_10, SFI,
    composite_layer::CompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀_all = RS_type 
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                composite_layer.ieJ₀⁺[:,1,n₁,n₀] = added_layer.t⁺⁺[:,:,n₁] * 
                            (composite_layer.ieJ₀⁺[:,1,n₁,n₀] + 
                            composite_layer.ieR⁺⁻[:,:,n₁,n₀] * added_layer.j₀⁻[:,1,n₀])

                composite_layer.ieJ₀⁻[:,1,n₁,n₀] = composite_layer.ieJ₀⁻[:,1,n₁,n₀] + 
                            composite_layer.ieT⁻⁻[:,:,n₁,n₀] * added_layer.j₀⁻[:,1,n₀]
            end
        end

        composite_layer.J₀⁺ = added_layer.j₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ (composite_layer.J₀⁺ .+ 
            composite_layer.R⁺⁻ ⊠ added_layer.j₀⁻)
        composite_layer.J₀⁻ = composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.j₀⁻
    end

    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        # Batched multiplication between added and composite
        if (n₁>0)
            composite_layer.ieT⁺⁺[:,:,n₁,n₀] = 
                    added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieT⁺⁺[:,:,n₁,n₀]
            composite_layer.ieT⁻⁻[:,:,n₁,n₀] = 
                composite_layer.ieT⁻⁻[:,:,n₁,n₀] * added_layer.t⁻⁻[:,:,n₀]
            composite_layer.ieR⁺⁻[:,:,n₁,n₀] = 
                added_layer.t⁺⁺[:,:,n₁] * composite_layer.ieR⁺⁻[:,:,n₁,n₀] * 
                added_layer.t⁻⁻[:,:,n₀]
        end
    end

    # Batched multiplication between added and composite
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(RS_type::RRS, ::ScatteringInterface_11, SFI,
                                composite_layer::Union{CompositeLayer, CompositeLayerRS}, 
                                added_layer::Union{AddedLayer,AddedLayerRS}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ = added_layer
    @unpack ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, ieJ₀⁺, ieJ₀⁻ = composite_layer
    
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    if SFI
        for Δn = 1:size(ieJ₀⁺,4)
            n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
            @inbounds @views composite_layer.ieJ₀⁻[:,:,n₁,Δn] = 
                    ieJ₀⁻[:,:,n₁,Δn] + 
                    T01_inv[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,:,n₁,Δn] +
                    added_layer.ieJ₀⁻[:,:,n₁,Δn]) +
                    (T01_inv[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
                    ieT⁻⁻[:,:,n₁,Δn]) ⊠
                    tmp_inv[:,:,n₀] ⊠ 
                    (added_layer.j₀⁻[:,:,n₀] + r⁻⁺[:,:,n₀] ⊠ J₀⁺[:,:,n₀]);
        end
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        composite_layer.J₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.j₀⁻) 
    end 
    for Δn = 1:size(ier⁻⁺,4)
        n₀, n₁ = get_n₀_n₁(ier⁻⁺,i_λ₁λ₀[Δn])
        @inbounds @views composite_layer.ieR⁻⁺[:,:,n₁,Δn] = ieR⁻⁺[:,:,n₁,Δn] +
            T01_inv[:,:,n₁] ⊠   
            (ier⁻⁺[:,:,n₁,Δn] ⊠ T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn]) +    
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) + 
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ r⁻⁺[:,:,n₀] ⊠ T⁺⁺[:,:,n₀]
        @inbounds @views composite_layer.ieT⁻⁻[:,:,n₁,Δn] = 
            T01_inv[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn] +  
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ t⁻⁻[:,:,n₀]
    end
    
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    composite_layer.R⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    composite_layer.T⁻⁻ .= T01_inv ⊠ t⁻⁻ #Suniti

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    if SFI
        for Δn = 1:size(ieJ₀⁺,4)
            n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
            @inbounds @views composite_layer.ieJ₀⁺[:,:,n₁,Δn] = 
                            added_layer.ieJ₀⁺[:,:,n₁,Δn] + 
                            T21_inv[:,:,n₁] ⊠ 
                            (ieJ₀⁺[:,:,n₁,Δn] + 
                            ieR⁺⁻[:,:,n₁,Δn] ⊠ added_layer.j₀⁻[:,:,n₀] +
                            R⁺⁻[:,:,n₁] ⊠ added_layer.ieJ₀⁻[:,:,n₁,Δn]) +
                            (T21_inv[:,:,n₁] ⊠ 
                            (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                            iet⁺⁺[:,:,n₁,Δn]) ⊠ 
                            tmp_inv[:,:,n₀] ⊠ (J₀⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₀] ⊠ added_layer.j₀⁻[:,:,n₀])
        end
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        composite_layer.J₀⁺ .= added_layer.j₀⁺ .+ T21_inv ⊠ (J₀⁺ + R⁺⁻ ⊠ added_layer.j₀⁻)
    end 
    for Δn = 1:size(ieJ₀⁺,4)
        n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
        
        @inbounds @views composite_layer.ieT⁺⁺[:,:,n₁,Δn] = 
                    T21_inv[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn] +
                    (T21_inv[:,:,n₁] ⊠ (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                    iet⁺⁺[:,:,n₁,Δn]) ⊠ tmp_inv[:,:,n₀] ⊠ T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

        @inbounds @views composite_layer.ieR⁺⁻[:,:,n₁,Δn] = 
                    ier⁺⁻[:,:,n₁,Δn] + 
                    T21_inv[:,:,n₁] ⊠ 
                    (ieR⁺⁻[:,:,n₁,Δn] ⊠ t⁻⁻[:,:,n₀] +
                    R⁺⁻[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn]) +
                    (T21_inv[:,:,n₁] ⊠ 
                    (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) + 
                    iet⁺⁺[:,:,n₁,Δn]) ⊠
                    tmp_inv[:,:,n₀] ⊠ R⁺⁻[:,:,n₀] ⊠ t⁻⁻[:,:,n₀]
    end
    
    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    composite_layer.T⁺⁺ .= T21_inv  ⊠ T⁺⁺ #Suniti
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    composite_layer.R⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti
    
end

function interaction_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
    ::ScatteringInterface_11, SFI,
    composite_layer::Union{CompositeLayer, CompositeLayerRS}, 
    added_layer::Union{AddedLayer,AddedLayerRS}, 
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
 
    @unpack i_λ₁λ₀_all = RS_type
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ = added_layer
    @unpack ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, ieJ₀⁺, ieJ₀⁻ = composite_layer
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;

    if SFI
        for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                composite_layer.ieJ₀⁻[:,1,n₁,1] = ieJ₀⁻[:,:,n₁,n₀] + 
                                        T01_inv[:,:,n₁] *
                                        (ier⁻⁺[:,:,n₁,n₀] * J₀⁺[:,:,n₀] + 
                                        r⁻⁺[:,:,n₁] * ieJ₀⁺[:,:,n₁,n₀] +
                                        added_layer.ieJ₀⁻[:,1,n₁,n₀]) +
                                        (T01_inv[:,:,n₁] * 
                                        (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + 
                                        r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                                        ieT⁻⁻[:,:,n₁,n₀]) *
                                        tmp_inv[:,:,n₀] * 
                                        (added_layer.j₀⁻[:,:,n₀] + r⁻⁺[:,:,n₀] * J₀⁺[:,:,n₀]) #Suniti: Eq 17 of Raman paper draft
            end
        end
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        composite_layer.J₀⁻[:] = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.j₀⁻) 
    end 

    for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            composite_layer.ieR⁻⁺[:,:,n₁,1] = ieR⁻⁺[:,:,n₁,n₀] +
                T01_inv[:,:,n₁] * 
                (ier⁻⁺[:,:,n₁,n₀] * T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀]) +    
                (T01_inv[:,:,n₁] * 
                (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) + 
                ieT⁻⁻[:,:,n₁,n₀]) * 
                tmp_inv[:,:,n₀] * r⁻⁺[:,:,n₀] * T⁺⁺[:,:,n₀] #Suniti: Eq 14 of Raman paper draft

            composite_layer.ieT⁻⁻[:,:,n₁,1] = T01_inv[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀] +  
                (T01_inv[:,:,n₁] * 
                (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                ieT⁻⁻[:,:,n₁,n₀]) * 
                tmp_inv[:,:,n₀] * t⁻⁻[:,:,n₀] #Suniti: Eq 13 of Raman paper draft
        end
    end

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    composite_layer.R⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    composite_layer.T⁻⁻ .= T01_inv ⊠ t⁻⁻ #Suniti

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all) #Δn in eachindex ieJ₁⁺[1,1,:,1]
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                composite_layer.ieJ₀⁺[:,1,n₁,n₀] = added_layer.ieJ₀⁺[:,:,n₁,n₀] + 
                        T21_inv[:,:,n₁] * 
                        (ieJ₀⁺[:,:,n₁,n₀] + 
                        ieR⁺⁻[:,:,n₁,n₀] * added_layer.j₀⁻[:,:,n₀] +
                        R⁺⁻[:,:,n₁] * added_layer.ieJ₀⁻[:,1,n₁,n₀]) +
                        (T21_inv[:,:,n₁] * 
                        (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                        R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                        iet⁺⁺[:,:,n₁,n₀]) * 
                        tmp_inv[:,:,n₀] * (J₀⁺[:,:,n₀] + 
                        R⁺⁻[:,:,n₀] * added_layer.j₀⁻[:,:,n₀])
            end
        end
    # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
    composite_layer.J₀⁺ .= added_layer.j₀⁺ .+ T21_inv ⊠ (J₀⁺ + R⁺⁻ ⊠ added_layer.j₀⁻)
    end 
    for Δn = 1:length(i_λ₁λ₀_all) #Δn in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            composite_layer.ieT⁺⁺[:,:,n₁,n₀] = T21_inv[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀] +
                    (T21_inv[:,:,n₁] * (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                    iet⁺⁺[:,:,n₁,n₀]) * tmp_inv[:,:,n₀] * T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

            composite_layer.ieR⁺⁻[:,:,n₁,n₀] = ier⁺⁻[:,:,n₁,n₀] + 
                    T21_inv[:,:,n₁] *
                    (ieR⁺⁻[:,:,n₁,n₀] * t⁻⁻[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀]) +
                    (T21_inv[:,:,n₁] * 
                    (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) + 
                    iet⁺⁺[:,:,n₁,n₀]) *
                    tmp_inv[:,:,n₀] * R⁺⁻[:,:,n₀] * t⁻⁻[:,:,n₀] #Suniti: Eq 15 of Raman paper draft
        end
    end
    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    composite_layer.T⁺⁺ .= T21_inv  ⊠ T⁺⁺ #Suniti
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    composite_layer.R⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti

end

"Compute interaction between composite and added layers"
function interaction!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::Union{CompositeLayer,CompositeLayerRS}, 
                        added_layer::Union{AddedLayer,AddedLayerRS},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    interaction_helper!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
    
end

function interaction!(RS_type::noRS, scattering_interface::AbstractScatteringInterface, SFI,
    composite_layer::Union{CompositeLayer,CompositeLayerRS}, 
    added_layer::Union{AddedLayer,AddedLayerRS},
    I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
    
    synchronize_if_gpu()

end