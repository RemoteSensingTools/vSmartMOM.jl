#=
 
This file contains RT interaction-related functions (noRS)
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper_ms!(::ScatteringInterface_00, SFI,
                                #composite_layer::CompositeLayerMS{FT}, 
                                #added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT2},
                                r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                                t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
                                j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                                R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                                T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                                J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT}
                                ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    if SFI
        #tmpJ₀⁺, tmpJ₀⁻ = similar(J₀⁺), similar(J₀⁺)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
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
                                r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                                t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
                                j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                                R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                                T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                                J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT}
                                ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        #J₀⁺, J₀⁻ = similar(composite_layer.J₀⁺), similar(composite_layer.J₀⁺)
        #J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻) 
        #J₀⁺ = added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺ 
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
                                r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                                t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
                                j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                                R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                                T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                                J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT}
                                ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    if SFI
        #tmpJ₀⁺, tmpJ₀⁻ = similar(J₀⁺), similar(J₀⁺)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
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
                                r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                                t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT},  
                                j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                                R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                                T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                                J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT}
                                ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    #@unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    #@unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer
    
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)
    
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    
    if SFI
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻ = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 
    end 
    
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    tmpR⁻⁺ = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ = T01_inv ⊠ t⁻⁻ #Suniti

    

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    if SFI
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        tmpJ₀⁺ = j₀⁺ .+ T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
    end 

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    tmpT⁺⁺ = T21_inv  ⊠ T⁺⁺ #Suniti
    
    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    tmpR⁺⁻ = r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti
    
    if SFI
        J₀⁺[:] = tmpJ₀⁺
        J₀⁻[:] = tmpJ₀⁻
    end
    R⁺⁻[:] = tmpR⁺⁻
    T⁻⁻[:] = tmpT⁻⁻
    R⁻⁺[:] = tmpR⁻⁺
    T⁺⁺[:] = tmpT⁺⁺
end

"Compute interaction between composite and added layers above the sensor"
function interaction_top!(ims::Int64, 
                        RS_type::Union{noRS, noRS_plus}, 
                        scattering_interface::AbstractScatteringInterface, 
                        SFI,
                        composite_layer::CompositeLayerMS{M}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2},
                        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2,M}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @show size(composite_layer.topT⁺⁺)
    #@unpack topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻ = composite_layer
    R⁻⁺ = arr_type(composite_layer.topR⁻⁺[ims]) 
    R⁺⁻ = arr_type(composite_layer.topR⁺⁻[ims]) 
                        
    T⁺⁺ = arr_type(composite_layer.topT⁺⁺[ims]) 
    T⁻⁻ = arr_type(composite_layer.topT⁻⁻[ims]) 
                        
    compJ₀⁺ = arr_type(composite_layer.topJ₀⁺[ims]) 
    compJ₀⁻ = arr_type(composite_layer.topJ₀⁻[ims])
    
    interaction_helper_ms!(scattering_interface, SFI, #composite_layer, added_layer, 
                        I_static,
                        r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                        R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, compJ₀⁺, compJ₀⁻);

    composite_layer.topR⁻⁺[ims][:] = Array(R⁻⁺)
    composite_layer.topR⁺⁻[ims][:] = Array(R⁺⁻)
                        
    composite_layer.topT⁺⁺[ims][:] = Array(T⁺⁺)
    composite_layer.topT⁻⁻[ims][:] = Array(T⁻⁻)
                        
    composite_layer.topJ₀⁺[ims][:] = Array(compJ₀⁺) 
    composite_layer.topJ₀⁻[ims][:] = Array(compJ₀⁻)

    synchronize_if_gpu()
    #@pack composite_layer = topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻   
end

"Compute interaction between composite and added layers above the sensor"
function interaction_bot!(ims::Int64, 
                        RS_type::Union{noRS, noRS_plus}, 
                        scattering_interface::AbstractScatteringInterface, 
                        SFI,
                        composite_layer::CompositeLayerMS{M}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}, arr_type) where {M,FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @show scattering_interface
    #@unpack botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻ = composite_layer
    R⁻⁺ = arr_type(composite_layer.botR⁻⁺[ims]) 
    R⁺⁻ = arr_type(composite_layer.botR⁺⁻[ims]) 
                        
    T⁺⁺ = arr_type(composite_layer.botT⁺⁺[ims]) 
    T⁻⁻ = arr_type(composite_layer.botT⁻⁻[ims]) 
                        
    compJ₀⁺ = arr_type(composite_layer.botJ₀⁺[ims]) 
    compJ₀⁻ = arr_type(composite_layer.botJ₀⁻[ims])
    
    interaction_helper_ms!(scattering_interface, SFI, #composite_layer, added_layer, 
                        I_static,
                        r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                        R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, compJ₀⁺, compJ₀⁻);
        
    composite_layer.botR⁻⁺[ims][:] = Array(R⁻⁺)
    composite_layer.botR⁺⁻[ims][:] = Array(R⁺⁻)
                        
    composite_layer.botT⁺⁺[ims][:] = Array(T⁺⁺)
    composite_layer.botT⁻⁻[ims][:] = Array(T⁻⁻)
                        
    composite_layer.botJ₀⁺[ims][:] = Array(compJ₀⁺) 
    composite_layer.botJ₀⁻[ims][:] = Array(compJ₀⁻)
    
    synchronize_if_gpu()
    #@pack composite_layer = botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻
    
end


#=============================#
#= Multi-sensor computations =#
#= with inelastic scattering =#
#=============================#

#=
 
This file contains RT interaction-related functions (noRS)
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper_ms!(RS_type, ::ScatteringInterface_00, SFI,
    #composite_layer::CompositeLayerMS{FT}, 
    #added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2},
    r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
    t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
    j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
    ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
    iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
    iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
    R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
    T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
    J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
    ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
    ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
    ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
    ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    # If SFI, interact source function in no scattering
    if SFI

        ieJ₀⁺[:] = 0.0 #ieJ₀⁺
        ieJ₀⁻[:] = 0.0 #ieJ₀⁻
        #tmpJ₀⁺, tmpJ₀⁻ = similar(J₀⁺), similar(J₀⁺)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
    end

    # Batched multiplication between added and composite
    T⁻⁻[:] = t⁻⁻ ⊠ T⁻⁻
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper_ms!(RS_type::RRS, ::ScatteringInterface_01, SFI,
    #composite_layer::CompositeLayerMS{FT}, 
    #added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2},
    r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
    t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
    j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
    ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
    iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
    iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
    R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
    T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
    J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
    ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
    ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
    ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
    ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type 
    if SFI
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                ieJ₀⁻[:,1,n₁,Δn] = 
                    T⁻⁻[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] + 
                    iej₀⁻[:,:,n₁,Δn]) 
                ieJ₀⁺[:,1,n₁,Δn] = 
                    iej₀⁺[:,:,n₁,Δn] + 
                    iet⁺⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀]
            end 
        end
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺ 
    end
    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            ieR⁻⁺[:,:,n₁,Δn] = 
                T⁻⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn] ⊠ 
                T⁺⁺[:,:,n₀]

            ieR⁺⁻[:,:,n₁,Δn] = ier⁺⁻[:,:,n₁,Δn]

            ieT⁺⁺[:,:,n₁,Δn] = 
                iet⁺⁺[:,:,n₁,Δn] ⊠ T⁺⁺[:,:,n₀]

            ieT⁻⁻[:,:,n₁,Δn] = 
                T⁻⁻[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn]    
        end
    end
    # Batched multiplication between added and composite
    R⁻⁺[:] = T⁻⁻ ⊠ r⁻⁺ ⊠ T⁺⁺
    R⁺⁻[:] = r⁺⁻
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
    T⁻⁻[:] = T⁻⁻ ⊠ t⁻⁻    
end

function interaction_helper_ms!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                ::ScatteringInterface_01, SFI,
                #composite_layer::CompositeLayerMS{FT}, 
                #added_layer::AddedLayer{FT}, 
                I_static::AbstractArray{FT2},
                r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
                j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
                iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
                iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
                R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
                ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
                ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
                ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
                ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀_all = RS_type 
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                ieJ₀⁻[:,1,n₁,n₀] = 
                    T⁻⁻[:,:,n₁] * 
                    (ier⁻⁺[:,:,n₁,n₀] * J₀⁺[:,:,n₀] + 
                    iej₀⁻[:,:,n₁,n₀]) 
                ieJ₀⁺[:,1,n₁,n₀] = 
                    iej₀⁺[:,:,n₁,n₀] + 
                    iet⁺⁺[:,:,n₁,n₀] * J₀⁺[:,:,n₀]
            end 
        end
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ J₀⁺ 
    end
    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            # Batched multiplication between added and composite
            ieR⁻⁺[:,:,n₁,n₀] = 
                T⁻⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀] * 
                T⁺⁺[:,:,n₀]

            ieR⁺⁻[:,:,n₁,n₀] = ier⁺⁻[:,:,n₁,n₀]

            ieT⁺⁺[:,:,n₁,n₀] = 
                iet⁺⁺[:,:,n₁,n₀] * T⁺⁺[:,:,n₀]

            ieT⁻⁻[:,:,n₁,n₀] = 
                T⁻⁻[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀]    
        end
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
function interaction_helper_ms!(RS_type::RRS, ::ScatteringInterface_10, SFI,
    #composite_layer::CompositeLayerMS{FT}, 
    #added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2},
    r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
    t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
    j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
    ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
    iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
    iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
    R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
    T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
    J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
    ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
    ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
    ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
    ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type 
    if SFI
        for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                
                ieJ₀⁺[:,1,n₁,Δn] = t⁺⁺[:,:,n₁] ⊠ 
                        (ieJ₀⁺[:,1,n₁,Δn] + 
                        ieR⁺⁻[:,:,n₁,Δn] ⊠ j₀⁻[:,1,n₀])
                ieJ₀⁻[:,1,n₁,Δn] = ieJ₀⁻[:,1,n₁,Δn] + 
                        ieT⁻⁻[:,:,n₁,Δn] ⊠ j₀⁻[:,1,n₀]
            end
        end        
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
    end
    for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₀  = n₁ + i_λ₁λ₀[Δn]
            # Batched multiplication between added and composite
            ieT⁺⁺[:,:,n₁,Δn] = 
                    t⁺⁺[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn]
            ieT⁻⁻[:,:,n₁,Δn] = 
                    ieT⁻⁻[:,:,n₁,Δn] ⊠ t⁻⁻[:,:,n₀]
            ieR⁺⁻[:,:,n₁,Δn] = 
                    t⁺⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn] ⊠ 
                    t⁻⁻[:,:,n₀]
        end
    end
    # Batched multiplication between added and composite
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
    T⁻⁻[:] = T⁻⁻ ⊠ t⁻⁻
    R⁺⁻[:] = t⁺⁺ ⊠ R⁺⁻ ⊠ t⁻⁻
end

function interaction_helper_ms!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                            ::ScatteringInterface_10, SFI,
                            #composite_layer::CompositeLayerMS{FT}, 
                            #added_layer::AddedLayer{FT}, 
                            I_static::AbstractArray{FT2},
                            r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                            t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT}, 
                            j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                            ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
                            iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
                            iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
                            R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                            T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                            J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
                            ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
                            ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
                            ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
                            ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀_all = RS_type 
    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                
                ieJ₀⁺[:,1,n₁,n₀] = t⁺⁺[:,:,n₁] * 
                        (ieJ₀⁺[:,1,n₁,n₀] + 
                        ieR⁺⁻[:,:,n₁,n₀] * j₀⁻[:,1,n₀])
                ieJ₀⁻[:,1,n₁,n₀] = ieJ₀⁻[:,1,n₁,n₀] + 
                        ieT⁻⁻[:,:,n₁,n₀] * j₀⁻[:,1,n₀]
            end
        end        
        J₀⁺[:] = j₀⁺ .+ t⁺⁺ ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
        J₀⁻[:] = J₀⁻ .+ T⁻⁻ ⊠ j₀⁻
    end
    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            # Batched multiplication between added and composite
            ieT⁺⁺[:,:,n₁,n₀] = 
                    t⁺⁺[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀]
            ieT⁻⁻[:,:,n₁,n₀] = 
                    ieT⁻⁻[:,:,n₁,n₀] * t⁻⁻[:,:,n₀]
            ieR⁺⁻[:,:,n₁,n₀] = 
                    t⁺⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀] * 
                    t⁻⁻[:,:,n₀]
        end
    end
    # Batched multiplication between added and composite
    T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
    T⁻⁻[:] = T⁻⁻ ⊠ t⁻⁻
    R⁺⁻[:] = t⁺⁺ ⊠ R⁺⁻ ⊠ t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper_ms!(RS_type::RRS, ::ScatteringInterface_11, SFI,
    #composite_layer::CompositeLayerMS{FT}, 
    #added_layer::AddedLayer{FT}, 
    I_static::AbstractArray{FT2},
    r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
    t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT},  
    j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
    ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
    iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
    iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
    R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
    T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
    J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
    ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
    ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
    ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
    ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type 
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
            @inbounds @views ieJ₀⁻[:,:,n₁,Δn] = 
                    ieJ₀⁻[:,:,n₁,Δn] + 
                    T01_inv[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,:,n₁,Δn] +
                    iej₀⁻[:,:,n₁,Δn]) +
                    (T01_inv[:,:,n₁] ⊠ 
                    (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
                    ieT⁻⁻[:,:,n₁,Δn]) ⊠
                    tmp_inv[:,:,n₀] ⊠ 
                    (j₀⁻[:,:,n₀] + r⁻⁺[:,:,n₀] ⊠ J₀⁺[:,:,n₀]);
        end
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        J₀⁻[:] = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 
    end 
    for Δn = 1:size(ier⁻⁺,4)
        n₀, n₁ = get_n₀_n₁(ier⁻⁺,i_λ₁λ₀[Δn])
        @inbounds @views ieR⁻⁺[:,:,n₁,Δn] = ieR⁻⁺[:,:,n₁,Δn] +
            T01_inv[:,:,n₁] ⊠   
            (ier⁻⁺[:,:,n₁,Δn] ⊠ T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn]) +    
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) + 
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ r⁻⁺[:,:,n₀] ⊠ T⁺⁺[:,:,n₀]
        @inbounds @views ieT⁻⁻[:,:,n₁,Δn] = 
            T01_inv[:,:,n₁] ⊠ iet⁻⁻[:,:,n₁,Δn] +  
            (T01_inv[:,:,n₁] ⊠ 
            (ier⁻⁺[:,:,n₁,Δn] ⊠ R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieR⁺⁻[:,:,n₁,Δn]) +
            ieT⁻⁻[:,:,n₁,Δn]) ⊠ 
            tmp_inv[:,:,n₀] ⊠ t⁻⁻[:,:,n₀]
    end
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    R⁻⁺[:] = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    T⁻⁻[:] = T01_inv ⊠ t⁻⁻ #Suniti



    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    if SFI
        for Δn = 1:size(ieJ₀⁺,4)
            n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
            @inbounds @views ieJ₀⁺[:,:,n₁,Δn] = 
                            iej₀⁺[:,:,n₁,Δn] + 
                            T21_inv[:,:,n₁] ⊠ 
                            (ieJ₀⁺[:,:,n₁,Δn] + 
                            ieR⁺⁻[:,:,n₁,Δn] ⊠ j₀⁻[:,:,n₀] +
                            R⁺⁻[:,:,n₁] ⊠ iej₀⁻[:,:,n₁,Δn]) +
                            (T21_inv[:,:,n₁] ⊠ 
                            (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                            iet⁺⁺[:,:,n₁,Δn]) ⊠ 
                            tmp_inv[:,:,n₀] ⊠ (J₀⁺[:,:,n₀] + 
                            R⁺⁻[:,:,n₀] ⊠ j₀⁻[:,:,n₀])
        end
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        J₀⁺[:] = j₀⁺ .+ T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
    end 
    for Δn = 1:size(ieJ₀⁺,4)
        n₀, n₁ = get_n₀_n₁(ieJ₀⁺,i_λ₁λ₀[Δn])
        
        @inbounds @views ieT⁺⁺[:,:,n₁,Δn] = 
                    T21_inv[:,:,n₁] ⊠ ieT⁺⁺[:,:,n₁,Δn] +
                    (T21_inv[:,:,n₁] ⊠ (ieR⁺⁻[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + 
                    R⁺⁻[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) +
                    iet⁺⁺[:,:,n₁,Δn]) ⊠ tmp_inv[:,:,n₀] ⊠ T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

        @inbounds @views ieR⁺⁻[:,:,n₁,Δn] = 
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
    T⁺⁺[:] = T21_inv  ⊠ T⁺⁺ #Suniti

    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    R⁺⁻[:] = r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti
end

function interaction_helper_ms!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                        ::ScatteringInterface_11, SFI,
                        #composite_layer::CompositeLayerMS{FT}, 
                        #added_layer::AddedLayer{FT}, 
                        I_static::AbstractArray{FT2},
                        r⁺⁻::AbstractArray{FT}, r⁻⁺::AbstractArray{FT}, 
                        t⁻⁻::AbstractArray{FT}, t⁺⁺::AbstractArray{FT},  
                        j₀⁺::AbstractArray{FT}, j₀⁻::AbstractArray{FT},
                        ier⁺⁻::AbstractArray{FT}, ier⁻⁺::AbstractArray{FT}, 
                        iet⁻⁻::AbstractArray{FT}, iet⁺⁺::AbstractArray{FT}, 
                        iej₀⁺::AbstractArray{FT}, iej₀⁻::AbstractArray{FT},
                        R⁻⁺::AbstractArray{FT}, R⁺⁻::AbstractArray{FT}, 
                        T⁺⁺::AbstractArray{FT}, T⁻⁻::AbstractArray{FT}, 
                        J₀⁺::AbstractArray{FT}, J₀⁻::AbstractArray{FT},
                        ieR⁻⁺::AbstractArray{FT}, ieR⁺⁻::AbstractArray{FT}, 
                        ieT⁺⁺::AbstractArray{FT}, ieT⁻⁻::AbstractArray{FT}, 
                        ieJ₀⁺::AbstractArray{FT}, ieJ₀⁻::AbstractArray{FT}
                        ) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀_all = RS_type 
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) #Suniti
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;

    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                @inbounds @views ieJ₀⁻[:,1,n₁,n₀] = 
                        ieJ₀⁻[:,1,n₁,n₀] + 
                        T01_inv[:,:,n₁] * 
                        (ier⁻⁺[:,:,n₁,n₀] * J₀⁺[:,1,n₀] + 
                        r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,n₀] +
                        iej₀⁻[:,1,n₁,n₀]) +
                        (T01_inv[:,:,n₁] * 
                        (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + 
                        r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                        ieT⁻⁻[:,:,n₁,n₀]) *
                        tmp_inv[:,:,n₀] * 
                        (j₀⁻[:,1,n₀] + r⁻⁺[:,:,n₀] * J₀⁺[:,1,n₀]);
            end
        end
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        J₀⁻[:] = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻) 
    end 
    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            @inbounds @views ieR⁻⁺[:,:,n₁,n₀] = ieR⁻⁺[:,:,n₁,n₀] +
                T01_inv[:,:,n₁] *   
                (ier⁻⁺[:,:,n₁,n₀] * T⁺⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀]) +    
                (T01_inv[:,:,n₁] * 
                (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) + 
                ieT⁻⁻[:,:,n₁,n₀]) * 
                tmp_inv[:,:,n₀] * r⁻⁺[:,:,n₀] * T⁺⁺[:,:,n₀]
            @inbounds @views ieT⁻⁻[:,:,n₁,n₀] = 
                T01_inv[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀] +  
                (T01_inv[:,:,n₁] * 
                (ier⁻⁺[:,:,n₁,n₀] * R⁺⁻[:,:,n₀] + r⁻⁺[:,:,n₁] * ieR⁺⁻[:,:,n₁,n₀]) +
                ieT⁻⁻[:,:,n₁,n₀]) * 
                tmp_inv[:,:,n₀] * t⁻⁻[:,:,n₀]
        end
    end
    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀ 
    R⁻⁺[:] = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺ #Suniti
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    T⁻⁻[:] = T01_inv ⊠ t⁻⁻ #Suniti



    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) #Suniti
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv

    if SFI
        for Δn = 1:length(i_λ₁λ₀_all)
            n₁ = i_λ₁λ₀_all[Δn]
            n₀ = 1
            if (n₁>0)
                @inbounds @views ieJ₀⁺[:,1,n₁,n₀] = 
                                iej₀⁺[:,1,n₁,n₀] + 
                                T21_inv[:,:,n₁] * 
                                (ieJ₀⁺[:,1,n₁,n₀] + 
                                ieR⁺⁻[:,:,n₁,n₀] * j₀⁻[:,1,n₀] +
                                R⁺⁻[:,:,n₁] * iej₀⁻[:,1,n₁,n₀]) +
                                (T21_inv[:,:,n₁] * 
                                (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                                R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                                iet⁺⁺[:,:,n₁,n₀]) * 
                                tmp_inv[:,:,n₀] * (J₀⁺[:,1,n₀] + 
                                R⁺⁻[:,:,n₀] * j₀⁻[:,1,n₀])
            end
        end
        # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ )
        J₀⁺[:] = j₀⁺ .+ T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
    end 
    for Δn = 1:length(i_λ₁λ₀_all)
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
        
            @inbounds @views ieT⁺⁺[:,:,n₁,n₀] = 
                        T21_inv[:,:,n₁] * ieT⁺⁺[:,:,n₁,n₀] +
                        (T21_inv[:,:,n₁] * (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                        R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) +
                        iet⁺⁺[:,:,n₁,n₀]) * tmp_inv[:,:,n₀] * T⁺⁺[:,:,n₀] #Suniti: Eq 12 of Raman paper draft

            @inbounds @views ieR⁺⁻[:,:,n₁,n₀] = 
                        ier⁺⁻[:,:,n₁,n₀] + 
                        T21_inv[:,:,n₁] * 
                        (ieR⁺⁻[:,:,n₁,n₀] * t⁻⁻[:,:,n₀] +
                        R⁺⁻[:,:,n₁] * iet⁻⁻[:,:,n₁,n₀]) +
                        (T21_inv[:,:,n₁] * 
                        (ieR⁺⁻[:,:,n₁,n₀] * r⁻⁺[:,:,n₀] + 
                        R⁺⁻[:,:,n₁] * ier⁻⁺[:,:,n₁,n₀]) + 
                        iet⁺⁺[:,:,n₁,n₀]) *
                        tmp_inv[:,:,n₀] * R⁺⁻[:,:,n₀] * t⁻⁻[:,:,n₀]
        end
    end
    
    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    T⁺⁺[:] = T21_inv  ⊠ T⁺⁺ #Suniti

    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    R⁺⁻[:] = r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻ #Suniti
end

"Compute interaction between composite and added layers above the sensor"
function interaction_top!(ims::Int64, 
RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, 
scattering_interface::AbstractScatteringInterface, 
SFI,
composite_layer::CompositeLayerMSRS{M}, 
added_layer::AddedLayerRS{FT},
I_static::AbstractArray{FT2},
arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2,M}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    @show size(composite_layer.topT⁺⁺)
    #@unpack topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻ = composite_layer
    R⁻⁺ = arr_type(composite_layer.topR⁻⁺[ims]) 
    R⁺⁻ = arr_type(composite_layer.topR⁺⁻[ims]) 

    T⁺⁺ = arr_type(composite_layer.topT⁺⁺[ims]) 
    T⁻⁻ = arr_type(composite_layer.topT⁻⁻[ims]) 

    compJ₀⁺ = arr_type(composite_layer.topJ₀⁺[ims]) 
    compJ₀⁻ = arr_type(composite_layer.topJ₀⁻[ims])

    ieR⁻⁺ = arr_type(composite_layer.topieR⁻⁺[ims]) 
    ieR⁺⁻ = arr_type(composite_layer.topieR⁺⁻[ims]) 

    ieT⁺⁺ = arr_type(composite_layer.topieT⁺⁺[ims]) 
    ieT⁻⁻ = arr_type(composite_layer.topieT⁻⁻[ims]) 

    compieJ₀⁺ = arr_type(composite_layer.topieJ₀⁺[ims]) 
    compieJ₀⁻ = arr_type(composite_layer.topieJ₀⁻[ims])

    interaction_helper_ms!(RS_type, scattering_interface, SFI, #composite_layer, added_layer, 
                            I_static,
                            r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                            ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻,
                            R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, compJ₀⁺, compJ₀⁻,
                            ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, compieJ₀⁺, compieJ₀⁻);

    composite_layer.topR⁻⁺[ims][:] = Array(R⁻⁺)
    composite_layer.topR⁺⁻[ims][:] = Array(R⁺⁻)

    composite_layer.topT⁺⁺[ims][:] = Array(T⁺⁺)
    composite_layer.topT⁻⁻[ims][:] = Array(T⁻⁻)

    composite_layer.topJ₀⁺[ims][:] = Array(compJ₀⁺) 
    composite_layer.topJ₀⁻[ims][:] = Array(compJ₀⁻)

    composite_layer.topieR⁻⁺[ims][:] = Array(ieR⁻⁺)
    composite_layer.topieR⁺⁻[ims][:] = Array(ieR⁺⁻)

    composite_layer.topieT⁺⁺[ims][:] = Array(ieT⁺⁺)
    composite_layer.topieT⁻⁻[ims][:] = Array(ieT⁻⁻)

    composite_layer.topieJ₀⁺[ims][:] = Array(compieJ₀⁺) 
    composite_layer.topieJ₀⁻[ims][:] = Array(compieJ₀⁻)

    synchronize_if_gpu()
    #@pack composite_layer = topR⁻⁺, topR⁺⁻, topT⁺⁺, topT⁻⁻, topJ₀⁺, topJ₀⁻   
end


function interaction_bot!(ims::Int64, 
                    RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, 
                    scattering_interface::AbstractScatteringInterface, 
                    SFI,
                    composite_layer::CompositeLayerMSRS{M}, 
                    added_layer::AddedLayerRS{FT},
                    I_static::AbstractArray{FT2},
                    arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2,M}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    @show scattering_interface
    #@unpack botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻ = composite_layer
    R⁻⁺ = arr_type(composite_layer.botR⁻⁺[ims]) 
    R⁺⁻ = arr_type(composite_layer.botR⁺⁻[ims]) 

    T⁺⁺ = arr_type(composite_layer.botT⁺⁺[ims]) 
    T⁻⁻ = arr_type(composite_layer.botT⁻⁻[ims]) 

    compJ₀⁺ = arr_type(composite_layer.botJ₀⁺[ims]) 
    compJ₀⁻ = arr_type(composite_layer.botJ₀⁻[ims])

    ieR⁻⁺ = arr_type(composite_layer.botieR⁻⁺[ims]) 
    ieR⁺⁻ = arr_type(composite_layer.botieR⁺⁻[ims]) 

    ieT⁺⁺ = arr_type(composite_layer.botieT⁺⁺[ims]) 
    ieT⁻⁻ = arr_type(composite_layer.botieT⁻⁻[ims]) 

    compieJ₀⁺ = arr_type(composite_layer.botieJ₀⁺[ims]) 
    compieJ₀⁻ = arr_type(composite_layer.botieJ₀⁻[ims])

    interaction_helper_ms!(RS_type, scattering_interface, SFI, #composite_layer, added_layer, 
                            I_static,
                            r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻,
                            ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻,
                            R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, compJ₀⁺, compJ₀⁻,
                            ieR⁻⁺, ieR⁺⁻, ieT⁺⁺, ieT⁻⁻, compieJ₀⁺, compieJ₀⁻);

    composite_layer.botR⁻⁺[ims][:] = Array(R⁻⁺)
    composite_layer.botR⁺⁻[ims][:] = Array(R⁺⁻)

    composite_layer.botT⁺⁺[ims][:] = Array(T⁺⁺)
    composite_layer.botT⁻⁻[ims][:] = Array(T⁻⁻)

    composite_layer.botJ₀⁺[ims][:] = Array(compJ₀⁺) 
    composite_layer.botJ₀⁻[ims][:] = Array(compJ₀⁻)

    composite_layer.botieR⁻⁺[ims][:] = Array(ieR⁻⁺)
    composite_layer.botieR⁺⁻[ims][:] = Array(ieR⁺⁻)

    composite_layer.botieT⁺⁺[ims][:] = Array(ieT⁺⁺)
    composite_layer.botieT⁻⁻[ims][:] = Array(ieT⁻⁻)

    composite_layer.botieJ₀⁺[ims][:] = Array(compieJ₀⁺) 
    composite_layer.botieJ₀⁻[ims][:] = Array(compieJ₀⁻)

    synchronize_if_gpu()
    #@pack composite_layer = botR⁻⁺, botR⁺⁻, botT⁺⁺, botT⁻⁻, botJ₀⁺, botJ₀⁻
end