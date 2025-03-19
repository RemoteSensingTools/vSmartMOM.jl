#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function interaction_helper!(::ScatteringInterface_00, SFI,
                                computed_layer_properties, 
                                computed_layer_properties_lin, 
                                composite_layer::CompositeLayer{FT}, 
                                composite_layer_lin::CompositeLayerLin{FT}, 
                                added_layer::AddedLayer{FT}, 
                                added_layer_lin::AddedLayerLin{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[1]
    
    # If SFI, interact source function in no scattering
    if SFI
        for iparam=1:Nparams 
            composite_layer_lin.J̇₀⁺[iparam,:] .= added_layer_lin.ap_J̇₀⁺[iparam,:] .+ 
                added_layer.t⁺⁺ ⊠ composite_layer_lin.J̇₀⁺[iparam,:] .+ 
                added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.J₀⁺
            composite_layer_lin.J̇₀⁻[iparam,:] .= composite_layer_lin.J̇₀⁻[iparam,:] .+ 
                composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_J̇₀⁻[iparam,:] .+ 
                composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ added_layer.J₀⁻
        end
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ added_layer.J₀⁻
    end
    # Batched multiplication between added and composite
    for iparam=1:Nparams 
        composite_layer_lin.Ṫ⁻⁻[iparam,:] = added_layer_lin.ap_ṫ⁻⁻[iparam,:] ⊠ composite_layer.T⁻⁻ .+
                                added_layer.t⁻⁻ ⊠ composite_layer_lin.Ṫ⁻⁻[iparam,:] 
        composite_layer_lin.Ṫ⁺⁺[iparam,:] = added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.T⁺⁺ .+
                                added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[iparam,:]
    end
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

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[1]
    if SFI
        #J₀⁺, J₀⁻ = similar(composite_layer.J₀⁺), similar(composite_layer.J₀⁺)
        #J₀⁻ = composite_layer.J₀⁻ .+ composite_layer.T⁻⁻ ⊠ (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻) 
        #J₀⁺ = added_layer.J₀⁺ .+ added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺ 
        for iparam=1:Nparams
            composite_layer_lin.J̇₀⁻[iparam,:] .= composite_layer_lin.J̇₀⁻[iparam,:] .+ 
                composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ 
                (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻) .+
                composite_layer.T⁻⁻ ⊠ 
                (added_layer_lin.ap_ṙ⁻⁺[iparam,:] ⊠ composite_layer.J₀⁺ .+ 
                added_layer.r⁻⁺ ⊠ composite_layer_lin.J̇₀⁺[iparam,:] .+ 
                added_layer_lin.ap_J̇₀⁻[iparam,:])
            composite_layer_lin.J̇₀⁺[iparam,:] .= added_layer_lin.ap_J̇₀⁺[iparam,:] .+ 
                added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.J₀⁺ .+
                added_layer.t⁺⁺ ⊠ composite_layer_lin.J̇₀⁺[iparam,:]  
        end
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ 
            (added_layer.r⁻⁺ ⊠ composite_layer.J₀⁺ .+ added_layer.J₀⁻)
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ composite_layer.J₀⁺         
    end

    # Batched multiplication between added and composite
    for iparam = 1:Nparams
        composite_layer_lin.Ṙ⁻⁺[iparam,:] = composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ added_layer.r⁻⁺ ⊠ composite_layer.T⁺⁺ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṙ⁻⁺[iparam,:] ⊠ composite_layer.T⁺⁺ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer.r⁻⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[iparam,:]
        composite_layer_lin.Ṙ⁺⁻[iparam,:] = added_layer_lin.ap_ṙ⁺⁻[iparam,:]
        composite_layer_lin.Ṫ⁺⁺[iparam,:] = added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.T⁺⁺ .+
                                    added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[iparam,:]
        composite_layer_lin.Ṫ⁻⁻[iparam,:] = composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ added_layer.t⁻⁻ .+
                                    composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[iparam,:]  
    end
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
                                composite_layer_lin::CompositeLayerLin{FT}, 
                                added_layer::AddedLayer{FT}, 
                                added_layer_lin::AddedLayerLin{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[1]
    if SFI
        for iparam=1:Nparams
            composite_layer_lin.J̇₀⁺[iparam,:] .= added_layer_lin.ap_J̇₀⁺[iparam,:] .+ 
                added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ 
                (composite_layer.J₀⁺ .+ composite_layer.R⁺⁻ ⊠ added_layer.J₀⁻) .+
                added_layer.t⁺⁺ ⊠ 
                (composite_layer_lin.J̇₀⁺[iparam,:] .+ 
                composite_layer_lin.Ṙ⁺⁻[iparam,:] ⊠ added_layer.J₀⁻ .+ 
                composite_layer.R⁺⁻ ⊠ added_layer_lin.ap_J̇₀⁻[iparam,:])
            composite_layer_lin.J̇₀⁻[iparam,:] .= composite_layer_lin.J̇₀⁻[iparam,:] .+ 
                composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ added_layer.J₀⁻ .+
                composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_J̇₀⁻[iparam,:] 
        end
        composite_layer.J₀⁺ .= added_layer.J₀⁺ .+ 
            added_layer.t⁺⁺ ⊠ 
            (composite_layer.J₀⁺ .+ composite_layer.R⁺⁻ ⊠ added_layer.J₀⁻)
        composite_layer.J₀⁻ .= composite_layer.J₀⁻ .+ 
            composite_layer.T⁻⁻ ⊠ added_layer.J₀⁻    
    end

    # Batched multiplication between added and composite
    for iparam=1:Nparams
        composite_layer_lin.Ṫ⁺⁺[iparam,:] = added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.T⁺⁺ .+
                                        added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṫ⁺⁺[iparam,:]
        composite_layer_lin.Ṫ⁻⁻[iparam,:] = composite_layer_lin.Ṫ⁻⁻[iparam,:] ⊠ added_layer.t⁻⁻ .+
                                        composite_layer.T⁻⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[iparam,:]
        composite_layer_lin.Ṙ⁺⁻[iparam,:] = added_layer_lin.ap_ṫ⁺⁺[iparam,:] ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻ .+
                                            added_layer.t⁺⁺ ⊠ composite_layer_lin.Ṙ⁺⁻[iparam,:] ⊠ added_layer.t⁻⁻ .+
                                            added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer_lin.ap_ṫ⁻⁻[iparam,:]
    end
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ ⊠ added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ ⊠ composite_layer.R⁺⁻ ⊠ added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT}, 
                                composite_layer_lin::CompositeLayerLin{FT}, 
                                added_layer::AddedLayer{FT}, 
                                added_layer_lin::AddedLayerLin{FT}, 
                                I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer #these are aliases to the respective struct elements  
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻ = composite_layer #these are aliases to the respective struct elements 
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺ = added_layer_lin #these are aliases to the respective struct elements  
    @unpack Ṙ⁻⁺, Ṙ⁺⁻, Ṫ⁺⁺, Ṫ⁻⁻, J̇₀⁺, J̇₀⁻ = composite_layer_lin #these are aliases to the respective struct elements 
    
    Nparams = size(composite_layer_lin.Ṫ⁻⁻)[1]
    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    tmp_inv = similar(t⁺⁺)
    tmp_inv_lin = similar(Ṫ⁺⁺)
    T01_inv_lin = similar(Ṫ⁺⁺)
    tmpṘ⁻⁺ = similar(Ṫ⁺⁺)
    tmpṪ⁻⁻ = similar(Ṫ⁺⁺)
    tmpap_J̇₀⁻ = similar(Ṫ⁺⁺)
    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹`
    @timeit "interaction inv1" batch_inv!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻) 
    # Temporary arrays:
    # T₁₂(I-R₀₁R₂₁)⁻¹
    T01_inv = T⁻⁻ ⊠ tmp_inv;
    for iparam=1:Nparams
        tmp_inv_lin[iparam,:] = tmp_inv ⊠ (ap_ṙ⁻⁺[iparam,:] ⊠ R⁺⁻ .+ r⁻⁺ ⊠ Ṙ⁺⁻[iparam,:]) ⊠ tmp_inv
        T01_inv_lin[iparam,:] = Ṫ⁻⁻[iparam,:] ⊠ tmp_inv .+ T⁻⁻ ⊠ tmp_inv_lin[iparam,:]
        # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
        tmpṘ⁻⁺[iparam,:] = Ṙ⁻⁺[iparam,:] .+ 
                        T01_inv_lin[iparam,:] ⊠ r⁻⁺ ⊠ T⁺⁺ .+
                        T01_inv ⊠ (ap_ṙ⁻⁺[iparam,:] ⊠ T⁺⁺ .+ r⁻⁺ ⊠ Ṫ⁺⁺[iparam])
    
        # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
        tmpṪ⁻⁻[iparam,:] = T01_inv_lin[iparam,:] ⊠ t⁻⁻ .+ T01_inv ⊠ ap_ṫ⁻⁻[iparam,:] 
    end
    
    if SFI
        #J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻)
        tmpJ₀⁻ = J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.J₀⁻) 
        for iparam=1:Nparams
            tmpap_J̇₀⁻[iparam,:] = ap_J̇₀⁻[iparam,:] .+ 
                T01_inv_lin[iparam,:] ⊠ (r⁻⁺ ⊠ J₀⁺ .+ added_layer.J₀⁻) .+
                T01_inv ⊠ (ap_ṙ⁻⁺[iparam,:] ⊠ J₀⁺ .+ r⁻⁺ ⊠ ap_J̇₀⁺[iparam,:] .+ added_layer_lin.ap_J̇₀⁻[iparam])  
        end
    end 

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    tmpR⁻⁺ = R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺
    
    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    tmpT⁻⁻ = T01_inv ⊠ t⁻⁻ 

    # Repeating for mirror-reflected directions
    T21_inv_lin = similar(Ṫ⁺⁺)
    tmpṘ⁺⁻ = similar(Ṫ⁺⁺)
    tmpṪ⁺⁺ = similar(Ṫ⁺⁺)
    tmpap_J̇₀⁺ = similar(Ṫ⁺⁺)
    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹`
    @timeit "interaction inv2" batch_inv!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺) 
    # T₂₁(I-R₀₁R₂₁)⁻¹
    T21_inv = t⁺⁺ ⊠ tmp_inv
    for iparam=1:Nparams
        tmp_inv_lin[iparam,:] = tmp_inv ⊠ (R⁺⁻ ⊠ ap_ṙ⁻⁺[iparam,:] .+ Ṙ⁺⁻[iparam,:] ⊠ r⁻⁺) ⊠ tmp_inv
        T21_inv_lin[iparam,:] = ap_ṫ⁺⁺[iparam,:] ⊠ tmp_inv .+ t⁺⁺ ⊠ tmp_inv_lin[iparam,:]

        # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
        tmpṪ⁺⁺[iparam,:] = T21_inv_lin[iparam,:] ⊠ T⁺⁺ .+ T21_inv ⊠ Ṫ⁺⁺[iparam,:] 
    
        # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
        tmpṘ⁺⁻[iparam,:] = ap_ṙ⁺⁻[iparam,:] .+ T21_inv_lin[iparam,:] ⊠ R⁺⁻ ⊠ t⁻⁻ .+ 
                                    T21_inv ⊠ (Ṙ⁺⁻[iparam,:] ⊠ t⁻⁻ .+ R⁺⁻ ⊠ ap_ṫ⁻⁻[iparam,:])  
    end
    if SFI
        for iparam=1:Nparams
            tmpap_J̇₀⁺[iparam,:] = added_layer_lin.ap_J̇₀⁺[iparam,:] .+ 
                T21_inv_lin[iparam,:] ⊠ (J₀⁺ .+ R⁺⁻ ⊠ added_layer.J₀⁻) .+
                T21_inv ⊠ (ap_J̇₀⁺[iparam,:] .+ 
                    Ṙ⁺⁻[iparam,:] ⊠ added_layer.J₀⁻ .+ 
                    R⁺⁻ ⊠ added_layer_lin.ap_J̇₀⁻[iparam,:])
        end
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
        
        composite_layer.ap_J̇₀⁺[:] = tmpap_J̇₀⁺
        composite_layer.ap_J̇₀⁻[:] = tmpap_J̇₀⁻
        
    end
    composite_layer.R⁺⁻[:] = tmpR⁺⁻
    composite_layer.T⁻⁻[:] = tmpT⁻⁻
    composite_layer.R⁻⁺[:] = tmpR⁻⁺
    composite_layer.T⁺⁺[:] = tmpT⁺⁺

    composite_layer.Ṙ⁺⁻[:] = tmpṘ⁺⁻
    composite_layer.Ṫ⁻⁻[:] = tmpṪ⁻⁻
    composite_layer.Ṙ⁻⁺[:] = tmpṘ⁻⁺
    composite_layer.Ṫ⁺⁺[:] = tmpṪ⁺⁺
end

"Compute interaction between composite and added layers"
function interaction!(scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        composite_layer_lin::CompositeLayerLin{FT}, 
                        added_layer::AddedLayer{FT},
                        added_layer_lin::AddedLayerLin{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    
    
    #@show A1[1,1,1], A2[1,1,1]
    interaction_helper!(scattering_interface, SFI, 
        composite_layer, composite_layer_lin, 
        added_layer, added_layer_lin, I_static)
    #A1 = Array(composite_layer.J₀⁻)
    #A2 = Array(composite_layer.J₀⁺)
    #@show A1[1,1,1], A2[1,1,1]
    synchronize_if_gpu()
end