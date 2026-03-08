#=
 
This file contains RT doubling-related functions
 
=#

"""
    $(FUNCTIONNAME)(pol_type, SFI, expk, ndoubl::Int, added_layer::AddedLayer, I_static::AbstractArray{FT}, 
                    architecture) where {FT}

Compute homogenous layer matrices from its elemental layer using Doubling 
"""
function doubling_helper!(pol_type, 
                          SFI, 
                          expk, expk_lin,
                          τ_sum, τ̇_sum, 
                          ndoubl::Int, 
                          #AMF,
                          quad_points::QuadPoints{FT}, 
                          added_layer::AddedLayer,
                          added_layer_lin::AddedLayerLin,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT}

    @unpack μ₀ = quad_points
    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺, ap_J̇₀⁺, ap_J̇₀⁻ = added_layer_lin
    # Device architecture
    dev = devi(architecture)
    arr_type = array_type(architecture)
    Nparams = size(expk_lin,1)
    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
    gp_refl_lin       = arr_type(zeros(Nparams, size(t⁺⁺)[1], size(t⁺⁺)[2], size(t⁺⁺)[3]))
    tt⁺⁺_gp_refl_lin  = arr_type(zeros(Nparams, size(t⁺⁺)[1], size(t⁺⁺)[2], size(t⁺⁺)[3]))
    if SFI
        # Dummy for source 
        J₁⁺ = similar(J₀⁺)
        ap_J̇₁⁺ = similar(ap_J̇₀⁺)
        # Dummy for J
        J₁⁻ = similar(J₀⁻)
        ap_J̇₁⁻ = similar(ap_J̇₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        #@show n, r⁻⁺, any(isnan, r⁻⁺)
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        #@show n, r⁻⁺, any(isnan, r⁻⁺)
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        #@show n, gp_refl, any(isnan, gp_refl)
        #@show n, t⁺⁺, any(isnan, t⁺⁺)
        
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl
        #@show n, tt⁺⁺_gp_refl, any(isnan, tt⁺⁺_gp_refl)
        
        for iparam = 1:Nparams
            @views gp_refl_lin[iparam,:,:,:] .= gp_refl ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ r⁻⁺ .+ r⁻⁺ ⊠ ap_ṙ⁻⁺[iparam,:,:,:]) ⊠ gp_refl 
            @views tt⁺⁺_gp_refl_lin[iparam,:,:,:] .= ap_ṫ⁺⁺[iparam,:,:,:] ⊠ gp_refl .+ t⁺⁺ ⊠ gp_refl_lin[iparam,:,:,:]
        end
        if SFI
            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            @views J₁⁺[:,1,:] = J₀⁺[:,1,:] .* reshape(expk, 1, :)
            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            @views J₁⁻[:,1,:] = J₀⁻[:,1,:] .* reshape(expk, 1, :)
            for iparam = 1:Nparams
                #if iparam == 1
                    @views ap_J̇₁⁺[iparam,:,1,:] .= ap_J̇₀⁺[iparam,:,1,:] .* reshape(expk, 1, :) .+ J₀⁺[:,1,:] .* reshape(expk_lin[iparam,:], 1, :)        
                    @views ap_J̇₁⁻[iparam,:,1,:] .= ap_J̇₀⁻[iparam,:,1,:] .* reshape(expk, 1, :) .+ J₀⁻[:,1,:] .* reshape(expk_lin[iparam,:], 1, :)
                    
                    @views expk_lin .= 2* reshape(expk, 1, length(expk)) .* expk_lin
                #else
                #    @views ap_J̇₁⁺[iparam,:,1,:] .= ap_J̇₀⁺[iparam,:,1,:] .* reshape(expk, 1, :)         
                #    @views ap_J̇₁⁻[iparam,:,1,:] .= ap_J̇₀⁻[iparam,:,1,:] .* reshape(expk, 1, :) 
                #end
                @views ap_J̇₀⁻[iparam,:,:,:] .= ap_J̇₀⁻[iparam,:,:,:] .+ 
                        (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (J₁⁻ .+ r⁻⁺ ⊠ J₀⁺)) .+
                        (tt⁺⁺_gp_refl ⊠ (ap_J̇₁⁻[iparam,:,:,:] .+ ap_ṙ⁻⁺[iparam,:,:,:] ⊠ J₀⁺ .+ r⁻⁺ ⊠ ap_J̇₀⁺[iparam,:,:,:]))  
                @views ap_J̇₀⁺[iparam,:,:,:] .= ap_J̇₁⁺[iparam,:,:,:] .+ 
                    (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (J₀⁺ .+ r⁻⁺ ⊠ J₁⁻)) .+
                    (tt⁺⁺_gp_refl ⊠ (ap_J̇₀⁺[iparam,:,:,:] .+ ap_ṙ⁻⁺[iparam, :,:,:] ⊠ J₁⁻ .+ r⁻⁺ ⊠ ap_J̇₁⁻[iparam, :,:,:]))
            end

            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁻[:] = J₀⁻ .+ (tt⁺⁺_gp_refl ⊠ (J₁⁻ .+ r⁻⁺ ⊠ J₀⁺)) 
            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁺[:] = J₁⁺ .+ (tt⁺⁺_gp_refl ⊠ (J₀⁺ .+ r⁻⁺ ⊠ J₁⁻))
            expk[:] = expk.^2
        end  

        for iparam = 1:Nparams
            ap_ṙ⁻⁺[iparam, :,:,:] .= ap_ṙ⁻⁺[iparam, :,:,:] .+ 
                        tt⁺⁺_gp_refl_lin[iparam, :,:,:] ⊠ r⁻⁺ ⊠ t⁺⁺ .+
                        tt⁺⁺_gp_refl ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ t⁺⁺ .+
                        r⁻⁺ ⊠ ap_ṫ⁺⁺[iparam, :,:,:])
            ap_ṫ⁺⁺[iparam, :,:,:]  = tt⁺⁺_gp_refl_lin[iparam, :,:,:] ⊠ t⁺⁺ .+ 
                        tt⁺⁺_gp_refl ⊠ ap_ṫ⁺⁺[iparam, :,:,:]
        end
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ .+ (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    
    # This needs to be moved to where the linearization of added layers with respect to all parameters is carried out because τ̇_sum is a function of all parameters and cannot be reduced to derivatives wrt just τ, ϖ and Z 
    # Move this out from elemental to after doubling (it is not necessary to consider this in elemental if Raman scattering is not involved)
    #if SFI
    #    J₀⁺[:, 1, :] .*= (exp.(-τ_sum[:]/μ₀))' #writing i_start:i_start to avoid scalar indexing errors with GPUArrays
    #    J₀⁻[:, 1, :] .*= (exp.(-τ_sum[:]/μ₀))'


    #    J̇₀⁺[1, :, 1, :] = J̇₀⁺[1, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' +
    #                        J₀⁺[:, 1, :] .* ((-1/μ₀) .* @view(τ̇_sum[1,:]))'
    #    J̇₀⁻[1, :, 1, :] = J̇₀⁻[1, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' +
    #                        J₀⁻[:, 1, :] .* ((-1/μ₀) .* @view(τ̇_sum[1,:]))'
    #    J̇₀⁺[2, :, 1, :] = J̇₀⁺[2, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
    #    J̇₀⁻[2, :, 1, :] = J̇₀⁻[2, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
    #    J̇₀⁺[3, :, 1, :] = J̇₀⁺[3, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
    #    J̇₀⁻[3, :, 1, :] = J̇₀⁻[3, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
    #end
    
    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻

    synchronize_if_gpu()

    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻,
                    added_layer_lin.ap_ṙ⁻⁺, added_layer_lin.ap_ṫ⁺⁺, added_layer_lin.ap_ṙ⁺⁻, added_layer_lin.ap_ṫ⁻⁻, architecture)

    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.J₀⁻, added_layer_lin.ap_J̇₀⁻, architecture)

    return nothing 

end

function doubling!(pol_type, SFI, expk, expk_lin,
                    τ_sum::AbstractArray,#{FT2,1}, #Suniti
                    τ̇_sum::AbstractArray,
                    ndoubl::Int, 
                    quad_points::QuadPoints,#{FT}, 
                    #AMF,
                    added_layer::AddedLayer,#{FT},
                    added_layer_lin::AddedLayerLin,
                    I_static::AbstractArray, 
                    architecture) #where {FT}

    doubling_helper!(pol_type, SFI, expk, expk_lin, 
        τ_sum, τ̇_sum,
        ndoubl, quad_points,
        #AMF,
        added_layer, added_layer_lin, I_static, architecture)
    synchronize_if_gpu()
end

# WARNING: make sure the linearized version does not clash with the Raman version
@kernel function apply_D!(n_stokes::Int,  
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if !(1<=i<=2) #(i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
        # Unroll the 1:3 indexing to avoid GPU kernel issues
        ṙ⁻⁺[1,iμ,jμ,n] = - ṙ⁻⁺[1,iμ, jμ,n]
        ṙ⁻⁺[2,iμ,jμ,n] = - ṙ⁻⁺[2,iμ, jμ,n]
        ṙ⁻⁺[3,iμ,jμ,n] = - ṙ⁻⁺[3,iμ, jμ,n]
    end
    
    #if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
    if (((1<=i<=2) & (1<=j<=2)) | (!(1<=i<=2) & !(1<=j<=2)))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
        # Unroll the 1:3 indexing to avoid GPU kernel issues
        ṙ⁺⁻[1,iμ,jμ,n] = ṙ⁻⁺[1,iμ,jμ,n]
        ṙ⁺⁻[2,iμ,jμ,n] = ṙ⁻⁺[2,iμ,jμ,n]
        ṙ⁺⁻[3,iμ,jμ,n] = ṙ⁻⁺[3,iμ,jμ,n]
        ṫ⁻⁻[1,iμ,jμ,n] = ṫ⁺⁺[1,iμ,jμ,n]
        ṫ⁻⁻[2,iμ,jμ,n] = ṫ⁺⁺[2,iμ,jμ,n]
        ṫ⁻⁻[3,iμ,jμ,n] = ṫ⁺⁺[3,iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
        # Unroll the 1:3 indexing to avoid GPU kernel issues
        ṙ⁺⁻[1,iμ,jμ,n] = - ṙ⁻⁺[1,iμ,jμ,n]
        ṙ⁺⁻[2,iμ,jμ,n] = - ṙ⁻⁺[2,iμ,jμ,n]
        ṙ⁺⁻[3,iμ,jμ,n] = - ṙ⁻⁺[3,iμ,jμ,n]
        ṫ⁻⁻[1,iμ,jμ,n] = - ṫ⁺⁺[1,iμ,jμ,n]
        ṫ⁻⁻[2,iμ,jμ,n] = - ṫ⁺⁺[2,iμ,jμ,n]
        ṫ⁻⁻[3,iμ,jμ,n] = - ṫ⁺⁺[3,iμ,jμ,n]
    end

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻, J̇₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    if !(1<=i<=2) #(i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
        # Unroll the 1:3 indexing to avoid GPU kernel issues
        J̇₀⁻[1, iμ, 1, n] = - J̇₀⁻[1, iμ, 1, n]
        J̇₀⁻[2, iμ, 1, n] = - J̇₀⁻[2, iμ, 1, n]
        J̇₀⁻[3, iμ, 1, n] = - J̇₀⁻[3, iμ, 1, n]
    end
end

function apply_D_matrix!(n_stokes::Int, 
        r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
        r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
        ṙ⁻⁺::AbstractArray{FT,4}, ṫ⁺⁺::AbstractArray{FT,4}, 
        ṙ⁺⁻::AbstractArray{FT,4}, ṫ⁻⁻::AbstractArray{FT,4}, architecture) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺  
        ṙ⁺⁻[:] = ṙ⁻⁺
        ṫ⁻⁻[:] = ṫ⁺⁺    
        return nothing
    else 
        device = devi(architecture)
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, 
                                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, 
                                ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻, 
                                ndrange=size(r⁻⁺));
        ##wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end

#=function apply_D_matrix!(n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        #wait(device, event);
        return nothing
    end
end=#

function apply_D_matrix_SFI!(n_stokes::Int, 
                    J₀⁻::AbstractArray{FT,3}, 
                    J̇₀⁻::AbstractArray{FT,4}, architecture) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture)
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, J̇₀⁻, ndrange=size(J₀⁻));
    ##wait(device, event);
    synchronize_if_gpu();
    nothing
end

#=
function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::Array{FT,3}) where {FT}
    
    n_stokes == 1 && return nothing

    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
    #wait(device, event);
    
    return nothing
end=#