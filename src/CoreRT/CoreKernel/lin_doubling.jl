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
                          expk, 
                          dexpk_fctr,
                          ndoubl::Int, 
                          added_layer::M,
                          lin_added_layer::M,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT,M}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer
    @unpack dr⁺⁻, dr⁻⁺, dt⁻⁻, dt⁺⁺, dj₀⁺, dj₀⁻ = lin_added_layer
    #@show typeof(expk), typeof(I_static)
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)

    dgp_refl      = similar(dt⁺⁺)
    dtt⁺⁺_gp_refl = similar(dt⁺⁺)
    
    # Dummy for source 
    j₁⁺ = similar(j₀⁺)
    dj₁⁺ = similar(dj₀⁺)
    # Dummy for J
    j₁⁻ = similar(j₀⁻)
    dj₁⁻ = similar(dj₀⁻)

    # Loop over number of doublings
    for n = 1:ndoubl
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        #@show typeof(gp_refl), typeof(I_static), typeof(I_static .- r⁻⁺), typeof(j₁⁺)
        @timeit "Batch Inv Doubling" batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl .= t⁺⁺ ⊠ gp_refl
        for ctr=1:3
            tmp = -(dr⁻⁺[ctr,:,:,:] ⊠ r⁻⁺ + r⁻⁺ ⊠ dr⁻⁺[ctr,:,:,:])
            dgp_refl[ctr,:,:,:] .= tmp ⊠ gp_refl ⊠ tmp
            dtt⁺⁺_gp_refl[ctr,:,:,:] .= dt⁺⁺[ctr,:,:,:] ⊠ gp_refl +
                                        t⁺⁺ ⊠ dgp_refl[ctr,:,:,:]
        end
        

        # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
        @inbounds @views j₁⁺[:,1,:] .= j₀⁺[:,1,:] .* expk'
        # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
        @inbounds @views j₁⁻[:,1,:] .= j₀⁻[:,1,:] .* expk'
        for ctr=1:3
            if (ctr==1)
                dj₁⁺[ctr,:,1,:] .= (dj₀⁺[ctr,:,1,:] + j₀⁺[:,1,:]*dexpk_fctr) .* expk'
                # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
                dj₁⁻[ctr,:,1,:] .= (dj₀⁻[ctr,:,1,:] + j₀⁻[:,1,:]*dexpk_fctr) .* expk'
            else
                dj₁⁺[ctr,:,1,:] .= (dj₀⁺[ctr,:,1,:]) .* expk'
                # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
                dj₁⁻[ctr,:,1,:] .= (dj₀⁻[ctr,:,1,:]) .* expk'
            end

            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            dj₀⁻[ctr,:,1,:] .= dj₀⁻[ctr,:,1,:] + 
                (dtt⁺⁺_gp_refl[ctr,:,:,:] ⊠ (j₁⁻ + r⁻⁺ ⊠ j₀⁺)) +
                (tt⁺⁺_gp_refl ⊠ (dj₁⁻[ctr,:,1,:] + dr⁻⁺[ctr,:,:,:] ⊠ j₀⁺ + r⁻⁺ ⊠ dj₀⁺[ctr,:,1,:]))
            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            dj₀⁺[ctr,:,1,:]  .= dj₁⁺[ctr,:,1,:] + 
                (dtt⁺⁺_gp_refl[ctr,:,:,:] ⊠ (j₀⁺ + r⁻⁺ ⊠ j₁⁻)) +
                (tt⁺⁺_gp_refl ⊠ (dj₀⁺[ctr,:,1,:] + dr⁻⁺[ctr,:,:,:] ⊠ j₁⁻ + r⁻⁺ ⊠ dj₁⁻[ctr,:,1,:]))
        end
    
        # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
        j₀⁻ .= j₀⁻ + (tt⁺⁺_gp_refl ⊠ (j₁⁻ + r⁻⁺ ⊠ j₀⁺)) 

        # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
        j₀⁺  .= j₁⁺ + (tt⁺⁺_gp_refl ⊠ (j₀⁺ + r⁻⁺ ⊠ j₁⁻))
        expk .= expk.^2
    
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺  .= r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺  .= tt⁺⁺_gp_refl ⊠ t⁺⁺
        for ctr=1:3
            # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
            dr⁻⁺[ctr,:,:,:]  .= dr⁻⁺[ctr,:,:,:] + 
                    (dtt⁺⁺_gp_refl[ctr,:,:,:] ⊠ r⁻⁺ ⊠ t⁺⁺) +
                    (tt⁺⁺_gp_refl ⊠ dr⁻⁺[ctr,:,:,:] ⊠ t⁺⁺) +
                    (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ dt⁺⁺[ctr,:,:,:])

            # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
            dt⁺⁺[ctr,:,:,:]  .= dtt⁺⁺_gp_refl[ctr,:,:,:] ⊠ t⁺⁺ +
                    tt⁺⁺_gp_refl ⊠ dt⁺⁺[ctr,:,:,:]
        end
    end
    synchronize_if_gpu()

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, dr⁻⁺, dt⁺⁺, dr⁺⁻, dt⁻⁻)

    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻
    apply_D_matrix_SFI!(pol_type.n, j₀⁻, dj₀⁻)

    return nothing 
end

function doubling!(pol_type, SFI, 
                    expk, dexpk_fctr,
                    ndoubl::Int, 
                    added_layer::AddedLayer{M},#{FT},
                    lin_added_layer::linAddedLayer{M},#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT,M}

    doubling_helper!(pol_type, SFI, 
                    expk, dexpk_fctr, 
                    ndoubl, 
                    added_layer, 
                    lin_added_layer,
                    I_static, architecture)
    synchronize_if_gpu()
end

@kernel function apply_D!(n_stokes::Int,  
        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
        dr⁻⁺, dt⁺⁺, dr⁺⁻, dt⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if (i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
    end
    
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
    end

    for ctr=1:3
        if (i > 2)
            dr⁻⁺[ctr,iμ,jμ,n] = - dr⁻⁺[ctr,iμ, jμ,n]
        end
        
        if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
            dr⁺⁻[ctr,iμ,jμ,n] = dr⁻⁺[ctr,iμ,jμ,n]
            dt⁻⁻[ctr,iμ,jμ,n] = dt⁺⁺[ctr,iμ,jμ,n]
        else
            dr⁺⁻[ctr,iμ,jμ,n] = - dr⁻⁺[ctr,iμ,jμ,n]
            dt⁻⁻[ctr,iμ,jμ,n] = - dt⁺⁺[ctr,iμ,jμ,n]
        end
    end

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻, dJ₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    if (i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
        dJ₀⁻[1:3, iμ, 1, n] .= - dJ₀⁻[1:3, iμ, 1, n] 
    end
end

function apply_D_matrix!(n_stokes::Int, 
    r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
    dr⁻⁺::AbstractArray{FT,4}, dt⁺⁺::AbstractArray{FT,4}, dr⁺⁻::AbstractArray{FT,4}, dt⁻⁻::AbstractArray{FT,4}) where {FT}
    if n_stokes == 1
        r⁺⁻ .= r⁻⁺
        t⁻⁻ .= t⁺⁺   
        dr⁺⁻ .= dr⁻⁺
        dt⁻⁻ .= dt⁺⁺   
        return nothing
    else 
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        applyD_kernel!(n_stokes, 
            r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
            dr⁻⁺, dt⁺⁺, dr⁺⁻, dt⁻⁻, 
            ndrange=size(r⁻⁺));
#        wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end


function apply_D_matrix_SFI!(n_stokes::Int, 
    J₀⁻::AbstractArray{FT,3},
    dJ₀⁻::AbstractArray{FT,4}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    applyD_kernel!(n_stokes, J₀⁻, dJ₀⁻, ndrange=size(J₀⁻));
   # wait(device, event);
    synchronize_if_gpu();
    nothing
end
