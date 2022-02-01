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
                          ndoubl::Int, 
                          added_layer::AddedLayer,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
 
    if SFI
        # Dummy for source 
        J₁⁺ = similar(J₀⁺)
        # Dummy for J
        J₁⁻ = similar(J₀⁻)

        # Dummy for source 
        ieJ₁⁺ = similar(ieJ₀⁺)
        # Dummy for J
        ieJ₁⁻ = similar(ieJ₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl

        if SFI

            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            J₁⁺[:,1,:] = J₀⁺[:,1,:] .* expk'
            ieJ₁⁺[:,1,:] = ieJ₀⁺[:,1,:] .* expk'

            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            J₁⁻[:,1,:] = J₀⁻[:,1,:] .* expk'
            ieJ₁⁻[:,1,:] = ieJ₀⁻[:,1,:] .* expk'

            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.17 in Raman paper draft)
            ieJ₀⁻[:,1,n₁,n₀] = ieJ₀⁻[:,1,n₁,n₀] + (tt⁺⁺_gp_refl[:,:,n₁]  
                * (ieJ₁⁻[:,1,n₁,n₀] + ier⁻⁺[:,:,n₁,n₀] * J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,n₀] 
                + (ier⁻⁺[:,:,n₁,n₀]*r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁]*ier⁻⁺[:,:,n₁,n₀])*gp_refl[:,:,n₀]*(J₁⁻[:,1,n₀] + r⁻⁺[:,:,n₀]*J₀⁺[:,1,n₀]))) 
                + iet⁻⁻[:,:,n₁,n₀]*gp_refl[:,:,n₀]*(J₁⁻[:,1,n₀] + r⁻⁺[:,:,n₀]*J₀⁺[:,1,n₀])#TODO

            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.16 in Raman paper draft)
            ieJ₀⁺[:,1,n₁,n₀] = ieJ₁⁺[:,1,n₁,n₀] + (tt⁺⁺_gp_refl[:,:,n₁] 
                * (ieJ₀⁺[:,1,n₁,n₀] + r⁻⁺[:,:,n₁] * ieJ₁⁻[:,1,n₁,n₀] + ier⁻⁺[:,:,n₁,n₀]*J₁⁻[:,1,n₀] 
                + (r⁻⁺[:,:,n₁]*ier⁻⁺[:,:,n₁,n₀] + ier⁻⁺[:,:,n₁,n₀]*r⁻⁺[:,:,n₀]) * gp_refl[:,:,n₀] 
                * (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀]))) + iet⁺⁺[:,:,n₁,n₀] * gp_refl[:,:,n₀] * (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀])
            
            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁻[:] = J₀⁻ + (tt⁺⁺_gp_refl ⊠ (J₁⁻ + r⁻⁺ ⊠ J₀⁺)) 

            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁺[:] = J₁⁺ + (tt⁺⁺_gp_refl ⊠ (J₀⁺ + r⁻⁺ ⊠ J₁⁻))
             
            expk[:] = expk.^2
        end  
        # (see Eqs.12 in Raman paper draft)
        iet⁺⁺[:] = t⁺⁺⊠gp_refl ⊠ (iet⁺⁺ + (ier⁻⁺⊠r⁻⁺ + r⁻⁺⊠ier⁻⁺) ⊠ gp_refl ⊠ t⁺⁺) + iet⁺⁺⊠gp_refl⊠t⁺⁺

        # (see Eqs.14 in Raman paper draft)
        ier⁻⁺[:] = ier⁻⁺ + t⁺⁺⊠gp_refl ⊠ r⁻⁺ ⊠ (iet⁺⁺ + (ier⁻⁺⊠r⁻⁺ + r⁻⁺⊠ier⁻⁺) ⊠ gp_refl ⊠ t⁺⁺) + iet⁺⁺⊠gp_refl⊠r⁻⁺⊠t⁺⁺ + t⁺⁺⊠gp_refl⊠ier⁻⁺⊠t⁺⁺

        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻

    synchronize_if_gpu()

    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)

    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.J₀⁻)

    return nothing 

end

function doubling!(pol_type, SFI, expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer,#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT}

    doubling_helper!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    synchronize_if_gpu()
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n, n₀  = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if (i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ, n]
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ, n]
        ier⁻⁺[iμ,jμ,n,n₀] = - ier⁻⁺[iμ, jμ, n, n₀]
        ier⁻⁺[iμ,jμ,n,n₀] = - ier⁻⁺[iμ, jμ, n, n₀]
    end
    
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
        ier⁺⁻[iμ,jμ,n,n₀] = ier⁻⁺[iμ,jμ,n,n₀]
        iet⁻⁻[iμ,jμ,n,n₀] = iet⁺⁺[iμ,jμ,n,n₀]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
        ier⁺⁻[iμ,jμ,n,n₀] = - ier⁻⁺[iμ,jμ,n,n₀]
        iet⁻⁻[iμ,jμ,n,n₀] = - iet⁺⁺[iμ,jμ,n,n₀]
    end

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻, ieJ₀⁻)
    iμ, _, n, n₀ = @index(Global, NTuple)
    i = mod(iμ, n_stokes)

    if (i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
        ieJ₀⁻[iμ, 1, n, n₀] = - ieJ₀⁻[iμ, 1, n, n₀] 
    end
end

#Suniti: is it possible to  use the same kernel for the 3D elastic and 4D inelastic terms or do we need to call two different kernels separately? 
function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}, ier⁻⁺::CuArray{FT,4}, iet⁺⁺::CuArray{FT,4}, ier⁺⁻::CuArray{FT,4}, iet⁻⁻::CuArray{FT,4}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺    
        ier⁺⁻[:] = ier⁻⁺
        iet⁻⁻[:] = iet⁺⁺  

        return nothing
    else 
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=size(ier⁻⁺)); #Suniti: is it possible to  use the same kernel for the 3D elastic and 4D inelastic terms or do we need to call two different kernels separately? 
        wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}, ier⁻⁺::CuArray{FT,4}, iet⁺⁺::CuArray{FT,4}, ier⁺⁻::CuArray{FT,4}, iet⁻⁻::CuArray{FT,4}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        ier⁺⁻[:] = ier⁻⁺
        iet⁻⁻[:] = iet⁺⁺  
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=size(r⁻⁺));
        wait(device, event);
        return nothing
    end
end

function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::CuArray{FT,3}, ieJ₀⁻::CuArray{FT,4}) where {FT}

    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻)) #Suniti: how to do this so that ieJ₀⁻ can also be included?
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ieJ₀⁻, ndrange=size(J₀⁻));
    wait(device, event);
    synchronize();
    
    return nothing
end
    
function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::Array{FT,3}, ieJ₀⁻::CuArray{FT,4}) where {FT}
    
    n_stokes == 1 && return nothing

    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ieJ₀⁻, ndrange=size(J₀⁻));
    wait(device, event);
    
    return nothing
end