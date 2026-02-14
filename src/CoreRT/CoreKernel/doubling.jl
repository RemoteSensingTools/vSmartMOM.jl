#=
 
This file contains RT doubling-related functions
 
=#

"""
    $(FUNCTIONNAME)(pol_type, SFI, expk, ndoubl::Int, added_layer::AddedLayer, I_static::AbstractArray{FT}, 
                    architecture; ws=nothing) where {FT}

Compute homogenous layer matrices from its elemental layer using Doubling.
If `ws::RTWorkspace` is provided, uses preallocated temporaries (critical for GPU performance).
"""
function doubling_helper!(pol_type, 
                          SFI, 
                          expk, 
                          ndoubl::Int, 
                          added_layer::Union{AddedLayer,AddedLayerRS},
                          I_static::AbstractArray{FT}, 
                          architecture;
                          ws::Union{RTWorkspace, Nothing}=nothing) where {FT}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Use workspace temporaries if available; otherwise allocate (CPU fallback)
    if ws !== nothing
        gp_refl      = ws.gp_refl
        tt⁺⁺_gp_refl = ws.tt_gp
        J₁⁺          = ws.J₁⁺
        J₁⁻          = ws.J₁⁻
    else
        gp_refl      = similar(t⁺⁺)
        tt⁺⁺_gp_refl = similar(t⁺⁺)
        J₁⁺          = similar(J₀⁺)
        J₁⁻          = similar(J₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        if ws !== nothing
            ws.tmp3d_a .= I_static .- r⁻⁺ ⊠ r⁻⁺
            batch_inv!(gp_refl, ws.tmp3d_a, ws)
        else
            batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        end
        tt⁺⁺_gp_refl .= t⁺⁺ ⊠ gp_refl

        if SFI

            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            @views J₁⁺[:,1,:] .= J₀⁺[:,1,:] .* expk'

            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            @views J₁⁻[:,1,:] .= J₀⁻[:,1,:] .* expk'

            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁻[:] = J₀⁻ .+ (tt⁺⁺_gp_refl ⊠ (J₁⁻ .+ r⁻⁺ ⊠ J₀⁺)) 

            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁺[:] = J₁⁺ .+ (tt⁺⁺_gp_refl ⊠ (J₀⁺ .+ r⁻⁺ ⊠ J₁⁻))
            expk .= expk.^2
        end  

        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ .+ (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

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

"""
    doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture; ws=nothing)

Public interface: doubling with optional workspace preallocation.
"""
function doubling!(pol_type, SFI, expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer,#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture;
                    ws::Union{RTWorkspace, Nothing}=nothing) where {FT}

    doubling_helper!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture; ws=ws)
    synchronize_if_gpu()
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if !(1<=i<=2) #(i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
    end
    
    #if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
    if (((1<=i<=2) & (1<=j<=2)) | (!(1<=i<=2) & !(1<=j<=2)))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
    end

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    if !(1<=i<=2) #(i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺    
        return nothing
    else 
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
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

function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
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