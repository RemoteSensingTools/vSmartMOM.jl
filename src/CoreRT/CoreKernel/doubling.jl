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
                          added_layer::M,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT,M}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, temp1, temp2, temp1_ptr,temp2_ptr = added_layer
    #@show typeof(expk), typeof(I_static)
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    #gp_refl      = temp1# similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
    #temp = similar(t⁺⁺)
    # Dummy for source 
    j₁⁺ = similar(j₀⁺)
    # Dummy for J
    j₁⁻  = similar(j₀⁻)
    #temp = similar(t⁺⁺)
    # Pointers to avoid memory allocation in CUBLAS routines
    #@timeit "Pointers" gp_ptrs   = CUBLAS.unsafe_strided_batch(gp_refl)
    #@timeit "Pointers" temp_ptrs = CUBLAS.unsafe_strided_batch(temp)
    # Loop over number of doublings
    for n = 1:ndoubl
        temp2 .= I_static .- r⁻⁺ ⊠ r⁻⁺
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        #@show typeof(gp_refl), typeof(I_static), typeof(I_static .- r⁻⁺), typeof(j₁⁺)
        @timeit "Batch Inv Doubling" batch_inv!(temp1, temp2,temp1_ptr, temp2_ptr)
        tt⁺⁺_gp_refl .= t⁺⁺ ⊠ temp1

        # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
        @inbounds @views j₁⁺[:,1,:] .= j₀⁺[:,1,:] .* expk'

        # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
        @inbounds @views j₁⁻[:,1,:] .= j₀⁻[:,1,:] .* expk'

        # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
        j₀⁻ .= j₀⁻ + (tt⁺⁺_gp_refl ⊠ (j₁⁻ + r⁻⁺ ⊠ j₀⁺)) 

        # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
        j₀⁺  .= j₁⁺ + (tt⁺⁺_gp_refl ⊠ (j₀⁺ + r⁻⁺ ⊠ j₁⁻))
        expk .= expk.^2
    
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺  .= r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺  .= tt⁺⁺_gp_refl ⊠ t⁺⁺
    end
    synchronize_if_gpu()

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻
    apply_D_matrix_SFI!(pol_type.n, j₀⁻)
#    CUBLAS.unsafe_free!(temp_ptrs);
#    CUBLAS.unsafe_free!(gp_ptrs);
    return nothing 
end

"""
    doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)

Double the elemental layer `ndoubl` times to build the full homogeneous-layer
reflectance, transmission, and source matrices stored in `added_layer`.

Delegates to [`doubling_helper!`](@ref), which iteratively applies the
adding equations (see Eqs. 8 in the Raman paper draft) and then restores
the `D`-matrix symmetry.  After completion a GPU synchronisation barrier
is issued.

# Arguments
- `pol_type`: polarization type (determines `D`-matrix structure)
- `SFI`: whether Source Function Integration is active
- `expk`: `exp(-dτ/μ₀)` attenuation factor (doubled each iteration)
- `ndoubl::Int`: number of doubling steps
- `added_layer`: [`AddedLayer`](@ref) whose `r`, `t`, `j` fields are updated in-place
- `I_static`: pre-allocated batched identity matrix
- `architecture`: CPU or GPU selector
"""
function doubling!(pol_type, SFI, 
                    expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer{M},#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT,M}

    doubling_helper!(pol_type, SFI, 
                    expk, ndoubl, added_layer, I_static, architecture)
    synchronize_if_gpu()
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
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

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    if (i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻ .= r⁻⁺
        t⁻⁻ .= t⁺⁺    
        return nothing
    else 
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end


function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
    #wait(device, event);
    synchronize_if_gpu();
    nothing
end
