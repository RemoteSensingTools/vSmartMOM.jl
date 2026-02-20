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

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, temp1, temp2, temp1_ptr, temp2_ptr,
       dbl_gp_refl, dbl_j₁⁺, dbl_j₁⁻) = added_layer
    dev = devi(architecture)

    ndoubl == 0 && return nothing
    
    @timeit "doubling_allocs" begin
    tt⁺⁺_gp_refl = dbl_gp_refl === nothing ? similar(t⁺⁺) : dbl_gp_refl
    j₁⁺ = dbl_j₁⁺ === nothing ? similar(j₀⁺) : dbl_j₁⁺
    j₁⁻ = dbl_j₁⁻ === nothing ? similar(j₀⁻) : dbl_j₁⁻
    end
    #temp = similar(t⁺⁺)
    # Pointers to avoid memory allocation in CUBLAS routines
    #@timeit "Pointers" gp_ptrs   = CUBLAS.unsafe_strided_batch(gp_refl)
    #@timeit "Pointers" temp_ptrs = CUBLAS.unsafe_strided_batch(temp)
    # Loop over number of doublings
    for n = 1:ndoubl
        @timeit "Batch Inv Doubling" compute_geometric_progression!(temp1, tt⁺⁺_gp_refl, r⁻⁺, t⁺⁺, I_static, temp2, temp1_ptr, temp2_ptr)
        @timeit "source_update" doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt⁺⁺_gp_refl, expk)
        @timeit "rt_update" doubling_rt_update!(r⁻⁺, t⁺⁺, tt⁺⁺_gp_refl, expk)
    end
    @timeit "sync_doubling" synchronize_if_gpu()

    @timeit "apply_D_matrix" begin
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    apply_D_matrix_SFI!(pol_type.n, j₀⁻)
    end
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

@inline function apply_D_matrix!(n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
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


@inline function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
    #wait(device, event);
    synchronize_if_gpu();
    nothing
end
