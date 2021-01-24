
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function rt_doubling_helper!(pol_type,ndoubl::Int, 
                            added_layer::AddedLayer{FT},
                            I_static::AbstractArray{FT}, 
                            architecture) where {FT}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    dev = devi(architecture)
    # @show FT
    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing

    # Used to store `inv(I - r⁻⁺ * r⁻⁺) * t⁺⁺`
    tmp_inv  = similar(t⁺⁺)
    # tmp_inv2 = similar(t⁺⁺)
    # tmp_inv2[:] = r⁻⁺ ⊠ r⁻⁺
    # Loop over each step
    for n = 1:ndoubl
        batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ r⁻⁺, t⁺⁺)   
        added_layer.r⁻⁺[:]  = r⁻⁺ + (t⁺⁺ ⊠ r⁻⁺ ⊠ tmp_inv)
        added_layer.t⁺⁺[:]  = t⁺⁺ ⊠ tmp_inv
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    synchronize()
    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)
    return nothing 
end

function rt_doubling!(pol_type,ndoubl::Int, 
                    added_layer::AddedLayer{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT}

    rt_doubling_helper!(pol_type, ndoubl, added_layer, I_static, architecture)
    synchronize()
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ - 1, n_stokes)
    j = mod(jμ - 1, n_stokes)

    if (i >= 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
    end
    
    if ((i <= 1) & (j <= 1)) | ((i >= 2) & (j >= 2))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        return nothing
    else 
        applyD_kernel! = apply_D!(KernelAbstractions.CUDADevice())
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        wait(KernelAbstractions.CUDADevice(), event);
        synchronize();
        return nothing
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        return nothing
    else 
        applyD_kernel! = apply_D!(KernelAbstractions.CPU())
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        wait(KernelAbstractions.CPU(), event);
        synchronize();
        return nothing
    end
end