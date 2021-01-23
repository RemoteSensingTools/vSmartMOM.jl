
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function rt_doubling_helper!(pol_type,ndoubl::Int, 
                            added_layer::AddedLayer,
                            D::AbstractArray{FT,3},
                            I_static::AbstractArray, 
                            architecture) where {FT}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    dev = devi(architecture)
    # # ToDo: Important output doubling applied to elemental layer, using same variables 
    # r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)

    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but 
    # using simplifications in eq. 29-32)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing

    # Used to store `inv(I - r⁻⁺ * r⁻⁺) * t⁺⁺`
    tmp_inv = similar(t⁺⁺)

    # Loop over each step
    for n = 1:ndoubl
        batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ r⁻⁺, t⁺⁺)   
        r⁻⁺[:]  = r⁻⁺ + (t⁺⁺ ⊠ r⁻⁺ ⊠ tmp_inv)
        t⁺⁺[:]  = t⁺⁺ ⊠ tmp_inv
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    event = apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    wait(event)
    # applyD_kernel! = apply_D!(dev)
    # @show pol_type.n
    # event = applyD_kernel!(r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=length(r⁻⁺));
    # wait(event)
    # synchronize()
    # @show r⁻⁺ ≈ r⁺⁻
    # @show t⁺⁺ ≈ t⁻⁻
    # 
    # r⁺⁻[:] = r⁻⁺
    # t⁻⁻[:] = t⁺⁺
    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # r⁻⁺[:] = D ⊠ r⁻⁺ ⊠ D

    # Using r⁺⁻ = Dr⁻⁺D
    # r⁺⁻[:] = D ⊠ r⁻⁺ ⊠ D
    
    # Using t⁻⁻ = Dt⁺⁺D
    # t⁻⁻[:] = D ⊠ t⁺⁺ ⊠ D

    return nothing 
end

function rt_doubling!(pol_type,ndoubl::Int, added_layer::AddedLayer,
                    D::AbstractArray{FT,3},
                    I_static::AbstractArray, 
                    architecture) where {FT}

    rt_doubling_helper!(pol_type, ndoubl, added_layer, D, I_static, architecture)
    synchronize()
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    if n_stokes == 1
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        i = mod(iμ - 1, n_stokes)
        j = mod(jμ - 1, n_stokes)
        # @show i,j
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
    @synchronize()
end

@kernel function apply_D!(r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    i, j, n = @index(Global, NTuple)
    r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
    t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}) where {FT}
    applyD_kernel! = apply_D!(KernelAbstractions.CUDADevice())
    applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
end

function apply_D_matrix!(r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}) where {FT}
    applyD_kernel! = apply_D!(KernelAbstractions.CPU())
    applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
end