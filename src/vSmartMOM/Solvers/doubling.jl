# <<Suniti>> it would be helpful to comment through this file. Thanks! 

# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps
function doubling_helper!(pol_type, SFI, expk, ndoubl::Int, 
                          added_layer::AddedLayer,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer

    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing

    #j0 = pol_type.n*(iμ0-1)+1
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
 
    if SFI
        # Dummy for source 
        J₁⁺ = similar(J₀⁺)
        # Dummy for J
        J₁⁻ = similar(J₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl
        # J₁⁺[:,1,:] =  J₀⁺[:,1,:] .* expk'
        if SFI
            J₁⁺[:,1,:] = J₀⁺[:,1,:] .* expk'
            J₁⁻[:,1,:] = J₀⁻[:,1,:] .* expk'
            #@show size(J₀⁺), size(J₁⁺), size((t⁺⁺ ⊠ gp_refl ⊠ (J₀⁺ .+ r⁻⁺ ⊠ J₁⁻)))
            J₀⁻[:] = J₀⁻ + (tt⁺⁺_gp_refl ⊠ (J₁⁻ + r⁻⁺ ⊠ J₀⁺)) 
            J₀⁺[:] = J₁⁺ + (tt⁺⁺_gp_refl ⊠ (J₀⁺ + r⁻⁺ ⊠ J₁⁻))
            expk[:] = expk.^2
        end  
        r⁻⁺[:]  = r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
        #@show r⁻⁺[1:3:end,28,1]./(J₀⁻[1:3:end,1,1]*0.005)#, r⁻⁺[1,28,1]

    end
    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻
    synchronize_if_gpu()
    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)
    if SFI
        apply_D_matrix_SFI!(pol_type.n, added_layer.J₀⁻)
    end
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻
    return nothing 
end

function doubling!(pol_type, SFI, expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer,#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT}

    doubling_helper!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    ### synchronize() #Suniti: does this convert the added layer to the current layer, so that added_layer.M = M?
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


function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺    
        
        return nothing
    else 
        applyD_kernel! = apply_D!(KernelAbstractions.CUDADevice())
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        wait(KernelAbstractions.CUDADevice(), event);
        synchronize_if_gpu();
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
        return nothing
    end
end

function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::CuArray{FT,3}) where {FT}
    if n_stokes == 1
        return nothing
    else 
        applyD_kernel! = apply_D_SFI!(KernelAbstractions.CUDADevice())
        event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(KernelAbstractions.CUDADevice(), event);
        synchronize();
        return nothing
    end
end
    
function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::Array{FT,3}) where {FT}
    if n_stokes == 1
        return nothing
    else 
        applyD_kernel! = apply_D_SFI!(KernelAbstractions.CPU())
        event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(KernelAbstractions.CPU(), event);
        return nothing
    end
end