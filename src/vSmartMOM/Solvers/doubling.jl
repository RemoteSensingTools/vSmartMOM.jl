
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function doubling_helper!(pol_type, SFI, expk, ndoubl::Int, 
                            added_layer::AddedLayer{FT},
                            I_static::AbstractArray{FT}, 
                            architecture) where {FT}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    dev = devi(architecture)
    # @show FT
    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    #j0 = pol_type.n*(iμ0-1)+1
    # Used to store `inv(I - r⁻⁺ * r⁻⁺) * t⁺⁺`
    tmp_inv  = similar(t⁺⁺)
    # tmp_inv2 = similar(t⁺⁺)
    # tmp_inv2[:] = r⁻⁺ ⊠ r⁻⁺
    #expk = exp.(-dτ_λ / qp_μN[j0])
    #assign dimensions (Nquad*pol_type.n, Nquad*pol_type.n,size(dτ_λ)) to E
    @show(added_layer.t⁺⁺[1,1,1],t⁺⁺[1,1,1])
    # Loop over each step
    for n = 1:ndoubl
        batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ r⁻⁺, t⁺⁺)   
        added_layer.r⁻⁺[:]  = r⁻⁺ + (t⁺⁺ ⊠ r⁻⁺ ⊠ tmp_inv)
        added_layer.t⁺⁺[:]  = t⁺⁺ ⊠ tmp_inv
        
        if SFI
            #nλ = @index(Global, NTuple)
            for nλ = 1:size(added_layer.J₀⁺,2)
                #E[:,:,nλ]=arr_type(Diagonal(repeat(expk[nλ]*pol_type.I₀,size(Nquad)))) #TODO: Pass exp(-dτ_λ ./ qp_μN[j0]) as an   argument to doubling and square it with every doubling iteration
                added_layer.J₀⁺[:,nλ] = (expk[nλ]*I_static + added_layer.t⁺⁺[:,:,nλ]) * J₀⁺[:,nλ] + added_layer.r⁻⁺[:,:,nλ] * J₀⁻[:,nλ]
                #@show size()
                #_a = added_layer.t⁺⁺[:,:,nλ] * J₀⁺[:,nλ] 
                #_b = (I_static + expk[nλ] * added_layer.t⁺⁺[:,:,nλ]) * J₀⁻[:,nλ]
                #@show size(_a), size(_b)
                #_c = _a + _b
                added_layer.J₀⁻[:,nλ] = added_layer.t⁺⁺[:,:,nλ] * J₀⁺[:,nλ] + (I_static + expk[nλ] * added_layer.t⁺⁺[:,:,nλ]) * J₀⁻[:,nλ]
            end
            #dτ_λ = 2*dτ_λ
            expk = expk .* expk
        end
        #@show(n,added_layer.t⁺⁺[1,1,1],t⁺⁺[1,1,1])
    end
    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻
    ### synchronize()
    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻, added_layer.J₀⁺, added_layer.J₀⁻)
    @pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻
    return nothing 
end

function doubling!(pol_type, SFI, expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT}

    doubling_helper!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    ### synchronize() #Suniti: does this convert the added layer to the current layer, so that added_layer.M = M?
end

@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, J₀⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if (i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
        if (SFI)
            J₀⁻[iμ, n] = - J₀⁻[iμ, n] 
        end
    end
    
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
    end


    #Suniti: Alternative version:
    #D = arr_type(Diagonal(repeat(pol_type.D, size(qp_μ)[1])))
    #r⁻⁺ = D.*r⁻⁺
    #J₀⁻ = D.*J₀⁻
    #r⁺⁻ = D.*r⁻⁺.*D
    #t⁻⁻ = D.*t⁺⁺.*D
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}, J₀⁺::CuArray{FT,2}, J₀⁻::CuArray{FT,2}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        if SFI
            J₀⁻ = J₀⁺
        end
        return nothing
    else 
        applyD_kernel! = apply_D!(KernelAbstractions.CUDADevice())
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, J₀⁻, ndrange=size(r⁻⁺));
        wait(KernelAbstractions.CUDADevice(), event);
        synchronize();
        return nothing
    end
end

function apply_D_matrix!(n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}, J₀⁺::Array{FT,2}, J₀⁻::Array{FT,2}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        if SFI
            J₀⁻ = J₀⁺
        end
        return nothing
    else 
        applyD_kernel! = apply_D!(KernelAbstractions.CPU())
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, J₀⁻, ndrange=size(r⁻⁺));
        wait(KernelAbstractions.CPU(), event);
        ### synchronize();
        return nothing
    end
end