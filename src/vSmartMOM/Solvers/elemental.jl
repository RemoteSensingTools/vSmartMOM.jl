# <<Suniti>> it would be helpful to comment through this file (and remove unwanted lines). Thanks! 

"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,#{FT2,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ_λ::AbstractArray{FT,1},   # ϖ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ::FT,                      # ϖ: single scattering albedo of elemental layer (no trace gas included)
                            Z⁺⁺::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺::AbstractArray{FT,2}, 
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::AddedLayer{FT}, 
                            I_static,
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
    arr_type = array_type(architecture)
    # @show FT
    # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    

    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))

    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))
    I₀_NquadN = arr_type(zeros(FT,size(qp_μN,1)));
    i_start   = pol_type.n*(iμ₀-1) + 1 
    i_end     = pol_type.n*iμ₀
    I₀_NquadN[iμ₀Nstart:i_end] = pol_type.I₀

    device = devi(architecture)

    # If in scattering mode:
    if scatter
   
        NquadN = length(qp_μN)

        # Needs explanation still, different weigths:
        wct0  = m == 0 ? FT(0.50) * ϖ * dτ     : FT(0.25) * ϖ * dτ
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        wct   = m == 0 ? FT(0.50) * ϖ * wt_μN  : FT(0.25) * ϖ * wt_μN
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4

        # Get the diagonal matrices first
        d_qp  = Diagonal(1 ./ qp_μN)
        d_wct = Diagonal(wct)

        # Calculate r⁻⁺ and t⁺⁺
        
        # Version 1: no absorption in batch mode (like before), need to separate these modes
        if maximum(dτ_λ) < 0.0001 
            #@show "Chose simple elemental"
            #@show typeof(τ_sum)
            r⁻⁺[:,:,:] .= d_qp * Z⁻⁺ * (d_wct * dτ)
            t⁺⁺[:,:,:] .= I_static - (d_qp * ((I_static - Z⁺⁺ * d_wct) * dτ))
            if SFI
                # Reminder: Add equation here what it does
                expk = exp.(-τ_sum/qp_μ[iμ₀])

                J₀⁺[:,1,:] .= ((d_qp * Z⁺⁺ * I₀_NquadN) * wct0) .* expk'
                J₀⁻[:,1,:] .= ((d_qp * Z⁻⁺ * I₀_NquadN) * wct0) .* expk'
              
            end
            #ii = pol_type.n*(iμ0-1)+1
            #@show 'A',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
    
        else
            #@show "Chose accurate elemental",  maximum(dτ_λ)   
            #@show('B')
            # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
            # This is not yet GPU ready as it has element wise operations (should work for CPU)
            
            kernel! = get_elem_rt!(device)
            event = kernel!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, qp_μN, wct, ndrange=size(r⁻⁺)); 
            wait(device, event)
            synchronize_if_gpu()

            if SFI
                kernel! = get_elem_rt_SFI!(device)
                event = kernel!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μN, ndoubl, wct02, pol_type.n, arr_type(pol_type.I₀), iμ₀, D, ndrange=size(J₀⁺))
                wait(device, event)
            end
            #ii = pol_type.n*(iμ0-1)+1
            #@show 'B',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
            synchronize_if_gpu()
        end
        kernel! = apply_D_elemental!(device)
        event = kernel!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        wait(device, event)
        synchronize_if_gpu()
        if SFI
            kernel! = apply_D_elemental_SFI!(device)
            event = kernel!(ndoubl, pol_type.n, J₀⁻, ndrange=size(J₀⁻));
            wait(device, event)
            synchronize_if_gpu()
        end
        #@show(r⁻⁺[1,1,1], t⁺⁺[1,1,1], J₀⁻[1,1], J₀⁺[1,1])
        
    else 
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
        #if SFI
        #    J₀⁺ = I₀_NquadN * exp(-τ_sum/qp_μ[iμ0])
        #end
    end    

    @pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

@kernel function get_elem_rt!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, qp_μN, wct2)
    i, j, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    #D = arr_type(Diagonal(repeat(pol_type.D, size(qp_μ4)[1]/pol_type.n))) #Suniti, #Chr: needs to be outside if using GPU
    if (wct2[j]>1.e-8) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = ϖ_λ[n] * Z⁻⁺[i,j] * (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))) * (wct2[j]) 
                    
        if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i])*(1 + ϖ_λ[n] * Z⁺⁺[i,i] * (dτ_λ[n] / qp_μN[i]) * wct2[i])
            else
                t⁺⁺[i,j,n] = 0.0
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = ϖ_λ[n] * Z⁺⁺[i,j] * (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j])) * wct2[j]
        end
    else
        r⁻⁺[i,j,n] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i]) #Suniti
        else
            t⁺⁺[i,j,n] = 0.0
        end
    end
    
end

@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μN, ndoubl, wct02, nStokes ,I₀, iμ0, D)
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    J₀⁺[i, 1, n]=0
    J₀⁻[i, 1, n]=0

    
    Z⁺⁺_I₀ = FT(0.0);
    Z⁻⁺_I₀ = FT(0.0);
    for ii = i_start:i_end
        Z⁺⁺_I₀ += Z⁺⁺[i,ii] * I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii] * I₀[ii-i_start+1] 
    end

    if (i>=i_start) && (i<=i_end)
        ctr = i-i_start+1
        #J₀⁺[i,n] = exp(-dτ_λ[n] / qp_μ4[i]) * pol_type.I₀[ctr]
        # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
        #J₀⁺[i, 1, n] = testCF * ϖ_λ[n] * (Z⁺⁺[i,i_start:i_end]'*pol_type.I₀) * (dτ_λ[n] / qp_μN[i]) * exp.(-dτ_λ[n] / qp_μN[i])
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (dτ_λ[n] / qp_μN[i]) * exp(-dτ_λ[n] / qp_μN[i])
    else
        #J₀⁺[i, 1, n] = testCF * ϖ_λ[n] * (Z⁺⁺[i,i_start:i_end]'*pol_type.I₀) * (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
    end
    # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
    #J₀⁻[i, 1, n] = ϖ_λ[n] * (Z⁻⁺[i,i_start:i_end]'*pol_type.I₀) * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * (1 - exp.(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))

    #J₀⁺[i, 1, n] *= testCF * exp(-τ_sum[n]/qp_μN[i_start])
    #J₀⁻[i, 1, n] *= testCF * exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #Suniti: define D = Diagonal{1,1,-1,-1,...Nquad times}
    end        
end

@kernel function apply_D_elemental!(ndoubl, pol_n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    i, j, n = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod(i, pol_n) #Suniti
        jj = mod(j, pol_n) #Suniti
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) #Suniti
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] #Suniti: added - sign
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] #Suniti: added - sign
        end
    else
        if mod(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
        end 
    end
end

@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if mod(i, pol_n) > 2
            J₀⁻[i, 1, n] = - J₀⁻[i, 1, n]
        end 
    end
end

#function elemental!(pol_type, SFI, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, 
#                              ndoubl, scatter, qp_μ, wt_μ, 
#                              added_layer::AddedLayer{FT}, 
#                              I_static,
#                              arr_type,
#                              architecture) where {FT}
#    
#    elemental_helper!(pol_type, SFI, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, added_layer, I_static, arr_type, architecture)
#    synchronize_if_gpu()
#end