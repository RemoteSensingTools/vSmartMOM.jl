"Elemental single-scattering layer"
function elemental_helper!(pol_type, SFI, iμ0,
                            τ_sum::AbstractArray{FT,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1}, 
                            dτ::FT, 
                            ϖ_λ::AbstractArray{FT,1}, 
                            ϖ::FT, 
                            Z⁺⁺::AbstractArray{FT,2}, 
                            Z⁻⁺::AbstractArray{FT,2}, 
                            m::Int, 
                            ndoubl::Int, 
                            scatter, 
                            qp_μ::AbstractArray{FT,1}, 
                            wt_μ::AbstractArray{FT,1}, 
                            added_layer::AddedLayer{FT}, 
                            I_static,
                            arr_type,
                            architecture) where {FT}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    # @show FT
    # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ (can be renamed to t⁺⁺, etc)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

    # dτ: optical depth of elemental layer
    # ϖ: single scattering albedo of elemental layer
    # bb: thermal source function at the upper boundary of the elemental layer
    # m: fourier moment
    # n: layer of which this is an elemental
    # ndoubl: number of doubling computations needed to progress from the elemental layer 
    #         to the full homogeneous layer n
    # scatter: flag indicating scattering

    Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)

    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)

    if scatter
        qp_μN = arr_type(reduce(vcat, (fill.(qp_μ, [pol_type.n]))))
        wt_μN = arr_type(reduce(vcat, (fill.(wt_μ, [pol_type.n]))))
        #for i=1:length(qp_μN)
        #    @show(i, qp_μN[i])
        #end
        NquadN = length(qp_μN)

        wct = m == 0 ? FT(0.50) *ϖ * wt_μN  : FT(0.25) * ϖ * wt_μN
        wct2 = m == 0 ? wt_μN  : wt_μN / 2
        # wct = m==0 ? 0.50 * 1 .* wt_μ4  : 0.25 .* 1 .* wt_μ4

        # Get the diagonal matrices first
        d_qp  = Diagonal(arr_type(1 ./ qp_μN))
        d_wct = Diagonal(arr_type(wct))

        # Calculate r⁻⁺ and t⁺⁺
        
        # Version 1: no absorption in batch mode (like before), need to separate these modes
        if maximum(dτ_λ) < 0.0001 

            r⁻⁺[:,:,:] .= d_qp * Z⁻⁺ * (d_wct * dτ)
            t⁺⁺[:,:,:] .= I_static - (d_qp * ((I_static - Z⁺⁺ * d_wct) * dτ))
        else    
        # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
        # This is not yet GPU ready as it has element wise operations (should work for CPU)
            
            kernel! = get_elem_rt!(device)
            event = kernel!(r⁻⁺, r⁺⁻, t⁺⁺, t⁻⁻, J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μN, wct2, ndoubl, pol_type, SFI, iμ0, D, ndrange=size(r⁻⁺)); #Suniti: what does this do? 
            
            wait(device, event)
            ### synchronize() # Check for CUDA here, only use with GPU!
        end
        #@show(r⁻⁺[1,1,1], t⁺⁺[1,1,1], J₀⁻[1,1], J₀⁺[1,1])
        
    else 
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@show(t⁺⁺[1,1,1,], added_layer.t⁺⁺[1,1,1])
    @pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

@kernel function get_elem_rt!(r⁻⁺, r⁺⁻, t⁺⁺, t⁻⁻, J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μ4, wct2, ndoubl, pol_type, SFI, iμ0, D)
    i, j, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    #D = arr_type(Diagonal(repeat(pol_type.D, size(qp_μ4)[1]/pol_type.n))) #Suniti, #Chr: needs to be outside if using GPU
    if (wct2[j]>1.e-8) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = ϖ_λ[n] * Z⁻⁺[i,j] * (qp_μ4[j] / (qp_μ4[i] + qp_μ4[j])) * (1 - exp.(-dτ_λ[n] * ((1 / qp_μ4[i]) + (1 / qp_μ4[j])))) * (wct2[j]) 
                    
        if (qp_μ4[i] == qp_μ4[j])

            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μ4[i]) + ϖ_λ[n] * Z⁺⁺[i,i] * (dτ_λ[n] / qp_μ4[i]) * exp.(-dτ_λ[n] / qp_μ4[i]) * wct2[i]
            else
                t⁺⁺[i,j,n] = 0.0
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = ϖ_λ[n] * Z⁺⁺[i,j] * (qp_μ4[j] / (qp_μ4[i] - qp_μ4[j])) * (exp(-dτ_λ[n] / qp_μ4[i]) - exp(-dτ_λ[n] / qp_μ4[j])) * wct2[j]
        end
    else
        r⁻⁺[i,j,n] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μ4[i]) #Suniti
        else
            t⁺⁺[i,j,n] = 0.0
        end
    end
    if ndoubl < 1
        ii = mod(i, pol_type.n) #Suniti
        jj = mod(j, pol_type.n) #Suniti
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) #Suniti
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] #Suniti: added - sign
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] #Suniti: added - sign
        end
    else
        if mod(i, pol_type.n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
        end 
    end
    if j==1
        if SFI
            J₀⁺[i,n]=0
             J₀⁻[i,n]=0
            #for j0=iμ0:(iμ0+pol_type.n-1)
            j0  = pol_type.n*(iμ0-1) + 1 
            #vj0 = j0:j0+pol_type.n-1    
            #Suniti: define τ_sum to be the sum of the optical thicknesses of all distinct homogeneous layers above the layer under consideration
            for ctr=0:pol_type.n-1
                if (qp_μ4[i] == qp_μ4[j0])
                    # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
                    J₀⁺[i,n] += exp(-dτ_λ[n] / qp_μ4[i]) + ϖ_λ[n] * (Z⁺⁺[i,j0+ctr]*pol_type.I₀[ctr+1]) * (dτ_λ[n] / qp_μ4[i]) * exp.(-dτ_λ[n] / qp_μ4[i]) * exp(-τ_sum[n]/qp_μ4[j0])
                else        
                    # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
                    # (𝑖 ≠ 𝑗)
                    J₀⁺[i,n] += ϖ_λ[n] * (Z⁺⁺[i,j0+ctr]*pol_type.I₀[ctr+1]) * (qp_μ4[j0] / (qp_μ4[i] - qp_μ4[j0])) * (exp(-dτ_λ[n] / qp_μ4[i]) - exp(-dτ_λ[n] / qp_μ4[j0])) * exp(-τ_sum[n]/qp_μ4[j0])
                end
                # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
                J₀⁻[i,n] += ϖ_λ[n] * (Z⁻⁺[i,j0+ctr]*pol_type.I₀[ctr+1]) * (qp_μ4[j0] / (qp_μ4[i] + qp_μ4[j0])) * (1 - exp.(-dτ_λ[n] * ((1 / qp_μ4[i]) + (1 / qp_μ4[j0])))) *exp(-τ_sum[n]/qp_μ4[j0])#*I0[i%4] Suniti 
            end 
            if ndoubl >= 1
                J₀⁻[i,n] = D[i,i]*J₀⁻[i,n] #Suniti: define D = Diagonal{1,1,-1,-1,...Nquad times}
            end
        end
    end
end

@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μ4, ndoubl, pol_type, iμ0, D)
    i, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?

    J₀⁺[i,n]=0
    J₀⁻[i,n]=0
    i_start  = pol_type.n*(iμ0-1) + 1 
    i_end    = pol_type.n*iμ0
    if (i>=i_start) && (i<=i_end)
        ctr = i-i_start+1
        J₀⁺[i,n] = exp(-dτ_λ[n] / qp_μ4[i]) * pol_type.I₀[ctr]
        # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
        J₀⁺[i,n] += ϖ_λ[n] * (Z⁺⁺[i,i_start:i_end]*pol_type.I₀) * (dτ_λ[n] / qp_μ4[i]) * exp.(-dτ_λ[n] / qp_μ4[i])
    else
        J₀⁺[i,n] += ϖ_λ[n] * (Z⁺⁺[i,i_start:i_end]*pol_type.I₀) * (qp_μ4[i_start] / (qp_μ4[i] - qp_μ4[i_start])) * (exp(-dτ_λ[n] / qp_μ4[i]) - exp(-dτ_λ[n] / qp_μ4[i_start]))
    end
    # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
    J₀⁻[i,n] += ϖ_λ[n] * (Z⁻⁺[i,i_start:i_end]*pol_type.I₀) * (qp_μ4[i_start] / (qp_μ4[i] + qp_μ4[i_start])) * (1 - exp.(-dτ_λ[n] * ((1 / qp_μ4[i]) + (1 / qp_μ4[i_start]))))

    J₀⁺[i,n] *= exp(-τ_sum[n]/qp_μ4[j0])
    J₀⁻[i,n] *= exp(-τ_sum[n]/qp_μ4[j0])
    
    if ndoubl >= 1
        J₀⁻[i,n] = D[i,i]*J₀⁻[i,n] #Suniti: define D = Diagonal{1,1,-1,-1,...Nquad times}
    end
        
    
end

function elemental!(pol_type, SFI, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                              ndoubl, scatter, qp_μ, wt_μ, 
                              added_layer::AddedLayer{FT}, 
                              I_static,
                              arr_type,
                              architecture) where {FT}
    
    elemental_helper!(pol_type, SFI, iμ0, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, added_layer, I_static, arr_type, architecture)
    ### synchronize()
end