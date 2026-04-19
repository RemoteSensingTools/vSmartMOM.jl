#=
 
This file contains RT elemental-related functions
 
=#
#=
"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,       #{FT2,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ_λ::AbstractArray{FT,1},   # ϖ_λ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ::FT,                      # ϖ: single scattering albedo of elemental layer (no trace gas absorption included)
                            Z⁺⁺::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺::AbstractArray{FT,2},   # Z matrix
                            τ̇_sum::AbstractArray{FT,2},       #{FT2,1}, #Suniti
                            dτ̇_λ::AbstractArray{FT,2},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ̇::AbstractArray{FT,1},                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ̇_λ::AbstractArray{FT,2},   # ϖ_λ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ̇::AbstractArray{FT,1},                      # ϖ: single scattering albedo of elemental layer (no trace gas absorption included)
                            Ż⁺⁺::AbstractArray{FT,3},   # Z matrix
                            Ż⁻⁺::AbstractArray{FT,3},   # Z matrix
                            F₀::AbstractArray{FT,2},    # Stokes vector of solar/stellar irradiance
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::AddedLayer{FT}, 
                            added_layer_lin::AddedLayerLin{FT}, 
                            I_static,
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ṙ⁺⁻, ṙ⁻⁺, ṫ⁻⁻, ṫ⁺⁺, J̇₀⁺, J̇₀⁻ = added_layer_lin
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
    #@unpack ϖ_Cabannes = RS_type
    arr_type = array_type(architecture)
    Nparams = size(Ż⁻⁺)[1]
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    # @show Array(τ_sum)[1], Array(dτ_λ)[1], Array(ϖ_λ)[1], Array(Z⁺⁺)[1,1]
    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    Ż⁺⁺_ = reshape(Ż⁺⁺, (Nparams, size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Ż⁻⁺_ = reshape(Ż⁻⁺, (Nparams, size(Z⁺⁺,1), size(Z⁺⁺,2),1))

    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))
    I₀_NquadN = arr_type(zeros(FT,size(qp_μN,1))); #incident irradiation
    i_end     = pol_type.n*iμ₀
    I₀_NquadN[iμ₀Nstart:i_end] = pol_type.I₀

    device = devi(architecture)

    # If in scattering mode:
    if scatter
   
        NquadN = length(qp_μN)

        # Needs explanation still, different weights: 
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct0  = m == 0 ? FT(0.50) * ϖ * dτ     : FT(0.25) * ϖ * dτ
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        wct   = m == 0 ? FT(0.50) * ϖ * wt_μN  : FT(0.25) * ϖ * wt_μN
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4
        wct0_lin = zeros(2)
        wct_lin  = zeros(2)
        wct0_lin[1] = m == 0 ? FT(0.50) * ϖ    : FT(0.25) * ϖ
        wct0_lin[2] = m == 0 ? FT(0.50) * dτ   : FT(0.25) * dτ
        wct_lin[1]  = 0
        wct_lin[2]  = m == 0 ? FT(0.50) * wt_μN: FT(0.25) * wt_μN

        # Get the diagonal matrices first
        d_qp  = Diagonal(1 ./ qp_μN)
        d_wct = Diagonal(wct)

        # Calculate r⁻⁺ and t⁺⁺
        
        # Version 1: no absorption in batch mode (initiation of a single scattering layer with no or low absorption)
        if false #maximum(dτ_λ) < 0.0001   
            # R⁻⁺₀₁(λ) = M⁻¹(0.5ϖₑ(λ)Z⁻⁺C)δ (See Eqs.7 in Raman paper draft)
            r⁻⁺[:,:,:]   .= d_qp * Z⁻⁺ * (d_wct * dτ)
            ṙ⁻⁺[:,:,:,1] .= d_qp * Z⁻⁺ * (d_wct + Diagonal(wct_lin[1]) * dτ)
            ṙ⁻⁺[:,:,:,2] .= d_qp * Z⁻⁺ * (Diagonal(wct_lin[2]) * dτ)
            ṙ⁻⁺[:,:,:,3] .= d_qp * d_wct * dτ
            # T⁺⁺₀₁(λ) = {I-M⁻¹[I - 0.5*ϖₑ(λ)Z⁺⁺C]}δ (See Eqs.7 in Raman paper draft)
            t⁺⁺[:,:,:] .= I_static - (d_qp * ((I_static - Z⁺⁺ * d_wct) * dτ))
            ṫ⁺⁺[:,:,:,1] .= (d_qp * (
                (Z⁺⁺ * Diagonal(wct_lin[1])) * dτ -
                (I_static - Z⁺⁺ * d_wct)))
            ṫ⁺⁺[:,:,:,2] .= (d_qp * (
                (Z⁺⁺ * Diagonal(wct_lin[2])) * dτ))
            ṫ⁺⁺[:,:,:,3] .=  d_qp * d_wct * dτ
            if SFI
                # Reminder: Add equation here what it does
                expk = exp.(-τ_sum/qp_μ[iμ₀]) #exp(-τ(z)/μ₀)
                # derivative with respect to dτ only
                expk_lin = exp.(-τ_sum/qp_μ[iμ₀]) * (-1/qp_μ[iμ₀]) 
                # J₀⁺ = 0.5[1+δ(m,0)]M⁻¹ϖₑ(λ)Z⁺⁺τI₀exp(-τ(z)/μ₀)
                J₀⁺[:,1,:]   .= (d_qp * Z⁺⁺ * I₀_NquadN * wct0) .* expk'
                J̇₀⁺[:, 1, :, 1] .= d_qp * Z⁺⁺ * I₀_NquadN * 
                                (wct0_lin[1] .* expk' + wct0 .* expk_lin')
                J̇₀⁺[:, 1, :, 2] .= (d_qp * Z⁺⁺ * I₀_NquadN * wct0_lin[2]) .* expk'
                J̇₀⁺[:, 1, :, 3] .= (d_qp * I₀_NquadN * wct0) .* expk'
                # J₀⁻ = 0.5[1+δ(m,0)]M⁻¹ϖₑ(λ)Z⁻⁺τI₀exp(-τ(z)/μ₀)
                J₀⁻[:,1,:]   .= (d_qp * Z⁻⁺ * I₀_NquadN * wct0) .* expk'
                J̇₀⁻[:, 1, :, 1] .= d_qp * Z⁻⁺ * I₀_NquadN * 
                                (wct0_lin[1] .* expk' + wct0 .* expk_lin')
                J̇₀⁻[:, 1, :, 2] .= (d_qp * Z⁻⁺ * I₀_NquadN * wct0_lin[2]) .* expk'
                J̇₀⁻[:, 1, :, 3] .= (d_qp * I₀_NquadN * wct0) .* expk'
            end
        else 
            # Version 2: More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
            # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
            kernel! = get_elem_rt!(device)
            event = kernel!(r⁻⁺, t⁺⁺, 
                ṙ⁻⁺, ṫ⁺⁺,
                ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, F₀,
                qp_μN, wct2, ndrange=size(r⁻⁺)); 
            #wait(device, event)
            synchronize_if_gpu()

            if SFI
                kernel! = get_elem_rt_SFI!(device)
                event = kernel!(J₀⁺, J₀⁻, 
                    J̇₀⁺, J̇₀⁻, 
                    ϖ_λ, dτ_λ, 
                    τ_sum, Z⁻⁺, Z⁺⁺, F₀, 
                    qp_μN, ndoubl, wct02, 
                    pol_type.n, arr_type(pol_type.I₀), iμ₀, D, ndrange=size(J₀⁺))
                #wait(device, event)
                synchronize_if_gpu()
            end
        end

        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, 
                                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                                ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)
        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, J₀⁻, J̇₀⁻)
        end      
    else 
        # Note: τ is not defined here
        # 01/16/25 Check why this is still tolerated (is this code still active?)
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
        ṫ⁺⁺[1, :] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}
        ṫ⁻⁻[1, :] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end
=#
"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                dτ::AbstractArray,
                F₀::AbstractArray,#{FT,2},    # Stokes vector of solar/stellar irradiance
                computed_layer_properties,
                m::Int,                     # m: fourier moment
                ndoubl::Int,                # ndoubl: number of doubling computations needed 
                scatter::Bool,              # scatter: flag indicating scattering
                quad_points::QuadPoints{FT}, # struct with quadrature points, weights, 
                added_layer::AddedLayer{FT}, 
                added_layer_lin::AddedLayerLin{FT}, 
                architecture) where {FT<:AbstractFloat}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ṙ⁺⁻, ṙ⁻⁺, ṫ⁻⁻, ṫ⁺⁺, J̇₀⁺, J̇₀⁻ = added_layer_lin
    @unpack qp_μ, iμ₀, wt_μN, qp_μN = quad_points
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    #@unpack ϖ_Cabannes = RS_type
    #@show architecture
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    qp_μN = arr_type(qp_μN)
    wt_μN = arr_type(wt_μN)
    #τ_sum = arr_type(τ_sum)
    #τ̇_sum = arr_type(τ̇_sum)
    I₀    = arr_type(pol_type.I₀)
    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    elem_fct = FT(1.0 / 2^ndoubl)

    device = devi(architecture)
    #@show typeof(ϖ),typeof(dτ),typeof(Z⁻⁺),typeof(Z⁺⁺) 
    #ϖ   = arr_type(ϖ);
    #dτ  = arr_type(dτ);
    #Z⁻⁺ = arr_type(Z⁻⁺);
    #Z⁺⁺ = arr_type(Z⁺⁺);
    #@show size(Z⁻⁺), size(ϖ)
    
    r⁻⁺ .= FT(0.0) 
    t⁺⁺ .= FT(0.0)
    ṙ⁻⁺ .= FT(0.0)
    ṫ⁺⁺ .= FT(0.0)
    J₀⁺ .= FT(0.0)
    J₀⁻ .= FT(0.0)
    J̇₀⁺ .= FT(0.0)
    J̇₀⁻ .= FT(0.0)

    # If in scattering mode:
    if scatter
   
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4
                        
        # More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
        # with absorption in batch mode, low tau_scatt but higher tau_total, needs exact equations
        kernel! = get_elem_rt!(device)
        #@show "Start event",   typeof(wct2)
        event = kernel!(r⁻⁺, t⁺⁺,
                    ṙ⁻⁺, ṫ⁺⁺, 
                    ϖ, dτ, Z⁻⁺, Z⁺⁺, 
                    qp_μN, wct2, elem_fct, ndrange=size(r⁻⁺)); 
        #@show t⁺⁺, any(isnan, t⁺⁺)

        #@show "Stop event"
        #wait(device, event)
        synchronize_if_gpu()

        if SFI
            kernel! = get_elem_rt_SFI!(device)
            #@show size(F₀)
            event = kernel!(J₀⁺, J₀⁻, 
                J̇₀⁺, J̇₀⁻, 
                ϖ, dτ, 
                #arr_type(τ_sum), arr_type(τ̇_sum), 
                Z⁻⁺, Z⁺⁺, 
                arr_type(F₀), 
                qp_μN, ndoubl, wct02, elem_fct,
                pol_type.n, I₀, iμ₀, D, ndrange=size(J₀⁺))
            #wait(device, event)
        end
        #ii = pol_type.n*(iμ0-1)+1
        #@show 'B',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
        synchronize_if_gpu()
        
        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n,         
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, J₀⁻, J̇₀⁻)
        end      
    else
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-dτ' ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-dτ' ./ qp_μN)}
        ṫ⁺⁺[:, 1] = Diagonal{exp(-dτ' ./ qp_μN).*reshape(-1 ./ qp_μN, length(qp_μN), 1)}*elem_fct
        ṫ⁻⁻[:, 1] = Diagonal{exp(-dτ' ./ qp_μN).*reshape(-1 ./ qp_μN, length(qp_μN), 1)}*elem_fct
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

@kernel function get_elem_rt!(r⁻⁺, t⁺⁺,
                        ṙ⁻⁺, ṫ⁺⁺, 
                        ϖ_λ, dτ_λ, 
                        Z⁻⁺, Z⁺⁺, 
                        qp_μN, wct, elem_fct) 
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    if (wct[j]>1.e-8) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dτ = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(1/μᵢ) ̇exp{-τ ̇(1/μᵢ + 1/μⱼ)}  ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dϖ = 𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dZ = ϖ ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        tmpF = (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))) 
        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * tmpF
            #Z⁻⁺[i,j] * 
            
        # derivative wrt τ_λ
        ṙ⁻⁺[i,j,n,1] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * elem_fct *
            (1/qp_μN[i]) * wct[j] * 
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))) 
        # derivative wrt ϖ
        ṙ⁻⁺[i,j,n,2] = Z⁻⁺[i,j,n2] * tmpF
            #Z⁻⁺[i,j] * 
            
        # derivative wrt Z
        ṙ⁻⁺[i,j,n,3] = ϖ_λ[n] * tmpF
                    
        #if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ}(1 + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ))) ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dτ_λ = (exp{-τ/μⱼ}/μᵢ)⋅(ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)⋅(1-τ/μⱼ)-1) ̇𝑤ⱼ  
            # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
        if qp_μN[i] == qp_μN[j]
            if i==j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = 
                    exp(-dτ_λ[n] / qp_μN[i]) * (-1 / qp_μN[i]) * 
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (dτ_λ[n] / qp_μN[i] - 1)) *
                    elem_fct
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = exp(-dτ_λ[n] / qp_μN[i]) *
                    (Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                    
                # derivative wrt Z
                ṫ⁺⁺[i,j,n,3] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[i])
            else
                tmpF = exp(-dτ_λ[n] / qp_μN[i]) *
                    (dτ_λ[n] / qp_μN[i]) * wct[j]

                t⁺⁺[i,j,n] = 
                    ϖ_λ[n] * Z⁺⁺[i,j,n2] * tmpF
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = 
                    exp(-dτ_λ[n] / qp_μN[i]) * (-1 / qp_μN[i]) * 
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * wct[j] * (dτ_λ[n] / qp_μN[i] - 1)) *
                    elem_fct
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = 
                    Z⁺⁺[i,j,n2] * tmpF   
                # derivative wrt Z
                ṫ⁺⁺[i,j,n,3] = 
                    ϖ_λ[n] * tmpF
            end
            #    # 𝐓⁺⁺(μᵢ, μⱼ) = (exp{-τ/μⱼ}(ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(τ/μᵢ))) ̇𝑤ⱼ        
            #    # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = (exp{-τ/μⱼ}⋅ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)/μᵢ)⋅(1 - τ/μⱼ) ̇𝑤ⱼ
            #    # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            #    # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            #    t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
            #        (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
            #    # derivative wrt τ_λ
            #    ṫ⁺⁺[i,j,n,1] = (exp(-dτ_λ[n] / qp_μN[j]) *
            #            ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) * 
            #            (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
            #    # derivative wrt ϖ_λ
            #    ṫ⁺⁺[i,j,n,2] = t⁺⁺[i,j,n] / ϖ_λ[n]
            #    # derivative wrt Z
            #    ṫ⁺⁺[i,j,n,3] = t⁺⁺[i,j,n] / Z⁺⁺[i,j,n2]
            #end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = -ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ}/μᵢ - exp{-τ/μⱼ}/μⱼ) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            tmpF = (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j]))
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * Z⁺⁺[i,j,n2] * tmpF
            #if n==1
            #    @show i, j, n, t⁺⁺[i,j,n], qp_μN[i], qp_μN[j]
            #end
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = -ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * elem_fct *
                (exp(-dτ_λ[n] / qp_μN[i])/ qp_μN[i] - 
                exp(-dτ_λ[n] / qp_μN[j])/ qp_μN[j]) 
            # derivative wrt ϖ_λ
            ṫ⁺⁺[i,j,n,2] = Z⁺⁺[i,j,n2] * tmpF
            # derivative wrt Z
            ṫ⁺⁺[i,j,n,3] = ϖ_λ[n] * tmpF
        end
    else
    
        r⁻⁺[i,j,n] = 0.0
        # derivative wrt τ_λ
        ṙ⁻⁺[i,j,n,1] = 0.0
        # derivative wrt ϖ
        ṙ⁻⁺[i,j,n,2] = 0.0
        # derivative wrt Z
        ṙ⁻⁺[i,j,n,3] = 0.0
                    
        #if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ}(1 + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ))) ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dτ_λ = (exp{-τ/μⱼ}/μᵢ)⋅(ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)⋅(1-τ/μⱼ)-1) ̇𝑤ⱼ  
            # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
        if i==j
            t⁺⁺[i,j,n] = 
                exp(-dτ_λ[n] / qp_μN[i]) 
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = 
                 t⁺⁺[i,j,n] * (-1 / qp_μN[i]) * 
                elem_fct
            # derivative wrt ϖ_λ
            ṫ⁺⁺[i,j,n,2] = 
                0.0  
            # derivative wrt Z
            ṫ⁺⁺[i,j,n,3] = 
                0.0
            #else
            #    # 𝐓⁺⁺(μᵢ, μⱼ) = (exp{-τ/μⱼ}(ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(τ/μᵢ))) ̇𝑤ⱼ        
            #    # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = (exp{-τ/μⱼ}⋅ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)/μᵢ)⋅(1 - τ/μⱼ) ̇𝑤ⱼ
            #    # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            #    # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            #    t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
            #        (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
            #    # derivative wrt τ_λ
            #    ṫ⁺⁺[i,j,n,1] = (exp(-dτ_λ[n] / qp_μN[j]) *
            #            ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) * 
            #            (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
            #    # derivative wrt ϖ_λ
            #    ṫ⁺⁺[i,j,n,2] = t⁺⁺[i,j,n] / ϖ_λ[n]
            #    # derivative wrt Z
            #    ṫ⁺⁺[i,j,n,3] = t⁺⁺[i,j,n] / Z⁺⁺[i,j,n2]
            #end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = -ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ}/μᵢ - exp{-τ/μⱼ}/μⱼ) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 0.0
                
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = 0.0
            # derivative wrt ϖ_λ
            ṫ⁺⁺[i,j,n,2] = 0.0
            # derivative wrt Z
            ṫ⁺⁺[i,j,n,3] = 0.0
        end
    end
    nothing
end

@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, 
                J̇₀⁺, J̇₀⁻, 
                ϖ_λ, dτ_λ, 
                #τ_sum, τ̇_sum, 
                Z⁻⁺, Z⁺⁺, F₀,
                qp_μN, ndoubl, wct02, elem_fct,
                nStokes,
                I₀, iμ0, D)
    
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    J₀⁺[i, 1, n]=0
    J₀⁻[i, 1, n]=0
    J̇₀⁺[i, 1, n, 1:3].=0
    J̇₀⁻[i, 1, n, 1:3].=0
    n2=1
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    #if scatter 
        Z⁺⁺_I₀ = FT(0.0);
        Z⁻⁺_I₀ = FT(0.0);
        
        for ii = i_start:i_end
            Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n] #I₀[ii-i_start+1]
            Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n] #I₀[ii-i_start+1] 
        end

        if (i>=i_start) && (i<=i_end)
            #ctr = i-i_start+1
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
            tmpF = wct02 * (dτ_λ[n] / qp_μN[i]) * exp(-dτ_λ[n] / qp_μN[i])
            J₀⁺[i, 1, n] = ϖ_λ[n] * Z⁺⁺_I₀ * tmpF
            # derivative wrt τ
            J̇₀⁺[i, 1, n, 1] = exp(-dτ_λ[n] / qp_μN[i]) * (1 / qp_μN[i]) * 
                    ϖ_λ[n] * Z⁺⁺_I₀ * (1 - dτ_λ[n] / qp_μN[i]) *
                    wct02 * elem_fct
            # derivative wrt ϖ
            J̇₀⁺[i, 1, n, 2] = Z⁺⁺_I₀ * tmpF
            # derivative wrt Z
            J̇₀⁺[i, 1, n, 3] = ϖ_λ[n] * F₀[1,n] * tmpF #Suniti: if the incident starlight were polarized, the third index of J₀ would be 3 for I, Q, and U (or 4 including V) components of F₀ instead of 1
        else
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
            tmpF = wct02 *  
                (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * 
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
            J₀⁺[i, 1, n] = ϖ_λ[n] * Z⁺⁺_I₀ * tmpF 

            # derivative wrt τ
            J̇₀⁺[i, 1, n, 1] = - wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * elem_fct *
                (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * 
                (exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i] - exp(-dτ_λ[n] / qp_μN[i_start]) / qp_μN[i_start])
            # derivative wrt ϖ
            J̇₀⁺[i, 1, n, 2] = Z⁺⁺_I₀ * tmpF
            # derivative wrt Z
            J̇₀⁺[i, 1, n, 3] = ϖ_λ[n] * F₀[1,n] * tmpF #Suniti: if the incident starlight were polarized, the third index of J₀ would be 3 for I, Q, and U (or 4 including V) components of F₀ instead of 1
        end
        #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
        tmpF = wct02 * 
            (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))
        J₀⁻[i, 1, n] = ϖ_λ[n] * Z⁻⁺_I₀ * tmpF
        # derivative wrt τ
        J̇₀⁻[i, 1, n, 1] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * elem_fct *
                (1 / qp_μN[i]) * 
                exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start])))
        # derivative wrt ϖ
        J̇₀⁻[i, 1, n, 2] = Z⁻⁺_I₀ * tmpF
        # derivative wrt Z
        J̇₀⁻[i, 1, n, 3] = ϖ_λ[n] * F₀[1,n] * tmpF #Suniti: if the incident starlight were polarized, the third index of J₀ would be 3 for I, Q, and U (or 4 including V) components of F₀ instead of 1
    #else
    #       J̇₀⁺[i, 1, n, 1:3] .= 0.0
    #       J̇₀⁻[i, 1, n, 1:3] .= 0.0
    #       J₀⁺[i, 1, n] = 0.0
    #       J₀⁻[i, 1, n] = 0.0
    #end
    ## Suniti: start here after lunch on 02/25/26    
    ## TODO: Move this out until after doubling (it is not necessary to consider this here already if Raman scattering is not involved)
    #J₀⁺[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])
    #J₀⁻[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])

    #J̇₀⁺[i, 1, n, 1] = J̇₀⁺[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start]) +
    #                    J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    #J̇₀⁻[i, 1, n, 1] = J̇₀⁻[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start]) +
    #                    J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    #J̇₀⁺[i, 1, n, 2] = J̇₀⁺[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
    #                    #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    #J̇₀⁻[i, 1, n, 2] = J̇₀⁻[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
    #                    #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    #J̇₀⁺[i, 1, n, 3] = J̇₀⁺[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
    #                    #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    #J̇₀⁻[i, 1, n, 3] = J̇₀⁻[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
    #                    #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #D = Diagonal{1,1,-1,-1,...Nquad times}
        J̇₀⁻[i, 1, n, 1] = D[i,i]*J̇₀⁻[i, 1, n, 1]
        J̇₀⁻[i, 1, n, 2] = D[i,i]*J̇₀⁻[i, 1, n, 2]
        J̇₀⁻[i, 1, n, 3] = D[i,i]*J̇₀⁻[i, 1, n, 3]
    end  
    #if (n==840||n==850)    
    #    @show i, n, J₀⁺[i, 1, n], J₀⁻[i, 1, n]      
    #end
    nothing
end

@kernel function apply_D_elemental!(ndoubl, pol_n, 
                                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                                ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)
    i, j, n = @index(Global, NTuple) #how best to do this for linearization? Is : okay, or should I use an iparam index?

    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        #if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
        if (((1<=ii<=2) & (1<=jj<=2)) | (!(1<=ii<=2) & !(1<=jj<=2))) 
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
            ṙ⁺⁻[i,j,n,1] = ṙ⁻⁺[i,j,n,1]
            ṙ⁺⁻[i,j,n,2] = ṙ⁻⁺[i,j,n,2]
            ṙ⁺⁻[i,j,n,3] = ṙ⁻⁺[i,j,n,3]
            ṫ⁻⁻[i,j,n,1] = ṫ⁺⁺[i,j,n,1]
            ṫ⁻⁻[i,j,n,2] = ṫ⁺⁺[i,j,n,2]
            ṫ⁻⁻[i,j,n,3] = ṫ⁺⁺[i,j,n,3]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] 
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] 
            ṙ⁺⁻[i,j,n,1] = -ṙ⁻⁺[i,j,n,1] 
            ṙ⁺⁻[i,j,n,2] = -ṙ⁻⁺[i,j,n,2] 
            ṙ⁺⁻[i,j,n,3] = -ṙ⁻⁺[i,j,n,3] 
            ṫ⁻⁻[i,j,n,1] = -ṫ⁺⁺[i,j,n,1] 
            ṫ⁻⁻[i,j,n,2] = -ṫ⁺⁺[i,j,n,2] 
            ṫ⁻⁻[i,j,n,3] = -ṫ⁺⁺[i,j,n,3] 
        end
    else
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
            ṙ⁻⁺[i,j,n,1] = - ṙ⁻⁺[i,j,n,1]
            ṙ⁻⁺[i,j,n,2] = - ṙ⁻⁺[i,j,n,2]
            ṙ⁻⁺[i,j,n,3] = - ṙ⁻⁺[i,j,n,3]
        end 
    end
    nothing
end

@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻, J̇₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            J₀⁻[i, 1, n] = - J₀⁻[i, 1, n]
            J̇₀⁻[i, 1, n, 1] = - J̇₀⁻[i, 1, n, 1]
            J̇₀⁻[i, 1, n, 2] = - J̇₀⁻[i, 1, n, 2]
            J̇₀⁻[i, 1, n, 3] = - J̇₀⁻[i, 1, n, 3]
        end 
    end
    nothing
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, 
                    r⁻⁺::AbstractArray{FT,3}, 
                    t⁺⁺::AbstractArray{FT,3}, 
                    r⁺⁻::AbstractArray{FT,3}, 
                    t⁻⁻::AbstractArray{FT,3},
                    ṙ⁻⁺::AbstractArray{FT,4}, 
                    ṫ⁺⁺::AbstractArray{FT,4}, 
                    ṙ⁺⁻::AbstractArray{FT,4}, 
                    ṫ⁻⁻::AbstractArray{FT,4}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        ṙ⁺⁻[:] = ṙ⁻⁺
        ṫ⁻⁻[:] = ṫ⁺⁺
        return nothing
    end
    device = devi(architecture(r⁻⁺))
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, 
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, 
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻, 
                        ndrange=size(r⁻⁺));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, 
                                J₀⁻::AbstractArray{FT,3},
                                J̇₀⁻::AbstractArray{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(architecture(J₀⁻))
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, J̇₀⁻, ndrange=size(J₀⁻));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end