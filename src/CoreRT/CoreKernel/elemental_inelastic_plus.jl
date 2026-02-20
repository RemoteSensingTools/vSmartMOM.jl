function getKernelDim(RS_type::RRS_plus,ier⁻⁺)
    return size(ier⁻⁺);
end

function getKernelDim(RS_type::Union{VS_0to1_plus, VS_1to0_plus},ier⁻⁺, i_λ₁λ₀)
    return (size(ier⁻⁺,1),size(ier⁻⁺,2), size(i_λ₁λ₀,1));
end

function getKernelDimSFI(RS_type::RRS_plus,ieJ₀⁻)
    return size(ieJ₀⁻);
end

function getKernelDimSFI(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                        ieJ₀⁻,i_λ₁λ₀)
    return (size(ieJ₀⁻,1), size(i_λ₁λ₀,1));
end

"Elemental single-scattering layer for RRS"
function elemental_inelastic!(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                            pol_type, SFI::Bool, 
                            τ_sum::AbstractArray{FT,1},
                            dτ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            ϖ::AbstractArray{FT,1},                     # dτ:   scattering optical depth of elemental layer (scalar)
                            Z⁺⁺_λ₁λ₀::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺_λ₁λ₀::AbstractArray{FT,2}, 
                            F₀::AbstractArray{FT,2},
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            I_static,
                            architecture) where {FT<:Real,FT2}

    (; ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻) = added_layer
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀) = quad_points
    arr_type = array_type(architecture)
    τ_sum = arr_type(τ_sum)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    
    # Later on, we can have Zs also vary with index, pretty easy here:
    #Z⁺⁺_ = reshape(Z⁺⁺_λ₁λ₀, (size(Z⁺⁺_λ₁λ₀,1), size(Z⁺⁺_λ₁λ₀,2),1))
    #Z⁻⁺_ = reshape(Z⁻⁺_λ₁λ₀, (size(Z⁺⁺_λ₁λ₀,1), size(Z⁺⁺_λ₁λ₀,2),1))

    D         = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))
    
    # If in scattering mode:
    if scatter
        # Needs explanation still, different weights: 
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)

        # Calculate r⁻⁺ and t⁺⁺
        #Version 2: More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
        # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
        ier⁻⁺.=0.0 
        iet⁺⁺.=0.0
        get_elem_rt!(RS_type, ier⁻⁺, iet⁺⁺, 
            dτ, ϖ, 
            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
            qp_μN, wct2)
        
        if SFI
            ieJ₀⁺.=0.0
            ieJ₀⁻.=0.0
            get_elem_rt_SFI!(RS_type, ieJ₀⁺, ieJ₀⁻, 
                τ_sum, dτ, ϖ, 
                Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                qp_μN, ndoubl,wct02, pol_type.n, 
                arr_type(pol_type.I₀), iμ₀, D);
        end
        
        # Apply D Matrix
        apply_D_matrix_elemental!(RS_type, ndoubl, pol_type.n, 
                                    ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
        if SFI
            apply_D_matrix_elemental_SFI!(RS_type, ndoubl, pol_type.n, 
                                            ieJ₀⁻)
        end
    else 
        # Note: τ is not defined here
        iet⁺⁺[:] = 0. #Diagonal{exp(-τ ./ qp_μN)}
        iet⁻⁻[:] = 0. #Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end


#Suniti: is there a way to pass information like ϖ_λ₁λ₀, i_λ₁λ₀, i_ref, etc. along with RS_type? So that they can be retrieved as RSS.ϖ_λ₁λ₀ for example?
# This one is only for RRS
# kernel wrapper:
function get_elem_rt!(RS_type::RRS_plus, 
                        ier⁻⁺, iet⁺⁺, 
                        dτ, ϖ, 
                        #Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                        qp_μN, wct2)
    (; fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref, iBand, bandSpecLim) = RS_type
    device = devi(architecture(ier⁻⁺))
    aType = array_type(architecture(ier⁻⁺))
    kernel! = get_elem_rt_RRS!(device)
    for iB in RS_type.iBand
        event = kernel!(fscattRayl[iB], 
                    aType(ϖ_λ₁λ₀[iB]), aType(i_λ₁λ₀[iB]), 
                    i_ref,
                    ier⁻⁺[:,:,bandSpecLim[iB],:], 
                    iet⁺⁺[:,:,bandSpecLim[iB],:], 
                    dτ[bandSpecLim[iB]], 
                    #ϖ[bandSpecLim[iB]], 
                    aType(Z⁻⁺_λ₁λ₀[:,:,bandSpecLim[iB]]), 
                    aType(Z⁺⁺_λ₁λ₀[:,:,bandSpecLim[iB]]), 
                    qp_μN, wct2, 
                    ndrange=getKernelDim(RS_type,ier⁻⁺[:,:,RS_type.bandSpecLim[iB],:])); 
        #wait(device, event);
        synchronize_if_gpu();
    end
end

function get_elem_rt!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
    ier⁻⁺, iet⁺⁺, 
    dτ, ϖ,
    Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
    qp_μN, wct2)
    (; fscattRayl, i_ref,
            ϖ_λ₁λ₀, i_λ₁λ₀, 
            ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2, 
            ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2, 
            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
            Z⁻⁺_λ₁λ₀_VS_n2, Z⁺⁺_λ₁λ₀_VS_n2,
            Z⁻⁺_λ₁λ₀_VS_o2, Z⁺⁺_λ₁λ₀_VS_o2) = RS_type
    device = devi(architecture(ier⁻⁺)) #change this 
    aType = array_type(architecture(ier⁻⁺)) 
    #RVRS
    kernel! = get_elem_rt_VS!(device)
    event = kernel!(aType(fscattRayl), 
        aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
        ier⁻⁺, iet⁺⁺, 
        dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), 
        qp_μN, wct2, 
        ndrange=getKernelDim(RS_type,ier⁻⁺,i_λ₁λ₀)); 
    synchronize_if_gpu();
    t_ier⁻⁺  = similar(ier⁻⁺)
    t_iet⁺⁺  = similar(ier⁻⁺)
    t_ier⁻⁺  .= 0.0 
    t_iet⁺⁺  .= 0.0 
    #VS - N2
    kernel! = get_elem_rt_VS!(device)
    event = kernel!(fscattRayl, 
        aType(ϖ_λ₁λ₀_VS_n2), aType(i_λ₁λ₀_VS_n2),
        t_ier⁻⁺, t_iet⁺⁺, 
        dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀_VS_n2), aType(Z⁺⁺_λ₁λ₀_VS_n2), 
        qp_μN, wct2, 
        ndrange=getKernelDim(RS_type,ier⁻⁺,i_λ₁λ₀_VS_n2)); 
    synchronize_if_gpu();
    ier⁻⁺ .+= t_ier⁻⁺
    iet⁺⁺ .+= t_iet⁺⁺
    t_ier⁻⁺  .= 0.0 
    t_iet⁺⁺  .= 0.0 
    #VS - O2
    kernel! = get_elem_rt_VS!(device)
    event = kernel!(aType(fscattRayl), 
        aType(ϖ_λ₁λ₀_VS_o2), aType(i_λ₁λ₀_VS_o2),
        t_ier⁻⁺, t_iet⁺⁺, 
        dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀_VS_o2), aType(Z⁺⁺_λ₁λ₀_VS_o2),
        qp_μN, wct2, 
        ndrange=getKernelDim(RS_type,ier⁻⁺, i_λ₁λ₀_VS_o2)); 
    synchronize_if_gpu();
    ier⁻⁺ .+= t_ier⁻⁺
    iet⁺⁺ .+= t_iet⁺⁺ 
    #t_iet⁺⁺  .= 0.0 
end


@kernel function get_elem_rt_VS!(fscattRayl, 
                            ϖ_λ₁λ₀, i_λ₁λ₀, 
                            #ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2, 
                            #ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2, 
                            ier⁻⁺, iet⁺⁺, 
                            dτ, ϖ,
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            #Z⁻⁺_λ₁λ₀_VS_n2, Z⁺⁺_λ₁λ₀_VS_n2, 
                            #Z⁻⁺_λ₁λ₀_VS_o2, Z⁺⁺_λ₁λ₀_VS_o2, 
                            qp_μN, wct2)
    i, j, Δn = @index(Global, NTuple) 
    #@unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, dτ₀, dτ₀_λ = RS_type 
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    #dτ₁ = 1 #dummy for now
    #Suniti: require that the incident wavelength is always the first element of 1:nSpec, and all the others belong to the same target VS band
    #Suniti: Then,
    n₀ = 1    
    n₁ = i_λ₁λ₀[Δn]  
    if n₁ >0
        if (wct2[j]>1.e-8)
            
            # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
            # optical thicknesses at wavelengths λ₀ and λ₁
            # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
            ier⁻⁺[i,j,n₁,1] += 
                    ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁻⁺_λ₁λ₀[i,j] * 
                    (1/( (qp_μN[i] / qp_μN[j]) + (dτ[n₁]/dτ[n₀]) )) * 
                    (1 - exp(-((dτ[n₁] / qp_μN[i]) + (dτ[n₀] / qp_μN[j])))) * wct2[j] 
            
            if (qp_μN[i] == qp_μN[j])
                # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
                #if i == j       
                if abs(dτ[n₀]-dτ[n₁])>1.e-8
                    iet⁺⁺[i,j,n₁,1] += 
                        ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * wct2[j] *
                        (exp(-dτ[n₁] / qp_μN[i]) - exp(-dτ[n₀] / qp_μN[j]))/
                        (1 - (dτ[n₁]/dτ[n₀]))  
                else    
                    iet⁺⁺[i,j,n₁,1] += 
                        (dτ[n₀] / qp_μN[i]) * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * 
                        Z⁺⁺_λ₁λ₀[i,j] * wct2[j] *
                        exp(-dτ[n₀] / qp_μN[j])   
                    end
                #else
                #    iet⁺⁺[i,j,n₁,1] += 0.0
                #end
            else
                # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
                # (𝑖 ≠ 𝑗)
                if (abs( (qp_μN[i]/qp_μN[j]) - (dτ[n₁]/dτ[n₀]) ) < 1.e-8)
                    iet⁺⁺[i,j,n₁,1] = 
                        (dτ[n₀]/qp_μN[i]) * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * 
                        Z⁺⁺_λ₁λ₀[i,j] * 
                        wct2[j] * exp(-dτ[n₀] / qp_μN[j])
                else
                    iet⁺⁺[i,j,n₁,1] += 
                            ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * 
                            (1 / ( (qp_μN[i]/qp_μN[j]) - (dτ[n₁]/dτ[n₀]) )) * 
                            wct2[j] * 
                            (exp(-dτ[n₁] / qp_μN[i]) - exp(-dτ[n₀] / qp_μN[j]))
                end
            end
        else
            ier⁻⁺[i,j,n₁,1] += 0.0
            if i==j
                iet⁺⁺[i,j,n₁,1] += 0.0
            else
                iet⁺⁺[i,j,n₁,1] += 0.0
            end
        end
    end
end

function get_elem_rt_SFI!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                        ieJ₀⁺, ieJ₀⁻, 
                        τ_sum, 
                        dτ, ϖ, 
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                        qp_μN, ndoubl,
                        wct02, nStokes,
                        I₀, iμ0,D)
    (; fscattRayl, i_ref,
    ϖ_λ₁λ₀, i_λ₁λ₀, 
    ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2, 
    ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2, 
    Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
    Z⁻⁺_λ₁λ₀_VS_n2, Z⁺⁺_λ₁λ₀_VS_n2,
    Z⁻⁺_λ₁λ₀_VS_o2, Z⁺⁺_λ₁λ₀_VS_o2) = RS_type

    device = devi(architecture(ieJ₀⁺))
    aType = array_type(architecture(ieJ₀⁺))
    kernel! = get_elem_rt_SFI_VS!(device)
    event = kernel!(aType(fscattRayl), 
        aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
        ieJ₀⁺, ieJ₀⁻, 
        τ_sum, dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), aType(F₀),
        qp_μN, ndoubl, wct02, nStokes, 
        I₀, iμ0, D, 
        ndrange=getKernelDimSFI(RS_type, ieJ₀⁻, i_λ₁λ₀)); #change this
    #wait(device, event)
    synchronize_if_gpu();

    t_ieJ₀⁺ = similar(ieJ₀⁻)
    t_ieJ₀⁻ = similar(ieJ₀⁻)
    t_ieJ₀⁻.= 0.0 
    t_ieJ₀⁺.= 0.0 
    
    event = kernel!(aType(fscattRayl), 
        aType(ϖ_λ₁λ₀_VS_n2), aType(i_λ₁λ₀_VS_n2), 
        t_ieJ₀⁺, t_ieJ₀⁻, 
        τ_sum, dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀_VS_n2), aType(Z⁺⁺_λ₁λ₀_VS_n2), aType(F₀),
        qp_μN, ndoubl, wct02, nStokes, 
        I₀, iμ0, D, 
        ndrange=getKernelDimSFI(RS_type, ieJ₀⁻, i_λ₁λ₀_VS_n2)); #change this
    #wait(device, event)
    synchronize_if_gpu();
    
    ieJ₀⁺ .+= t_ieJ₀⁺
    ieJ₀⁻ .+= t_ieJ₀⁻
    t_ieJ₀⁻ .= 0.0 
    t_ieJ₀⁺ .= 0.0 
    event = kernel!(fscattRayl, 
        aType(ϖ_λ₁λ₀_VS_o2), aType(i_λ₁λ₀_VS_o2), 
        t_ieJ₀⁺, t_ieJ₀⁻, 
        τ_sum, dτ, ϖ, 
        aType(Z⁻⁺_λ₁λ₀_VS_o2), aType(Z⁺⁺_λ₁λ₀_VS_o2), aType(F₀),
        qp_μN, ndoubl, wct02, nStokes, 
        I₀, iμ0, D, 
        ndrange=getKernelDimSFI(RS_type, ieJ₀⁻, i_λ₁λ₀_VS_o2)); #change this
    #wait(device, event)
    synchronize_if_gpu();
        
    ieJ₀⁺ .+= t_ieJ₀⁺
    ieJ₀⁻ .+= t_ieJ₀⁻ 
end

#  TODO: Nov 30, 2021
@kernel function get_elem_rt_SFI_VS!(fscattRayl,
                            ϖ_λ₁λ₀, i_λ₁λ₀, 
                            ieJ₀⁺, ieJ₀⁻, 
                            τ_sum, dτ, ϖ, 
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                            qp_μN, ndoubl,
                            wct02, nStokes, 
                            I₀, iμ0, D)
    
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0

    i, Δn = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 

    #Suniti: require that the incident wavelength is always the first element of 1:nSpec, and all the others belong to the same target VS band
    #Suniti: Then,
    n₀ = 1    
    n₁ = i_λ₁λ₀[Δn]  
      
    FT = eltype(I₀)
    
    if n₁>0
        ieJ₀⁺[i, 1, n₁, 1]=0
        ieJ₀⁻[i, 1, n₁, 1]=0
        Z⁺⁺_I₀ = FT(0.0);
        Z⁻⁺_I₀ = FT(0.0);
        for ii = i_start:i_end
            Z⁺⁺_I₀ += Z⁺⁺_λ₁λ₀[i,ii] * F₀[ii-i_start+1, n₀] #I₀[ii-i_start+1]
            Z⁻⁺_I₀ += Z⁻⁺_λ₁λ₀[i,ii] * F₀[ii-i_start+1, n₀] #I₀[ii-i_start+1] 
        end
        
        if (i_start ≤ i ≤ i_end)
            #ctr = i-i_start+1
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
            if abs(dτ[n₀]-dτ[n₁])>1.e-8
                ieJ₀⁺[i, 1, n₁, 1] = 
                        ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_I₀ * wct02 *
                        (exp(-dτ[n₁] / qp_μN[i]) - exp(-dτ[n₀] / qp_μN[i])) /
                        (1-(dτ[n₁]/dτ[n₀])) 
                        
            else
                ieJ₀⁺[i, 1, n₁, 1] = 
                        (dτ[n₀]/qp_μN[i]) *  wct02 * ϖ_λ₁λ₀[Δn] * 
                        fscattRayl[n₀] * Z⁺⁺_I₀ * 
                        exp(-dτ[n₀] / qp_μN[i])
            end
        else
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
            if (abs( (qp_μN[i]/qp_μN[i_start]) - (dτ[n₁]/dτ[n₀]) ) < 1.e-8)
                ieJ₀⁺[i, 1, n₁, 1] = 
                    (dτ[n₀]/qp_μN[i]) * wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * 
                    Z⁺⁺_I₀ * 
                    exp(-dτ[n₀] / qp_μN[i_start])
            else
            ieJ₀⁺[i, 1, n₁, 1] = 
                        wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_I₀ * 
                        (1 /( (qp_μN[i]/qp_μN[i_start]) - (dτ[n₁]/dτ[n₀]) ) ) * 
                        (exp(-dτ[n₁] / qp_μN[i]) - exp(-dτ[n₀] / qp_μN[i_start]))  
            end
        end
        #TODO
        #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]                    
        ieJ₀⁻[i, 1, n₁, 1] = 
                    wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁻⁺_I₀ * 
                    (1/( (qp_μN[i] / qp_μN[i_start]) + (dτ[n₁]/dτ[n₀]) )) *
                    (1 - exp(-( (dτ[n₁] / qp_μN[i]) + (dτ[n₀] / qp_μN[i_start]) ) ))  

        ieJ₀⁺[i, 1, n₁, 1] *= exp(-τ_sum[n₀]/qp_μN[i_start])
        ieJ₀⁻[i, 1, n₁, 1] *= exp(-τ_sum[n₀]/qp_μN[i_start])

        if ndoubl >= 1
            ieJ₀⁻[i, 1, n₁, 1] = D[i,i]*ieJ₀⁻[i, 1, n₁, 1] #D = Diagonal{1,1,-1,-1,...Nquad times}
        end     
    end    
end

#  TODO: Nov 30, 2021
function get_elem_rt_SFI!(RS_type::RRS_plus, 
                        ieJ₀⁺, ieJ₀⁻, 
                        τ_sum, dτ_λ, ϖ_λ, 
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                        qp_μN, ndoubl,
                        wct02, nStokes,
                        I₀, iμ0,D)
    (; fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref) = RS_type
    device = devi(architecture(ieJ₀⁺))
    aType = array_type(architecture(ieJ₀⁺))
    kernel! = get_elem_rt_SFI_RRS!(device)
    event = kernel!(aType(fscattRayl), aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
                i_ref, ieJ₀⁺, ieJ₀⁻, 
                τ_sum, dτ_λ, ϖ_λ, 
                aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), aType(F₀),
                qp_μN, ndoubl, wct02, nStokes, 
                I₀, iμ0, D, 
                ndrange=getKernelDimSFI(RS_type,ieJ₀⁻));
    #wait(device, event)
    synchronize_if_gpu();
end

# only for RRS