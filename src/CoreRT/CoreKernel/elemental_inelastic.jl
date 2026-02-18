function getKernelDim(RS_type::RRS,ier⁻⁺)
    return size(ier⁻⁺);
end

function getKernelDim(RS_type::Union{VS_0to1, VS_1to0},ier⁻⁺)
    return (size(ier⁻⁺,1),size(ier⁻⁺,2), size(RS_type.i_λ₁λ₀));
end

function getKernelDimSFI(RS_type::RRS,ieJ₀⁻)
    return size(ieJ₀⁻);
end

function getKernelDimSFI(RS_type::Union{VS_0to1, VS_1to0},ieJ₀⁻)
    return (size(ieJ₀⁻,1),size(ieJ₀⁻,2), size(RS_type.i_λ₁λ₀));
end

"Elemental single-scattering layer for RRS"
function elemental_inelastic!(RS_type::Union{RRS, RRS_plus},
                            pol_type, SFI::Bool, 
                            τ_sum::AbstractArray{FT,1},
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            ϖ_λ::AbstractArray{FT,1},                     # dτ:   scattering optical depth of elemental layer (scalar)
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

    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
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
        #@show RS_type
        get_elem_rt!(RS_type, ier⁻⁺, iet⁺⁺, 
            dτ_λ, ϖ_λ, Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, qp_μN, wct2)
        
        if SFI
            get_elem_rt_SFI!(RS_type, ieJ₀⁺, ieJ₀⁻, 
                τ_sum, dτ_λ, ϖ_λ, Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀,
                F₀, 
                qp_μN, ndoubl,wct02, pol_type.n, 
                arr_type(pol_type.I₀), iμ₀, D);
        end
        # Apply D Matrix
        apply_D_matrix_elemental!(RS_type, ndoubl, pol_type.n, 
                                    ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
        #println("Apply D matrix done")
        if SFI
            #@show "here 1"
            apply_D_matrix_elemental_SFI!(RS_type, 
                                        ndoubl, 
                                        pol_type.n, 
                                        ieJ₀⁻)
            #@show "here 2"
        end
        #println("Apply D matrix SFI done")      
    else 
        # Note: τ is not defined here
        iet⁺⁺[:] = 0.0 #Diagonal{exp(-τ ./ qp_μN)}
        iet⁻⁻[:] = 0.0 #Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

#Suniti: is there a way to pass information like ϖ_λ₁λ₀, i_λ₁λ₀, i_ref, etc. along with RS_type? So that they can be retrieved as RSS.ϖ_λ₁λ₀ for example?
# This one is only for RRS
@kernel function get_elem_rt_RRS!(fscattRayl, 
                            ϖ_λ₁λ₀, i_λ₁λ₀, i_ref,
                            ier⁻⁺, iet⁺⁺, 
                            dτ_λ, ϖ_λ,
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, wct2)

    i, j, n₁, Δn = @index(Global, NTuple)
    
    nMax = length(dτ_λ) 
    # n₁ covers the full range of wavelengths, while n₀ = n₁+Δn only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    n₀  = n₁ + i_λ₁λ₀[Δn]
    ier⁻⁺[i,j,n₁,Δn]=0
    iet⁺⁺[i,j,n₁,Δn]=0

    
    if (1 ≤ n₀ ≤ nMax) & (wct2[j]>1.e-8) 

        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,Δn] = 
            fscattRayl[n₀] * ϖ_λ₁λ₀[Δn] * Z⁻⁺_λ₁λ₀[i,j] * 
            (1/( (qp_μN[i] / qp_μN[j]) + (dτ_λ[n₁]/dτ_λ[n₀]) )) * 
            (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * wct2[j] 


        if (qp_μN[i] == qp_μN[j])
            # @show i,j
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            #if i == j       
                if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-8
                    iet⁺⁺[i,j,n₁,Δn] = 
                        ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * wct2[j] *
                        (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j]))/
                        (1 - (dτ_λ[n₁]/dτ_λ[n₀]))   
                    
                else    
                    iet⁺⁺[i,j,n₁,Δn] =  
                        (dτ_λ[n₀]/ qp_μN[i]) * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * 
                        Z⁺⁺_λ₁λ₀[i,j] * wct2[j] *
                        exp(-dτ_λ[n₀] / qp_μN[j])

                end
            #else
            #    iet⁺⁺[i,j,n₁,Δn] =  0.0
            #end
        else
            #@show  qp_μN[i], qp_μN[j]  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)

            if (abs( (qp_μN[i]/qp_μN[j]) - (dτ_λ[n₁]/dτ_λ[n₀]) ) < 1.e-8)
                iet⁺⁺[i,j,n₁,Δn] = 
                (dτ_λ[n₀]/qp_μN[i]) * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * 
                wct2[j] * exp(-dτ_λ[n₀] / qp_μN[j])
            else
                iet⁺⁺[i,j,n₁,Δn] = 
                ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * 
                (1 / ( (qp_μN[i]/qp_μN[j]) - (dτ_λ[n₁]/dτ_λ[n₀]) )) * 
                wct2[j] * 
                (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j]))
            end
        end
        
    else
        ier⁻⁺[i,j,n₁,Δn] = 0.0
        if i==j
            iet⁺⁺[i,j,n₁,Δn] = 0.0
        else
            iet⁺⁺[i,j,n₁,Δn] = 0.0
        end
    end

end

# kernel wrapper:
function get_elem_rt!(RS_type::RRS, 
                        ier⁻⁺, iet⁺⁺, 
                        dτ_λ, ϖ_λ,
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                        qp_μN, wct2)
        @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
        device = devi(architecture(ier⁻⁺))
        aType = array_type(architecture(ier⁻⁺))
        kernel! = get_elem_rt_RRS!(device)
       #@show typeof(i_ref), typeof(ϖ_λ₁λ₀)
       # @show typeof(dτ_λ), typeof(ϖ_λ)
       # @show typeof(qp_μN), typeof(wct2)
       # @show typeof(ier⁻⁺), typeof(iet⁺⁺)
        #@show typeof(fscattRayl), (fscattRayl[1:10]), typeof(aType(ϖ_λ₁λ₀))
        #exit()
        #sleep(10)
       
            #@show minimum(ϖ_λ₁λ₀), maximum(ϖ_λ₁λ₀)
            
        event = kernel!(aType(fscattRayl), 
                    aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
                    i_ref,
                    ier⁻⁺, iet⁺⁺, 
                    dτ_λ, ϖ_λ,
                    aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), 
                    qp_μN, wct2, 
                    ndrange=getKernelDim(RS_type,ier⁻⁺)); 
        #wait(device, event);
        synchronize_if_gpu();
        #for j=1:1:length(qp_μN)
        #    @show minimum(iet⁺⁺[1:3:end,j,200,50]), minimum(ier⁻⁺[1:3:end,j,200,50]) 
        #    @show maximum(iet⁺⁺[1:3:end,j,200,50]), maximum(ier⁻⁺[1:3:end,j,200,50]) 
        #end
        #for j=1:3:length(qp_μN)
        #    @show minimum(iet⁺⁺[1:3:end,j,200,50]), minimum(ier⁻⁺[1:3:end,j,200,50]) 
        #    @show maximum(iet⁺⁺[1:3:end,j,200,50]), maximum(ier⁻⁺[1:3:end,j,200,50]) 
        #end
end

function get_elem_rt!(RS_type::Union{VS_0to1, VS_1to0}, 
    ier⁻⁺, iet⁺⁺, 
    dτ_λ, ϖ_λ,
    Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
    qp_μN, wct2)
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
    device = devi(architecture(ier⁻⁺))
    aType = array_type(architecture(ier⁻⁺))
    kernel! = get_elem_rt_VS!(device)
    #@show typeof(Z⁻⁺_λ₁λ₀), typeof(Z⁺⁺_λ₁λ₀), typeof(ϖ_λ₁λ₀), typeof(i_λ₁λ₀), typeof(i_ref)
    event = kernel!(aType(fscattRayl), 
        aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
        i_ref,
        ier⁻⁺, iet⁺⁺, 
        dτ_λ, ϖ_λ,
        aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), 
        qp_μN, wct2, 
        ndrange=getKernelDim(RS_type,ier⁻⁺)); 
    #wait(device, event);
    synchronize_if_gpu();
end


@kernel function get_elem_rt_VS!(fscattRayl,
                            ϖ_λ₁λ₀, i_λ₁λ₀, i_ref,
                            ier⁻⁺, iet⁺⁺, 
                            dτ_λ, ϖ_λ,
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, wct2)
    i, j, Δn = @index(Global, NTuple) 
    #@unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, dτ₀, dτ₀_λ = RS_type 
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    #dτ₁ = 1 #dummy for now
    #Suniti: require that the incident wavelength is always the first element of 1:nSpec, and all the others belong to the same target VS band
    #Suniti: Then,
    n₀ = 1    
    n₁ = n₀ + i_λ₁λ₀[Δn]  
    if (wct2[j]>1.e-8) 
        
        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,1] = 
                ϖ_λ₁λ₀[Δn] * ϖ_λ[n₀] * fscattRayl[n₀] * Z⁻⁺_λ₁λ₀[i,j] * 
                (1/( (qp_μN[i] / qp_μN[j]) + (dτ_λ[n₁]/dτ_λ[n₀]) )) * 
                (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * wct2[j] 
                    
        if (qp_μN[i] == qp_μN[j])
            # @show i,j
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j       
                if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                    iet⁺⁺[i,j,n₁,1] = 
                        ϖ_λ₁λ₀[Δn] * ϖ_λ[n₀] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
                        (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[i]))/
                        (1 - (dτ_λ[n₁]/dτ_λ[n₀]))  
                else    
                    iet⁺⁺[i,j,n₁,1] = 
                        (dτ_λ[n₀]/ qp_μN[i]) * ϖ_λ₁λ₀[Δn] * ϖ_λ[n₀] * fscattRayl[n₀] * 
                        Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
                        exp(-dτ_λ[n₀] / qp_μN[i])   
                end
            else
                iet⁺⁺[i,j,n₁,1] = 0.0
            end
        else
            #@show  qp_μN[i], qp_μN[j]  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            iet⁺⁺[i,j,n₁,1] = 
                    ϖ_λ₁λ₀[Δn] * ϖ_λ[n₀] * fscattRayl[n₀] * Z⁺⁺_λ₁λ₀[i,j] * 
                    (1 / ( (qp_μN[i]/qp_μN[j]) - (dτ_λ[n₁]/dτ_λ[n₀]) )) * 
                    wct2[j] * 
                    (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j]))
        end
    else
        ier⁻⁺[i,j,n₁,1] = 0.0
        if i==j
            iet⁺⁺[i,j,n₁,1] = 0.0
        else
            iet⁺⁺[i,j,n₁,1] = 0.0
        end
    end
end

function get_elem_rt_SFI!(RS_type::Union{VS_0to1, VS_1to0}, 
                        ieJ₀⁺, ieJ₀⁻, 
                        τ_sum, dτ_λ, ϖ_λ, 
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                        qp_μN, ndoubl,
                        wct02, nStokes,
                        I₀, iμ0,D)
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
    #@show fscattRayl
    device = devi(architecture(ieJ₀⁺))
    aType = array_type(architecture(ieJ₀⁺))
    kernel! = get_elem_rt_SFI_VS!(device)
    #@show typeof(ieJ₀⁺), typeof(τ_sum), typeof(dτ_λ),typeof(wct02), typeof(qp_μN), typeof(dτ_λ) 
    event = kernel!(fscattRayl, aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
    i_ref, ieJ₀⁺, ieJ₀⁻, 
    τ_sum, dτ_λ, 
    aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), aType(F₀),
    qp_μN, ndoubl, wct02, nStokes, 
    I₀, iμ0, D, 
    ndrange=getKernelDimSFI(RS_type,ieJ₀⁻));
    #wait(device, event)
    synchronize_if_gpu();
end

#  TODO: Nov 30, 2021
function get_elem_rt_SFI!(RS_type::RRS, 
                        ieJ₀⁺, ieJ₀⁻, 
                        τ_sum, dτ_λ, ϖ_λ, 
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                        F₀,
                        qp_μN, ndoubl,
                        wct02, nStokes,
                        I₀, iμ0,D)
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
    #@show fscattRayl
    device  = devi(architecture(ieJ₀⁺))
    aType   = array_type(architecture(ieJ₀⁺))
    kernel! = get_elem_rt_SFI_RRS!(device)
    #@show typeof(ieJ₀⁺), typeof(τ_sum), typeof(dτ_λ),typeof(wct02), typeof(qp_μN), typeof(dτ_λ) 
    event = kernel!(aType(fscattRayl), aType(ϖ_λ₁λ₀), aType(i_λ₁λ₀), 
                i_ref, ieJ₀⁺, ieJ₀⁻, 
                τ_sum, dτ_λ, ϖ_λ,
                aType(Z⁻⁺_λ₁λ₀), aType(Z⁺⁺_λ₁λ₀), 
                aType(F₀),
                qp_μN, ndoubl, wct02, nStokes, 
                I₀, iμ0, D, 
                ndrange=getKernelDimSFI(RS_type,ieJ₀⁻));
    
    #wait(device, event)
    synchronize_if_gpu();
    #@show minimum(ieJ₀⁺[1:3:end,1,200,50]), minimum(ieJ₀⁻[1:3:end,1,200,50]) 
    #@show maximum(ieJ₀⁺[1:3:end,1,200,50]), maximum(ieJ₀⁻[1:3:end,1,200,50]) 
end

# only for RRS
@kernel function get_elem_rt_SFI_RRS!(fscattRayl, 
                            ϖ_λ₁λ₀, i_λ₁λ₀, i_ref, 
                            ieJ₀⁺, ieJ₀⁻, 
                            τ_sum, dτ_λ, ϖ_λ,
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, F₀,
                            qp_μN, ndoubl,
                            wct02, nStokes,
                            I₀, iμ0, D)

    # 
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    nMax = length(dτ_λ)
    i, _, n₁, Δn = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    n₀  = n₁ + i_λ₁λ₀[Δn]
    #i_ϖ = i_ref + i_λ₁λ₀[Δn]
    ieJ₀⁺[i, 1, n₁, Δn]=0
    ieJ₀⁻[i, 1, n₁, Δn]=0        
    FT = eltype(I₀)
    if (1 ≤ n₀ ≤ nMax)
         
        Z⁺⁺_I₀ = FT(0.0);
        Z⁻⁺_I₀ = FT(0.0);
        for ii = i_start:i_end
            Z⁺⁺_I₀ += Z⁺⁺_λ₁λ₀[i,ii] * F₀[ii-i_start+1,n₀] #I₀[ii-i_start+1]
            Z⁻⁺_I₀ += Z⁻⁺_λ₁λ₀[i,ii] * F₀[ii-i_start+1,n₀] #I₀[ii-i_start+1] 
        end  
        if (i_start ≤ i ≤ i_end)
            #ctr = i-i_start+1
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
            if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-8
                ieJ₀⁺[i, 1, n₁, Δn] = 
                        ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_I₀ * wct02 *
                        (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[i])) /
                        (1 - (dτ_λ[n₁]/dτ_λ[n₀])) 
            else
                ieJ₀⁺[i, 1, n₁, Δn] = 
                        (dτ_λ[n₀]/ qp_μN[i]) * wct02 * ϖ_λ₁λ₀[Δn] * 
                        fscattRayl[n₀] * 
                        Z⁺⁺_I₀ * 
                        exp(-dτ_λ[n₀] / qp_μN[i])
            end
        else
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
            if (abs( (qp_μN[i]/qp_μN[i_start]) - (dτ_λ[n₁]/dτ_λ[n₀]) ) < 1.e-8)
                ieJ₀⁺[i, 1, n₁, Δn] = 
                (dτ_λ[n₀]/qp_μN[i]) * wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_I₀ * 
                exp(-dτ_λ[n₀] / qp_μN[i_start])
            else
                ieJ₀⁺[i, 1, n₁, Δn] = 
                    wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁺⁺_I₀ * 
                    (1 /( (qp_μN[i]/qp_μN[i_start]) - (dτ_λ[n₁]/dτ_λ[n₀]) ) ) * 
                    (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[i_start]))
            end
        end
        
        #TODO
        #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]                    
        ieJ₀⁻[i, 1, n₁, Δn] = wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl[n₀] * Z⁻⁺_I₀ * 
                (1/( (qp_μN[i] / qp_μN[i_start]) + (dτ_λ[n₁]/dτ_λ[n₀]) )) *
                (1 - exp(-( (dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[i_start]) ) ))  
        
        ieJ₀⁺[i, 1, n₁, Δn] *= exp(-τ_sum[n₀]/qp_μN[i_start]) #correct this to include n₀ap
        ieJ₀⁻[i, 1, n₁, Δn] *= exp(-τ_sum[n₀]/qp_μN[i_start]) 
    end
    if ndoubl >= 1 #double check to make sure this isnt repeated using apply_D
        ieJ₀⁻[i, 1, n₁, Δn] = D[i,i] * ieJ₀⁻[i, 1, n₁, Δn] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end    
    #if ((n₀==840||n₀==850)&&(i==3))       
    #    @show i, n₀, n₁, Δn, ieJ₀⁺[i, 1, n₁, Δn], ieJ₀⁻[i, 1, n₁, Δn]
    #end
    #@show i, n₁, Δn
end

@kernel function apply_D_elemental_RRS!(ndoubl, pol_n, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
    i, j, n₁, Δn = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            ier⁺⁻[i, j, n₁, Δn] = ier⁻⁺[i, j, n₁, Δn]
            iet⁻⁻[i, j, n₁, Δn] = iet⁺⁺[i, j ,n₁, Δn]
        else
            ier⁺⁻[i, j, n₁, Δn] = -ier⁻⁺[i, j, n₁, Δn] 
            iet⁻⁻[i, j, n₁, Δn] = -iet⁺⁺[i, j, n₁, Δn] 
        end
    else
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            ier⁻⁺[i, j, n₁, Δn] = - ier⁻⁺[i, j, n₁, Δn]
        end 
    end
end

@kernel function apply_D_elemental_VS!(ndoubl, 
                                pol_n, 
                                i_λ₁λ₀, 
                                ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)

    i, j, Δn = @index(Global, NTuple)
    n₁ = i_λ₁λ₀[Δn]
    if n₁>0
        if ndoubl < 1
            ii = mod(i, pol_n) 
            jj = mod(j, pol_n) 
            #if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            if (((1<=ii<=2) & (1<=jj<= 2)) | (!(1<=ii<=2) & !(1<=jj<=2))) 
                ier⁺⁻[i, j, n₁, 1] = ier⁻⁺[i, j, n₁, 1]
                iet⁻⁻[i, j, n₁, 1] = iet⁺⁺[i, j ,n₁, 1]
            else
                ier⁺⁻[i, j, n₁, 1] = -ier⁻⁺[i, j, n₁, 1] 
                iet⁻⁻[i, j, n₁, 1] = -iet⁺⁺[i, j, n₁, 1] 
            end
        else
            if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
                ier⁻⁺[i, j, n₁, 1] = - ier⁻⁺[i, j, n₁, 1] 
            end 
        end
    end
end

@kernel function apply_D_elemental_SFI_RRS!(ndoubl, pol_n, ieJ₀⁻)
    i, _, n₁, Δn = @index(Global, NTuple)
          
    if ndoubl>1
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            ieJ₀⁻[i, 1, n₁, Δn] = - ieJ₀⁻[i, 1, n₁, Δn] #this assumes an unpolarized source
        end 
    end
end

@kernel function apply_D_elemental_SFI_VS!(ndoubl, 
        pol_n, 
        i_λ₁λ₀, 
        ieJ₀⁻)
    i, Δn = @index(Global, NTuple)
    #@unpack i_λ₁λ₀ = RS_type
    
    n₁ = i_λ₁λ₀[Δn]
    
    if ndoubl>1
        if (n₁>0)
            if !(1<=mod(i, pol_n)<=2)
                ieJ₀⁻[i, 1, n₁, 1] = - ieJ₀⁻[i, 1, n₁, 1]
            end 
        end
    end
end


function apply_D_matrix_elemental!(RS_type::Union{RRS, RRS_plus}, ndoubl::Int, n_stokes::Int, 
                                    ier⁻⁺::AbstractArray{FT,4}, 
                                    iet⁺⁺::AbstractArray{FT,4}, 
                                    ier⁺⁻::AbstractArray{FT,4}, 
                                    iet⁻⁻::AbstractArray{FT,4}) where {FT}
    device = devi(architecture(ier⁻⁺))
    applyD_kernel! = apply_D_elemental_RRS!(device)
    event = applyD_kernel!(ndoubl,
        n_stokes, 
        ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, 
        ndrange=size(ier⁻⁺));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
                        ndoubl::Int, n_stokes::Int, 
                        ier⁻⁺::AbstractArray{FT,4}, 
                        iet⁺⁺::AbstractArray{FT,4}, 
                        ier⁺⁻::AbstractArray{FT,4}, 
                        iet⁻⁻::AbstractArray{FT,4}) where {FT}
    
    device = devi(architecture(ier⁻⁺))
    applyD_kernel! = apply_D_elemental_VS!(device)
    event = applyD_kernel!(ndoubl,
                    n_stokes, RS_type.i_λ₁λ₀_all,
                    ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, 
                    ndrange=getKernelDim(RS_type,ier⁻⁺,RS_type.i_λ₁λ₀_all));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental_SFI!(RS_type::Union{RRS, RRS_plus},
        ndoubl::Int, n_stokes::Int, ieJ₀⁻::AbstractArray{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        #@show "here 1.1"
        device = devi(architecture(ieJ₀⁻))
        #@show "here 1.2"
        applyD_kernel! = apply_D_elemental_SFI_RRS!(device)
        #@show "here 1.3", RS_type
        event = applyD_kernel!(ndoubl,
                                n_stokes, 
                                ieJ₀⁻, 
                                ndrange=size(ieJ₀⁻));
        #@show "here 1.4"
        #wait(device, event);
        synchronize();
        return nothing
    end
end

function apply_D_matrix_elemental_SFI!(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                    ndoubl::Int, n_stokes::Int, ieJ₀⁻::AbstractArray{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        #@show "here 1.1"
        device = devi(architecture(ieJ₀⁻))
        #@show "here 1.2"
        applyD_kernel! = apply_D_elemental_SFI_VS!(device)
        #@show "here 1.3", RS_type
        event = applyD_kernel!(ndoubl,
                            n_stokes, 
                            RS_type.i_λ₁λ₀_all,    
                            ieJ₀⁻, 
                            ndrange = getKernelDimSFI(RS_type,ieJ₀⁻,RS_type.i_λ₁λ₀_all));
                            #ndrange=size(ieJ₀⁻));
        #@show "here 1.4"
        #wait(device, event);
        synchronize();
        return nothing
    end
end