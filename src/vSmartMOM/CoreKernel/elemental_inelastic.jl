#=
 
This file contains RT elemental-related functions
 
=#
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
function elemental_inelastic!(RS_type,
                            pol_type, SFI::Bool, 
                            τ_sum::AbstractArray{FT,1},
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            Z⁺⁺_λ₁λ₀::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺_λ₁λ₀::AbstractArray{FT,2}, 
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            I_static,
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    
    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺_λ₁λ₀, (size(Z⁺⁺_λ₁λ₀,1), size(Z⁺⁺_λ₁λ₀,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺_λ₁λ₀, (size(Z⁺⁺_λ₁λ₀,1), size(Z⁺⁺_λ₁λ₀,2),1))

    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))
    I₀_NquadN = arr_type(zeros(FT,size(qp_μN,1))); #incident irradiation
    i_start   = pol_type.n*(iμ₀-1) + 1 
    i_end     = pol_type.n*iμ₀
    I₀_NquadN[iμ₀Nstart:i_end] = pol_type.I₀

    device = devi(architecture)

    # If in scattering mode:
    if scatter
   
        NquadN = length(qp_μN)

        # Needs explanation still, different weights: 
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        # scalars
        #wct0  = m == 0 ? FT(0.50) * ϖ * dτ     : FT(0.25) * ϖ * dτ 
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        # vectors
        #wct   = m == 0 ? FT(0.50) * ϖ * wt_μN  : FT(0.25) * ϖ * wt_μN
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4

        # Get the diagonal matrices first
        #d_qp  = Diagonal(1 ./ qp_μN)
        #d_wct = Diagonal(wct2)

        # Calculate r⁻⁺ and t⁺⁺
        #Version 2: More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
        # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
        
        kernel! = get_elem_rt!(device)
        @show getKernelDim(RS_type,ier⁻⁺)
        @show size(qp_μN), size(dτ_λ), size(Z⁻⁺_λ₁λ₀)
        @unpack ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
        @show size(ϖ_λ₁λ₀), size(i_λ₁λ₀), size(i_ref)
        event = kernel!(RS_type, 
                        ier⁻⁺, iet⁺⁺, 
                        dτ, dτ_λ, 
                        Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                        qp_μN, wct2, 
                        ndrange=getKernelDim(RS_type,ier⁻⁺)); 

        wait(device, event)
        synchronize_if_gpu()

        if SFI
            kernel! = get_elem_rt_SFI!(device)
            event = kernel!(RS_type, 
                            ieJ₀⁺, ieJ₀⁻, 
                            τ_sum, dτ, dτ_λ, 
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, ndoubl,wct02, 
                            pol_type.n, 
                            arr_type(pol_type.I₀), 
                            iμ₀, D, 
                            ndrange=getKernelDimSFI(RS_type,ieJ₀⁻));
            wait(device, event)
        end
            #ii = pol_type.n*(iμ0-1)+1
            #@show 'B',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
            synchronize_if_gpu()

        # Apply D Matrix
        apply_D_matrix_elemental!(RS_type, ndoubl, pol_type.n, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(RS_type, ndoubl, pol_type.n, ieJ₀⁻)
        end      
    else 
        # Note: τ is not defined here
        iet⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        iet⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end
#Suniti: is there a way to pass information like ϖ_λ₁λ₀, i_λ₁λ₀, i_ref, etc. along with RS_type? So that they can be retrieved as RSS.ϖ_λ₁λ₀ for example?
@kernel function get_elem_rt!(RS_type::RRS, 
                            ier⁻⁺, iet⁺⁺, 
                            dτ, dτ_λ, 
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, wct2)

    i, j, n₁, Δn = @index(Global, NTuple)
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type
    nMax = length(dτ_λ) 
    # n₁ covers the full range of wavelengths, while n₀ = n₁+Δn only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    n₀  = n₁ + i_λ₁λ₀[Δn]
    i_ϖ = i_ref + i_λ₁λ₀[Δn]
    #@show   n₀ , i_ϖ 
    if (wct2[j]>1.e-8) & (1 ≤ n₀ ≤ nMax)

        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,Δn] = fscattRayl * ϖ_λ₁λ₀[Δn] * Z⁻⁺_λ₁λ₀[i,j] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * 
            (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * wct2[j] 
                    
        if (qp_μN[i] == qp_μN[j])
            # @show i,j
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j       
                if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                    iet⁺⁺[i,j,n₁,Δn] = 
                        ϖ_λ₁λ₀[Δn] * fscattRayl * dτ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
                        ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) 
                        
                else    
                    iet⁺⁺[i,j,n₁,Δn] = 
                        ϖ_λ₁λ₀[Δn] * fscattRayl * dτ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] *
                        exp(-dτ_λ[n₀] / qp_μN[j])/ qp_μN[j]
                end
            else
                iet⁺⁺[i,j,n₁,Δn] = 0.0
            end
        else
            #@show  qp_μN[i], qp_μN[j]  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            iet⁺⁺[i,j,n₁,Δn] = 
                ϖ_λ₁λ₀[Δn] * fscattRayl * Z⁺⁺_λ₁λ₀[i,j] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct2[j] * 
                (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j]))
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

@kernel function get_elem_rt!(RS_type::Union{VS_0to1, VS_1to0}, 
                            ier⁻⁺, iet⁺⁺, 
                            dτ₁, dτ_λ, 
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, wct2)
    i, j, Δn = @index(Global, NTuple) 
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, dτ₀, dτ₀_λ = RS_type 
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    
    n₁ = i_λ₁λ₀[Δn]  
    if (wct2[j]>1.e-8) 

        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,1] = 
                ϖ_λ₁λ₀[n₁] * fscattRayl * (dτ₀/dτ₁) * Z⁻⁺_λ₁λ₀[i,j] * 
                (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[j]*dτ₁)) * wct2[j] * 
                (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ₀_λ / qp_μN[j])))) 
                    
        if (qp_μN[i] == qp_μN[j])
            # @show i,j
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j       
                if abs(dτ₀_λ-dτ_λ[n₁])>1.e-6
                    iet⁺⁺[i,j,n₁,1] = 
                            ((exp(-dτ₀_λ / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ₀_λ)) * 
                            ϖ_λ₁λ₀[n₁] * fscattRayl * dτ₀ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i]
                else    
                    iet⁺⁺[i,j,n₁,1] = 
                            ϖ_λ₁λ₀[n₁] * fscattRayl * dτ₀ * Z⁺⁺_λ₁λ₀[i,i] * wct2[i] * 
                            exp(-dτ₀_λ / qp_μN[j])/ qp_μN[j]
                end
            else
                iet⁺⁺[i,j,n₁,1] = 0.0
            end
        else
            #@show  qp_μN[i], qp_μN[j]  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            iet⁺⁺[i,j,n₁,1] = 
                    ϖ_λ₁λ₀[n₁] * fscattRayl * (dτ₀/dτ₁) * Z⁺⁺_λ₁λ₀[i,j] * 
                    (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[j]*dτ₁)) * wct2[j] * 
                    (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ₀_λ / qp_μN[j]))
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

#  TODO: Nov 30, 2021
@kernel function get_elem_rt_SFI!(RS_type::Union{VS_0to1, VS_1to0}, 
                            ieJ₀⁺, ieJ₀⁻, 
                            τ_sum, dτ₁, dτ_λ, 
                            Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
                            qp_μN, ndoubl,wct02, 
                            nStokes, I₀, iμ0,D)
    
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, dτ₀, dτ₀_λ = RS_type 

    i, _, Δn = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    n₁ =  i_λ₁λ₀[Δn]
      
    #if (wct2[j]>1.e-8) 
    
        FT = eltype(I₀)
        ieJ₀⁺[i, 1, n₁, 1]=0
        ieJ₀⁻[i, 1, n₁, 1]=0
    
        Z⁺⁺_I₀ = FT(0.0);
        Z⁻⁺_I₀ = FT(0.0);
        for ii = i_start:i_end
            Z⁺⁺_I₀ += Z⁺⁺_λ₁λ₀[i,ii] * I₀[ii-i_start+1]
            Z⁻⁺_I₀ += Z⁻⁺_λ₁λ₀[i,ii] * I₀[ii-i_start+1] 
        end
    
    if (i>=i_start) && (i<=i_end)
        #ctr = i-i_start+1
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        if abs(dτ₀_λ-dτ_λ[n₁])>1.e-6
            ieJ₀⁺[i, 1, n₁, 1] = 
                    ((exp(-dτ₀_λ / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ₀_λ)) * 
                    ϖ_λ₁λ₀[n₁] * fscattRayl * dτ₀ * Z⁺⁺_I₀ * wct02
        else
            ieJ₀⁺[i, 1, n₁, 1] = 
                    wct02 * ϖ_λ₁λ₀[n₁] * fscattRayl * Z⁺⁺_I₀ * 
                    (dτ₀ / qp_μN[j]) * exp(-dτ₀_λ / qp_μN[j])
        end
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        ieJ₀⁺[i, 1, n₁, 1] = 
            wct02 * ϖ_λ₁λ₀[n₁] * fscattRayl * (dτ₀/dτ₁) * Z⁺⁺_I₀ * 
            (qp_μN[i_start]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[i_start]*dτ₁)) * 
            (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ₀_λ / qp_μN[i_start]))
    end
    #TODO
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]                    
    ieJ₀⁻[i, 1, n₁, 1] = 
            wct02 * ϖ_λ₁λ₀[n₁] * fscattRayl * (dτ₀/dτ₁) * Z⁻⁺_I₀ * 
            (qp_μN[i_start]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[i_start]*dτ₁)) * 
            (1 - exp(-( (dτ_λ[n₁] / qp_μN[i]) + (dτ₀_λ / qp_μN[i_start]) )))

    ieJ₀⁺[i, 1, n₁, 1] *= exp(-τ_sum[n]/qp_μN[i_start])
    ieJ₀⁻[i, 1, n₁, 1] *= exp(-τ_sum[n]/qp_μN[i_start])

    if ndoubl >= 1
        ieJ₀⁻[i, 1, n₁, 1] = D[i,i]*ieJ₀⁻[i, 1, n₁, 1] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end        
end

#  TODO: Nov 30, 2021
@kernel function get_elem_rt_SFI!(RS_type::RRS, 
    ieJ₀⁺, ieJ₀⁻, 
    τ_sum, dτ, dτ_λ, 
    Z⁻⁺_λ₁λ₀, Z⁺⁺_λ₁λ₀, 
    qp_μN, ndoubl,
    wct02, nStokes,
    I₀, iμ0,D)

    @unpack fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, i_ref = RS_type 
    
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    nMax = length(dτ_λ)
    i, _, n₁, Δn = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    n₀  = n₁ + i_λ₁λ₀[Δn]
    i_ϖ = i_ref + i_λ₁λ₀[Δn]     
    FT = eltype(I₀)
    if (1 ≤ n₀ ≤ nMax)
        ieJ₀⁺[i, 1, n₁, Δn]=0
        ieJ₀⁻[i, 1, n₁, Δn]=0    
        Z⁺⁺_I₀ = FT(0.0);
        Z⁻⁺_I₀ = FT(0.0);
        for ii = i_start:i_end
            Z⁺⁺_I₀ += Z⁺⁺_λ₁λ₀[i,ii] * I₀[ii-i_start+1]
            Z⁻⁺_I₀ += Z⁻⁺_λ₁λ₀[i,ii] * I₀[ii-i_start+1] 
        end  
        if (i_start ≤ i ≤ i_end)
            #ctr = i-i_start+1
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
            if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                ieJ₀⁺[i, 1, n₁, Δn] = 
                        ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) * 
                        ϖ_λ₁λ₀[Δn] * fscattRayl * dτ * Z⁺⁺_I₀ * wct02
            else
                ieJ₀⁺[i, 1, n₁, Δn] = 
                        wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl * Z⁺⁺_I₀ * 
                        (dτ / qp_μN[i_start]) * exp(-dτ_λ[n₀] / qp_μN[i_start])
            end
        else
            # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
            ieJ₀⁺[i, 1, n₁, Δn] = 
                    wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl * Z⁺⁺_I₀ * 
                    (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * 
                    (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[i_start]))
        end
        #TODO
        #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]                    
        ieJ₀⁻[i, 1, n₁, Δn] = wct02 * ϖ_λ₁λ₀[Δn] * fscattRayl * Z⁻⁺_I₀ * 
                (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * 
                (1 - exp(-( (dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[i_start]) )))  
        ieJ₀⁺[i, 1, n₁, Δn] *= exp(-τ_sum[n₀]/qp_μN[i_start]) #correct this to include n₀ap
        ieJ₀⁻[i, 1, n₁, Δn] *= exp(-τ_sum[n₀]/qp_μN[i_start]) 
    end
    if ndoubl >= 1 #double check to make sure this isnt repeated using apply_D
        ieJ₀⁻[i, 1, n₁, Δn] = D[i,i] * ieJ₀⁻[i, 1, n₁, Δn] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end           
end

@kernel function apply_D_elemental!(RS_type::RRS, ndoubl, pol_n, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
    i, j, n₁, n₀ = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            ier⁺⁻[i, j, n₁, n₀] = ier⁻⁺[i, j, n₁, n₀]
            iet⁻⁻[i, j, n₁, n₀] = iet⁺⁺[i, j ,n₁, n₀]
        else
            ier⁺⁻[i, j, n₁, n₀] = -ier⁻⁺[i, j, n₁, n₀] 
            iet⁻⁻[i, j, n₁, n₀] = -iet⁺⁺[i, j, n₁, n₀] 
        end
    else
        if mod(i, pol_n) > 2
            ier⁻⁺[i, j, n₁, n₀] = - ier⁻⁺[i, j, n₁, n₀]
        end 
    end
end

@kernel function apply_D_elemental_SFI!(RS_type::RRS, ndoubl, pol_n, ieJ₀⁻)
    i, _, n₁, n₀ = @index(Global, NTuple)
    
    if ndoubl>1
        if mod(i, pol_n) > 2
            ieJ₀⁻[i, 1, n₁, n₀] = - ieJ₀⁻[i, 1, n₁, n₀]
        end 
    end
end

@kernel function apply_D_elemental!(RS_type::Union{VS_0to1, VS_1to0}, 
                        ndoubl, pol_n, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)

    i, j, Δn = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type
    
    n₁ = i_λ₁λ₀[Δn]
    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            ier⁺⁻[i, j, n₁, 1] = ier⁻⁺[i, j, n₁, 1]
            iet⁻⁻[i, j, n₁, 1] = iet⁺⁺[i, j ,n₁, 1]
        else
            ier⁺⁻[i, j, n₁, 1] = -ier⁻⁺[i, j, n₁, 1] 
            iet⁻⁻[i, j, n₁, 1] = -iet⁺⁺[i, j, n₁, 1] 
        end
    else
        if mod(i, pol_n) > 2
            ier⁻⁺[i, j, n₁, 1] = - ier⁻⁺[i, j, n₁, 1]
        end 
    end
end

@kernel function apply_D_elemental_SFI!(RS_type::Union{VS_0to1, VS_1to0}, 
                        ndoubl, pol_n, ieJ₀⁻)
    i, _, Δn = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type
    
    n₁ = i_λ₁λ₀[Δn]
    if ndoubl>1
        if mod(i, pol_n) > 2
            ieJ₀⁻[i, 1, n₁, 1] = - ieJ₀⁻[i, 1, n₁, 1]
        end 
    end
end

function apply_D_matrix_elemental!(RS_type::RRS, ndoubl::Int, n_stokes::Int, 
                                    ier⁻⁺::CuArray{FT,4}, 
                                    iet⁺⁺::CuArray{FT,4}, 
                                    ier⁺⁻::CuArray{FT,4}, 
                                    iet⁻⁻::CuArray{FT,4}) where {FT}
    device = devi(Architectures.GPU())
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(RS_type, ndoubl,n_stokes, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=size(ier⁻⁺));
    wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental!(RS_type::RRS, ndoubl::Int, n_stokes::Int, 
                                    ier⁻⁺::Array{FT,4}, 
                                    iet⁺⁺::Array{FT,4}, 
                                    ier⁺⁻::Array{FT,4}, 
                                    iet⁻⁻::Array{FT,4}) where {FT}
    device = devi(Architectures.CPU())
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(RS_type, ndoubl,n_stokes, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=size(ier⁻⁺));
    wait(device, event);
    return nothing
end

function apply_D_matrix_elemental_SFI!(RS_type::RRS,ndoubl::Int, n_stokes::Int, J₀⁻::CuArray{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(Architectures.GPU())
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(RS_type::RRS,ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(device, event);
        synchronize();
        return nothing
    end
end
    
function apply_D_matrix_elemental_SFI!(RS_type::RRS,ndoubl::Int, n_stokes::Int, J₀⁻::Array{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(RS_type::RRS,ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(device, event);
        return nothing
    end
end