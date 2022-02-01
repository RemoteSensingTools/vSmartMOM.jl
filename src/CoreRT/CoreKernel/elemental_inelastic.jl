#=
 
This file contains RT elemental-related functions
 
=#

"Elemental single-scattering layer for RRS"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,#{FT2,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ_λ::AbstractArray{FT,1},   # ϖ_λ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ::FT,                      # ϖ: single scattering albedo of elemental layer (no trace gas absorption included)
                                                        # Rayleigh_XS/(Raman_XS+Rayleigh_XS)
                            ϖ_λ₀λ₁::AbstractArray{FT,2},# Raman_XS/(Raman_XS+Rayleigh_XS)
                            Z⁺⁺::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺::AbstractArray{FT,2}, 
                            Z⁺⁺_λ₀λ₁::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺_λ₀λ₁::AbstractArray{FT,2}, 
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::AddedLayer{FT}, 
                            I_static,
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    
    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺_λ₀λ₁, (size(Z⁺⁺_λ₀λ₁,1), size(Z⁺⁺_λ₀λ₁,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺_λ₀λ₁, (size(Z⁺⁺_λ₀λ₁,1), size(Z⁺⁺_λ₀λ₁,2),1))

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
        wct0  = m == 0 ? FT(0.50) * ϖ * dτ     : FT(0.25) * ϖ * dτ 
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        # vectors
        wct   = m == 0 ? FT(0.50) * ϖ * wt_μN  : FT(0.25) * ϖ * wt_μN
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4

        # Get the diagonal matrices first
        d_qp  = Diagonal(1 ./ qp_μN)
        d_wct = Diagonal(wct)

        # Calculate r⁻⁺ and t⁺⁺
        
        # Version 1: no absorption in batch mode (initiation of a single scattering layer with no or low absorption)
        if false #maximum(dτ_λ) < 0.0001   
            # R⁻⁺₀₁(λ) = M⁻¹(0.5ϖₑ(λ)Z⁻⁺C)δ (See Eqs.7 in Raman paper draft)
            r⁻⁺[:,:,:] .= d_qp * Z⁻⁺ * (d_wct * dτ)
            # T⁺⁺₀₁(λ) = {I-M⁻¹[I - 0.5*ϖₑ(λ)Z⁺⁺C]}δ (See Eqs.7 in Raman paper draft)
            t⁺⁺[:,:,:] .= I_static - (d_qp * ((I_static - Z⁺⁺ * d_wct) * dτ))
            if SFI
                # Reminder: Add equation here what it does
                expk = exp.(-τ_sum/qp_μ[iμ₀]) #exp(-τ(z)/μ₀)
                # J₀⁺ = 0.5[1+δ(m,0)]M⁻¹ϖₑ(λ)Z⁺⁺τI₀exp(-τ(z)/μ₀)
                J₀⁺[:,1,:] .= ((d_qp * Z⁺⁺ * I₀_NquadN) * wct0) .* expk'
                # J₀⁻ = 0.5[1+δ(m,0)]M⁻¹ϖₑ(λ)Z⁻⁺τI₀exp(-τ(z)/μ₀)
                J₀⁻[:,1,:] .= ((d_qp * Z⁻⁺ * I₀_NquadN) * wct0) .* expk'
              
            end
        else 
            #Version 2: More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
            # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
            kernel! = get_elem_rt!(device)
            event = kernel!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2, ndrange=size(ier⁻⁺)); 
            wait(device, event)
            synchronize_if_gpu()

            if SFI
                kernel! = get_elem_rt_SFI!(device)
                event = kernel!(ieJ₀⁺, ieJ₀⁻, ϖ_λ, ϖ_λ₀λ₁, dτ_λ, τ_sum, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, ndoubl, wct02, pol_type.n, arr_type(pol_type.I₀), iμ₀, D, ndrange=size(J₀⁺))
                wait(device, event)
            end
            #ii = pol_type.n*(iμ0-1)+1
            #@show 'B',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
            synchronize_if_gpu()
        end
        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, ieJ₀⁻)
        end      
    else 
        # Note: τ is not defined here
        iet⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        iet⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

@kernel function get_elem_rt!(ier⁻⁺, iet⁺⁺, ϖ_λ, ϖ_λ₀λ₁, dτ₀, dτ₁, dτ_λ, Z⁻⁺_λ₀λ₁, Z⁺⁺_λ₀λ₁, qp_μN, wct2)
    i, j, n₁, n₀ = @index(Global, NTuple) 
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    if (wct2[j]>1.e-8) 
        # dτ₀, dτ₁ are the purely scattering (elastic+inelastic) molecular elemental 
        # optical thicknesses at wavelengths λ₀ and λ₁
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        ier⁻⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁻⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[j]*dτ₁)) * (1 - exp(-((dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[j])))) * (wct2[j]) 
                    
        if (qp_μN[i] == qp_μN[j])
            # @show i,j
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j       
                if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
                    iet⁺⁺[i,j,n₁,n₀] = ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) * ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i]
                else    
                    iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_λ₀λ₁[i,i] * wct2[i] * exp(-dτ_λ[n₀] / qp_μN[j])/ qp_μN[j]
                end
            else
                iet⁺⁺[i,j,n₁,n₀] = 0.0
            end
        else
            #@show  qp_μN[i], qp_μN[j]  
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            iet⁺⁺[i,j,n₁,n₀] = ϖ_λ₀λ₁[n₁,n₀] * (dτ₀/dτ₁) * Z⁺⁺_λ₀λ₁[i,j] * (qp_μN[j]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[j]*dτ₁)) * (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[j])) * wct2[j]
        end
    else
        ier⁻⁺[i,j,n₁,n₀] = 0.0
        if i==j
            iet⁺⁺[i,j,n₁,n₀] = 0.0
        else
            iet⁺⁺[i,j,n₁,n₀] = 0.0
        end
    end
end

#  TODO: Nov 30, 2021
@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, qp_μN, ndoubl, wct02, nStokes ,I₀, iμ0, D)
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n₁, n₀ = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    # let n₁ cover the full range of wavelengths, while n₀ only includes wavelengths at intervals 
    # that contribute significantly enough to inelastic scattering, so that n₀≪n₁ 
    FT = eltype(I₀)
    J₀⁺[i, 1, n₁, n₀]=0
    J₀⁻[i, 1, n₁, n₀]=0

    
    Z⁺⁺_I₀ = FT(0.0);
    Z⁻⁺_I₀ = FT(0.0);
    for ii = i_start:i_end
        Z⁺⁺_I₀ += Z⁺⁺_λ₀λ₁[i,ii] * I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺_λ₀λ₁[i,ii] * I₀[ii-i_start+1] 
    end
    
    if (i>=i_start) && (i<=i_end)
        #ctr = i-i_start+1
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        if abs(dτ_λ[n₀]-dτ_λ[n₁])>1.e-6
            J₀⁺[i, 1, n₁, n₀] = ((exp(-dτ_λ[n₀] / qp_μN[i]) - exp(-dτ_λ[n₁] / qp_μN[i]))/(dτ_λ[n₁]-dτ_λ[n₀])) * ϖ_λ₀λ₁[n₁,n₀] * dτ₀ * Z⁺⁺_I₀ * wct02
        else
            J₀⁺[i, 1, n₁, n₀] = wct02 * ϖ_λ₁λ₀[n₁, n₀] * Z⁺⁺_I₀ * (dτ₀[n] / qp_μN[j]) * exp(-dτ_λ[n₀] / qp_μN[j])
        end
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        J₀⁺[i, 1, n₁, n₀] = wct02 * ϖ_λ₁λ₀[n₁, n₀] * (dτ₀/dτ₁) * Z⁺⁺_I₀ * (qp_μN[i_start]*dτ₁ / (qp_μN[i]*dτ₀ - qp_μN[i_start]*dτ₁)) * (exp(-dτ_λ[n₁] / qp_μN[i]) - exp(-dτ_λ[n₀] / qp_μN[i_start]))
    end
    #TODO
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]                    
    J₀⁻[i, 1, n₁, n₀] = wct02 * ϖ_λ₁λ₀[n₁, n₀] * (dτ₀/dτ₁) * Z⁻⁺_I₀ * (qp_μN[i_start]*dτ₁ / (qp_μN[i]*dτ₀ + qp_μN[i_start]*dτ₁)) * (1 - exp(-( (dτ_λ[n₁] / qp_μN[i]) + (dτ_λ[n₀] / qp_μN[i_start]) )))

    J₀⁺[i, 1, n₁, n₀] *= exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁻[i, 1, n₁, n₀] *= exp(-τ_sum[n]/qp_μN[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n₁, n₀] = D[i,i]*J₀⁻[i, 1, n₁, n₀] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end        
end

@kernel function apply_D_elemental!(ndoubl, pol_n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    i, j, n = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            r⁺⁻[i, j, n₁, n₀] = r⁻⁺[i, j, n₁, n₀]
            t⁻⁻[i, j, n₁, n₀] = t⁺⁺[i, j ,n₁, n₀]
        else
            r⁺⁻[i, j, n₁, n₀] = -r⁻⁺[i, j, n₁, n₀] 
            t⁻⁻[i, j, n₁, n₀] = -t⁺⁺[i, j, n₁, n₀] 
        end
    else
        if mod(i, pol_n) > 2
            r⁻⁺[i, j, n₁, n₀] = - r⁻⁺[i, j, n₁, n₀]
        end 
    end
end

@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if mod(i, pol_n) > 2
            J₀⁻[i, 1, n₁, n₀] = - J₀⁻[i, 1, n₁, n₀]
        end 
    end
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}) where {FT}
    device = devi(Architectures.GPU())
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
    wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}) where {FT}
    device = devi(Architectures.CPU())
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
    wait(device, event);
    return nothing
end

function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, J₀⁻::CuArray{FT,3}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(Architectures.GPU())
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(device, event);
        synchronize();
        return nothing
    end
end
    
function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, J₀⁻::Array{FT,3}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        wait(device, event);
        return nothing
    end
end
