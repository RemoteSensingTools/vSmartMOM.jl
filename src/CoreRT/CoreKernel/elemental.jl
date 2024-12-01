#=
 
This file contains RT elemental-related functions
 
=#

"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,       #{FT2,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ_λ::AbstractArray{FT,1},   # ϖ_λ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ::FT,                      # ϖ: single scattering albedo of elemental layer (no trace gas absorption included)
                            Z⁺⁺::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺::AbstractArray{FT,2},   # Z matrix
                            F₀::AbstractArray{FT,2},    # Stokes vector of solar/stellar irradiance
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            I_static,
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀ = quad_points
    #@unpack ϖ_Cabannes = RS_type
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    # @show Array(τ_sum)[1], Array(dτ_λ)[1], Array(ϖ_λ)[1], Array(Z⁺⁺)[1,1]
    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))

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
            # Version 2: More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
            # Version 2: with absorption in batch mode, low tau_scatt but higher tau_total, needs different equations
            kernel! = get_elem_rt!(device)
            event = kernel!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, F₀,
                qp_μN, wct2, ndrange=size(r⁻⁺)); 
            #wait(device, event)
            synchronize_if_gpu()

            if SFI
                kernel! = get_elem_rt_SFI!(device)
                event = kernel!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, F₀, qp_μN, ndoubl, wct02, pol_type.n, arr_type(pol_type.I₀), iμ₀, D, ndrange=size(J₀⁺))
                #wait(device, event)
                synchronize_if_gpu()
            end
        end

        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, J₀⁻)
        end      
    else 
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,#{FT2,1}, #Suniti
                            dτ::AbstractArray,
                            F₀::AbstractArray,#{FT,2},    # Stokes vector of solar/stellar irradiance
                            computed_layer_properties,
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack qp_μ, iμ₀, wt_μN, qp_μN = quad_points
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    #@unpack ϖ_Cabannes = RS_type
    #@show architecture
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    qp_μN = arr_type(qp_μN)
    wt_μN = arr_type(wt_μN)
    τ_sum = arr_type(τ_sum)
    I₀    = arr_type(pol_type.I₀)
    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)
    #@show typeof(ϖ),typeof(dτ),typeof(Z⁻⁺),typeof(Z⁺⁺) 
    #ϖ   = arr_type(ϖ);
    #dτ  = arr_type(dτ);
    #Z⁻⁺ = arr_type(Z⁻⁺);
    #Z⁺⁺ = arr_type(Z⁺⁺);
    #@show size(Z⁻⁺), size(ϖ)
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
        event = kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, qp_μN, wct2, ndrange=size(r⁻⁺)); 
        #@show "Stop event"
        #wait(device, event)
        synchronize_if_gpu()

        if SFI
            kernel! = get_elem_rt_SFI!(device)
            #@show size(F₀)
            event = kernel!(J₀⁺, J₀⁻, ϖ, dτ, arr_type(τ_sum), Z⁻⁺, Z⁺⁺, 
            arr_type(F₀), qp_μN, ndoubl, wct02, pol_type.n, I₀, iμ₀, D, ndrange=size(J₀⁺))
            #wait(device, event)
        end
        #ii = pol_type.n*(iμ0-1)+1
        #@show 'B',iμ0,  r⁻⁺[1,ii,1]/(J₀⁻[1,1,1]*wt_μ[iμ0]), r⁻⁺[1,ii,1], J₀⁻[1,1,1]*wt_μ[iμ0], J₀⁺[1,1,1]*wt_μ[iμ0]
        synchronize_if_gpu()
        
        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, J₀⁻)
        end      
    else
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

@kernel function get_elem_rt!(r⁻⁺, t⁺⁺, 
                        ϖ_λ, dτ_λ, 
                        Z⁻⁺, Z⁺⁺, 
                        qp_μN, wct) 
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    if (wct[j]>1.e-8) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            #Z⁻⁺[i,j] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))) 
                    
        if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                    #(1 + ϖ_λ[n] * Z⁺⁺[i,i] * (dτ_λ[n] / qp_μN[i]) * wct[i])
            else
                #t⁺⁺[i,j,n] = 0.0
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
                    
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                #Z⁺⁺[i,j] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j])) 
        end
    else
        r⁻⁺[i,j,n] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i]) #Suniti
        else
            t⁺⁺[i,j,n] = 0.0
        end
    end
    nothing
end

@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, 
                ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, F₀,
                qp_μN, ndoubl, wct02, nStokes,
                I₀, iμ0, D)
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    J₀⁺[i, 1, n]=0
    J₀⁻[i, 1, n]=0
    n2=1
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    Z⁺⁺_I₀ = FT(0.0);
    Z⁻⁺_I₀ = FT(0.0);
    
    for ii = i_start:i_end
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n2] #I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n2] #I₀[ii-i_start+1] 
    end

    if (i>=i_start) && (i<=i_end)
        ctr = i-i_start+1
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (dτ_λ[n] / qp_μN[i]) * exp(-dτ_λ[n] / qp_μN[i])
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
    end
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))

    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end  
    #if (n==840||n==850)    
    #    @show i, n, J₀⁺[i, 1, n], J₀⁻[i, 1, n]      
    #end
    nothing
end

@kernel function apply_D_elemental!(ndoubl, pol_n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    i, j, n = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod(i, pol_n) 
        jj = mod(j, pol_n) 
        #if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
        if (((1<=ii<=2) & (1<=jj<=2)) | (!(1<=ii<=2) & !(1<=jj<=2))) 
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] 
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] 
        end
    else
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
        end 
    end
    nothing
end

@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if !(1<=mod(i, pol_n)<=2) #mod(i, pol_n) > 2
            J₀⁻[i, 1, n] = - J₀⁻[i, 1, n]
        end 
    end
    nothing
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
    device = devi(architecture(r⁻⁺))
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(architecture(J₀⁻))
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end