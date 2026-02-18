#=
 
This file contains RT elemental-related functions
 
=#
"Elemental single-scattering layer"
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,#{FT2,1}, #Suniti
                            dτ::AbstractArray,
                            computed_layer_properties::CoreDirectionalScatteringOpticalProperties,
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2,M}
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer
    @unpack qp_μ, iμ₀, wt_μN, qp_μN = quad_points
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺, G = computed_layer_properties
    
    arr_type = array_type(architecture)
    
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    I₀    = arr_type(pol_type.I₀)
    D     = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)
    #@show maximum(Array(ϖ)), maximum(Array(dτ))
    # If in scattering mode:
    if scatter
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = m == 0 ? FT(0.50)              : FT(0.25)
        wct2  = m == 0 ? wt_μN/2               : wt_μN/4
 
        # More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
        # with absorption in batch mode, low tau_scatt but higher tau_total, needs exact equations
        kernel! = get_canopy_elem_rt!(device)
        event = kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, G, Z⁻⁺, Z⁺⁺, qp_μN, wct2, ndrange=size(r⁻⁺)); 
        #wait(device, event)
        synchronize_if_gpu()
        #@show G
        # SFI part
        kernel! = get_canopy_elem_rt_SFI!(device)
        event = kernel!(j₀⁺, j₀⁻, ϖ, dτ, arr_type(τ_sum), G, Z⁻⁺, Z⁺⁺, qp_μN, ndoubl, wct02, pol_type.n, I₀, iμ₀, D, ndrange=size(j₀⁺))
        #wait(device, event)
        synchronize_if_gpu()
        
        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        # apply D matrix for SFI
        apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, j₀⁻)   
    else
        # Note: τ is not defined here
        t⁺⁺ .= Diagonal{exp.(-τ*G ./ qp_μN)}
        t⁻⁻ .= Diagonal{exp.(-τ*G ./ qp_μN)}
    end    
end

@kernel function get_canopy_elem_rt!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, G, Z⁻⁺, Z⁺⁺, μ, wct) 
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    if (wct[j]>1.e-8) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ

        r⁻⁺[i,j,n] = 
            ϖ_λ[n] *  Z⁻⁺[i,j,n2] * 
            (μ[j] / (μ[i]*G[j] + μ[j]*G[i])) * wct[j] * 
            (1 - exp(-dτ_λ[n] * ((G[i] / μ[i]) + (G[j] / μ[j]))))
                      
        if (μ[i] == μ[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n]*G[i] / μ[i]) *
                    (1 + ϖ_λ[n]  * Z⁺⁺[i,i,n2] * (dτ_λ[n]  / μ[i]) * wct[i])
            else
                t⁺⁺[i,j,n] = 0.0
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n]  * Z⁺⁺[i,j,n2] * 
                (μ[j] / (μ[i]*G[j] - μ[j]*G[i])) * wct[j] * 
                (exp(-dτ_λ[n] * G[i] / μ[i]) - exp(-dτ_λ[n] * G[j] / μ[j]))
                #(exp(-dτ_λ[n] * G[j] / μ[j]) - exp(-dτ_λ[n] * G[i] / μ[i]))  
        end
    else
        r⁻⁺[i,j,n] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] * G[i] / μ[i]) #Suniti
        else
            t⁺⁺[i,j,n] = 0.0
        end
    end
    nothing
end

@kernel function get_canopy_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, G, Z⁻⁺, Z⁺⁺, μ, ndoubl, wct02, nStokes ,I₀, iμ0, D)
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
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * I₀[ii-i_start+1] 
    end

    if (i>=i_start) && (i<=i_end)
        ctr = i-i_start+1
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        # 1.54 in Fell
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (G[i] * dτ_λ[n] / μ[i]) * exp(-dτ_λ[n] *  G[i] / μ[i])
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        # 1.53 in Fell; 2.14 in Myneni Book 
        J₀⁺[i, 1, n] = 
        wct02 * ϖ_λ[n]  *  Z⁺⁺_I₀ * 
        (μ[i_start] / (μ[i]*G[i_start] - μ[i_start]*G[i])) * 
        (exp(-dτ_λ[n] * G[i] / μ[i]) - exp(-dτ_λ[n] * G[i_start] / μ[i_start]))
        #(exp(-dτ_λ[n] * G[i_start] / μ[i_start]) - exp(-dτ_λ[n] * G[i] / μ[i]))
    end
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
    # 1.52 in Fell
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] *  Z⁻⁺_I₀ * 
            (μ[i_start] / (μ[i]*G[i_start] + μ[i_start]*G[i])) *
            (1 - exp(-dτ_λ[n] * ((G[i] / μ[i]) + (G[i_start] / μ[i_start]))))
             
        #(1 - exp(-(dτ_λ[n] * (G[i_start] * μ[i] + G[i] * μ[i_start]))/(μ[i_start] * μ[i])))
        
    # Multiply with incoming:
    #G is now included in tau_sum already!
    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/μ[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/μ[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end  
    #if (n==840||n==850)    
    #    @show i, n, J₀⁺[i, 1, n], J₀⁻[i, 1, n]      
    #end
    nothing
end