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
                            architecture) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer
    (; qp_μ, iμ₀, wt_μN, qp_μN) = quad_points
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, G) = computed_layer_properties
    
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
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)
 
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

"""
    get_canopy_elem_rt!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, G, Z⁻⁺, Z⁺⁺, μ, wct)

KernelAbstractions elemental R/T kernel for directional canopy scattering.
Each workitem owns one `(i, j, n)` matrix element and evaluates the same
finite-δ single-scattering formulas as the elastic elemental kernel, with
directional path-length factors `G` included in the optical-depth exponents
and stream denominators.
"""
@kernel function get_canopy_elem_rt!(r⁻⁺, t⁺⁺, @Const(ϖ_λ), @Const(dτ_λ),
                                     @Const(G), @Const(Z⁻⁺), @Const(Z⁺⁺),
                                     @Const(μ), @Const(wct))
    FT = eltype(r⁻⁺)
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    if (wct[j] > rt_weight_tol(eltype(wct)))
        # CanopyOptics normalizes Z by ϖ * G(μ_in), so multiply by the
        # incoming-stream G[j] here to recover the physical Γ kernel.

        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * G[j] * Z⁻⁺[i,j,n2] * 
            (μ[j] / (μ[i]*G[j] + μ[j]*G[i])) * wct[j] * 
            -expm1(-dτ_λ[n] * ((G[i] / μ[i]) + (G[j] / μ[j])))
                      
        if (μ[i] == μ[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n]*G[i] / μ[i]) *
                    (1 + ϖ_λ[n] * G[i] * Z⁺⁺[i,i,n2] *
                         (dτ_λ[n]  / μ[i]) * wct[i])
            else
                t⁺⁺[i,j,n] = zero(FT)
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * G[j] * Z⁺⁺[i,j,n2] * 
                (μ[j] / (μ[i]*G[j] - μ[j]*G[i])) * wct[j] * 
                expdiff_neg(dτ_λ[n] * G[i] / μ[i], dτ_λ[n] * G[j] / μ[j])
                #(exp(-dτ_λ[n] * G[j] / μ[j]) - exp(-dτ_λ[n] * G[i] / μ[i]))  
        end
    else
        r⁻⁺[i,j,n] = zero(FT)
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] * G[i] / μ[i]) #Suniti
        else
            t⁺⁺[i,j,n] = zero(FT)
        end
    end
    nothing
end

"""
    get_canopy_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, G, Z⁻⁺, Z⁺⁺, μ,
                            ndoubl, wct02, nStokes, I₀, iμ0, D)

KernelAbstractions source-function kernel for canopy elemental layers. Each
workitem computes the direct-beam `Z * I₀` contractions for one stream and
wavelength, applies the canopy path factor `G` in the finite-δ source
formulas, multiplies by the direct-beam attenuation above the layer, and
applies the upwelling D-matrix sign when required.
"""
@kernel function get_canopy_elem_rt_SFI!(J₀⁺, J₀⁻, @Const(ϖ_λ), @Const(dτ_λ),
                                         @Const(τ_sum), @Const(G), @Const(Z⁻⁺),
                                         @Const(Z⁺⁺), @Const(μ), ndoubl, wct02,
                                         nStokes, @Const(I₀), iμ0, @Const(D))
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    J₀⁺[i, 1, n] = zero(FT)
    J₀⁻[i, 1, n] = zero(FT)
    n2=1
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    Z⁺⁺_I₀ = zero(FT);
    Z⁻⁺_I₀ = zero(FT);
    
    for ii = i_start:i_end
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * I₀[ii-i_start+1]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * I₀[ii-i_start+1] 
    end

    if (i >= i_start) & (i <= i_end)
        ctr = i-i_start+1
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        # 1.54 in Fell
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (G[i] * dτ_λ[n] / μ[i]) * exp(-dτ_λ[n] *  G[i] / μ[i])
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        # 1.53 in Fell; 2.14 in Myneni Book 
        J₀⁺[i, 1, n] = 
        wct02 * ϖ_λ[n] * G[i_start] * Z⁺⁺_I₀ * 
        (μ[i_start] / (μ[i]*G[i_start] - μ[i_start]*G[i])) * 
        expdiff_neg(dτ_λ[n] * G[i] / μ[i], dτ_λ[n] * G[i_start] / μ[i_start])
        #(exp(-dτ_λ[n] * G[i_start] / μ[i_start]) - exp(-dτ_λ[n] * G[i] / μ[i]))
    end
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
    # 1.52 in Fell
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * G[i_start] * Z⁻⁺_I₀ * 
            (μ[i_start] / (μ[i]*G[i_start] + μ[i_start]*G[i])) *
            -expm1(-dτ_λ[n] * ((G[i] / μ[i]) + (G[i_start] / μ[i_start])))
             
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
