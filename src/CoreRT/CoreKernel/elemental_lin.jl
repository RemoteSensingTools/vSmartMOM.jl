#=
 
This file contains RT elemental-related functions
 
=#
"""
    elemental!(pol_type, SFI, τ_sum, τ̇_sum, dτ, F₀, computed_layer_properties,
               m, ndoubl, scatter, quad_points, added_layer, added_layer_lin, architecture)

Compute the elemental (single-scattering) layer reflection and transmission matrices
and their derivatives with respect to the 3 core optical parameters ``(\\tau, \\varpi, \\mathbf{Z})``.

The elemental layer has optical depth ``d\\tau = \\tau / 2^{n_d}`` where ``n_d`` is the number
of doublings needed to build the full layer.

# Single-scattering formulas (Sanghavi & Stephens 2013, Eqs. 19–20)

**Reflection matrix:**
```math
\\mathbf{r}^{-+}(\\mu_i, \\mu_j) = \\varpi \\, \\mathbf{Z}^{-+}(\\mu_i, \\mu_j) 
  \\frac{\\mu_j}{\\mu_i + \\mu_j} \\left(1 - e^{-d\\tau(1/\\mu_i + 1/\\mu_j)}\\right) w_j
```

**Transmission matrix:**
```math
\\mathbf{t}^{++}(\\mu_i, \\mu_j) = \\delta_{ij} e^{-d\\tau/\\mu_i} + 
  \\varpi \\, \\mathbf{Z}^{++}(\\mu_i, \\mu_j) 
  \\frac{\\mu_j}{\\mu_i - \\mu_j} \\left(e^{-d\\tau/\\mu_i} - e^{-d\\tau/\\mu_j}\\right) w_j
```

# Core derivatives (3 per matrix element)
For each matrix ``\\mathbf{M} \\in \\{\\mathbf{r}^{-+}, \\mathbf{t}^{++}, \\ldots\\}``:
- ``\\dot{\\mathbf{M}}[1]``: ``\\partial \\mathbf{M}/\\partial(d\\tau)`` — optical depth derivative
- ``\\dot{\\mathbf{M}}[2]``: ``\\partial \\mathbf{M}/\\partial\\varpi`` — single-scattering albedo derivative
- ``\\dot{\\mathbf{M}}[3]``: ``\\partial \\mathbf{M}/\\partial\\mathbf{Z}`` — phase matrix derivative

When `SFI=true`, the source function vectors ``\\mathbf{j}_0^+, \\mathbf{j}_0^-`` and their
derivatives are also computed for the solar beam contribution.

# Arguments
- `pol_type`: Polarization type (I, IQU, or IQUV).
- `SFI::Bool`: Whether to compute Source Function Integration terms.
- `τ_sum`: Cumulative optical depth above this layer `[nSpec]`.
- `τ̇_sum`: Derivative of cumulative τ w.r.t. parameters `[Nparams × nSpec]`.
- `dτ`: Elemental optical depth ``\\tau/2^{n_d}`` `[nSpec]`.
- `F₀`: Solar irradiance Stokes vector `[nStokes × nSpec]`.
- `computed_layer_properties`: Forward optical properties ``(\\tau, \\varpi, \\mathbf{Z}^{++}, \\mathbf{Z}^{-+})``.
- `m::Int`: Fourier moment index.
- `ndoubl::Int`: Number of doublings.
- `scatter::Bool`: Whether the layer scatters.
- `quad_points`: Quadrature points and weights.
- `added_layer`: Output: forward RT matrices (modified in-place).
- `added_layer_lin`: Output: linearized RT matrices (modified in-place).
- `architecture`: CPU or GPU.
"""
function elemental!(pol_type, SFI::Bool,
                τ_sum::AbstractArray,#{FT2,1}, #Suniti
                τ̇_sum::AbstractArray,
                dτ::AbstractArray,
                F₀::AbstractArray,#{FT,2},    # Stokes vector of solar/stellar irradiance
                computed_layer_properties,
                computed_layer_properties_lin,
                m::Int,                     # m: fourier moment
                ndoubl::Int,                # ndoubl: number of doubling computations needed 
                scatter::Bool,              # scatter: flag indicating scattering
                quad_points::QuadPoints{FT}, # struct with quadrature points, weights, 
                added_layer::AddedLayer{FT}, 
                added_layer_lin::AddedLayerLin{FT}, 
                architecture) where {FT<:AbstractFloat}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer
    (; ṙ⁺⁻, ṙ⁻⁺, ṫ⁻⁻, ṫ⁺⁺, J̇₀⁺, J̇₀⁻) = added_layer_lin
    (; qp_μ, iμ₀, wt_μN, qp_μN) = quad_points
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    (; τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺) = computed_layer_properties_lin

    arr_type = array_type(architecture)
    qp_μN = arr_type(qp_μN)
    wt_μN = arr_type(wt_μN)
    τ_sum = arr_type(τ_sum)
    τ̇_sum = arr_type(τ̇_sum)
    I₀    = arr_type(pol_type.I₀)
    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)

    # Chain-rule inputs (convert to device arrays)
    nparams = size(τ̇, 2)   # τ̇ is [nSpec, Nparams] — Nparams in last dim
    dτ̇_dev = arr_type(τ̇ ./ FT(2^ndoubl))   # elemental τ̇
    ϖ̇_dev  = arr_type(ϖ̇)
    Ż⁻⁺_dev = arr_type(Ż⁻⁺)
    Ż⁺⁺_dev = arr_type(Ż⁺⁺)

    # If in scattering mode:
    if scatter
   
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)
        # Zero forward, 3-core, AND ap_ arrays
        r⁻⁺ .= zero(FT)
        t⁺⁺ .= zero(FT)
        ṙ⁻⁺ .= zero(FT)
        ṫ⁺⁺ .= zero(FT)
        j₀⁺ .= zero(FT)
        j₀⁻ .= zero(FT)
        J̇₀⁺ .= zero(FT)
        J̇₀⁻ .= zero(FT)
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)
        added_layer_lin.ap_ṙ⁻⁺ .= zero(FT)
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        added_layer_lin.ap_J̇₀⁻ .= zero(FT)

        # Fused elemental + chain rule kernel
        kernel! = get_elem_rt_fused!(device)
        event = kernel!(r⁻⁺, t⁺⁺,
                    ṙ⁻⁺, ṫ⁺⁺,
                    added_layer_lin.ap_ṙ⁻⁺, added_layer_lin.ap_ṫ⁺⁺,
                    added_layer_lin.ap_ṙ⁺⁻, added_layer_lin.ap_ṫ⁻⁻,
                    ϖ, dτ, Z⁻⁺, Z⁺⁺,
                    dτ̇_dev, ϖ̇_dev, Ż⁻⁺_dev, Ż⁺⁺_dev,
                    qp_μN, wct2,
                    nparams, ndoubl, pol_type.n,
                    ndrange=size(r⁻⁺))
        synchronize_if_gpu()

        if SFI
            # Fused SFI + chain rule + Bug 22 fix kernel
            kernel! = get_elem_rt_SFI_fused!(device)
            event = kernel!(j₀⁺, j₀⁻,
                J̇₀⁺, J̇₀⁻,
                added_layer_lin.ap_J̇₀⁺, added_layer_lin.ap_J̇₀⁻,
                ϖ, dτ,
                τ_sum, τ̇_sum,
                Z⁻⁺, Z⁺⁺,
                arr_type(F₀),
                dτ̇_dev, ϖ̇_dev, Ż⁻⁺_dev, Ż⁺⁺_dev,
                qp_μN, ndoubl, wct02,
                pol_type.n, I₀, iμ₀, D, nparams,
                ndrange=size(j₀⁺))
        end
        synchronize_if_gpu()

        # Apply D Matrix to forward quantities (fused kernel handles derivative D internally)
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        # SFI D-matrix already applied inside fused kernel
        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, j₀⁻)
        end
    else
        # No scattering: zero ap_ arrays, set transmission only
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)
        added_layer_lin.ap_ṙ⁻⁺ .= zero(FT)
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        added_layer_lin.ap_J̇₀⁻ .= zero(FT)

        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
        ṫ⁺⁺[:, :, :, 1] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}
        ṫ⁻⁻[:, :, :, 1] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}

        # Chain rule for no-scatter: ap_ṫ = ṫ[1]*dτ̇ (only τ derivative matters)
        nspec_here = size(τ, 1)
        for iparam = 1:nparams
            for iλ = 1:nspec_here
                @views added_layer_lin.ap_ṫ⁺⁺[:,:,iλ,iparam] .= ṫ⁺⁺[:,:,iλ,1] .* dτ̇_dev[iλ,iparam]
                @views added_layer_lin.ap_ṫ⁻⁻[:,:,iλ,iparam] .= ṫ⁻⁻[:,:,iλ,1] .* dτ̇_dev[iλ,iparam]
            end
        end
    end
end

@kernel function get_elem_rt!(r⁻⁺, t⁺⁺,
                        ṙ⁻⁺, ṫ⁺⁺, 
                        ϖ_λ, dτ_λ, 
                        Z⁻⁺, Z⁺⁺, 
                        qp_μN, wct) 
    FT = eltype(r⁻⁺)
    n2 = 1
    i, j, n = @index(Global, NTuple) 
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    if (wct[j] > eps(FT)) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dτ = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(1/μᵢ) ̇exp{-τ ̇(1/μᵢ + 1/μⱼ)}  ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dϖ = 𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        # d𝐑⁻⁺(μᵢ, μⱼ)/dZ = ϖ ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            #Z⁻⁺[i,j] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))) 
        # derivative wrt τ_λ
        ṙ⁻⁺[i,j,n,1] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            (1/qp_μN[i]) * wct[j] * 
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))) 
        # derivative wrt ϖ
        ṙ⁻⁺[i,j,n,2] = ϖ_λ[n] == 0 ? FT(0) : r⁻⁺[i, j, n] / ϖ_λ[n]
        # derivative wrt Z
        # derivative wrt Z: direct formula avoids 0/0 when Z=0
        ṙ⁻⁺[i,j,n,3] = ϖ_λ[n] * 
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))))
                    
        if (qp_μN[i] == qp_μN[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ}(1 + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ))) ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dτ_λ = (exp{-τ/μⱼ}/μᵢ)⋅(ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)⋅(1-τ/μⱼ)-1) ̇𝑤ⱼ  
            # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = 
                    exp(-dτ_λ[n] / qp_μN[i]) * (1 / qp_μN[i]) *
                    (-1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (1 - dτ_λ[n] / qp_μN[i]))
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i]    
                # derivative wrt Z
                ṫ⁺⁺[i,j,n,3] = 
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[i]
            else
                # 𝐓⁺⁺(μᵢ, μⱼ) = (exp{-τ/μⱼ}(ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(τ/μᵢ))) ̇𝑤ⱼ        
                # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = (exp{-τ/μⱼ}⋅ϖ ̇𝐙⁺⁺(μᵢ, μᵢ)/μᵢ)⋅(1 - τ/μⱼ) ̇𝑤ⱼ
                # d𝐓⁺⁺(μᵢ, μᵢ)/dϖ_λ = 𝐙⁺⁺(μᵢ, μᵢ)⋅(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
                # d𝐓⁺⁺(μᵢ, μᵢ)/dZ   = ϖ ̇(τ/μᵢ) ̇exp{-τ/μᵢ} ̇𝑤ᵢ
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
                # derivative wrt τ_λ
                ṫ⁺⁺[i,j,n,1] = (exp(-dτ_λ[n] / qp_μN[j]) *
                        ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) * 
                        (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
                # derivative wrt ϖ_λ
                ṫ⁺⁺[i,j,n,2] = ϖ_λ[n] == 0 ? FT(0) : t⁺⁺[i, j, n] / ϖ_λ[n]
                # derivative wrt Z
                # derivative wrt Z: direct formula avoids 0/0
                ṫ⁺⁺[i,j,n,3] = exp(-dτ_λ[n] / qp_μN[j]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[j]
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # d𝐓⁺⁺(μᵢ, μⱼ)/dτ_λ = -ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ}/μᵢ - exp{-τ/μⱼ}/μⱼ) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                #Z⁺⁺[i,j] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j])) 
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = -ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i])/ qp_μN[i] - 
                exp(-dτ_λ[n] / qp_μN[j])/ qp_μN[j]) 
            # derivative wrt ϖ_λ
            ṫ⁺⁺[i,j,n,2] = ϖ_λ[n] == 0 ? FT(0) : t⁺⁺[i, j, n] / ϖ_λ[n]
            # derivative wrt Z
            # derivative wrt Z: direct formula avoids 0/0
            ṫ⁺⁺[i,j,n,3] = ϖ_λ[n] * 
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] * 
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j]))
        end
    else
        #r⁻⁺[i,j,n] = 0.0
        #ṙ⁻⁺[i,j,n,:] = 0.0
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i]) #Suniti
            # derivative wrt τ_λ
            ṫ⁺⁺[i,j,n,1] = -exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i]
        #else
        #    t⁺⁺[i,j,n] = 0.0
            # derivative wrt τ_λ
        #    ṫ⁺⁺[i,j,n,1] = 0.0
        end
        # derivative wrt ϖ_λ
        #ṫ⁺⁺[i,j,n,2] = 0.0
        # derivative wrt Z
        #ṫ⁺⁺[i,j,n,3] = 0.0
    end
    nothing
end

@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, 
                J̇₀⁺, J̇₀⁻, 
                ϖ_λ, dτ_λ, 
                τ_sum, τ̇_sum, 
                Z⁻⁺, Z⁺⁺, F₀,
                qp_μN, ndoubl, wct02, nStokes,
                I₀, iμ0, D)
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    #J₀⁺[i, 1, n]=0
    #J₀⁻[i, 1, n]=0
    #J̇₀⁺[i, 1, n, 1:3]=0
    #J̇₀⁻[i, 1, n, 1:3]=0
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
        # derivative wrt τ
        J̇₀⁺[i, 1, n, 1] = J₀⁺[i, 1, n]*(1/dτ_λ[n] - 1/qp_μN[i])
        # derivative wrt ϖ
        J̇₀⁺[i, 1, n, 2] = ϖ_λ[n] == 0 ? FT(0) : J₀⁺[i, 1, n] / ϖ_λ[n]
        # derivative wrt Z (safe division: 0/0 → 0)
        J̇₀⁺[i, 1, n, 3] = Z⁺⁺_I₀ == 0 ? FT(0) : J₀⁺[i, 1, n] / Z⁺⁺_I₀
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * 
            (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
        # derivative wrt τ
        J̇₀⁺[i, 1, n, 1] = - wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * 
            (exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i] - exp(-dτ_λ[n] / qp_μN[i_start]) / qp_μN[i_start])
        # derivative wrt ϖ
        J̇₀⁺[i, 1, n, 2] = ϖ_λ[n] == 0 ? FT(0) : J₀⁺[i, 1, n] / ϖ_λ[n]
        # derivative wrt Z (safe division: 0/0 → 0)
        J̇₀⁺[i, 1, n, 3] = Z⁺⁺_I₀ == 0 ? FT(0) : J₀⁺[i, 1, n] / Z⁺⁺_I₀
    end
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * 
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))
    # derivative wrt τ
    J̇₀⁻[i, 1, n, 1] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) * 
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))) *
            ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))
    # derivative wrt ϖ
    J̇₀⁻[i, 1, n, 2] = ϖ_λ[n] == 0 ? FT(0) : J₀⁻[i, 1, n] / ϖ_λ[n]
    # derivative wrt Z (safe division: 0/0 → 0)
    J̇₀⁻[i, 1, n, 3] = Z⁻⁺_I₀ == 0 ? FT(0) : J₀⁻[i, 1, n] / Z⁻⁺_I₀

    # TODO: Move this out until after doubling (it is not necessary to consider this here already if Raman scattering is not involved)
    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/qp_μN[i_start])

    # Bug 22 fix: Remove τ̇_sum[1,n] contribution from core derivative.
    # The τ̇_sum beam attenuation derivative is per-physical-parameter and must be
    # added AFTER the chain rule (in rt_kernel!), not here in the 3-core framework.
    # Old code used τ̇_sum[1,n] which only captured parameter 1's contribution.
    J̇₀⁺[i, 1, n, 1] = J̇₀⁺[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 1] = J̇₀⁻[i, 1, n, 1]*exp(-τ_sum[n]/qp_μN[i_start])
    J̇₀⁺[i, 1, n, 2] = J̇₀⁺[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 2] = J̇₀⁻[i, 1, n, 2]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁺[i, 1, n, 3] = J̇₀⁺[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁺[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])
    J̇₀⁻[i, 1, n, 3] = J̇₀⁻[i, 1, n, 3]*exp(-τ_sum[n]/qp_μN[i_start]) #+
                        #J₀⁻[i, 1, n] * (-τ̇_sum[1,n]/qp_μN[i_start])


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

# ============================================================================
# Fused kernels: combine elemental RT + chain rule in a single pass.
# These eliminate the separate lin_added_layer_all_params! call by computing
# ap_ṙ⁻⁺, ap_ṫ⁺⁺ (and optionally ap_ṙ⁺⁻, ap_ṫ⁻⁻) directly from local
# 3-core scalar intermediates, avoiding ~12 full-array reads.
# ============================================================================

"""
    get_elem_rt_fused!(...)

Fused elemental R/T kernel: computes forward r⁻⁺, t⁺⁺ and their per-parameter
derivatives ap_ṙ⁻⁺, ap_ṫ⁺⁺ (and ap_ṙ⁺⁻, ap_ṫ⁻⁻ for ndoubl < 1) in a single pass.

The 3-core derivatives (ṙ⁻⁺[1:3], ṫ⁺⁺[1:3]) are kept as local scalars and used
directly for the chain rule, then also written to their arrays for backward
compatibility with the 3-core doubling path.
"""
@kernel function get_elem_rt_fused!(r⁻⁺, t⁺⁺,
                        ṙ⁻⁺, ṫ⁺⁺,
                        ap_ṙ⁻⁺, ap_ṫ⁺⁺, ap_ṙ⁺⁻, ap_ṫ⁻⁻,
                        ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺,
                        dτ̇, ϖ̇, Ż⁻⁺, Ż⁺⁺_lin,
                        qp_μN, wct,
                        nparams, ndoubl, pol_n)
    FT = eltype(r⁻⁺)
    i, j, n = @index(Global, NTuple)
    n2 = 1
    if size(Z⁻⁺, 3) > 1
        n2 = n
    end
    n2_lin = 1
    if size(Ż⁻⁺, 3) > 1
        n2_lin = n
    end

    # D-matrix Stokes signs
    i_stokes = mod(i, pol_n)
    j_stokes = mod(j, pol_n)
    i12 = (pol_n == 1) | ((1 <= i_stokes) & (i_stokes <= 2))
    j12 = (pol_n == 1) | ((1 <= j_stokes) & (j_stokes <= 2))
    same_block = (i12 & j12) | (!i12 & !j12)
    d_sign = ifelse(same_block, one(FT), -one(FT))
    di = ifelse(i12, one(FT), -one(FT))
    dj = ifelse(j12, one(FT), -one(FT))

    # R⁻⁺ row-sign correction for elemental D-matrix (ndoubl >= 1 negates Stokes 3,4 rows)
    sign_r = ifelse((ndoubl >= 1) & !i12, -one(FT), one(FT))

    # Local 3-core derivative scalars
    ṙ_tau = FT(0); ṙ_w = FT(0); ṙ_Z = FT(0)
    ṫ_tau = FT(0); ṫ_w = FT(0); ṫ_Z = FT(0)

    if (wct[j] > eps(FT))
        # ---- R⁻⁺(μᵢ, μⱼ) ----
        r⁻⁺[i,j,n] =
            ϖ_λ[n] * Z⁻⁺[i,j,n2] *
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] *
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))))

        ṙ_tau = ϖ_λ[n] * Z⁻⁺[i,j,n2] *
            (1/qp_μN[i]) * wct[j] *
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j])))
        ṙ_w = ϖ_λ[n] == 0 ? FT(0) : r⁻⁺[i,j,n] / ϖ_λ[n]
        ṙ_Z = ϖ_λ[n] *
            (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] *
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[j]))))

        # Write 3-core (backward compat)
        ṙ⁻⁺[i,j,n,1] = ṙ_tau
        ṙ⁻⁺[i,j,n,2] = ṙ_w
        ṙ⁻⁺[i,j,n,3] = ṙ_Z

        # ---- T⁺⁺(μᵢ, μⱼ) ----
        if (qp_μN[i] == qp_μN[j])
            if i == j
                t⁺⁺[i,j,n] =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i])
                ṫ_tau =
                    exp(-dτ_λ[n] / qp_μN[i]) * (1 / qp_μN[i]) *
                    (-1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (1 - dτ_λ[n] / qp_μN[i]))
                ṫ_w =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    Z⁺⁺[i,i,n2] * (dτ_λ[n] / qp_μN[i]) * wct[i]
                ṫ_Z =
                    exp(-dτ_λ[n] / qp_μN[i]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[i]
            else
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / qp_μN[i]) * wct[j])
                ṫ_tau = (exp(-dτ_λ[n] / qp_μN[j]) *
                        ϖ_λ[n] * Z⁺⁺[i,j,n2] / qp_μN[i]) *
                        (1 - dτ_λ[n] / qp_μN[j]) * wct[j]
                ṫ_w = ϖ_λ[n] == 0 ? FT(0) : t⁺⁺[i,j,n] / ϖ_λ[n]
                ṫ_Z = exp(-dτ_λ[n] / qp_μN[j]) *
                    ϖ_λ[n] * (dτ_λ[n] / qp_μN[i]) * wct[j]
            end
        else
            t⁺⁺[i,j,n] =
                ϖ_λ[n] * Z⁺⁺[i,j,n2] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j]))
            ṫ_tau = -ϖ_λ[n] * Z⁺⁺[i,j,n2] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                (exp(-dτ_λ[n] / qp_μN[i])/ qp_μN[i] -
                exp(-dτ_λ[n] / qp_μN[j])/ qp_μN[j])
            ṫ_w = ϖ_λ[n] == 0 ? FT(0) : t⁺⁺[i,j,n] / ϖ_λ[n]
            ṫ_Z = ϖ_λ[n] *
                (qp_μN[j] / (qp_μN[i] - qp_μN[j])) * wct[j] *
                (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[j]))
        end

        # Write 3-core (backward compat)
        ṫ⁺⁺[i,j,n,1] = ṫ_tau
        ṫ⁺⁺[i,j,n,2] = ṫ_w
        ṫ⁺⁺[i,j,n,3] = ṫ_Z

        # ---- Fused chain rule: ap_ = ṙ_tau*dτ̇ + ṙ_w*ϖ̇ + ṙ_Z*Ż ----
        for iparam = 1:nparams
            val_r = ṙ_tau * dτ̇[n,iparam] + ṙ_w * ϖ̇[n,iparam] + ṙ_Z * Ż⁻⁺[i,j,n2_lin,iparam]
            val_t = ṫ_tau * dτ̇[n,iparam] + ṫ_w * ϖ̇[n,iparam] + ṫ_Z * Ż⁺⁺_lin[i,j,n2_lin,iparam]

            ap_ṙ⁻⁺[i,j,n,iparam] = sign_r * val_r
            ap_ṫ⁺⁺[i,j,n,iparam] = val_t

            if ndoubl < 1
                # For ndoubl < 1 (no doubling): compute ⁺⁻ and ⁻⁻ via D-matrix
                # ṙ⁺⁻ = d_sign * (ṙ_tau*dτ̇ + ṙ_w*ϖ̇ + ṙ_Z * D·Ż⁻⁺·D)
                # where D·Ż·D at (i,j) = di*dj*Ż
                ap_ṙ⁺⁻[i,j,n,iparam] = d_sign * (ṙ_tau * dτ̇[n,iparam] + ṙ_w * ϖ̇[n,iparam] + ṙ_Z * di * dj * Ż⁻⁺[i,j,n2_lin,iparam])
                ap_ṫ⁻⁻[i,j,n,iparam] = d_sign * (ṫ_tau * dτ̇[n,iparam] + ṫ_w * ϖ̇[n,iparam] + ṫ_Z * di * dj * Ż⁺⁺_lin[i,j,n2_lin,iparam])
            end
        end
    else
        # No scattering weight: only diagonal transmission
        if i == j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / qp_μN[i])
            ṫ_tau = -exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i]
            ṫ⁺⁺[i,j,n,1] = ṫ_tau
            for iparam = 1:nparams
                val_t = ṫ_tau * dτ̇[n,iparam]
                ap_ṫ⁺⁺[i,j,n,iparam] = val_t
                if ndoubl < 1
                    # diagonal: same_block=true so d_sign=1
                    ap_ṫ⁻⁻[i,j,n,iparam] = val_t
                end
            end
        end
    end
    nothing
end

"""
    get_elem_rt_SFI_fused!(...)

Fused SFI source kernel: computes J₀⁺, J₀⁻ and their per-parameter derivatives
ap_J̇₀⁺, ap_J̇₀⁻ in a single pass, including the Bug 22 beam attenuation fix.

Eliminates the separate chain-rule pass for SFI terms and the per-parameter
τ̇_sum correction loop in rt_kernel!.
"""
@kernel function get_elem_rt_SFI_fused!(J₀⁺, J₀⁻,
                J̇₀⁺, J̇₀⁻,
                ap_J̇₀⁺, ap_J̇₀⁻,
                ϖ_λ, dτ_λ,
                τ_sum, τ̇_sum,
                Z⁻⁺, Z⁺⁺, F₀,
                dτ̇, ϖ̇, Ż⁻⁺, Ż⁺⁺_lin,
                qp_μN, ndoubl, wct02, nStokes,
                I₀, iμ0, D, nparams)
    i_start  = nStokes*(iμ0-1) + 1
    i_end    = nStokes*iμ0

    i, _, n = @index(Global, NTuple)
    FT = eltype(I₀)

    n2 = 1
    if size(Z⁻⁺, 3) > 1
        n2 = n
    end
    n2_lin = 1
    if size(Ż⁻⁺, 3) > 1
        n2_lin = n
    end

    # Forward Z·I₀ products
    Z⁺⁺_I₀ = FT(0.0)
    Z⁻⁺_I₀ = FT(0.0)
    for ii = i_start:i_end
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n2]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n2]
    end

    # ---- J₀⁺ and 3-core scalars ----
    J̇⁺_tau = FT(0); J̇⁺_w = FT(0); J̇⁺_Z = FT(0)

    if (i>=i_start) && (i<=i_end)
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (dτ_λ[n] / qp_μN[i]) * exp(-dτ_λ[n] / qp_μN[i])
        J̇⁺_tau = J₀⁺[i, 1, n]*(1/dτ_λ[n] - 1/qp_μN[i])
        J̇⁺_w = ϖ_λ[n] == 0 ? FT(0) : J₀⁺[i, 1, n] / ϖ_λ[n]
        J̇⁺_Z = Z⁺⁺_I₀ == 0 ? FT(0) : J₀⁺[i, 1, n] / Z⁺⁺_I₀
    else
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ *
            (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) * (exp(-dτ_λ[n] / qp_μN[i]) - exp(-dτ_λ[n] / qp_μN[i_start]))
        J̇⁺_tau = - wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] - qp_μN[i_start])) *
            (exp(-dτ_λ[n] / qp_μN[i]) / qp_μN[i] - exp(-dτ_λ[n] / qp_μN[i_start]) / qp_μN[i_start])
        J̇⁺_w = ϖ_λ[n] == 0 ? FT(0) : J₀⁺[i, 1, n] / ϖ_λ[n]
        J̇⁺_Z = Z⁺⁺_I₀ == 0 ? FT(0) : J₀⁺[i, 1, n] / Z⁺⁺_I₀
    end

    # ---- J₀⁻ and 3-core scalars ----
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) *
            (1 - exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))))
    J̇⁻_tau = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (qp_μN[i_start] / (qp_μN[i] + qp_μN[i_start])) *
            exp(-dτ_λ[n] * ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))) *
            ((1 / qp_μN[i]) + (1 / qp_μN[i_start]))
    J̇⁻_w = ϖ_λ[n] == 0 ? FT(0) : J₀⁻[i, 1, n] / ϖ_λ[n]
    J̇⁻_Z = Z⁻⁺_I₀ == 0 ? FT(0) : J₀⁻[i, 1, n] / Z⁻⁺_I₀

    # ---- Apply beam attenuation exp(-τ_sum/μ₀) ----
    beam_atten = exp(-τ_sum[n]/qp_μN[i_start])
    J₀⁺[i, 1, n] *= beam_atten
    J₀⁻[i, 1, n] *= beam_atten
    J̇⁺_tau *= beam_atten
    J̇⁺_w   *= beam_atten
    J̇⁺_Z   *= beam_atten
    J̇⁻_tau *= beam_atten
    J̇⁻_w   *= beam_atten
    J̇⁻_Z   *= beam_atten

    # Write 3-core arrays (backward compat)
    J̇₀⁺[i, 1, n, 1] = J̇⁺_tau
    J̇₀⁺[i, 1, n, 2] = J̇⁺_w
    J̇₀⁺[i, 1, n, 3] = J̇⁺_Z
    J̇₀⁻[i, 1, n, 1] = J̇⁻_tau
    J̇₀⁻[i, 1, n, 2] = J̇⁻_w
    J̇₀⁻[i, 1, n, 3] = J̇⁻_Z

    # ---- D-matrix for J₀⁻ (ndoubl >= 1) ----
    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n]
        J̇⁻_tau = D[i,i]*J̇⁻_tau
        J̇⁻_w   = D[i,i]*J̇⁻_w
        J̇⁻_Z   = D[i,i]*J̇⁻_Z
        J̇₀⁻[i, 1, n, 1] = J̇⁻_tau
        J̇₀⁻[i, 1, n, 2] = J̇⁻_w
        J̇₀⁻[i, 1, n, 3] = J̇⁻_Z
    end

    # ---- Fused chain rule + Bug 22 fix ----
    for iparam = 1:nparams
        # Compute Ż·I₀ dot products for this parameter
        Ż⁺⁺_I₀_p = FT(0)
        Ż⁻⁺_I₀_p = FT(0)
        for ii = i_start:i_end
            Ż⁺⁺_I₀_p += Ż⁺⁺_lin[i, ii, n2_lin, iparam] * F₀[ii-i_start+1, n2]
            Ż⁻⁺_I₀_p += Ż⁻⁺[i, ii, n2_lin, iparam] * F₀[ii-i_start+1, n2]
        end

        # Chain rule: ap_J̇ = J̇_tau*dτ̇ + J̇_w*ϖ̇ + J̇_Z*Ż_I₀
        ap_J̇₀⁺[i, 1, n, iparam] = J̇⁺_tau * dτ̇[n,iparam] + J̇⁺_w * ϖ̇[n,iparam] + J̇⁺_Z * Ż⁺⁺_I₀_p
        ap_J̇₀⁻[i, 1, n, iparam] = J̇⁻_tau * dτ̇[n,iparam] + J̇⁻_w * ϖ̇[n,iparam] + J̇⁻_Z * Ż⁻⁺_I₀_p

        # Bug 22 fix: per-parameter τ̇_sum beam attenuation derivative
        # d(exp(-τ_sum/μ₀))/dp_j * J₀ = -τ̇_sum[j]/μ₀ * J₀
        ap_J̇₀⁺[i, 1, n, iparam] += J₀⁺[i, 1, n] * (-τ̇_sum[n, iparam] / qp_μN[i_start])
        ap_J̇₀⁻[i, 1, n, iparam] += J₀⁻[i, 1, n] * (-τ̇_sum[n, iparam] / qp_μN[i_start])
    end

    nothing
end

@kernel function apply_D_elemental!(ndoubl, pol_n, 
                                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                                ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)
    i, j, n = @index(Global, NTuple) #how best to do this for linearization? Is : okay, or should I use an iparam index?

    if ndoubl < 1
        ii = mod1(i, pol_n)
        jj = mod1(j, pol_n)
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
            ṙ⁺⁻[i,j,n,1] = ṙ⁻⁺[i,j,n,1]
            ṙ⁺⁻[i,j,n,2] = ṙ⁻⁺[i,j,n,2]
            ṙ⁺⁻[i,j,n,3] = ṙ⁻⁺[i,j,n,3]
            ṫ⁻⁻[i,j,n,1] = ṫ⁺⁺[i,j,n,1]
            ṫ⁻⁻[i,j,n,2] = ṫ⁺⁺[i,j,n,2]
            ṫ⁻⁻[i,j,n,3] = ṫ⁺⁺[i,j,n,3]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] 
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] 
            ṙ⁺⁻[i,j,n,1] = -ṙ⁻⁺[i,j,n,1] 
            ṙ⁺⁻[i,j,n,2] = -ṙ⁻⁺[i,j,n,2] 
            ṙ⁺⁻[i,j,n,3] = -ṙ⁻⁺[i,j,n,3] 
            ṫ⁻⁻[i,j,n,1] = -ṫ⁺⁺[i,j,n,1] 
            ṫ⁻⁻[i,j,n,2] = -ṫ⁺⁺[i,j,n,2] 
            ṫ⁻⁻[i,j,n,3] = -ṫ⁺⁺[i,j,n,3] 
        end
    else
        if mod1(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
            ṙ⁻⁺[i,j,n,1] = - ṙ⁻⁺[i,j,n,1]
            ṙ⁻⁺[i,j,n,2] = - ṙ⁻⁺[i,j,n,2]
            ṙ⁻⁺[i,j,n,3] = - ṙ⁻⁺[i,j,n,3]
        end 
    end
    nothing
end

@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻, J̇₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if mod1(i, pol_n) > 2
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
                    ṙ⁻⁺::AbstractArray{FT,4}, 
                    ṫ⁺⁺::AbstractArray{FT,4}, 
                    ṙ⁺⁻::AbstractArray{FT,4}, 
                    ṫ⁻⁻::AbstractArray{FT,4}) where {FT}
    device = devi(architecture(r⁻⁺))
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, 
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, 
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻, 
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
