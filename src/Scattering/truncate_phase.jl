#=

This file specifies how to truncate the AerosolOptics struct, given the
truncation type. All `truncate_phase` methods share the contract:

    truncate_phase(method::AbstractTruncationType, aero::AerosolOptics; kwargs...)
        -> AerosolOptics

The returned `AerosolOptics` carries the truncated Greek coefficients
and the `fᵗ = 1 - c₀` retained-fraction parameter; downstream pipeline
code applies the τ/ω rescaling per Sanghavi & Stephens 2015 Eq. 8.

=#

"""
    truncate_phase(::NoTruncation, aero::AerosolOptics; kwargs...) -> AerosolOptics

Identity passthrough. Returns an `AerosolOptics` with the same Greek
coefficients, ω̃ and k, and `fᵗ = 0` (the `f_tr → 0` limit of Sanghavi
& Stephens 2015 Eq. 8 — the truncation-modified `τ*`, `ω*`, `Z*`
collapse to the originals).

Note: raw Mie outputs initialise `fᵗ = 1` as a "untruncated yet"
sentinel — passing them through unchanged would let downstream
`delta_m_forward` interpret the 1 as "everything is in the forward
peak" and silently zero out aerosol scattering. We return `fᵗ = 0`
so the rescaling is a true no-op.
"""
truncate_phase(::NoTruncation, aero::AerosolOptics{FT}; kwargs...) where {FT} =
    AerosolOptics(greek_coefs = aero.greek_coefs, ω̃ = aero.ω̃,
                  k = aero.k, fᵗ = zero(FT), derivs = aero.derivs)

"""
    truncate_phase_lowconf(::NoTruncation, aero::AerosolOptics; kwargs...) -> AerosolOptics

Identity passthrough; matches [`truncate_phase`](@ref) for `NoTruncation`.
Same `fᵗ = 0` reset.
"""
truncate_phase_lowconf(::NoTruncation, aero::AerosolOptics; kwargs...) =
    truncate_phase(NoTruncation(), aero; kwargs...)


@doc raw"""
    truncate_phase_lowconf(mod::δBGE, aero::AerosolOptics; reportFit=false) -> AerosolOptics

Legacy/low-confidence δ-BGE truncation variant.

Fits truncated coefficients outside the forward exclusion cone (`Δ_angle`) and
rescales by retained scattering fraction ``c_0``. The returned truncation factor
is:

```math
f^t = 1 - c_0.
```
"""
function truncate_phase_lowconf(mod::δBGE, aero::AerosolOptics{FT}; reportFit=false) where {FT}
    (; greek_coefs, ω̃, k) = aero
    (; α, β, γ, δ, ϵ, ζ) = greek_coefs
    (; l_max, Δ_angle) = mod


    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre(length(β));

    # Reconstruct phase matrix elements:
    scattering_matrix, P, P² = reconstruct_phase(greek_coefs, μ; returnLeg=true)

    (; f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄) = scattering_matrix

    # Find elements that exclude the peak (if wanted!)
    iμ = findall(x -> x < cosd(Δ_angle), μ)

    # Prefactor for P2:
    fac = zeros(FT,l_max);
    for l = 2:l_max - 1
        fac[l + 1] = sqrt(FT(1) / ( ( l - FT(1)) * l * (l + FT(1)) * (l + FT(2)) ));
    end

    # Create subsets (for Ax=y weighted least-squares fits):
    y₁₁ = view(f₁₁, iμ)
    y₁₂ = view(f₁₂, iμ)
    y₃₄ = view(f₃₄, iμ)
    A   = view(P, iμ, 1:l_max)
    B   = fac[3:end]' .* view(P², iμ, 3:l_max)

    # Weights (also avoid division by 0)
    minY = zeros(length(iμ)) .+ FT(1e-8);
    W₁₁ = Diagonal(w_μ[iμ] ./ max(abs.(y₁₁), minY));
    W₁₂ = Diagonal(w_μ[iμ] ./ max(abs.(y₁₂), minY));
    W₃₄ = Diagonal(w_μ[iμ] ./ max(abs.(y₃₄), minY));
    # W₁₂ = Diagonal(w_μ[iμ]);
    # W₃₄ = Diagonal(w_μ[iμ]);
    # Julia backslash operator for least squares (just like Matlab);
    cl = ((W₁₁ * A) \ (W₁₁ * y₁₁))   # B in δ-BGR (β)
    γᵗ = similar(cl); γᵗ[1:2] .=0
    ϵᵗ = similar(cl); ϵᵗ[1:2] .=0
    γᵗ[3:end] = ((W₁₂ * B) \ (W₁₂ * y₁₂))   # G in δ-BGE (γ)
    ϵᵗ[3:end] = ((W₃₄ * B) \ (W₃₄ * y₃₄))   # E in δ-BGE (ϵ)
    
    if reportFit
        println("Errors in δ-BGE fits:")
        mod_y = convert.(FT, A * cl)
        mod_γ = convert.(FT, B * γᵗ[3:end])
        mod_ϵ = convert.(FT, B * ϵᵗ[3:end])
        @show StatsBase.rmsd(mod_y, y₁₁; normalize=true)
        @show StatsBase.rmsd(mod_γ, y₁₂; normalize=true)
        @show StatsBase.rmsd(mod_ϵ, y₃₄; normalize=true)
    end

    # Integrate truncated function for later renormalization (here: fraction that IS still scattered):
    c₀ = FT(cl[1]) # ( w_μ' * (P[:,1:l_max] * cl) ) / 2
    
    # Compute truncated greek coefficients:
    βᵗ = cl / c₀                                    # Eq. 38a, B in δ-BGR (β)
    δᵗ = (δ[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38b, derived from β
    αᵗ = (α[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38c, derived from β
    ζᵗ = (ζ[1:l_max] .- (β[1:l_max] .- cl)) / c₀    # Eq. 38d, derived from β

    # Truncated Greek coefficients only — ω̃ and k pass through. The
    # τ / ω rescaling per Sanghavi & Stephens 2015 Eq. 8 is applied
    # later in the pipeline by `delta_m_forward` (see
    # CoreRT/LayerOpticalProperties/delta_m_truncation.jl): given
    # `(τ, ω̃, fᵗ)` it returns `(τ_mod, ϖ_mod)` with the proper
    # `(1 − fᵗ·ω̃)` and `(1−fᵗ)·ω̃/(1−fᵗ·ω̃)` factors. Re-applying
    # them here would double-count.
    greek_coefs = GreekCoefs(αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ)
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, fᵗ=(FT(1) - c₀))
end

@doc raw"""
    truncate_phase(mod::δBGE, aero::AerosolOptics; reportFit=false) -> AerosolOptics

Apply δ-BGE truncation to aerosol Greek coefficients.

The method removes/approximates the forward peak using a least-squares fit over
angles outside `Δ_angle`, then renormalizes with retained scattering fraction
``c_0``:

```math
\beta^t = \frac{c}{c_0},\qquad
\delta^t,\alpha^t,\zeta^t \text{ adjusted consistently from } \beta^t,
\qquad
f^t = 1-c_0.
```

Returns a new [`AerosolOptics`](@ref) with truncated coefficients and updated
`fᵗ`.
"""
function truncate_phase(mod::δBGE, aero::AerosolOptics{FT}; reportFit=false) where {FT}
    (; greek_coefs, ω̃, k) = aero
    (; α, β, γ, δ, ϵ, ζ) = greek_coefs
    (; l_max, Δ_angle) = mod

    l_tr = l_max
    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre(length(β));

    # Reconstruct phase matrix elements:
    scattering_matrix, P, P² = reconstruct_phase(greek_coefs, μ; returnLeg=true)

    (; f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄) = scattering_matrix

    # Find elements that exclude the peak (if wanted!)
    iμ = findall(x -> x < cosd(Δ_angle), μ)

    # Prefactor for P2:
    fac = zeros(FT,l_tr);
    for l = 2:l_tr - 1
        fac[l + 1] = sqrt(FT(1) / ( ( l - FT(1)) * l * (l + FT(1)) * (l + FT(2)) ));
    end

    # Create subsets (for Ax=y weighted least-squares fits):
    y₁₁ = view(f₁₁, iμ)
    y₁₂ = view(f₁₂, iμ)
    y₃₄ = view(f₃₄, iμ)
    
    #= for β
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ Pᵢ(μₖ)Pⱼ(μₖ)/f₁₁²(μₖ), xᵢ=cᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ Pᵢ(μₖ)/f₁₁(μₖ)
    =#   
    A = zeros(l_tr, l_tr)
    x = zeros(l_tr)
    b = zeros(l_tr)

    for i = 1:l_tr
        b[i] = sum(w_μ.*P[:,i]./f₁₁)
        A[i,i] = sum(w_μ.*(P[:,i]./f₁₁).^2)
        for j = i+1:l_tr
            A[i,j] = sum(w_μ.*P[:,i].*P[:,j]./(f₁₁.^2))
            A[j,i] = A[i,j]
        end
    end
    cl = A\b # Julia backslash operator for SVD (just like Matlab);
    # B in δ-BGR (β)
    if reportFit
        println("Errors in δ-BGE fits:")
        mod_y = convert.(FT, A * cl)
        @show StatsBase.rmsd(mod_y, y₁₁; normalize=true)
    end

    #= for γ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)facⱼP²ⱼ(μₖ)/f₁₂²(μₖ), xᵢ=gᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)/f₁₂(μₖ)
    =#  
    A = zeros(l_tr, l_tr)
    x = zeros(l_tr)
    b = zeros(l_tr)

    for i = 3:l_tr
        b[i] = fac[i]*sum(w_μ.*P²[:,i]./f₁₂)
        A[i,i] = (fac[i])^2*sum(w_μ.*(P²[:,i]./f₁₂).^2)
        for j = i+1:l_tr
            A[i,j] = fac[i]*fac[j]*sum(w_μ.*P²[:,i].*P²[:,j]./(f₁₂.^2))
            A[j,i] = A[i,j]
        end
    end
    γᵗ = similar(cl); γᵗ[1:2] .=0
    γᵗ[3:end] = A[3:end,3:end] \ b[3:end]   # G in δ-BGE (γ)

    if reportFit
        println("Errors in δ-BGE fits:")
        mod_γ = convert.(FT, B * γᵗ[3:end])
        @show StatsBase.rmsd(mod_γ, y₁₂; normalize=true)
    end
    
    #= for ϵ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)facⱼP²ⱼ(μₖ)/f₁₂²(μₖ), xᵢ=eᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)/f₃₄(μₖ)
    =#  
    A = zeros(l_tr, l_tr)
    x = zeros(l_tr)
    b = zeros(l_tr)

    for i = 3:l_tr
        b[i] = fac[i]*sum(w_μ.*P²[:,i]./f₃₄)
        A[i,i] = (fac[i])^2*sum(w_μ.*(P²[:,i]./f₃₄).^2)
        for j = i+1:l_tr
            A[i,j] = fac[i]*fac[j]*sum(w_μ.*P²[:,i].*P²[:,j]./(f₃₄.^2))
            A[j,i] = A[i,j]
        end
    end
    
    ϵᵗ = similar(cl); ϵᵗ[1:2] .=0
    ϵᵗ[3:end] = A[3:end,3:end] \ b[3:end]   # E in δ-BGE (ϵ)
    
    if reportFit
        println("Errors in δ-BGE fits:")
        mod_ϵ = convert.(FT, B * ϵᵗ[3:end])
        @show StatsBase.rmsd(mod_ϵ, y₃₄; normalize=true)
    end

    # Integrate truncated function for later renormalization (here: fraction that IS still scattered):
    c₀ = FT(cl[1]) # ( w_μ' * (P[:,1:l_max] * cl) ) / 2
    
    # Compute truncated greek coefficients:
    βᵗ = cl / c₀                                    # Eq. 38a, B in δ-BGR (β)
    δᵗ = (δ[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38b, derived from β
    αᵗ = (α[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38c, derived from β
    ζᵗ = (ζ[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38d, derived from β

    # Truncated Greek coefficients only — ω̃ and k pass through. The
    # τ / ω rescaling per Sanghavi & Stephens 2015 Eq. 8 is applied
    # later in the pipeline by `delta_m_forward` (see
    # CoreRT/LayerOpticalProperties/delta_m_truncation.jl): given
    # `(τ, ω̃, fᵗ)` it returns `(τ_mod, ϖ_mod)` with the proper
    # `(1 − fᵗ·ω̃)` and `(1−fᵗ)·ω̃/(1−fᵗ·ω̃)` factors. Re-applying
    # them here would double-count.
    greek_coefs = GreekCoefs(αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ)
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, fᵗ=(FT(1) - c₀))
end
