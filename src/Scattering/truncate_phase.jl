#=
 
This file specifies how to truncate the AerosolOptics struct, given the truncation type
 
=#

"""
    $(FUNCTIONNAME)(mod::δBGE, aero::AerosolOptics))
    
Returns the truncated aerosol optical properties as [`AerosolOptics`](@ref) 
- `mod` a [`δBGE`](@ref) struct that defines the truncation order (new length of greek parameters) and exclusion angle
- `aero` a [`AerosolOptics`](@ref) set of aerosol optical properties that is to be truncated
"""
function truncate_phase_lowconf(mod::δBGE, aero::AerosolOptics{FT}; reportFit=false) where {FT}
    @unpack greek_coefs, ω̃, k = aero
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
    @unpack l_max, Δ_angle =  mod


    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre(length(β));

    # Reconstruct phase matrix elements:
    scattering_matrix, P, P² = reconstruct_phase(greek_coefs, μ; returnLeg=true)

    @unpack f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = scattering_matrix

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

    # Adjust scattering and extinction cross section!
    greek_coefs = GreekCoefs(αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ  )
  
    # C_sca  = (ω̃ * k);
    # C_scaᵗ = C_sca * c₀; 
    # C_ext  = k - (C_sca - C_scaᵗ);
    #@show typeof(ω̃), typeof(k),typeof(c₀)
    # return AerosolOptics(greek_coefs = greek_coefs, ω̃=C_scaᵗ / C_ext, k=C_ext, fᵗ = 1-c₀) 
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, fᵗ=(FT(1) - c₀))
end

"""
$(FUNCTIONNAME)(mod::δBGE, aero::AerosolOptics))
    
Returns the truncated aerosol optical properties as [`AerosolOptics`](@ref) 
- `mod` a [`δBGE`](@ref) struct that defines the truncation order (new length of greek parameters) and exclusion angle
- `aero` a [`AerosolOptics`](@ref) set of aerosol optical properties that is to be truncated
"""
function truncate_phase(mod::δBGE, aero::AerosolOptics{FT}; reportFit=false) where {FT}
    @unpack greek_coefs, ω̃, k = aero
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
    @unpack l_max, Δ_angle =  mod

    l_tr = l_max
    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre(length(β));
    μ = convert.(FT, μ); w_μ = convert.(FT, w_μ);
    #show μ
    # Reconstruct phase matrix elements:
    scattering_matrix, P, P² = reconstruct_phase(greek_coefs, μ; returnLeg=true)

    @unpack f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = scattering_matrix

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
    A = zeros(FT,l_tr, l_tr)
    x = zeros(FT,l_tr)
    b = zeros(FT,l_tr)

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
    A = zeros(FT,l_tr, l_tr)
    x = zeros(FT,l_tr)
    b = zeros(FT,l_tr)

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

    # Adjust scattering and extinction cross section!
    greek_coefs = GreekCoefs(αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ  )
  
    # C_sca  = (ω̃ * k);
    # C_scaᵗ = C_sca * c₀; 
    # C_ext  = k - (C_sca - C_scaᵗ);
    #@show typeof(ω̃), typeof(k),typeof(c₀)
    # return AerosolOptics(greek_coefs = greek_coefs, ω̃=C_scaᵗ / C_ext, k=C_ext, fᵗ = 1-c₀) 
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, fᵗ=(FT(1) - c₀))
end