#=
 
This file specifies how to truncate the AerosolOptics struct, given the truncation type
 
=#

#=
 
This file specifies how to truncate the AerosolOptics struct, given the truncation type
 
=#

"""
$(FUNCTIONNAME)(mod::δBGE, aero::AerosolOptics))
    
Returns the truncated aerosol optical properties as [`AerosolOptics`](@ref) 
- `mod` a [`δBGE`](@ref) struct that defines the truncation order (new length of greek parameters) and exclusion angle
- `aero` a [`AerosolOptics`](@ref) set of aerosol optical properties that is to be truncated
"""
function truncate_phase(mod::δBGE, aero::AerosolOptics{FT}, lin::dAerosolOptics{FT}; reportFit=false) where {FT}
    @unpack greek_coefs, ω̃, k = aero
    @unpack d_greek_coefs, dω̃, dk = lin   
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
    #@unpack dα, dβ, dγ, dδ, dϵ, dζ = d_greek_coefs
    @unpack l_max, Δ_angle =  mod


    l_tr = l_max
    # Obtain Gauss-Legendre quadrature points and weights for phase function:
    μ, w_μ = gausslegendre(Int((length(β)-1)/2));

    # Reconstruct phase matrix elements:
    scattering_matrix, dscattering_matrix, P, P² = reconstruct_phase(greek_coefs, d_greek_coefs, μ; returnLeg=true)

    @unpack f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = scattering_matrix
    @unpack df₁₁, df₁₂, df₂₂, df₃₃, df₃₄, df₄₄ = dscattering_matrix

    # Find elements that exclude the peak (if wanted!)
    iμ = findall(x -> x < cosd(Δ_angle), μ)

    # Prefactor for P2:
    fac = zeros(FT,l_tr);
    for l = 2:l_tr - 1
        fac[l + 1] = sqrt(FT(1) / ( ( l - FT(1)) * l * (l + FT(1)) * (l + FT(2)) ));
    end

    # Create subsets (for Ax=y weighted least-squares fits):
    #y₁₁ = view(f₁₁, iμ)
    #y₁₂ = view(f₁₂, iμ)
    #y₃₄ = view(f₃₄, iμ)
    
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
        @show StatsBase.rmsd(mod_y, b; normalize=true)
    end

    #= for dβ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ Pᵢ(μₖ)Pⱼ(μₖ)/f₁₁²(μₖ), xᵢ=c'ᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ ((2/f₁₁(μₖ))∑ₗcₗPₗ(μₖ) - 1) Pᵢ(μₖ)f₁₁'(μₖ)/(f₁₁(μₖ))^2
    =#   
    #A = zeros(l_tr, l_tr)
    x = zeros(l_tr)
    b = zeros(l_tr)
    dcl = zeros(4, l_tr)

    for ctr=1:4
        for i = 1:l_tr
            b[i] = sum(((2. / f₁₁).*sum(cl.*(P[:,1:l_tr])', dims=1) .- 1).*w_μ.*P[:,i].*df₁₁[ctr,:]./f₁₁.^2)
        end
        dcl[ctr,:] = A\b # Julia backslash operator for SVD (just like Matlab);
        # B in δ-BGR (β)
        if reportFit
            println("Errors in δ-BGE fits:")
            mod_y = convert.(FT, A * dcl[ctr,:])
            @show StatsBase.rmsd(mod_y, b; normalize=true)
        end
    end


    #= for γ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)facⱼP²ⱼ(μₖ)/f₁₂²(μₖ), xᵢ=γᵢ (as in Sanghavi & Stephens 2015), and
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
        @show StatsBase.rmsd(mod_γ, b; normalize=true)
    end

    #= for dγ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)facⱼP²ⱼ(μₖ)/f₁₂²(μₖ), xᵢ=γ'ᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)f₁₂'(μₖ)/f₁₂²(μₖ) {2∑ₗP²ₗ(μₖ)γₗ/f₁₂(μₖ) - 1}
    =#  
    x = zeros(l_tr)
    b = zeros(l_tr)
    dγᵗ = zeros(4, l_tr)
    for ctr=1:4
        for i = 3:l_tr
            b[i] = fac[i]*sum(((2. / f₁₂).*sum(γᵗ.*(P²[:,1:l_tr])', dims=1) .- 1).*w_μ.*P²[:,i].*df₁₂[ctr,:]./f₁₂.^2)
        end
        dγᵗ[ctr,1:2] .=0
        dγᵗ[ctr,3:end] = A[3:end,3:end] \ b[3:end]   # G in δ-BGE (γ)
        if reportFit
            println("Errors in δ-BGE fits:")
            mod_γ = convert.(FT, A * dγᵗ[ctr,:])
            @show StatsBase.rmsd(mod_γ, b; normalize=true)
        end
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
        mod_ϵ = convert.(FT, A * ϵᵗ)
        @show StatsBase.rmsd(mod_ϵ, b; normalize=true)
    end

    #= for dϵ
       Ax=b, where
       Aᵢⱼ = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)facⱼP²ⱼ(μₖ)/f₃₄²(μₖ), xᵢ=ϵ'ᵢ (as in Sanghavi & Stephens 2015), and
       bᵢ  = ∑ₖ w_μₖ facᵢP²ᵢ(μₖ)f₃₄'(μₖ)/f₃₄²(μₖ) {2∑ₗP²ₗ(μₖ)ϵₗ/f₃₄(μₖ) - 1}
    =#  
    x = zeros(l_tr)
    b = zeros(l_tr)
    dϵᵗ = zeros(4, l_tr)
    for ctr=1:4
        for i = 3:l_tr
            b[i] = fac[i]*sum(((2. / f₃₄).*sum(γᵗ.*(P²[:,1:l_tr])', dims=1) .- 1).*w_μ.*P²[:,i].*df₃₄[ctr,:]./f₃₄.^2)
        end
        dϵᵗ[ctr,1:2] .=0
        dϵᵗ[ctr,3:end] = A[3:end,3:end] \ b[3:end]   # G in δ-BGE (γ)
        if reportFit
            println("Errors in δ-BGE fits:")
            mod_ϵ = convert.(FT, A * dϵᵗ[ctr,:])
            @show StatsBase.rmsd(mod_ϵ, b; normalize=true)
        end
    end

    # Integrate truncated function for later renormalization (here: fraction that IS still scattered):
    c₀ = FT(cl[1]) # ( w_μ' * (P[:,1:l_max] * cl) ) / 2
    
    # Compute truncated greek coefficients:
    βᵗ = cl / c₀                                    # Eq. 38a, B in δ-BGR (β)
    δᵗ = (δ[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38b, derived from β
    αᵗ = (α[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38c, derived from β
    ζᵗ = (ζ[1:l_tr] .- (β[1:l_tr] .- cl)) / c₀    # Eq. 38d, derived from β

    dβᵗ = zeros(4, l_tr)
    dδᵗ = zeros(4, l_tr)
    dαᵗ = zeros(4, l_tr)
    dζᵗ = zeros(4, l_tr)
    for ctr = 1:4
        #dβᵗ[ctr,1] = 0.
        for i = 1:l_tr
            dβᵗ[ctr,:] = (dcl[ctr,:] - βᵗ*dcl[ctr,1])/c₀
            dδᵗ[ctr,:] = (d_greek_coefs[ctr].δ[1:l_tr] .- (d_greek_coefs[ctr].β[1:l_tr] .- dcl[ctr,:]) - δᵗ*dcl[ctr,1]) / c₀    # Eq. 38b, derived from β
            dαᵗ[ctr,:] = (d_greek_coefs[ctr].α[1:l_tr] .- (d_greek_coefs[ctr].β[1:l_tr] .- dcl[ctr,:]) - αᵗ*dcl[ctr,1]) / c₀    # Eq. 38c, derived from β
            dζᵗ[ctr,:] = (d_greek_coefs[ctr].ζ[1:l_tr] .- (d_greek_coefs[ctr].β[1:l_tr] .- dcl[ctr,:]) - ζᵗ*dcl[ctr,1]) / c₀    # Eq. 38d, derived from β
        end
    end 
    # Adjust scattering and extinction cross section!
    greek_coefs = GreekCoefs(αᵗ, βᵗ, γᵗ, δᵗ, ϵᵗ, ζᵗ)
    d_greek_coefs = [GreekCoefs(dαᵗ[i,:], dβᵗ[i,:], dγᵗ[i,:], dδᵗ[i,:], dϵᵗ[i,:], dζᵗ[i,:]) for i=1:4]
    dfᵗ = zeros(4)
    dfᵗ = (FT(1) .- dcl[:,1])
    dk_ref=zeros(4)
    # C_sca  = (ω̃ * k);
    # C_scaᵗ = C_sca * c₀; 
    # C_ext  = k - (C_sca - C_scaᵗ);
    #@show typeof(ω̃), typeof(k),typeof(c₀)
    # return AerosolOptics(greek_coefs = greek_coefs, ω̃=C_scaᵗ / C_ext, k=C_ext, fᵗ = 1-c₀) 
    return AerosolOptics(greek_coefs=greek_coefs, ω̃=ω̃, k=k, k_ref=FT(0), fᵗ=(FT(1) - c₀)), 
        dAerosolOptics(d_greek_coefs, dω̃, dk, dk_ref, dfᵗ)
end

