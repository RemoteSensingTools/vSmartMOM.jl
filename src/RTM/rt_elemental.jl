"Elemental single-scattering layer"
# function rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
#     # ToDo: Main output is r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
#     # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM

#     #dτ: optical depth of elemental layer
#     #ϖ: single scattering albedo of elemental layer
#     #bb: thermal source function at the upper boundary of the elemental layer
#     #m: fourier moment
#     #n: layer of which this is an elemental
#     #ndoubl: number of doubling computations needed to progress from the elemental layer to the full homogeneous layer n
#     #scatter: flag indicating scattering

#     # Create temporary matrices
#     I_static = one(similar(Z⁻⁺))
#     aux1 = similar(Z⁻⁺)

#     if scatter
#         #TODO: import vector containing quadrature cosines qp_μ of length Nquad4
#         #TODO: import vector containing quadrature weights wt_μ of length Nquad4
#         #TODO: construct composite, post-truncation dτ=τ/2^{ndoubl} , ϖ, p⁺⁺, p⁻⁺ matrices and import them here
        
#         qp_μ4 = reduce(vcat, (fill.(qp_μ,[pol_type.n])))
#         wt_μ4 = reduce(vcat, (fill.(wt_μ,[pol_type.n])))
#         Nquadn = length(qp_μ4)

#         wct = m==0 ? 0.50 * ϖ * wt_μ4  : 0.25 * ϖ * wt_μ4


#         # The following section performs: 
#         # r⁻⁺[:] = Diagonal(1 ./ qp_μ4) * Z⁻⁺ * Diagonal(wct) * dτ
#         # t⁺⁺[:] = I - (Diagonal(1 ./ qp_μ4) * (I - Z⁺⁺ * Diagonal(wct)) * dτ)

#         # Get the diagonal matrices first
#         d_qp = Diagonal(1 ./ qp_μ4)
#         d_wct = Diagonal(wct)

#         # Calculate r⁻⁺
#         mul!(aux1, d_qp, Z⁻⁺)        # Diagonal(1 ./ qp_μ4) * Z⁻⁺
#         mul!(r⁻⁺, aux1, d_wct * dτ)  # r⁻⁺ = (Diagonal(1 ./ qp_μ4) * Z⁻⁺) * Diagonal(wct) * dτ

#         # Calculate t⁺⁺
#         mul!(aux1, Z⁺⁺, d_wct)      # Z⁺⁺ * Diagonal(wct)
#         @. aux1 = I_static - aux1   # (I - Z⁺⁺ * Diagonal(wct))
#         rmul!(aux1, dτ)             # (I - Z⁺⁺ * Diagonal(wct)) * dτ
#         mul!(aux1, d_qp, aux1)      # (Diagonal(1 ./ qp_μ4) * (I - Z⁺⁺ * Diagonal(wct)) * dτ)
#         @. t⁺⁺ = I_static - aux1    # t⁺⁺ = I - ⤴

#         #test = I - Diagonal(1 ./ qp_μ4) * dτ
#         #@show t⁺⁺[1],  test[1], 1 ./ qp_μ4[1]
        
#         if ndoubl<1
#             mul!(aux1, D, r⁻⁺)
#             mul!(r⁺⁻, aux1, D)

#             mul!(aux1, D, t⁺⁺)
#             mul!(t⁻⁻, aux1, D)
#             #=for iμ = 1:Nquad4, jμ = 1:Nquad4
#                 # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
#                 i=mod(iμ-1,4)
#                 j=mod(jμ-1,4)
#                 if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
#                     r⁺⁻[iμ,jμ] = r⁻⁺[iμ,jμ]
#                     t⁻⁻[iμ,jμ] = t⁺⁺[iμ,jμ]
#                 else
#                     r⁺⁻[iμ,jμ] = - r⁻⁺[iμ,jμ]
#                     t⁻⁻[iμ,jμ] = - t⁺⁺[iμ,jμ]
#                 end
#             end =#  
#         else
#             #For doubling, transform R->DR, where D = Diagonal{1,1,-1,-1}
#             mul!(aux1, D, r⁻⁺)
#             @. r⁻⁺ = aux1
#             #= for iμ = 1:Nquad4;
#                 i=mod(iμ-1,4)    
#                 if (i>=2)
#                     r⁻⁺[iμ,:] = - r⁻⁺[iμ,:]
#                 end
#             end =#  
#         end
#     else
#         t⁺⁺[:] = Diagonal{exp(-τ./qp_μ4)}
#         t⁻⁻[:] = Diagonal{exp(-τ./qp_μ4)}
#     end 
#     return nothing 
# end

function rt_elemental_helper!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                       ndoubl::Int, qp_μ, wt_μ, 
                       r⁻⁺::AbstractArray{FT,3}, 
                       t⁺⁺::AbstractArray{FT,3}, 
                       r⁺⁻::AbstractArray{FT,3}, 
                       t⁻⁻::AbstractArray{FT,3}, 
                       D::AbstractArray{Float32,2}, 
                       I_static::AbstractArray) where {FT}

    r⁻⁺_ = r⁻⁺
    t⁺⁺_ = t⁺⁺
    r⁺⁻_ = r⁺⁻
    t⁻⁻_ = t⁻⁻


    Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)

    aux1 = similar(Z⁻⁺_)


    qp_μ4 = reduce(vcat, (fill.(qp_μ,[pol_type.n])))
    wt_μ4 = reduce(vcat, (fill.(wt_μ,[pol_type.n])))

    D_ = repeat(D, 1, 1, 1)

    Nquadn = length(qp_μ4)

    wct = m==0 ? 0.50 * ϖ * wt_μ4  : 0.25 * ϖ * wt_μ4

    d_qp = Array(Diagonal(1 ./ qp_μ4)) 
    d_wct = Array(Diagonal(wct))

    aux1 = d_qp ⊠ Z⁻⁺_
    r⁻⁺_ = aux1 ⊠ (d_wct * dτ)

    t⁺⁺_ = I_static .- (d_qp ⊠ ((I_static .- Z⁺⁺_ ⊠ d_wct) * dτ))

    if ndoubl<1
        r⁺⁻_ = D_ ⊠ r⁻⁺_ ⊠ D_
        t⁻⁻_ = D_ ⊠ t⁺⁺_ ⊠ D_
    else
        r⁻⁺_ = D_ ⊠ r⁻⁺_
    end

    @. r⁻⁺ = r⁻⁺_[:,:,1]
    @. t⁺⁺ = t⁺⁺_[:,:,1]
    @. r⁺⁻ = r⁺⁻_[:,:,1]
    @. t⁻⁻ = t⁻⁻_[:,:,1]
end

function rt_elemental!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, 
                              ndoubl::Int, qp_μ, wt_μ, 
                              r⁻⁺::AbstractArray{FT,3}, 
                              t⁺⁺::AbstractArray{FT,3}, 
                              r⁺⁻::AbstractArray{FT,3}, 
                              t⁻⁻::AbstractArray{FT,3}, 
                              D::AbstractArray{Float32,2}, 
                              I_static::AbstractArray) where {FT}

    rt_elemental_helper!(pol_type, dτ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, qp_μ, wt_μ, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    synchronize()
end