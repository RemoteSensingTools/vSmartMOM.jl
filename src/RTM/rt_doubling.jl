
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function rt_doubling_helper!(ndoubl::Int, 
                             r⁻⁺::AbstractArray{FT,3}, 
                             t⁺⁺::AbstractArray{FT,3}, 
                             r⁺⁻::AbstractArray{FT,3}, 
                             t⁻⁻::AbstractArray{FT,3},
                             D::AbstractArray{Float32,2},
                             I_static::AbstractArray) where {FT}

    # # ToDo: Important output doubling applied to elemental layer, using same variables 
    # r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)

    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but 
    # using simplifications in eq. 29-32)

    r⁻⁺_ = r⁻⁺ # repeat(r⁻⁺, 1, 1, 1)
    t⁺⁺_ = t⁺⁺ # repeat(t⁺⁺, 1, 1, 1)
    r⁺⁻_ = r⁺⁻ # repeat(r⁻⁺, 1, 1, 1)
    t⁻⁻_ = t⁻⁻ # repeat(r⁻⁺, 1, 1, 1)
    D_ = repeat(D, 1, 1, 1)

    aux1 = similar(t⁺⁺_)
    aux2 = similar(t⁺⁺_)
    aux3 = similar(t⁺⁺_)

    if (ndoubl == 0)
        return nothing
    end

    # Loop over each step
    for n = 1:ndoubl

        # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺
        aux1 = I_static .- r⁻⁺_ ⊠ r⁻⁺_     # (I - r⁻⁺ * r⁻⁺)      
        batch_solve!(aux2, aux1, t⁺⁺_)   # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺

        # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1
        aux1 = r⁻⁺_ ⊠ aux2               # r⁻⁺ * M1
        aux3 = t⁺⁺_ ⊠ aux1               # t⁺⁺ * r⁻⁺ * M1
        r⁻⁺_  = r⁻⁺_ + aux3               # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1

        # t⁺⁺[:] = t⁺⁺ * M1 
        aux1 = t⁺⁺_ ⊠ aux2           # t⁺⁺ * M1 
        t⁺⁺_  = aux1                   # t⁺⁺[:] = t⁺⁺ * M1 

    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    r⁻⁺_ = D_ ⊠ r⁻⁺_ ⊠ D_
    # Using r⁺⁻ = Dr⁻⁺D
    r⁺⁻_ = D_ ⊠ r⁻⁺_ ⊠ D_
    # Using t⁻⁻ = Dt⁺⁺D
    t⁻⁻_ = D_ ⊠ t⁺⁺_ ⊠ D_

    @. r⁻⁺ = r⁻⁺_[:,:,1]
    @. t⁺⁺ = t⁺⁺_[:,:,1]
    @. r⁺⁻ = r⁺⁻_[:,:,1]
    @. t⁻⁻ = t⁻⁻_[:,:,1]

    return nothing 
end

function rt_doubling!(ndoubl::Int, 
                      r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
                      r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
                      D::AbstractArray{Float32,2}, I_static::AbstractArray) where {FT}

    rt_doubling_helper!(ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D, I_static)
    synchronize()
end

"minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
#The starting optical depth dτ is also determined from its maximum possible value, τ"
function doubling_number(dτ_max, τ_end) # check if τ_end can be replaced by τ_end*ϖ for absorbing atmospheres
    FT = eltype(dτ_max)
    # @show FT, eltype(τ_end), dτ_max, τ_end
    # minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
    # The starting optical depth dτ is also determined from its maximum possible value, dτ_max
    if τ_end <= dτ_max
        dτ = τ_end
        ndoubl = 0
        return dτ, ndoubl
    else
        q1 = log10(2.0)
        q2 = log10(dτ_max)
        q3 = log10(τ_end)
        tlimit = (q3 - q2) / q1
        nlimit = floor(Int, tlimit)
        diff = tlimit - nlimit
        if diff < eps(FT)
            dτ = dτ_max
            ndoubl = nlimit
        else
            ndoubl = nlimit + 1       
            x = q3 - q1 * ndoubl
            dτ = 10.0^x
        end 
        return dτ, ndoubl
    end
end

# function rt_doubling!(dτ, τ_total, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, D)
#     # # ToDo: Important output doubling applied to elemental layer, using same variables r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)
#     # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but using simplifications in eq. 29-32)
#     Nquadn = size(r⁻⁺, 1)
#     if (ndoubl == 0)
#         @assert (τ_total == dτ * 2^ndoubl)
#         return nothing 
#     end
#     τ_total = dτ

#     # Create temporary matrices
#     I_static = one(similar(t⁺⁺))
#     aux1 = similar(t⁺⁺)
#     aux2 = similar(t⁺⁺)
#     aux3 = similar(t⁺⁺)

#     for n = 1:ndoubl

#         # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺
#         mul!(aux1, r⁻⁺, r⁻⁺)            # r⁻⁺ * r⁻⁺
#         @. aux1 = I_static - aux1       # (I - r⁻⁺ * r⁻⁺)
#         ldiv!(aux2, qr!(aux1), t⁺⁺)     # M1 = (I - r⁻⁺ * r⁻⁺) \ t⁺⁺

#         # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1
#         mul!(aux1, r⁻⁺, aux2)           # r⁻⁺ * M1
#         mul!(aux3, t⁺⁺, aux1)           # t⁺⁺ * r⁻⁺ * M1
#         @. r⁻⁺ = r⁻⁺ + aux3             # r⁻⁺[:] = r⁻⁺ + t⁺⁺ * r⁻⁺ * M1

#         # t⁺⁺[:] = t⁺⁺ * M1 
#         mul!(aux1, t⁺⁺, aux2)           # t⁺⁺ * M1 
#         @. t⁺⁺ = aux1                   # t⁺⁺[:] = t⁺⁺ * M1 

#         τ_total = 2 * τ_total
#     end
#     # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
#     mul!(aux1, D, r⁻⁺)
#     @. r⁻⁺ = aux1

#     # Using r⁺⁻ = Dr⁻⁺D
#     mul!(aux1, D, r⁻⁺)
#     mul!(r⁺⁻, aux1, D)

#     # Using t⁻⁻ = Dt⁺⁺D
#     mul!(aux1, D, t⁺⁺)
#     mul!(t⁻⁻, aux1, D)
#     #= for iμ = 1:Nquad4, jμ = 1:Nquad4
#         # That "4" and Nquad4 needs to be dynamic, coming from the PolType struct.
#         i=mod(iμ-1,4)
#         j=mod(jμ-1,4)
#         #@show i,j
#         if (i>=2)
#             r⁻⁺[iμ,jμ] = - r⁻⁺[iμ, jμ]
#         end
#         if ((i<=1)&(j<=1)) | ((i>=2)&(j>=2))
#             r⁺⁻[iμ,jμ] = r⁻⁺[iμ,jμ]
#             t⁻⁻[iμ,jμ] = t⁺⁺[iμ,jμ]
#         else
#             r⁺⁻[iμ,jμ] = - r⁻⁺[iμ,jμ]
#             t⁻⁻[iμ,jμ] = - t⁺⁺[iμ,jμ]
#         end
#     end =#

#     @assert (τ_total == dτ * 2^ndoubl)
#     return nothing 
# end