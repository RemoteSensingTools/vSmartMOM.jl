
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function rt_doubling_helper!(ndoubl::Int, 
                             r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
                             r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
                             D::AbstractArray{FT,3},
                             I_static::AbstractArray) where {FT}

    # # ToDo: Important output doubling applied to elemental layer, using same variables 
    # r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)

    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but 
    # using simplifications in eq. 29-32)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing

    println("running!")

    # Used to store `inv(I - r⁻⁺ * r⁻⁺) * t⁺⁺`
    tmp_inv = similar(t⁺⁺)

    # Loop over each step
    for n = 1:ndoubl
        batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ r⁻⁺, t⁺⁺)   
        r⁻⁺[:]  = r⁻⁺ + (t⁺⁺ ⊠ r⁻⁺ ⊠ tmp_inv)
        t⁺⁺[:]  = t⁺⁺ ⊠ tmp_inv
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    r⁻⁺[:] = D ⊠ r⁻⁺ ⊠ D

    # Using r⁺⁻ = Dr⁻⁺D
    r⁺⁻[:] = D ⊠ r⁻⁺ ⊠ D
    
    # Using t⁻⁻ = Dt⁺⁺D
    t⁻⁻[:] = D ⊠ t⁺⁺ ⊠ D

    return nothing 
end

function rt_doubling!(ndoubl::Int, 
                      r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
                      r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
                      D::AbstractArray{FT,3}, I_static::AbstractArray) where {FT}

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