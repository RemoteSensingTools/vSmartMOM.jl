
# Prototype doubling methods, compute homogenous layer matrices from its elemental layer in 
# `ndoubl` doubling steps

function rt_doubling_helper!(ndoubl::Int, added_layer::AddedLayer,
                             D::AbstractArray{FT,3},
                             I_static::AbstractArray) where {FT}

    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer

    # # ToDo: Important output doubling applied to elemental layer, using same variables 
    # r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ (can be renamed to t⁺⁺, etc)

    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM but 
    # using simplifications in eq. 29-32)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing

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

function rt_doubling!(ndoubl::Int, added_layer::AddedLayer,
                      D::AbstractArray{FT,3}, I_static::AbstractArray) where {FT}

    rt_doubling_helper!(ndoubl, added_layer, D, I_static)
    synchronize()
end
