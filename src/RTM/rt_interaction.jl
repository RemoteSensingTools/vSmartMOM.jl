"Simulates the full atmosphere from n distinct homogeneous layers"
function rt_interaction_helper!(kn::Int,
                                R⁻⁺::AbstractArray{FT,3}, 
                                T⁺⁺::AbstractArray{FT,3}, 
                                R⁺⁻::AbstractArray{FT,3}, 
                                T⁻⁻::AbstractArray{FT,3}, 
                                r⁻⁺::AbstractArray{FT,3}, 
                                t⁺⁺::AbstractArray{FT,3}, 
                                r⁺⁻::AbstractArray{FT,3}, 
                                t⁻⁻::AbstractArray{FT,3},
                                I_static::AbstractArray) where {FT}
    
    # ToDo: Important output from this routine is R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ (can be renamed to 𝐓⁻⁻, etc later)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)

    # kn = 1: no scattering in either the added layer or composite layer.
    # kn = 2: composite layer has no scattering but added layer does.
    # kn = 3: composite layer has scattering but added layer does not.
    # kn = 4: both composite layer and added layer have scattering.
    
    # ----------------

    if kn == 1

        # No scattering in either the added layer or the composite layer.
        T⁻⁻[:] = t⁻⁻ ⊠ T⁻⁻
        T⁺⁺[:] = t⁺⁺ ⊠ T⁺⁺
        return nothing

    elseif kn == 2

        # No scattering in inhomogeneous composite layer.
        # Scattering in homogeneous layer, added to bottom of the composite layer.
        # Produces a new, scattering composite layer.
        R⁻⁺[:] = T⁻⁻ * r⁻⁺ * T⁺⁺
        R⁺⁻[:] = r⁺⁻
        T⁺⁺[:] = t⁺⁺ * T⁺⁺
        T⁻⁻[:] = T⁻⁻ * t⁻⁻
        return nothing 

    elseif kn == 3

        # Scattering in inhomogeneous composite layer.
        # no scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        T⁺⁺[:] = t⁺⁺ * T⁺⁺
        T⁻⁻[:] = T⁻⁻ * t⁻⁻
        R⁺⁻[:] = t⁺⁺ * R⁺⁻ * t⁻⁻
        return nothing 

    elseif kn == 4

        # Scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.

        # Used to store `inv(I - R⁺⁻ * r⁻⁺) * T⁺⁺`
        tmp_inv = similar(t⁺⁺)

        # Compute and store `inv(I - R⁺⁻ * r⁻⁺) * T⁺⁺`
        batch_solve!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺, T⁺⁺)

        R⁻⁺[:] = R⁻⁺ + (T⁻⁻ ⊠ r⁻⁺ ⊠ tmp_inv)
        T⁺⁺[:] = t⁺⁺ ⊠ tmp_inv

        # Repeating for mirror-reflected directions

        # Compute and store `inv(I - r⁻⁺ * R⁺⁻) * t⁻⁻`
        batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻, t⁻⁻)

        R⁺⁻[:] = r⁺⁻ + t⁺⁺ ⊠ R⁺⁻ ⊠ tmp_inv
        T⁻⁻[:] = T⁺⁺ ⊠ tmp_inv

        return nothing
        
    else 
        error("kn is ($kn), must be in (1, 2, 3, 4)")
    end

end

function rt_interaction!(kn::Int,
                         R⁻⁺::AbstractArray{FT,3}, T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
                         I_static::AbstractArray) where {FT}

    rt_interaction_helper!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, I_static)
    synchronize()

end