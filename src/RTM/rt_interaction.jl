"Simulates the full atmosphere from n distinct homogeneous layers"
function rt_interaction!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
                             
    # ToDo: Important output from this routine is R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ (can be renamed to 𝐓⁻⁻, etc later)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)
    Nquadn = size(r⁻⁺, 1)
    # kn = 1: no scattering in either the added layer or composite layer.
    # kn = 2: composite layer has no scattering but added layer does.
    # kn = 3: composite layer has scattering but added layer does not.
    # kn = 4: both composite layer and added layer have scattering.

    # Create temporary matrices
    I_static = one(similar(R⁺⁻))
    aux1 = similar(R⁺⁻)
    aux2 = similar(R⁺⁻)
    aux3 = similar(R⁺⁻)

    if kn==1
        # No scattering in either the added layer or the composite layer.
        T⁻⁻ = t⁻⁻ * T⁻⁻
        T⁺⁺ = t⁺⁺ * T⁺⁺
        
        return nothing 
    elseif kn==2
        # No scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is added 
        # to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        M1=T⁻⁻
        M2=T⁺⁺
        R⁻⁺[:] = M1 * r⁻⁺ * M2
        R⁺⁻[:] = r⁺⁻
        T⁺⁺[:] = t⁺⁺ * M2
        T⁻⁻[:] = M1 * t⁻⁻
        return nothing 
    elseif kn==3
        # Scattering in inhomogeneous composite layer.
        # no scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        T⁺⁺[:] = t⁺⁺ * T⁺⁺
        T⁻⁻[:] = T⁻⁻ * t⁻⁻
        R⁺⁻[:] = t⁺⁺ * R⁺⁻ * t⁻⁻
        return nothing 
    elseif kn==4
        # Scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.

        # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
        mul!(aux1, R⁺⁻, r⁻⁺)        # R⁺⁻ * r⁻⁺
        @. aux1 = I_static - aux1   # (I - R⁺⁻ * r⁻⁺)
        ldiv!(aux2, qr!(aux1), T⁺⁺) # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺

        # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        mul!(aux1, r⁻⁺, aux2)   # r⁻⁺ * M1
        mul!(aux3, T⁻⁻, aux1)   # T⁻⁻ * r⁻⁺ * M1
        @. R⁻⁺ = R⁻⁺ + aux3     # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        
        # t_T⁺⁺ = t⁺⁺ * M1
        mul!(T⁺⁺, t⁺⁺, aux2)

        # Repeating for mirror-reflected directions

        # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
        mul!(aux1, r⁻⁺, R⁺⁻)        # r⁻⁺ * R⁺⁻
        @. aux1 = I_static - aux1   # (I - r⁻⁺ * R⁺⁻)
        ldiv!(aux2, qr!(aux1), t⁻⁻) # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

        # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
        mul!(aux3, R⁺⁻, aux2)   # R⁺⁻ * M1
        mul!(aux1, t⁺⁺, aux3)   # t⁺⁺ * R⁺⁻ * M1
        @. R⁺⁻ = r⁺⁻ + aux1     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

        # t_T⁻⁻ = T⁺⁺ * M1
        mul!(T⁻⁻, T⁺⁺, aux2)
                 
        return nothing 
    end
end