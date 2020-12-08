"Simulates the full atmosphere from n distinct homogeneous layers"
function rt_interaction(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    # ToDo: Important output from this routine is R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ (can be renamed to 𝐓⁻⁻, etc later)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)
    Nquad4 = size(r⁻⁺, 1)
    # kn = 1: no scattering in either the added layer or composite layer.
    # kn = 2: composite layer has no scattering but added layer does.
    # kn = 3: composite layer has scattering but added layer does not.
    # kn = 4: both composite layer and added layer have scattering.
    if kn==1
        # No scattering in either the added layer or the composite layer.
        for iμ in 1:Nquad4
            T⁻⁻[iμ, iμ] = t⁻⁻[iμ, iμ]*T⁻⁻[iμ, iμ]
            T⁺⁺[iμ, iμ] = t⁺⁺[iμ, iμ]*T⁺⁺[iμ, iμ]
        end

        # diag_ind = diagind(T⁻⁻)[1:Nquad4]
        # T⁻⁻[diag_ind] = t⁻⁻[diag_ind] .* T⁻⁻[diag_ind]
        # T⁺⁺[diag_ind] = t⁺⁺[diag_ind] .* T⁺⁺[diag_ind]

        return R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻
    elseif kn==2
        # No scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is added 
        # to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        M1=T⁻⁻
        M2=T⁺⁺
        R⁻⁺ = M1 * r⁻⁺ * M2
        R⁺⁻ = r⁺⁻
        T⁺⁺ = t⁺⁺ * M2
        T⁻⁻ = M1 * t⁻⁻
        return R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻
    elseif kn==3
        # Scattering in inhomogeneous composite layer.
        # no scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        T⁺⁺ = t⁺⁺ * T⁺⁺
        T⁻⁻ = T⁻⁻ * t⁻⁻
        R⁺⁻ = t⁺⁺ * R⁺⁻ * t⁻⁻
        return R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻
    elseif kn==4
        # Scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        # M1 = inv(I - R⁺⁻ * r⁻⁺)
        M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
        t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1 # M1 * T⁺⁺
        t_T⁺⁺ = t⁺⁺ * M1 # M1 * T⁺⁺

        #repeating for mirror-reflected directions
        # M1 = inv(I - r⁻⁺ * R⁺⁻)
        M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
        t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1 # M1 * t⁻⁻
        t_T⁻⁻ = T⁺⁺ * M1 # M1 * t⁻⁻

        return t_R⁻⁺, t_T⁺⁺, t_R⁺⁻, t_T⁻⁻
    end
end