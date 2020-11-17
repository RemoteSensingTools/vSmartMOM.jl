"Simulates the full atmosphere from n distinct homogeneous layers"
function rt_interaction!(mo::Stokes_IQUV, kn, n, m, T_pp, T_mm, R_mp, R_pm, t_pp, t_mm, r_mp, r_pm)
    # ToDo: Important output from this routine is R_mp, R_pm, T_pp, T_mm (can be renamed to 𝐓⁻⁻, etc later)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)

    # kn = 1: no scattering in either the added layer or composite layer.
    # kn = 2: composite layer has no scattering but added layer does.
    # kn = 3: composite layer has scattering but added layer does not.
    # kn = 4: both composite layer and added layer have scattering.
    if kn==1
        # No scattering in either the added layer or the composite layer.
        for iμ in 1:Nquad4
            T_mm[iμ, iμ] = t_mm[iμ, iμ]*T_mm[iμ, iμ]
            T_pp[iμ, iμ] = t_pp[iμ, iμ]*T_pp[iμ, iμ]
        end
        return nothing
    elseif kn==2
        # No scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is added 
        # to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        M1=T_mm
        M2=T_pp
        R_mp = M1 * r_mp * M2
        R_pm = r_pm
        T_pp = t_pp * M2
        T_mm = M1 * t_mm
        return nothing
    elseif kn==3
        # Scattering in inhomogeneous composite layer.
        # no scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        T_pp = t_pp * T_pp
        T_mm = T_mm * t_mm
        R_pm = t_pp * R_pm * t_mm
        return nothing
    elseif kn==4
        # Scattering in inhomogeneous composite layer.
        # scattering in homogeneous layer which is 
        # added to the bottom of the composite layer.
        # Produces a new, scattering composite layer.
        M1 = inv(I - R_pm * r_mp)
        t_R_mp = R_mp + T_mm * r_mp * M1 * T_pp
        t_T_pp = t_pp * M1 * T_pp
        #repeating for mirror-reflected directions
        M1 = inv(I - r_mp * R_pm)
        t_R_pm = r_pm + t_pp * R_pm * M1 * t_mm
        t_T_mm = T_pp * M1 * t_mm
        return nothing
    end
end