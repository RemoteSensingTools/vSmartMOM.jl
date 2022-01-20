

"""
    $(FUNCTIONNAME)(RS_type, λ, grid_in)

Returns the Raman SSA per layer at reference wavelength `λ` from nearby source wavelengths for RRS and a (currently) single incident wavelength for VRS/RVRS
(N₂,O₂ atmosphere, i.e. terrestrial)

Input: 
    - `RS_type` Raman scattering type (RRS/RVRS/VRS)
    - `λ` wavelength in `[μm]`
    - `grid_in` wavenumber grid with equidistant gridpoints
"""
function getRamanSSProp(RS_type::noRS, λ, grid_in)

    return nothing
end

function getRamanSSProp!(RS_type::Union{VS_0to1,VS_1to0}, λ, grid_in)
    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    λ_scatt = 2.e7/(grid_in[1]+grid_in[end])     
    #determine Rayleigh scattering cross-section at mean scattered wavelength
    atmo_σ_Rayl_scatt = compute_optical_Rayl!(λ_scatt, n2, o2)
    # determine Rayleigh scattering cross-section at incident wavelength λ 
    atmo_σ_Rayl = compute_optical_Rayl!(λ, n2, o2)
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
    index_VRSgrid_out, atmo_σ_VRS, index_RVRSgrid_out, atmo_σ_RVRS = 
        compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
    RS_type.ϖ_VRS = atmo_σ_VRS/atmo_σ_Rayl
    RS_type.i_VRS = index_VRSgrid_out
    RS_type.ϖ_RVRS = atmo_σ_RVRS/atmo_σ_Rayl
    RS_type.i_RVRS = index_RVRSgrid_out
    RS_type.k_Rayl_scatt = atmo_σ_Rayl_scatt/atmo_σ_Rayl
    RS_type.n_Raman = 1;
    return nothing
end
#function getRamanLayerSSA(RS_type::VRS_1to0, T, λ, grid_in)
#    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at single monochromatic wavelength λ of the spectral band (assumed constant throughout the band)
#    atmo_σ_Rayl = compute_optical_Rayl!(λ, n2, o2)
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
#    index_VRSgrid_out, atmo_σ_VRS_1to0, index_RVRSgrid_out, atmo_σ_RVRS_1to0 = 
#        compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
#    ϖ_VRS = atmo_σ_VRS_1to0/atmo_σ_Rayl
#    i_VRS = index_VRSgrid_out
#    ϖ_RVRS = atmo_σ_RVRS_1to0/atmo_σ_Rayl
#    i_RVRS = index_RVRSgrid_out
#    return ϖ_RVRS, i_RVRS, ϖ_VRS, i_VRS
#end

function getRamanSSProp!(RS_type::RRS, λ, grid_in) 
    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at central wavelength λ of the spectral band (assumed constant throughout the band)
    atmo_σ_Rayl = compute_optical_Rayl!(λ, n2, o2)
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
    index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
    RS_type.ϖ_λ₁λ₀ = atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
    RS_type.i_λ₁λ₀ = index_raman_grid[end:-1:1]
    RS_type.n_Raman = length(RS_type.ϖ_λ₁λ₀)
    return nothing
end

