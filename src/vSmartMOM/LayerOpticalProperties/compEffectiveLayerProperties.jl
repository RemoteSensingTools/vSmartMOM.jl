function constructCoreOpticalProperties(RS_type, iBand, iz, m, τRayl, τAer, τ_abs, aerosol_optics, greek_rayleigh, pol_type, qp_μ)

    # Quadrature points:
    μ = Array(qp_μ)
    # Number of Aerosols:
    nAero = size(τAer[iBand],1)
    arr_type = Array
    # Rayleigh Z matrix:
    Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, greek_rayleigh, m, arr_type = arr_type);

    rayl = CoreScatteringOpticalProperties(τRayl[iBand][iz], RS_type.ϖ_Cabannes, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺)
    #rayl = CoreScatteringOpticalProperties(τRayl[iBand][iz], [RS_type.ϖ_Cabannes], [Rayl𝐙⁺⁺], [Rayl𝐙⁻⁺])
    #iz = 1
    # Create Core Optical Properties of all aerosols combined:
    #Scattering.compute_Z_moments(pol_type, μ, aerosol_optics[iBand][i].greek_coefs, m, arr_type=arr_type)...)
    aer = sum([createAero(τAer[iBand][i,iz], aerosol_optics[iBand][i], Scattering.compute_Z_moments(pol_type, μ, aerosol_optics[iBand][i].greek_coefs, m, arr_type=arr_type)...) for i=1:nAero])
    return rayl + aer
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end