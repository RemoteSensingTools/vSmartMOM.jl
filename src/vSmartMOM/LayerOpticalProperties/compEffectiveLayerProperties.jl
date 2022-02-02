function constructCoreOpticalProperties(RS_type, iBand, m, model)
    @unpack τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh,  = model
    @assert iBand ≤ length(τ_rayl) "iBand exceeded number of bands"
    
    pol_type = model.params.polarization_type
    # Do this in CPU space only first:
    arr_type = Array
    # Quadrature points:
    μ = Array(model.quad_points.qp_μ )
    # Number of Aerosols:
    nAero = size(τ_aer[iBand],1)
    
    # Rayleigh Z matrix:
    Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, greek_rayleigh, m, arr_type = arr_type);

    # Create Rayleight Core properties per layer
    rayl = CoreScatteringOpticalProperties.(τ_rayl[iBand][:], [RS_type.ϖ_Cabannes], [Rayl𝐙⁺⁺], [Rayl𝐙⁻⁺])
    
    # Initiate combined properties with rayleigh
    combo = rayl

    # Loop over all aerosol types:
    for i=1:nAero
        # Precomute Z matrices per type (constant per layer)
        AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, aerosol_optics[iBand][i].greek_coefs, m, arr_type=arr_type)
        # Generate Core optical properties for Aerosols i
        aer   = createAero.(τ_aer[iBand][i,:], [aerosol_optics[iBand][i]], [AerZ⁺⁺], [AerZ⁻⁺])
        # Mix with previous Core Optical Properties
        combo = combo .+ aer
    end

    # Somewhere here we can add canopy later as well!
    ###

    # fScattRayleigh:
    fScattRayleigh = [rayl[i].τ * rayl[i].ϖ / combo[i].τ for i=1:length(combo)]

    # Create Core Optical Properties merged with trace gas absorptions:
    combo = combo .+ [CoreAbsorptionOpticalProperties(τ_abs[iBand][:,i]) for i=1:length(combo)]
    #aType = array_type(model.params.architecture)
    #combo2 = [CoreScatteringOpticalProperties(aType(combo[i].τ),aType(combo[i].ϖ), aType(combo[i].Z⁺⁺), aType(combo[i].Z⁻⁺)) for i in eachindex(combo)]
    # Need to check how to convert to GPU later as well!
    return combo,fScattRayleigh
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end