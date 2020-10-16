
function phase_shorthand(x)
    println(x[1], x[2], x[3], x[4])

    μ  = x[1]                # Log mean radius
    σ  = x[2]               # Log stddev of radius
    r_max = 30.0            # Maximum radius
    nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
    nᵣ = x[3]                # Real part of refractive index
    nᵢ = x[4]                # Imag part of refractive index

    size_distribution = LogNormal(log(μ), log(σ))

    # Create the aerosol
    aero = make_univariate_aerosol(size_distribution, r_max, nquad_radius, nᵣ, nᵢ)

    ## 
    ## STEP 2: Create the Mie Calculations model
    ## 

    λ = 0.55   # Incident wavelength
    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10)

    wigner_file_path = "/home/rjeyaram/RadiativeTransfer/src/PhaseFunction/wigner_values.jld"
    model_NAI2 = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_file_path)

    aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);

    return [aerosol_optics_NAI2.greek_coefs.α 
            aerosol_optics_NAI2.greek_coefs.β 
            aerosol_optics_NAI2.greek_coefs.γ 
            aerosol_optics_NAI2.greek_coefs.δ 
            aerosol_optics_NAI2.greek_coefs.ϵ 
            aerosol_optics_NAI2.greek_coefs.ζ
            aerosol_optics_NAI2.ω̃
            aerosol_optics_NAI2.k]
end