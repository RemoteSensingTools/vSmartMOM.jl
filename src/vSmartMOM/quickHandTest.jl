using Parameters

pol_type = model.params.polarization_type;           # Polarization type (IQUV)
obs_geom = model.obs_geom::ObsGeometry;  #, # Solar Zenith, Viewing Zenith, Viewing Azimuthal 
τRayl = model.τRayl;        # Rayleigh optical depth 
    #nAer,                 # Number of aerosol species 
τAer =  model.τAer ;                # Aerosol optical depth and single-scattering albedo
qp_μ =model.qp_μ;
wt_μ =   model.wt_μ  ;      # Quadrature points and weights
Ltrunc =  3  ;            # Trunction length for legendre terms
aerosol_optics = model.aerosol_optics;       # AerosolOptics (greek_coefs, ω̃, k, fᵗ)
GreekRayleigh =  model.greek_rayleigh ;       # Greek coefficients of Rayleigh Phase Function
τ_abs =    model.τ_abs  ;           # nSpec x N
architecture = model.params.architecture;

@unpack obs_alt, sza, vza, vaz = obs_geom   # Observational geometry properties
    FT = eltype(sza)                    # Get the float-type to use
    Nz = length(τRayl)                  # Number of vertical slices
    nSpec = size(τ_abs, 1)              # Number of spectral points
    μ0 = cosd(sza)                      # μ0 defined as cos(θ); θ = sza
    iμ0 = vSmartMOM.nearest_point(qp_μ, μ0)       # Find the closest point to μ0 in qp_μ
    arr_type = array_type(architecture)

    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively
    R = zeros(FT, length(vza), pol_type.n, nSpec)
    T = zeros(FT, length(vza), pol_type.n, nSpec)
    R_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    T_SFI = zeros(FT, length(vza), pol_type.n, nSpec)

    # Copy qp_μ "pol_type.n" times
    qp_μN = arr_type(reshape(transpose(repeat(qp_μ, 1, pol_type.n)),pol_type.n*size(qp_μ)[1],1))
    #for i = 1:length(qp_μN)
    #   @show(i,qp_μN[i]) 
    #end
    println("Processing on: ", architecture)
    println("With FT: ", FT)

    #= 
    Loop over number of truncation terms =#
    SFI = true
