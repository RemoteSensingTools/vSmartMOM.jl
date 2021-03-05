



function default_parameters(FT::DataType=Float64)

    λ = FT(0.770) 
    depol = FT(0.0)
    l_trunc = 20 
    Δ_angle = FT(2.0)
    polarization_type = Stokes_IQU{FT}()

    vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
    vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
    sza = FT(60.)

    obs_alt = FT(1000.0)

    μ            = FT(1.3) #characteristic radius [μm]
    σ            = FT(2.0) #characteristic width
    r_max        = FT(30.0) #maximum radius in distribution [μm]
    nquad_radius = 2500
    nᵣ           = FT(1.3)
    nᵢ           = FT(0.00000001)

    # Aerosol vertical distribution profiles
    p₀          = FT(90000) # Pressure peak [Pa]
    σp          = FT(5000.) # Pressure peak width [Pa]

    file = "/Users/sanghavi/data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
    timeIndex = 2

    lat = 34.1377;
    lon = -118.1253;

    profile_reduction_n = 20; 

    decomp_type = NAI2()

    quadrature_type = RadauQuad() 

    architecture = default_architecture;

    SFI = 1 #Suniti: 0:= DNI, 1:= SFI
    grid_start = FT(1e7 / 774)
    grid_end = FT(1e7 / 757)
    grid_n = 100

    broadening_function = Voigt()
    wing_cutoff = 100
    CEF=HumlicekWeidemann32SDErrorFunction()
    vmr = FT(0.21)

    max_m = 3

    return vSmartMOM_Parameters(FT,
                                λ,
                                depol,
                                l_trunc,
                                Δ_angle,
                                max_m,
                                polarization_type,
                                obs_alt,
                                sza,
                                vza,
                                vaz,
                                μ,
                                σ,
                                r_max,
                                nquad_radius,
                                nᵣ,
                                nᵢ,
                                p₀,
                                σp,
                                file,
                                timeIndex,

                                lat,
                                lon,

                                profile_reduction_n,
                                grid_start,
                                grid_end,
                                grid_n,

                                broadening_function,
                                wing_cutoff,
                                CEF,
                                vmr,
                                decomp_type,
                                quadrature_type,
                                architecture)

end


function default_model(params::vSmartMOM_Parameters)

    
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)
    obs_geom = ObsGeometry(params.obs_alt, params.sza, params.vza, params.vaz)

    size_distribution = LogNormal(log(params.μ), log(params.σ))

    # Create the aerosols (needs to be generalized through loops):
    aerosol = make_univariate_aerosol(size_distribution, params.r_max, params.nquad_radius, params.nᵣ, params.nᵢ)

    # Read profile (and generate dry/wet VCDs per layer)
    profile_hr = vSmartMOM.read_atmos_profile(params.file, params.lat, params.lon, params.timeIndex);
    profile = vSmartMOM.reduce_profile(params.profile_reduction_n, profile_hr);

    mie_model = make_mie_model(params.decomp_type, aerosol, params.λ, params.polarization_type, truncation_type)
    aerosol_optics = compute_aerosol_optical_properties(mie_model, params.float_type);

    # Truncate:
    aerosol_optics_trunc = Scattering.truncate_phase(truncation_type, aerosol_optics; reportFit=true)

    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)

    Nquad, qp_μ, wt_μ = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom);

    τRayl =  getRayleighLayerOptProp(profile.psurf / 100, params.λ, params.depol, profile.vcd_dry);
    ϖRayl = ones(params.float_type, length(τRayl));

    # Compute Naer aerosol optical thickness profiles
    τAer = params.float_type(0.0) * vSmartMOM.getAerosolLayerOptProp(1.0, params.p₀, params.σp, profile.p_levels)
    @show τAer, sum(τAer)
    # Can be done with arbitrary length later:
    # τAer = FT(0.2) * τAer_1; # [τAer_1 τAer_2]
    ϖAer = params.float_type[aerosol_optics.ω̃]; 
    # fᵗ   = FT[aerosol_optics_trunc_aero1.fᵗ]; # [aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];

    # aerosol_optics = [aerosol_optics_trunc_aero1] # 


    grid = range(params.grid_start, params.grid_end, length=params.grid_n);
    τ_abs = zeros(params.float_type, length(grid), length(profile.p));

    hitran_data = read_hitran(artifact("O2"), iso=1)
    absorption_model = make_hitran_model(hitran_data, 
                                         params.broadening_function, 
                                         wing_cutoff = params.wing_cutoff, 
                                         CEF = params.CEF, 
                                         architecture = params.architecture, 
                                         vmr = params.vmr)
    
    τ_abs = zeros(params.float_type, length(grid), length(profile.p));
    
    compute_absorption_profile!(τ_abs, absorption_model, grid, profile);

    return vSmartMOM_Model(params, truncation_type, mie_model, aerosol_optics, aerosol_optics_trunc, greek_rayleigh, Nquad, qp_μ, wt_μ, τRayl, ϖRayl, τAer, ϖAer, obs_geom, aerosol, profile, absorption_model,τ_abs)

end