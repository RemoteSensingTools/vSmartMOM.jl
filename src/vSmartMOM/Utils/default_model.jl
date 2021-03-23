



function default_parameters(FT::DataType=Float64)

    # Note, we will need many parameters in Arrays, as we will have `nAer` aerosol types (1 in default)
    # as well as nBands (we can model multiple spectral bands)
    λ_band = FT[0.770] 
    λ_ref = FT(0.770) 
    depol = FT(0.0)
    l_trunc = 20 
    Δ_angle = FT(0.0)
    polarization_type = Stokes_IQU{FT}()

    vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
    vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
    sza = FT(60.)

    obs_alt = FT(1000.0)
    nAer = 1 #Number of aerosol species

    τAer_ref     = FT[0.5] #AOTs at reference wavelength
    μ            = FT[1.3] #characteristic radius [μm]
    σ            = FT[2.0] #characteristic width
    r_max        = FT(30.0) #maximum radius in distribution [μm] #baseline setting
    nquad_radius = 2500 #baseline setting
    nᵣ           = FT[1.3] #TODO: make this a 2D matrix were each species and band has an independent entry
    nᵢ           = FT[0.00000001] #TODO: same as above

    # Aerosol vertical distribution profiles
    p₀          = FT[90000] # Pressure peak [Pa]
    σp          = FT[5000.] # Pressure peak width [Pa]

    # Note: We should change the default profile to an ASCII file and have a simple ascii reader...
    file = "/Users/sanghavi/data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
    timeIndex = 2

    lat = 34.1377;
    lon = -118.1253;

    profile_reduction_n = 20; 

    decomp_type = NAI2()

    quadrature_type = RadauQuad() 

    architecture = default_architecture;

    SFI = 1 #Suniti: 0:= DNI, 1:= SFI
    spec_grid_start = FT[(1e7 / 777)]
    spec_grid_end = FT[(1e7 / 747)]
    spec_grid_n = Integer[2]

    broadening_function = Voigt()
    wing_cutoff = 100
    CEF=HumlicekWeidemann32SDErrorFunction()
    vmr = FT(0.21)
    vmr = FT(0.00)

    # This is only a hard-coded quick fix, eventually we need to compute after how many Fourier 
    # components it converges and then stop (and just set a convergence criterion)
    max_m = l_trunc #temp

    return vSmartMOM_Parameters(FT,
                                λ_band, λ_ref,
                                depol,
                                l_trunc,
                                Δ_angle,
                                max_m,
                                polarization_type,
                                obs_alt,
                                sza,
                                vza,
                                vaz,
                                nAer, 
                                τAer_ref,
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
                                spec_grid_start,
                                spec_grid_end,
                                spec_grid_n,

                                broadening_function,
                                wing_cutoff,
                                CEF,
                                vmr,
                                decomp_type,
                                quadrature_type,
                                architecture)

end


function default_model(params::vSmartMOM_Parameters)
    @unpack λ_band = params
    nBands = length(λ_band)
    obs_geom = ObsGeometry(params.obs_alt, params.sza, params.vza, params.vaz)
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)
    Nquad, qp_μ, wt_μ = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom);
    # Read profile (and generate dry/wet VCDs per layer)
    profile_hr = vSmartMOM.read_atmos_profile(params.file, params.lat, params.lon, params.timeIndex);
    profile = vSmartMOM.reduce_profile(params.profile_reduction_n, profile_hr);

    #Rayleigh
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τRayl = [zeros(params.float_type, length(profile.p)) for i=1:nBands];
    #ϖRayl = similar(ϖRayl)

    # Create empty array (can be changed later)
    # spec_grid[iBand][iSpec]
    spec_grid = [zeros(params.float_type,n) for n in params.spec_grid_n];

    # τ_abs[iBand][iSpec,iZ]
    τ_abs     = [zeros(params.float_type,n, length(profile.p)) for n in params.spec_grid_n]
    # Loop over all bands:
    for ib=1:length(λ_band)
        # Compute Rayleigh properties per layer for `ib` band center 
        τRayl[ib]   = getRayleighLayerOptProp(profile.psurf / 100, λ_band[ib], params.depol, profile.vcd_dry);
        # Define spectral grid per band
        spec_grid[ib] = range(params.spec_grid_start[ib], params.spec_grid_end[ib], length=params.spec_grid_n[ib]);
        
        # This has to be more flexible in the future, so that we can have different trace gases per band, Newton's work 
        # can be checked for that. (ideally a dictionary that defines gases to be used and profiles per band)
        hitran_data = read_hitran(artifact("O2"), iso=1)
        absorption_model = make_hitran_model(hitran_data, 
                                         params.broadening_function, 
                                         wing_cutoff = params.wing_cutoff, 
                                         CEF = params.CEF, 
                                         architecture = params.architecture, 
                                         vmr = params.vmr)
        #τ_abs = zeros(params.float_type, length(λ), length(spec_grid), length(profile.p)); 
        #compute_absorption_profile!(τ_abs[ib], absorption_model, range(params.spec_grid_start[ib], params.spec_grid_end[ib], length=params.spec_grid_n[ib]), profile);
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (params.nAer)) for i=1:nBands];

    # τAer[iBand][iAer,iProfile]
    τAer = [zeros(params.float_type, params.nAer, length(profile.p)) for i=1:nBands];

    # Loop over aerosol type
    for iaer=1:params.nAer
        # Create Aerosol size distribution for each aerosol species
        size_distribution = LogNormal(log(params.μ[iaer]), log(params.σ[iaer]))
        # Create a univariate aerosol distribution
        aerosol = make_univariate_aerosol(size_distribution, params.r_max, params.nquad_radius, params.nᵣ[iaer], params.nᵢ[iaer]) #Suniti: why is the refractive index needed here?
        # Create the aerosol extinction cross-section at the reference wavelength:
        #@show params
        mie_model      = make_mie_model(params.decomp_type, aerosol, params.λ_ref, params.polarization_type, truncation_type)       
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)
        # Loop over bands
        for ib=1:nBands
            # Create the aerosols:
            mie_model      = make_mie_model(params.decomp_type, aerosol, params.λ_band[ib], params.polarization_type, truncation_type)
            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, params.float_type);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            @show iaer,ib
            aerosol_optics[ib][iaer] = Scattering.truncate_phase(truncation_type, aerosol_optics_raw; reportFit=false)
            # Compute nAer aerosol optical thickness profiles
            τAer[ib][iaer,:] = params.τAer_ref[iaer] * (aerosol_optics[ib][iaer].k/k_ref) * vSmartMOM.getAerosolLayerOptProp(1.0, params.p₀[iaer], params.σp[iaer], profile.p)
            # @show τAer, sum(τAer)
            # Can be done with arbitrary length later:
            # τAer = FT(0.2) * τAer_1; # [τAer_1 τAer_2]
            #ϖAer = params.float_type[aerosol_optics[iaer,ib].ω̃]; 
            # fᵗ   = FT[aerosol_optics_trunc_aero1.fᵗ]; # [aerosol_optics_trunc_aero1.fᵗ aerosol_optics_trunc_aero2.fᵗ];
            # aerosol_optics = [aerosol_optics_trunc_aero1] #
        end 
    end
    
    # Just return band 1 for now:
    iBand = 1;
    return vSmartMOM_Model(params, aerosol_optics[iBand],  greek_rayleigh, Nquad, qp_μ, wt_μ, τ_abs[iBand], τRayl[iBand],τAer[iBand],  obs_geom,profile)

end