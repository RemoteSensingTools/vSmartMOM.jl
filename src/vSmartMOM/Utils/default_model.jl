
function default_parameters(FT::DataType=Float64)

    # Note, we will need many parameters in Arrays, as we will have `nAer` aerosol types (1 in default)
    # as well as nBands (we can model multiple spectral bands)
    λ_band = FT[0.770] 
    λ_ref = FT(0.770) 
    depol = FT(0.0)
    l_trunc = 30 
    Δ_angle = FT(0.0)
    polarization_type = Stokes_IQU{FT}()

    # Define Surface here (per band, we just use 1 for now!)
    BRDF = LambertianSurfaceScalar(0.35)
    BRDF_per_band = [BRDF];

    
    #vza = FT[60., 45., 30., 15., 0., 15., 30., 45., 60.]
    #vaz = FT[180., 180., 180., 180., 0., 0., 0., 0., 0.]
    vza = FT[10.]
    vaz = FT[0.]
    sza = FT(60.)

    obs_alt = FT(1000.0)
    nAer = 1 #Number of aerosol species

    # These are arrays by aerosol type!
    τAer_ref     = FT[0.12]  # AOTs at reference wavelength
    μ            = FT[1.3]  # characteristic radius [μm]
    σ            = FT[2.0]  # characteristic width
    nᵣ           = FT[1.3]  # TODO: make this a 2D matrix were each species and band has an independent entry
    nᵢ           = FT[0.00000001] # TODO: same as above
    # Aerosol vertical distribution profiles (per aerosol type)
    p₀          = FT[90000] # Pressure peak [Pa]
    σp          = FT[5000.] # Pressure peak width [Pa]

    # More global setting!
    r_max        = FT(30.0) # maximum radius in distribution [μm] #baseline setting
    nquad_radius = 2500     # baseline setting for quadrature points in size distribution

    

    # Note: We should change the default profile to an ASCII file and have a simple ascii reader...
    #file = "/Users/cfranken/data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
    file = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4"
    timeIndex = 2

    lat = FT(34.1377);
    lon = FT(-118.1253);

    profile_reduction_n = 20; 

    decomp_type = NAI2()

    #quadrature_type = RadauQuad() 
    quadrature_type = GaussQuadFullSphere()

    architecture = default_architecture;

    SFI = 1 #Suniti: 0:= DNI, 1:= SFI
    spec_grid_start = FT[(1e7 / 777)]
    spec_grid_end   = FT[(1e7 / 757)]
    spec_grid_n     = Integer[15000]

    broadening_function = Voigt()
    wing_cutoff = 100
    CEF=HumlicekWeidemann32SDErrorFunction()
    vmr = FT(0.21)
    #vmr = FT(0.00)

    # This is only a hard-coded quick fix, eventually we need to compute after how many Fourier 
    # components it converges and then stop (and just set a convergence criterion)
    #max_m = l_trunc #temp
    max_m = 3 
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
                                BRDF_per_band,
                                decomp_type,
                                quadrature_type,
                                #BRDF_type,
                                #a_surf,
                                architecture)

end

function parameters_from_yaml(file_path)

    params_dict = YAML.load_file(file_path)
    broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
    CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))

    polarization_type = eval(Meta.parse(params_dict["scattering"]["polarization_type"]))
    decomp_type = eval(Meta.parse(params_dict["scattering"]["decomp_type"]))

    quadrature_type = eval(Meta.parse(params_dict["vSmartMOM"]["quadrature_type"]))

    architecture = eval(Meta.parse(params_dict["vSmartMOM"]["architecture"]))

    FT = eval(Meta.parse(params_dict["vSmartMOM"]["float_type"]))

    BRDF_per_band = map(x -> eval(Meta.parse(x)), params_dict["surface"]["BRDFs"]) 

    lengths = convert.(Integer, map(x -> length(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"]))

    profile_path = params_dict["atmospheric_profile"]["file"]

    if endswith(profile_path, ".yaml")
        time_index = nothing
        lat = nothing
        lon = nothing
    else
        time_index = params_dict["atmospheric_profile"]["time_index"]
        lat = params_dict["atmospheric_profile"]["lat"]
        lon = params_dict["atmospheric_profile"]["lon"]
    end

    return vSmartMOM_Parameters(FT,
                                convert.(FT, params_dict["scattering"]["λ"]),
                                FT(params_dict["scattering"]["λ_ref"]),
                                FT(params_dict["scattering"]["depol"]),
                                params_dict["truncation_type"]["l_trunc"],
                                FT(params_dict["truncation_type"]["Δ_angle"]),
                                params_dict["vSmartMOM"]["max_m"],
                                polarization_type,
                                FT(params_dict["geometry"]["obs_alt"]),
                                FT(params_dict["geometry"]["sza"]),
                                convert.(FT, params_dict["geometry"]["vza"]),
                                convert.(FT, params_dict["geometry"]["vaz"]),
                                length(params_dict["scattering"]["aerosols"]), 
                                convert.(FT, map(x -> x["τ_ref"], params_dict["scattering"]["aerosols"])),
                                # FT(params_dict["scattering"]["aerosols"][1]["τ_ref"]),
                                convert.(FT, map(x -> x["μ"], params_dict["scattering"]["aerosols"])),
                                # FT(params_dict["scattering"]["aerosols"][1]["μ"]),
                                convert.(FT, map(x -> x["σ"], params_dict["scattering"]["aerosols"])),
                                # FT(params_dict["scattering"]["aerosols"][1]["σ"]),
                                FT(params_dict["scattering"]["r_max"]),
                                params_dict["scattering"]["nquad_radius"],
                                convert.(FT, map(x -> x["nᵣ"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["nᵢ"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["p₀"], params_dict["scattering"]["aerosols"])),
                                convert.(FT, map(x -> x["σp"], params_dict["scattering"]["aerosols"])),
                                profile_path,
                                time_index,
                                lat,
                                lon,
                                params_dict["atmospheric_profile"]["profile_reduction"],
                                convert.(FT, map(x -> first(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"])),
                                convert.(FT, map(x -> last(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"])),
                                convert(Array{Integer}, convert.(Integer, map(x -> length(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"]))),
                                # Integer[10000],
                                # reshape(lengths, length(lengths), 1),
                                # convert.(Int, map(x -> length(collect(eval(Meta.parse(x)))), params_dict["absorption"]["spec_bands"])),
                                # FT(params_dict["absorption"]["ν_max"]),
                                # FT(params_dict["absorption"]["ν_step"]),
                                broadening_function,
                                params_dict["absorption"]["wing_cutoff"],
                                CEF,
                                FT(params_dict["absorption"]["vmr"]),
                                BRDF_per_band,
                                decomp_type,
                                quadrature_type,
                                architecture)
end

function parameters_from_json(file_path, FT::DataType=Float64)

    JSON_obj = JSON.parsefile(file_path)
    # validate_JSON_parameters(JSON_obj)

    polarization_type = JSON_obj["polarization_type"] == "Stokes_I" ? Stokes_I{FT}() : 
                        JSON_obj["polarization_type"] == "Stokes_IQU" ? Stokes_IQU{FT}() : 
                        JSON_obj["polarization_type"] == "Stokes_IQUV" ? Stokes_IQUV{FT}() : Stokes_IQU{FT}()

    decomp_type = JSON_obj["decomp_type"] == "NAI2" ? NAI2() : 
                  JSON_obj["decomp_type"] == "PCW" ? PCW() : NAI2()

    quadrature_type = JSON_obj["quadrature_type"] == "RadauQuad" ? RadauQuad() : 
                      JSON_obj["quadrature_type"] == "GaussQuadHemisphere" ? GaussQuadHemisphere() : 
                      JSON_obj["quadrature_type"] == "GaussQuadFullSphere" ? GaussQuadFullSphere() : RadauQuad()

    architecture = JSON_obj["decomp_type"] == "Default" ? default_architecture : 
                   JSON_obj["decomp_type"] == "GPU" ? GPU : 
                   JSON_obj["decomp_type"] == "CPU" ? CPU : default_architecture

    broadening_function = JSON_obj["broadening_function"] == "Voigt" ? Voigt() : 
                          JSON_obj["broadening_function"] == "Doppler" ? Doppler() : 
                          JSON_obj["broadening_function"] == "Lorentz" ? Lorentz() : Voigt()

    CEF = JSON_obj["CEF"] == "HumlicekWeidemann32SDErrorFunction" ? HumlicekWeidemann32SDErrorFunction() : HumlicekWeidemann32SDErrorFunction()

    return vSmartMOM_Parameters(FT,
                                FT(JSON_obj["λ"]),
                                FT(JSON_obj["depol"]),
                                JSON_obj["l_trunc"],
                                FT(JSON_obj["Δ_angle"]),
                                JSON_obj["max_m"],
                                polarization_type,
                                FT(JSON_obj["obs_alt"]),
                                FT(JSON_obj["sza"]),
                                convert.(FT, JSON_obj["vza"]),
                                convert.(FT, JSON_obj["vaz"]),
                                FT(JSON_obj["μ"]),
                                FT(JSON_obj["σ"]),
                                FT(JSON_obj["r_max"]),
                                JSON_obj["nquad_radius"],
                                FT(JSON_obj["nᵣ"]),
                                FT(JSON_obj["nᵢ"]),
                                FT(JSON_obj["p₀"]),
                                FT(JSON_obj["σp"]),
                                JSON_obj["file"],
                                JSON_obj["timeIndex"],
                                JSON_obj["lat"],
                                JSON_obj["lon"],
                                JSON_obj["profile_reduction_n"],
                                FT(JSON_obj["grid_start"]),
                                FT(JSON_obj["grid_end"]),
                                JSON_obj["grid_n"],
                                broadening_function,
                                JSON_obj["wing_cutoff"],
                                CEF,
                                FT(JSON_obj["vmr"]),
                                decomp_type,
                                quadrature_type,
                                architecture)

end


function model_from_parameters(params::vSmartMOM_Parameters)

    
    @unpack λ_band = params
    nBands = length(λ_band)

    # Create observation geometry:
    obs_geom = ObsGeometry(params.obs_alt, params.sza, params.vza, params.vaz)

    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    #Nquad, qp_μ, wt_μ = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom);
    quadPoints = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))
    # Read profile (and generate dry/wet VCDs per layer)
    if isnothing(params.timeIndex)
        profile_hr = vSmartMOM.read_atmos_profile(params.file);
    else
        profile_hr = vSmartMOM.read_atmos_profile(params.file, params.lat, params.lon, params.timeIndex);
    end
    
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
        compute_absorption_profile!(τ_abs[ib], absorption_model, range(params.spec_grid_start[ib], params.spec_grid_end[ib], length=params.spec_grid_n[ib]), profile);
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (params.nAer)) for i=1:nBands];

    FT2 = typeof(params.τAer_ref[1])
    # τAer[iBand][iAer,iProfile]
    τAer = [zeros(FT2, params.nAer, length(profile.p)) for i=1:nBands];

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
            aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            @show iaer,ib
            aerosol_optics[ib][iaer] = Scattering.truncate_phase(truncation_type, aerosol_optics_raw; reportFit=false)
            # Compute nAer aerosol optical thickness profiles
            τAer[ib][iaer,:] = params.τAer_ref[iaer] * (aerosol_optics[ib][iaer].k/k_ref) * vSmartMOM.getAerosolLayerOptProp(1.0, params.p₀[iaer], params.σp[iaer], profile.p)
            
        
        end 
    end
    
    # Just return band 1 for now:
    iBand = 1;
    return vSmartMOM_Model(params, aerosol_optics[iBand],  greek_rayleigh, quadPoints, τ_abs[iBand], τRayl[iBand],τAer[iBand], obs_geom, profile, params.brdf[iBand])#, BRDF_type, a_surf[iBand])

end