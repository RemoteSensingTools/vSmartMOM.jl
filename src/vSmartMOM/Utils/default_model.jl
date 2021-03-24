
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

    μ            = FT(1.3)
    σ            = FT(2.0)
    r_max        = FT(30.0)
    nquad_radius = 2500
    nᵣ           = FT(1.3)
    nᵢ           = FT(0.00000001)

    # Aerosol vertical distribution profiles
    p₀          = FT(90000) # Pressure peak [Pa]
    σp          = FT(5000.) # Pressure peak width [Pa]

    file = "/net/fluo/data1/ftp/XYZT_ESE156/Data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4" 
    timeIndex = 2

    lat = 34.1377;
    lon = -118.1253;

    profile_reduction_n = 20; 

    decomp_type = NAI2()

    quadrature_type = RadauQuad() 

    architecture = default_architecture;

    grid_start = FT(1e7 / 774)
    grid_end = FT(1e7 / 757)
    grid_n = 10000

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

function parameters_from_yaml(file_path)

    params_dict = YAML.load_file(file_path)
    broadening_function = eval(Meta.parse(params_dict["absorption"]["broadening"]))
    CEF = eval(Meta.parse(params_dict["absorption"]["CEF"]))

    polarization_type = eval(Meta.parse(params_dict["scattering"]["polarization_type"]))
    decomp_type = eval(Meta.parse(params_dict["scattering"]["decomp_type"]))

    quadrature_type = eval(Meta.parse(params_dict["vSmartMOM"]["quadrature_type"]))

    architecture = eval(Meta.parse(params_dict["vSmartMOM"]["architecture"]))

    FT = eval(Meta.parse(params_dict["vSmartMOM"]["float_type"]))

    return vSmartMOM_Parameters(FT,
                                FT(params_dict["scattering"]["λ"]),
                                FT(params_dict["scattering"]["depol"]),
                                params_dict["truncation_type"]["l_trunc"],
                                FT(params_dict["truncation_type"]["Δ_angle"]),
                                params_dict["vSmartMOM"]["max_m"],
                                polarization_type,
                                FT(params_dict["geometry"]["obs_alt"]),
                                FT(params_dict["geometry"]["sza"]),
                                convert.(FT, params_dict["geometry"]["vza"]),
                                convert.(FT, params_dict["geometry"]["vaz"]),
                                FT(params_dict["scattering"]["aerosols"][1]["μ"]),
                                FT(params_dict["scattering"]["aerosols"][1]["σ"]),
                                FT(params_dict["scattering"]["aerosols"][1]["r_max"]),
                                params_dict["scattering"]["aerosols"][1]["nquad_radius"],
                                FT(params_dict["scattering"]["nᵣ"]),
                                FT(params_dict["scattering"]["nᵢ"]),
                                FT(params_dict["scattering"]["p₀"]),
                                FT(params_dict["scattering"]["σp"]),
                                params_dict["atmospheric_profile"]["file"],
                                params_dict["atmospheric_profile"]["time_index"],
                                params_dict["atmospheric_profile"]["lat"],
                                params_dict["atmospheric_profile"]["lon"],
                                params_dict["atmospheric_profile"]["profile_reduction"],
                                FT(params_dict["absorption"]["ν_min"]),
                                FT(params_dict["absorption"]["ν_max"]),
                                FT(params_dict["absorption"]["ν_step"]),
                                broadening_function,
                                params_dict["absorption"]["wing_cutoff"],
                                CEF,
                                FT(params_dict["absorption"]["vmr"]),
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
    τAer = params.float_type(0.2) * vSmartMOM.getAerosolLayerOptProp(1.0, params.p₀, params.σp, profile.p_levels)

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