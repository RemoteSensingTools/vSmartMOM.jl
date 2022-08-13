using Interpolations
using DelimitedFiles
using Printf


# Read Spectral Response Functions:
sr = readdlm("test/rami/Sentinel2A_ILS.dat")
n_desert      = readdlm("test/rami/refractive_aero_desert.txt")
n_continental = readdlm("test/rami/refractive_aero_continental.txt")

const sr_2 = sr[:,3]
const sr_3 = sr[:,4]
const sr_4 = sr[:,5]
const sr_8A = sr[:,10]
const sr_11 = sr[:,13]
const sr_12 = sr[:,14]
const wl_sr = sr[:,1]

# Some LUT for convenience:
const sentinel_band_to_index = Dict([("2" , 1),    # Blue
                                     ("3" , 2),    # Green
                                     ("4" , 3),    # Red
                                     ("8A", 4),    # Narrow NIR 
                                     ("11", 5),  # SWIR
                                     ("12", 6)]) # SWIR

const scenario_yaml_I = Dict([("AtmosphereType.RAYLEIGH"           , "RamiNoGasI.yaml"),    
                            ("AtmosphereType.ABSORBING"            , "RamiGasI.yaml"),      
                            ("AtmosphereType.SCATTERING_ABSORBING" , "RamiGasI.yaml"),      
                            ("AtmosphereType.AEROSOLS"             , "RamiNoGasI.yaml"),     
                            ("AtmosphereType.SCATTERING_AEROSOLS"  , "RamiNoGasI.yaml"),    
                            ("AtmosphereType.ABSORBING_AEROSOLS"   , "RamiGasI.yaml"),
                            ("AtmosphereType.COMPLETE"             , "RamiGasI.yaml" )]) # SWIR

const noRayleighType   = ("AtmosphereType.ABSORBING", "AtmosphereType.AEROSOLS", "AtmosphereType.ABSORBING_AEROSOLS")
const noAbsorptionType = ("AtmosphereType.RAYLEIGH", "AtmosphereType.AEROSOLS", "AtmosphereType.SCATTERING_AEROSOLS")


const sentinel_band_to_LUT = Dict([("2" , 2),    # Blue
                            ("3" , 3),    # Green
                            ("4" , 4),    # Red
                            ("8A", 9),    # Narrow NIR 
                            ("11", 12),  # SWIR
                            ("12", 13)]) # SWIR


# Add aerosols (if required by RAMI, otherwise do nothing)
function add_aerosols!(rami_atmosphere, params, band, n_desert, n_continental)# TODO We still need to check the two fractions of coarse and fine! Right now, only one is
    if length(rami_atmosphere["aerosols"]) > 0# 
        if startswith(rami_atmosphere["aerosols"][1]["name"], "D")
            @show "Desert Aerosol"
            μ_fine   = 0.0478666
            σ_fine   = 1.87411
            μ_coarse = 0.604127
            σ_coarse = 1.75172
            # Scaled so that sum is 1
            n_coarse = 0.0033219635
            
            # These refractive indices depend on the band. 
            # TODO: We need to pull those from the table actually, as they are band dependent 
            # (in vSmartMOM, we so far only have one ni/nr per type).
            # See https://rami-benchmark.jrc.ec.europa.eu/_www/RAMI4ATM/down/RAMI4ATM_aerosols_v1.0/refractive_index/continental.txt
            # Need to generalize
            
            # TODO: Make it wavelength dependent (and reference scenario fixed)
            #nᵣ = 1.477538814814815
            #nᵢ = 0.004342592592592592
            @show band
            iN = sentinel_band_to_LUT[band]
            nᵣ = n_desert[iN,2]
            #nᵣ = 1.4220380000000001 
            nᵢ = n_desert[iN,3]
            n_ref = n_desert[3,2] - im*n_desert[3,3]
            @show nᵣ, nᵢ

        elseif startswith(rami_atmosphere["aerosols"][1]["name"], "C")
            @show "Continental Aerosol"
            #TODO: Change, not just one refractice index
            μ_fine   = 0.0807989
            σ_fine   = 1.50180
            μ_coarse = 0.682651
            σ_coarse = 2.10400
            # Scaled so that sum is 1
            n_coarse   = 0.00046374026257

            # TODO: Make it wavelength dependent (and reference scenario fixed)
            @show band[1]
            iN = sentinel_band_to_LUT[band]
            nᵣ = n_continental[iN,2]
            #nᵣ = 1.4220380000000001 
            nᵢ = n_continental[iN,3]
            n_ref = n_continental[3,2] - im*n_continental[3,3]
            @show nᵣ, nᵢ
           

        else
            println("weird aerosol here")
        end

        # TODO: Suniti, please double-check

        # size_distribution = LogNormal(log(μ), log(σ))
        size_distribution = MixtureModel(LogNormal,
                            [(log(μ_fine), log(σ_fine)), (log(μ_coarse), log(σ_coarse))],
                            [1-n_coarse, n_coarse])

        # Create the aerosol(s)
        aero = vSmartMOM.CoreRT.Aerosol(size_distribution, nᵣ, nᵢ)

        # Need to do the pressure peak more carefully
        # TODO: We need to just put it into 1-2 layers...
        RT_aerosol = vSmartMOM.CoreRT.RT_Aerosol(aero, rami_atmosphere["aerosols"][1]["tau_550"], Uniform(795.0,1013.0))

        # Assemble scattering parameters
        scattering_params = vSmartMOM.CoreRT.ScatteringParameters([RT_aerosol], 20.0, 1000, 0.550, n_ref, vSmartMOM.Scattering.NAI2())

        params.scattering_params = scattering_params
        @show params.scattering_params.rt_aerosols[1].τ_ref
    end
end

# Scale gas concentrations (if required by RAMI, otherwise do nothing)
function scale_gases!(rami_atmosphere, params)
    conc = rami_atmosphere["concentrations"]
    if !isempty(conc)
        # Reference standrad values kg/m2
        ref_h2o = 14.274
        ref_o3  = 0.746e-2
        h2o = conc["H2O"]["value"]
        o3  = conc["O3"]["value"]
        params.absorption_params.vmr["O3"]  .*= h2o/ref_h2o
        params.absorption_params.vmr["H2O"] .*= o3/ref_o3
        @show h2o/ref_h2o, o3/ref_o3
    end
end

# Convolve with Sentinel response function (input is in wavenumbers, convolution in nm space)
function convolve_2_sentinel(ν, BRF, band)
    ex = Meta.parse("sentinelILS = sr_" * band)
    eval(ex)
    srr = LinearInterpolation(wl_sr, sentinelILS);
    d1,d2,d3 = size(BRF)

    wl_in = 1e7./ν;
    weight = srr.(wl_in)
    weight = weight/sum(weight)
    BRF_out = zeros(d1)
    for i in eachindex(BRF_out)
        BRF_out[i] = weight' * BRF[i,1,:]
    end
    return BRF_out
end

function getParams(scenario)
    # Get Band (always just one!)
    band              = scenario["measures"][1]["bands"][1]
    rami_atmosphere   = scenario["atmosphere"]
    atm_type          = rami_atmosphere["atmosphere_type"]
    if atm_type ∈ noAbsorptionType
        yamlFile = "test/rami/RamiNoGasI.yaml"
        @info "Using " * yamlFile * " as input"
        params   = vSmartMOM.parameters_from_yaml(yamlFile)
        # just select the right band:
        #@show band, sentinel_band_to_index[band]
        params.spec_bands = [params.spec_bands[sentinel_band_to_index[band]]]
        @info "Mean wl is " * string(mean(1e7./params.spec_bands[1]))
    else
        # For gases, we have a separate YAML per band (easier for LUT)
        yamlFile = "test/rami/RamiGasI_" * band * ".yaml"
        @info "Using " * yamlFile * " as input"
        params   = vSmartMOM.parameters_from_yaml(yamlFile)
    end
    params.l_trunc = 40
    return params
end

# Determin max amount of fourier moments (to be extended):
function assignMoment!(params,atm_type)

    if atm_type == "AtmosphereType.RAYLEIGH"
        params.max_m = 3
    end
end

function setGeometry!(scenario, params)
    # q precomputed from vmr (needs to be that way for now)
    params.q = [1.3683608600937127e-7, 1.6171537681996363e-7, 1.92814491391246e-7, 2.301334304285797e-7, 2.923316659193965e-7, 4.322777129669201e-7, 6.779608309202629e-7, 1.0511505175455402e-6, 1.5176378638869645e-6, 1.9685758834495292e-6, 2.3946347910935506e-6, 2.783374862836421e-6, 3.0632678282313282e-6, 3.212544115378091e-6, 3.259192960666205e-6, 3.259192960666205e-6, 3.228093730180227e-6, 3.165895272735097e-6, 3.103696819992401e-6, 3.063267828231329e-6, 3.022838838457035e-6, 2.9699701625192205e-6, 2.8953320375674745e-6, 2.8020343909002165e-6, 2.7149565968900444e-6, 2.643428415848055e-6, 2.5719002410250296e-6, 2.5034819926657405e-6, 2.450613350103426e-6, 2.4101843904364436e-6, 2.388414951438685e-6, 2.388414951438685e-6, 2.425733990073239e-6, 2.783374862836421e-6, 3.3991395124012653e-6, 5.233996504028426e-6, 9.329785263484815e-6, 1.7166886651582844e-5, 3.299681516351052e-5, 7.090902170019672e-5, 0.00016328651918899218, 0.0002920724597408198, 0.00046568540880824054, 0.0007233721407928284, 0.0011078736971087291, 0.0016623702032229144, 0.0024324309812018903, 0.003334347912200177, 0.004309152659531596]
    rami_illumination = scenario["illumination"]

    # ########################################
    params.sza = rami_illumination["sza"]["value"]

    # Viewing zenith angles (#TODO: this needs to be looped, can use hardcoded values!)
    #@show size(rami_measures)
    #for ra in rami_measures
    #    @show ra["measure_type"], ra["height"]
    #end

    vza_start = 1.0  # rami_measures[1]["vza_start"]["value"]
    vza_end   = 75.0 # rami_measures[1]["vza_end"]["value"]
    vza_step  = 2.0  # rami_measures[1]["vza_step"]["value"]
    vzas = collect(vza_start:vza_step:vza_end)
    params.vza  = [reverse(vzas); vzas; reverse(vzas); vzas]
    params.vaz  = [repeat([0.0], length(vzas)); repeat([180.0], length(vzas));repeat([-90.0], length(vzas));repeat([90.0], length(vzas)) ]
end

# Set surface type (Lambertian for now!)
function setSurface!(scenario, params)
    rami_surface      = scenario["surface"]
    
    if rami_surface["name"] ∈ ["WHI", "BLA", "LAM"]
        params.brdf = [vSmartMOM.LambertianSurfaceScalar(rami_surface["surface_parameters"]["reflectance"])]
    elseif rami_surface["name"] == "RPV"
        p = rami_surface["surface_parameters"]
        params.brdf = [vSmartMOM.CoreRT.rpvSurfaceScalar(p["rho_0"][1],p["rho_c"][1], p["k"][1],p["theta"][1] )]
    else
        @assert rami_surface["name"] in ["WHI", "BLA", "LAM", "RPV"] && startswith(scenario["name"], "HOM00") "Currently only supporting Lambertian (HOM00) and RPV surfaces"
    end
    @info "Surface type: " * rami_surface["name"]
end
 

function save_toa_results(BRF, model, name, folder)
    fileNamePP = folder * name * "-brfpp_" * "vSmartMOM.mes" 
    open(fileNamePP, "w") do io
        @printf io "%4d    %4d   %.6f \n" 76	5	-1.000000
        for i=1:76
            @printf io "%.6f    %.6f   %.6f     %.6f    %.6f \n" deg2rad(model.obs_geom.sza) deg2rad(model.obs_geom.vza[i]) deg2rad(model.obs_geom.vaz[i]) BRF[i] -1
        end
    end;

    fileNameOP = folder * name * "-brfop_" * "vSmartMOM.mes" 
    open(fileNameOP, "w") do io
        @printf io "%4d    %4d   %.6f \n" 76	5	-1.000000
        for i=77:152
            @printf io "%.6f    %.6f   %.6f     %.6f    %.6f \n" deg2rad(model.obs_geom.sza) deg2rad(model.obs_geom.vza[i]) deg2rad(model.obs_geom.vaz[i]) BRF[i] -1
        end
    end;
end



# Get q from vmr:
# q = -1/((dry_mass - dry_mass./vmr_h2o)/wet_mass - 1)