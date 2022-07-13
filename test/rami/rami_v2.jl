using JSON
using vSmartMOM
using DelimitedFiles
using Distributions
                                  
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


## 

# Read in the standard atmosphere

# atmosphere ascii file location (MAYBE add as artifact)
rami_atmos_fn = "test/rami/RAMI4ATM_AFGLUSstandard_ap_v1.0.txt"
US_std_atmos = readdlm(rami_atmos_fn, ' ')


## 

# Starts with HOM00_LAM, HOM00_BLA, or HOM00_WHI for lambertian surfaces
# Starts with HOM00_RPV, HOM00_RLI for Anisotropic surfaces
# Starts with HOM45_LAM, HOM25_LAM, or HOM35_LAM for homogeneous discrete with isotropic background

function produce_rami_results(experiment_name::String;
                              rami_json::String = "test/rami/RAMI4ATM_experiments_v1.0.json")

    # Get rami scenario from experiment name
    all_scenarios = JSON.parsefile(rami_json);
    filter!(x->x["name"] == experiment_name, all_scenarios);
    @assert length(all_scenarios) ≤ 1 "Multiple matching experiment names!"
    @assert length(all_scenarios) ≥ 1 "Experiment name not found in JSON input"
    curr_scenario = all_scenarios[1]["observations"][1];
    band = curr_scenario["measures"][1]["bands"]

    # There are 7 rami categories to be mindful of: 
    # atmosphere, name, canopy, illumination, measures, time, surface

    rami_atmosphere   = curr_scenario["atmosphere"]
    rami_name         = curr_scenario["name"]
    rami_canopy       = curr_scenario["canopy"]
    rami_illumination = curr_scenario["illumination"]
    rami_measures     = curr_scenario["measures"]
    rami_time         = curr_scenario["time"]
    rami_surface      = curr_scenario["surface"]

    conc              = rami_atmosphere["concentrations"]
    atm_type          = rami_atmosphere["atmosphere_type"]
    @show atm_type
    
    # Get a default set of Rayleigh params
    default_params = "test/rami/" * scenario_yaml_I[atm_type]
    @show default_params
    params = vSmartMOM.parameters_from_yaml(default_params)

    
    # ########################################
    # Change params fields to match rami input
    # ########################################

    
    # ########################################
    # 1. Atmosphere
    # ########################################

    # aerosols 
    add_aerosols!(rami_atmosphere, params)
    @show params.scattering_params.rt_aerosols[1].τ_ref

    # TODO: Fix q from vmr (needs to be that way)
    params.q = [1.3683608600937127e-7, 1.6171537681996363e-7, 1.92814491391246e-7, 2.301334304285797e-7, 2.923316659193965e-7, 4.322777129669201e-7, 6.779608309202629e-7, 1.0511505175455402e-6, 1.5176378638869645e-6, 1.9685758834495292e-6, 2.3946347910935506e-6, 2.783374862836421e-6, 3.0632678282313282e-6, 3.212544115378091e-6, 3.259192960666205e-6, 3.259192960666205e-6, 3.228093730180227e-6, 3.165895272735097e-6, 3.103696819992401e-6, 3.063267828231329e-6, 3.022838838457035e-6, 2.9699701625192205e-6, 2.8953320375674745e-6, 2.8020343909002165e-6, 2.7149565968900444e-6, 2.643428415848055e-6, 2.5719002410250296e-6, 2.5034819926657405e-6, 2.450613350103426e-6, 2.4101843904364436e-6, 2.388414951438685e-6, 2.388414951438685e-6, 2.425733990073239e-6, 2.783374862836421e-6, 3.3991395124012653e-6, 5.233996504028426e-6, 9.329785263484815e-6, 1.7166886651582844e-5, 3.299681516351052e-5, 7.090902170019672e-5, 0.00016328651918899218, 0.0002920724597408198, 0.00046568540880824054, 0.0007233721407928284, 0.0011078736971087291, 0.0016623702032229144, 0.0024324309812018903, 0.003334347912200177, 0.004309152659531596]

    # ########################################
    # 2. Name
    # ########################################

    # ########################################
    # 3. Illumination
    # ########################################
    params.sza = rami_illumination["sza"]["value"]

    # ########################################
    # 4. Measures
    # ########################################
    # Spectral bands
    iBand        = sentinel_band_to_index[rami_measures[1]["bands"][1]]
    params.spec_bands = [params.spec_bands[iBand]]
    @show params.spec_bands

    # Viewing zenith angles (#TODO: this needs to be looped, can use hardcoded values!)
    @show size(rami_measures)
    for ra in rami_measures
        @show ra["measure_type"], ra["height"]
    end
    vza_start = 1.0  # rami_measures[1]["vza_start"]["value"]
    vza_end   = 75.0 # rami_measures[1]["vza_end"]["value"]
    vza_step  = 2.0  # rami_measures[1]["vza_step"]["value"]
    vzas = collect(vza_start:vza_step:vza_end)
    params.vza  = [vzas; vzas; vzas; vzas]


    # Viewing azimuth (TODO: Check definitions!)
    # vaa = rami_measures[1]["delta_vaa"]["value"]
    params.vaz  = [repeat([0.0], length(vzas)); repeat([90.0], length(vzas));repeat([180.0], length(vzas));repeat([270.0], length(vzas)) ]
    @show params.vaz
    
    #@show params.spec_bands
    # ########################################
    # 5. Time
    # ########################################

    # ########################################
    # 6. Surface
    # ########################################
    @assert rami_surface["name"] in ["WHI", "BLA", "LAM"] && startswith(curr_scenario["name"], "HOM00") "Currently only supporting Lambertian (HOM00) surfaces"
    params.brdf = [vSmartMOM.LambertianSurfaceScalar(rami_surface["surface_parameters"]["reflectance"])]
    params.brdf = [vSmartMOM.LambertianSurfaceScalar(0.0)]
    @show 
    model = model_from_parameters(params)

    return rt_run(model), model
end

R, models = produce_rami_results("HOM00_WHI_SD2S_M02_z30a000")

##

#=
experiment_name = "HOM00_WHI_SD2S_M02_z30a000"
rami_json = "/home/rjeyaram/vSmartMOM/test/rami/RAMI4ATM_experiments_v1.0.json"
all_scenarios = JSON.parsefile(rami_json)
filter(x->x["name"] == experiment_name, all_scenarios)
rayleigh_scenario = all_scenarios[1]["observations"][1]
params = vSmartMOM.parameters_from_yaml("/home/rjeyaram/vSmartMOM/test/test_parameters/PureRayleighParameters.yaml")
=#

# Rayleigh: HOM00_WHI_S00S_M02_z30a000
# Rayleigh + Scattering: HOM00_WHI_SD2S_M02_z30a000
# Absorbing: HOM00_WHI_A00S_M02_z30a000
# Scattering + Absorbing: HOM00_WHI_AD2S_M02_z30a000