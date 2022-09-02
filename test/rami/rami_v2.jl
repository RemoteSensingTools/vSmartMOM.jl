using JSON
using vSmartMOM
using DelimitedFiles
using Distributions

include("test/rami/rami_tools.jl")

# Read in the standard atmosphere

# atmosphere ascii file location (MAYBE add as artifact)
# rami_atmos_fn = "test/rami/RAMI4ATM_AFGLUSstandard_ap_v1.0.txt"
# US_std_atmos = readdlm(rami_atmos_fn, ' ')

# Starts with HOM00_LAM, HOM00_BLA, or HOM00_WHI for lambertian surfaces
# Starts with HOM00_RPV, HOM00_RLI for Anisotropic surfaces
# Starts with HOM45_LAM, HOM25_LAM, or HOM35_LAM for homogeneous discrete with isotropic background

function produce_rami_results(experiment_name::String;
                        all_scenarios = all_scenarios,
                        n_desert=n_desert,n_continental=n_continental, sr=sr )

    # Get rami scenario from experiment name
    # all_scenarios = JSON.parsefile(rami_json);
    filt_scenario = filter(x->x["name"] == experiment_name, all_scenarios);
    @assert length(filt_scenario) ≤ 1 "Multiple matching experiment names!"
    @assert length(filt_scenario) ≥ 1 "Experiment name not found in JSON input"
    curr_scenario = filt_scenario[1]["observations"][1];
    #name = curr_scenario["name"]
    # Band to use    
    band = curr_scenario["measures"][1]["bands"][1]

    rami_atmosphere   = curr_scenario["atmosphere"]
    atm_type = rami_atmosphere["atmosphere_type"]
    @info "Atmosphere type is " * atm_type 
    
    params = getParams(curr_scenario)
    
    # Add Aerosols: 
    add_aerosols!(rami_atmosphere, params, band, n_desert, n_continental)
    
    # Scale trace gases:
    scale_gases!(rami_atmosphere, params)

    # Set viewing geometries:
    setGeometry!(curr_scenario, params)

    # Set surface parameters:
    setSurface!(curr_scenario, params)

    # Create vSmartMOM model:
    model = model_from_parameters(params)

    # Turning off Rayleigh scattering for non Rayleigh types (has to be done in model, not params):
    if atm_type ∈ noRayleighType
        @info "Turning off Rayleigh ", atm_type
        model.τ_rayl[1] .= 0.0
    end
    #  model.τ_rayl[1] .= 0.00000000001
    ########################################################

    # Run model (can think about including the BOA and hemispheric data here as well)
    #R = rt_run(model)
    R, _, _, _, hdr, bhr_uw, bhr_dw = rt_run(model)
    # Convolve results:
    BRF = convolve_2_sentinel(model.params.spec_bands[1], R[1], band)
    save_toa_results(BRF, model, experiment_name, "/home/cfranken/rami/")
    return BRF, R, hdr, bhr_uw, bhr_dw, model
    #return BRF, R, model

    # Can do postprocessing here (saving data, convolution, etc)
end

# Get all scenarios
rami_json = "test/rami/RAMI4ATM_experiments_v1.0.json";
all_scenarios = JSON.parsefile(rami_json);

# R, models = produce_rami_results("HOM00_WHI_SD2S_M02_z30a000")
#R, models = produce_rami_results("HOM00_LAM_ED2D_M12_z30a000")
#R, models = produce_rami_results("HOM00_WHI_A00S_M02_z30a000")
#BRF, R, model = produce_rami_results("HOM00_WHI_A00S_M12_z30a000")
#white lambertian
BRF1, R1, hdr1, bhr_uw1, bhr_dw1, model = produce_rami_results("HOM00_LAM_S00S_M03_z30a000")#("HOM00_WHI_EC6L_M03_z30a000")
#BRF1, R, model = produce_rami_results("HOM00_WHI_SD2S_M02_z30a000")
#RPV
BRF2, R2, hdr2, bhr_uw2, bhr_dw2, model = produce_rami_results("HOM00_RPV_S00S_M03_z30a000")#("HOM00_RPV_EC6L_M03_z30a000")
# for RossLi
BRF3, R3, hdr3, bhr_uw3, bhr_dw3, model = produce_rami_results("HOM00_RLI_S00S_M03_z30a000")#("HOM00_RLI_EC6L_M03_z30a000") 

##

#for i = 1:10:2000
for i = 1000:1050
    println( i, " ",  all_scenarios[i]["observations"][1]["atmosphere"]["atmosphere_type"], " ", all_scenarios[i]["name"])
end

#=
for scenario in all_scenarios
    @info "Processing " * scenario["name"]
    try
        BRF, R,mods = produce_rami_results(scenario["name"])
    catch e
        @info e
    end
end
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