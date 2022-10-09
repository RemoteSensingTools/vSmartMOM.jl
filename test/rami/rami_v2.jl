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



# Get all scenarios
rami_json = "test/rami/RAMI4ATM_experiments_v1.0.json";
all_scenarios = JSON.parsefile(rami_json);

# R, models = produce_rami_results("HOM00_WHI_SD2S_M02_z30a000")
#R, models = produce_rami_resplults("HOM00_LAM_ED2D_M12_z30a000")
#R, models = produce_rami_results("HOM00_WHI_A00S_M02_z30a000")
#BRF, R, model = produce_rami_results("HOM00_WHI_A00S_M12_z30a000")
#white lambertian
BRF1, R1, hdrf1, bhr1, model = produce_rami_results("HOM00_LAM_S00S_M03_z30a000")#("HOM00_WHI_EC6L_M03_z30a000")
BRF25, R11, hdrf11, bhr11, model1 = produce_rami_results("HOM25_LAM_S00S_M03_z30a000")
BRF35, R11, hdrf11, bhr11, model1 = produce_rami_results("HOM35_LAM_S00S_M03_z30a000")
BRF45, R11, hdrf11, bhr11, model1 = produce_rami_results("HOM45_LAM_S00S_M03_z30a000")
#BRF1, R, model = produce_rami_results("HOM00_WHI_SD2S_M02_z30a000")
#RPV
BRF2, R2, hdrf2, bhr2, model = produce_rami_results("HOM00_RPV_S00S_M03_z30a000")#("HOM00_RPV_EC6L_M03_z30a000")
# for RossLi
BRF3, R3, hdrf3, bhr3, model = produce_rami_results("HOM00_RLI_S00S_M03_z30a000")#("HOM00_RLI_EC6L_M03_z30a000") 

##

#for i = 1:10:3000
for i = 1000:10:2050
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