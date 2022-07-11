using JSON
using vSmartMOM
using DelimitedFiles
using Distributions

# TODO: Add proper band limits and resolutions for gases:
const sentinel_band_to_wn = Dict([("2" , 1e7 ./ reverse(collect(455.0:0.1:534.0))),    # Blue
                                  ("3" , 1e7 ./ [559.8 559.9]),    # Green
                                  ("4" , 1e7 ./ [664.6 664.7]),    # Red
                                  ("8A", 1e7 ./ [864.7 864.8]),    # Narrow NIR 
                                  ("11", 1e7 ./ [1613.7 1613.8]),  # SWIR
                                  ("12", 1e7 ./ [2202.4 2202.5])]) # SWIR 


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
                              rami_json::String = "test/rami/RAMI4ATM_experiments_v1.0.json",
                              default_params::String = "test/rami/RamiGas.yaml")

    # Get rami scenario from experiment name
    all_scenarios = JSON.parsefile(rami_json);
    filter!(x->x["name"] == experiment_name, all_scenarios);
    @assert length(all_scenarios) ≤ 1 "Multiple matching experiment names!"
    @assert length(all_scenarios) ≥ 1 "Experiment name not found in JSON input"
    curr_scenario = all_scenarios[1]["observations"][1];

    # Get a default set of Rayleigh params
    params = vSmartMOM.parameters_from_yaml(default_params)

    # ########################################
    # Change params fields to match rami input
    # ########################################

    # There are 7 rami categories to be mindful of: 
    # atmosphere, name, canopy, illumination, measures, time, surface

    rami_atmosphere = curr_scenario["atmosphere"]
    rami_name = curr_scenario["name"]
    rami_canopy = curr_scenario["canopy"]
    rami_illumination = curr_scenario["illumination"]
    rami_measures = curr_scenario["measures"]
    rami_time = curr_scenario["time"]
    rami_surface = curr_scenario["surface"]

    # ########################################
    # 1. Atmosphere
    # ########################################

    # aerosols 
    if length(rami_atmosphere["aerosols"]) > 0
        # TODO We still need to check the two fractions of coarse and fine! Right now, only one is
        # 
        if startswith(rami_atmosphere["aerosols"][1]["name"], "D")

            μ_fine   = 0.0478666
            σ_fine   = 1.87411
            μ_coarse = 0.604127
            σ_coarse = 1.75172
            n_coarse = 0.00332189
            # These refractive indices depend on the band. 
            # TODO: We need to pull those from the table actually, as they are band dependent 
            # (in vSmartMOM, we so far only have one ni/nr per type).
            # See https://rami-benchmark.jrc.ec.europa.eu/_www/RAMI4ATM/down/RAMI4ATM_aerosols_v1.0/refractive_index/continental.txt
            # Need to generalize
            nᵣ = 1.4434925925925925
            nᵢ = 0.0015797

        elseif startswith(rami_atmosphere["aerosols"][1]["name"], "C")

            #TODO: Change, not just one refractice index
            μ_fine   = 0.0807989
            σ_fine   = 1.50180
            μ_coarse = 0.682651
            σ_coarse = 2.10400
            n_coarse   = 0.00046373
            nᵣ = 1.477538814814815
            nᵢ = 0.004342592592592592

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
        RT_aerosol = vSmartMOM.CoreRT.RT_Aerosol(aero, 0.550, 89880.0, 5000.0)

        # Assemble scattering parameters
        scattering_params = vSmartMOM.CoreRT.ScatteringParameters([RT_aerosol], 30.0, 2500, 0.550, vSmartMOM.Scattering.NAI2())

        params.scattering_params = scattering_params

    end

    # p/T profiles 
    params.p = reverse(US_std_atmos[:,2])
    params.T = reverse(US_std_atmos[:,3])
    params.T = (params.T[1:end-1] + params.T[2:end])/ 2
    params.q = zero(params.T)

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

    # Viewing zenith 
    vza_start = rami_measures[1]["vza_start"]["value"]
    vza_end   = rami_measures[1]["vza_end"]["value"]
    vza_step  = rami_measures[1]["vza_step"]["value"]
    params.vza  = collect(vza_start:vza_step:vza_end)

    # Viewing azimuth (TODO: Check definitions!)
    vaa = rami_measures[1]["delta_vaa"]["value"]
    params.vaz  = repeat([vaa], length(params.vza))

    # Spectral bands
    # Can expand to run multiple bands, just one for testing
    spec_bands        = sentinel_band_to_wn[rami_measures[1]["bands"][1]]
    params.spec_bands = [spec_bands]
    #@show params.spec_bands
    # ########################################
    # 5. Time
    # ########################################

    # ########################################
    # 6. Surface
    # ########################################
    @assert rami_surface["name"] in ["WHI", "BLA", "LAM"] && startswith(curr_scenario["name"], "HOM00") "Currently only supporting Lambertian (HOM00) surfaces"
    params.brdf = [vSmartMOM.LambertianSurfaceScalar(rami_surface["surface_parameters"]["reflectance"])]
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