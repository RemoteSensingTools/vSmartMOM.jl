using JSON
using vSmartMOM
using DelimitedFiles
using Distributions

include("rami_tools.jl")

# Get all scenarios
rami_json = "test/rami/RAMI4ATM_experiments_v1.0.json";
all_scenarios = JSON.parsefile(rami_json);
folder = "/home/cfranken/rami_IQUV/"

for scenario in all_scenarios
    @info "Processing " * scenario["name"]
    try

        file = folder * scenario["name"] * "-brfpp_" * "vSmartMOM-JPL.mes"
        if !isfile(file)
            #@show file
            if startswith(scenario["name"],"HOM00")
                BRF, R, hdrf, bhr, model = produce_rami_results(scenario["name"], polarized=true, folder=folder)
            end
        end

    catch e
        @info e
    end
end
