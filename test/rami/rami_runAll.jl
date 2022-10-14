using JSON
using vSmartMOM
using DelimitedFiles
using Distributions

include("rami_tools.jl")

# Get all scenarios
rami_json = "test/rami/RAMI4ATM_experiments_v1.0.json";
all_scenarios = JSON.parsefile(rami_json);

for scenario in all_scenarios
    @info "Processing " * scenario["name"]
    try
        BRF, R, hdrf, bhr, model = produce_rami_results(scenario["name"])
    catch e
        @info e
    end
end
