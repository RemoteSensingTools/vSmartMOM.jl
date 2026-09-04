#!/usr/bin/env julia

"""Generate and validate the static 80-state bottom-layer CO2 campaign tables."""

include(joinpath(@__DIR__, "bottom_layer_truth_common.jl"))
using .BottomLayerTruthCommon

function main()
    paths = write_campaign_tables!()
    states = read_bottom_states(paths.state_file)
    @info("bottom-layer truth tables ready",
          states=length(states),
          state_table=paths.state_file,
          sha256=state_table_sha256(paths.state_file))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
