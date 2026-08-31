using OCORamanWorkflow

run = get(ENV, "OCORAMAN_RUN_ID", run_id("o2a_smoke"))
outdir = get(ENV, "GRID_OUTDIR", joinpath(output_root(), run, "grid_O2A"))
yaml = get(ENV, "GRID_YAML",
    joinpath(repo_root(), "test", "test_parameters", "O2_parameters2_SIF_grid.yaml"))

ENV["GRID_YAML"] = yaml
ENV["GRID_OUTDIR"] = outdir
ENV["SZA"] = get(ENV, "SZA", "19")
ENV["ALBEDO"] = get(ENV, "ALBEDO", "0.0")
ENV["SIF_ON"] = get(ENV, "SIF_ON", "1")

write_manifest(outdir;
    script=@__FILE__,
    yaml=yaml,
    extra=Dict(
        "SZA" => ENV["SZA"],
        "ALBEDO" => ENV["ALBEDO"],
        "SIF_ON" => ENV["SIF_ON"],
        "driver" => "test/benchmarks/creategrid_O2Aband_RamanSIF.jl",
    ))

println("OCORaman O2A smoke output")
println("  YAML:   ", ENV["GRID_YAML"])
println("  outdir: ", ENV["GRID_OUTDIR"])
println("  SZA:    ", ENV["SZA"])
println("  albedo: ", ENV["ALBEDO"])
println("  SIF_ON: ", ENV["SIF_ON"])

include(joinpath(repo_root(), "test", "benchmarks", "creategrid_O2Aband_RamanSIF.jl"))

try
    print_grid_summary(latest_jld2(outdir))
catch err
    @warn "Could not summarize smoke output" exception=(err, catch_backtrace())
end
