using OCORamanWorkflow

function parse_float_list(value::AbstractString)
    isempty(strip(value)) && return Float64[]
    return parse.(Float64, split(value, ","))
end

function env_string_pairs(pairs::Pair{String,String}...)
    env = copy(ENV)
    for (k, v) in pairs
        env[k] = v
    end
    return env
end

sza_list = parse_float_list(get(ENV, "OCORAMAN_SZA_LIST", "0,30,50,70"))
albedo_list = parse_float_list(get(ENV, "OCORAMAN_ALBEDO_LIST", "0.0,0.1,0.3,0.5"))
yaml = get(ENV, "GRID_YAML",
    joinpath(repo_root(), "test", "test_parameters", "O2_parameters2_SIF_grid.yaml"))
batch_id = get(ENV, "OCORAMAN_RUN_ID", run_id("o2a_grid"))
batch_dir = get(ENV, "OCORAMAN_BATCH_DIR", joinpath(output_root(), batch_id))
dry_run = get(ENV, "DRY_RUN", "1") != "0"
smoke_script = joinpath(@__DIR__, "02_o2a_smoke_output.jl")

write_manifest(batch_dir;
    script=@__FILE__,
    yaml=yaml,
    extra=Dict(
        "SZA" => sza_list,
        "ALBEDO" => albedo_list,
        "DRY_RUN" => dry_run,
        "child_script" => smoke_script,
    ))

println("OCORaman O2A grid runner")
println("  batch_dir: ", batch_dir)
println("  yaml:      ", yaml)
println("  dry_run:   ", dry_run)
println("  cases:     ", length(sza_list) * length(albedo_list))

for sza in sza_list, albedo in albedo_list
    case_tag = "sza$(replace(string(sza), "." => "p"))_alb$(replace(string(albedo), "." => "p"))"
    case_dir = joinpath(batch_dir, case_tag)
    child_env = env_string_pairs(
        "OCORAMAN_RUN_ID" => "$(batch_id)_$(case_tag)",
        "GRID_YAML" => yaml,
        "GRID_OUTDIR" => case_dir,
        "SZA" => string(sza),
        "ALBEDO" => string(albedo),
        "SIF_ON" => get(ENV, "SIF_ON", "1"),
    )
    cmd = `$(Base.julia_cmd()) --project=$(workflow_root()) $smoke_script`
    println("  ", dry_run ? "[dry-run] " : "[run] ", case_tag)
    println("    GRID_OUTDIR=", case_dir)
    if !dry_run
        run(setenv(cmd, child_env))
    end
end
