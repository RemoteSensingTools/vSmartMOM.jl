# Generate the FixedPressureGrid golden file used by test_vertical_grid_psurf.jl.
# Run once with the UNMODIFIED code (commit 7acab570, before vertical-grid
# conventions existed), from a copy of its tree whose test/ holds this file's
# fixture:  julia --project=. reference/generate_vertical_grid_golden.jl OUT
# The default convention must reproduce these arrays bit for bit.
isdefined(@__MODULE__, :_vg_parameters) ||
    include(joinpath(@__DIR__, "..", "vertical_grid_psurf_fixture.jl"))

function golden_arrays()
    arrays = Pair{String,Array{Float64}}[]
    for kind in (:line, :absco)
        params = _vg_parameters(Float64, kind; water=true, aerosols=true)
        # No vertical_grid keyword: this is the pre-convention API.
        model, lin = model_from_parameters(LinMode(), params;
            compute_aerosol_microphysics_jacobians=false)
        push!(arrays, "$(kind)_tau_abs_psurf" => lin.τ̇_abs_psurf[1],
                      "$(kind)_tau_rayl_psurf" => lin.τ̇_rayl_psurf[1],
                      "$(kind)_tau_aer_psurf" => lin.τ̇_aer_psurf[1])
        result = rt_run(model, lin, 2, size(lin.τ̇_abs[1], 1), 1)
        push!(arrays, "$(kind)_toa" => Array(result.toa),
                      "$(kind)_toa_jacobian" => Array(result.toa_jacobian))
    end
    return arrays
end

function write_golden(path, arrays)
    open(path, "w") do io
        println(io, "# FixedPressureGrid golden values from vSmartMOM 7acab570 (Julia $(VERSION), ",
                Sys.ARCH, "); name dims : hex(Float64 bits)")
        for (name, a) in arrays
            println(io, name, " ", join(size(a), "x"), " : ",
                    join(string.(reinterpret(UInt64, vec(a)); base=16), " "))
        end
    end
end

function read_golden(path)
    result = Dict{String,Array{Float64}}()
    for line in eachline(path)
        startswith(line, "#") && continue
        head, data = split(line, " : ")
        name, dims = split(head)
        shape = Tuple(parse.(Int, split(dims, "x")))
        bits = parse.(UInt64, split(data); base=16)
        result[name] = reshape(reinterpret(Float64, bits), shape)
    end
    return result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    write_golden(ARGS[1], golden_arrays())
end
