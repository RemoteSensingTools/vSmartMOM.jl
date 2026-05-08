# =============================================================================
# benchmark_natraj_cds.jl -- full Natraj et al. (2009) CDS Rayleigh benchmark.
#
# Reference:
#   Natraj, V., Li, K.-F., & Yung, Y. L. 2009, ApJ, 691, 1909-1920.
#   DOI: 10.1088/0004-637X/691/2/1909
#
# Run from the test environment:
#   cd test
#   julia --project=. benchmarks/benchmark_natraj_cds.jl --quick
#   julia --project=. benchmarks/benchmark_natraj_cds.jl
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using Printf
using Statistics

const BENCH_DIR = @__DIR__
const YAML_PATH = joinpath(BENCH_DIR, "natraj.yaml")
const TABLE_DIR = joinpath(BENCH_DIR, "natraj_cds", "CDS")

const CDS_TAUS = [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.0]
const CDS_TAU_NAMES = Dict(0.02 => "0.02", 0.05 => "0.05", 0.1 => "0.1",
                           0.15 => "0.15", 0.25 => "0.25", 0.5 => "0.5",
                           1.0 => "1")
const CDS_ALBEDOS = [0.0, 0.25, 0.8]
const CDS_MU0S = [0.1, 0.2, 0.4, 0.6, 0.8, 0.92, 1.0]
const CDS_MUS = [0.02, 0.06, 0.1, 0.16, 0.2, 0.28, 0.32, 0.4,
                 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.0]
const CDS_PHIS = collect(0.0:30.0:180.0)
const OBS_VZA = vcat((acosd.(CDS_MUS) for _ in CDS_PHIS)...)
const OBS_VAZ = vcat((fill(phi, length(CDS_MUS)) for phi in CDS_PHIS)...)

const STOKES_INDEX = Dict(:I => 1, :Q => 2, :U => 3)
const DIR_NAME = Dict(:up => "UP", :dn => "DN")
const LAMBERTIAN = vSmartMOM.CoreRT.LambertianSurfaceScalar

mutable struct Options
    taus::Vector{Float64}
    albedos::Vector{Float64}
    mu0s::Vector{Float64}
    stokes::Vector{Symbol}
    directions::Vector{Symbol}
    rel_floor::Float64
    quick::Bool
    dry_run::Bool
    verbose::Bool
    csv_path::Union{Nothing,String}
end

Options() = Options(copy(CDS_TAUS), copy(CDS_ALBEDOS), copy(CDS_MU0S),
                    [:I, :Q, :U], [:up, :dn], 0.01, false, false, false,
                    nothing)

qkey(x::Real) = round(Int, Float64(x) * 1_000_000)
table_key(stokes, direction, tau, albedo, mu0, mu) =
    (stokes, direction, qkey(tau), qkey(albedo), qkey(mu0), qkey(mu))
obs_index(iphi, imu) = (iphi - 1) * length(CDS_MUS) + imu

function usage()
    println("""
    Usage:
      julia --project=. benchmarks/benchmark_natraj_cds.jl [options]

    Options:
      --quick                 Run tau=0.5, albedo=0, mu0=0.2 only.
      --dry-run               Parse tables and print selected work only.
      --verbose               Do not suppress model/rt_run output.
      --tau=0.5,1             Select optical thickness values.
      --albedo=0,0.25         Select Lambertian surface albedos.
      --mu0=0.2,0.8           Select incident cosines.
      --stokes=I,Q,U          Select Stokes components.
      --directions=up,dn      Select upwelling TOA and/or downwelling BOA.
      --rel-floor=0.01        Relative stats ignore |truth| below this.
      --csv=path              Write scalar comparison rows to a CSV file.
      -h, --help              Show this message.
    """)
end

function split_value_arg(args, i)
    arg = args[i]
    if occursin("=", arg)
        lhs, rhs = split(arg, "="; limit = 2)
        return lhs, rhs, i
    end
    i == length(args) && error("missing value for $(arg)")
    return arg, args[i + 1], i + 1
end

parse_float_list(s) =
    [parse(Float64, strip(x)) for x in split(s, ",") if !isempty(strip(x))]

function parse_stokes_list(s)
    out = Symbol[]
    for token in split(uppercase(replace(s, "," => " ")))
        if token == "IQU"
            append!(out, [:I, :Q, :U])
        elseif token in ("I", "Q", "U")
            push!(out, Symbol(token))
        else
            error("unsupported Stokes component: $token")
        end
    end
    return unique(out)
end

function parse_direction_list(s)
    out = Symbol[]
    for token in split(lowercase(replace(s, "," => " ")))
        if token in ("up", "upwelling")
            push!(out, :up)
        elseif token in ("dn", "down", "downwelling")
            push!(out, :dn)
        else
            error("unsupported direction: $token")
        end
    end
    return unique(out)
end

function validate_subset(values, allowed, label)
    for value in values
        any(isapprox(value, x; atol = 1e-10, rtol = 0) for x in allowed) ||
            error("$label=$value is not in the CDS grid")
    end
    return values
end

function parse_options(args = ARGS)
    opts = Options()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            usage()
            exit(0)
        elseif arg == "--quick"
            opts.quick = true
        elseif arg == "--dry-run"
            opts.dry_run = true
        elseif arg == "--verbose"
            opts.verbose = true
        elseif startswith(arg, "--tau")
            _, value, i = split_value_arg(args, i)
            opts.taus = validate_subset(parse_float_list(value), CDS_TAUS, "tau")
        elseif startswith(arg, "--albedo")
            _, value, i = split_value_arg(args, i)
            opts.albedos = validate_subset(parse_float_list(value), CDS_ALBEDOS, "albedo")
        elseif startswith(arg, "--mu0")
            _, value, i = split_value_arg(args, i)
            opts.mu0s = validate_subset(parse_float_list(value), CDS_MU0S, "mu0")
        elseif startswith(arg, "--stokes")
            _, value, i = split_value_arg(args, i)
            opts.stokes = parse_stokes_list(value)
        elseif startswith(arg, "--direction")
            _, value, i = split_value_arg(args, i)
            opts.directions = parse_direction_list(value)
        elseif startswith(arg, "--rel-floor")
            _, value, i = split_value_arg(args, i)
            opts.rel_floor = parse(Float64, value)
        elseif startswith(arg, "--csv")
            _, value, i = split_value_arg(args, i)
            opts.csv_path = value
        else
            error("unknown argument: $arg")
        end
        i += 1
    end

    if opts.quick
        opts.taus = [0.5]
        opts.albedos = [0.0]
        opts.mu0s = [0.2]
    end
    isempty(opts.stokes) && error("no Stokes components selected")
    isempty(opts.directions) && error("no directions selected")
    return opts
end

function tau_name(tau)
    for x in CDS_TAUS
        isapprox(tau, x; atol = 1e-10, rtol = 0) && return CDS_TAU_NAMES[x]
    end
    error("tau=$tau is not in the CDS grid")
end

function parse_cds_file!(tables, stokes, direction, tau)
    path = joinpath(TABLE_DIR, "$(stokes)_$(DIR_NAME[direction])_TAU_$(tau_name(tau))")
    isfile(path) || error("missing CDS table: $path")

    current_albedo = NaN
    open(path, "r") do io
        for line in eachline(io)
            text = strip(line)
            if startswith(lowercase(text), "albedo")
                current_albedo = parse(Float64, strip(split(text, "="; limit = 2)[2]))
                continue
            end

            parts = split(text)
            length(parts) == 9 || continue
            nums = Float64[]
            ok = true
            for part in parts
                value = tryparse(Float64, part)
                if value === nothing
                    ok = false
                    break
                end
                push!(nums, value)
            end
            ok || continue
            isnan(current_albedo) &&
                error("encountered data row before albedo in $path")

            mu0 = nums[1]
            mu = nums[2]
            tables[table_key(stokes, direction, tau, current_albedo, mu0, mu)] =
                nums[3:end]
        end
    end
    return tables
end

function load_reference_tables(opts)
    tables = Dict{Tuple{Symbol,Symbol,Int,Int,Int,Int},Vector{Float64}}()
    for stokes in opts.stokes, direction in opts.directions, tau in opts.taus
        parse_cds_file!(tables, stokes, direction, tau)
    end
    return tables
end

function maybe_quiet(f, verbose)
    verbose && return f()
    result = nothing
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            result = f()
        end
    end
    return result
end

function build_model(base_params, taus, albedo, mu0, verbose)
    maybe_quiet(verbose) do
        params = deepcopy(base_params)
        params.spec_bands = [fill(1e7 / 360.0, length(taus))]
        params.vza = OBS_VZA
        params.vaz = OBS_VAZ
        params.sza = acosd(mu0)
        params.brdf = [LAMBERTIAN(albedo)]

        model = model_from_parameters(params)
        model.τ_rayl[1] .= reshape(taus, :, 1)
        return model
    end
end

run_rt(model, verbose) = maybe_quiet(verbose) do
    CoreRT.rt_run(model, i_band = 1)
end

function direct_beam_correction(stokes, direction, tau, mu0, mu)
    # CDS downwelling I is diffuse-only. CoreRT.T includes the unscattered
    # direct beam at the solar ordinate; after multiplying by pi this is
    # 0.5 * exp(-tau / mu0). Q and U are unaffected.
    if direction == :dn && stokes == :I &&
       isapprox(mu, mu0; atol = 1e-12, rtol = 0)
        return 0.5 * exp(-tau / mu0)
    end
    return 0.0
end

function compare_case!(results, tables, opts, taus, albedo, mu0, up, dn)
    arrays = Dict(:up => π .* up, :dn => π .* dn)
    for (ispec, tau) in enumerate(taus)
        for direction in opts.directions
            radiance = arrays[direction]
            for stokes in opts.stokes
                istokes = STOKES_INDEX[stokes]
                for (imu, mu) in enumerate(CDS_MUS)
                    truth_values = tables[table_key(stokes, direction, tau, albedo, mu0, mu)]
                    for (iphi, phi) in enumerate(CDS_PHIS)
                        correction = direct_beam_correction(stokes, direction, tau, mu0, mu)
                        modeled = radiance[obs_index(iphi, imu), istokes, ispec] -
                                  correction
                        truth = truth_values[iphi]
                        abs_err = abs(modeled - truth)
                        rel_err = iszero(truth) ? NaN : abs_err / abs(truth)
                        push!(results, (;
                            stokes, direction, tau, albedo, mu0, mu, phi,
                            truth, modeled, abs_err, rel_err,
                            direct_correction = correction,
                        ))
                    end
                end
            end
        end
    end
    return results
end

finite_values(rows, field) =
    [Float64(getproperty(row, field)) for row in rows
     if isfinite(Float64(getproperty(row, field)))]

function stats_tuple(values)
    isempty(values) && return (; n = 0, mean = NaN, median = NaN, max = NaN)
    return (; n = length(values), mean = mean(values), median = median(values),
            max = maximum(values))
end

function summarize_group(label, rows, opts)
    abs_stats = stats_tuple(finite_values(rows, :abs_err))
    rel_rows = [row for row in rows
                if isfinite(row.rel_err) && abs(row.truth) >= opts.rel_floor]
    rel_stats = stats_tuple(finite_values(rel_rows, :rel_err))
    @printf("%-11s %7d %12.4e %12.4e %12.4e %7d %12.4e %12.4e %12.4e\n",
            label, abs_stats.n, abs_stats.mean, abs_stats.median,
            abs_stats.max, rel_stats.n, rel_stats.mean, rel_stats.median,
            rel_stats.max)
end

function summarize(results, opts, elapsed_s)
    println("\n=== Natraj et al. (2009) CDS Rayleigh benchmark ===")
    @printf("Compared %d scalar table values in %.1f s\n", length(results), elapsed_s)
    @printf("Relative stats ignore |truth| < %.3g; absolute stats use all rows.\n",
            opts.rel_floor)
    println("Downwelling I is compared after removing the direct-beam ordinate.")
    println()
    @printf("%-11s %7s %12s %12s %12s %7s %12s %12s %12s\n",
            "group", "n_abs", "mean_abs", "median_abs", "max_abs",
            "n_rel", "mean_rel", "median_rel", "max_rel")
    summarize_group("overall", results, opts)
    for stokes in opts.stokes, direction in opts.directions
        rows = [row for row in results
                if row.stokes == stokes && row.direction == direction]
        summarize_group("$(stokes)_$(uppercase(string(direction)))", rows, opts)
    end
    return nothing
end

function write_csv(path, results)
    dir = dirname(abspath(path))
    isdir(dir) || mkpath(dir)
    open(path, "w") do io
        println(io, "stokes,direction,tau,albedo,mu0,mu,phi,truth,modeled,abs_err,rel_err,direct_correction")
        for row in results
            @printf(io, "%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                    string(row.stokes), uppercase(string(row.direction)),
                    row.tau, row.albedo, row.mu0, row.mu, row.phi, row.truth,
                    row.modeled, row.abs_err, row.rel_err, row.direct_correction)
        end
    end
    return path
end

function run_benchmark(opts = parse_options())
    tables = load_reference_tables(opts)
    n_runs = length(opts.albedos) * length(opts.mu0s)
    n_scalars = length(opts.taus) * n_runs *
                length(opts.directions) * length(opts.stokes) *
                length(CDS_MUS) * length(CDS_PHIS)

    println("Selected $n_runs RT runs and $n_scalars scalar CDS comparisons.")
    println("Tau values are batched as spectral points: $(opts.taus)")
    println("Tables: $TABLE_DIR")
    println("YAML:   $YAML_PATH")
    if opts.dry_run
        println("Dry run only; no RT model was executed.")
        return NamedTuple[]
    end

    base_params = parameters_from_yaml(YAML_PATH)
    results = NamedTuple[]
    started = time()
    case_i = 0
    for albedo in opts.albedos, mu0 in opts.mu0s
        case_i += 1
        @printf("[%3d/%3d] albedo=%-4.2g mu0=%-4.2g n_tau=%d\n",
                case_i, n_runs, albedo, mu0, length(opts.taus))
        model = build_model(base_params, opts.taus, albedo, mu0, opts.verbose)
        up, dn = run_rt(model, opts.verbose)
        compare_case!(results, tables, opts, opts.taus, albedo, mu0, up, dn)
    end

    elapsed_s = time() - started
    summarize(results, opts, elapsed_s)
    if opts.csv_path !== nothing
        write_csv(opts.csv_path, results)
        println("\nWrote CSV: $(opts.csv_path)")
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
end
