#!/usr/bin/env julia

# Chunked O2 A-band Raman/Cabannes/Rayleigh LUT driver.
#
# Purpose:
#   Same science workflow as createRamanLUT_O2A.jl, but intended for smaller
#   failure/restart units. Each run writes only the selected SZA/albedo subset
#   to its own NetCDF file and records the selected full-grid indices.
#
# What this automates:
#   - psurf-dependent p/T/q profile truncation/interpolation
#   - SIF-off Raman LUT generation through SurfaceSIF = 0
#   - Float32 GPU execution
#   - two broad O2 A-band chunks with 250 cm^-1 Raman shoulders
#   - unshouldered NetCDF output only
#
# Default grid:
#   rho       = 0:0.05:1
#   SZA       = 14 values with equally spaced cos(SZA), max SZA = 70 deg
#   VZA/VAZ   = 0/0
#   psurf     = 1000, 750, 500 hPa
#   SIF       = off only
#
# NetCDF variables:
#   stokes_rayleigh(psurf, sif, albedo, sza, stokes, wn)
#   stokes_cabannes(psurf, sif, albedo, sza, stokes, wn)
#   stokes_rrs(psurf, sif, albedo, sza, stokes, wn)
#   F0(wn)
#   SIF0(sif, wn)
#   tau_rayl(psurf, wn)
#   tau_abs(psurf, wn)
#   tau_aer(psurf, wn) = 0 for this no-aerosol LUT
#
# Chunk controls:
#   SZA_IDXS     = one-based full-grid SZA indices, e.g. 1:4 or 1,7,14
#   ALBEDO_IDXS  = one-based full-grid albedo indices, e.g. all or 1,7,21
#   FULL_VZAS    = comma/range list, or "sza" to use the full SZA grid
#   VZA_IDXS     = one-based selected VZA indices from FULL_VZAS
#   RAZS         = relative azimuth grid, e.g. 0:15:345
#   CHUNK_TAG    = optional output-name tag; otherwise built from indices
#
# Examples:
#   # Inspect axes, spectral sizes, and timing estimate only:
#   CUDA_DEVICE=1 DRY_RUN=1 julia --project=. \
#       workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl
#
#   # Example SZA chunks for one pressure slice:
#   CUDA_DEVICE=1 PSURFS=1000 SZA_IDXS=1:4 ALBEDO_IDXS=all \
#       julia --project=. workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl
#   CUDA_DEVICE=1 PSURFS=1000 SZA_IDXS=5:9 ALBEDO_IDXS=all \
#       julia --project=. workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl
#   CUDA_DEVICE=1 PSURFS=1000 SZA_IDXS=10:14 ALBEDO_IDXS=all \
#       julia --project=. workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl

using Dates
using Interpolations
using NCDatasets
using Printf
using Statistics
using YAML

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel

const FT = Float32

const REPO_ROOT = pkgdir(vSmartMOM)
const SCRIPT_PATH = abspath(get(ENV, "SCRIPT_PATH_OVERRIDE", @__FILE__))
const RUN_STARTED = now()
const CUDA_IMPORTED = Ref(false)

parse_float_list(text::AbstractString) =
    [parse(Float64, strip(x)) for x in split(text, ","; keepempty = false)]

function parse_string_list(text::AbstractString)
    token = strip(text)
    lowercase(token) in ("", "none", "[]") && return String[]
    vals = String[]
    for raw in split(token, ","; keepempty = false)
        part = strip(raw)
        isempty(part) || push!(vals, part)
    end
    return vals
end

function parse_path_list(text::AbstractString)
    vals = parse_string_list(text)
    return vals
end

function parse_float_grid(text::AbstractString)
    vals = Float64[]
    for raw in split(strip(text), ","; keepempty = false)
        part = strip(raw)
        if occursin(":", part)
            pieces = split(part, ":")
            length(pieces) in (2, 3) ||
                error("Float range '$part' must be start:stop or start:step:stop")
            nums = parse.(Float64, pieces)
            r = length(nums) == 2 ? (nums[1]:nums[2]) : (nums[1]:nums[2]:nums[3])
            append!(vals, collect(r))
        else
            push!(vals, parse(Float64, part))
        end
    end
    isempty(vals) && error("Float grid '$text' selected no values")
    return vals
end

function parse_bool_list(text::AbstractString)
    vals = Bool[]
    for item in split(text, ","; keepempty = false)
        token = lowercase(strip(item))
        if token in ("1", "true", "t", "yes", "y", "on", "sif", "sif_on")
            push!(vals, true)
        elseif token in ("0", "false", "f", "no", "n", "off", "nosif", "sif_off")
            push!(vals, false)
        else
            error("Cannot parse boolean token '$item' in SIF_STATES")
        end
    end
    return vals
end

function parse_index_selector(text::AbstractString, n::Integer, name::AbstractString)
    token = lowercase(strip(text))
    token in ("", "all", "*") && return collect(1:n)

    idxs = Int[]
    for raw in split(text, ","; keepempty = false)
        part = strip(raw)
        if occursin(":", part)
            pieces = split(part, ":")
            length(pieces) in (2, 3) ||
                error("$name range '$part' must be start:stop or start:step:stop")
            vals = parse.(Int, pieces)
            r = length(vals) == 2 ? (vals[1]:vals[2]) : (vals[1]:vals[2]:vals[3])
            append!(idxs, collect(r))
        else
            push!(idxs, parse(Int, part))
        end
    end
    isempty(idxs) && error("$name selected no indices")
    any(i -> i < 1 || i > n, idxs) &&
        error("$name contains index outside 1:$n: $(join(idxs, ","))")
    unique_idxs = unique(idxs)
    length(unique_idxs) == length(idxs) ||
        error("$name contains duplicate indices: $(join(idxs, ","))")
    return unique_idxs
end

function index_tag(idxs::AbstractVector{<:Integer}, nfull::Integer)
    length(idxs) == nfull && all(idxs .== collect(1:nfull)) && return "all"
    if !isempty(idxs) && all(diff(idxs) .== 1)
        return @sprintf("%03d-%03d", first(idxs), last(idxs))
    end
    return join((@sprintf("%03d", i) for i in idxs), "p")
end

function sza_axis_equal_mu(n::Integer, sza_max::Real)
    n > 1 || return [0.0]
    sza_max <= 70 || error("SZA_MAX_DEG must be <= 70; got $sza_max")
    mu = collect(range(cosd(Float64(sza_max)), stop = 1.0, length = n))
    return reverse(acosd.(mu))
end

function vza_axis_from_env(text::AbstractString, full_sza_axis::AbstractVector{<:Real})
    token = lowercase(strip(text))
    token in ("sza", "sza_grid", "same_as_sza") && return Float64.(full_sza_axis)
    return parse_float_grid(text)
end

const DRY_RUN = get(ENV, "DRY_RUN", "0") == "1"
const FORCE = get(ENV, "FORCE", "0") == "1"
const RESUME = get(ENV, "RESUME", "0") == "1"
const ARCH = uppercase(get(ENV, "ARCH", "GPU"))
const CUDA_DEVICE = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
FORCE && RESUME && error("Use either FORCE=1 or RESUME=1, not both.")

const BAND_TAG = get(ENV, "BAND_TAG", "O2A")
const BAND_TITLE = get(ENV, "BAND_TITLE", "O2 A-band")
const BAND_OUTPUT_PREFIX = get(ENV, "BAND_OUTPUT_PREFIX",
    "RamanLUT_$(BAND_TAG)_float32_gpu")

const PARAM_YAML = get(ENV, "PARAM_YAML",
    joinpath(REPO_ROOT, "workflows", "OCORaman", "config", "paramsRamanLUT_O2A.yaml"))
const PROFILE_SOURCE_NOTE =
    "atmospheric_profile.p/T/q copied from ParamsEMIT_MODTRANcomp_newLUT.yaml; " *
    "AFGL 1986 Midlatitude Winter, model/table 1c, on the EMIT/MODTRAN " *
    "20 hPa grid; H2O VMR is derived internally from specific humidity q"

const HITRAN_LUT_DIR = get(ENV, "VSMARTMOM_HITRAN_LUT_DIR",
    "/home/sanghavi/data/HITRAN_LUTs")
const USE_O2_ABSCO = get(ENV, "USE_O2_ABSCO", "1") == "1"
const FIXED_MOLECULES = parse_string_list(get(ENV, "FIXED_MOLECULES", "O2"))
const VARIABLE_MOLECULES = parse_string_list(get(ENV, "VARIABLE_MOLECULES", ""))

const O2_ABSCO = get(ENV, "O2_ABSCO",
    "/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2")
const DEFAULT_H2O_LUT = begin
    joinpath(HITRAN_LUT_DIR, "H2O.jld2")
end
const H2O_LUT = get(ENV, "H2O_LUT", isfile(DEFAULT_H2O_LUT) ? DEFAULT_H2O_LUT : "")
const CIA_FILES = parse_path_list(get(ENV, "CIA_FILES", ""))
const MTCKD_FILE = get(ENV, "MTCKD_FILE", "")

const SOLAR_OUT = get(ENV, "SOLAR_OUT",
    joinpath(homedir(), "Raman_misc", "workflows", "worktrees",
             "uni_vSmartMOM_sanghavi_2025-03-18", "src", "SolarModel", "solar.out"))

const WN_STEP = parse(Float64, get(ENV, "WN_STEP", "0.1"))
const SHOULDER_CM = parse(Float64, get(ENV, "SHOULDER_CM", "250.0"))
const WL_EDGES_NM = parse_float_list(get(ENV, "WL_EDGES_NM", "785,762.5,740"))
const NBANDS = length(WL_EDGES_NM) - 1

const PSURF_AXIS = parse_float_list(get(ENV, "PSURFS", "1000,750,500"))
const SIF_AXIS = parse_bool_list(get(ENV, "SIF_STATES", "0"))
const FULL_ALBEDO_AXIS = parse_float_list(get(ENV, "FULL_ALBEDOS",
    get(ENV, "ALBEDOS", join(string.(0.0:0.05:1.0), ","))))
const NSZA = parse(Int, get(ENV, "NSZA", "14"))
const SZA_MAX_DEG = parse(Float64, get(ENV, "SZA_MAX_DEG", "70.0"))
const FULL_SZA_AXIS = sza_axis_equal_mu(NSZA, SZA_MAX_DEG)
const FULL_VZA_AXIS = vza_axis_from_env(get(ENV, "FULL_VZAS", "0.0"), FULL_SZA_AXIS)
const ALBEDO_IDXS = parse_index_selector(get(ENV, "ALBEDO_IDXS", "all"),
    length(FULL_ALBEDO_AXIS), "ALBEDO_IDXS")
const SZA_IDXS = parse_index_selector(get(ENV, "SZA_IDXS", "all"),
    length(FULL_SZA_AXIS), "SZA_IDXS")
const VZA_IDXS = parse_index_selector(get(ENV, "VZA_IDXS", "all"),
    length(FULL_VZA_AXIS), "VZA_IDXS")
const ALBEDO_AXIS = FULL_ALBEDO_AXIS[ALBEDO_IDXS]
const SZA_AXIS = FULL_SZA_AXIS[SZA_IDXS]
const MU0_AXIS = cosd.(SZA_AXIS)
const VZA_AXIS = FULL_VZA_AXIS[VZA_IDXS]
const RAZ_AXIS = parse_float_grid(get(ENV, "RAZS", "0.0"))
const OBS_VZA_AXIS = repeat(VZA_AXIS, inner = length(RAZ_AXIS))
const OBS_VAZ_AXIS = repeat(RAZ_AXIS, outer = length(VZA_AXIS))

const PROFILE_REDUCTION = parse(Int, get(ENV, "PROFILE_REDUCTION", "12"))
const NSTREAMS = parse(Int, get(ENV, "NSTREAMS", "3"))

const REF_NSPEC = parse(Float64, get(ENV, "REFERENCE_NSPEC", "12747"))
const REF_SECONDS_PER_SCENE = parse(Float64, get(ENV, "REFERENCE_SECONDS_PER_SCENE", "3650"))

const OUTDIR = get(ENV, "OUTDIR", joinpath(homedir(), "data", "RamanSIFgrid"))
const PSURF_TAG = length(PSURF_AXIS) == 1 ? @sprintf("psurf%04.0fhpa", PSURF_AXIS[1]) : "psurf_grid"
const VIEW_TAG = (length(FULL_VZA_AXIS) == 1 && only(FULL_VZA_AXIS) == 0.0 &&
                  length(VZA_IDXS) == 1 && only(VZA_IDXS) == 1 &&
                  length(RAZ_AXIS) == 1 && only(RAZ_AXIS) == 0.0) ? "" :
    "_vza$(index_tag(VZA_IDXS, length(FULL_VZA_AXIS)))_razn$(lpad(length(RAZ_AXIS), 3, '0'))"
const DEFAULT_CHUNK_TAG =
    "sza$(index_tag(SZA_IDXS, length(FULL_SZA_AXIS)))_alb$(index_tag(ALBEDO_IDXS, length(FULL_ALBEDO_AXIS)))$(VIEW_TAG)"
const CHUNK_TAG = get(ENV, "CHUNK_TAG", DEFAULT_CHUNK_TAG)
const OUT_NC = get(ENV, "OUT_NC",
    joinpath(OUTDIR, "$(BAND_OUTPUT_PREFIX)_$(PSURF_TAG)_$(CHUNK_TAG).nc"))

label_float(x::Real; digits::Int = 2) =
    replace(@sprintf("%.*f", digits, Float64(x)), "." => "p")

yaml_array(xs) = "[" * join(string.(xs), ", ") * "]"

yaml_string_array(xs::AbstractVector{<:AbstractString}) =
    isempty(xs) ? "[]" : "[" * join(xs, ", ") * "]"

function yaml_value(x)
    if x isa AbstractVector
        return "[" * join(string.(x), ", ") * "]"
    end
    return string(x)
end

function lut_path_for_molecule(mol::AbstractString)
    env_key = uppercase(mol) * "_LUT"
    haskey(ENV, env_key) && return ENV[env_key]
    uppercase(mol) == "O2" && USE_O2_ABSCO && return O2_ABSCO
    return joinpath(HITRAN_LUT_DIR, "$(mol).jld2")
end

function absorption_lut_paths()
    paths = String[]
    seen = Set{String}()
    for mol in vcat(FIXED_MOLECULES, VARIABLE_MOLECULES)
        uppercase(mol) == "H2O" &&
            error("Do not list H2O in FIXED_MOLECULES/VARIABLE_MOLECULES; q plus H2O_LUT activates water vapor.")
        path = lut_path_for_molecule(mol)
        if !(path in seen)
            push!(paths, path)
            push!(seen, path)
        end
    end
    if !isempty(H2O_LUT) && !(H2O_LUT in seen)
        push!(paths, H2O_LUT)
    end
    return paths
end

function validate_absorption_files!()
    for path in absorption_lut_paths()
        isfile(path) || error("Absorption LUT not found: $path")
    end
    for path in CIA_FILES
        isfile(path) || error("CIA file not found: $path")
    end
    !isempty(MTCKD_FILE) && !isfile(MTCKD_FILE) &&
        error("MTCKD_FILE not found: $MTCKD_FILE")
    return nothing
end

function absorption_note()
    gas = isempty(vcat(FIXED_MOLECULES, VARIABLE_MOLECULES)) ?
        "no listed line absorbers" :
        join(vcat(FIXED_MOLECULES, VARIABLE_MOLECULES), ", ")
    h2o = isempty(H2O_LUT) ? "" : "; H2O is q-driven"
    cia = isempty(CIA_FILES) ? "" : "; CIA enabled"
    mtckd = isempty(MTCKD_FILE) ? "" : "; MT_CKD H2O continuum enabled"
    return "tau_abs is gas absorption from $(gas)$(h2o)$(cia)$(mtckd); no aerosol implies tau_aer=0, not tau_abs=0."
end

function ensure_cuda()
    if !CUDA_IMPORTED[]
        @eval import CUDA
        CUDA_IMPORTED[] = true
    end
    return getfield(@__MODULE__, :CUDA)
end

function warmup_cuda_linalg!(CUDA)
    A = CUDA.CuArray(FT[1 0; 0 1])
    B = CUDA.CuArray(FT[1 2; 3 4])
    C = A * B
    CUDA.synchronize()
    @info "CUDA/cuBLAS warmup" check = Array(C)[1, 1]
    A = B = C = nothing
    GC.gc()
    CUDA.reclaim()
    return nothing
end

struct ProfileCase
    psurf::Float64
    p::Vector{Float64}
    T::Vector{Float64}
    q::Vector{Float64}
end

struct BandDef
    solve_start::Float64
    solve_stop::Float64
    output_start::Float64
    output_stop::Float64
    solve_wn::Vector{FT}
    keep::Vector{Int}
    out_range::UnitRange{Int}
end

function load_profile_from_yaml(path::AbstractString)
    isfile(path) || error("Profile YAML not found: $path")
    cfg = YAML.load_file(path)
    prof = cfg["atmospheric_profile"]
    p = Float64.(prof["p"])
    T = Float64.(prof["T"])
    q = Float64.(prof["q"])
    length(p) == length(T) + 1 ||
        error("Profile p length $(length(p)) must equal T length + 1 ($(length(T)+1))")
    length(T) == length(q) ||
        error("Profile T/q length mismatch: $(length(T)) vs $(length(q))")
    return (; p, T, q)
end

function load_vmr_from_yaml(path::AbstractString)
    cfg = YAML.load_file(path)
    abs_cfg = get(cfg, "absorption", Dict{String, Any}())
    raw_vmr = get(abs_cfg, "vmr", Dict{String, Any}())
    vmr = Dict{String, Any}()
    for (k, v) in raw_vmr
        vmr[string(k)] = v
    end
    get!(vmr, "O2", 0.21)
    return vmr
end

function interp_at_p(profile::AbstractVector, p_centers::AbstractVector, p_target::Real)
    k = searchsortedlast(p_centers, p_target)
    k < 1 && return profile[1]
    k >= length(p_centers) && return profile[end]
    lp1, lp2 = log(p_centers[k]), log(p_centers[k + 1])
    frac = (log(p_target) - lp1) / (lp2 - lp1)
    return profile[k] + frac * (profile[k + 1] - profile[k])
end

function truncate_profile(base, psurf::Real)
    p0, T0, q0 = base.p, base.T, base.q
    p_surf = Float64(psurf)
    p0[1] < p_surf <= p0[end] ||
        error("Requested psurf=$p_surf hPa outside profile range $(p0[1])--$(p0[end]) hPa")

    k = searchsortedlast(p0, p_surf - eps(p_surf))
    if k >= length(p0)
        p = copy(p0)
        p[end] = p_surf
        return ProfileCase(p_surf, p, copy(T0), copy(q0))
    end

    p = vcat(p0[1:k], p_surf)
    n_layers = k
    p_centers = [(p0[i] + p0[i + 1]) / 2 for i in 1:(length(p0) - 1)]

    function trunc_layer(prof)
        out = similar(prof, n_layers)
        n_layers > 1 && (out[1:(n_layers - 1)] .= prof[1:(n_layers - 1)])
        out[n_layers] = interp_at_p(prof, p_centers, p_surf)
        return out
    end

    return ProfileCase(p_surf, p, trunc_layer(T0), trunc_layer(q0))
end

function build_band_defs()
    NBANDS >= 1 || error("Need at least one band")
    defs = BandDef[]
    output_count = 0
    previous_output_stop = -Inf
    for ib in 1:NBANDS
        wl_hi = WL_EDGES_NM[ib]
        wl_lo = WL_EDGES_NM[ib + 1]
        wl_hi > wl_lo || error("WL_EDGES_NM must be descending wavelengths; got $wl_hi then $wl_lo")

        output_start = 1e7 / wl_hi
        output_stop = 1e7 / wl_lo
        solve_start = output_start - SHOULDER_CM
        solve_stop = output_stop + SHOULDER_CM
        solve_wn = collect(FT(solve_start):FT(WN_STEP):FT(solve_stop))

        keep = findall(ν -> (ν > previous_output_stop + 10eps(FT(output_stop))) &&
                            (ν >= FT(output_start) - 10eps(FT(output_start))) &&
                            (ν <= FT(output_stop) + 10eps(FT(output_stop))),
                       solve_wn)
        isempty(keep) && error("Band $ib has no unshouldered output points")

        out_range = (output_count + 1):(output_count + length(keep))
        output_count += length(keep)
        push!(defs, BandDef(solve_start, solve_stop, output_start, output_stop,
                            solve_wn, keep, out_range))
        previous_output_stop = output_stop
    end
    return defs
end

function write_case_yaml(path::AbstractString, profile::ProfileCase, sza::Real, band_defs, vmr)
    mkpath(dirname(path))
    lut_entries = ["\"$(path)\"" for path in absorption_lut_paths()]
    lut_entry = "[" * join(lut_entries, ", ") * "]"
    active_molecules = unique(vcat(FIXED_MOLECULES, VARIABLE_MOLECULES))

    open(path, "w") do io
        println(io, "# Auto-generated by $(basename(SCRIPT_PATH))")
        println(io, "# LUT tag: $BAND_TAG")
        println(io, "# LUT title: $BAND_TITLE")
        println(io, "# Source profile: $PROFILE_SOURCE_NOTE")
        println(io, "radiative_transfer:")
        println(io, "  spec_bands:")
        for b in band_defs
            println(io, "    - $(b.solve_start):$(WN_STEP):$(b.solve_stop)")
        end
        println(io, "  surface:")
        for _ in band_defs
            println(io, "    - LambertianSurfaceScalar{Float32}(0.0)")
        end
        println(io, "  polarization_type: Stokes_IQU()")
        println(io, "  nstreams: $(NSTREAMS)")
        println(io, "  truncation: auto")
        println(io, "  depol: -1")
        println(io, "  float_type: Float32")
        println(io, "  architecture: $(ARCH == "GPU" ? "Architectures.GPU()" : "Architectures.CPU()")")
        println(io)
        println(io, "geometry:")
        println(io, "  sza: $(Float64(sza))")
        println(io, "  vza: $(yaml_array(OBS_VZA_AXIS))")
        println(io, "  vaz: $(yaml_array(OBS_VAZ_AXIS))")
        println(io, "  obs_alt: [0]")
        println(io)
        println(io, "atmospheric_profile:")
        println(io, "  T: $(yaml_array(profile.T))")
        println(io, "  p: $(yaml_array(profile.p))")
        println(io, "  q: $(yaml_array(profile.q))")
        println(io, "  profile_reduction: $(PROFILE_REDUCTION)")
        println(io)
        println(io, "absorption:")
        println(io, "  fixed_molecules:")
        for _ in band_defs
            println(io, "    - $(yaml_string_array(FIXED_MOLECULES))")
        end
        println(io, "  variable_molecules:")
        for _ in band_defs
            println(io, "    - $(yaml_string_array(VARIABLE_MOLECULES))")
        end
        println(io, "  LUTfiles:")
        for _ in band_defs
            println(io, "    - $lut_entry")
        end
        println(io, "  vmr:")
        for mol in active_molecules
            haskey(vmr, mol) || error("No VMR entry for molecule $mol in PARAM_YAML=$PARAM_YAML")
            println(io, "    $mol: $(yaml_value(vmr[mol]))")
        end
        println(io, "  broadening: Voigt()")
        println(io, "  CEF: HumlicekWeidemann32SDErrorFunction()")
        println(io, "  wing_cutoff: 150")
        if !isempty(CIA_FILES)
            println(io, "  cia_files:")
            for path in CIA_FILES
                println(io, "    - \"$(path)\"")
            end
        end
        if !isempty(MTCKD_FILE)
            println(io, "  mtckd_file: \"$(MTCKD_FILE)\"")
        end
    end
    return path
end

function configure_architecture!()
    ARCH in ("GPU", "CPU") || error("ARCH must be GPU or CPU; got $ARCH")
    if ARCH == "GPU"
        CUDA = ensure_cuda()
        CUDA.functional() || error("CUDA is not functional, but ARCH=GPU")
        CUDA.device!(CUDA_DEVICE)
        CUDA.allowscalar(false)
        warmup_cuda_linalg!(CUDA)
        @info "CUDA device" name = CUDA.name(CUDA.device()) device = CUDA_DEVICE
    end
    return nothing
end

function with_lambertian_albedo(model::CoreRT.RTModel, albedo::Real)
    surfaces = [CoreRT.LambertianSurfaceScalar{FT}(FT(albedo))
                for _ in 1:length(CoreRT.get_surfaces(model))]
    return CoreRT.RTModel(
        model.architecture,
        model.solver,
        model.numerics,
        model.geometry,
        model.quad_points,
        model.atmosphere,
        model.optics,
        surfaces,
        model.sources,
    )
end

function solar_interpolator()
    if isfile(SOLAR_OUT)
        solar = solar_transmission_from_file(SOLAR_OUT)
        @info "Using solar transmission" SOLAR_OUT size = size(solar)
        return LinearInterpolation(FT.(solar[4:end, 1]), FT.(solar[4:end, 2]),
                                   extrapolation_bc = Line())
    end
    @warn "SOLAR_OUT not found; using unit solar transmission" SOLAR_OUT
    return ν -> one(FT)
end

function solar_F0(n_pol::Int, ν, solar_T)
    P = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-05 * π))
    F_vec = FT.(solar_T.(ν)) .* P
    F0 = zeros(FT, n_pol, length(ν))
    F0[1, :] .= F_vec
    return F0
end

function sif_matrix(n_pol::Int, ν, sif_on::Bool)
    SIF0 = zeros(FT, n_pol, length(ν))
    if sif_on
        ν_sif, jSIF = load_sif_spectrum()
        build_sif_source(SIF0, collect(ν), ν_sif, jSIF; pol_component = 1)
    end
    return SIF0
end

function make_rrs(model, i_band::Int, F0)
    n_pol = CoreRT.polarization_type(model).n
    ν = CoreRT.get_spec_bands(model)[i_band]
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    ν_center = FT(0.5) * (first(ν) + last(ν))
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν_center, effT)
    n_model_bands = length(CoreRT.get_spec_bands(model))

    rs = InelasticScattering.RRS(
        n2 = n2,
        o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = ones(FT, n_model_bands),
        ϖ_Cabannes = ones(FT, n_model_bands),
        ϖ_λ₁λ₀ = zeros(FT, 1),
        i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = argmin(abs.(ν .- mean(ν))),
        n_Raman = 0,
        iBand = [i_band],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, n_pol, length(ν)),
    )
    CoreRT.getRamanSSProp!(rs, 1e7 / mean(ν), ν)
    rs.n_Raman > 0 || error(
        "Band $i_band produced no usable Raman offsets " *
        "(nSpec=$(length(ν)), solve range=$(first(ν)):$(last(ν)) cm^-1, " *
        "WN_STEP=$WN_STEP, SHOULDER_CM=$SHOULDER_CM). " *
        "Use a wider/finer smoke grid; the production grid is WN_STEP=0.1 " *
        "with 250 cm^-1 shoulders."
    )
    rs.ϖ_Cabannes .= FT.(model.ϖ_Cabannes)
    return rs
end

function make_nors(model, i_band::Int, F0)
    n_pol = CoreRT.polarization_type(model).n
    n_model_bands = length(CoreRT.get_spec_bands(model))
    return InelasticScattering.noRS{FT}(
        fscattRayl = FT[1],
        ϖ_Cabannes = ones(FT, n_model_bands),
        bandSpecLim = UnitRange{Int}[],
        iBand = [i_band],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, n_pol, size(F0, 2)),
    )
end

function assert_rayleigh_spectral(model, i_band::Int)
    τ_col = vec(sum(Array(model.τ_rayl[i_band]), dims = 2))
    span = maximum(τ_col) - minimum(τ_col)
    if !(span > 0)
        error("Rayleigh optical depth appears spectrally constant in band $i_band")
    end
    return span
end

function prepare_output_axes(band_defs, solar_T)
    n_out = maximum(last(b.out_range) for b in band_defs)
    wn = zeros(FT, n_out)
    F0 = zeros(FT, n_out)
    SIF0 = zeros(FT, length(SIF_AXIS), n_out)
    band_index = zeros(Int32, n_out)
    n_pol = 3

    for (ib, b) in enumerate(band_defs)
        ν = b.solve_wn
        F_band = solar_F0(n_pol, ν, solar_T)
        for (isif, sif_on) in enumerate(SIF_AXIS)
            SIF_band = sif_matrix(n_pol, ν, sif_on)
            SIF0[isif, b.out_range] .= SIF_band[1, b.keep]
        end
        wn[b.out_range] .= ν[b.keep]
        F0[b.out_range] .= F_band[1, b.keep]
        band_index[b.out_range] .= Int32(ib)
    end
    return (; wn, wavelength_nm = FT.(1e7 ./ wn), F0, SIF0, band_index)
end

function init_netcdf(path::AbstractString, axes, profiles, band_defs)
    mkpath(dirname(path))
    if isfile(path)
        FORCE || error("Output exists: $path (set FORCE=1 to overwrite)")
        rm(path)
    end

    ds = NCDataset(path, "c")
    ds.attrib["title"] = "$(BAND_TITLE) Rayleigh/Cabannes/RRS Float32 GPU LUT"
    ds.attrib["source"] = "vSmartMOM.jl"
    ds.attrib["created"] = string(RUN_STARTED)
    ds.attrib["script"] = abspath(SCRIPT_PATH)
    ds.attrib["git_commit"] = try
        readchomp(`git -C $(REPO_ROOT) rev-parse HEAD`)
    catch
        "unknown"
    end
    ds.attrib["param_yaml"] = PARAM_YAML
    ds.attrib["band_tag"] = BAND_TAG
    ds.attrib["band_title"] = BAND_TITLE
    ds.attrib["profile_source_note"] = PROFILE_SOURCE_NOTE
    ds.attrib["fixed_molecules"] = join(FIXED_MOLECULES, ",")
    ds.attrib["variable_molecules"] = join(VARIABLE_MOLECULES, ",")
    ds.attrib["absorber_lut_files"] = join(absorption_lut_paths(), ",")
    ds.attrib["use_o2_absco"] = string(USE_O2_ABSCO)
    ds.attrib["o2_absco"] = USE_O2_ABSCO ? O2_ABSCO : "not used"
    ds.attrib["h2o_lut"] = isempty(H2O_LUT) ? "none; fallback artifact path" : H2O_LUT
    ds.attrib["cia_files"] = isempty(CIA_FILES) ? "none" : join(CIA_FILES, ",")
    ds.attrib["mtckd_file"] = isempty(MTCKD_FILE) ? "none" : MTCKD_FILE
    ds.attrib["solar_out"] = isfile(SOLAR_OUT) ? SOLAR_OUT : "unit solar transmission fallback"
    ds.attrib["float_type"] = "Float32"
    ds.attrib["architecture"] = ARCH
    ds.attrib["cuda_device"] = ARCH == "GPU" ? CUDA_DEVICE : "not used"
    ds.attrib["full_vza_axis_deg"] = join(FULL_VZA_AXIS, ",")
    ds.attrib["selected_vza_indices_1based"] = join(VZA_IDXS, ",")
    ds.attrib["selected_vza_axis_deg"] = join(VZA_AXIS, ",")
    ds.attrib["raz_axis_deg"] = join(RAZ_AXIS, ",")
    ds.attrib["camera_pairing"] =
        "Generated YAML repeats each selected VZA for every RAz; output is reshaped as vza x raz."
    ds.attrib["shoulder_cm-1"] = SHOULDER_CM
    ds.attrib["wn_step_cm-1"] = WN_STEP
    ds.attrib["profile_reduction"] = PROFILE_REDUCTION
    ds.attrib["chunked_lut"] = "true"
    ds.attrib["chunk_tag"] = CHUNK_TAG
    ds.attrib["full_albedo_axis"] = join(FULL_ALBEDO_AXIS, ",")
    ds.attrib["selected_albedo_indices_1based"] = join(ALBEDO_IDXS, ",")
    ds.attrib["full_sza_axis_deg"] = join(FULL_SZA_AXIS, ",")
    ds.attrib["selected_sza_indices_1based"] = join(SZA_IDXS, ",")
    ds.attrib["chunk_merge_note"] =
        "This file contains only selected albedo/sza coordinates; merge by coordinate " *
        "values or by the *_full_index variables, then validate the full grid."
    ds.attrib["note_tau_abs"] = absorption_note()

    n_psurf = length(PSURF_AXIS)
    n_sif = length(SIF_AXIS)
    n_alb = length(ALBEDO_AXIS)
    n_sza = length(SZA_AXIS)
    n_vza = length(VZA_AXIS)
    n_raz = length(RAZ_AXIS)
    n_stokes = 3
    n_wn = length(axes.wn)
    max_nlevel = maximum(length(p.p) for p in profiles)
    max_nlayer = maximum(length(p.T) for p in profiles)

    defDim(ds, "psurf", n_psurf)
    defDim(ds, "sif", n_sif)
    defDim(ds, "albedo", n_alb)
    defDim(ds, "sza", n_sza)
    defDim(ds, "vza", n_vza)
    defDim(ds, "raz", n_raz)
    defDim(ds, "full_albedo", length(FULL_ALBEDO_AXIS))
    defDim(ds, "full_sza", length(FULL_SZA_AXIS))
    defDim(ds, "full_vza", length(FULL_VZA_AXIS))
    defDim(ds, "stokes", n_stokes)
    defDim(ds, "wn", n_wn)
    defDim(ds, "band", length(band_defs))
    defDim(ds, "profile_level", max_nlevel)
    defDim(ds, "profile_layer", max_nlayer)

    v = defVar(ds, "psurf", Float32, ("psurf",))
    v.attrib["units"] = "hPa"
    v[:] = FT.(PSURF_AXIS)

    v = defVar(ds, "sif_on", Int32, ("sif",))
    v.attrib["description"] = "0 = SurfaceSIF off, 1 = SurfaceSIF on"
    v[:] = Int32.(SIF_AXIS)

    v = defVar(ds, "albedo", Float32, ("albedo",))
    v[:] = FT.(ALBEDO_AXIS)
    v = defVar(ds, "albedo_full_index", Int32, ("albedo",))
    v.attrib["description"] = "1-based index of each selected albedo in full_albedo"
    v[:] = Int32.(ALBEDO_IDXS)
    v = defVar(ds, "full_albedo", Float32, ("full_albedo",))
    v.attrib["description"] = "Full intended albedo axis before chunk selection"
    v[:] = FT.(FULL_ALBEDO_AXIS)

    v = defVar(ds, "sza", Float32, ("sza",))
    v.attrib["units"] = "degree"
    v[:] = FT.(SZA_AXIS)
    v = defVar(ds, "sza_full_index", Int32, ("sza",))
    v.attrib["description"] = "1-based index of each selected SZA in full_sza"
    v[:] = Int32.(SZA_IDXS)
    v = defVar(ds, "full_sza", Float32, ("full_sza",))
    v.attrib["units"] = "degree"
    v.attrib["description"] = "Full intended SZA axis before chunk selection"
    v[:] = FT.(FULL_SZA_AXIS)

    v = defVar(ds, "mu0", Float32, ("sza",))
    v.attrib["description"] = "cos(sza), equally spaced by construction"
    v[:] = FT.(MU0_AXIS)

    v = defVar(ds, "vza", Float32, ("vza",))
    v.attrib["units"] = "degree"
    v[:] = FT.(VZA_AXIS)
    v = defVar(ds, "vza_full_index", Int32, ("vza",))
    v.attrib["description"] = "1-based index of each selected VZA in full_vza"
    v[:] = Int32.(VZA_IDXS)
    v = defVar(ds, "full_vza", Float32, ("full_vza",))
    v.attrib["units"] = "degree"
    v.attrib["description"] = "Full intended VZA axis before chunk selection"
    v[:] = FT.(FULL_VZA_AXIS)

    v = defVar(ds, "raz", Float32, ("raz",))
    v.attrib["units"] = "degree"
    v.attrib["description"] = "Unique relative azimuth grid; 360 deg is intentionally omitted when using 0:15:345"
    v[:] = FT.(RAZ_AXIS)

    v = defVar(ds, "stokes_index", Int32, ("stokes",))
    v.attrib["description"] = "1=I, 2=Q, 3=U"
    v[:] = Int32[1, 2, 3]

    v = defVar(ds, "wn", Float32, ("wn",))
    v.attrib["units"] = "cm-1"
    v[:] = axes.wn

    v = defVar(ds, "wavelength_nm", Float32, ("wn",))
    v.attrib["units"] = "nm"
    v[:] = axes.wavelength_nm

    v = defVar(ds, "band_index", Int32, ("wn",))
    v.attrib["description"] = "1-based solve band that produced this output wavenumber"
    v[:] = axes.band_index

    v = defVar(ds, "F0", Float32, ("wn",))
    v.attrib["description"] = "Solar source flux used in SolarBeam, Stokes I component"
    v[:] = axes.F0

    v = defVar(ds, "SIF0", Float32, ("sif", "wn"))
    v.attrib["description"] = "SurfaceSIF source term, Stokes I component"
    v[:, :] = axes.SIF0

    v = defVar(ds, "band_solve_start", Float32, ("band",))
    v.attrib["units"] = "cm-1"
    v[:] = FT.([b.solve_start for b in band_defs])
    v = defVar(ds, "band_solve_stop", Float32, ("band",))
    v.attrib["units"] = "cm-1"
    v[:] = FT.([b.solve_stop for b in band_defs])
    v = defVar(ds, "band_output_start", Float32, ("band",))
    v.attrib["units"] = "cm-1"
    v[:] = FT.([b.output_start for b in band_defs])
    v = defVar(ds, "band_output_stop", Float32, ("band",))
    v.attrib["units"] = "cm-1"
    v[:] = FT.([b.output_stop for b in band_defs])

    profile_p = fill(FT(NaN), n_psurf, max_nlevel)
    profile_T = fill(FT(NaN), n_psurf, max_nlayer)
    profile_q = fill(FT(NaN), n_psurf, max_nlayer)
    nlevel = zeros(Int32, n_psurf)
    nlayer = zeros(Int32, n_psurf)
    for (ip, pcase) in enumerate(profiles)
        profile_p[ip, 1:length(pcase.p)] .= FT.(pcase.p)
        profile_T[ip, 1:length(pcase.T)] .= FT.(pcase.T)
        profile_q[ip, 1:length(pcase.q)] .= FT.(pcase.q)
        nlevel[ip] = Int32(length(pcase.p))
        nlayer[ip] = Int32(length(pcase.T))
    end
    v = defVar(ds, "profile_p", Float32, ("psurf", "profile_level"))
    v.attrib["units"] = "hPa"
    v[:, :] = profile_p
    v = defVar(ds, "profile_T", Float32, ("psurf", "profile_layer"))
    v.attrib["units"] = "K"
    v[:, :] = profile_T
    v = defVar(ds, "profile_q", Float32, ("psurf", "profile_layer"))
    v.attrib["units"] = "kg kg-1"
    v[:, :] = profile_q
    v = defVar(ds, "profile_nlevel", Int32, ("psurf",)); v[:] = nlevel
    v = defVar(ds, "profile_nlayer", Int32, ("psurf",)); v[:] = nlayer

    dims = ("psurf", "sif", "albedo", "sza", "vza", "raz", "stokes", "wn")
    rayl = defVar(ds, "stokes_rayleigh", Float32, dims)
    cab = defVar(ds, "stokes_cabannes", Float32, dims)
    rrs = defVar(ds, "stokes_rrs", Float32, dims)
    rayl.attrib["description"] = "TOA Stokes vector from noRS/Rayleigh run"
    cab.attrib["description"] = "TOA elastic Stokes vector from RRS/Cabannes run"
    rrs.attrib["description"] = "TOA inelastic rotational Raman Stokes vector from RRS run"

    tau_rayl = defVar(ds, "tau_rayl", Float32, ("psurf", "wn"))
    tau_rayl.attrib["description"] = "Column Rayleigh optical depth, summed over layers"
    tau_abs = defVar(ds, "tau_abs", Float32, ("psurf", "wn"))
    tau_abs.attrib["description"] = "Column gas absorption optical depth, summed over layers"
    tau_aer = defVar(ds, "tau_aer", Float32, ("psurf", "wn"))
    tau_aer.attrib["description"] = "Column aerosol optical depth; zero in this no-aerosol LUT"
    tau_aer[:, :] = zeros(FT, n_psurf, n_wn)

    varpi_cab = defVar(ds, "varpi_cabannes", Float32, ("psurf", "band"))
    varpi_cab.attrib["description"] = "Band-center Cabannes single-scattering albedo from model"
    effT = defVar(ds, "effective_temperature", Float32, ("psurf",))
    effT.attrib["description"] = "Dry-column weighted temperature after profile reduction"

    return ds, (; rayl, cab, rrs, tau_rayl, tau_abs, varpi_cab, effT)
end

function validate_resume_axis(ds, name::AbstractString, expected; atol = 1f-4)
    actual = collect(ds[name][:])
    length(actual) == length(expected) || error(
        "RESUME axis '$name' length mismatch: file has $(length(actual)), expected $(length(expected))")
    all(isapprox.(Float64.(actual), Float64.(expected); atol = atol, rtol = 0)) || error(
        "RESUME axis '$name' values do not match current run setup")
    return nothing
end

function validate_resume_int_axis(ds, name::AbstractString, expected)
    haskey(ds, name) || error("RESUME expected chunk variable '$name' in output file")
    actual = collect(ds[name][:])
    length(actual) == length(expected) || error(
        "RESUME axis '$name' length mismatch: file has $(length(actual)), expected $(length(expected))")
    collect(Int.(actual)) == collect(Int.(expected)) || error(
        "RESUME axis '$name' values do not match current chunk selection")
    return nothing
end

function open_resume_netcdf(path::AbstractString, axes, profiles)
    isfile(path) || error("RESUME=1 requested, but output does not exist: $path")
    ds = NCDataset(path, "a")
    try
        validate_resume_axis(ds, "psurf", FT.([p.psurf for p in profiles]))
        validate_resume_axis(ds, "albedo", FT.(ALBEDO_AXIS))
        validate_resume_axis(ds, "sza", FT.(SZA_AXIS))
        validate_resume_axis(ds, "vza", FT.(VZA_AXIS))
        validate_resume_axis(ds, "raz", FT.(RAZ_AXIS))
        validate_resume_int_axis(ds, "albedo_full_index", ALBEDO_IDXS)
        validate_resume_int_axis(ds, "sza_full_index", SZA_IDXS)
        validate_resume_int_axis(ds, "vza_full_index", VZA_IDXS)
        validate_resume_axis(ds, "wn", axes.wn; atol = 5f-3)
        vars = (;
            rayl = ds["stokes_rayleigh"],
            cab = ds["stokes_cabannes"],
            rrs = ds["stokes_rrs"],
            tau_rayl = ds["tau_rayl"],
            tau_abs = ds["tau_abs"],
            varpi_cab = ds["varpi_cabannes"],
            effT = ds["effective_temperature"],
        )
        return ds, vars
    catch
        close(ds)
        rethrow()
    end
end

function scene_complete(vars, ip::Int, isif::Int, ia::Int, isza::Int)
    fill_limit = FT(1e20)
    for v in (vars.rayl, vars.cab, vars.rrs)
        vals = Array(v[ip, isif, ia, isza, :, :, :, :])
        physical = isfinite.(vals) .& (abs.(vals) .< fill_limit)
        all(physical) || return false
        maximum(abs.(vals)) > zero(FT) || return false
    end
    return true
end

function print_plan(band_defs, profiles, axes)
    solve_points = sum(length(b.solve_wn) for b in band_defs)
    scene_count = length(PSURF_AXIS) * length(SIF_AXIS) * length(ALBEDO_AXIS) * length(SZA_AXIS)
    scale = solve_points / REF_NSPEC
    est_seconds = scene_count * REF_SECONDS_PER_SCENE * scale

    println("$BAND_TAG Raman LUT plan")
    println("  output NetCDF: $OUT_NC")
    println("  Raman parameter YAML: $PARAM_YAML")
    println("  band title: $BAND_TITLE")
    println("  architecture: $ARCH")
    ARCH == "GPU" && println("  CUDA device: $CUDA_DEVICE")
    println("  profile note: $PROFILE_SOURCE_NOTE")
    println("  psurf axis: $(join(PSURF_AXIS, ", ")) hPa")
    println("  SIF states: $(join(SIF_AXIS, ", "))")
    println("  chunk tag: $CHUNK_TAG")
    println("  full albedo count: $(length(FULL_ALBEDO_AXIS))")
    println("  selected albedo indices: $(join(ALBEDO_IDXS, ", "))")
    println("  selected albedo values: $(join(ALBEDO_AXIS, ", "))")
    dmu = length(MU0_AXIS) > 1 ? abs(MU0_AXIS[2] - MU0_AXIS[1]) : 0.0
    full_dmu = length(FULL_SZA_AXIS) > 1 ? abs(cosd(FULL_SZA_AXIS[2]) - cosd(FULL_SZA_AXIS[1])) : 0.0
    println("  full SZA count: $(length(FULL_SZA_AXIS)); max=$(maximum(FULL_SZA_AXIS)) deg; |full mu0 spacing|=$(full_dmu)")
    println("  selected SZA indices: $(join(SZA_IDXS, ", "))")
    println("  selected SZA values: $(join(SZA_AXIS, ", ")) deg")
    println("  selected SZA count: $(length(SZA_AXIS)); selected |mu0 spacing|=$(dmu)")
    println("  full VZA count: $(length(FULL_VZA_AXIS))")
    println("  selected VZA indices: $(join(VZA_IDXS, ", "))")
    println("  selected VZA values: $(join(VZA_AXIS, ", ")) deg")
    println("  RAz values: $(join(RAZ_AXIS, ", ")) deg")
    println("  output camera grid: $(length(VZA_AXIS)) VZA x $(length(RAZ_AXIS)) RAz = $(length(OBS_VZA_AXIS)) camera pairs")
    println("  bands: $(length(band_defs)); shoulder=$(SHOULDER_CM) cm^-1; step=$(WN_STEP) cm^-1")
    for (ib, b) in enumerate(band_defs)
        @printf("    band %d solve %.3f:%.3f:%.3f  output %.3f--%.3f  nsolve=%d nout=%d\n",
                ib, b.solve_start, WN_STEP, b.solve_stop,
                b.output_start, b.output_stop, length(b.solve_wn), length(b.keep))
    end
    println("  output wavelengths: $(minimum(axes.wavelength_nm))--$(maximum(axes.wavelength_nm)) nm")
    println("  output spectral points: $(length(axes.wn))")
    println("  scene count: $scene_count")
    @printf("  crude estimate: %.2f h total for selected PSURFS (ref %.0f s/scene at %.0f solve pts, scale %.2f)\n",
            est_seconds / 3600, REF_SECONDS_PER_SCENE, REF_NSPEC, scale)
    for ps in PSURF_AXIS
        n_ps_scene = length(SIF_AXIS) * length(ALBEDO_AXIS) * length(SZA_AXIS)
        @printf("    psurf %.0f hPa slice: %.2f h\n", ps,
                n_ps_scene * REF_SECONDS_PER_SCENE * scale / 3600)
    end
    println("  chunking note: run separate chunks on separate GPUs/servers; avoid parallel chunks on one saturated GPU")
    println("  fixed molecules: $(isempty(FIXED_MOLECULES) ? "none" : join(FIXED_MOLECULES, ", "))")
    println("  variable molecules: $(isempty(VARIABLE_MOLECULES) ? "none" : join(VARIABLE_MOLECULES, ", "))")
    println("  absorption LUTs: $(join(absorption_lut_paths(), ", "))")
    println("  O2 ABSCO: $(USE_O2_ABSCO ? O2_ABSCO : "not used")")
    println("  H2O LUT: $(isempty(H2O_LUT) ? "none/fallback" : H2O_LUT)")
    println("  CIA files: $(isempty(CIA_FILES) ? "none" : join(CIA_FILES, ", "))")
    println("  MT_CKD file: $(isempty(MTCKD_FILE) ? "none" : MTCKD_FILE)")
    println("  solar: $(isfile(SOLAR_OUT) ? SOLAR_OUT : "unit fallback")")
    println("  DRY_RUN=$DRY_RUN")
    println("  RESUME=$RESUME")
end

function run!()
    length(WL_EDGES_NM) == NBANDS + 1 || error("Bad WL_EDGES_NM")
    NBANDS <= 2 || error("This workflow is intended to use no more than 2 bands; got $NBANDS")
    isfile(PARAM_YAML) || error("PARAM_YAML not found: $PARAM_YAML")
    validate_absorption_files!()

    band_defs = build_band_defs()
    base_profile = load_profile_from_yaml(PARAM_YAML)
    vmr = load_vmr_from_yaml(PARAM_YAML)
    profiles = [truncate_profile(base_profile, ps) for ps in PSURF_AXIS]
    solar_T = solar_interpolator()
    axes = prepare_output_axes(band_defs, solar_T)
    print_plan(band_defs, profiles, axes)
    DRY_RUN && return nothing

    configure_architecture!()
    ds, vars = RESUME ? open_resume_netcdf(OUT_NC, axes, profiles) :
                        init_netcdf(OUT_NC, axes, profiles, band_defs)
    config_dir = joinpath(dirname(OUT_NC), "config_" * splitext(basename(OUT_NC))[1])
    mkpath(config_dir)

    try
        for (ip, profile) in enumerate(profiles)
            for (isza, sza) in enumerate(SZA_AXIS)
                yaml_path = joinpath(config_dir,
                    "$(lowercase(BAND_TAG))_psurf$(label_float(profile.psurf; digits=0))_sza$(label_float(sza; digits=3)).yaml")
                write_case_yaml(yaml_path, profile, sza, band_defs, vmr)
                params = parameters_from_yaml(yaml_path)
                params.float_type === Float32 ||
                    error("Generated YAML did not resolve to Float32")
                model0 = model_from_parameters(params)

                if isza == 1
                    vars.varpi_cab[ip, :] = FT.(model0.ϖ_Cabannes)
                    vars.effT[ip] = FT((model0.profile.vcd_dry' * model0.profile.T) / sum(model0.profile.vcd_dry))
                    for (ib, b) in enumerate(band_defs)
                        assert_rayleigh_spectral(model0, ib)
                        vars.tau_rayl[ip, b.out_range] =
                            FT.(vec(sum(Array(model0.τ_rayl[ib]), dims = 2))[b.keep])
                        vars.tau_abs[ip, b.out_range] =
                            FT.(vec(sum(Array(model0.τ_abs[ib]), dims = 2))[b.keep])
                    end
                end

                for (ia, albedo) in enumerate(ALBEDO_AXIS)
                    model = with_lambertian_albedo(model0, albedo)
                    for (isif, sif_on) in enumerate(SIF_AXIS)
                        if RESUME && scene_complete(vars, ip, isif, ia, isza)
                            @info "Skipping completed RT scene" psurf = profile.psurf sza albedo sif_on
                            continue
                        end
                        for (ib, b) in enumerate(band_defs)
                            ν = CoreRT.get_spec_bands(model)[ib]
                            n_pol = CoreRT.polarization_type(model).n
                            F0 = solar_F0(n_pol, ν, solar_T)
                            SIF0 = sif_matrix(n_pol, ν, sif_on)
                            sources = CoreRT.SolarBeam(F₀ = F0) +
                                      CoreRT.SurfaceSIF(SIF₀ = SIF0)
                            rs_rrs = make_rrs(model, ib, F0)
                            rs_nors = make_nors(model, ib, F0)

                            @info "RT scene" psurf = profile.psurf sza albedo sif_on ib NBANDS nSpec = length(ν) nRaman = rs_rrs.n_Raman
                            t0 = time()
                            rrs_result = CoreRT.rt_run_test(rs_rrs, model, ib; sources)
                            nors_result = CoreRT.rt_run_test(rs_nors, model, ib; sources)
                            dt = round(time() - t0; digits = 2)
                            @info "RT scene done" psurf = profile.psurf sza albedo sif_on ib seconds = dt

                            Rcab = Array(rrs_result[1])
                            Rrrs = Array(rrs_result[3])
                            Rray = Array(nors_result[1])
                            keep = b.keep
                            out = b.out_range
                            expected_cameras = length(VZA_AXIS) * length(RAZ_AXIS)
                            size(Rcab, 1) == expected_cameras ||
                                error("Expected $expected_cameras camera rows, got $(size(Rcab, 1))")
                            for iv in eachindex(VZA_AXIS), ir in eachindex(RAZ_AXIS)
                                icam = (iv - 1) * length(RAZ_AXIS) + ir
                                vars.cab[ip, isif, ia, isza, iv, ir, :, out] =
                                    FT.(Rcab[icam, 1:3, keep])
                                vars.rrs[ip, isif, ia, isza, iv, ir, :, out] =
                                    FT.(Rrrs[icam, 1:3, keep])
                                vars.rayl[ip, isif, ia, isza, iv, ir, :, out] =
                                    FT.(Rray[icam, 1:3, keep])
                            end
                        end
                    end
                    if ARCH == "GPU"
                        CUDA = ensure_cuda()
                        GC.gc()
                        CUDA.reclaim()
                    end
                    NCDatasets.sync(ds)
                end
            end
        end
    finally
        close(ds)
    end

    @info "Wrote NetCDF" OUT_NC
    return OUT_NC
end

if abspath(PROGRAM_FILE) == @__FILE__
    run!()
end
