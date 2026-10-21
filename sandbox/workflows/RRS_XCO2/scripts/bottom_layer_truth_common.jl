module BottomLayerTruthCommon

using DelimitedFiles
using Printf
using SHA
using UUIDs

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common

export BottomTruthState, SURFACES, AEROSOL_CASES, SIF_CASES,
       SIF_DEFINITION_VERSION, SIF_ANGULAR_INTEGRAL_760,
       SIF_RADIANCE_760, SIF760_NATIVE, MSIF_NATIVE,
       BOTTOM_CO2_PPM, BACKGROUND_CO2_PPM, BOTTOM_LAYER_INDEX,
       BOTTOM_DRY_COLUMN_FRACTION, CAMPAIGN_ROOT, TRUTH_ROOT, STATE_FILE,
       FULL_COLUMN_TRUTH_ROOT, FULL_COLUMN_STATE_FILE,
       state_index, old_control_index, paired_no_sif_index,
       column_xco2_ppm, read_bottom_states, write_campaign_tables!,
       state_table_sha256

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CAMPAIGN_ROOT = get(ENV, "BOTTOM_XCO2_CAMPAIGN_ROOT",
    joinpath(ROOT, "bottom_layer_XCO2_retrievals"))
const TRUTH_ROOT = joinpath(CAMPAIGN_ROOT, "truth")
const STATE_FILE = joinpath(TRUTH_ROOT, "true_states.dat")
const FULL_COLUMN_TRUTH_ROOT = get(ENV, "FULL_COLUMN_TRUTH_ROOT",
    joinpath(ROOT, "truth_map"))
const FULL_COLUMN_STATE_FILE = joinpath(FULL_COLUMN_TRUTH_ROOT,
                                        "true_states.dat")

const SURFACES = ("urban", "rural", "desert", "forest")
const AEROSOL_CASES = ("none", "aod760_0p28")
const SIF_CASES = ("off", RRSXCO2Common.SIF_CASE_ON)
const SIF_DEFINITION_VERSION = RRSXCO2Common.SIF_DEFINITION_VERSION
const SIF_ANGULAR_INTEGRAL_760 =
    RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760
const SIF_RADIANCE_760 = RRSXCO2Common.SIF_RADIANCE_760
const _CAMPAIGN_SIF_STATE = RRSXCO2Common.campaign_sif_state()
const SIF760_NATIVE = _CAMPAIGN_SIF_STATE.SIF760
const MSIF_NATIVE = _CAMPAIGN_SIF_STATE.mSIF
const BOTTOM_CO2_PPM = (360, 380, 400, 420, 440)
const BACKGROUND_CO2_PPM = 400.0
const BOTTOM_LAYER_INDEX = 16

# Derived from the exact Float32 production p/T/q profile after truncation to
# 1000 hPa and reduction to 16 layers.  The profile is ordered TOA -> BOA, so
# this is vcd_dry[end] / sum(vcd_dry).  Keep the extra digits: the state table
# records the nominal column average implied by the production dry-air weights.
# The solver stores VMRs as Float32, so its realized weighted mean can differ
# from this intended-ppm value by less than 1e-5 ppm.
const BOTTOM_DRY_COLUMN_FRACTION = 0.06229780735737046

struct BottomTruthState
    index::Int
    surface_index::Int
    surface::String
    aerosol_index::Int
    aerosol_case::String
    sif_index::Int
    sif_case::String
    bottom_co2_index::Int
    background_co2_ppm::Float64
    bottom_layer_index::Int
    bottom_co2_ppm::Float64
    xco2_ppm::Float64
    psurf_hpa::Float64
    sza_deg::Float64
    vza_deg::Float64
    relative_azimuth_deg::Float64
    aod550::NTuple{3,Float64}
    sif_angular_integral760::Float64
    sif760::Float64
    msif::Float64
    coeff::NTuple{3,NTuple{3,Float64}}
end

"State index with surface -> aerosol -> SIF -> bottom-CO2 nesting."
state_index(surface_index::Integer, aerosol_index::Integer,
            sif_index::Integer, bottom_co2_index::Integer) =
    20 * (surface_index - 1) + 10 * (aerosol_index - 1) +
     5 * (sif_index - 1) + bottom_co2_index

"Matching uniform-400-ppm state in the accepted 64-scene truth map."
old_control_index(surface_index::Integer, aerosol_index::Integer,
                  sif_index::Integer) =
    16 * (surface_index - 1) + 8 * (aerosol_index - 1) +
     4 * (sif_index - 1) + 2

old_control_index(state::BottomTruthState) = old_control_index(
    state.surface_index, state.aerosol_index, state.sif_index)

"New no-SIF state with the same surface, aerosol, and bottom CO2 profile."
paired_no_sif_index(state::BottomTruthState) = state_index(
    state.surface_index, state.aerosol_index, 1, state.bottom_co2_index)

column_xco2_ppm(bottom_co2_ppm::Real) =
    BACKGROUND_CO2_PPM + BOTTOM_DRY_COLUMN_FRACTION *
    (Float64(bottom_co2_ppm) - BACKGROUND_CO2_PPM)

function _full_column_rows()
    isfile(FULL_COLUMN_STATE_FILE) || error(
        "missing accepted full-column state table: $FULL_COLUMN_STATE_FILE")
    raw = readdlm(FULL_COLUMN_STATE_FILE; comments=true)
    size(raw) == (64, 27) || error(
        "expected a 64x27 accepted state table, found $(size(raw))")
    result = Dict{Int,Vector{Any}}()
    for row in eachrow(raw)
        values = Any[item for item in row]
        index = Int(values[1])
        haskey(result, index) && error(
            "duplicate state $index in $FULL_COLUMN_STATE_FILE")
        result[index] = values
    end
    length(result) == 64 || error("accepted state table is incomplete")
    return result
end

function _write_state_table(path)
    old_states = _full_column_rows()
    open(path, "w") do io
        println(io, "# Bottom-layer CO2 truth table; state order: surface -> aerosol -> SIF -> bottom-layer CO2")
        println(io, "# Layers 1:15 are 400 ppm; layer 16 (TOA-to-BOA ordering) takes the tabulated bottom_co2_ppm.")
        println(io, "# xco2_ppm is the nominal dry-column average, not the bottom-layer VMR; bottom dry-column fraction = $(BOTTOM_DRY_COLUMN_FRACTION).")
        println(io, "# SIF definition v$(SIF_DEFINITION_VERSION): 2pi*L_lambda(760 nm)=0.5 mW m^-2 nm^-1; every isotropic BOA stream has L_lambda=0.5/(2pi) per sr.")
        println(io, "# index surface_index surface aerosol_index aerosol_case sif_index sif_case bottom_co2_index background_co2_ppm bottom_layer_index bottom_co2_ppm xco2_ppm psurf_hpa sza_deg vza_deg relative_azimuth_deg sulfate_aod550 organic_aod550 stratospheric_aod550 sif_angular_integral760 SIF760 mSIF o2a_P0 o2a_P1 o2a_P2 weak_P0 weak_P1 weak_P2 strong_P0 strong_P1 strong_P2")
        for (isurf, surface) in enumerate(SURFACES),
            (iaer, aerosol) in enumerate(AEROSOL_CASES),
            (isif, sif_case) in enumerate(SIF_CASES),
            (ico2, bottom_co2) in enumerate(BOTTOM_CO2_PPM)
            index = state_index(isurf, iaer, isif, ico2)
            old = old_states[old_control_index(isurf, iaer, isif)]
            Int(old[2]) == isurf && String(old[3]) == surface &&
                Int(old[4]) == iaer && String(old[5]) == aerosol &&
                Int(old[6]) == isif && String(old[7]) == sif_case &&
                Int(old[9]) == 400 || error(
                    "accepted control metadata is inconsistent for new state $index")
            aod = Tuple(Float64.(old[13:15]))
            sc = Tuple(Float64.(old[16:18]))
            c1 = Tuple(Float64.(old[19:21]))
            c2 = Tuple(Float64.(old[22:24]))
            c3 = Tuple(Float64.(old[25:27]))
            @printf(io,
                "%03d %d %s %d %s %d %s %d %.1f %d %.1f %.12f 1000.0 30.0 0.0 0.0 %.12f %.12f %.12f %.8e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                index, isurf, surface, iaer, aerosol, isif, sif_case,
                ico2, BACKGROUND_CO2_PPM, BOTTOM_LAYER_INDEX,
                Float64(bottom_co2), column_xco2_ppm(bottom_co2), aod...,
                sc..., c1..., c2..., c3...)
        end
    end
    return path
end

function _write_reuse_map(path)
    open(path, "w") do io
        println(io, "# Accepted uniform-400 full-column source -> new bottom-layer 400-ppm control")
        println(io, "# old_state new_state surface aerosol_case sif_case")
        for (isurf, surface) in enumerate(SURFACES),
            (iaer, aerosol) in enumerate(AEROSOL_CASES),
            (isif, sif_case) in enumerate(SIF_CASES)
            old = old_control_index(isurf, iaer, isif)
            new = state_index(isurf, iaer, isif, 3)
            @printf(io, "%03d %03d %s %s %s\n",
                    old, new, surface, aerosol, sif_case)
        end
    end
    return path
end

function _write_component_catalog(path)
    shared = joinpath(FULL_COLUMN_TRUTH_ROOT, "scene_components.dat")
    isfile(shared) || error("missing shared component catalog: $shared")
    open(path, "w") do io
        println(io, "# Bottom-layer-XCO2 campaign component catalog")
        println(io, "# The aerosol and surface definitions below are copied from the accepted full-column campaign; relative documentation paths are rebased for this directory.")
        println(io, "# For exact component reuse, true_states.dat is authoritative: its shared aerosol, SIF, and surface inputs are cloned from the corrected full-column true_states.dat rows.")
        println(io, "# SIF-on convention: 2pi*L_lambda(760 nm)=0.5 mW m^-2 nm^-1, hence L_lambda(760)=0.5/(2pi) per sr; this is not a wavelength-integrated spectrum.")
        println(io, "# The copied catalog may show additional descriptive digits; do not use those digits to regenerate solver inputs.")
        println(io)
        println(io, "[CO2_PROFILES]")
        println(io, "# case_index background_ppm bottom_layer_index bottom_ppm column_XCO2_ppm")
        for (index, bottom) in enumerate(BOTTOM_CO2_PPM)
            @printf(io, "%d %.1f %d %.1f %.12f\n", index,
                    BACKGROUND_CO2_PPM, BOTTOM_LAYER_INDEX, Float64(bottom),
                    column_xco2_ppm(bottom))
        end
        println(io, "# Layers are ordered TOA-to-BOA. Layers 1:15 remain at 400 ppm; only layer 16 changes.")
        println(io, "# Current layer-16 dry-column fraction: $(BOTTOM_DRY_COLUMN_FRACTION)")
        println(io)
        shared_text = read(shared, String)
        shared_text = replace(
            shared_text,
            "../truth_map_aerosols/" => "../../truth_map_aerosols/",
            "../surface_albedos/" => "../../surface_albedos/",
        )
        write(io, shared_text)
    end
    return path
end

function _publish_text!(writer::Function, destination::AbstractString)
    mkpath(dirname(destination))
    temporary = destination * ".tmp.$(uuid4())"
    ispath(temporary) && error("temporary output already exists: $temporary")
    writer(temporary)
    if isfile(destination)
        if read(destination) == read(temporary)
            rm(temporary)
            return destination
        end
        rm(temporary)
        error("refusing to replace a different campaign table: $destination")
    end
    try
        mv(temporary, destination)
    catch error_value
        # Two disjoint host launchers may bootstrap the static tables at the
        # same instant.  Treat a concurrently published byte-identical file as
        # success, but never replace a different one.
        if isfile(destination) && read(destination) == read(temporary)
            rm(temporary)
            return destination
        end
        rethrow(error_value)
    end
    return destination
end

function write_campaign_tables!(root::AbstractString=TRUTH_ROOT)
    state_file = joinpath(root, "true_states.dat")
    reuse_file = joinpath(root, "control_reuse_map.dat")
    component_file = joinpath(root, "scene_components.dat")
    _publish_text!(_write_state_table, state_file)
    _publish_text!(_write_reuse_map, reuse_file)
    _publish_text!(_write_component_catalog, component_file)
    return (; state_file, reuse_file, component_file)
end

function _state_from_row(row)
    return BottomTruthState(
        Int(row[1]), Int(row[2]), String(row[3]), Int(row[4]), String(row[5]),
        Int(row[6]), String(row[7]), Int(row[8]), Float64(row[9]), Int(row[10]),
        Float64(row[11]), Float64(row[12]), Float64(row[13]), Float64(row[14]),
        Float64(row[15]), Float64(row[16]),
        (Float64(row[17]), Float64(row[18]), Float64(row[19])),
        Float64(row[20]), Float64(row[21]), Float64(row[22]),
        ((Float64(row[23]), Float64(row[24]), Float64(row[25])),
         (Float64(row[26]), Float64(row[27]), Float64(row[28])),
         (Float64(row[29]), Float64(row[30]), Float64(row[31]))),
    )
end

function _validate_states(states::Vector{BottomTruthState})
    length(states) == 80 || error("expected 80 bottom-layer truth states, found $(length(states))")
    for (expected, state) in enumerate(states)
        state.index == expected || error(
            "truth table row $expected carries state index $(state.index)")
        state.index == state_index(state.surface_index, state.aerosol_index,
                                   state.sif_index, state.bottom_co2_index) ||
            error("state $(state.index) violates the campaign indexing rule")
        state.surface == SURFACES[state.surface_index] || error(
            "state $(state.index) has inconsistent surface metadata")
        state.aerosol_case == AEROSOL_CASES[state.aerosol_index] || error(
            "state $(state.index) has inconsistent aerosol metadata")
        state.sif_case == SIF_CASES[state.sif_index] || error(
            "state $(state.index) has inconsistent SIF metadata")
        state.bottom_co2_ppm == BOTTOM_CO2_PPM[state.bottom_co2_index] ||
            error("state $(state.index) has inconsistent bottom-layer CO2")
        state.background_co2_ppm == BACKGROUND_CO2_PPM || error(
            "state $(state.index) has a non-canonical CO2 background")
        state.bottom_layer_index == BOTTOM_LAYER_INDEX || error(
            "state $(state.index) has a non-canonical bottom-layer index")
        isapprox(state.xco2_ppm, column_xco2_ppm(state.bottom_co2_ppm);
                 atol=5e-10, rtol=0) || error(
            "state $(state.index) has an inconsistent column XCO2")
        state.psurf_hpa == 1000.0 && state.sza_deg == 30.0 &&
            state.vza_deg == 0.0 && state.relative_azimuth_deg == 0.0 ||
            error("state $(state.index) has non-canonical pressure or geometry")
        if state.sif_index == 1
            all(iszero, (state.sif_angular_integral760,
                         state.sif760, state.msif)) || error(
                "state $(state.index) has nonzero SIF in the off case")
        else
            state.sif_angular_integral760 == SIF_ANGULAR_INTEGRAL_760 ||
                error("state $(state.index) has the wrong SIF angular integral")
            isapprox(state.sif760, SIF760_NATIVE; atol=5e-15, rtol=0) ||
                error("state $(state.index) has the wrong SIF760")
            isapprox(state.msif, MSIF_NATIVE; atol=5e-16, rtol=0) ||
                error("state $(state.index) has the wrong mSIF")
            isapprox(2π * state.sif760 * 1.0e7 / 760.0^2,
                     SIF_ANGULAR_INTEGRAL_760; atol=5e-13, rtol=0) ||
                error("state $(state.index) violates 2pi*L_lambda(760)=0.5")
        end
    end
    return states
end

function read_bottom_states(path::AbstractString=STATE_FILE)
    isfile(path) || error(
        "missing $path; run generate_bottom_layer_truth_states.jl first")
    raw = readdlm(path; comments=true)
    size(raw, 2) == 31 || error(
        "expected 31 columns in $path, found $(size(raw, 2))")
    return _validate_states([_state_from_row(row) for row in eachrow(raw)])
end

state_table_sha256(path::AbstractString=STATE_FILE) = bytes2hex(sha256(read(path)))

end # module
