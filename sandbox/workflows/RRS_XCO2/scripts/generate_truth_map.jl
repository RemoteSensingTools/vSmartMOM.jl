#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using CUDA
using DataInterpolations
using Dates
using DelimitedFiles
using Distributions
using NCDatasets
using Printf
using Statistics
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TRUTH_ROOT = joinpath(ROOT, "truth_map")
const OUT = get(ENV, "TRUTH_OUT", TRUTH_ROOT)
const STATE_FILE = joinpath(TRUTH_ROOT, "true_states.dat")
const DEVICE = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
const FIRST_STATE = parse(Int, get(ENV, "FIRST_STATE", "1"))
const LAST_STATE = parse(Int, get(ENV, "LAST_STATE", "64"))
const AEROSOL_CASE_FILTER = lowercase(get(ENV, "AEROSOL_CASE_FILTER", "all"))
const SIF_CASE_FILTER = lowercase(get(ENV, "SIF_CASE_FILTER", "all"))
const FORCE = get(ENV, "FORCE", "0") == "1"
const NLAYERS = parse(Int, get(ENV, "NLAYERS", "16"))
const AEROSOL_NSTREAMS = parse(Int, get(ENV, "AEROSOL_NSTREAMS", "9"))
const FT = get(ENV, "RRS_XCO2_FLOAT_TYPE", "Float32") == "Float64" ?
    Float64 : Float32
const SZA_DEG = FT(parse(Float64, get(ENV, "TRUTH_SZA_DEG", "30")))
const VZA_DEG = FT(parse(Float64, get(ENV, "TRUTH_VZA_DEG", "0")))
const RELATIVE_AZIMUTH_DEG = FT(parse(Float64,
    get(ENV, "TRUTH_RELATIVE_AZIMUTH_DEG", "0")))
const SHOULDER_CM = FT(parse(Float64, get(ENV, "RAMAN_SHOULDER_CM", "234")))
const STEP = FT(0.1)
const SOLAR_OUT = get(ENV, "SOLAR_OUT",
    joinpath(homedir(), "Raman_misc", "workflows", "worktrees",
             "uni_vSmartMOM_sanghavi_2025-03-18", "src", "SolarModel", "solar.out"))
const RAMAN_CONFIG = joinpath(pkgdir(vSmartMOM), "workflows", "OCORaman",
                              "config", "paramsRamanLUT_O2A.yaml")
const SURFACES = ("urban", "rural", "desert", "forest")
const XCO2_PPM = (380, 400, 420, 440)

AEROSOL_CASE_FILTER in ("all", "none", "aerosol") || error(
    "AEROSOL_CASE_FILTER must be all, none, or aerosol")
SIF_CASE_FILTER in ("all", "off", "on") || error(
    "SIF_CASE_FILTER must be all, off, or on")
zero(FT) <= SZA_DEG < FT(90) || error("TRUTH_SZA_DEG must be in [0, 90)")
zero(FT) <= VZA_DEG < FT(90) || error("TRUTH_VZA_DEG must be in [0, 90)")
zero(FT) <= RELATIVE_AZIMUTH_DEG <= FT(360) || error(
    "TRUTH_RELATIVE_AZIMUTH_DEG must be in [0, 360]")

struct TruthState
    index::Int
    surface_index::Int
    surface::String
    aerosol_index::Int
    aerosol_case::String
    sif_index::Int
    sif_case::String
    xco2_index::Int
    xco2_ppm::Int
    aod550::NTuple{3,Float64}
    sif_total::Float64
    sif760::Float64
    msif::Float64
    coeff::NTuple{3,NTuple{3,Float64}}
end

function read_states(path=STATE_FILE)
    raw = readdlm(path; comments=true)
    states = TruthState[]
    for r in eachrow(raw)
        push!(states, TruthState(Int(r[1]), Int(r[2]), String(r[3]), Int(r[4]),
            String(r[5]), Int(r[6]), String(r[7]), Int(r[8]), Int(r[9]),
            (Float64(r[13]), Float64(r[14]), Float64(r[15])), Float64(r[16]),
            Float64(r[17]), Float64(r[18]),
            ((Float64(r[19]), Float64(r[20]), Float64(r[21])),
             (Float64(r[22]), Float64(r[23]), Float64(r[24])),
             (Float64(r[25]), Float64(r[26]), Float64(r[27])))))
    end
    length(states) == 64 || error("Expected 64 truth states, found $(length(states))")
    return states
end

function selected_states(path=STATE_FILE)
    states = read_states(path)
    1 <= FIRST_STATE <= LAST_STATE <= length(states) || error(
        "invalid state interval $FIRST_STATE:$LAST_STATE")
    selected = states[FIRST_STATE:LAST_STATE]
    if AEROSOL_CASE_FILTER == "none"
        filter!(state -> !any(>(0), state.aod550), selected)
    elseif AEROSOL_CASE_FILTER == "aerosol"
        filter!(state -> any(>(0), state.aod550), selected)
    end
    if SIF_CASE_FILTER == "off"
        filter!(state -> iszero(state.sif_total), selected)
    elseif SIF_CASE_FILTER == "on"
        filter!(state -> !iszero(state.sif_total), selected)
    end
    isempty(selected) && error(
        "state selection is empty after AEROSOL_CASE_FILTER=" *
        "$AEROSOL_CASE_FILTER and SIF_CASE_FILTER=$SIF_CASE_FILTER")
    return selected
end

function truncate_profile!(params, psurf=1000.0)
    return RRSXCO2Common.truncate_profile_to_surface!(params, psurf)
end

function select_band!(params, iband, ν, surface)
    params.spec_bands = [FT.(ν)]
    params.brdf = [surface]
    ap = params.absorption_params
    ap.fixed_molecules = [ap.fixed_molecules[iband]]
    ap.variable_molecules = [ap.variable_molecules[iband]]
    ap.luts = [ap.luts[iband]]
    ap.h2o_lut = [ap.h2o_lut[iband]]
    if length(ap.cia_files) > 1
        ap.cia_files = [ap.cia_files[iband]]
        ap.cia_reference_codes = [ap.cia_reference_codes[iband]]
        ap.cia_negative_policies = [ap.cia_negative_policies[iband]]
    end
    return params
end

function set_common!(params, state)
    params.float_type = FT
    params.architecture = GPU()
    params.sza = SZA_DEG
    params.vza = FT[VZA_DEG]
    params.vaz = FT[RELATIVE_AZIMUTH_DEG]
    aerosol_active = any(>(0), state.aod550)
    if aerosol_active
        # Requested production aerosol resolution. δBGE's l_max is an
        # inclusive truncation parameter. The stream cap follows the public
        # 2N-1 contract; Nquad remains a derived operator dimension and is not
        # set directly. AEROSOL_NSTREAMS permits isolated memory benchmarks
        # without changing the nine-stream production default.
        params.nstreams = AEROSOL_NSTREAMS
        params.l_trunc = 2 * AEROSOL_NSTREAMS - 1
        params.stream_l_cap = 2 * AEROSOL_NSTREAMS - 1
        params.legacy_l_cap_override = nothing
        params.truncation = vSmartMOM.Scattering.δBGE{FT}(16, zero(FT))
    else
        # Rayleigh has exact Fourier support through m=2; retain the existing
        # five-stream resolution for aerosol-free scenes.
        params.nstreams = 5
        params.l_trunc = 9
        params.stream_l_cap = 9
        params.legacy_l_cap_override = nothing
        params.truncation = vSmartMOM.Scattering.NoTruncation()
    end
    params.absorption_params.vmr["CO2"] = FT(state.xco2_ppm * 1e-6)
    RRSXCO2Common.prepare_shared_profile!(params; psurf=1000.0, nlayers=NLAYERS)
    for (aer, τ) in zip(params.scattering_params.rt_aerosols, state.aod550)
        aer.τ_ref = FT(τ)
    end
    return params
end

"Convert the experiment aerosol configuration to the solver precision."
function scattering_as_ft32(sp)
    convert_dist(d::LogNormal) = LogNormal(FT.(Distributions.params(d))...)
    convert_dist(d::Normal) = Normal(FT.(Distributions.params(d))...)
    aerosols = CoreRT.RT_Aerosol{FT}[]
    for rt in sp.rt_aerosols
        mie = rt.aerosol
        mie32 = vSmartMOM.Scattering.Aerosol(
            convert_dist(mie.size_distribution), FT(mie.nᵣ), FT(mie.nᵢ))
        profile32 = convert_dist(rt.profile)
        if rt.phase_function === nothing
            push!(aerosols, CoreRT.RT_Aerosol(mie32, FT(rt.τ_ref), profile32))
        else
            push!(aerosols, CoreRT.RT_Aerosol(
                mie32, FT(rt.τ_ref), profile32, rt.phase_function; ϖ=FT(rt.ϖ)))
        end
    end
    return CoreRT.ScatteringParameters(
        aerosols, FT(sp.r_max), sp.nquad_radius, FT(sp.λ_ref),
        Complex{FT}(sp.n_ref), sp.decomp_type)
end

"""
    band_surface(coeff, band_ν, solve_ν)

Evaluate the three-term Legendre albedo model in the coordinate of the complete
output band, then return its values on a (possibly chunked or shoulder-expanded)
solve grid.  Defining `x = -1:1` from `band_ν` *before* chunking is essential:
`LambertianSurfaceLegendre` defines that coordinate from the grid it receives,
so constructing one independently for every chunk would reset the surface
polynomial at every chunk boundary.

O₂ Raman shoulder points lie outside `band_ν`; evaluating the same polynomial
there provides the intended smooth extrapolation without changing its
normalization.  An explicit spectrum prevents the surface implementation from
introducing any second, local normalization.
"""
function band_surface(coeff, band_ν, solve_ν)
    length(band_ν) > 1 || error("the complete surface band needs at least two points")
    span = last(band_ν) - first(band_ν)
    iszero(span) && error("the complete surface band has zero wavenumber span")
    x = 2 .* (solve_ν .- first(band_ν)) ./ span .- 1
    ρ = coeff[1] .+ coeff[2] .* x .+
        coeff[3] .* (3 .* x.^2 .- 1) ./ 2
    return CoreRT.LambertianSurfaceSpectrum(FT.(ρ))
end

# Compatibility name for small analysis scripts.  The second argument must be
# the complete band, not the local retained core of a spectral chunk.
transformed_surface(coeff, band_ν, solve_ν) = band_surface(coeff, band_ν, solve_ν)

function output_grids()
    return RRSXCO2Common.basis_grids(FT)
end

function o2_solve_grid(output_ν)
    return RRSXCO2Common.raman_solve_grid(
        output_ν; shoulder_cm=SHOULDER_CM, step_cm=STEP)
end

function solar_interpolator()
    isfile(SOLAR_OUT) || error("Required high-resolution solar spectrum missing: $SOLAR_OUT")
    solar = solar_transmission_from_file(SOLAR_OUT)
    # DataInterpolations uses the (dependent value, coordinate) argument order.
    return LinearInterpolation(FT.(solar[4:end, 2]), FT.(solar[4:end, 1]);
                               extrapolation=ExtrapolationType.Linear)
end

function source_set(ν, sif_on, solar_T)
    nPol = 3
    P = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-5 * π))
    F0 = zeros(FT, nPol, length(ν))
    F0[1, :] .= FT.(solar_T.(ν)) .* P
    SIF0 = zeros(FT, nPol, length(ν))
    if sif_on
        state = sif_reference_state(total_sif=0.5, reference_wavelength_nm=760)
        build_sif_source(SIF0, collect(ν), state.ν, FT(π) .* state.spectrum)
    end
    return F0, SolarBeam(F₀=F0) + SurfaceSIF(SIF₀=SIF0)
end

function make_rrs(model, F0)
    ν = CoreRT.get_spec_bands(model)[1]
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(mean(ν), effT)
    rs = InelasticScattering.RRS(n2=n2, o2=o2,
        greek_raman=InelasticScattering.GreekCoefs(FT[1],FT[1],FT[1],FT[1],FT[1],FT[1]),
        fscattRayl=FT[1], ϖ_Cabannes=FT[1], ϖ_λ₁λ₀=zeros(FT,1),
        i_λ₁λ₀=zeros(Int,1), Z⁻⁺_λ₁λ₀=zeros(FT,1,1), Z⁺⁺_λ₁λ₀=zeros(FT,1,1),
        i_ref=argmin(abs.(ν .- mean(ν))), n_Raman=0, iBand=[1],
        F₀=copy(F0), SIF₀=zeros(FT,3,length(ν)))
    CoreRT.getRamanSSProp!(rs, 1e7/mean(ν), ν)
    rs.ϖ_Cabannes .= FT.(model.ϖ_Cabannes)
    return rs
end

make_nors(F0) = InelasticScattering.noRS{FT}(fscattRayl=FT[1], ϖ_Cabannes=FT[1],
    bandSpecLim=UnitRange{Int}[], iBand=[1], F₀=copy(F0), SIF₀=zeros(FT,3,size(F0,2)))

function toa3(result, keep=Colon())
    a = Array(result)
    return FT.(a[1, 1:3, keep])
end

function build_o2(state, band_ν, solve_ν)
    # O2, p/T/q, and aerosols come from the same ABSCO-backed configuration as
    # the retrieval forward model. The historical Raman YAML is no longer an
    # independent source of atmospheric state or spectroscopy.
    p = load_parameters(; float_type=FT, architecture=:GPU,
                        nstreams=AEROSOL_NSTREAMS)
    set_common!(p, state)
    surface_anchor = RRSXCO2Common.surface_basis_grids(FT)[1]
    surface = band_surface(state.coeff[1], surface_anchor, solve_ν)
    select_band!(p, 1, solve_ν, surface)
    # Keep the direct solar direction outside the diffuse angular operator.
    # Forward RRS uses the dedicated Z₀/R₀/T₀ source columns, so embedding
    # μ₀ here only inflates the dense Raman workspace without changing the
    # radiative-transfer result.
    return model_from_parameters(p; external_solar=true,
                                 aerosol_anchor_bands=[band_ν])
end

function build_co2(state, original_band, band_ν, solve_ν)
    p = load_parameters(; float_type=FT, architecture=:GPU,
                        nstreams=AEROSOL_NSTREAMS)
    set_common!(p, state)
    surface_anchor = RRSXCO2Common.surface_basis_grids(FT)[original_band]
    surface = band_surface(state.coeff[original_band], surface_anchor, solve_ν)
    select_band!(p, original_band, solve_ν, surface)
    return model_from_parameters(p; external_solar=true,
                                 aerosol_anchor_bands=[band_ν])
end

function simulate_o2(state, grids, solar_T; rayleigh_core_only::Bool=true)
    bandν, solveν, keep = grids
    model = build_o2(state, bandν, solveν)
    F0, sources = source_set(solveν, state.sif_total > 0, solar_T)
    rrs_result = CoreRT.rt_run_toa(make_rrs(model,F0), model;
                                   i_band=1, sources)
    cabannes = toa3(rrs_result.elastic, keep)
    rrs = toa3(rrs_result.inelastic, keep)

    if rayleigh_core_only
        # Raman shoulders are source wavelengths for the inelastic solve.
        # Pure Rayleigh/noRS is monochromatic, so solving those discarded
        # samples repeats the full Fourier/layer calculation without changing
        # any retained radiance. Build the noRS model directly on the exact
        # retained core instead. `bandν` remains the complete spectral anchor
        # for aerosol and surface interpolation in both models.
        coreν = solveν[keep]
        model = nothing
        rrs_result = nothing
        F0 = nothing
        sources = nothing
        GC.gc()
        CUDA.reclaim()
        ray_model = build_o2(state, bandν, coreν)
        _, ray_sources = source_set(coreν, state.sif_total > 0, solar_T)
        ray_result = CoreRT.rt_run_toa(ray_model; i_band=1, sources=ray_sources)
        rayleigh = toa3(ray_result)
    else
        ray_result = CoreRT.rt_run_toa(model; i_band=1, sources)
        rayleigh = toa3(ray_result, keep)
    end
    return (; rayleigh, cabannes, rrs)
end

function simulate_co2(state, iband, bandν, solveν, solar_T)
    model = build_co2(state, iband, bandν, solveν)
    F0, sources = source_set(solveν, false, solar_T)
    result = CoreRT.rt_run_toa(model; i_band=1, sources)
    return toa3(result)
end

# Unchunked convenience method: the solve grid is the complete band.
simulate_co2(state, iband, ν, solar_T) = simulate_co2(state, iband, ν, ν, solar_T)

function write_wavelengths(grids)
    path = joinpath(OUT, "sim_wavelength.nc")
    NCDataset(path,"c") do ds
        RRSXCO2Common.write_absco_provenance!(ds.attrib)
        RRSXCO2Common.write_fourier_convergence_provenance!(ds.attrib)
        ds.attrib["spectral_grid"] = "uniform 0.1 cm-1 basis grids; wavelength is derived as 1e7/wavenumber"
        ds.attrib["o2a_raman_shoulder_cm-1"] = SHOULDER_CM
        ds.attrib["sza_deg"] = SZA_DEG
        ds.attrib["vza_deg"] = VZA_DEG
        ds.attrib["relative_azimuth_deg"] = RELATIVE_AZIMUTH_DEG
        ds.attrib["strong_co2_short_shoulder_points"] = Int32(8)
        ds.attrib["strong_co2_convolution_support_sigma"] = 6.0
        for (name,ν) in zip(("o2a","weak_co2","strong_co2"),grids)
            defDim(ds,name,length(ν))
            v=defVar(ds,"$(name)_wavelength",Float64,(name,)); v.attrib["units"]="nm"; v[:]=1e7 ./ ν
            w=defVar(ds,"$(name)_wavenumber",Float64,(name,)); w.attrib["units"]="cm-1"; w[:]=ν
        end
    end
end

function write_scene(state, grids, o2, weak, strong)
    path = joinpath(OUT,@sprintf("hiressim_%03d.nc",state.index))
    NCDataset(path,"c") do ds
        RRSXCO2Common.write_absco_provenance!(ds.attrib)
        defDim(ds,"stokes",3)
        for (name,ν) in zip(("o2a","weak_co2","strong_co2"),grids); defDim(ds,name,length(ν)); end
        for (name,data) in (("radiance_rayleigh_o2a",o2.rayleigh),
                            ("radiance_cabannes_o2a",o2.cabannes),
                            ("radiance_rrs_o2a",o2.rrs))
            v=defVar(ds,name,Float32,("stokes","o2a")); v.attrib["units"]="mW m-2 sr-1 (cm-1)-1"; v[:,:]=data
        end
        for (name,band,data) in (("radiance_rayleigh_weak_co2","weak_co2",weak),
                                 ("radiance_rayleigh_strong_co2","strong_co2",strong))
            v=defVar(ds,name,Float32,("stokes",band)); v.attrib["units"]="mW m-2 sr-1 (cm-1)-1"; v[:,:]=data
        end
        ds.attrib["state_index"] = Int32(state.index)
        ds.attrib["surface"] = state.surface; ds.attrib["aerosol_case"] = state.aerosol_case
        ds.attrib["sif_case"] = state.sif_case; ds.attrib["xco2_ppm"] = state.xco2_ppm
        ds.attrib["psurf_hpa"] = 1000.0
        ds.attrib["sza_deg"] = SZA_DEG
        ds.attrib["vza_deg"] = VZA_DEG
        ds.attrib["relative_azimuth_deg"] = RELATIVE_AZIMUTH_DEG
        ds.attrib["strong_co2_short_shoulder_points"] = Int32(8)
        ds.attrib["strong_co2_convolution_support_sigma"] = 6.0
        ds.attrib["atmospheric_layers"] = Int32(NLAYERS)
        ds.attrib["aod550_sulfate"] = state.aod550[1]
        ds.attrib["aod550_organic_carbon"] = state.aod550[2]
        ds.attrib["aod550_stratospheric"] = state.aod550[3]
        ds.attrib["sif_total_mW_m-2_sr-1"] = state.sif_total
        ds.attrib["source_state_table"] = "true_states.dat"
        ds.attrib["created"] = string(now())
        # Written only after every spectral variable has been populated, so
        # downstream processing can reject an interrupted partial file.
        ds.attrib["simulation_complete"] = Int32(1)
        ds.attrib["completed"] = string(now())
    end
    return path
end

function completed_scene(path)
    isfile(path) || return false
    return NCDataset(path, "r") do dataset
        haskey(dataset.attrib, "simulation_complete") &&
            Int(dataset.attrib["simulation_complete"]) == 1
    end
end

function main()
    mkpath(OUT)
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(DEVICE)
    states = selected_states()
    output = output_grids()
    solveν, keep = o2_solve_grid(output[1])
    write_wavelengths(output)
    solar_T = solar_interpolator()
    o2cache, weakcache, strongcache = Dict(), Dict(), Dict()
    for state in states
        path = joinpath(OUT,@sprintf("hiressim_%03d.nc",state.index))
        if isfile(path) && !FORCE
            completed_scene(path) || error(
                "existing scene is not marked complete: $path; inspect it, " *
                "then rerun with FORCE=1 to replace it")
            @info "skip existing" path
            continue
        end
        okey=(state.surface_index,state.aerosol_index,state.sif_index)
        ckey=(state.surface_index,state.aerosol_index,state.xco2_index)
        o2 = get!(o2cache,okey) do
            @info "simulate O2A" state.index okey nsolve=length(solveν)
            simulate_o2(state,(output[1],solveν,keep),solar_T)
        end
        weak = get!(weakcache,ckey) do
            @info "simulate weak CO2" state.index ckey
            simulate_co2(state,2,output[2],output[2],solar_T)
        end
        strong = get!(strongcache,ckey) do
            @info "simulate strong CO2" state.index ckey
            simulate_co2(state,3,output[3],output[3],solar_T)
        end
        write_scene(state,output,o2,weak,strong)
        CUDA.synchronize(); GC.gc(); CUDA.reclaim()
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
