module RRSXCO2Common

using DataInterpolations
using YAML
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel: planck_spectrum_wn, solar_transmission_from_file

export WORKFLOW_ROOT, ROOT, CONFIG, BAND_NAMES, ABSCO_VERSION, configure_absco_luts!,
       configure_luts!, absco_lut_paths, load_parameters, wavelengths_nm,
       truncate_profile_to_surface!, materialize_profile!, prepare_shared_profile!,
       surface_basis_grids, basis_grids, raman_solve_grid,
       solar_interpolator, sources_for_band,
       write_absco_provenance!, write_fourier_convergence_provenance!

# Source lives under sandbox/workflows on the surface-split branch, whereas
# large generated products can remain in the historical top-level RRS_XCO2
# tree.  Keep those locations independent so changing branches never implies
# moving or regenerating truth/retrieval data.
const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, ".."))
const ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const CONFIG = normpath(get(
    ENV, "RRS_XCO2_CONFIG",
    joinpath(WORKFLOW_ROOT, "config", "oco_grass_3aerosol.yaml")))
const BAND_NAMES = ("o2a", "weak_co2", "strong_co2")
const ABSCO_VERSION = "5.2"
const DEFAULT_SURFACE_PRESSURE_HPA = 1000.0
const DEFAULT_PROFILE_LAYERS = 16
const SOLAR_OUT = get(ENV, "SOLAR_OUT",
    joinpath(homedir(), "Raman_misc", "workflows", "worktrees",
             "uni_vSmartMOM_sanghavi_2025-03-18", "src", "SolarModel", "solar.out"))
const _SOLAR_INTERPOLATORS = Dict{DataType,Any}()

"Return the spectroscopy tables used by the three-band experiment."
function absco_lut_paths()
    lut_dir = get(ENV, "VSMARTMOM_ABSCO_DIR",
        "/net/fluo/data1/ABSCO_CS_Database/v5.2_final")
    hitran_lut_dir = get(ENV, "VSMARTMOM_HITRAN_LUT_DIR",
        joinpath(homedir(), "data", "HITRAN_LUTs"))
    defaults = (
        O2_ABSCO_LUT = "o2_v52_v2.jld2",
        O2_H2O_LUT = joinpath(hitran_lut_dir, "H2O.jld2"),
        WEAK_H2O_ABSCO_LUT = "wh2o_v52.jld2",
        WEAK_CO2_ABSCO_LUT = "wco2_v52.jld2",
        STRONG_H2O_ABSCO_LUT = "sh2o_v52.jld2",
        STRONG_CO2_ABSCO_LUT = "sco2_v52.jld2",
    )
    return NamedTuple{keys(defaults)}(Tuple(
        get(ENV, String(key), key == :O2_H2O_LUT ? filename : joinpath(lut_dir, filename))
        for (key, filename) in pairs(defaults)))
end

"Configure and validate the spectroscopy tables used by truth and retrieval."
function configure_absco_luts!()
    paths = absco_lut_paths()
    for (key, path) in pairs(paths)
        haskey(ENV, String(key)) || (ENV[String(key)] = path)
        isfile(path) || error(
            "Missing $key at $path. Set $key or VSMARTMOM_ABSCO_DIR.")
    end
    return paths
end

# Compatibility name used by existing RRS_XCO2 scripts. It now means ABSCO
# only; silently selecting the historical HITRAN-built LUTs would invalidate
# truth/retrieval closure.
configure_luts!() = configure_absco_luts!()

"Write the experiment's spectroscopy/profile provenance into NetCDF attributes."
function write_absco_provenance!(attributes)
    paths = absco_lut_paths()
    attributes["spectroscopy_database"] = "ABSCO"
    attributes["spectroscopy_version"] = ABSCO_VERSION
    attributes["o2_absco_lut"] = abspath(paths.O2_ABSCO_LUT)
    attributes["o2_h2o_lut"] = abspath(paths.O2_H2O_LUT)
    attributes["weak_h2o_absco_lut"] = abspath(paths.WEAK_H2O_ABSCO_LUT)
    attributes["weak_co2_absco_lut"] = abspath(paths.WEAK_CO2_ABSCO_LUT)
    attributes["strong_h2o_absco_lut"] = abspath(paths.STRONG_H2O_ABSCO_LUT)
    attributes["strong_co2_absco_lut"] = abspath(paths.STRONG_CO2_ABSCO_LUT)
    attributes["h2o_line_absorption_by_band"] =
        "o2a=rebuilt HITRAN LUT (matches archived truth); weak_co2=ABSCO; strong_co2=ABSCO"
    attributes["o2_vmr"] = 0.21
    attributes["atmospheric_profile_configuration"] = abspath(CONFIG)
    attributes["profile_preparation"] =
        "truncate native p/T/q at 1000 hPa, then materialize one common 16-layer profile"
    return attributes
end

"Write the Fourier-loop policy selected by the shared experiment YAML."
function write_fourier_convergence_provenance!(attributes)
    input = YAML.load_file(CONFIG)
    numerics = get(input["radiative_transfer"], "numerics", Dict{Any,Any}())
    selection = lowercase(string(get(numerics, "fourier_convergence", "all")))
    attributes["fourier_convergence"] = selection
    if selection in ("intensity", "stokes", "iqu")
        attributes["fourier_convergence_tolerance"] =
            Float64(get(numerics, "fourier_tolerance", 1e-5))
        attributes["fourier_convergence_minimum_m"] =
            Int32(get(numerics, "fourier_min_m", 3))
        attributes["fourier_convergence_guard_rule"] =
            "test only for m > min(fourier_min_m-1, solver m_max)"
        attributes["fourier_convergence_consecutive_moments"] =
            Int32(get(numerics, "fourier_n_consecutive", 2))
        attributes["fourier_convergence_observable"] = selection == "intensity" ?
            "forward Stokes I at every output wavelength and view" :
            "every forward Stokes component at every output wavelength and view"
    end
    return attributes
end

"Load the shared three-band parameters at a requested precision/backend."
function load_parameters(; float_type::DataType=Float64,
                           architecture::Symbol=:GPU,
                           nstreams::Union{Nothing,Integer}=nothing,
                           config_path::AbstractString=CONFIG)
    architecture in (:CPU, :GPU) || throw(ArgumentError(
        "architecture must be :CPU or :GPU"))
    float_type in (Float32, Float64) || throw(ArgumentError(
        "float_type must be Float32 or Float64"))
    configure_absco_luts!()
    input = YAML.load_file(config_path)
    rt = input["radiative_transfer"]
    rt["float_type"] = string(float_type)
    rt["architecture"] = architecture == :GPU ? "GPU()" : "CPU()"
    nstreams === nothing || (rt["nstreams"] = Int(nstreams))
    params = parameters_from_dict(input)
    length(params.scattering_params.rt_aerosols) == 3 ||
        error("Expected exactly three aerosol species")
    length(params.spec_bands) == 3 || error("Expected exactly three OCO bands")
    return params
end

"Canonical coefficient-definition grids before supplemental ILS support nodes."
function surface_basis_grids(::Type{FT}=Float32) where {FT<:AbstractFloat}
    step = FT(0.1)
    o2a = collect(FT(1e7 / 773):step:FT(1e7 / 757))
    weak = collect(FT(1e7 / 1622):step:FT(1e7 / 1589))
    strong = collect(FT(1e7 / 2084):step:FT(1e7 / 2042))
    return (o2a, weak, strong)
end

"Canonical 0.1-cm⁻¹ solve grids, including strong-band ILS support."
function basis_grids(::Type{FT}=Float32) where {FT<:AbstractFloat}
    o2a, weak, strong_base = surface_basis_grids(FT)
    step = FT(0.1)
    strong_shoulder = FT.(Float64(last(strong_base)) .+
                           Float64(step) .* (1:8))
    return (o2a, weak, vcat(strong_base, strong_shoulder))
end

"""
    raman_solve_grid(core; shoulder_cm=234, step_cm=0.1)

Surround a canonical retained O₂ A-band grid with Raman source shoulders.
Only the shoulders are constructed: `core` itself is inserted verbatim, and
the returned `keep` range therefore satisfies `solve[keep] == core`
bit-for-bit. noRS calculations should use `core` directly and must not pay for
these solve-only Raman wavelengths.

This deliberately avoids rebuilding the retained samples as part of a longer
floating-point range. A different range origin can move internal Float32
nodes by one ULP even when both ranges nominally use the same 0.1 cm⁻¹ step.
"""
function raman_solve_grid(core::AbstractVector{FT};
                          shoulder_cm::Real=234,
                          step_cm::Real=0.1) where {FT<:AbstractFloat}
    isempty(core) && throw(ArgumentError("O2 core grid must be nonempty"))
    shoulder = FT(shoulder_cm)
    step = FT(step_cm)
    shoulder >= zero(FT) || throw(ArgumentError(
        "Raman shoulder width must be nonnegative"))
    step > zero(FT) || throw(ArgumentError("spectral step must be positive"))
    nshoulder = round(Int, shoulder / step)
    isapprox(FT(nshoulder) * step, shoulder;
             atol=2eps(FT) * max(abs(shoulder), one(FT)), rtol=0) ||
        throw(ArgumentError(
            "Raman shoulder width must be an integer multiple of the spectral step"))

    offsets = step .* FT.(1:nshoulder)
    left = reverse(FT(first(core)) .- offsets)
    right = FT(last(core)) .+ offsets
    solve = vcat(left, FT.(core), right)
    keep = (nshoulder + 1):(nshoulder + length(core))
    solve[keep] == FT.(core) || error(
        "internal error: Raman solve grid displaced its retained O2 core")
    return solve, keep
end

"Truncate the native p/T/q grid at a requested surface pressure."
function truncate_profile_to_surface!(params,
                                      psurf::Real=DEFAULT_SURFACE_PRESSURE_HPA)
    p0, T0, q0 = Float64.(params.p), Float64.(params.T), Float64.(params.q)
    first(p0) < psurf <= last(p0) || throw(ArgumentError(
        "surface pressure $psurf hPa lies outside the profile"))
    k = searchsortedlast(p0, psurf - eps(Float64(psurf)))
    pcenter = (p0[1:end-1] .+ p0[2:end]) ./ 2
    function at_pressure(values, pressure)
        hi = searchsortedfirst(pcenter, pressure)
        hi <= 1 && return values[1]
        hi > length(values) && return values[end]
        lo = hi - 1
        weight = (log(pressure) - log(pcenter[lo])) /
                 (log(pcenter[hi]) - log(pcenter[lo]))
        return values[lo] + weight * (values[hi] - values[lo])
    end
    FT = params.float_type
    params.p = FT.(vcat(p0[1:k], psurf))
    params.T = FT.(vcat(T0[1:k-1], at_pressure(T0, psurf)))
    params.q = FT.(vcat(q0[1:k-1], at_pressure(q0, psurf)))
    return params
end

"Materialize profile reduction once so every model consumes identical p/T/q arrays."
function materialize_profile!(params, nlayers::Integer=DEFAULT_PROFILE_LAYERS)
    nlayers > 0 || throw(ArgumentError("nlayers must be positive"))
    FT = params.float_type
    vmr = params.absorption_params === nothing ? Dict() :
          params.absorption_params.vmr
    obs_alt = params.obs_alt isa Real ? FT(params.obs_alt) :
              Vector{FT}(params.obs_alt)
    profile, _ = CoreRT.prepare_observer_profile(
        Vector{FT}(params.T), Vector{FT}(params.p), Vector{FT}(params.q),
        vmr, obs_alt, nlayers)
    length(profile.T) == nlayers || error(
        "profile reduction produced $(length(profile.T)) layers, expected $nlayers")
    params.T = copy(profile.T)
    params.p = copy(profile.p_half)
    params.q = copy(profile.q)
    params.profile_reduction_n = -1
    return params
end

"Apply the shared surface-pressure truncation and vertical reduction exactly once."
function prepare_shared_profile!(params;
                                 psurf::Real=DEFAULT_SURFACE_PRESSURE_HPA,
                                 nlayers::Integer=DEFAULT_PROFILE_LAYERS)
    truncate_profile_to_surface!(params, psurf)
    return materialize_profile!(params, nlayers)
end

wavelengths_nm(ν) = 1.0e7 ./ collect(ν)

"""High-resolution solar transmission used by the truth-map simulations."""
function solar_interpolator(::Type{FT}) where {FT<:AbstractFloat}
    return get!(_SOLAR_INTERPOLATORS, FT) do
        isfile(SOLAR_OUT) || error(
            "Required high-resolution solar spectrum missing: $SOLAR_OUT")
        solar = solar_transmission_from_file(SOLAR_OUT)
        # Keep this construction identical to generate_truth_map.jl. The first
        # three rows are file metadata, and DataInterpolations takes dependent
        # values before coordinates.
        LinearInterpolation(FT.(solar[4:end, 2]), FT.(solar[4:end, 1]);
                            extrapolation=ExtrapolationType.Linear)
    end
end

"""
Physical solar beam for every band; retrievable SIF is present only in O2 A.

Pass `solar_T=solar_interpolator(FT)` when Fraunhofer structure is required.
The default remains the historical Planck-only source used by standalone
benchmarks; the retrieval adapter opts into the truth-map spectrum explicitly.
"""
function sources_for_band(params, iband; SIF760=nothing, mSIF=nothing,
                          solar_T=nothing)
    FT = params.float_type
    ν = FT.(params.spec_bands[iband])
    nPol, nSpec = params.polarization_type.n, length(ν)
    planck = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .*
                 FT(2.1629e-5 * π))
    solar = solar_T === nothing ? planck : FT.(solar_T.(ν)) .* planck
    F₀ = zeros(FT, nPol, nSpec)
    @views F₀[1, :] .= solar
    beam = CoreRT.SolarBeam(F₀=F₀)
    iband == 1 || return beam
    state = sif_reference_state(total_sif=0.5, reference_wavelength_nm=760)
    SIF760 === nothing && (SIF760 = state.SIF760)
    mSIF === nothing && (mSIF = state.mSIF)
    return beam + CoreRT.SurfaceSIF(SIF760=FT(SIF760), mSIF=FT(mSIF),
        wavenumber_cm1=ν)
end

end
