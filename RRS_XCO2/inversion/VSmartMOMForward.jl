module VSmartMOMForward

using CUDA
using Distributions
using LinearAlgebra: dot
using NCDatasets
using vSmartMOM
using vSmartMOM.CoreRT

include(joinpath(@__DIR__, "instrument", "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "..", "scripts", "common.jl"))
using .SyntheticOCO2
using .RRSXCO2Common
using ..OptimalEstimation: ForwardEvaluation

export OCOForwardEvaluator,
       prepare_retrieval_parameters,
       apply_retrieval_state!,
       set_fixed_upper_co2_ppm!,
       column_averaged_co2_ppm,
       evaluate_oco_forward,
       surface_coefficient_transform,
       tau_ref_per_aod760

const ACTIVE_STATE_COUNT = 30
const FIXED_UPPER_CO2_LAYERS = 1:4
const ACTIVE_CO2_LAYERS = 5:16
const CO2_RANGE = 2:13
const LOG_AOD_RANGE = 14:16
const LOG_HEIGHT_RANGE = 17:19
const SURFACE_RANGE = 20:28
const SIF_RANGE = 29:30
const DEFAULT_COEFFICIENT_PATH = joinpath(
    @__DIR__, "instrument", "representative_stokes_coefficients.nc")
const DEFAULT_COMPONENT_PATH = joinpath(
    @__DIR__, "..", "truth_map", "scene_components.dat")
"""
Map truth-state Legendre coefficients into a retrieval solve-grid basis.

Truth-map surface polynomials are normalized on each complete output band
before any supplemental convolution points are added. The retrieval strong
band includes eight such points in its solve grid. Because the historical
`LambertianSurfaceLegendre` basis normalizes over spectral index, extending
that grid would otherwise silently change the physical surface. This affine
three-coefficient transform preserves the truth-map polynomial exactly; the
same matrix is applied to the analytic surface Jacobian by the evaluator.
"""
function surface_coefficient_transform(params, iband::Integer)
    1 <= iband <= length(RRSXCO2Common.BAND_NAMES) ||
        throw(BoundsError(RRSXCO2Common.BAND_NAMES, iband))
    FT = params.float_type
    # This is deliberately the same canonical basis used by the truth-map
    # builder. Keeping the definition in common.jl prevents an ILS-support
    # change from silently changing the physical surface polynomial in only
    # one side of the retrieval.
    base = RRSXCO2Common.surface_basis_grids(FT)[iband]
    solve = params.spec_bands[iband]
    length(base) > 1 && length(solve) > 1 || error(
        "surface coefficient transform requires multi-point spectral grids")
    base_span = last(base) - first(base)
    solve_span = last(solve) - first(solve)
    α = solve_span / base_span
    β = FT(2) * (first(solve) - first(base)) / base_span + α - one(FT)
    return FT[
        1 β (3β^2 - 1 + α^2) / 2
        0 α 3α*β
        0 0 α^2
    ]
end

function _validate_retrieval_convolution_support(params; support_sigma=6)
    for (iband, spec) in enumerate(BAND_SPECS)
        extension = required_wavenumber_extensions(
            collect(params.spec_bands[iband]), spec;
            step_cm=0.1, support_sigma)
        isempty(extension.short) && isempty(extension.long) || error(
            "retrieval $(spec.name) grid lacks Gaussian convolution support: " *
            "short_nodes=$(length(extension.short)) " *
            "long_nodes=$(length(extension.long))")
    end
    return params
end

"""Parse the fixed-microphysics conversion tau_ref(550 nm)/AOD(760 nm)."""
function tau_ref_per_aod760(path::AbstractString=DEFAULT_COMPONENT_PATH)
    isfile(path) || throw(ArgumentError("missing aerosol component catalog: $path"))
    for line in eachline(path)
        fields = split(strip(line))
        isempty(fields) && continue
        first(fields) == "aod760_0p28" || continue
        length(fields) >= 9 || error("malformed aod760_0p28 row in $path")
        aod550 = parse.(Float64, fields[2:4])
        aod760 = parse.(Float64, fields[6:8])
        all(>(0), aod760) || error("AOD760 conversion values must be positive")
        return aod550 ./ aod760
    end
    error("aod760_0p28 conversion row was not found in $path")
end

function _truncate_profile!(params, psurf=1000.0)
    return RRSXCO2Common.truncate_profile_to_surface!(params, psurf)
end

"""
Construct the fixed 16-layer, three-band retrieval template.

The template is parsed directly at the requested float type and architecture;
changing those fields after parsing would leave nested aerosol and surface
objects at the wrong precision. The final 16-layer profile becomes the input
grid so layer-resolved CO2 state entries map one-to-one without interpolation.
"""
function prepare_retrieval_parameters(;
        architecture::Symbol=:GPU,
        float_type::DataType=Float32,
        nstreams::Int=9,
        config_path::AbstractString=RRSXCO2Common.CONFIG)
    architecture in (:CPU, :GPU) || throw(ArgumentError(
        "retrieval architecture must be :CPU or :GPU"))
    float_type in (Float32, Float64) || throw(ArgumentError(
        "retrieval float type must be Float32 or Float64"))
    nstreams == 9 || throw(ArgumentError(
        "the first retrieval suite is fixed to nstreams=9"))
    params = RRSXCO2Common.load_parameters(;
        float_type, architecture, nstreams, config_path)
    params.spec_bands = [collect(grid) for grid in
                         RRSXCO2Common.basis_grids(float_type)]
    FT = params.float_type
    params.sza = FT(30)
    params.vza = FT[0]
    params.vaz = FT[0]
    RRSXCO2Common.prepare_shared_profile!(params; psurf=1000.0, nlayers=16)
    params.absorption_params.vmr["CO2"] = fill(FT(400e-6), 16)
    return _validate_retrieval_convolution_support(params)
end

mutable struct OCOForwardEvaluator
    base_parameters
    coefficients::Dict{Symbol,Vector{Float64}}
    tau_ref_scale::Vector{Float64}
    fixed_upper_co2_vmr::Float64
    solar_transmission
    synchronize_backend::Function
end

function OCOForwardEvaluator(;
        architecture::Symbol=:GPU,
        float_type::DataType=Float32,
        nstreams::Int=9,
        fixed_upper_co2_ppm::Real=400.0,
        coefficient_path::AbstractString=DEFAULT_COEFFICIENT_PATH)
    isfinite(fixed_upper_co2_ppm) && 0 < fixed_upper_co2_ppm < 10_000 ||
        throw(ArgumentError("fixed upper-atmosphere CO2 must lie in (0,10000) ppm"))
    params = prepare_retrieval_parameters(;
        architecture, float_type, nstreams)
    synchronize_backend = () -> nothing
    if architecture == :GPU
        CUDA.functional() || error("CUDA is not functional on this host")
        device = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
        CUDA.device!(device)
        synchronize_backend = CUDA.synchronize
    end
    coefficients = read_representative_coefficients(coefficient_path)
    solar_transmission = RRSXCO2Common.solar_interpolator(float_type)
    return OCOForwardEvaluator(
        params, coefficients, tau_ref_per_aod760(),
        Float64(fixed_upper_co2_ppm) * 1e-6, solar_transmission,
        synchronize_backend)
end

"""
    set_fixed_upper_co2_ppm!(evaluator, truth_xco2_ppm)

Set the non-retrieved CO2 VMR in layers 1:4 to the uniform truth-scene value.
The synthetic truth map varies its full 16-layer profile uniformly between
380 and 440 ppm, while these four layers have zero prior variance and are
absent from the 30-element active state. Configuring them per scene therefore
restores representability without adding retrieval coordinates.
"""
function set_fixed_upper_co2_ppm!(evaluator::OCOForwardEvaluator,
                                  truth_xco2_ppm::Real)
    isfinite(truth_xco2_ppm) && 0 < truth_xco2_ppm < 10_000 ||
        throw(ArgumentError("fixed upper-atmosphere CO2 must lie in (0,10000) ppm"))
    evaluator.fixed_upper_co2_vmr = Float64(truth_xco2_ppm) * 1e-6
    return evaluator
end

"""
    column_averaged_co2_ppm(evaluator, state)

Return the dry-air-column-weighted CO2 volume-mixing ratio in ppm for one
active retrieval state,

```math
X_{CO_2} = 10^6
    \\frac{\\sum_z q_{CO_2,z} N_{dry,z}}{\\sum_z N_{dry,z}}.
```

The diagnostic uses the same pressure interfaces, specific-humidity profile,
and hydrostatic dry-air columns as the retrieval forward model. Consequently,
the layer weights respond to retrieved surface pressure. The non-active upper
CO2 layers use `evaluator.fixed_upper_co2_vmr`, which is set from the current
truth scene before each retrieval.
"""
function column_averaged_co2_ppm(evaluator::OCOForwardEvaluator,
                                  state::AbstractVector)
    length(state) == ACTIVE_STATE_COUNT || throw(DimensionMismatch(
        "OCO retrieval state must contain $ACTIVE_STATE_COUNT entries"))
    params = evaluator.base_parameters
    FT = params.float_type

    p_half = copy(params.p)
    p_half[end] = FT(state[1])
    p_half[end] > p_half[end - 1] || throw(ArgumentError(
        "surface pressure $(state[1]) hPa is not below the final layer top " *
        "$(p_half[end - 1]) hPa"))

    co2 = fill(FT(evaluator.fixed_upper_co2_vmr), 16)
    co2[ACTIVE_CO2_LAYERS] .= FT.(state[CO2_RANGE])
    all(x -> isfinite(x) && 0 < x < 0.01, co2) || throw(ArgumentError(
        "CO2 VMR trial state is outside the physical guard range (0,0.01)"))

    obs_alt = params.obs_alt isa Real ?
        FT(params.obs_alt) : convert(Vector{FT}, params.obs_alt)
    profile, _ = CoreRT.prepare_observer_profile(
        convert(Vector{FT}, params.T), convert(Vector{FT}, p_half),
        convert(Vector{FT}, params.q), Dict("CO2" => co2), obs_alt, -1)
    weights = Float64.(profile.vcd_dry)
    vmr = Float64.(profile.vmr["CO2"])
    return 1.0e6 * dot(vmr, weights) / sum(weights)
end

"""Map the 30 retrieval coordinates into a fresh parameter template."""
function apply_retrieval_state!(params, state::AbstractVector,
                                tau_ref_scale::AbstractVector;
                                fixed_upper_co2_vmr::Real=400e-6)
    length(state) == ACTIVE_STATE_COUNT || throw(DimensionMismatch(
        "OCO retrieval state must contain $ACTIVE_STATE_COUNT entries"))
    length(tau_ref_scale) == 3 || throw(DimensionMismatch(
        "three AOD760-to-tau_ref scale factors are required"))
    FT = params.float_type
    psurf = FT(state[1])
    psurf > params.p[end-1] || throw(ArgumentError(
        "surface pressure $psurf hPa is not below the final layer top " *
        "$(params.p[end-1]) hPa"))
    params.p[end] = psurf

    isfinite(fixed_upper_co2_vmr) && 0 < fixed_upper_co2_vmr < 0.01 ||
        throw(ArgumentError(
            "fixed upper-atmosphere CO2 VMR must lie in (0,0.01)"))
    co2 = fill(FT(fixed_upper_co2_vmr), 16)
    co2[ACTIVE_CO2_LAYERS] .= FT.(state[CO2_RANGE])
    all(x -> 0 < x < 0.01, co2) || throw(ArgumentError(
        "CO2 VMR trial state is outside the physical guard range (0,0.01)"))
    params.absorption_params.vmr["CO2"] = co2

    aod760 = exp.(Float64.(state[LOG_AOD_RANGE]))
    heights = exp.(Float64.(state[LOG_HEIGHT_RANGE]))
    all(isfinite, aod760) && all(isfinite, heights) || throw(ArgumentError(
        "log-aerosol coordinate overflowed"))
    for iaerosol in 1:3
        aerosol = params.scattering_params.rt_aerosols[iaerosol]
        aerosol.τ_ref = FT(aod760[iaerosol] * tau_ref_scale[iaerosol])
        sigma0 = aerosol.profile.σ
        aerosol.profile = LogNormal(log(FT(heights[iaerosol])), sigma0)
    end

    for iband in 1:3
        first_index = first(SURFACE_RANGE) + 3 * (iband - 1)
        coefficients = FT.(state[first_index:first_index + 2])
        native_coefficients =
            surface_coefficient_transform(params, iband) * coefficients
        params.brdf[iband] =
            CoreRT.LambertianSurfaceLegendre(native_coefficients)
    end
    return (; SIF760=FT(state[first(SIF_RANGE)]),
            mSIF=FT(state[last(SIF_RANGE)]),
            tau_ref=FT.(aod760 .* tau_ref_scale),
            aerosol_height=FT.(heights))
end

function _canonical_toa(result)
    toa = Array(result.toa)
    jacobian_local = Array(result.toa_jacobian)
    ndims(toa) == 3 || error("unexpected TOA radiance shape $(size(toa))")
    ndims(jacobian_local) == 4 || error(
        "unexpected TOA Jacobian shape $(size(jacobian_local))")
    # ObserverRTResult uses (VZA, Stokes, spectral[, parameter]); the
    # synthetic OCO operator uses (Stokes, spectral[, parameter]).
    size(toa, 1) == 1 || error(
        "retrieval evaluator requires exactly one VZA; received $(size(toa, 1))")
    size(jacobian_local, 1) == 1 || error(
        "retrieval Jacobian requires exactly one VZA; received " *
        "$(size(jacobian_local, 1))")
    return Array(@view(toa[1, :, :])), jacobian_local
end

function _assert_finite_band_output(stokes, jacobian, layout, band_name)
    all(isfinite, stokes) || error(
        "$band_name forward Stokes spectrum contains " *
        "$(count(!isfinite, stokes)) non-finite entries")
    bad_columns = String[]
    names = layout.global_parameter_names
    for iparam in axes(jacobian, 3)
        column = @view jacobian[:, :, iparam]
        nbad = count(!isfinite, column)
        iszero(nbad) || push!(bad_columns, "$(names[iparam])=$nbad")
    end
    isempty(bad_columns) || error(
        "$band_name high-resolution Jacobian contains non-finite entries: " *
        join(bad_columns, ", "))
    return nothing
end

"""Evaluate noRS vSmartMOM and its selected Jacobian through the OCO operator."""
function evaluate_oco_forward(evaluator::OCOForwardEvaluator,
                              state::AbstractVector)
    params = deepcopy(evaluator.base_parameters)
    physical = apply_retrieval_state!(
        params, state, evaluator.tau_ref_scale;
        fixed_upper_co2_vmr=evaluator.fixed_upper_co2_vmr)
    evaluator.synchronize_backend()
    model_pair = nothing
    model_build_seconds = @elapsed begin
        model_pair = model_from_parameters(
            OCO_RRS_synth(), params; external_solar=true)
        evaluator.synchronize_backend()
    end
    model, lin_model = model_pair

    band_measurements = Vector{Vector{Float64}}(undef, 3)
    band_jacobians = Vector{Matrix{Float64}}(undef, 3)
    rt_seconds = 0.0
    instrument_seconds = 0.0
    for (iband, spec) in enumerate(BAND_SPECS)
        sources = RRSXCO2Common.sources_for_band(
            params, iband; SIF760=physical.SIF760, mSIF=physical.mSIF,
            solar_T=evaluator.solar_transmission)
        result = nothing
        rt_seconds += @elapsed begin
            result = rt_run_lin(model, lin_model; i_band=iband, sources)
            evaluator.synchronize_backend()
        end
        stokes, local_jacobian = _canonical_toa(result)
        global_jacobian = Array(globalize_jacobian(
            local_jacobian, result.layout))
        size(global_jacobian, 1) == 1 || error(
            "retrieval evaluator requires exactly one VZA; received " *
            "$(size(global_jacobian, 1))")
        jacobian = Array(@view(global_jacobian[1, :, :, :]))
        _assert_finite_band_output(stokes, jacobian, result.layout, spec.name)
        wavelength_nm = 1.0e7 ./ Float64.(params.spec_bands[iband])
        instrument_seconds += @elapsed begin
            band_measurements[iband] = process_stokes_spectrum(
                wavelength_nm, stokes, evaluator.coefficients[spec.name], spec)
            band_jacobians[iband] = process_stokes_jacobian(
                wavelength_nm, jacobian, evaluator.coefficients[spec.name], spec)
            first_surface = first(SURFACE_RANGE) + 3 * (iband - 1)
            surface_columns = first_surface:first_surface + 2
            band_jacobians[iband][:, surface_columns] =
                band_jacobians[iband][:, surface_columns] *
                Float64.(surface_coefficient_transform(params, iband))
        end
    end

    measurement = vcat(band_measurements...)
    jacobian = vcat(band_jacobians...)
    # The selected vSmartMOM basis is native (tau_ref,z0); the OE basis uses
    # log(AOD760),log(z0). At fixed microphysics dF/dlog(AOD760)=tau_ref*dF/dtau_ref.
    jacobian[:, LOG_AOD_RANGE] .*= reshape(Float64.(physical.tau_ref), 1, :)
    jacobian[:, LOG_HEIGHT_RANGE] .*= reshape(
        Float64.(physical.aerosol_height), 1, :)
    lengths = length.(band_measurements)
    stops = cumsum(lengths)
    starts = vcat(1, stops[1:end-1] .+ 1)
    ranges = UnitRange{Int}[first:last for (first, last) in zip(starts, stops)]
    timing = (model_build_seconds=model_build_seconds,
              rt_linearized_seconds=rt_seconds,
              instrument_seconds=instrument_seconds,
              forward_seconds=model_build_seconds + rt_seconds,
              jacobian_seconds=rt_seconds)
    return ForwardEvaluation(measurement, jacobian, ranges; timing)
end

(evaluator::OCOForwardEvaluator)(state) = evaluate_oco_forward(evaluator, state)

end # module VSmartMOMForward
