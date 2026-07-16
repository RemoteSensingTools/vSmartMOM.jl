using Test
using vSmartMOM
using vSmartMOM.CoreRT
using YAML

const _MS_CONFIG = joinpath(@__DIR__, "..", "config", "quickstart.yaml")

function _height_model(obs_alt; reduction=-1)
    params = parameters_from_yaml(_MS_CONFIG)
    params.obs_alt = obs_alt
    params.profile_reduction_n = reduction
    return model_from_parameters(params)
end

function _direct_height_model(::Type{FT}, obs_alt;
                              vza=[60.0], vaz=fill(180.0, length(vza)),
                              polarization=vSmartMOM.Scattering.Stokes_I(),
                              architecture=vSmartMOM.Architectures.CPU()) where {FT}
    params = parameters_from_yaml(_MS_CONFIG)
    params.float_type = FT
    params.architecture = architecture
    params.obs_alt = obs_alt
    params.vza = Float64.(vza)
    params.vaz = Float64.(vaz)
    params.polarization_type = polarization
    return model_from_parameters(params)
end

"Direct-beam attenuation above one quickstart observer interface."
function _expected_unscattered(model, boundary_index, F₀)
    FT = CoreRT.float_type(model)
    τ_above = sum(@view model.τ_abs[1][1, 1:boundary_index]) +
              sum(@view model.τ_rayl[1][1, 1:boundary_index])
    return FT(0.5 / π) * F₀ * exp(-τ_above / model.quad_points.μ₀)
end

@testset "observer-height profile framing" begin
    model = _height_model([0.0, 5.0])
    geom = model.obs_geom
    z_half = CoreRT.half_level_altitudes(model.profile)

    @test geom.include_toa
    @test geom.include_boa
    @test geom.sensor_altitudes == [5.0]
    @test length(geom.sensor_levels) == 1
    @test z_half[only(geom.sensor_levels) + 1] ≈ 5.0 atol=1e-7
    @test length(model.profile.T) == 3  # original two layers split at H
    @test length(model.profile.p_half) == 4
    @test all(diff(model.profile.p_half) .> 0)
    @test all(model.profile.Δz .> 0)
    @test size(model.τ_rayl[1], 2) == length(model.profile.T)
    @test size(model.τ_abs[1], 2) == length(model.profile.T)

    two_levels = _height_model([3.0, 8.0])
    @test two_levels.obs_geom.sensor_altitudes == [8.0, 3.0]
    @test CoreRT.half_level_altitudes(two_levels.profile)[
        two_levels.obs_geom.sensor_levels .+ 1] ≈ [8.0, 3.0] atol=1e-7

    mixed = _height_model([5.0, 5.0, 1.0e6])
    @test mixed.obs_geom.include_toa
    @test !mixed.obs_geom.include_boa
    @test mixed.obs_geom.sensor_altitudes == [5.0]  # duplicate collapsed

    # A request already present on the input grid must reuse that interface,
    # not create a near-duplicate through exp(log(p)) roundoff.
    base_model = _height_model([0.0])
    base_z_half = CoreRT.half_level_altitudes(base_model.profile)
    existing_height = base_z_half[2]
    existing = _height_model(existing_height)
    @test length(existing.profile.T) == length(base_model.profile.T)
    @test existing.obs_geom.sensor_altitudes == [existing_height]
    @test existing.obs_geom.sensor_levels == [1]

    # Machine-equivalent values at the top interface resolve to TOA rather
    # than constructing a numerically zero-thickness top layer.
    near_toa = _height_model(prevfloat(first(base_z_half)))
    @test near_toa.obs_geom.include_toa
    @test !near_toa.obs_geom.include_boa
    @test isempty(near_toa.obs_geom.sensor_levels)

    # Legacy profiles may define q/VMR on the N+1 pressure interfaces. They
    # are normalized to the N layer centers in log pressure before H framing.
    legacy_profile, _ = CoreRT.prepare_observer_profile(
        Float64[250, 275], Float64[100, 500, 1000],
        Float64[0.0, 0.01, 0.02], Dict("X" => [1.0, 2.0, 4.0]),
        Float64[0], -1)
    @test length(legacy_profile.q) == 2
    @test length(legacy_profile.vmr["X"]) == 2
    @test legacy_profile.q == CoreRT._layer_centered_input(
        "specific humidity", Float64[0.0, 0.01, 0.02],
        Float64[100, 500, 1000])
    @test legacy_profile.vmr["X"] == CoreRT._layer_centered_input(
        "VMR X", Float64[1.0, 2.0, 4.0], Float64[100, 500, 1000])

    interface_params = parameters_from_yaml(_MS_CONFIG)
    interface_params.obs_alt = 5.0
    interface_params.q = [0.0, 0.01, 0.02]
    interface_model = model_from_parameters(interface_params)
    @test length(interface_model.profile.q) == 3
    @test CoreRT.half_level_altitudes(interface_model.profile)[
        only(interface_model.obs_geom.sensor_levels) + 1] ≈ 5.0 atol=1e-7

    # Exercise an N+1 trace-gas profile through parsing, model construction,
    # H insertion, and a BatchContext update without requiring HITRAN data.
    vmr_config = YAML.load_file(_MS_CONFIG)
    vmr_config["absorption"] = Dict(
        "fixed_molecules" => [String[]],
        "variable_molecules" => [String[]],
        "vmr" => Dict("X" => [1.0, 2.0, 4.0]),
        "broadening" => "Voigt()",
        "CEF" => "HumlicekWeidemann32SDErrorFunction()",
        "wing_cutoff" => 25,
    )
    vmr_params = parameters_from_dict(vmr_config)
    vmr_params.obs_alt = 5.0
    vmr_model = model_from_parameters(vmr_params)
    input_p_full = (vmr_params.p[1:end-1] .+ vmr_params.p[2:end]) ./ 2
    X_layer = CoreRT._layer_centered_input(
        "VMR X", [1.0, 2.0, 4.0], vmr_params.p)
    @test vmr_model.profile.vmr["X"] ≈ CoreRT._interp_layer_field(
        input_p_full, X_layer, vmr_model.profile.p_full)

    vmr_ctx = BatchContext(vmr_params)
    q_interfaces = [0.0, 0.01, 0.02]
    X_interfaces = [1.5, 2.5, 4.5]
    update_model!(vmr_ctx;
                  q=q_interfaces, vmr=Dict("X" => X_interfaces))
    X_updated = CoreRT._layer_centered_input(
        "VMR X", X_interfaces, vmr_params.p)
    @test vmr_ctx.model.profile.vmr["X"] ≈ CoreRT._interp_layer_field(
        input_p_full, X_updated, vmr_ctx.model.profile.p_full)
    @test vmr_ctx.current_q == CoreRT._layer_centered_input(
        "specific humidity", q_interfaces, vmr_params.p)
    @test CoreRT.half_level_altitudes(vmr_ctx.model.profile)[
        only(vmr_ctx.model.obs_geom.sensor_levels) + 1] ≈ 5.0 atol=1e-7
end

@testset "observer-height CIA VMR reframing" begin
    mktempdir() do tmpdir
        # One constant-temperature X-X CIA block spanning the quickstart's
        # single spectral point. Keep the fixed-width HITRAN header explicit
        # so this regression is self-contained and does not need CIA data.
        cia_path = joinpath(tmpdir, "X-X.cia")
        header = rpad("X-X", 20) *
                 lpad("12986.000", 10) *
                 lpad("12988.000", 10) *
                 lpad("2", 7) *
                 lpad("250.0", 7) *
                 lpad("1.000e-43", 10)
        open(cia_path, "w") do io
            println(io, header)
            println(io, "12986.000 1.000e-43")
            println(io, "12988.000 1.000e-43")
        end

        config = YAML.load_file(_MS_CONFIG)
        config["absorption"] = Dict(
            "fixed_molecules" => [String[]],
            "variable_molecules" => [String[]],
            # Deliberately use the original N-layer shape. Inserting H creates
            # N+1 layers, so a stale use of the input dict would bounds-error.
            "vmr" => Dict("X" => [0.1, 0.2]),
            "broadening" => "Voigt()",
            "CEF" => "HumlicekWeidemann32SDErrorFunction()",
            "wing_cutoff" => 25,
            "cia_files" => [cia_path],
        )
        params = parameters_from_dict(config)
        params.obs_alt = 5.0
        @test length(params.absorption_params.vmr["X"]) == 2

        model = model_from_parameters(params)
        @test length(model.profile.T) == 3
        @test length(model.profile.vmr["X"]) == 3

        table = vSmartMOM.Absorption.load_cia_table(
            cia_path, params.spec_bands[1])
        expected = zeros(Float64, length(params.spec_bands[1]),
                         length(model.profile.T))
        vSmartMOM.Absorption.compute_τ_cia!(
            expected, table, model.profile, model.profile.vmr)

        @test all(expected .> 0)
        @test model.τ_abs[1] ≈ expected rtol=1e-14
    end
end

@testset "observer-height anchored reduction" begin
    params = parameters_from_yaml(_MS_CONFIG)
    params.obs_alt = [2.0, 6.0]
    params.profile_reduction_n = 1

    model = @test_logs (:warn, r"increasing layer count") model_from_parameters(params)
    @test length(model.profile.T) == 3  # N=2 interior H_i requires N+1 layers
    @test model.obs_geom.sensor_altitudes == [6.0, 2.0]
    @test CoreRT.half_level_altitudes(model.profile)[
        model.obs_geom.sensor_levels .+ 1] ≈ [6.0, 2.0] atol=1e-7

    # K interior interfaces still require K+1 physical layers when the
    # requested reduction target equals K.
    equal_target = deepcopy(params)
    equal_target.profile_reduction_n = 2
    equal_model = @test_logs (:warn, r"increasing layer count") model_from_parameters(equal_target)
    @test length(equal_model.profile.T) == 3

    # A reduction request larger than the naturally H-framed profile must not
    # expand it. The two-layer quickstart column needs only one split at 5 km.
    unexpanded = _height_model(5.0; reduction=20)
    @test length(unexpanded.profile.T) == 3
    @test CoreRT.half_level_altitudes(unexpanded.profile)[
        only(unexpanded.obs_geom.sensor_levels) + 1] ≈ 5.0 atol=1e-7

    @test_throws ArgumentError CoreRT.prepare_observer_profile(
        Float64[250, 275], Float64[100, 500, 1000], zeros(2), Dict(),
        Float64[2, 6], 0)
    @test_throws ArgumentError CoreRT.prepare_observer_profile(
        Float64[250, 275], Float64[100, 500, 1000], zeros(2), Dict(),
        Float64[2, 6], -2)
end

@testset "observer-height vertical source reframing" begin
    params = parameters_from_yaml(_MS_CONFIG)
    params.obs_alt = 5.0

    # B is supplied on the original two-layer input grid. Model construction
    # must interpolate it to the three-layer grid created by inserting H.
    B_input = reshape([1.0, 3.0], 2, 1)
    sources = SolarBeam() + ThermalEmission(B_layer=B_input)
    model = model_from_parameters(params; sources=sources)
    thermal = only(filter(source -> source isa ThermalEmission,
                          model.sources.sources))
    @test size(thermal.B_layer) == (length(model.profile.T), 1) == (3, 1)

    input_p_full = (params.p[1:end-1] .+ params.p[2:end]) ./ 2
    expected = CoreRT._interp_layer_field(
        input_p_full, vec(B_input), model.profile.p_full)
    @test vec(thermal.B_layer) ≈ expected

    # Batch updates use the same anchored builder and interpolate from a stable
    # copy of the original source values, so they match a fresh scene build and
    # do not accumulate interpolation drift.
    ctx = BatchContext(params; sources=sources)
    params_new = deepcopy(params)
    params_new.T .+= 2.0
    params_new.p .*= 1.01
    fresh = model_from_parameters(params_new; sources=sources)
    update_model!(ctx; T=params_new.T, p_half=params_new.p, q=params_new.q)

    thermal_ctx = only(filter(source -> source isa ThermalEmission,
                              ctx.model.sources.sources))
    thermal_fresh = only(filter(source -> source isa ThermalEmission,
                                fresh.sources.sources))
    @test thermal_ctx.B_layer ≈ thermal_fresh.B_layer rtol=0 atol=0
    @test ctx.model.profile.p_full ≈ fresh.profile.p_full rtol=0 atol=0
    @test ctx.model.geometry.sensor_levels == fresh.geometry.sensor_levels
    @test ctx.model.geometry.toa_altitude ≈ fresh.geometry.toa_altitude rtol=0 atol=0
    @test CoreRT.half_level_altitudes(ctx.model.profile)[
        only(ctx.model.geometry.sensor_levels) + 1] ≈ 5.0 atol=1e-7

    raw_thermal = only(filter(source -> source isa ThermalEmission,
                              ctx.input_sources.sources))
    @test raw_thermal.B_layer == B_input
end

@testset "observer-height output convention" begin
    endpoints_model = _height_model([0.0])
    endpoint_column = CoreRT._rt_run_column(
        vSmartMOM.InelasticScattering.noRS{Float64}(), endpoints_model, [1])
    endpoints = rt_run(endpoints_model)
    @test endpoints isa ObserverRTResult
    @test isempty(endpoints.levels)
    @test endpoints.toa == endpoint_column[1]
    @test endpoints.boa == endpoint_column[2]
    @test endpoints[1] === endpoints.toa
    @test endpoints[end] === endpoints.bhr_dw
    @test Tuple(endpoints)[1:2] == (endpoints.toa, endpoints.boa)
    R, T = endpoints
    @test R == endpoints.toa
    @test T == endpoints.boa

    boa = rt_run(_height_model(0.0))
    @test boa.toa === nothing
    @test boa.boa !== nothing
    @test isempty(boa.levels)

    toa_model = _height_model(1.0e6)
    toa = rt_run(toa_model)
    @test toa.toa !== nothing
    @test toa.boa === nothing
    @test isempty(toa.levels)

    interior_model = _height_model(5.0)
    interior = rt_run(interior_model)
    @test interior.toa === nothing
    @test interior.boa === nothing
    @test length(interior.levels) == 1
    @test interior.levels[1].height_km == 5.0
    @test size(interior.levels[1].upwelling) == (1, 1, 1)
    @test size(interior.levels[1].downwelling) == (1, 1, 1)

    # Exporting full-column stream state must not silently erase the separately
    # requested interior result. The dedicated stream API remains usable on the
    # same H-framed model without performing that second observer solve.
    callback_hits = Ref(0)
    with_callback = CoreRT.rt_run(
        vSmartMOM.InelasticScattering.noRS{Float64}(), interior_model, [1];
        streams_callback = _ -> (callback_hits[] += 1))
    @test callback_hits[] > 0
    @test only(with_callback.levels).upwelling ≈ only(interior.levels).upwelling
    @test only(with_callback.levels).downwelling ≈ only(interior.levels).downwelling
    streams = CoreRT.rt_run_streams(interior_model)
    @test !isempty(streams.J⁻_per_m)
    @test size(streams.τ_total, 2) == length(interior_model.profile.T)

    arbitrary = rt_run(_height_model([3.0, 8.0]))
    @test arbitrary.toa === nothing
    @test arbitrary.boa === nothing
    @test getfield.(arbitrary.levels, :height_km) == [8.0, 3.0]
    @test all(level -> size(level.upwelling) == (1, 1, 1), arbitrary.levels)
    @test all(level -> size(level.downwelling) == (1, 1, 1), arbitrary.levels)

    combined_model = _height_model([0.0, 5.0])
    combined_column = CoreRT._rt_run_column(
        vSmartMOM.InelasticScattering.noRS{Float64}(), combined_model, [1])
    combined = rt_run(combined_model)
    @test combined.toa ≈ combined_column[1]
    @test combined.boa ≈ combined_column[2]
    @test length(combined.levels) == 1
    @test combined.levels[1].height_km == 5.0
    @test combined.levels[1].upwelling ≈ interior.levels[1].upwelling
    @test combined.levels[1].downwelling ≈ interior.levels[1].downwelling
end

@testset "full-column-only solver observer semantics" begin
    interior_model = _height_model(5.0)

    ss_error = try
        CoreRT.rt_run_ss(interior_model; i_band=1)
        nothing
    catch err
        err
    end
    @test ss_error isa ArgumentError
    @test sprint(showerror, ss_error) ==
          "ArgumentError: rt_run_ss does not support interior-height radiances; " *
          "use rt_run(model) for obs_alt requests inside the atmosphere"

    lin_error = try
        rt_run(interior_model, nothing, 0, 0, 0; i_band=1)
        nothing
    catch err
        err
    end
    @test lin_error isa ArgumentError
    @test sprint(showerror, lin_error) ==
          "ArgumentError: linearized rt_run does not support interior-height radiances; " *
          "use rt_run(model) for obs_alt requests inside the atmosphere"

    lin_params = parameters_from_yaml(_MS_CONFIG)
    lin_params.obs_alt = 0.0
    lin_model, jac_model = model_from_parameters(LinMode(), lin_params)
    lin_boa = rt_run(lin_model, jac_model, 0, 0, 1; i_band=1)
    @test lin_boa[1] === nothing
    @test lin_boa[2] !== nothing
    @test lin_boa[3] === nothing
    @test lin_boa[4] !== nothing

    # Single-scatter still computes a full column internally, but its public
    # slots honor the same resolved scalar endpoint convention as rt_run.
    boa_ss = CoreRT.rt_run_ss(_height_model(0.0); i_band=1)
    @test boa_ss[1] === nothing
    @test boa_ss[2] !== nothing
    @test boa_ss[3] === nothing
    @test boa_ss[4] !== nothing
    @test boa_ss[5] === nothing
    @test boa_ss[6] !== nothing

    toa_ss = CoreRT.rt_run_ss(_height_model(1.0e6); i_band=1)
    @test toa_ss[1] !== nothing
    @test toa_ss[2] === nothing
    @test toa_ss[3] !== nothing
    @test toa_ss[4] === nothing
    @test toa_ss[5] !== nothing
    @test toa_ss[6] === nothing
end

@testset "multisensor explicit solar carrier" begin
    model = _height_model(5.0)
    rs = vSmartMOM.InelasticScattering.noRS{Float64}()
    rs.F₀ = fill(3.0, 1, 1)
    CoreRT._prepare_multisensor_F₀!(rs, model, [1], nothing)
    @test rs.F₀ == fill(3.0, 1, 1)
end

@testset "solar-aligned unscattered observer radiance" begin
    # The direct carrier follows the solver's ordinate convention: it is
    # present whenever a requested VZA maps to the solar quadrature node,
    # independent of VAZ. A nearby but separately inserted VZA node must not
    # receive any of the delta-function beam.
    vza = [60.0, 60.0, 60.0, 45.0, 59.9]
    vaz = [0.0, 90.0, 180.0, 180.0, 180.0]

    for FT in (Float64, Float32)
        model = _direct_height_model(FT, [0.0, 1.0e-6]; vza=vza, vaz=vaz)
        F₀ = reshape(FT[3], 1, 1)
        result = rt_run(model; sources=SolarBeam(F₀=F₀))
        level = only(result.levels)
        direct = level.unscattered_downwelling
        expected = _expected_unscattered(model, level.boundary_index, F₀[1])
        direct_rtol = FT === Float64 ? 5e-13 : 5e-6

        @test eltype(direct) === FT
        @test size(direct) == size(level.downwelling) == (length(vza), 1, 1)
        @test vec(direct[1:3, 1, 1]) ≈ fill(expected, 3) rtol=direct_rtol
        @test direct[1, 1, 1] == direct[2, 1, 1] == direct[3, 1, 1]
        @test all(iszero, @view direct[4:5, :, :])

        # At one millimetre above BOA, diffuse plus separately reported
        # unscattered radiance approaches the legacy total BOA radiance.
        total_down = level.downwelling[1, 1, 1] + direct[1, 1, 1]
        boundary_rtol = FT === Float64 ? 1e-8 : 2e-5
        @test total_downwelling(level) == level.downwelling .+ direct
        @test total_down ≈ result.boa[1, 1, 1] rtol=boundary_rtol

        # The upper-boundary limit independently verifies attenuation and the
        # existing upwelling solution. Float32 uses a 10 cm offset so the
        # requested interface remains distinct at a ~17 km TOA altitude.
        toa_offset = FT === Float64 ? 1.0e-6 : 1.0e-4
        near_toa_model = _direct_height_model(
            FT, [0.0, Float64(model.obs_geom.toa_altitude) - toa_offset])
        near_toa_result = rt_run(
            near_toa_model; sources=SolarBeam(F₀=F₀))
        near_toa = only(near_toa_result.levels)
        near_toa_direct = near_toa.unscattered_downwelling[1, 1, 1]
        expected_near_toa = _expected_unscattered(
            near_toa_model, near_toa.boundary_index, F₀[1])

        @test near_toa_direct ≈ expected_near_toa rtol=direct_rtol
        @test near_toa_direct ≈ F₀[1] / (2π) rtol=(FT === Float64 ? 1e-8 : 2e-6)
        @test near_toa.upwelling ≈ near_toa_result.toa rtol=(FT === Float64 ? 1e-8 : 2e-6)
    end

    # An explicit source spectrum scales the unscattered output, rather than
    # silently falling back to the historical unit Stokes-I carrier.
    scale_model = _direct_height_model(Float64, 5.0)
    unit = rt_run(scale_model; sources=SolarBeam(F₀=ones(1, 1)))
    triple = rt_run(scale_model; sources=SolarBeam(F₀=fill(3.0, 1, 1)))
    @test only(triple.levels).unscattered_downwelling ≈
          3 .* only(unit.levels).unscattered_downwelling rtol=5e-14

    # Preserve the established m=0 Stokes selector used by endpoint output:
    # incident I/Q survive while U/V do not contribute to this Fourier slot.
    polarized_model = _direct_height_model(
        Float64, 5.0; polarization=vSmartMOM.Scattering.Stokes_IQUV())
    polarized_F₀ = reshape([2.0, -0.5, 0.25, -0.125], 4, 1)
    polarized = only(rt_run(
        polarized_model; sources=SolarBeam(F₀=polarized_F₀)).levels)
    attenuation = _expected_unscattered(
        polarized_model, polarized.boundary_index, 1.0)
    @test vec(polarized.unscattered_downwelling[1, :, 1]) ≈
          attenuation .* [2.0, -0.5, 0.0, 0.0] rtol=5e-13 atol=1e-15

    # Preserve source compatibility for downstream callers that construct a
    # LevelRadiance with the original six arguments.
    legacy_down = fill(2.0, 1, 1, 1)
    legacy_level = LevelRadiance(
        5.0, 1, ones(1, 1, 1), legacy_down, nothing, nothing)
    @test all(iszero, legacy_level.unscattered_downwelling)
    @test total_downwelling(legacy_level) == legacy_down
end

@testset "observer-height display" begin
    model = _height_model([0.0, 5.0])
    model_text = sprint(show, MIME("text/plain"), model)
    params = parameters_from_yaml(_MS_CONFIG)
    params.obs_alt = [0.0, 5.0]
    params_text = sprint(show, MIME("text/plain"), params)
    result_text = sprint(show, MIME("text/plain"), rt_run(model))

    @test occursin("outputs=TOA, 5.0 km, BOA", model_text)
    @test occursin("observer: TOA + BOA + requested [5.0] km", params_text)
    @test occursin("interior levels: 1", result_text)
    @test occursin("5.0 km", result_text)
end
