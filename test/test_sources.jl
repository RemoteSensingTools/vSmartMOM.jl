# =========================================================================
# v0.6 source-term refactor — Phase 1 + Phase 2 regression tests
#
# Phase 1: AbstractSource type vocabulary (composition, NoSource identity,
#          AD-mode traits).
# Phase 2: SolarBeam + PreparedSolarBeam + prepare_source seam, plus the
#          rt_run(_lin/_ss) `sources=` kwarg routed through the legacy
#          RS_type.F₀ channel — verified bit-equal against the unchanged
#          default path.
#
# Plan:  ~/.claude/plans/gpt-also-had-some-velvety-whale.md
# =========================================================================

using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.InelasticScattering
using Test

@testset "v0.6 sources — vocabulary + SolarBeam routing" begin

    @testset "Phase 1: AbstractSource composition" begin
        struct _StubSrc <: AbstractSource end
        a = _StubSrc(); b = _StubSrc(); c = _StubSrc()

        # NoSource is the additive identity
        @test (NoSource() + NoSource()) === NoSource()
        @test (NoSource() + a) === a
        @test (a + NoSource()) === a

        # Concrete + concrete → SourceSet of length 2
        s_ab = a + b
        @test s_ab isa SourceSet
        @test length(s_ab) == 2

        # Associativity: flattened (no nesting)
        @test length((a + b) + c) == 3
        @test length(a + (b + c)) == 3

        # SourceSet × SourceSet flattens
        @test length(s_ab + s_ab) == 4

        # NoSource × SourceSet routes through identity
        @test (NoSource() + s_ab) === s_ab
        @test (s_ab + NoSource()) === s_ab

        # Iteration / indexing
        n = 0
        for src in s_ab
            n += 1
            @test src isa _StubSrc
        end
        @test n == 2
        @test s_ab[1] isa _StubSrc
        @test firstindex(s_ab) == 1
        @test lastindex(s_ab) == 2

        # AD-mode traits
        @test source_ad_mode(NoSource()) isa NoSourceJacobian
        @test source_ad_mode(a) isa AnalyticSourceJacobian
        @test source_ad_mode(s_ab) isa AnalyticSourceJacobian
    end

    @testset "Phase 2: SolarBeam construction + prepare_source" begin
        sb_default = SolarBeam()
        @test sb_default isa SolarBeam
        @test sb_default.F₀ === nothing
        @test sb_default.sza === nothing

        sb_sza = SolarBeam(sza = 35.0)
        @test sb_sza.sza == 35.0

        F₀_in = [1.0 2.0 3.0; 0.0 0.0 0.0]
        sb_F₀ = SolarBeam(F₀ = F₀_in)
        @test sb_F₀.F₀ === F₀_in

        # AD mode trait
        @test source_ad_mode(sb_default) isa AnalyticSourceJacobian
        @test source_ad_mode(prepare_source(sb_default, Float64, 3, 5, identity)) isa
              AnalyticSourceJacobian

        # SolarBeam + SolarBeam composes
        combo = sb_default + sb_sza
        @test combo isa SourceSet
        @test length(combo) == 2
    end

    @testset "Phase 2: prepare_source materialises F₀ correctly" begin
        # default ⇒ unit Stokes-I
        prepared = prepare_source(SolarBeam(), Float64, 3, 5, identity)
        @test prepared isa PreparedSolarBeam
        @test size(prepared.F₀) == (3, 5)
        @test all(prepared.F₀[1, :] .== 1.0)
        @test all(prepared.F₀[2, :] .== 0.0)

        # custom F₀ round-trips
        F₀_custom = randn(3, 5)
        prepared2 = prepare_source(SolarBeam(F₀ = F₀_custom), Float64, 3, 5, identity)
        @test prepared2.F₀ ≈ F₀_custom

        # FT precision conversion
        prepared3 = prepare_source(SolarBeam(F₀ = F₀_custom), Float32, 3, 5, identity)
        @test eltype(prepared3.F₀) === Float32
        @test prepared3.F₀ ≈ Float32.(F₀_custom)

        # Shape mismatch errors out (rather than broadcasting silently)
        sb_bad = SolarBeam(F₀ = [1.0 2.0; 3.0 4.0])
        @test_throws ErrorException prepare_source(sb_bad, Float64, 3, 5, identity)
    end

    @testset "Phase 2: prepare_sources walks SourceSet, NoSource round-trips" begin
        @test prepare_sources(NoSource(), Float64, 3, 5, identity) === NoSource()

        prep_one = prepare_sources(SolarBeam(), Float64, 3, 5, identity)
        @test prep_one isa SourceSet
        @test length(prep_one) == 1
        @test prep_one.sources[1] isa PreparedSolarBeam

        prep_two = prepare_sources(SolarBeam() + SolarBeam(sza = 30.0), Float64, 3, 5, identity)
        @test prep_two isa SourceSet
        @test length(prep_two) == 2
        @test all(s -> s isa PreparedSolarBeam, prep_two.sources)
    end

    @testset "Phase 2: extract_solar_F₀ + has_solar_beam helpers" begin
        prepared = prepare_source(SolarBeam(), Float64, 3, 5, identity)
        prep_set = prepare_sources(SolarBeam() + SolarBeam(sza = 30.0), Float64, 3, 5, identity)

        # has_solar_beam dispatch
        @test CoreRT.has_solar_beam(prepared)
        @test CoreRT.has_solar_beam(prep_set)
        @test !CoreRT.has_solar_beam(NoSource())

        # extract_solar_F₀ pulls from the first SolarBeam in the set
        F₀_extracted = CoreRT.extract_solar_F₀(prep_set, Float64, 3, 5, identity)
        @test F₀_extracted === prep_set.sources[1].F₀

        # extract_solar_F₀(NoSource) ⇒ zeros (no solar)
        F₀_zero = CoreRT.extract_solar_F₀(NoSource(), Float64, 3, 5, identity)
        @test size(F₀_zero) == (3, 5)
        @test all(F₀_zero .== 0.0)
    end

    @testset "Phase 2c: sources lives on RTModel (model.sources field)" begin
        # Model field defaults to a SolarBeam (matches rt_run's historical
        # unit-Stokes-I default).
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = vSmartMOM.Architectures.CPU()

        model_default = model_from_parameters(params)
        @test model_default.sources isa SolarBeam

        # `model_from_parameters(...; sources=...)` overrides the default.
        F₀_custom = nothing  # SolarBeam() with default keeps unit Stokes I
        custom = SolarBeam()
        model_explicit = model_from_parameters(params; sources = custom)
        @test model_explicit.sources === custom

        # Composed SourceSet round-trips into the model
        combo = SolarBeam() + SolarBeam(sza = 30.0)
        model_combo = model_from_parameters(params; sources = combo)
        @test model_combo.sources isa SourceSet
        @test length(model_combo.sources) == 2
    end

    @testset "Phase 2: rt_run(model) bit-equal to rt_run(model; sources=SolarBeam())" begin
        # End-to-end bit-equality check on a small CPU forward run. Default
        # path (sources=nothing) should match the explicit unit-Stokes-I
        # SolarBeam to the last bit, since both flow F₀ = ones(1, :) through
        # the same kernels.
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = vSmartMOM.Architectures.CPU()
        model = model_from_parameters(params)

        R_default = rt_run(model)[1]                                      # legacy
        R_explicit = rt_run(model; sources = SolarBeam())[1]              # default-equivalent
        R_explicit2 = rt_run(model; sources = SolarBeam(F₀ = nothing))[1] # same as above

        @test R_default == R_explicit
        @test R_default == R_explicit2

        # SourceSet wrapper around a default SolarBeam — same answer.
        R_in_set = rt_run(model; sources = SourceSet((SolarBeam(),)))[1]
        @test R_default == R_in_set

        # NoSource ⇒ F₀ all zero ⇒ TOA reflectance is zero (no incident solar).
        R_no_source = rt_run(model; sources = NoSource())[1]
        @test all(iszero, R_no_source)

        # A user-supplied F₀ scaled by 2× scales the linear-in-F₀ TOA
        # reflectance the same way (single-scattering and propagation are
        # both linear in F₀, even with multi-bounce surface coupling).
        FT = CoreRT.float_type(model)
        pol_n = CoreRT.polarization_type(model).n
        nSpec = sum(size(model.τ_abs[i], 1) for i in 1:length(model.τ_abs))
        F₀_2x = zeros(FT, pol_n, nSpec)
        F₀_2x[1, :] .= FT(2)
        R_2x = rt_run(model; sources = SolarBeam(F₀ = F₀_2x))[1]
        @test R_2x ≈ FT(2) .* R_default rtol = 1e-12
    end
end
