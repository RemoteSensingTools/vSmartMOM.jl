# =================================================================
# Type Stability Tests for vSmartMOM
# Verifies that core types and operations preserve Float32/Float64
# without unwanted promotion. Uses @inferred for key helper functions.
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using Test

println("="^60)
println("Type Stability Tests")
println("="^60)

@testset "Type Stability" begin

    @testset "CoreScatteringOpticalProperties construction" begin
        # Minimal arrays for inline test (no YAML)
        n = 2
        nSpec = 1

        # Float64
        τ64 = [0.1]
        ϖ64 = [0.99]
        Z⁺⁺64 = zeros(Float64, n, n, nSpec)
        Z⁺⁺64[1, 1, 1] = 1.0
        Z⁻⁺64 = zeros(Float64, n, n, nSpec)
        Z⁻⁺64[1, 1, 1] = 0.5

        opt64 = CoreRT.CoreScatteringOpticalProperties(τ=τ64, ϖ=ϖ64, Z⁺⁺=Z⁺⁺64, Z⁻⁺=Z⁻⁺64)
        @test opt64 isa CoreRT.CoreScatteringOpticalProperties
        @test eltype(opt64.τ) == Float64
        @test eltype(opt64.ϖ) == Float64
        @test eltype(opt64.Z⁺⁺) == Float64
        @test eltype(opt64.Z⁻⁺) == Float64

        # Float32
        τ32 = Float32[0.1f0]
        ϖ32 = Float32[0.99f0]
        Z⁺⁺32 = zeros(Float32, n, n, nSpec)
        Z⁺⁺32[1, 1, 1] = 1.0f0
        Z⁻⁺32 = zeros(Float32, n, n, nSpec)
        Z⁻⁺32[1, 1, 1] = 0.5f0

        opt32 = CoreRT.CoreScatteringOpticalProperties(τ=τ32, ϖ=ϖ32, Z⁺⁺=Z⁺⁺32, Z⁻⁺=Z⁻⁺32)
        @test opt32 isa CoreRT.CoreScatteringOpticalProperties
        @test eltype(opt32.τ) == Float32
        @test eltype(opt32.ϖ) == Float32
        @test eltype(opt32.Z⁺⁺) == Float32
        @test eltype(opt32.Z⁻⁺) == Float32
    end

    @testset "CoreScatteringOpticalProperties + and * preserve float type" begin
        n = 2
        nSpec = 1

        # Float32 pair
        τ32 = Float32[0.1f0]
        ϖ32 = Float32[0.99f0]
        Z⁺⁺32 = zeros(Float32, n, n, nSpec)
        Z⁺⁺32[1, 1, 1] = 1.0f0
        Z⁻⁺32 = zeros(Float32, n, n, nSpec)
        Z⁻⁺32[1, 1, 1] = 0.5f0

        opt_a = CoreRT.CoreScatteringOpticalProperties(τ=τ32, ϖ=ϖ32, Z⁺⁺=Z⁺⁺32, Z⁻⁺=Z⁻⁺32)
        opt_b = CoreRT.CoreScatteringOpticalProperties(τ=Float32[0.05f0], ϖ=Float32[0.98f0], Z⁺⁺=copy(Z⁺⁺32), Z⁻⁺=copy(Z⁻⁺32))

        # + should preserve Float32 (no promotion to Float64)
        sum_opt = opt_a + opt_b
        @test eltype(sum_opt.τ) == Float32
        @test eltype(sum_opt.ϖ) == Float32
        @test eltype(sum_opt.Z⁺⁺) == Float32
        @test eltype(sum_opt.Z⁻⁺) == Float32

        # * (vertical concatenation) should preserve Float32
        cat_opt = opt_a * opt_b
        @test eltype(cat_opt.τ) == Float32
        @test eltype(cat_opt.ϖ) == Float32
        @test eltype(cat_opt.Z⁺⁺) == Float32
        @test eltype(cat_opt.Z⁻⁺) == Float32

        # Float64 pair - same checks
        τ64 = [0.1]
        ϖ64 = [0.99]
        Z⁺⁺64 = zeros(Float64, n, n, nSpec)
        Z⁺⁺64[1, 1, 1] = 1.0
        Z⁻⁺64 = zeros(Float64, n, n, nSpec)
        Z⁻⁺64[1, 1, 1] = 0.5

        opt_a64 = CoreRT.CoreScatteringOpticalProperties(τ=τ64, ϖ=ϖ64, Z⁺⁺=Z⁺⁺64, Z⁻⁺=Z⁻⁺64)
        opt_b64 = CoreRT.CoreScatteringOpticalProperties(τ=[0.05], ϖ=[0.98], Z⁺⁺=copy(Z⁺⁺64), Z⁻⁺=copy(Z⁻⁺64))

        sum_opt64 = opt_a64 + opt_b64
        @test eltype(sum_opt64.τ) == Float64
        @test eltype(sum_opt64.ϖ) == Float64

        cat_opt64 = opt_a64 * opt_b64
        @test eltype(cat_opt64.τ) == Float64
        @test eltype(cat_opt64.ϖ) == Float64

        # Vacuum layers are valid: zero scattering plus zero absorption should
        # remain finite and non-scattering, not create NaN SSA through 0 / 0.
        zero_opt = CoreRT.CoreScatteringOpticalProperties(
            τ = zeros(Float64, nSpec),
            ϖ = zeros(Float64, nSpec),
            Z⁺⁺ = copy(Z⁺⁺64),
            Z⁻⁺ = copy(Z⁻⁺64),
        )
        zero_abs = CoreRT.CoreAbsorptionOpticalProperties(τ = zeros(Float64, nSpec))
        vacuum_opt = zero_opt + zero_abs
        @test all(vacuum_opt.τ .== 0)
        @test all(vacuum_opt.ϖ .== 0)
        @test all(isfinite.(vacuum_opt.ϖ))
    end

    @testset "Helper functions type-stable (@inferred)" begin
        # get_indices(iμ, pol_type) -> (st_iμ, istart, iend)
        pol_I = Scattering.Stokes_I()  # Float64 default
        pol_I32 = Scattering.Stokes_I(D=Float32[1.0], I₀=Float32[1.0])
        @test @inferred CoreRT.get_indices(1, pol_I) === (0, 1, 1)
        @test @inferred CoreRT.get_indices(2, pol_I) === (1, 2, 2)
        @test @inferred CoreRT.get_indices(1, pol_I32) === (0, 1, 1)

        pol_IQU = Scattering.Stokes_IQU()  # Float64 default
        @test @inferred CoreRT.get_indices(1, pol_IQU) === (0, 1, 3)
        @test @inferred CoreRT.get_indices(2, pol_IQU) === (3, 4, 6)

        # doubling_number(dτ_max, τ_end) -> (dτ, ndoubl)
        @test @inferred CoreRT.doubling_number(Float64(0.01), Float64(0.005)) === (0.005, 0)
        # Float32: verify correct behavior (may have type instability from 10.0^x in implementation)
        dτ32, ndoubl32 = CoreRT.doubling_number(Float32(0.01), Float32(0.005))
        @test dτ32 === Float32(0.005)
        @test ndoubl32 === 0
        dτ, ndoubl = CoreRT.doubling_number(Float32(0.01), Float32(1.0))
        @test dτ isa Float32
        @test ndoubl isa Int

        # nearest_point(f_array, f) -> Int (or CartesianIndex for multi-dim)
        f_arr32 = Float32[0.1f0, 0.2f0, 0.3f0, 0.4f0]
        @test @inferred CoreRT.nearest_point(f_arr32, 0.25f0) isa Integer
        f_arr64 = [0.1, 0.2, 0.3, 0.4]
        @test @inferred CoreRT.nearest_point(f_arr64, 0.25) isa Integer
    end

    @testset "Raman noRS keyword constructors infer float type" begin
        rs_default = vSmartMOM.InelasticScattering.noRS()
        @test rs_default isa vSmartMOM.InelasticScattering.noRS{Float64}
        @test eltype(rs_default.F₀) === Float64

        rs32 = vSmartMOM.InelasticScattering.noRS(
            fscattRayl = Float32[1],
            ϖ_Cabannes = Float32[1],
            bandSpecLim = UnitRange{Int}[],
            iBand = [1],
            F₀ = zeros(Float32, 1, 2),
        )
        @test rs32 isa vSmartMOM.InelasticScattering.noRS{Float32}
        @test eltype(rs32.F₀) === Float32

        rs_plus = vSmartMOM.InelasticScattering.noRS_plus(ϖ_Cabannes = Float32(1))
        @test rs_plus isa vSmartMOM.InelasticScattering.noRS_plus{Float32}
        @test rs_plus.ϖ_Cabannes === Float32(1)
    end

    @testset "AddedLayer and CompositeLayer with Float32" begin
        RS_type = vSmartMOM.InelasticScattering.noRS{Float32}()
        FT = Float32
        arr_type = Array
        dims = (2, 2)
        nSpec = 1

        added = CoreRT.make_added_layer(RS_type, FT, arr_type, dims, nSpec)
        @test added isa CoreRT.AddedLayer
        @test eltype(added.r⁻⁺) == Float32
        @test eltype(added.t⁺⁺) == Float32
        @test eltype(added.j₀⁺) == Float32

        composite = CoreRT.make_composite_layer(RS_type, FT, arr_type, dims, nSpec)
        @test composite isa CoreRT.CompositeLayer
        @test eltype(composite.R⁻⁺) == Float32
        @test eltype(composite.T⁺⁺) == Float32
        @test eltype(composite.J₀⁺) == Float32
    end

    @testset "AddedLayer and CompositeLayer with Float64" begin
        RS_type = vSmartMOM.InelasticScattering.noRS{Float64}()
        FT = Float64
        arr_type = Array
        dims = (2, 2)
        nSpec = 1

        added = CoreRT.make_added_layer(RS_type, FT, arr_type, dims, nSpec)
        @test added isa CoreRT.AddedLayer
        @test eltype(added.r⁻⁺) == Float64

        composite = CoreRT.make_composite_layer(RS_type, FT, arr_type, dims, nSpec)
        @test composite isa CoreRT.CompositeLayer
        @test eltype(composite.R⁻⁺) == Float64
    end
end

println("\n" * "="^60)
println("Type stability tests complete.")
println("="^60)
