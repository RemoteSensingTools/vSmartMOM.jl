using Test
using vSmartMOM
using vSmartMOM.StandaloneSS

function _as_vector(x, n_spec, ::Type{FT}) where {FT}
    x isa Number && return fill(convert(FT, x), n_spec)
    return convert(Vector{FT}, collect(x))
end

_ref_phase(::RayleighSSContributor, cosΘ) = 0.75 * (1 + cosΘ^2)
function _ref_phase(c::HGAerosolSSContributor, cosΘ)
    g = c.g
    (1 - g^2) / (1 + g^2 - 2g * cosΘ)^1.5
end
_ref_phase(::AbsorptionSSContributor, cosΘ) = zero(cosΘ)

_ref_ssa(::RayleighSSContributor) = 1
_ref_ssa(c::HGAerosolSSContributor) = c.ϖ
_ref_ssa(::AbsorptionSSContributor) = 0

function _ref_cos_scatter(μ₀, μv, Δϕ)
    -μ₀ * μv + sqrt(max(0, 1 - μ₀^2)) * sqrt(max(0, 1 - μv^2)) * cos(Δϕ)
end

function _ref_outputs(config; paths=:paths_1_2)
    contributors = config.contributors
    n_layers, n_spec = size(first(contributors).τ)
    n_geom = length(config.geometry.μv)
    I0 = _as_vector(config.I0, n_spec, typeof(config.geometry.μ₀))
    albedo = _as_vector(config.surface.albedo, n_spec, typeof(config.geometry.μ₀))

    τ_total_layer = zeros(typeof(config.geometry.μ₀), n_layers, n_spec)
    τ_scat_layer = zeros(typeof(config.geometry.μ₀), n_layers, n_spec)
    weighted_phase = zeros(typeof(config.geometry.μ₀), n_geom, n_layers, n_spec)

    for c in contributors
        for s in 1:n_spec, iz in 1:n_layers
            τ = c.τ[iz, s]
            ϖ = _ref_ssa(c)
            τ_total_layer[iz, s] += τ
            τ_scat_layer[iz, s] += τ * ϖ
            for iv in 1:n_geom
                cosΘ = _ref_cos_scatter(config.geometry.μ₀,
                                        config.geometry.μv[iv],
                                        config.geometry.Δϕ[iv])
                weighted_phase[iv, iz, s] += τ * ϖ * _ref_phase(c, cosΘ)
            end
        end
    end

    τ_cum = zeros(typeof(config.geometry.μ₀), n_layers + 1, n_spec)
    for s in 1:n_spec, iz in 1:n_layers
        τ_cum[iz + 1, s] = τ_cum[iz, s] + τ_total_layer[iz, s]
    end

    path1 = zeros(typeof(config.geometry.μ₀), n_geom, 1, n_spec)
    path2 = zeros(typeof(config.geometry.μ₀), n_geom, 1, n_spec)

    if paths in (:path1, :paths_1_2)
        for iv in 1:n_geom, s in 1:n_spec
            μ₀ = config.geometry.μ₀
            μv = config.geometry.μv[iv]
            a = 1 / μ₀ + 1 / μv
            prefactor = I0[s] / (4π * μv * a)
            for iz in 1:n_layers
                τ_scat = τ_scat_layer[iz, s]
                τ_total = τ_total_layer[iz, s]
                τ_total == 0 && continue
                τ_scat == 0 && continue
                ϖ_eff = τ_scat / τ_total
                P_eff = weighted_phase[iv, iz, s] / τ_scat
                layer_factor = exp(-τ_cum[iz, s] * a) -
                               exp(-τ_cum[iz + 1, s] * a)
                path1[iv, 1, s] += prefactor * ϖ_eff * P_eff * layer_factor
            end
        end
    end

    if paths in (:path2, :paths_1_2)
        for iv in 1:n_geom, s in 1:n_spec
            μ₀ = config.geometry.μ₀
            μv = config.geometry.μv[iv]
            τ = τ_cum[end, s]
            path2[iv, 1, s] =
                (μ₀ * I0[s] * albedo[s] / π) * exp(-τ / μ₀) * exp(-τ / μv)
        end
    end

    (; path1, path2, total=path1 .+ path2)
end

function _mixed_config(::Type{FT}) where {FT}
    geometry = SSGeometry(μ₀=FT(0.82),
                          μv=FT[0.35, 0.74],
                          Δϕ=FT[0.2, 1.4])
    surface = LambertianSSSurface(albedo=FT[0.2, 0.35])
    rayleigh = RayleighSSContributor(τ=FT[0.08 0.03; 0.11 0.05])
    aerosol = HGAerosolSSContributor(g=FT(0.72), ϖ=FT(0.91),
                                     τ=FT[0.04 0.02; 0.07 0.01])
    absorption = AbsorptionSSContributor(τ=FT[0.01 0.05; 0.02 0.04])
    ExactSSConfig(geometry=geometry, surface=surface,
                  contributors=(rayleigh, aerosol, absorption),
                  I0=FT[1.0, 0.8])
end

@testset "StandaloneSS Phase 1" begin
    @testset "Top-level exports" begin
        @test run_exact_ss === vSmartMOM.run_exact_ss
        @test StandaloneSS === vSmartMOM.StandaloneSS
    end

    @testset "Rayleigh path 1 and Lambertian path 2" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.8), μv=FT[0.4, 0.7], Δϕ=FT[0.0, 0.9])
        surface = LambertianSSSurface(albedo=FT(0.3))
        rayleigh = RayleighSSContributor(τ=FT[0.05 0.07; 0.1 0.04])
        config = ExactSSConfig(geometry=geometry, surface=surface,
                               contributors=(rayleigh,), I0=FT[1.0, 0.7])

        result = run_exact_ss(config; paths=:paths_1_2)
        expected = _ref_outputs(config; paths=:paths_1_2)
        @test size(result.total) == (2, 1, 2)
        @test result.path1 ≈ expected.path1 rtol=1e-12 atol=1e-14
        @test result.path2 ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test result.total ≈ expected.total rtol=1e-12 atol=1e-14
        @test result.metadata.n_stokes == 1
        @test result.quadrature_info.kernel_backend == :KernelAbstractionsCPU
    end

    @testset "Mixed τϖ-weighted phase mixing" begin
        config = _mixed_config(Float64)
        result = run_exact_ss(config)
        expected = _ref_outputs(config)
        @test result.path1 ≈ expected.path1 rtol=1e-12 atol=1e-14
        @test result.path2 ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test all(isfinite.(result.total))
        @test all(result.total .> 0)
    end

    @testset "Path selection" begin
        config = _mixed_config(Float64)
        only1 = run_exact_ss(config; paths=:path1)
        only2 = run_exact_ss(config; paths=:path2)
        both = run_exact_ss(config; paths=:paths_1_2)
        @test only1.path2 == zero(only1.path2)
        @test only2.path1 == zero(only2.path1)
        @test both.total ≈ only1.total .+ only2.total rtol=1e-12
    end

    @testset "Black surface and pure absorption" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.7), μv=FT[0.5], Δϕ=FT[0.3])
        black = LambertianSSSurface(albedo=zero(FT))
        rayleigh = RayleighSSContributor(τ=reshape(FT[0.1, 0.2], 2, 1))
        black_config = ExactSSConfig(geometry=geometry, surface=black,
                                     contributors=(rayleigh,), I0=FT[1.0])
        black_result = run_exact_ss(black_config)
        @test black_result.path2 == zero(black_result.path2)
        @test black_result.path1[1, 1, 1] > 0

        absorption = AbsorptionSSContributor(τ=reshape(FT[0.2, 0.3], 2, 1))
        surface = LambertianSSSurface(albedo=FT(0.4))
        abs_config = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(absorption,), I0=FT[1.0])
        abs_result = run_exact_ss(abs_config)
        abs_expected = _ref_outputs(abs_config)
        @test abs_result.path1 == zero(abs_result.path1)
        @test abs_result.path2 ≈ abs_expected.path2 rtol=1e-12 atol=1e-14
        @test abs_result.path2[1, 1, 1] > 0
    end

    @testset "Float32" begin
        config = _mixed_config(Float32)
        result = run_exact_ss(config)
        expected = _ref_outputs(config)
        @test eltype(result.total) == Float32
        @test all(isfinite.(result.total))
        @test result.path1 ≈ expected.path1 rtol=5e-6 atol=1e-7
        @test result.path2 ≈ expected.path2 rtol=5e-6 atol=1e-7
    end

    @testset "Validation errors" begin
        FT = Float64
        good_surface = LambertianSSSurface(albedo=FT(0.2))
        c = RayleighSSContributor(τ=FT[0.1;;])
        bad_μ₀ = ExactSSConfig(
            geometry=SSGeometry(μ₀=FT(-0.1), μv=FT[0.5], Δϕ=FT[0.0]),
            surface=good_surface, contributors=(c,), I0=FT[1.0])
        bad_μv = ExactSSConfig(
            geometry=SSGeometry(μ₀=FT(0.8), μv=FT[0.0], Δϕ=FT[0.0]),
            surface=good_surface, contributors=(c,), I0=FT[1.0])
        good = ExactSSConfig(
            geometry=SSGeometry(μ₀=FT(0.8), μv=FT[0.5], Δϕ=FT[0.0]),
            surface=good_surface, contributors=(c,), I0=FT[1.0])

        @test_throws ArgumentError run_exact_ss(bad_μ₀)
        @test_throws ArgumentError run_exact_ss(bad_μv)
        @test_throws ArgumentError run_exact_ss(good; paths=:all)
        @test_throws ArgumentError run_exact_ss(ExactSSConfig(
            geometry=good.geometry, surface=good.surface,
            contributors=good.contributors, I0=good.I0, n_stokes=3))
    end

    @testset "Required-Nquad diagnostics" begin
        config = _mixed_config(Float64)
        rayleigh, aerosol, absorption = config.contributors

        @test determine_required_l_from_moments([1.0]) == 0
        @test determine_required_l_from_moments([0.0, 0.0, 0.0]) == 0
        @test determine_required_l_aerosol(rayleigh) == 2
        @test determine_required_l_aerosol(absorption) == 0
        @test determine_required_l_aerosol(aerosol; target_relative_error=1e-4) ==
              ceil(Int, log(1e-4) / log(abs(aerosol.g)))
        @test determine_required_l_aerosol(
            HGAerosolSSContributor(g=0.0, ϖ=1.0, τ=[1.0;;])) == 0

        @test determine_required_nbrdf(config.surface) == 1
        @test determine_required_nbrdf_coxmunk(0.0; target_relative_error=1e-4) ==
              clamp(ceil(Int, log(1 / 1e-4) / sqrt(0.003)), 16, 200)
        @test determine_required_nbrdf_coxmunk(1000.0; target_relative_error=1e-4) == 16

        stream_diag = determine_required_nstreams(config.contributors, config.surface)
        @test stream_diag.L_aerosol ==
              determine_required_l_aerosol(aerosol; target_relative_error=1e-4)
        @test stream_diag.Nstreams == ceil(Int, (stream_diag.L_aerosol + 1) / 2)
        @test stream_diag.NSTREAMS_BRDF == 1
        @test stream_diag.surface_kind == :LambertianSSSurface

        inner = determine_required_nquad_inner([0.1, 1.5], config.contributors)
        @test 8 <= inner <= 64
        @test inner >= determine_required_nquad_inner(0.0, (absorption,))

        combined = determine_required_nquad(config)
        τ_total = zeros(2)
        for c in config.contributors
            τ_total .+= vec(sum(c.τ; dims=1))
        end
        @test combined.Nstreams == stream_diag.Nstreams
        @test combined.NSTREAMS_BRDF == stream_diag.NSTREAMS_BRDF
        @test combined.inner_nquad == determine_required_nquad_inner(τ_total,
                                                                     config.contributors)

        @test_throws ArgumentError determine_required_l_from_moments(Float64[])
        @test_throws ArgumentError determine_required_l_from_moments([1.0]; target_relative_error=0.0)
        @test_throws ArgumentError determine_required_l_aerosol(
            HGAerosolSSContributor(g=1.0, ϖ=1.0, τ=[1.0;;]))
        @test_throws ArgumentError determine_required_nbrdf_coxmunk(-1.0)
        @test_throws ArgumentError determine_required_nquad_inner(-0.1, config.contributors)
    end
end
