using Test
using vSmartMOM
using vSmartMOM.StandaloneSS
using ForwardDiff
using LinearAlgebra

const _STANDALONE_CUDA_LOADED = Ref(false)
try
    @eval using CUDA
    _STANDALONE_CUDA_LOADED[] = true
catch
end

_standalone_cuda_functional() =
    _STANDALONE_CUDA_LOADED[] && CUDA.functional() &&
    vSmartMOM.Architectures.has_cuda()

_want_ref_path1(paths) = paths in (:path1, :paths_1_2, :all, :all_four)
_want_ref_path2(paths) = paths in (:path2, :paths_1_2, :all, :all_four)
_want_ref_path3(paths) = paths in (:path3, :all, :all_four)
_want_ref_path4(paths) = paths in (:path4, :all, :all_four)

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

function _ref_gauss_legendre_01(n, ::Type{FT}) where {FT}
    β = [FT(k) / sqrt(FT(4k^2 - 1)) for k in 1:(n - 1)]
    T = SymTridiagonal(zeros(FT, n), β)
    eig = eigen(T)
    x = FT(0.5) .* (eig.values .+ one(FT))
    w = FT(0.5) .* (FT(2) .* eig.vectors[1, :].^2)
    return x, w
end

function _ref_azimuthal_average_phase(c, μa, μb, n_phi)
    a = μa * μb
    b = sqrt(max(0, 1 - μa^2)) * sqrt(max(0, 1 - μb^2))
    s = zero(promote_type(typeof(μa), typeof(μb)))
    for k in 0:(n_phi - 1)
        ϕ = 2π * k / n_phi
        cosΘ = clamp(a + b * cos(ϕ), -1, 1)
        s += _ref_phase(c, cosΘ)
    end
    return s / n_phi
end

function _ref_effective_azimuthal_phase(contributors, τ_scat_layer, iz, ispec,
                                        μa, μb, n_phi)
    τ_scat = τ_scat_layer[iz, ispec]
    τ_scat == 0 && return zero(τ_scat)
    weighted = zero(τ_scat)
    for c in contributors
        scat_weight = c.τ[iz, ispec] * _ref_ssa(c)
        weighted += scat_weight * _ref_azimuthal_average_phase(c, μa, μb, n_phi)
    end
    return weighted / τ_scat
end

function _ref_τ_integral(τ_top, τ_bot, τ_total, μ_first, μ_second)
    b = 1 / μ_first - 1 / μ_second
    f_top = exp(-τ_top / μ_first - (τ_total - τ_top) / μ_second)
    f_bot = exp(-τ_bot / μ_first - (τ_total - τ_bot) / μ_second)
    abs(b) < 1e-10 && return 0.5 * (f_top + f_bot) * (τ_bot - τ_top)
    return (f_top - f_bot) / b
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

    FT = typeof(config.geometry.μ₀)
    path1 = zeros(FT, n_geom, 1, n_spec)
    path2 = zeros(FT, n_geom, 1, n_spec)
    path3 = zeros(FT, n_geom, 1, n_spec)
    path4 = zeros(FT, n_geom, 1, n_spec)

    if _want_ref_path1(paths)
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

    if _want_ref_path2(paths)
        for iv in 1:n_geom, s in 1:n_spec
            μ₀ = config.geometry.μ₀
            μv = config.geometry.μv[iv]
            τ = τ_cum[end, s]
            path2[iv, 1, s] =
                (μ₀ * I0[s] * albedo[s] / π) * exp(-τ / μ₀) * exp(-τ / μv)
        end
    end

    if _want_ref_path3(paths) || _want_ref_path4(paths)
        μ_nodes, μ_weights = _ref_gauss_legendre_01(config.inner_nquad, FT)
        for iv in 1:n_geom, s in 1:n_spec
            μ₀ = config.geometry.μ₀
            μv = config.geometry.μv[iv]
            τ_total = τ_cum[end, s]

            if _want_ref_path3(paths)
                F_surface = zero(FT)
                for iz in 1:n_layers
                    τ_scat = τ_scat_layer[iz, s]
                    τ_total_layer[iz, s] == 0 && continue
                    τ_scat == 0 && continue
                    ϖ_eff = τ_scat / τ_total_layer[iz, s]
                    layer_sum = zero(FT)
                    for k in eachindex(μ_nodes)
                        μd = μ_nodes[k]
                        P̄ = _ref_effective_azimuthal_phase(
                            contributors, τ_scat_layer, iz, s, μd, μ₀,
                            config.azimuth_nquad)
                        τ_int = _ref_τ_integral(τ_cum[iz, s], τ_cum[iz + 1, s],
                                                τ_total, μ₀, μd)
                        layer_sum += μ_weights[k] * P̄ * τ_int * FT(0.5)
                    end
                    F_surface += ϖ_eff * I0[s] * layer_sum
                end
                path3[iv, 1, s] =
                    (albedo[s] / π) * F_surface * exp(-τ_total / μv)
            end

            if _want_ref_path4(paths)
                L_surface = (albedo[s] / π) * μ₀ * I0[s] * exp(-τ_total / μ₀)
                for iz in 1:n_layers
                    τ_scat = τ_scat_layer[iz, s]
                    τ_total_layer[iz, s] == 0 && continue
                    τ_scat == 0 && continue
                    ϖ_eff = τ_scat / τ_total_layer[iz, s]
                    layer_sum = zero(FT)
                    for k in eachindex(μ_nodes)
                        μu = μ_nodes[k]
                        P̄ = _ref_effective_azimuthal_phase(
                            contributors, τ_scat_layer, iz, s, μu, μv,
                            config.azimuth_nquad)
                        τ_int = _ref_τ_integral(τ_cum[iz, s], τ_cum[iz + 1, s],
                                                τ_total, μv, μu)
                        layer_sum += μ_weights[k] * P̄ * τ_int * FT(0.5)
                    end
                    path4[iv, 1, s] += L_surface * ϖ_eff * layer_sum / μv
                end
            end
        end
    end

    (; path1, path2, path3, path4, total=path1 .+ path2 .+ path3 .+ path4)
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

function _benchmark_config(::Type{FT}; n_geom::Int = 8, n_spec::Int = 8,
                           n_layers::Int = 4) where {FT}
    geometry = SSGeometry(
        μ₀ = FT(0.82),
        μv = collect(range(FT(0.3), FT(0.85); length=n_geom)),
        Δϕ = collect(range(FT(0.1), FT(2.4); length=n_geom)))
    surface = LambertianSSSurface(albedo=fill(FT(0.22), n_spec))
    rayleigh = RayleighSSContributor(τ=fill(FT(0.03), n_layers, n_spec))
    aerosol = HGAerosolSSContributor(g=FT(0.7), ϖ=FT(0.9),
                                     τ=fill(FT(0.01), n_layers, n_spec))
    absorption = AbsorptionSSContributor(τ=fill(FT(0.02), n_layers, n_spec))
    return ExactSSConfig(geometry=geometry, surface=surface,
                         contributors=(rayleigh, aerosol, absorption),
                         I0=fill(one(FT), n_spec), inner_nquad=8,
                         azimuth_nquad=16)
end

function _average_elapsed(f, n::Int)
    total = 0.0
    for _ in 1:n
        GC.gc(false)
        total += @elapsed f()
    end
    return total / n
end

@testset "StandaloneSS Phase 1" begin
    @testset "Top-level exports" begin
        @test run_exact_ss === vSmartMOM.run_exact_ss
        @test run_exact_ss_with_jacobians === vSmartMOM.run_exact_ss_with_jacobians
        @test StandaloneSS === vSmartMOM.StandaloneSS
        @test CoxMunkSSSurface === vSmartMOM.CoxMunkSSSurface
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
        @test result.quadrature_info.inner_quadrature == config.inner_nquad
    end

    @testset "Cox-Munk path 2" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.78), μv=FT[0.42, 0.71],
                              Δϕ=FT[0.15, 1.2])
        n_water = Complex{FT}(FT(1.34), FT(1e-8))
        surface = CoxMunkSSSurface(wind_speed=FT(5.0), n_water=n_water,
                                   include_whitecaps=false, shadowing=true)
        absorption = AbsorptionSSContributor(τ=reshape(FT[0.02, 0.06], 2, 1))
        config = ExactSSConfig(geometry=geometry, surface=surface,
                               contributors=(absorption,), I0=FT[1.0])

        result = run_exact_ss(config; paths=:path2)
        expected = zeros(FT, 2, 1, 1)
        core_surface = vSmartMOM.CoreRT.CoxMunkSurface{FT}(
            wind_speed=surface.wind_speed,
            n_water=n_water,
            whitecap_albedo=FT(0.22),
            include_whitecaps=false,
            shadowing=true)
        τ = sum(absorption.τ[:, 1])
        for iv in eachindex(geometry.μv)
            M = vSmartMOM.CoreRT.coxmunk_brdf_mueller(
                core_surface, 1, geometry.μv[iv], geometry.μ₀,
                geometry.Δϕ[iv]; n_water=n_water)
            expected[iv, 1, 1] = geometry.μ₀ * M[1, 1] *
                                 exp(-τ / geometry.μ₀) *
                                 exp(-τ / geometry.μv[iv])
        end

        @test result.path1 == zero(result.path1)
        @test result.path2 ≈ expected rtol=1e-12 atol=1e-14
        @test result.total ≈ expected rtol=1e-12 atol=1e-14
        @test_throws ArgumentError run_exact_ss(config; paths=:all)
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

    @testset "Lambertian paths 3 and 4" begin
        config = _mixed_config(Float64)
        result = run_exact_ss(config; paths=:all)
        expected = _ref_outputs(config; paths=:all)
        @test result.path1 ≈ expected.path1 rtol=1e-12 atol=1e-14
        @test result.path2 ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test result.path3 ≈ expected.path3 rtol=1e-12 atol=1e-14
        @test result.path4 ≈ expected.path4 rtol=1e-12 atol=1e-14
        @test result.total ≈ expected.total rtol=1e-12 atol=1e-14
        @test all(result.path3 .> 0)
        @test all(result.path4 .> 0)
        @test result.quadrature_info.inner_quadrature == config.inner_nquad
        @test result.quadrature_info.azimuth_quadrature == config.azimuth_nquad
    end

    @testset "Path selection" begin
        config = _mixed_config(Float64)
        only1 = run_exact_ss(config; paths=:path1)
        only2 = run_exact_ss(config; paths=:path2)
        only3 = run_exact_ss(config; paths=:path3)
        only4 = run_exact_ss(config; paths=:path4)
        both = run_exact_ss(config; paths=:paths_1_2)
        all_paths = run_exact_ss(config; paths=:all)
        @test only1.path2 == zero(only1.path2)
        @test only2.path1 == zero(only2.path1)
        @test only3.path1 == zero(only3.path1)
        @test only3.path2 == zero(only3.path2)
        @test only3.path4 == zero(only3.path4)
        @test only4.path1 == zero(only4.path1)
        @test only4.path2 == zero(only4.path2)
        @test only4.path3 == zero(only4.path3)
        @test both.total ≈ only1.total .+ only2.total rtol=1e-12
        @test all_paths.total ≈
              only1.total .+ only2.total .+ only3.total .+ only4.total rtol=1e-12
        @test all_paths.path3 ≈ only3.path3 rtol=1e-12
        @test all_paths.path4 ≈ only4.path4 rtol=1e-12
    end

    @testset "Type inference and path parity" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.79), μv=FT[0.46], Δϕ=FT[0.4])
        surface = LambertianSSSurface(albedo=FT(0.27))
        rayleigh = RayleighSSContributor(τ=reshape(FT[0.05, 0.08], 2, 1))
        config = ExactSSConfig(geometry=geometry, surface=surface,
                               contributors=(rayleigh,), I0=FT[1.0],
                               inner_nquad=4, azimuth_nquad=8)

        optics = @inferred vSmartMOM.StandaloneSS._precompute_optics(config)
        ρ = @inferred vSmartMOM.StandaloneSS._precompute_surface_brdf(
            config.surface, config.geometry, 1, FT)
        μ_nodes, _ = vSmartMOM.StandaloneSS._gauss_legendre_01(config.inner_nquad, FT)
        reference_μ₀ = fill(config.geometry.μ₀, length(config.geometry.μv))
        P̄3, P̄4 = @inferred vSmartMOM.StandaloneSS._precompute_azimuthal_phase_pair(
            config, optics.τ_scat_layer, μ_nodes, reference_μ₀, config.geometry.μv)

        @test eltype(optics.τ_cum) == FT
        @test eltype(ρ) == FT
        @test size(P̄3) == size(P̄4) == (1, 2, 1, config.inner_nquad)

        all_paths = @inferred run_exact_ss(config; paths=:all)
        split3 = run_exact_ss(config; paths=:path3)
        split4 = run_exact_ss(config; paths=:path4)
        @test all_paths.path3 ≈ split3.path3 rtol=1e-12 atol=1e-14
        @test all_paths.path4 ≈ split4.path4 rtol=1e-12 atol=1e-14

        dualτ = ForwardDiff.Dual{Nothing}(FT(0.05), one(FT))
        dual_config = ExactSSConfig(
            geometry=geometry,
            surface=surface,
            contributors=(RayleighSSContributor(
                τ=reshape([dualτ, zero(dualτ) + FT(0.08)], 2, 1)),),
            I0=FT[1.0])
        dual_result = @inferred run_exact_ss(dual_config; paths=:path1)
        @test eltype(dual_result.total) <: ForwardDiff.Dual
        @test isfinite(ForwardDiff.value(dual_result.total[1, 1, 1]))
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
        black_all = run_exact_ss(black_config; paths=:all)
        @test black_all.path3 == zero(black_all.path3)
        @test black_all.path4 == zero(black_all.path4)
        @test black_result.path1[1, 1, 1] > 0

        absorption = AbsorptionSSContributor(τ=reshape(FT[0.2, 0.3], 2, 1))
        surface = LambertianSSSurface(albedo=FT(0.4))
        abs_config = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(absorption,), I0=FT[1.0])
        abs_result = run_exact_ss(abs_config)
        abs_expected = _ref_outputs(abs_config)
        @test abs_result.path1 == zero(abs_result.path1)
        abs_all = run_exact_ss(abs_config; paths=:all)
        @test abs_all.path3 == zero(abs_all.path3)
        @test abs_all.path4 == zero(abs_all.path4)
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

    @testset "ForwardDiff compatibility" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.76), μv=FT[0.48], Δϕ=FT[0.35])

        fτ(x) = begin
            surface = LambertianSSSurface(albedo=FT(0.31))
            absorption = AbsorptionSSContributor(τ=reshape([x[1], FT(0.04)], 2, 1))
            config = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(absorption,), I0=FT[1.0])
            vec(run_exact_ss(config; paths=:path2).total)
        end
        Jτ = ForwardDiff.jacobian(fτ, FT[0.18])
        baseτ = only(fτ(FT[0.18]))
        atten_slope = inv(geometry.μ₀) + inv(only(geometry.μv))
        @test Jτ[1, 1] ≈ -baseτ * atten_slope rtol=1e-12 atol=1e-14

        fw(u) = begin
            n_water = Complex{FT}(FT(1.34), FT(1e-8))
            surface = CoxMunkSSSurface(wind_speed=u[1], n_water=n_water,
                                       include_whitecaps=false, shadowing=true)
            absorption = AbsorptionSSContributor(τ=reshape(FT[0.07, 0.03], 2, 1))
            config = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(absorption,), I0=FT[1.0])
            vec(run_exact_ss(config; paths=:path2).total)
        end
        Jw = ForwardDiff.jacobian(fw, FT[4.0])
        h = FT(1e-5)
        fd = (only(fw(FT[4.0 + h])) - only(fw(FT[4.0 - h]))) / (2h)
        @test Jw[1, 1] ≈ fd rtol=2e-5 atol=1e-8
    end

    @testset "Handcoded f2 Jacobians for paths 1 and 2" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.81), μv=FT[0.47], Δϕ=FT[0.25])

        path1_config_from_τ(x) = ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=zero(FT)),
            contributors=(RayleighSSContributor(τ=reshape(x, 2, 1)),),
            I0=FT[1.0])
        τ_rayleigh = FT[0.07, 0.04]
        f1(x) = vec(run_exact_ss(path1_config_from_τ(x); paths=:path1).total)
        J1 = ForwardDiff.jacobian(f1, τ_rayleigh)
        jac1 = run_exact_ss_with_jacobians(path1_config_from_τ(τ_rayleigh);
                                           paths=:path1).jacobians
        @test vec(jac1.path1.τ_layer[1, 1, 1, :]) ≈ vec(J1) rtol=1e-12 atol=1e-14
        @test all(jac1.path1.ϖ_eff .>= 0)
        @test all(jac1.path1.P_eff .>= 0)

        path2_config_from_τ(x) = ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=FT(0.33)),
            contributors=(AbsorptionSSContributor(τ=reshape(x, 2, 1)),),
            I0=FT[1.0])
        τ_abs = FT[0.09, 0.03]
        f2(x) = vec(run_exact_ss(path2_config_from_τ(x); paths=:path2).total)
        J2 = ForwardDiff.jacobian(f2, τ_abs)
        jac2 = run_exact_ss_with_jacobians(path2_config_from_τ(τ_abs);
                                           paths=:path2).jacobians
        @test vec(jac2.path2.τ_layer[1, 1, 1, :]) ≈ vec(J2) rtol=1e-12 atol=1e-14

        falbedo(a) = only(run_exact_ss(ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=a),
            contributors=(AbsorptionSSContributor(τ=reshape(τ_abs, 2, 1)),),
            I0=FT[1.0]); paths=:path2).total)
        d_albedo = ForwardDiff.derivative(falbedo, FT(0.33))
        @test jac2.path2.surface_brdf[1, 1, 1] / π ≈ d_albedo rtol=1e-12 atol=1e-14
        @test jac2.total.τ_layer ≈ jac2.path2.τ_layer
    end

    @testset "CUDA equivalency and speed smoke" begin
        if !_standalone_cuda_functional()
            @test_skip "CUDA not available or not functional; skipping StandaloneSS GPU smoke"
        else
            FT = Float64
            cpu_config = _benchmark_config(FT; n_geom=32, n_spec=32, n_layers=6)
            gpu_config = ExactSSConfig(
                geometry=cpu_config.geometry,
                surface=cpu_config.surface,
                contributors=cpu_config.contributors,
                I0=cpu_config.I0,
                inner_nquad=cpu_config.inner_nquad,
                azimuth_nquad=cpu_config.azimuth_nquad,
                architecture=vSmartMOM.GPU())

            cpu_result = run_exact_ss(cpu_config; paths=:all)
            gpu_result = run_exact_ss(gpu_config; paths=:all)
            CUDA.synchronize()

            @test gpu_result.quadrature_info.kernel_backend == :CUDA
            @test gpu_result.total isa CUDA.CuArray
            for field in (:total, :path1, :path2, :path3, :path4)
                @test Array(getfield(gpu_result, field)) ≈
                      getfield(cpu_result, field) rtol=5e-12 atol=5e-14
            end

            cpu_time = _average_elapsed(3) do
                run_exact_ss(cpu_config; paths=:all)
                nothing
            end
            gpu_time = _average_elapsed(3) do
                run_exact_ss(gpu_config; paths=:all)
                CUDA.synchronize()
                nothing
            end
            @info "StandaloneSS CPU/GPU timing smoke" cpu_time gpu_time speedup=cpu_time / gpu_time
            @test isfinite(cpu_time) && cpu_time > 0
            @test isfinite(gpu_time) && gpu_time > 0
        end
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
        @test_throws ArgumentError run_exact_ss(good; paths=:not_a_path)
        @test_throws ArgumentError run_exact_ss(ExactSSConfig(
            geometry=good.geometry, surface=good.surface,
            contributors=good.contributors, I0=good.I0, n_stokes=3))
        @test_throws ArgumentError run_exact_ss(ExactSSConfig(
            geometry=good.geometry, surface=good.surface,
            contributors=good.contributors, I0=good.I0, inner_nquad=0))
        @test_throws ArgumentError run_exact_ss(ExactSSConfig(
            geometry=good.geometry, surface=good.surface,
            contributors=good.contributors, I0=good.I0, azimuth_nquad=0))
        @test_throws ArgumentError run_exact_ss(ExactSSConfig(
            geometry=good.geometry,
            surface=CoxMunkSSSurface(wind_speed=FT(4.0)),
            contributors=good.contributors, I0=good.I0); paths=:path3)
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
        coxmunk = CoxMunkSSSurface(wind_speed=4.0)
        @test determine_required_nbrdf(coxmunk) ==
              determine_required_nbrdf_coxmunk(coxmunk.wind_speed)
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
