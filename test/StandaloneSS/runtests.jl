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

function _ref_phase(c::RayleighSSContributor, cosΘ)
    dpl_p = (1 - c.depol) / (1 + c.depol / 2)
    P₂ = (3 * cosΘ^2 - 1) / 2
    return 1 + 0.5 * dpl_p * P₂
end
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

function _ref_rayleigh_phase_first_column(μ₀, μv, Δϕ, depol, n_stokes)
    cosΘ = _ref_cos_scatter(μ₀, μv, Δϕ)
    dpl_p = (1 - depol) / (1 + depol / 2)
    f11 = 1 + 0.5 * dpl_p * ((3 * cosΘ^2 - 1) / 2)
    f12 = 0.75 * dpl_p * max(0, 1 - cosΘ^2)
    p = zeros(typeof(f11), n_stokes)
    p[1] = f11
    n_stokes == 1 && return p

    sinΘ² = max(0, 1 - cosΘ^2)
    if sinΘ² <= eps(typeof(f11))
        cos2χ = one(f11)
        sin2χ = zero(f11)
    else
        sinΘ = sqrt(sinΘ²)
        sin₀ = sqrt(max(0, 1 - μ₀^2))
        sinᵥ = sqrt(max(0, 1 - μv^2))
        cosχ = (μ₀ * sinᵥ + μv * sin₀ * cos(Δϕ)) / sinΘ
        sinχ = sin₀ * sin(Δϕ) / sinΘ
        cos2χ = cosχ^2 - sinχ^2
        sin2χ = 2 * sinχ * cosχ
    end
    p[2] = f12 * cos2χ
    n_stokes >= 3 && (p[3] = f12 * sin2χ)
    return p
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
        @test chain_rule_combine_dP === vSmartMOM.chain_rule_combine_dP
        @test chain_rule_combine_dτ === vSmartMOM.chain_rule_combine_dτ
        @test chain_rule_combine_dϖ === vSmartMOM.chain_rule_combine_dϖ
        @test chain_rule_combine_surface_brdf ===
              vSmartMOM.chain_rule_combine_surface_brdf
        @test exact_ss_config_from_model === vSmartMOM.exact_ss_config_from_model
        @test SSMeasurementSelector === vSmartMOM.SSMeasurementSelector
        @test GreekCoefsSSContributor === vSmartMOM.GreekCoefsSSContributor
        @test selected_measurements === vSmartMOM.selected_measurements
        @test selected_measurement_jacobian ===
              vSmartMOM.selected_measurement_jacobian
        @test StandaloneSS === vSmartMOM.StandaloneSS
        @test CoxMunkSSSurface === vSmartMOM.CoxMunkSSSurface
        @test vSmartMOM.HenyeyGreensteinPhaseFunction ===
              vSmartMOM.Scattering.HenyeyGreensteinPhaseFunction
        @test vSmartMOM.SyntheticPolarizedHenyeyGreensteinPhaseFunction ===
              vSmartMOM.Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction
    end

    @testset "RTModel adapter" begin
        cfg = Dict(
            "radiative_transfer" => Dict(
                "spec_bands" => ["[13000.0]"],
                "surface" => ["LambertianSurfaceScalar(0.2)"],
                "quadrature_type" => "GaussQuadHemisphere()",
                "polarization_type" => "Stokes_IQU()",
                "max_m" => 3,
                "Δ_angle" => 0.0,
                "l_trunc" => 8,
                "truncation" => "NoTruncation(l_max=8)",
                "depol" => 0.0,
                "float_type" => "Float64",
                "architecture" => "CPU()",
            ),
            "geometry" => Dict(
                "sza" => 35.0,
                "vza" => [20.0, 50.0],
                "vaz" => [0.0, 90.0],
                "obs_alt" => 1000.0,
            ),
            "atmospheric_profile" => Dict(
                "T" => [280.0],
                "p" => [100.0, 1000.0],
                "profile_reduction" => -1,
            ),
            "scattering" => Dict(
                "aerosols" => [Dict(
                    "τ_ref" => 0.02,
                    "phase_function" => "SyntheticPolarizedHenyeyGreensteinPhaseFunction(g=0.45, polarization_fraction=0.6)",
                    "ssa" => 0.92,
                    "p₀" => 550.0,
                    "σp" => 200.0,
                )],
                "r_max" => 1.0,
                "nquad_radius" => 16,
                "λ_ref" => 0.55,
                "decomp_type" => "NAI2()",
            ),
        )

        params = parameters_from_dict(cfg)
        model = model_from_parameters(params)
        config = exact_ss_config_from_model(model)
        @test config.geometry.μ₀ ≈ cosd(params.sza)
        @test config.surface isa LambertianSSSurface
        @test length(config.contributors) == 2
        @test config.contributors[1] isa GreekCoefsSSContributor
        @test config.contributors[2] isa GreekCoefsSSContributor
        @test config.polarization_type isa vSmartMOM.Scattering.Stokes_IQU{Float64}

        result = run_exact_ss(config; paths=:paths_1_2)
        model_result = run_exact_ss(model; paths=:paths_1_2)
        @test model_result.total ≈ result.total rtol=1e-12 atol=1e-14
        @test size(result.total) == (2, 3, 1)
        @test all(isfinite.(result.total))
        @test any(abs.(result.path1[:, 2:3, :]) .> 0)

        scalar_cfg = deepcopy(cfg)
        scalar_cfg["radiative_transfer"]["polarization_type"] = "Stokes_I()"
        scalar_model = model_from_parameters(parameters_from_dict(scalar_cfg))
        scalar_all = run_exact_ss(scalar_model; paths=:all)
        @test size(scalar_all.total) == (2, 1, 1)
        @test all(isfinite.(scalar_all.total))
        @test any(scalar_all.path3 .!= 0)
        @test any(scalar_all.path4 .!= 0)

        cfg_trunc = deepcopy(cfg)
        delete!(cfg_trunc["radiative_transfer"], "truncation")
        truncated_model = model_from_parameters(parameters_from_dict(cfg_trunc))
        @test_throws ArgumentError exact_ss_config_from_model(truncated_model)

        cox_cfg = deepcopy(cfg)
        cox_cfg["radiative_transfer"]["surface"] = ["CoxMunkSurface(wind_speed=5.0)"]
        delete!(cox_cfg, "scattering")
        cox_model = model_from_parameters(parameters_from_dict(cox_cfg))
        cox_config = exact_ss_config_from_model(cox_model)
        @test cox_config.surface isa CoxMunkSSSurface
        cox_result = run_exact_ss(cox_config; paths=:path2)
        cox_model_result = run_exact_ss(cox_model; paths=:path2)
        @test cox_model_result.path2 ≈ cox_result.path2 rtol=1e-12 atol=1e-14
        @test size(cox_result.path2) == (2, 3, 1)
        @test all(isfinite.(cox_result.path2))
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

        pol_config = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(rayleigh,), I0=FT[1.0, 0.7],
                                   polarization_type=vSmartMOM.Scattering.Stokes_I{FT}())
        pol_result = run_exact_ss(pol_config; paths=:paths_1_2)
        @test pol_result.total ≈ result.total rtol=1e-12 atol=1e-14
        @test pol_result.metadata.n_stokes == 1

        iq_config = ExactSSConfig(
            geometry=geometry,
            surface=surface,
            contributors=(rayleigh,),
            I0=FT[1.0, 0.7],
            polarization_type=vSmartMOM.Scattering.Stokes_IQ{FT}())
        iq_path1 = run_exact_ss(iq_config; paths=:path1)

        vector_config = ExactSSConfig(
            geometry=geometry,
            surface=surface,
            contributors=(rayleigh,),
            I0=FT[1.0, 0.7],
            polarization_type=vSmartMOM.Scattering.Stokes_IQU{FT}())
        vector_path1 = run_exact_ss(vector_config; paths=:path1)
        expected_path1 = zeros(FT, 2, 3, 2)
        for iv in eachindex(geometry.μv), ispec in 1:2
            phase_col = _ref_rayleigh_phase_first_column(
                geometry.μ₀, geometry.μv[iv], geometry.Δϕ[iv],
                rayleigh.depol, 3)
            expected_path1[iv, :, ispec] .=
                expected.path1[iv, 1, ispec] .* phase_col ./ phase_col[1]
        end
        @test size(vector_path1.path1) == (2, 3, 2)
        @test vector_path1.path1 ≈ expected_path1 rtol=1e-12 atol=1e-14
        @test vector_path1.total ≈ expected_path1 rtol=1e-12 atol=1e-14
        @test vector_path1.path2 == zero(vector_path1.path2)
        @test size(iq_path1.path1) == (2, 2, 2)
        @test iq_path1.path1 ≈ expected_path1[:, 1:2, :] rtol=1e-12 atol=1e-14

        vector_path2 = run_exact_ss(vector_config; paths=:path2)
        @test size(vector_path2.path2) == (2, 3, 2)
        @test vector_path2.path2[:, 1:1, :] ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test vector_path2.path2[:, 2:3, :] == zero(vector_path2.path2[:, 2:3, :])
        iq_path2 = run_exact_ss(iq_config; paths=:path2)
        @test iq_path2.path2[:, 1:1, :] ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test iq_path2.path2[:, 2:2, :] == zero(iq_path2.path2[:, 2:2, :])
        vector_paths12 = run_exact_ss(vector_config; paths=:paths_1_2)
        expected_paths12 = copy(expected_path1)
        expected_paths12[:, 1:1, :] .+= expected.path2
        @test vector_paths12.path1 ≈ expected_path1 rtol=1e-12 atol=1e-14
        @test vector_paths12.path2[:, 1:1, :] ≈ expected.path2 rtol=1e-12 atol=1e-14
        @test vector_paths12.path2[:, 2:3, :] == zero(vector_paths12.path2[:, 2:3, :])
        @test vector_paths12.total ≈ expected_paths12 rtol=1e-12 atol=1e-14
        @test_throws ArgumentError run_exact_ss_with_jacobians(
            vector_config; paths=:path1)
        @test_throws ArgumentError run_exact_ss_with_jacobians(
            vector_config; paths=:path2)
    end

    @testset "Greek-coefficient aerosol vector path 1" begin
        FT = Float64
        phase = vSmartMOM.Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction(
            g=FT(0.45), polarization_fraction=FT(0.6))
        greek = vSmartMOM.Scattering.greek_coefficients(phase; l_max=24,
                                                        nquad=80)
        optics = vSmartMOM.Scattering.analytic_aerosol_optics(
            phase; single_scattering_albedo=FT(0.92),
            extinction_cross_section=FT(1.7), l_max=24, nquad=80)
        @test length(greek.β) == 24
        @test optics.greek_coefs ≈ greek
        @test optics.ω̃ == FT(0.92)
        @test optics.k == FT(1.7)

        geometry = SSGeometry(μ₀=FT(0.82),
                              μv=FT[0.51, 0.69],
                              Δϕ=FT[0.4, 1.1])
        aerosol = GreekCoefsSSContributor(greek_coefs=greek, ϖ=FT(0.92),
                                          τ=reshape(FT[0.03, 0.04], 2, 1))
        absorption = AbsorptionSSContributor(
            τ=reshape(FT[0.01, 0.02], 2, 1))
        scalar_config = ExactSSConfig(
            geometry=geometry, surface=LambertianSSSurface(albedo=FT(0.0)),
            contributors=(aerosol, absorption), I0=FT(1.0))
        vector_config = ExactSSConfig(
            geometry=geometry, surface=scalar_config.surface,
            contributors=scalar_config.contributors, I0=scalar_config.I0,
            polarization_type=vSmartMOM.Scattering.Stokes_IQU{FT}())

        scalar = run_exact_ss(scalar_config; paths=:path1)
        vector = run_exact_ss(vector_config; paths=:path1)
        expected = zeros(FT, 2, 3, 1)
        for iv in eachindex(geometry.μv)
            phase_col = vSmartMOM.Scattering.phase_matrix_first_column(
                greek, geometry.μ₀, geometry.μv[iv], geometry.Δϕ[iv],
                Val(3))
            expected[iv, :, 1] .= scalar.path1[iv, 1, 1] .*
                                  collect(phase_col) ./ phase_col[1]
        end
        @test vector.path1 ≈ expected rtol=1e-10 atol=1e-13
        @test vector.path1[:, 1:1, :] ≈ scalar.path1 rtol=1e-10 atol=1e-13
        @test any(abs.(vector.path1[:, 2:3, :]) .> 0)

        vector4_config = ExactSSConfig(
            geometry=geometry, surface=scalar_config.surface,
            contributors=scalar_config.contributors, I0=scalar_config.I0,
            polarization_type=vSmartMOM.Scattering.Stokes_IQUV{FT}())
        vector4 = run_exact_ss(vector4_config; paths=:path1)
        expected4 = zeros(FT, 2, 4, 1)
        for iv in eachindex(geometry.μv)
            phase_col = vSmartMOM.Scattering.phase_matrix_first_column(
                greek, geometry.μ₀, geometry.μv[iv], geometry.Δϕ[iv],
                Val(4))
            expected4[iv, :, 1] .= scalar.path1[iv, 1, 1] .*
                                   collect(phase_col) ./ phase_col[1]
        end
        @test vector4.path1 ≈ expected4 rtol=1e-10 atol=1e-13
        @test vector4.path1[:, 1:1, :] ≈ scalar.path1 rtol=1e-10 atol=1e-13
        @test vector4.path1[:, 4:4, :] == zero(vector4.path1[:, 4:4, :])
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

        vector_config = ExactSSConfig(
            geometry=geometry, surface=surface,
            contributors=(absorption,), I0=FT[1.0],
            polarization_type=vSmartMOM.Scattering.Stokes_IQUV{FT}())
        vector_result = run_exact_ss(vector_config; paths=:path2)
        vector_expected = zeros(FT, 2, 4, 1)
        for iv in eachindex(geometry.μv)
            M = vSmartMOM.CoreRT.coxmunk_brdf_mueller(
                core_surface, 4, geometry.μv[iv], geometry.μ₀,
                geometry.Δϕ[iv]; n_water=n_water)
            vector_expected[iv, :, 1] .= geometry.μ₀ .* M[:, 1] .*
                                         exp(-τ / geometry.μ₀) .*
                                         exp(-τ / geometry.μv[iv])
        end
        @test vector_result.path2 ≈ vector_expected rtol=1e-12 atol=1e-14

        iq_config = ExactSSConfig(
            geometry=geometry, surface=surface,
            contributors=(absorption,), I0=FT[1.0],
            polarization_type=vSmartMOM.Scattering.Stokes_IQ{FT}())
        iq_result = run_exact_ss(iq_config; paths=:path2)
        iq_expected = zeros(FT, 2, 2, 1)
        for iv in eachindex(geometry.μv)
            M = vSmartMOM.CoreRT.coxmunk_brdf_mueller(
                core_surface, 2, geometry.μv[iv], geometry.μ₀,
                geometry.Δϕ[iv]; n_water=n_water)
            iq_expected[iv, :, 1] .= geometry.μ₀ .* M[:, 1] .*
                                     exp(-τ / geometry.μ₀) .*
                                     exp(-τ / geometry.μv[iv])
        end
        @test iq_result.path2 ≈ iq_expected rtol=1e-12 atol=1e-14
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

    @testset "Rayleigh azimuth average is analytic" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.82),
                              μv=FT[0.43, 0.77],
                              Δϕ=FT[0.1, 1.3])
        surface = LambertianSSSurface(albedo=FT(0.24))
        rayleigh = RayleighSSContributor(τ=reshape(FT[0.04, 0.07], 2, 1))
        base = ExactSSConfig(geometry=geometry, surface=surface,
                             contributors=(rayleigh,), I0=FT[1.0],
                             inner_nquad=4, azimuth_nquad=1)
        refined = ExactSSConfig(geometry=geometry, surface=surface,
                                contributors=(rayleigh,), I0=FT[1.0],
                                inner_nquad=4, azimuth_nquad=32)
        coarse_result = run_exact_ss(base; paths=:all)
        refined_result = run_exact_ss(refined; paths=:all)
        @test coarse_result.path3 ≈ refined_result.path3 rtol=1e-12 atol=1e-14
        @test coarse_result.path4 ≈ refined_result.path4 rtol=1e-12 atol=1e-14

        depol = FT(0.028)
        depol_rayleigh = RayleighSSContributor(
            τ=reshape(FT[0.04, 0.07], 2, 1), depol=depol)
        cosΘ = FT(0.31)
        dpl_p = (one(FT) - depol) / (one(FT) + depol / FT(2))
        expected_phase = one(FT) + FT(0.5) * dpl_p *
                         ((FT(3) * cosΘ^2 - one(FT)) / FT(2))
        @test exact_phase_function(depol_rayleigh, cosΘ) ≈ expected_phase

        depol_base = ExactSSConfig(geometry=geometry, surface=surface,
                                   contributors=(depol_rayleigh,),
                                   I0=FT[1.0], inner_nquad=4,
                                   azimuth_nquad=1)
        depol_refined = ExactSSConfig(geometry=geometry, surface=surface,
                                      contributors=(depol_rayleigh,),
                                      I0=FT[1.0], inner_nquad=4,
                                      azimuth_nquad=32)
        depol_coarse = run_exact_ss(depol_base; paths=:all)
        depol_fine = run_exact_ss(depol_refined; paths=:all)
        @test depol_coarse.path3 ≈ depol_fine.path3 rtol=1e-12 atol=1e-14
        @test depol_coarse.path4 ≈ depol_fine.path4 rtol=1e-12 atol=1e-14
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

    @testset "Measurement selector" begin
        config = _mixed_config(Float64)
        result = run_exact_ss(config; paths=:all)

        default_selector = SSMeasurementSelector()
        @test selected_measurements(result, default_selector) ==
              vec(result.total[:, :, :])
        @test selected_measurements(result) ==
              selected_measurements(result, default_selector)

        one_geometry = SSMeasurementSelector(geometry_indices=2)
        @test selected_measurements(result, one_geometry) ==
              vec(result.total[2:2, 1:1, :])

        multiple_geometries = SSMeasurementSelector(geometry_indices=[2, 1])
        @test selected_measurements(result, multiple_geometries) ==
              vec(result.total[[2, 1], 1:1, :])

        spectral_subset = SSMeasurementSelector(spectral_indices=2)
        @test selected_measurements(result, spectral_subset) ==
              vec(result.total[:, 1:1, 2:2])

        path_selector = SSMeasurementSelector(paths=(:path1, :path2),
                                              geometry_indices=1:2,
                                              spectral_indices=1)
        @test selected_measurements(result, path_selector) ==
              vcat(vec(result.path1[:, 1:1, 1:1]),
                   vec(result.path2[:, 1:1, 1:1]))

        J4 = reshape(collect(1.0:18.0), 2, 1, 3, 3)
        jac_selector = SSMeasurementSelector(geometry_indices=2,
                                             spectral_indices=[1, 3])
        @test selected_measurement_jacobian(J4, jac_selector) ==
              reshape(J4[2:2, 1:1, [1, 3], :], :, 3)

        polarized_result = (; total = reshape(collect(1.0:12.0), 2, 3, 2))
        @test selected_measurements(polarized_result) ==
              vec(polarized_result.total)
        stokes_i_selector = SSMeasurementSelector(stokes_indices=1)
        @test selected_measurements(polarized_result, stokes_i_selector) ==
              vec(polarized_result.total[:, 1:1, :])

        J_polarized = reshape(collect(1.0:48.0), 2, 3, 2, 4)
        @test selected_measurement_jacobian(J_polarized) ==
              reshape(J_polarized, :, 4)
        stokes_subset = SSMeasurementSelector(stokes_indices=[3, 1])
        @test selected_measurement_jacobian(J_polarized, stokes_subset) ==
              reshape(J_polarized[:, [3, 1], :, :], :, 4)

        @test_throws ArgumentError SSMeasurementSelector(paths=())
        @test_throws ArgumentError SSMeasurementSelector(paths=(:bogus,))
        @test_throws ArgumentError selected_measurement_jacobian(
            J4, SSMeasurementSelector(paths=(:total, :path1)))
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
        path1_with_jac = run_exact_ss_with_jacobians(
            path1_config_from_τ(τ_rayleigh); paths=:path1)
        jac1 = path1_with_jac.jacobians
        @test path1_with_jac.measurements ==
              selected_measurements(path1_with_jac, SSMeasurementSelector())
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
        path2_selector = SSMeasurementSelector(paths=:path2)
        path2_with_jac = run_exact_ss_with_jacobians(
            path2_config_from_τ(τ_abs); paths=:path2, selector=path2_selector)
        jac2 = path2_with_jac.jacobians
        @test path2_with_jac.measurement_selector == path2_selector
        @test path2_with_jac.measurements ==
              selected_measurements(path2_with_jac, path2_selector)
        @test vec(jac2.path2.τ_layer[1, 1, 1, :]) ≈ vec(J2) rtol=1e-12 atol=1e-14

        falbedo(a) = only(run_exact_ss(ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=a),
            contributors=(AbsorptionSSContributor(τ=reshape(τ_abs, 2, 1)),),
            I0=FT[1.0]); paths=:path2).total)
        d_albedo = ForwardDiff.derivative(falbedo, FT(0.33))
        @test jac2.path2.surface_brdf[1, 1, 1] / π ≈ d_albedo rtol=1e-12 atol=1e-14
        @test jac2.total.τ_layer ≈ jac2.path2.τ_layer
        @test_throws ArgumentError run_exact_ss_with_jacobians(
            path1_config_from_τ(τ_rayleigh); paths=:path1,
            selector=SSMeasurementSelector(paths=:path2))
    end

    @testset "f2 τ chain-rule contraction" begin
        FT = Float64
        geometry = SSGeometry(μ₀=FT(0.81), μv=FT[0.47], Δϕ=FT[0.25])

        path1_baseτ = FT[0.07, 0.04]
        path1_config_from_p(p) = ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=zero(FT)),
            contributors=(RayleighSSContributor(
                τ=reshape(path1_baseτ .* p, 2, 1)),),
            I0=FT[1.0])
        p0 = ones(FT, 2)
        f1p(p) = vec(run_exact_ss(path1_config_from_p(p); paths=:path1).total)
        J1p = ForwardDiff.jacobian(f1p, p0)
        jac1 = run_exact_ss_with_jacobians(path1_config_from_p(p0);
                                           paths=:path1).jacobians
        dτ1_dp = zeros(FT, 2, 1, 2)
        dτ1_dp[1, 1, 1] = path1_baseτ[1]
        dτ1_dp[2, 1, 2] = path1_baseτ[2]
        J1_chain = @inferred chain_rule_combine_dτ(jac1.path1.τ_layer, dτ1_dp)
        @test size(J1_chain) == (1, 1, 1, 2)
        @test vec(J1_chain) ≈ vec(J1p) rtol=1e-12 atol=1e-14
        path1_selector = SSMeasurementSelector(paths=:path1)
        J1_selected = @inferred chain_rule_combine_dτ(
            jac1.path1.τ_layer, dτ1_dp, path1_selector)
        @test size(J1_selected) == (1, 2)
        @test J1_selected ≈ reshape(J1_chain, :, 2) rtol=1e-12 atol=1e-14
        @test_throws ArgumentError chain_rule_combine_dτ(
            jac1.path1.τ_layer, dτ1_dp,
            SSMeasurementSelector(paths=(:path1, :path2)))

        path2_baseτ = FT[0.09, 0.03]
        path2_config_from_p(p) = ExactSSConfig(
            geometry=geometry,
            surface=LambertianSSSurface(albedo=FT(0.33)),
            contributors=(AbsorptionSSContributor(
                τ=reshape(path2_baseτ .* p, 2, 1)),),
            I0=FT[1.0])
        f2p(p) = vec(run_exact_ss(path2_config_from_p(p); paths=:path2).total)
        J2p = ForwardDiff.jacobian(f2p, p0)
        jac2 = run_exact_ss_with_jacobians(path2_config_from_p(p0);
                                           paths=:path2).jacobians
        dτ2_dp = zeros(FT, 2, 1, 2)
        dτ2_dp[1, 1, 1] = path2_baseτ[1]
        dτ2_dp[2, 1, 2] = path2_baseτ[2]
        J2_chain = @inferred chain_rule_combine_dτ(jac2.path2.τ_layer, dτ2_dp)
        @test vec(J2_chain) ≈ vec(J2p) rtol=1e-12 atol=1e-14

        @test_throws ArgumentError chain_rule_combine_dτ(
            zeros(FT, 1, 1, 1, 2), zeros(FT, 3, 1, 1))
        @test_throws ArgumentError chain_rule_combine_dτ(
            zeros(FT, 1, 1, 1, 2), zeros(FT, 2, 2, 1))
        @test_throws ArgumentError chain_rule_combine_dτ(
            zeros(FT, 1, 1, 1), zeros(FT, 1, 1))

        hg_geometry = SSGeometry(μ₀=FT(0.77),
                                 μv=FT[0.43, 0.69],
                                 Δϕ=FT[0.15, 1.0])
        τ_hg = FT[0.05 0.08; 0.03 0.04]
        g0 = FT(0.35)
        ϖ0 = FT(0.82)
        hg_contributor(g, ϖ) = begin
            T = promote_type(typeof(g), typeof(ϖ))
            HGAerosolSSContributor(g=convert(T, g), ϖ=convert(T, ϖ),
                                   τ=τ_hg)
        end
        hg_config(g, ϖ) = ExactSSConfig(
            geometry=hg_geometry,
            surface=LambertianSSSurface(albedo=zero(FT)),
            contributors=(hg_contributor(g, ϖ),),
            I0=FT[1.0, 0.9])
        hg_selector = SSMeasurementSelector(paths=:path1,
                                            geometry_indices=[2, 1],
                                            spectral_indices=1:2)
        hg_jac = run_exact_ss_with_jacobians(hg_config(g0, ϖ0);
                                             paths=:path1,
                                             selector=hg_selector).jacobians

        fϖ(p) = selected_measurements(
            run_exact_ss(hg_config(g0, p[1]); paths=:path1),
            hg_selector)
        Jϖ = ForwardDiff.jacobian(fϖ, FT[ϖ0])
        dϖ_dp = ones(FT, 2, 2, 1)
        Jϖ_chain = @inferred chain_rule_combine_dϖ(
            hg_jac.path1.ϖ_eff, dϖ_dp, hg_selector)
        @test Jϖ_chain ≈ Jϖ rtol=1e-12 atol=1e-14

        fg(p) = selected_measurements(
            run_exact_ss(hg_config(p[1], ϖ0); paths=:path1),
            hg_selector)
        Jg = ForwardDiff.jacobian(fg, FT[g0])
        dP_dp = zeros(FT, 2, 2, 2, 1)
        for iv in eachindex(hg_geometry.μv)
            cosΘ = _ref_cos_scatter(hg_geometry.μ₀, hg_geometry.μv[iv],
                                    hg_geometry.Δϕ[iv])
            dP_dg = ForwardDiff.derivative(
                g -> exact_phase_function(hg_contributor(g, ϖ0), cosΘ),
                g0)
            dP_dp[iv, :, :, 1] .= dP_dg
        end
        Jg_chain = @inferred chain_rule_combine_dP(
            hg_jac.path1.P_eff, dP_dp, hg_selector)
        @test Jg_chain ≈ Jg rtol=1e-12 atol=1e-14

        @test_throws ArgumentError chain_rule_combine_dϖ(
            zeros(FT, 1, 1, 1, 2), zeros(FT, 3, 1, 1))
        @test_throws ArgumentError chain_rule_combine_dP(
            zeros(FT, 1, 1, 1, 2), zeros(FT, 2, 2, 1, 1))
        @test_throws ArgumentError chain_rule_combine_dP(
            zeros(FT, 1, 1, 1), zeros(FT, 1, 1, 1))

        surface_geometry = SSGeometry(μ₀=FT(0.79),
                                      μv=FT[0.41, 0.73],
                                      Δϕ=FT[0.2, 1.1])
        surface_base_albedo = FT[0.25, 0.41]
        surface_abs = AbsorptionSSContributor(τ=FT[0.06 0.02; 0.03 0.05])
        surface_config_from_p(p) = ExactSSConfig(
            geometry=surface_geometry,
            surface=LambertianSSSurface(albedo=surface_base_albedo .* p),
            contributors=(surface_abs,),
            I0=FT[1.0, 0.8])
        fp(p) = vec(run_exact_ss(surface_config_from_p(p); paths=:path2).total)
        Jp = ForwardDiff.jacobian(fp, ones(FT, 2))
        surface_jac = run_exact_ss_with_jacobians(
            surface_config_from_p(ones(FT, 2)); paths=:path2).jacobians
        dρ_dp = zeros(FT, 2, 2, 2)
        dρ_dp[:, 1, 1] .= surface_base_albedo[1] / π
        dρ_dp[:, 2, 2] .= surface_base_albedo[2] / π
        J_surface_chain = @inferred chain_rule_combine_surface_brdf(
            surface_jac.path2.surface_brdf, dρ_dp)
        @test size(J_surface_chain) == (2, 1, 2, 2)
        for ip in 1:2
            @test vec(J_surface_chain[:, :, :, ip]) ≈ Jp[:, ip] rtol=1e-12 atol=1e-14
        end
        surface_selector = SSMeasurementSelector(paths=:path2,
                                                 geometry_indices=2,
                                                 spectral_indices=1:2)
        fp_selected(p) = selected_measurements(
            run_exact_ss(surface_config_from_p(p); paths=:path2),
            surface_selector)
        Jp_selected = ForwardDiff.jacobian(fp_selected, ones(FT, 2))
        J_surface_selected = @inferred chain_rule_combine_surface_brdf(
            surface_jac.path2.surface_brdf, dρ_dp, surface_selector)
        @test size(J_surface_selected) == (2, 2)
        @test J_surface_selected ≈ Jp_selected rtol=1e-12 atol=1e-14
        @test_throws ArgumentError chain_rule_combine_surface_brdf(
            zeros(FT, 1, 1, 1), zeros(FT, 2, 1, 1))
        @test_throws ArgumentError chain_rule_combine_surface_brdf(
            zeros(FT, 1, 1, 1), zeros(FT, 1, 2, 1))
        @test_throws ArgumentError chain_rule_combine_surface_brdf(
            zeros(FT, 1, 1), zeros(FT, 1, 1, 1))
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

            cpu32_config = _benchmark_config(Float32; n_geom=8, n_spec=16,
                                             n_layers=4)
            gpu32_config = ExactSSConfig(
                geometry=cpu32_config.geometry,
                surface=cpu32_config.surface,
                contributors=cpu32_config.contributors,
                I0=cpu32_config.I0,
                inner_nquad=cpu32_config.inner_nquad,
                azimuth_nquad=cpu32_config.azimuth_nquad,
                architecture=vSmartMOM.GPU())
            cpu32_result = run_exact_ss(cpu32_config; paths=:all)
            gpu32_result = run_exact_ss(gpu32_config; paths=:all)
            CUDA.synchronize()
            @test eltype(gpu32_result.total) == Float32
            for field in (:total, :path1, :path2, :path3, :path4)
                @test Array(getfield(gpu32_result, field)) ≈
                      getfield(cpu32_result, field) rtol=5e-5 atol=1e-6
            end

            phase = vSmartMOM.Scattering.SyntheticPolarizedHenyeyGreensteinPhaseFunction(
                g=FT(0.35), polarization_fraction=FT(0.4))
            greek = vSmartMOM.Scattering.greek_coefficients(
                phase; l_max=16, nquad=48)
            greek_cpu = ExactSSConfig(
                geometry=SSGeometry(μ₀=FT(0.81),
                                    μv=FT[0.42, 0.68],
                                    Δϕ=FT[0.2, 1.0]),
                surface=LambertianSSSurface(albedo=FT(0.24)),
                contributors=(GreekCoefsSSContributor(
                                  greek_coefs=greek, ϖ=FT(0.9),
                                  τ=reshape(FT[0.02, 0.03], 2, 1)),
                              AbsorptionSSContributor(
                                  τ=reshape(FT[0.01, 0.01], 2, 1))),
                I0=FT[1.0], inner_nquad=8, azimuth_nquad=16)
            greek_gpu = ExactSSConfig(
                geometry=greek_cpu.geometry,
                surface=greek_cpu.surface,
                contributors=greek_cpu.contributors,
                I0=greek_cpu.I0,
                inner_nquad=greek_cpu.inner_nquad,
                azimuth_nquad=greek_cpu.azimuth_nquad,
                architecture=vSmartMOM.GPU())
            greek_cpu_result = run_exact_ss(greek_cpu; paths=:all)
            greek_gpu_result = run_exact_ss(greek_gpu; paths=:all)
            CUDA.synchronize()
            for field in (:total, :path1, :path2, :path3, :path4)
                @test Array(getfield(greek_gpu_result, field)) ≈
                      getfield(greek_cpu_result, field) rtol=5e-12 atol=5e-14
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
            contributors=good.contributors, I0=good.I0, n_stokes=3);
            paths=:all)
        @test size(vSmartMOM.StandaloneSS._precompute_optics(
            Val(3), good).P_eff) == (1, 3, 1, 1)
        hg_vector = ExactSSConfig(
            geometry=good.geometry, surface=good.surface,
            contributors=(HGAerosolSSContributor(g=FT(0.4), ϖ=FT(0.9),
                                                 τ=FT[0.1;;]),),
            I0=FT[1.0], n_stokes=3)
        @test_throws ArgumentError run_exact_ss(hg_vector; paths=:path1)
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
