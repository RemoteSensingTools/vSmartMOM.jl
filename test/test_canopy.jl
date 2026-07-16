using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.IO
using vSmartMOM.InelasticScattering
using CanopyOptics
using LinearAlgebra
using Statistics
using Test

const CO = CanopyOptics

function _canopy_test_quad_points(; FT=Float64, l_trunc=15, pol_type=CoreRT.Stokes_I())
    obs = CoreRT.ObsGeometry(FT(30), FT[0, 30], FT[0, 0], FT[0])
    return CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), l_trunc, obs,
                                 pol_type, Array)
end

function _spherical_reference_Z(scatter, μ; n_φ=96)
    FT = eltype(μ)
    ω = FT(scatter.R + scatter.T)
    φ_az = collect(range(zero(FT), FT(π); length=n_φ + 1))[1:end-1]
    dφ = FT(π) / n_φ
    w_azi = (FT(2) / FT(π)) * dφ

    nμ = length(μ)
    Zpp = zeros(FT, nμ, nμ)
    Zmp = zeros(FT, nμ, nμ)
    for j in 1:nμ, i in 1:nμ
        Ωin = CO.dirVector_μ(μ[j], zero(FT))
        Γdown = zero(FT)
        Γup = zero(FT)
        for φ in φ_az
            Γdown += CO.compute_Γ_isotropic(scatter, Ωin, CO.dirVector_μ( μ[i], φ)) * w_azi
            Γup   += CO.compute_Γ_isotropic(scatter, Ωin, CO.dirVector_μ(-μ[i], φ)) * w_azi
        end
        scale = FT(2) / (ω * FT(0.5)) # spherical LAD has G(μ) ≡ 0.5
        Zpp[i, j] = scale * Γdown
        Zmp[i, j] = scale * Γup
    end
    return Zpp, Zmp
end

@testset "CanopySurface construction and YAML" begin
    soil = LambertianSurfaceScalar(0.1)
    canopy = CanopySurface(; soil, LAI=3.0,
                            canopy_clumping=CO.ConstantClumping(Ω=0.75))

    @test canopy.LAI == 3.0
    @test canopy.n_layers == 1
    @test canopy.soil === soil
    @test canopy.canopy_quadrature isa CO.CanopyQuadrature
    @test canopy.canopy_clumping isa CO.ConstantClumping{Float64}
    @test canopy._cache === nothing

    params = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    @test params.brdf[1] isa CanopySurface
    @test params.brdf[1].soil isa LambertianSurfaceScalar
    @test params.brdf[1].canopy_quadrature.n_azimuth == CO.CanopyQuadrature().n_azimuth
end

@testset "Analytic canopy Z invariants" begin
    FT = Float64
    μ, w = CO.gauleg(8, zero(FT), one(FT))
    μ = collect(μ)
    w = collect(w)
    scatter = CO.BiLambertianCanopyScattering(R=FT(0.5), T=FT(0.5))
    LAD = CO.spherical_leaves(FT)
    quadrature = CO.CanopyQuadrature(n_leaf=32, n_azimuth=32)

    Zpp, Zmp = CO.compute_Z_matrices_aniso_analytic(scatter, μ, LAD, 2;
                                                    quadrature)
    @test all(isfinite.(Zpp))
    @test all(isfinite.(Zmp))
    @test minimum(Zpp[:, :, 1]) ≥ -1e-12
    @test minimum(Zmp[:, :, 1]) ≥ -1e-12

    fluxes = [sum(w .* (Zpp[:, j, 1] .+ Zmp[:, j, 1])) for j in eachindex(μ)]
    @test all(isapprox.(fluxes, FT(2); rtol=0.03, atol=0.03))

    Zpp_ref, Zmp_ref = _spherical_reference_Z(scatter, μ; n_φ=96)
    @test norm(Zpp[:, :, 1] - Zpp_ref) / norm(Zpp_ref) < 5e-3
    @test norm(Zmp[:, :, 1] - Zmp_ref) / norm(Zmp_ref) < 5e-3
end

@testset "Spectral canopy single-point coarse grid" begin
    qp = _canopy_test_quad_points(l_trunc=9)
    soil = LambertianSurfaceScalar(0.1)
    canopy = CanopySurface(; soil, LAI=2.0,
                            leaf_reflectance=[0.40, 0.45, 0.48],
                            leaf_transmittance=[0.05, 0.08, 0.10],
                            leaf_optics_grid=[760.0, 770.0, 780.0],
                            grid_unit=:nm)
    R_interp, T_interp, wn_coarse, Zpp, Zmp =
        CoreRT._build_spectral_canopy_cache(canopy, collect(qp.qp_μ), collect(qp.wt_μ),
                                            [1e7 / 770.0], 1, Float64, 1)

    @test length(R_interp) == 1
    @test length(T_interp) == 1
    @test length(wn_coarse) ≥ 2
    @test size(Zpp, 3) == length(wn_coarse)
    @test size(Zmp, 3) == length(wn_coarse)
end

@testset "CanopySurface Stokes Z expansion" begin
    pol = CoreRT.Stokes_IQUV{Float64}()
    qp = _canopy_test_quad_points(l_trunc=7, pol_type=pol)
    μ = collect(qp.qp_μ)
    soil = LambertianSurfaceScalar(0.1)
    diffuse = CO.BiLambertianCanopyScattering(R=0.40, T=0.10)
    specular = CO.SpecularCanopyScattering(nᵣ=1.5, κ=0.2)
    quadrature = CO.CanopyQuadrature(n_leaf=16, n_azimuth=8)

    diffuse_canopy = CanopySurface(; soil, LAI=1.0, canopy_scattering=diffuse,
                                    canopy_quadrature=quadrature)
    # Phase B: fourth arg is now `m_max` (order). m_max=1 ⇒ 2 moments,
    # matching the original (count=2) test expectation at line 117.
    Zpp_d, Zmp_d = CoreRT._compute_canopy_Z_stack(
        diffuse_canopy, diffuse, μ, 1, Float64, pol.n)

    @test size(Zpp_d) == (length(μ) * pol.n, length(μ) * pol.n, 2)
    for si in 1:pol.n, sj in 1:pol.n
        (si, sj) == (1, 1) && continue
        @test all(iszero, Zpp_d[si:pol.n:end, sj:pol.n:end, :])
        @test all(iszero, Zmp_d[si:pol.n:end, sj:pol.n:end, :])
    end

    composite = diffuse + specular
    specular_canopy = CanopySurface(; soil, LAI=1.0, canopy_scattering=composite,
                                     canopy_quadrature=quadrature)
    # Phase B: fourth arg is `m_max` (order). m_max=1 ⇒ 2 moments.
    Zpp_c, Zmp_c = CoreRT._compute_canopy_Z_stack(
        specular_canopy, composite, μ, 1, Float64, pol.n)

    @test size(Zpp_c) == size(Zpp_d)
    @test all(isfinite, Zpp_c)
    @test all(isfinite, Zmp_c)
    @test any(abs.(Zpp_c[2:pol.n:end, 1:pol.n:end, :]) .> 1e-12) ||
          any(abs.(Zmp_c[2:pol.n:end, 1:pol.n:end, :]) .> 1e-12)
end

@testset "CanopySurface forward RT smoke" begin
    params = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()

    model = model_from_parameters(params)
    @test model.surfaces[1] isa CanopySurface

    R = rt_run(model; i_band=1)[1]
    @test ndims(R) == 3
    @test all(isfinite.(R))
    @test all(R[:, 1, :] .> 0)
end

@testset "CanopySurface bright-soil coupling stays positive" begin
    ν_800nm = 1e7 / 800.0
    params = parameters_from_dict(Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => ["[$ν_800nm]"],
            "surface" => ["LambertianSurfaceScalar(0.0)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 4,
            "Δ_angle" => 2.0,
            "l_trunc" => 20,
            "depol" => -1.0,
            "float_type" => "Float64",
            "architecture" => "CPU()",
        ),
        "geometry" => Dict{String,Any}(
            "sza" => 30.0,
            "vza" => [30.0, 0.0, 30.0],
            "vaz" => [180.0, 0.0, 0.0],
            "obs_alt" => [0.0],
        ),
        "atmospheric_profile" => Dict{String,Any}(
            "T" => [285.0],
            "p" => [1012.99, 1013.0],
            "profile_reduction" => -1,
        ),
    ))

    params.brdf[1] = CanopySurface(;
        soil = LambertianSurfaceScalar(1.0),
        LAI = 4.0,
        n_layers = 1,
        LAD = CO.planophile_leaves2(Float64),
        leaf_reflectance = 0.4424617419246772,
        leaf_transmittance = 0.4745523417367442,
    )

    R = rt_run(model_from_parameters(params); i_band=1)[1]
    @test all(isfinite.(R))
    @test minimum(R[:, 1, 1]) >= -1e-10
    @test maximum(R[:, 1, 1]) < 1.0
end

@testset "CanopySurface linearized soil-albedo Jacobian axis regression" begin
    # Build a minimal canopy RT problem: 1 spectral point, 6 streams, Stokes_I.
    # With no aerosols and no gases, Nparams = n_surface = 1 and iparam = 1.
    # The correct layout for AddedLayerLin is [nμ, nμ, nSpec, Nparams] so the
    # soil-albedo derivative must land in [:,:,:,1], NOT [1,:,:,:].
    # Before the fix, iparam was used as the FIRST axis, corrupting a row of
    # the nμ dimension and leaving the actual parameter slot zero.

    FT = Float64
    pol_type = CoreRT.Stokes_I{FT}()
    n_stokes = pol_type.n  # = 1
    arch = vSmartMOM.Architectures.CPU()
    arr_type = Array

    qp = CoreRT.rt_set_streams(CoreRT.RadauQuad(), 4, CoreRT.ObsGeometry(FT(30), FT[30.0], FT[0.0], FT[0]), pol_type, arr_type)
    nμ   = size(qp.qp_μ, 1) * n_stokes
    nSpec = 1

    soil = LambertianSurfaceScalar(FT(0.15))
    canopy = CanopySurface(; soil, LAI=FT(2.0), n_layers=1,
                            leaf_reflectance=FT(0.4),
                            leaf_transmittance=FT(0.05))

    # Allocate AddedLayer and AddedLayerLin (Nparams=1, matching one surface param)
    Nparams = 1
    dims = (nμ, nμ)

    al, all_lin = CoreRT.make_added_layer(LinMode(), noRS{FT}(), FT, arr_type,
                                          Nparams, dims, nSpec)

    # Zero optical depth (τ̇_sum shape: [nSpec, Nparams] — unused by canopy lin path)
    τ_sum  = zeros(FT, nSpec)
    τ̇_sum  = zeros(FT, nSpec, Nparams)

    CoreRT.create_surface_layer!(noRS{FT}(), canopy, al, all_lin, 1,
                                 false, 0, pol_type, qp, τ_sum, τ̇_sum,
                                 zeros(FT, nSpec), arch)

    # The soil-albedo Jacobian must be in the LAST axis (iparam=1 → [:,:,:,1])
    deriv_correct = all_lin.ap_ṙ⁻⁺[:, :, :, 1]
    # Before the fix, [1,:,:,:] would hold the (incorrectly placed) data instead
    deriv_wrong_slot = all_lin.ap_ṙ⁻⁺[1, :, :, :]  # old axis — should be zero or differ

    # The correct parameter slot must be non-trivially populated (canopy with
    # soil albedo 0.15 has a nonzero soil-albedo Jacobian)
    @test any(abs.(deriv_correct) .> 1e-10)

    # Verify the wrong slot (first row of nμ dim) is NOT the sole nonzero holder.
    # If the fix is correct, deriv_correct sums to strictly more than deriv_wrong_slot
    # summed over a full [:,:,:,1] slice, because the correct slot has a full nμ×nμ
    # matrix while the wrong slot is just one row.
    @test sum(abs.(deriv_correct)) > sum(abs.(deriv_wrong_slot))

    # Shape invariant: parameter axis is last
    @test size(all_lin.ap_ṙ⁻⁺) == (nμ, nμ, nSpec, Nparams)
end
