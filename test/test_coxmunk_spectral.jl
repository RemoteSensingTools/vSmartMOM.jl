using Test, LinearAlgebra, YAML
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using vSmartMOM.StandaloneSS: exact_ss_config_from_model

# Cox–Munk water refractive index (decision 9 of the v2.2 decision log):
# Segelstein (1981) table, interpolated linearly in wavelength at every
# spectral point; a vector `n_water` is used element-wise.
const CMS = vSmartMOM.CoreRT

"Quickstart-based Cox–Munk model on `ν` (cm⁻¹)."
function _coxmunk_spectral_params(ν; surface = "CoxMunkSurface(wind_speed=5.0)",
                                  float_type = "Float64", polarization = "Stokes_IQU()")
    d = YAML.load_file(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
    d["radiative_transfer"]["surface"] = [surface]
    d["radiative_transfer"]["polarization_type"] = polarization
    d["radiative_transfer"]["float_type"] = float_type
    d["radiative_transfer"]["spec_bands"] = ["[" * join(string.(ν), ", ") * "]"]
    d["geometry"]["sza"] = 30.0
    d["geometry"]["vza"] = [30.0, 45.0]
    d["geometry"]["vaz"] = [0.0, 35.0]
    return read_parameters(d)
end

@testset "Segelstein (1981) water refractive index" begin
    # Rows of the stored table are returned exactly.
    i = findfirst(λ -> isapprox(λ, 549.54086; rtol = 1e-12), CMS._WATER_λ_NM)
    @test (CMS._WATER_N_REAL[i], CMS._WATER_K_IMAG[i]) == (1.335972, 2.4422501e-9)
    @test CMS.water_refractive_index(CMS._WATER_λ_NM[i]) == complex(1.335972, 2.4422501e-9)
    # Second redistribution (omlc.org: λ and k to 4 significant digits).
    for (λ_nm, n, k) in ((549.5, 1.335972, 2.442e-9), (760.3, 1.326909, 1.580e-7),
                         (1600.0, 1.309642, 9.347e-5), (2061.0, 1.293973, 6.206e-4))
        row = argmin(abs.(CMS._WATER_λ_NM .- λ_nm))
        @test abs(CMS._WATER_λ_NM[row] / λ_nm - 1) < 5e-4
        m = CMS.water_refractive_index(CMS._WATER_λ_NM[row])
        @test real(m) == n
        @test imag(m) ≈ k rtol = 1e-3
    end
    # Linear in wavelength between rows, for n and k.
    λ₁, λ₂ = CMS._WATER_λ_NM[300], CMS._WATER_λ_NM[301]
    m₁, m₂ = CMS.water_refractive_index(λ₁), CMS.water_refractive_index(λ₂)
    @test CMS.water_refractive_index(0.25λ₁ + 0.75λ₂) ≈ 0.25m₁ + 0.75m₂ rtol = 1e-14
    # The whole solar range is covered; beyond the table is an error.
    @test first(CMS._WATER_λ_NM) < 200 && last(CMS._WATER_λ_NM) > 4000
    @test all(λ -> isfinite(CMS.water_refractive_index(λ)), 200.0:7.0:4000.0)
    @test_throws ArgumentError CMS.water_refractive_index(10.0)
    @test_throws ArgumentError CMS.water_refractive_index(2e10)
end

@testset "Per-spectral-point n_water resolution" begin
    ν = collect(range(1e7 / 775, 1e7 / 757, length = 5))
    s = CoxMunkSurface(wind_speed = 5.0)
    r = CMS._resolve_spectral_surface(s, ν)
    @test r.n_water == ComplexF64.(CMS.water_refractive_index.(1e7 ./ ν))
    @test (r.wind_speed, r.whitecap_albedo, r.include_whitecaps, r.shadowing) ==
          (s.wind_speed, s.whitecap_albedo, s.include_whitecaps, s.shadowing)
    s32 = CoxMunkSurface{Float32}(wind_speed = 5f0)
    @test eltype(CMS._resolve_spectral_surface(s32, ν).n_water) === ComplexF32
    # A scalar keeps its meaning; a vector is checked against the grid.
    scalar = CoxMunkSurface(wind_speed = 5.0, n_water = 1.33 + 0im)
    @test CMS._resolve_spectral_surface(scalar, ν) === scalar
    @test_throws ArgumentError CMS._resolve_spectral_surface(
        CoxMunkSurface(wind_speed = 5.0, n_water = r.n_water[1:3]), ν)
    lamb = LambertianSurfaceScalar(0.1)
    @test CMS._resolve_spectral_surface(lamb, ν) === lamb
    # A vector with equal entries is the scalar case.
    @test CMS._spectral_n_water(CoxMunkSurface(wind_speed = 5.0,
                                               n_water = fill(1.33 + 0im, 5)), 5) === 1.33 + 0im
    # Without a wavelength, `nothing` has no value.
    @test_throws ArgumentError CMS._spectral_n_water(s, 5)
    @test_throws ArgumentError CMS.coxmunk_brdf_mueller(s, 3, 0.8, 0.7, 0.3)
    @test_throws ArgumentError CMS.reflectance(r, Stokes_I(), [0.5, 0.8], 0)
    @test CMS.coxmunk_brdf_mueller(r, 3, 0.8, 0.7, 0.3; n_water = r.n_water[2]) ==
          CMS.coxmunk_brdf_mueller(CoxMunkSurface(wind_speed = 5.0, n_water = r.n_water[2]),
                                   3, 0.8, 0.7, 0.3)
end

@testset "Operator knots in n_water" begin
    # Constant index: one knot, no interpolation.
    knots, lo, hi, t = CMS._n_water_knots(fill(1.33 + 1e-8im, 50))
    @test knots == [1] && all(==(1), lo) && all(==(1), hi) && all(iszero, t)
    # Collinear indices with a small total change: the two end points.
    line = [1.33 + 0.0002s / 99 + (1e-7 + 1e-9s) * im for s in 0:99]
    knots, lo, hi, t = CMS._n_water_knots(line)
    @test knots == [1, 100]
    @test t ≈ (0:99) ./ 99 rtol = 1e-12
    # Knots reproduce their own index; the interpolation reproduces each
    # point's Fresnel elements (the knot criterion).
    ν = collect(range(1e7 / 2083, 1e7 / 2042, length = 4000))   # 2.06 µm band
    n = ComplexF64.(CMS.water_refractive_index.(1e7 ./ ν))
    knots, lo, hi, t = CMS._n_water_knots(n)
    @test first(knots) == 1 && last(knots) == length(n) && issorted(knots)
    @test 2 < length(knots) < 40
    @test all(t[knots] .∈ Ref((0.0, 1.0)))
    worst = 0.0
    for c in CMS._N_WATER_PROXY_COS, s in eachindex(n)
        fa, fe = CMS._fresnel_elements(n[knots[lo[s]]], c), CMS._fresnel_elements(n[knots[hi[s]]], c)
        fs = CMS._fresnel_elements(n[s], c)
        worst = max(worst, maximum(abs.((1 - t[s]) .* fa .+ t[s] .* fe .- fs)) / fs[1])
    end
    @test worst <= 1.0001 * CMS._N_WATER_TOL
    # An index close to 1 (where the reflectance is most sensitive) needs more knots.
    near1 = [1.01 + 0.0003s / 99 + 0im for s in 0:99]
    @test length(first(CMS._n_water_knots(near1))) > 2
    # Reviewer counterexamples: n ≈ 2 (no curvature at normal incidence) and a
    # real part below 1 (critical angle): the middle point must be a knot.
    @test first(CMS._n_water_knots(ComplexF64[1.99, 2.0, 2.01])) == [1, 2, 3]
    n₀ = 0.9936288716722868
    @test first(CMS._n_water_knots(ComplexF64[n₀ - 1e-6, n₀, n₀ + 1e-6])) == [1, 2, 3]
    # The plan is reused for equal values and rebuilt when any value changes,
    # also in place and for long vectors (whose `hash` samples few elements).
    plan = CMS._n_water_knot_plan(n)
    @test CMS._n_water_knot_plan(copy(n)) === plan
    long = [1.33 + 1e-9s + 0im for s in 1:50_000]
    long_plan = CMS._n_water_knot_plan(long)
    @test long_plan == CMS._n_water_knots(long)
    long[3] = 0.99 + 0im
    @test CMS._n_water_knot_plan(long) == CMS._n_water_knots(long) != long_plan
end

"""
Worst relative error of the knot-interpolated Fourier reflectance for the
indices `n` (or the Segelstein index on a band around `λ_μm`, 1 % wide,
`n_spec` points), against the reflectance at each point's own index, at the
points farthest from their knots.
"""
function _knot_operator_error(λ_μm::Real, ::Type{FT}; n_spec = 3000) where {FT}
    ν = collect(range(1e4 / (1.005λ_μm), 1e4 / (0.995λ_μm), length = n_spec))
    return _knot_operator_error(Complex{FT}.(CMS.water_refractive_index.(1e7 ./ ν)))
end

function _knot_operator_error(n::AbstractVector{Complex{FT}}) where {FT}
    knots, lo, hi, t = CMS._n_water_knots(n)
    surf = CoxMunkSurface{FT}(wind_speed = FT(3), n_water = n)
    μ = FT[0.12, 0.4, 0.7, 0.97]
    pol = Stokes_IQU{FT}()
    worst = 0.0
    for m in (0, 2)
        Rk = [CMS.reflectance(surf, pol, μ, m; n_water = n[k]) for k in knots]
        for q in 1:(length(knots) - 1)
            knots[q + 1] - knots[q] > 1 || continue
            s = knots[q] + argmin(abs.(t[(knots[q] + 1):(knots[q + 1] - 1)] .- 0.5))
            R_point = CMS.reflectance(CoxMunkSurface{Float64}(wind_speed = 3.0), Stokes_IQU(),
                                      Float64.(μ), m; n_water = ComplexF64(n[s]))
            R_interp = (1 - t[s]) * Rk[lo[s]] + t[s] * Rk[hi[s]]
            worst = max(worst, maximum(abs.(R_interp .- R_point)) / maximum(abs.(R_point)))
        end
    end
    return worst, length(knots)
end

@testset "Knot interpolation accuracy, 0.2–4 µm" begin
    bands = (0.25, 0.45, 0.55, 0.76, 0.87, 1.24, 1.38, 1.6, 1.88, 2.06, 2.3, 2.6, 2.75,
             2.8, 2.94, 3.2, 3.9)
    errors64 = [_knot_operator_error(λ, Float64) for λ in bands]
    @info "Cox–Munk knot interpolation (Float64): worst relative operator error, knots" bands errors64
    @test maximum(first.(errors64)) < 1e-7
    errors32 = [_knot_operator_error(λ, Float32; n_spec = 1000) for λ in (0.55, 1.6, 2.06, 2.8)]
    @info "Cox–Munk knot interpolation (Float32)" errors32
    @test maximum(first.(errors32)) < 1e-6
    # User vectors far from water: high, low and absorbing indices.
    others = ([2.0 + 0.01x + 0im for x in range(-1, 1, length = 201)],
              [1.5 + 0.02x + 0im for x in range(0, 1, length = 201)],
              [1.02 + 0.005x + 0im for x in range(0, 1, length = 201)],
              [complex(1.3 + 0.02x, 0.1 + 0.05x) for x in range(0, 1, length = 201)])
    errors_other = [_knot_operator_error(ComplexF64.(n)) for n in others]
    @info "Cox–Munk knot interpolation (other indices)" errors_other
    @test maximum(first.(errors_other)) < 1e-7
end

@testset "Spectral Cox–Munk RT" begin
    # Seven points across the O₂ A band: knots and interpolated points.
    ν = collect(range(1e7 / 775, 1e7 / 757, length = 7))
    p = _coxmunk_spectral_params(ν)
    model, lin = model_from_parameters(LinMode(), p)
    νg = model.atmosphere.spec_bands[1]
    n_point = ComplexF64.(CMS.water_refractive_index.(1e7 ./ νg))
    knots = first(CMS._n_water_knots(n_point))
    @test 1 < length(knots) < length(νg)          # some points are interpolated
    R = rt_run(model).toa

    # Same radiances as one run per point with that point's scalar index.
    R_point = similar(R)
    original = model.surfaces[1]
    for s in eachindex(νg)
        model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0, n_water = n_point[s])
        R_point[:, :, s] = rt_run(model).toa[:, :, s]
    end
    model.surfaces[1] = original
    scale = maximum(abs, R_point[:, 1, :])
    @test maximum(abs.(R[:, :, knots] .- R_point[:, :, knots])) / scale < 1e-12
    @test maximum(abs.(R .- R_point)) / scale < 1e-7

    # A vector n_water is used element-wise; equal entries are the scalar case.
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0, n_water = n_point)
    @test rt_run(model).toa == R
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0, n_water = fill(n_point[3], 7))
    R_uniform = rt_run(model).toa
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0, n_water = n_point[3])
    @test rt_run(model).toa == R_uniform
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0, n_water = n_point[1:6])
    @test_throws ArgumentError rt_run(model)
    model.surfaces[1] = original

    # Atmosphere/surface split and linearized solver use the same operator.
    cache = rt_run_atmosphere(model)
    @test rt_run_surface(cache, original)[1] == R
    ngas = size(lin.τ̇_abs[1], 1)
    L = rt_run(model, lin, 0, ngas, 1)
    @test maximum(abs.(L[1] .- R)) / scale < 1e-13
    h = 1e-4
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0 + h); plus = rt_run(model).toa
    model.surfaces[1] = CoxMunkSurface(wind_speed = 5.0 - h); minus = rt_run(model).toa
    model.surfaces[1] = original
    wind = first(CoreRT.surface_range(L.layout))
    fd = (plus - minus) / (2h)
    @test maximum(abs.(L[3][:, :, :, wind] .- fd)) / maximum(abs, fd) < 1e-6

    # Purely absorbing atmosphere: the TOA field is the direct surface path,
    # exact at every point (knot or interpolated) with that point's index.
    fill!(model.τ_rayl[1], 0)
    fill!(model.τ_abs[1], 0.2 / size(model.τ_abs[1], 2))
    actual = rt_run(model).toa
    μ0 = cosd(p.sza)
    for iv in eachindex(p.vza), s in eachindex(νg)
        μv = cosd(p.vza[iv])
        M = CMS.coxmunk_brdf_mueller(original, 3, μ0, μv, deg2rad(p.vaz[iv]);
                                     n_water = n_point[s])
        @test actual[iv, :, s] ≈ μ0 * M[:, 1] * exp(-0.2 * (inv(μ0) + inv(μv))) rtol = 2e-10 atol = 2e-14
    end

    # StandaloneSS adapter sees the same per-point index.
    cfg = exact_ss_config_from_model(model_from_parameters(p))
    @test cfg.surface.n_water == n_point
end

@testset "Spectral Cox–Munk RT in Float32" begin
    ν = collect(range(1e7 / 775, 1e7 / 757, length = 5))
    R64 = rt_run(model_from_parameters(_coxmunk_spectral_params(ν))).toa
    R32 = rt_run(model_from_parameters(_coxmunk_spectral_params(ν; float_type = "Float32"))).toa
    @test eltype(R32) === Float32
    @test maximum(abs.(R32 .- R64)) / maximum(abs, R64[:, 1, :]) < 1e-4
end

@testset "Canopy over a Cox–Munk soil" begin
    # The soil index is resolved per spectral point in every solver that
    # builds the canopy: same results as the explicit per-point vector.
    ν = collect(range(1e7 / 775, 1e7 / 757, length = 3))
    p = _coxmunk_spectral_params(ν; polarization = "Stokes_I()")
    n_point = ComplexF64.(CMS.water_refractive_index.(1e7 ./ ν))
    canopy(soil) = CanopySurface(; soil, LAI = 1.0, n_layers = 1,
                                 leaf_reflectance = 0.4, leaf_transmittance = 0.05)
    results = map((CoxMunkSurface(wind_speed = 5.0),
                   CoxMunkSurface(wind_speed = 5.0, n_water = n_point))) do soil
        p.brdf[1] = canopy(soil)
        model = model_from_parameters(p)
        forward = rt_run(model).toa
        single = rt_run_ss(model)[1]
        p.brdf[1] = canopy(soil)
        model, lin = model_from_parameters(LinMode(), p)
        linearized = rt_run(model, lin, 0, size(lin.τ̇_abs[1], 1), 1)[1]
        (forward, single, linearized)
    end
    @test all(isfinite, results[1][1])
    @test results[1] == results[2]
end
