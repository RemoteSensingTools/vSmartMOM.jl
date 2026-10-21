# Fixed radius grid for the log-normal Mie paths (v2.2 decision 11).
#
# The radius nodes depend only on `r_min`, `r_max` and `nquad_radius`, never on
# the size-distribution parameters. The analytic size tangents (which act on
# the normalized quadrature weights only) are then the exact derivatives of the
# discretized optics at any node count. With the former quantile-based lower
# bound, `r_min = quantile(dist, 1e-8)`, the nodes moved with μ_logr/σ_logr and
# central differences of an under-resolved quadrature disagreed with the
# tangents. For the aerosol below at the JacobianTestFast radius settings the
# release audit (BJ-4) found size-tangent errors of 4–14 % in ω̃̇ and up to 7×
# in the aerosol optical-depth tangent.

using Test, Logging, YAML
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using Distributions: LogNormal

const _SC = vSmartMOM.Scattering

# The release audit's aerosol (rₘ = 0.2 µm, σ_g = 1.6, n = 1.3 − 0.005i) on the
# JacobianTestFast radius settings: r_max = 50 µm, 100 nodes.
const _FG_r_max = 50.0
const _FG_nquad = 100
_fg_aerosol(μ=log(0.2), σ=log(1.6); nᵣ=1.3, nᵢ=0.005) = Aerosol(LogNormal(μ, σ), nᵣ, nᵢ)
_fg_model(a, λ; kw...) = make_mie_model(NAI2(), a, λ, Stokes_IQU(), NoTruncation(),
                                        _FG_r_max, _FG_nquad; kw...)

@testset "fixed radius grid: NAI2, linearized NAI2, PCW" begin

    @testset "nodes do not depend on the size distribution" begin
        m₁ = _fg_model(_fg_aerosol(), 0.77)
        m₂ = _fg_model(_fg_aerosol(log(1.2), log(2.0)), 0.77)
        @test m₁.r_min == _SC.DEFAULT_R_MIN == 1e-3
        r₁, w₁ = _SC.radius_quadrature(m₁)
        r₂, w₂ = _SC.radius_quadrature(m₂)
        @test r₁ == r₂ && w₁ == w₂
        @test (r₁, w₁) == _SC.gauleg_log(_FG_nquad, 1e-3, _FG_r_max; norm=false)
        @test 1e-3 < minimum(r₁) && maximum(r₁) < _FG_r_max

        m₃ = _fg_model(_fg_aerosol(), 0.77; r_min=0.01)
        @test m₃.r_min == 0.01
        @test minimum(first(_SC.radius_quadrature(m₃))) > 0.01
        @test_throws ArgumentError _fg_model(_fg_aerosol(), 0.77; r_min=_FG_r_max)
        @test_throws ArgumentError _fg_model(_fg_aerosol(), 0.77; r_min=0.0)

        # Float32 models carry a Float32 bound; the other constructors default it.
        a32 = Aerosol(LogNormal(log(0.2f0), log(1.6f0)), 1.3f0, 0.005f0)
        m32 = make_mie_model(NAI2(), a32, 0.77f0, Stokes_IQU(), NoTruncation(), 50.0f0, 100)
        @test m32.r_min === Float32(_SC.DEFAULT_R_MIN)
        o32 = compute_aerosol_optical_properties(m32)
        @test o32.k isa Float32
        @test o32.k ≈ compute_aerosol_optical_properties(m₁).k rtol=1e-5
        mk = MieModel(computation_type=NAI2(), aerosol=_fg_aerosol(), λ=0.77,
                      polarization_type=Stokes_IQU(), truncation_type=NoTruncation(),
                      r_max=_FG_r_max, nquad_radius=_FG_nquad)
        @test mk.r_min == _SC.DEFAULT_R_MIN
        @test MieModel(computation_type=NAI2(), aerosol=_fg_aerosol(), λ=0.77,
                       polarization_type=Stokes_IQU(), truncation_type=NoTruncation(),
                       r_max=_FG_r_max, nquad_radius=_FG_nquad, r_min=0.02).r_min == 0.02
        mp = MieModel(NAI2(), _fg_aerosol(), 0.77, Stokes_IQU(), NoTruncation(), _FG_r_max,
                      _FG_nquad, zeros(1, 1, 1), zeros(1, 1, 1), vSmartMOM.Architectures.CPU(),
                      nothing)
        @test mp.r_min == _SC.DEFAULT_R_MIN
        @test compute_aerosol_optical_properties(mp).k == compute_aerosol_optical_properties(m₁).k
    end

    @testset "optical-property size tangents equal central differences ($_FG_nquad nodes)" begin
        a₀ = _fg_aerosol()
        μ₀, σ₀ = a₀.size_distribution.μ, a₀.size_distribution.σ
        h = 1e-5
        for λ in (0.55, 0.77)
            optics, lin_optics = compute_aerosol_optical_properties(LinMode(), _fg_model(a₀, λ), Float64)
            k_ref, k̇_ref = compute_ref_aerosol_extinction(LinMode(), _fg_model(a₀, λ), Float64)
            @test k_ref == optics.k
            for (slot, δμ, δσ) in ((3, h, 0.0), (4, 0.0, h))
                p = compute_aerosol_optical_properties(_fg_model(_fg_aerosol(μ₀ + δμ, σ₀ + δσ), λ))
                m = compute_aerosol_optical_properties(_fg_model(_fg_aerosol(μ₀ - δμ, σ₀ - δσ), λ))
                fd_k = (p.k - m.k) / 2h
                fd_ϖ = (p.ω̃ - m.ω̃) / 2h
                L = min(length(p.greek_coefs.β), length(m.greek_coefs.β))
                @test lin_optics.k̇[slot] ≈ fd_k rtol=1e-7
                @test k̇_ref[slot] ≈ fd_k rtol=1e-7
                @test lin_optics.ω̃̇[slot] ≈ fd_ϖ rtol=1e-6
                for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    fd = (getproperty(p.greek_coefs, f)[1:L] .-
                          getproperty(m.greek_coefs, f)[1:L]) ./ 2h
                    an = getproperty(lin_optics.lin_greek_coefs, Symbol(f, '̇'))[slot, 1:L]
                    @test maximum(abs, an .- fd) ≤ 1e-6 * maximum(abs, fd)
                end
            end
        end
    end

    @testset "τ_aer and TOA size Jacobians at JacobianTestFast settings" begin
        # Gas absorption does not touch the radius grid; dropping it keeps the
        # five model builds fast.
        cfg = YAML.load_file(joinpath(@__DIR__, "test_parameters", "JacobianTestFast.yaml"))
        delete!(cfg, "absorption")
        aerosol = only(cfg["scattering"]["aerosols"])
        aerosol["μ"] = 0.2; aerosol["σ"] = 1.6; aerosol["nᵢ"] = 0.005
        @test cfg["scattering"]["r_max"] == _FG_r_max
        @test cfg["scattering"]["nquad_radius"] == _FG_nquad
        params = with_logger(NullLogger()) do
            read_parameters(cfg)
        end
        @test params.scattering_params.r_min == _SC.DEFAULT_R_MIN
        with_logger(NullLogger()) do
            model, lin = model_from_parameters(LinMode(), params)
            result = rt_run(model, lin, 1, size(lin.τ̇_abs[1], 1), 1)
            layout = result.layout
            h = 1e-5
            for (slot, δμ, δσ) in ((4, h, 0.0), (5, 0.0, h))
                perturbed(s) = begin
                    p = deepcopy(params)
                    a = only(p.scattering_params.rt_aerosols).aerosol
                    d = a.size_distribution
                    a.size_distribution = LogNormal(d.μ + s * δμ, d.σ + s * δσ)
                    m = model_from_parameters(p)
                    (τ=m.τ_aer[1], toa=rt_run(m).toa)
                end
                plus, minus = perturbed(1), perturbed(-1)
                # τ_aer = τ_ref k(λ)/k_ref nearly cancels near λ_ref, so its
                # size tangent is the most sensitive check of node motion.
                fd_τ = (plus.τ .- minus.τ) ./ 2h
                an_τ = lin.τ̇_aer[1][:, slot, :, :]
                @test maximum(abs, fd_τ) > 0
                @test maximum(abs, an_τ .- fd_τ) ≤ 1e-6 * maximum(abs, fd_τ)
                fd = (plus.toa .- minus.toa) ./ 2h
                an = result.toa_jacobian[:, :, :, CoreRT.aerosol_range(layout, 1)[slot]]
                @test maximum(abs, fd) > 0
                @test maximum(abs, an .- fd) ≤ 1e-6 * maximum(abs, fd)
            end
        end
    end

    @testset "scattering.r_min configuration" begin
        base = YAML.load_file(joinpath(@__DIR__, "test_parameters", "JacobianTestFast.yaml"))
        delete!(base, "absorption")
        aerosol = only(base["scattering"]["aerosols"])
        aerosol["μ"] = 0.2; aerosol["σ"] = 1.6; aerosol["nᵢ"] = 0.005
        parse_cfg(c) = with_logger(NullLogger()) do
            read_parameters(c)
        end
        @test parse_cfg(base).scattering_params.r_min === _SC.DEFAULT_R_MIN
        r_min = 0.02   # below 5e-7 of this mode's particle number; a coarser grid than the default
        cfg = deepcopy(base); cfg["scattering"]["r_min"] = r_min
        params = parse_cfg(cfg)
        @test params.scattering_params.r_min == r_min
        for bad in (0.0, -1e-3, 50.0, 80.0, "1e-3", true)
            c = deepcopy(base); c["scattering"]["r_min"] = bad
            @test_throws ArgumentError parse_cfg(c)
        end

        # The configured bound reaches every Mie call of the model builders,
        # the BatchContext reference-extinction cache and its microphysics update.
        sp = params.scattering_params
        a = only(sp.rt_aerosols).aerosol
        mie(aer, λ; kw...) = make_mie_model(NAI2(), aer, λ, params.polarization_type,
                                            NoTruncation(), sp.r_max, sp.nquad_radius; kw...)
        ref_aerosol(aer) = Aerosol(aer.size_distribution, real(sp.n_ref), -imag(sp.n_ref))
        λ₁ = 1e4 / params.spec_bands[1][1]   # first anchor wavelength of the band
        with_logger(NullLogger()) do
            expected = compute_aerosol_optical_properties(mie(a, λ₁; r_min))
            default = compute_aerosol_optical_properties(mie(a, λ₁))
            @test abs(expected.k / default.k - 1) > 1e-6   # the two grids differ visibly
            expected_lin = last(compute_aerosol_optical_properties(LinMode(), mie(a, λ₁; r_min), Float64))
            forward = model_from_parameters(params)
            model, lin = model_from_parameters(LinMode(), params)
            @test first(forward.aerosol_optics[1][1].k) ≈ expected.k rtol=1e-12
            @test first(model.aerosol_optics[1][1].k) ≈ expected.k rtol=1e-12
            @test lin.lin_aerosol_optics[1][1].k̇[:, 1] ≈ expected_lin.k̇ rtol=1e-12
            @test model.τ_aer[1] ≈ forward.τ_aer[1] rtol=1e-12

            ctx = CoreRT.BatchContext(params)
            k_ref = compute_ref_aerosol_extinction(mie(ref_aerosol(a), sp.λ_ref; r_min), Float64)
            @test only(ctx.k_ref) ≈ k_ref rtol=1e-14
            new_aerosol = Aerosol(LogNormal(log(0.3), log(1.5)), 1.4, 0.002)
            update_aerosol_microphysics!(ctx, 1, new_aerosol)
            fresh_params = deepcopy(params)
            only(fresh_params.scattering_params.rt_aerosols).aerosol = new_aerosol
            fresh = model_from_parameters(fresh_params)
            @test ctx.model.τ_aer[1] ≈ fresh.τ_aer[1] rtol=1e-12
            @test first(ctx.model.aerosol_optics[1][1].k) ≈ first(fresh.aerosol_optics[1][1].k) rtol=1e-12
            @test only(ctx.k_ref) ≈ compute_ref_aerosol_extinction(
                mie(ref_aerosol(new_aerosol), sp.λ_ref; r_min), Float64) rtol=1e-14
        end

        # The six-argument ScatteringParameters constructor keeps working.
        old = CoreRT.ScatteringParameters(sp.rt_aerosols, sp.r_max, sp.nquad_radius,
                                          sp.λ_ref, sp.n_ref, sp.decomp_type)
        @test old.r_min === _SC.DEFAULT_R_MIN
    end

    @testset "PCW uses the same fixed nodes" begin
        # Coarse grid: any difference in node placement would show at ≫ 1e-10.
        a = _fg_aerosol(log(0.3), log(1.5))
        λ, r_max, nq = 0.77, 2.0, 16
        nai2 = compute_aerosol_optical_properties(
            make_mie_model(NAI2(), a, λ, Stokes_IQUV(), NoTruncation(), r_max, nq))
        N_max = _SC.get_n_max(2π * r_max / λ)
        wigner_A, wigner_B = compute_wigner_values(N_max)
        pcw = with_logger(NullLogger()) do
            compute_aerosol_optical_properties(
                make_mie_model(PCW(), a, λ, Stokes_IQUV(), NoTruncation(), r_max, nq,
                               wigner_A, wigner_B))
        end
        @test pcw.k ≈ nai2.k rtol=1e-12
        @test pcw.ω̃ ≈ nai2.ω̃ rtol=1e-12
    end
end
