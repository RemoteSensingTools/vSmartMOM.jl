# =========================================================================
# RT entry-point guards and consistency fixes
# =========================================================================
#
# Regression tests for solve-boundary defects: a reused RS object keeping a
# driver-written F₀; HDR/BHR (and rt_run_ss hem_T) double-counting the direct
# beam when μ₀ is a weighted node; TMS with stream exports and with combined
# bands; the Cox–Munk glint residual on multisensor TOA output; U/V-polarized
# incident beams; ThermalEmission on paths that ignore it or with a wrong
# layer count; and linearized runs on surfaces without a linearized layer.

module RTEntryGuardsTests

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.CoreRT: TMSCorrection, NoSSCorrection, RTNumericalParameters
using vSmartMOM.InelasticScattering
using Distributions: LogNormal
using Statistics
using Test
using YAML

const _QS = joinpath(@__DIR__, "..", "config", "quickstart.yaml")

"Quickstart-derived parameters: pure Rayleigh, two layers, nstreams=3."
function _qs_params(; surface = "LambertianSurfaceScalar(0.15)", pol = "Stokes_I()",
                    sza = 60.0, vza = [60.0, 30.0], vaz = [180.0, 90.0],
                    obs_alt = [0], bands = ["[12987.0]"], nstreams = 3,
                    quadrature = nothing, ss_correction = NoSSCorrection(),
                    float_type = "Float64", aerosol = false)
    d = YAML.load_file(_QS)
    rt = d["radiative_transfer"]
    rt["spec_bands"] = bands
    rt["surface"] = fill(surface, length(bands))
    rt["polarization_type"] = pol
    rt["nstreams"] = nstreams
    rt["float_type"] = float_type
    if aerosol   # one truncated (δBGE) aerosol mode, so TMS is active
        rt["truncation"] = "auto"
        rt["greek_beta_cutoff"] = nothing
        d["scattering"] = Dict(
            "aerosols" => [Dict("τ_ref" => 0.3, "μ" => 0.6, "σ" => 1.6,
                                "nᵣ" => 1.45, "nᵢ" => 0.001,
                                "p₀" => 500.0, "σp" => 100.0)],
            "r_max" => 20.0, "nquad_radius" => 50, "λ_ref" => 0.770,
            "decomp_type" => "NAI2()")
    end
    quadrature === nothing || (rt["quadrature_type"] = quadrature)
    g = d["geometry"]
    g["sza"] = sza; g["vza"] = vza; g["vaz"] = vaz; g["obs_alt"] = obs_alt
    p = read_parameters(d)
    p.numerics = RTNumericalParameters{p.float_type}(ss_correction = ss_correction)
    return p
end

_qs_model(; sources = SolarBeam(), kw...) =
    model_from_parameters(_qs_params(; kw...); sources)

_nspec(model) = size(model.τ_abs[1], 1)
_nlay(model) = size(model.τ_abs[1], 2)
_beam(model, scale) = (F = zeros(CoreRT.polarization_type(model).n, _nspec(model));
                       F[1, :] .= scale; F)

@testset "reused RS objects do not keep a driver-written F₀" begin
    model = _qs_model()
    rs = InelasticScattering.noRS{Float64}()
    ref = CoreRT.rt_run(InelasticScattering.noRS{Float64}(), model, 1).toa

    @test CoreRT.rt_run(rs, model, 1; sources = SolarBeam(F₀ = _beam(model, 2.0))).toa ≈ 2ref
    @test CoreRT.rt_run(rs, model, 1).toa == ref          # back to model.sources

    # model.sources is honoured after an earlier explicit-source call.
    model2 = _qs_model(sources = SolarBeam(F₀ = _beam(model, 2.0)))
    rs2 = InelasticScattering.noRS{Float64}()
    CoreRT.rt_run(rs2, model2, 1; sources = SolarBeam(F₀ = _beam(model, 3.0)))
    @test CoreRT.rt_run(rs2, model2, 1).toa ≈ 2ref

    # Legacy presets keep working: an assigned carrier, and an in-place edit
    # of the carrier a previous solve wrote.
    rs3 = InelasticScattering.noRS{Float64}()
    rs3.F₀ = _beam(model, 3.0)
    @test CoreRT.rt_run(rs3, model, 1).toa ≈ 3ref
    rs4 = InelasticScattering.noRS{Float64}()
    CoreRT.rt_run(rs4, model, 1)
    rs4.F₀ .*= 4
    @test CoreRT.rt_run(rs4, model, 1).toa ≈ 4ref

    # Carrier precision differing from the model (assignment converts).
    model32 = model_from_parameters(_qs_params(float_type = "Float32"))
    ref32 = CoreRT.rt_run(InelasticScattering.noRS{Float64}(), model32, 1).toa
    rs64 = InelasticScattering.noRS{Float64}()
    CoreRT.rt_run(rs64, model32, 1; sources = SolarBeam(F₀ = _beam(model32, 2.0)))
    @test CoreRT.rt_run(rs64, model32, 1).toa == ref32

    # Single-scatter driver.
    ss_ref = CoreRT.rt_run_ss(InelasticScattering.noRS{Float64}(), model, 1)[1]
    rs5 = InelasticScattering.noRS{Float64}()
    CoreRT.rt_run_ss(rs5, model, 1; sources = SolarBeam(F₀ = _beam(model, 2.0)))
    @test CoreRT.rt_run_ss(rs5, model, 1)[1] == ss_ref

    # Interior-height (multisensor) driver.
    model_ms = _qs_model(obs_alt = [0, 5.0])
    ms_ref = CoreRT.rt_run(InelasticScattering.noRS{Float64}(), model_ms, 1)
    rs6 = InelasticScattering.noRS{Float64}()
    CoreRT.rt_run(rs6, model_ms, 1; sources = SolarBeam(F₀ = _beam(model, 2.0)))
    ms = CoreRT.rt_run(rs6, model_ms, 1)
    @test ms.toa == ms_ref.toa
    @test only(ms.levels).upwelling == only(ms_ref.levels).upwelling

    # Linearized driver.
    fwd, lin = model_from_parameters(LinMode(), _qs_params())
    ngas = size(lin.τ̇_abs[1], 1)
    lin_ref = rt_run(InelasticScattering.noRS{Float64}(), fwd, lin, 0, ngas, 1, 1)
    rs7 = InelasticScattering.noRS{Float64}()
    rt_run(rs7, fwd, lin, 0, ngas, 1, 1; sources = SolarBeam(F₀ = _beam(fwd, 2.0)))
    lin_again = rt_run(rs7, fwd, lin, 0, ngas, 1, 1)
    @test lin_again.toa == lin_ref.toa
    @test lin_again.toa_jacobian == lin_ref.toa_jacobian
end

@testset "Raman carriers: RRS and VS objects" begin
    p = parameters_from_yaml(joinpath(@__DIR__, "test_parameters",
                                      "Phase1b_RRS_761-764nm.yaml"))
    p.spec_bands = [(1e7 / 763.1):0.5:(1e7 / 762.9)]
    model = model_from_parameters(p)
    FT = CoreRT.float_type(model)
    ν = model.atmosphere.spec_bands[1]
    ν̃ = mean(ν)
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    npol = CoreRT.polarization_type(model).n
    unit = zeros(FT, npol, length(ν)); unit[1, :] .= 1
    function make_rrs()
        rs = InelasticScattering.RRS(
            n2 = n2, o2 = o2,
            greek_raman = InelasticScattering.GreekCoefs(
                [one(FT)], [one(FT)], [one(FT)], [one(FT)], [one(FT)], [one(FT)]),
            fscattRayl = [one(FT)], ϖ_Cabannes = [one(FT)],
            ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
            Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
            i_ref = argmin(abs.(ν .- ν̃)), n_Raman = 0,
            F₀ = zeros(FT, 1, 1), SIF₀ = zeros(FT, npol, length(ν)))
        CoreRT.getRamanSSProp!(rs, 1e7 / ν̃, ν)
        rs
    end

    ref = CoreRT.rt_run(make_rrs(), model, 1)
    rs = make_rrs()
    CoreRT.rt_run(rs, model, 1; sources = SolarBeam(F₀ = 2 .* unit))
    again = CoreRT.rt_run(rs, model, 1)
    @test again.toa == ref.toa
    @test again.inelastic_toa == ref.inelastic_toa

    # Thermal emission is not injected by the Raman kernels.
    th = ThermalEmission(B_layer = fill(1.0, _nlay(model), length(ν)))
    @test_throws ArgumentError CoreRT.rt_run(make_rrs(), model, 1;
                                             sources = SolarBeam() + th)

    # Carrier resolution shared by every driver, on a VS object.
    vs = InelasticScattering.VS_0to1_plus{FT}(n2 = n2, o2 = o2)
    CoreRT._prepare_multisensor_F₀!(vs, model, [1], SolarBeam(F₀ = 2 .* unit))
    @test vs.F₀ == 2 .* unit
    CoreRT._prepare_multisensor_F₀!(vs, model, [1], nothing)
    @test vs.F₀ == unit
    vs.F₀ = 3 .* unit
    CoreRT._prepare_multisensor_F₀!(vs, model, [1], nothing)
    @test vs.F₀ == 3 .* unit
end

@testset "HDR/BHR count the direct beam once" begin
    # Vacuum over a Lambertian surface: HDR = μ₀A/π, BHR-down = μ₀,
    # BHR-up = Aμ₀. nstreams=3 Gauss puts μ = 0.5 on a weighted node, so
    # SZA = 60° is on-node; Radau always has μ₀ on a weighted node.
    A = 0.3
    for (quad, sza) in (("GaussLegQuad()", 60.0), ("GaussLegQuad()", acosd(0.5 + 1e-9)),
                        ("RadauQuad()", 30.0))
        model = _qs_model(surface = "LambertianSurfaceScalar($A)", sza = sza,
                          quadrature = quad)
        model.τ_rayl[1] .= 0
        model.τ_abs[1] .= 0
        r = rt_run(model)
        μ₀ = cosd(sza)
        @test r.hdr[:, 1, 1] ≈ fill(μ₀ * A / π, 2) rtol = 1e-12
        @test r.bhr_dw[1] ≈ μ₀ rtol = 1e-12
        @test r.bhr_uw[1] ≈ A * μ₀ rtol = 1e-12
    end

    # Scattering atmosphere: continuous across the node.
    on  = rt_run(_qs_model(sza = 60.0))
    off = rt_run(_qs_model(sza = acosd(0.5 + 1e-9)))
    @test on.hdr ≈ off.hdr rtol = 1e-7
    @test on.bhr_dw ≈ off.bhr_dw rtol = 1e-7
    @test on.bhr_uw ≈ off.bhr_uw rtol = 1e-7

    # rt_run_ss hem_T is the diffuse downwelling flux for every scheme.
    ss_on  = rt_run_ss(_qs_model(sza = 60.0))
    ss_off = rt_run_ss(_qs_model(sza = acosd(0.5 + 1e-9)))
    @test ss_on[6] ≈ ss_off[6] rtol = 1e-7
    @test ss_on[5] ≈ ss_off[5] rtol = 1e-7
    vac = _qs_model(sza = 60.0)
    vac.τ_rayl[1] .= 0
    @test all(iszero, rt_run_ss(vac)[6])
end

@testset "TMS: stream exports rejected, combined bands per band" begin
    model = _qs_model(ss_correction = TMSCorrection())
    @test_throws ArgumentError rt_run_streams(model)
    @test_throws ArgumentError rt_run(model; streams_callback = _ -> nothing)

    # Pure Rayleigh: exact and truncated single scattering coincide, so TMS
    # must leave every band of a combined solve unchanged.
    bands = ["[12987.0]", "[18000.0, 18001.0]"]
    on  = _qs_model(bands = bands, ss_correction = TMSCorrection())
    off = _qs_model(bands = bands)
    R_on  = CoreRT._rt_run_column(InelasticScattering.noRS{Float64}(), on, [1, 2])[1]
    R_off = CoreRT._rt_run_column(InelasticScattering.noRS{Float64}(), off, [1, 2])[1]
    @test R_on ≈ R_off rtol = 1e-10

    # Truncated aerosol: the correction of a combined two-band solve equals
    # that of the separate single-band solves, band by band.
    col(m, b) = CoreRT._rt_run_column(InelasticScattering.noRS{Float64}(), m, b)[1]
    aer = (bands = ["[12987.0]", "[18000.0]"], aerosol = true, sza = 40.0,
           vza = [11.0, 50.0], vaz = [0.0, 120.0])
    on  = _qs_model(; aer..., ss_correction = TMSCorrection())
    off = _qs_model(; aer...)
    ΔR  = col(on, [1, 2]) .- col(off, [1, 2])
    ΔR1 = col(on, [1]) .- col(off, [1])
    ΔR2 = col(on, [2]) .- col(off, [2])
    @test maximum(abs, ΔR) > 1e-3 * maximum(abs, col(off, [1, 2]))   # TMS active
    @test ΔR ≈ cat(ΔR1, ΔR2; dims = 3) rtol = 1e-10
end

@testset "μ₀ / VZA within rounding of a quadrature node" begin
    # A cosine a few ulps from a node must behave as on the node (merged;
    # the node moved onto the exact μ₀) and continuously with a clearly
    # separate ordinate.
    sza_ulp = 60.0
    while cosd(sza_ulp) == 0.5
        sza_ulp = prevfloat(sza_ulp)
    end
    @test cosd(sza_ulp) == nextfloat(0.5)
    for quad in ("GaussLegQuad()", "RadauQuad()")
        runs = map((60.0, sza_ulp, acosd(0.5 + 1e-9))) do sza
            m = _qs_model(sza = sza, quadrature = quad)
            (m, rt_run(m), rt_run_ss(m))
        end
        (m_on, r_on, s_on), (m_ulp, r_ulp, s_ulp), (_, r_off, s_off) = runs
        @test m_ulp.quad_points.Nquad == m_on.quad_points.Nquad
        @test m_ulp.quad_points.iμ₀ == m_on.quad_points.iμ₀
        # μ₀ stays exact and is itself the ordinate at iμ₀ (Gauss moves the
        # coinciding node onto it; Radau builds its nodes around μ₀ and merges
        # the coinciding VZA).
        @test m_ulp.quad_points.μ₀ == cosd(sza_ulp)
        @test Array(m_ulp.quad_points.qp_μ)[m_ulp.quad_points.iμ₀] == cosd(sza_ulp)
        for (a, b) in ((r_on, r_ulp), (r_on, r_off))
            @test a.toa ≈ b.toa rtol = 1e-7
            @test a.hdr ≈ b.hdr rtol = 1e-7
        end
        @test s_ulp[6] ≈ s_on[6] rtol = 1e-7
        @test s_off[6] ≈ s_on[6] rtol = 1e-7
    end

    # Clustered ordinates: μ₀ and a VZA a few ulps apart near a node. No two
    # distinct ordinates of the final grid may be numerically coincident.
    e = eps(1.0)
    clustered = _qs_model(sza = acosd(0.5 + 3e), vza = [acosd(0.5 + 6e)], vaz = [180.0])
    μ_grid = sort(Array(clustered.quad_points.qp_μ))
    @test minimum(diff(μ_grid)) > 4e
    @test rt_run(clustered).toa ≈
          rt_run(_qs_model(sza = 60.0, vza = [60.0], vaz = [180.0])).toa rtol = 1e-7

    # A VZA whose cosine is a rounding error away from a Gauss node.
    node = Array(_qs_model(sza = 30.0).quad_points.qp_μ)[1]
    vza_near = acosd(node)
    @test cosd(vza_near) != node
    near = _qs_model(sza = 30.0, vza = [vza_near])
    far  = _qs_model(sza = 30.0, vza = [acosd(node + 1e-9)])
    @test near.quad_points.Nquad == far.quad_points.Nquad - 1   # merged
    @test rt_run(near).toa ≈ rt_run(far).toa rtol = 1e-7
end

@testset "multisensor TOA keeps the Cox–Munk glint residual" begin
    model = _qs_model(surface = "CoxMunkSurface(wind_speed=2.0)", pol = "Stokes_IQU()",
                      nstreams = 6, sza = 30.0, vza = [30.0, 30.0, 25.0],
                      vaz = [0.0, 180.0, 180.0], obs_alt = [0, 5.0])
    ms = rt_run(model)
    column = CoreRT._rt_run_column(InelasticScattering.noRS{Float64}(), model, [1])
    @test ms.toa ≈ column[1] rtol = 1e-12

    # Float32 model with Float64 aerosol inputs: the optics (and τ) are
    # promoted, the glint correction runs in the radiance precision.
    p32 = _qs_params(float_type = "Float32", aerosol = true,
                     surface = "CoxMunkSurface(wind_speed=2.0)", pol = "Stokes_IQU()",
                     sza = 30.0, vza = [30.0, 25.0], vaz = [0.0, 180.0], obs_alt = [0, 5.0])
    p32.scattering_params.rt_aerosols[1].aerosol =
        vSmartMOM.Scattering.Aerosol(LogNormal(log(0.6), log(1.6)), 1.45, 0.001)
    model32 = model_from_parameters(p32)
    ms32 = rt_run(model32)
    column32 = CoreRT._rt_run_column(InelasticScattering.noRS{Float32}(), model32, [1])
    @test eltype(ms32.toa) === Float32
    @test ms32.toa ≈ column32[1] rtol = 1e-5
end

@testset "U/V-polarized incident beams are rejected" begin
    model = _qs_model(pol = "Stokes_IQU()")
    nspec = _nspec(model)
    FU = zeros(3, nspec); FU[1, :] .= 1; FU[3, :] .= 0.1
    FQ = zeros(3, nspec); FQ[1, :] .= 1; FQ[2, :] .= 0.1
    @test_throws ArgumentError rt_run(model; sources = SolarBeam(F₀ = FU))
    @test_throws ArgumentError rt_run_ss(model; sources = SolarBeam(F₀ = FU))
    @test all(isfinite, rt_run(model; sources = SolarBeam(F₀ = FQ)).toa)
    rs = InelasticScattering.noRS{Float64}()
    rs.F₀ = copy(FU)
    @test_throws ArgumentError CoreRT.rt_run(rs, model, 1)
    model_ms = _qs_model(pol = "Stokes_IQU()", obs_alt = [0, 5.0])
    @test_throws ArgumentError rt_run(model_ms; sources = SolarBeam(F₀ = FU))
    fwd, lin = model_from_parameters(LinMode(), _qs_params(pol = "Stokes_IQU()"))
    @test_throws ArgumentError rt_run(fwd, lin, 0, size(lin.τ̇_abs[1], 1), 1;
                                      sources = SolarBeam(F₀ = FU))
    model_v = _qs_model(pol = "Stokes_IQUV()")
    FV = zeros(4, nspec); FV[1, :] .= 1; FV[4, :] .= 0.1
    @test_throws ArgumentError rt_run(model_v; sources = SolarBeam(F₀ = FV))
end

@testset "ThermalEmission: unsupported paths and layer count" begin
    model = _qs_model()
    nspec, nlay = _nspec(model), _nlay(model)
    th = ThermalEmission(B_layer = fill(1.0, nlay, nspec))
    @test_throws ArgumentError rt_run_ss(model; sources = th)
    @test all(isfinite, rt_run_ss(model; sources = ThermalEmission())[1])

    @test_throws ArgumentError rt_run(model;
        sources = ThermalEmission(B_layer = fill(1.0, nlay + 1, nspec)))
    @test_throws ArgumentError rt_run(model;
        sources = SolarBeam() + ThermalEmission(B_layer = fill(1.0, 1, nspec)))
    @test rt_run(model; sources = ThermalEmission(B_layer = zeros(1, nspec))).toa ==
          rt_run(model; sources = NoSource()).toa
    @test all(isfinite, rt_run(model; sources = th).toa)
end

@testset "linearized RT rejects surfaces without a linearized layer" begin
    for (surface, name) in (("rpvSurfaceScalar(0.1, 0.5, 1.0, -0.2)", "rpvSurfaceScalar"),
                            ("LambertianSurfaceSpectrum([0.2])", "LambertianSurfaceSpectrum"))
        fwd, lin = model_from_parameters(LinMode(), _qs_params(surface = surface))
        err = try
            rt_run(fwd, lin, 0, size(lin.τ̇_abs[1], 1), 1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(name, sprint(showerror, err))
    end
    @test supports_linearized_surface(CoreRT.LambertianSurfaceScalar(0.1))
    @test supports_linearized_surface(CoreRT.CoxMunkSurface(wind_speed = 2.0))
    @test !supports_linearized_surface(CoreRT.RossLiSurfaceScalar(0.1, 0.1, 0.1))
end

end # module RTEntryGuardsTests
