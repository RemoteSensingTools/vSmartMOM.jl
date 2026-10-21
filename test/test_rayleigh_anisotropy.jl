# Automatic (`depol: -1`) Rayleigh depolarization and the tensor dispersion
# (v2.2 decision 13c). `params.rayleigh_anisotropy` selects the N₂/O₂
# polarizability tensor behind greek_rayleigh and τ_rayl:
# DispersiveAnisotropy() (default) or ConstantAnisotropy() (the Raman-module
# constants, the v2.2 pre-release behaviour). Cabannes quantities always use the Raman-module
# constants.

using Test
using vSmartMOM

const _RA_IS = vSmartMOM.InelasticScattering
const _RA_CRT = vSmartMOM.CoreRT

_ra_params(file = "PureRayleighParameters.yaml") =
    parameters_from_yaml(joinpath("test_parameters", file))

function _ra_dict(; extra...)
    rt = Dict{String,Any}(
        "spec_bands" => ["[13000.0 13000.1]", "[6200.0 6200.1]"],
        "surface" => ["LambertianSurfaceScalar(0.1)", "LambertianSurfaceScalar(0.1)"],
        "polarization_type" => "Stokes_IQU()", "nstreams" => 3,
        "truncation" => "NoTruncation()", "depol" => -1,
        "float_type" => "Float64", "architecture" => "CPU()")
    for (k, v) in extra
        rt[string(k)] = v
    end
    return Dict{String,Any}(
        "radiative_transfer" => rt,
        "geometry" => Dict("sza" => 30.0, "vza" => [0.0], "vaz" => [0.0],
                           "obs_alt" => [0]),
        "atmospheric_profile" => Dict(
            "T" => [220.0, 240.0, 270.0, 288.0],
            "p" => [1.0, 100.0, 400.0, 700.0, 1013.0],
            "profile_reduction" => -1))
end

# The pre-release rule, written out: Rayleigh depolarization from the
# Raman-module constants at the band centre and column-weighted temperature.
function _ra_constant_depol(λₘ, T)
    n2, o2 = _RA_IS.getRamanAtmoConstants(oftype(λₘ, 1.0e7) / λₘ, T)
    γ, _ = _RA_IS.compute_γ_air_Rayleigh!(λₘ, n2, o2)
    return 2γ / (1 + γ)
end

@testset "rayleigh_anisotropy parsing" begin
    @test parameters_from_dict(_ra_dict()).rayleigh_anisotropy === DispersiveAnisotropy()
    @test _ra_params().rayleigh_anisotropy === DispersiveAnisotropy()
    for spec in ("ConstantAnisotropy()", "ConstantAnisotropy")
        @test parameters_from_dict(_ra_dict(rayleigh_anisotropy = spec)).rayleigh_anisotropy ===
              ConstantAnisotropy()
    end
    @test parameters_from_dict(_ra_dict(rayleigh_anisotropy = "DispersiveAnisotropy()")).rayleigh_anisotropy ===
          DispersiveAnisotropy()
    @test_throws Exception parameters_from_dict(_ra_dict(rayleigh_anisotropy = "Dispersive()"))

    # Positional construction with the pre-release field list still works and
    # selects the default.
    p = _ra_params()
    old_fields = [getfield(p, f) for f in fieldnames(typeof(p))[1:end-1]]
    for q in (vSmartMOM.CoreRT.vSmartMOM_Parameters(old_fields...),
              typeof(p)(old_fields...))
        @test q isa typeof(p)
        @test q.rayleigh_anisotropy === DispersiveAnisotropy()
        @test q.depol == p.depol && q.spec_bands == p.spec_bands
    end
    p32 = parameters_from_dict(_ra_dict(float_type = "Float32"))
    q32 = typeof(p32)([getfield(p32, f) for f in fieldnames(typeof(p32))[1:end-1]]...)
    @test q32 isa vSmartMOM.CoreRT.vSmartMOM_Parameters{Float32}
    @test q32.rayleigh_anisotropy === DispersiveAnisotropy()
    @test_throws MethodError vSmartMOM.CoreRT.vSmartMOM_Parameters(old_fields[1:end-1]...)
end

@testset "depol: -1 Rayleigh depolarization follows rayleigh_anisotropy" begin
    # Float64; Float32 parsed natively (Float32 spectral grids); and
    # Float64-parsed parameters switched to Float32 (mixed-precision grids).
    for (FT, base) in ((Float64, parameters_from_dict(_ra_dict())),
                       (Float32, parameters_from_dict(_ra_dict(float_type = "Float32"))),
                       (Float32, let b = parameters_from_dict(_ra_dict())
                                     b.float_type = Float32; b end))
        models = Dict{Any,Any}()
        for a in (DispersiveAnisotropy(), ConstantAnisotropy())
            params = deepcopy(base)
            params.rayleigh_anisotropy = a
            model = model_from_parameters(params)
            models[a] = model
            profile = model.atmosphere.profile
            T = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry)
            for (i, band) in enumerate(params.spec_bands)
                νₘ = FT(0.5) * (band[1] + band[end])
                λₘ = FT(1.0e7) / νₘ
                depol = a isa ConstantAnisotropy ? _ra_constant_depol(λₘ, FT(T)) :
                        first(_RA_IS.rayleigh_air_depolarization(a, λₘ, FT(T)))
                greek = vSmartMOM.Scattering.get_greek_rayleigh(FT(depol))
                @test model.optics.rayleigh.greek_rayleigh[i].β == greek.β
                @test model.optics.rayleigh.greek_rayleigh[i].δ == greek.δ
                τ = _RA_CRT.getRayleighLayerOptProp(profile.p_half[end],
                        FT.(1e4 ./ band), FT(depol), profile.vcd_dry)
                @test model.optics.τ_rayl[i] == τ
            end
        end
        d, c = models[DispersiveAnisotropy()], models[ConstantAnisotropy()]
        # Cabannes quantities stay on the Raman-module constants.
        @test d.optics.rayleigh.ϖ_Cabannes == c.optics.rayleigh.ϖ_Cabannes
        @test all(x.β == y.β for (x, y) in zip(d.optics.rayleigh.greek_cabannes,
                                               c.optics.rayleigh.greek_cabannes))
        # The dispersive anisotropy lowers ρ and τ_rayl in the NIR
        # (769 nm: about -6 % and -0.3 %; 1613 nm: about -9 % and -0.4 %).
        for i in 1:2
            ratio = sum(d.optics.τ_rayl[i]) / sum(c.optics.τ_rayl[i])
            @test FT(0.99) < ratio < 1
        end
        @test d.optics.rayleigh.greek_rayleigh[2].β[3] > c.optics.rayleigh.greek_rayleigh[2].β[3]
    end

    # An explicit depol ≥ 0 does not depend on the anisotropy.
    params = parameters_from_dict(_ra_dict(depol = 0.03))
    mc = deepcopy(params); mc.rayleigh_anisotropy = ConstantAnisotropy()
    m1, m2 = model_from_parameters(params), model_from_parameters(mc)
    @test m1.optics.τ_rayl == m2.optics.τ_rayl
    @test m1.optics.rayleigh.greek_rayleigh[1].β == m2.optics.rayleigh.greek_rayleigh[1].β
end

@testset "linearized and batch builders use the same anisotropy" begin
    for a in (DispersiveAnisotropy(), ConstantAnisotropy())
        p32 = parameters_from_dict(_ra_dict(float_type = "Float32"))
        p32.rayleigh_anisotropy = a
        fwd32 = model_from_parameters(p32)
        model32, _ = model_from_parameters(LinMode(), p32)
        @test eltype(model32.optics.τ_rayl[2]) === Float32
        @test model32.optics.τ_rayl == fwd32.optics.τ_rayl
        @test all(x.β == y.β for (x, y) in zip(model32.optics.rayleigh.greek_rayleigh,
                                               fwd32.optics.rayleigh.greek_rayleigh))

        params = _ra_params("JacobianTestRayleigh.yaml")
        params.rayleigh_anisotropy = a
        fwd = model_from_parameters(params)
        model, lin = model_from_parameters(LinMode(), params)
        @test model.optics.τ_rayl == fwd.optics.τ_rayl
        @test model.optics.rayleigh.greek_rayleigh[1].β == fwd.optics.rayleigh.greek_rayleigh[1].β

        ctx = BatchContext(params)
        newT = params.T .+ 5.0
        update_model!(ctx; T = newT)
        fresh_params = deepcopy(params); fresh_params.T = newT
        fresh = model_from_parameters(fresh_params)
        @test ctx.model.optics.τ_rayl[1] == fresh.optics.τ_rayl[1]
        @test ctx.model.optics.rayleigh.greek_rayleigh[1].β ==
              fresh.optics.rayleigh.greek_rayleigh[1].β
    end
end
