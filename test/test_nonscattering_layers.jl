# =========================================================================
# Non-scattering layers (τ·ϖ ≤ 2eps) inside a scattering column.
#
# A layer whose scattering optical depth vanishes must be optically
# indistinguishable from (a) the same layer with a vanishingly small
# scattering optical depth and (b) the column with that layer merged away.
# These checks pin three defects of the non-scattering shortcut:
#
#   * stale downward sources: the reused added layer kept `j₀⁺` and the
#     per-source slots from the previous layer / Fourier moment;
#   * surface interface: the surface was added with the LAST atmospheric
#     layer's scattering-interface trait instead of (column scatters,
#     surface reflects), dropping surface–atmosphere coupling;
#   * thermal emission: a non-scattering layer emitted nothing (and leaked
#     stale thermal slots), so a pure absorber gave TOA = 0 instead of
#     B(1 − e^{−τ/μ}).
#
# The column (`rt_run`), multisensor (interior `obs_alt`) and linearized
# paths are all exercised.
# =========================================================================

using Test
using Logging
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering: Stokes_I, Stokes_IQU
using vSmartMOM.SolarModel: planck_spectrum_wn

const _NS_YAML = joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml")

_ns_quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Warn))

function _ns_params(n; FT=Float64, pol=Stokes_IQU(), albedo=0.0, obs_alt=[0.0],
                    spec_bands=nothing)
    p = parameters_from_yaml(_NS_YAML)
    p.architecture = vSmartMOM.Architectures.CPU()
    p.float_type = FT
    p.polarization_type = pol
    spec_bands === nothing || (p.spec_bands = spec_bands)
    p.sza = 30.0
    p.vza = [0.0, 30.0, 60.0]
    p.vaz = [0.0, 90.0, 180.0]
    p.T = fill(250.0, n)
    p.q = zeros(n)
    p.p = collect(range(100.0, 1000.0, length=n + 1))
    p.brdf = [LambertianSurfaceScalar(albedo)]
    p.obs_alt = obs_alt
    return p
end

# Overwrite the layer optics with a prescribed Rayleigh / absorption column.
function _ns_set_optics!(model, τ_rayl, τ_abs)
    FT = CoreRT.float_type(model)
    for iz in eachindex(τ_rayl)
        model.τ_rayl[1][:, iz] .= FT(τ_rayl[iz])
        model.τ_abs[1][:, iz] .= FT(τ_abs[iz])
    end
    return model
end

function _ns_model(τ_rayl; τ_abs=zero(τ_rayl), kwargs...)
    p = _ns_params(length(τ_rayl); kwargs...)
    return _ns_set_optics!(model_from_parameters(p; external_solar=false), τ_rayl, τ_abs)
end

_ns_run(τ_rayl; sources=nothing, kwargs...) = _ns_quiet() do
    model = _ns_model(τ_rayl; kwargs...)
    sources === nothing ? rt_run(model) : rt_run(model; sources=sources(model))
end

_ns_rel(a, b) = maximum(abs.(a .- b)) / maximum(abs.(b))

@testset "Non-scattering layers" begin

    # ---------------------------------------------------------------------
    @testset "column: τ_layer = 0 ≡ 1e-12 ≡ merged (albedo $albedo)" for albedo in (0.0, 0.3)
        cases = (
            # (with zero layer, with 1e-12 layer, merged column)
            ([0.1, 0.0, 0.1], [0.1, 1e-12, 0.1], [0.1, 0.1]),   # interior
            ([0.1, 0.0],      [0.1, 1e-12],      [0.1]),        # bottom layer
            ([0.0, 0.1],      [1e-12, 0.1],      [0.1]),        # top layer
            ([0.1, 0.0, 0.0], [0.1, 1e-12, 1e-12], [0.1]),      # two at the bottom
        )
        for (τ0, τε, τm) in cases
            r0 = _ns_run(τ0; albedo)
            rε = _ns_run(τε; albedo)
            rm = _ns_run(τm; albedo)
            # TOA upwelling and BOA downwelling, all three Stokes components.
            @test r0.toa ≈ rε.toa rtol=1e-9
            @test r0.boa ≈ rε.boa rtol=1e-9
            @test r0.toa ≈ rm.toa rtol=1e-12
            @test r0.boa ≈ rm.boa rtol=1e-12
            # Polarization survives (Q is non-trivial for Rayleigh).
            @test maximum(abs, r0.toa[:, 2, :]) > 1e-3
        end
        # A purely absorbing interior layer (no merge partner).
        τ_abs = [0.0, 0.05, 0.0]
        r0 = _ns_run([0.1, 0.0, 0.1]; τ_abs, albedo)
        rε = _ns_run([0.1, 1e-12, 0.1]; τ_abs, albedo)
        @test r0.toa ≈ rε.toa rtol=1e-9
        @test r0.boa ≈ rε.boa rtol=1e-9
    end

    @testset "column: fully non-scattering atmosphere over a Lambertian surface" begin
        # Every layer (and therefore the composite) is non-scattering: the
        # surface must still be added as a reflecting layer.
        r0 = _ns_run([0.0, 0.0]; τ_abs=[0.05, 0.05], albedo=0.3)
        rε = _ns_run([1e-12, 1e-12]; τ_abs=[0.05, 0.05], albedo=0.3)
        @test r0.toa ≈ rε.toa rtol=1e-9
        @test r0.boa ≈ rε.boa rtol=1e-9
        @test minimum(r0.toa[:, 1, :]) > 0
    end

    # ---------------------------------------------------------------------
    @testset "multisensor: interior observers around a zero layer (albedo $albedo)" for albedo in (0.0, 0.3)
        # Heights of the existing interfaces, so the requested observers
        # snap to them and the layer count is unchanged.
        z = _ns_quiet() do
            CoreRT.half_level_altitudes(_ns_model([0.1, 0.0, 0.1, 0.05]).profile)
        end
        # Observers directly above and directly below the empty layer 2.
        obs = [0.0, z[2], z[3]]
        ms(τ) = _ns_quiet() do
            rt_run(_ns_model(τ; albedo, obs_alt=obs))
        end
        r0 = ms([0.1, 0.0, 0.1, 0.05])
        rε = ms([0.1, 1e-12, 0.1, 0.05])
        @test length(r0.levels) == 2
        @test [l.boundary_index for l in r0.levels] == [1, 2]
        @test r0.toa ≈ rε.toa rtol=1e-9
        @test r0.boa ≈ rε.boa rtol=1e-9
        for (l0, lε) in zip(r0.levels, rε.levels)
            @test l0.upwelling ≈ lε.upwelling rtol=1e-9
            @test l0.downwelling ≈ lε.downwelling rtol=1e-9
        end
        # The empty layer is invisible: both of its faces see the same field.
        @test r0.levels[1].upwelling ≈ r0.levels[2].upwelling rtol=1e-12
        @test r0.levels[1].downwelling ≈ r0.levels[2].downwelling rtol=1e-12
        # Multisensor endpoints agree with the column solver and the merged column.
        rc = _ns_run([0.1, 0.0, 0.1, 0.05]; albedo)
        rm = _ns_run([0.1, 0.1, 0.05]; albedo)
        @test r0.toa ≈ rc.toa rtol=1e-10
        @test r0.boa ≈ rc.boa rtol=1e-10
        @test r0.toa ≈ rm.toa rtol=1e-10
        @test r0.boa ≈ rm.boa rtol=1e-10

        # Non-scattering bottom layer below an interior observer (surface step).
        zb = _ns_quiet() do
            CoreRT.half_level_altitudes(_ns_model([0.1, 0.1, 0.0]).profile)
        end
        msb(τ) = _ns_quiet() do
            rt_run(_ns_model(τ; albedo, obs_alt=[0.0, zb[2]]))
        end
        b0 = msb([0.1, 0.1, 0.0])
        bε = msb([0.1, 0.1, 1e-12])
        @test b0.toa ≈ bε.toa rtol=1e-9
        @test b0.boa ≈ bε.boa rtol=1e-9
        @test only(b0.levels).upwelling ≈ only(bε.levels).upwelling rtol=1e-9
        @test only(b0.levels).downwelling ≈ only(bε.levels).downwelling rtol=1e-9
        @test b0.toa ≈ _ns_run([0.1, 0.1]; albedo).toa rtol=1e-10
    end

    # ---------------------------------------------------------------------
    @testset "thermal: isothermal pure absorber B(1 − e^{−τ/μ})" begin
        T_iso = 250.0
        thermal(model) = ThermalEmission(fill(T_iso, length(model.profile.T)),
                                         collect(CoreRT.get_spec_bands(model)[1]))
        for τ_abs in ([0.2, 0.3, 0.1], [40.0, 40.0, 40.0])
            r0 = _ns_run([0.0, 0.0, 0.0]; τ_abs, pol=Stokes_I(), sources=thermal)
            rε = _ns_run([1e-12, 1e-12, 1e-12]; τ_abs, pol=Stokes_I(), sources=thermal)
            ν = collect(Float64, _ns_quiet(() -> CoreRT.get_spec_bands(
                _ns_model([0.0, 0.0, 0.0]))[1]))
            B = planck_spectrum_wn(T_iso, ν)
            μ = cosd.([0.0, 30.0, 60.0])
            expected = [B[n] * (1 - exp(-sum(τ_abs) / μ[i]))
                        for i in eachindex(μ), n in eachindex(B)]
            # TOA upwelling and BOA downwelling of a black-surface slab.
            @test r0.toa[:, 1, :] ≈ expected rtol=1e-10
            @test r0.boa[:, 1, :] ≈ expected rtol=1e-10
            @test rε.toa ≈ r0.toa rtol=1e-9
        end
        # Over a reflecting surface the downwelling emission is reflected:
        # TOA(μ) = B(1 − e^{−τ/μ}) + e^{−τ/μ}·a·B·(1 − 2E₃(τ)), with
        # 1 − 2E₃(τ) = 2∫₀¹ μ(1 − e^{−τ/μ})dμ (no surface emission here).
        a, τ_abs = 0.3, [0.4, 0.4, 0.2]
        τ_tot = sum(τ_abs)
        r0 = _ns_run([0.0, 0.0, 0.0]; τ_abs, albedo=a, pol=Stokes_I(), sources=thermal)
        ν = collect(Float64, _ns_quiet(() -> CoreRT.get_spec_bands(
            _ns_model([0.0, 0.0, 0.0]))[1]))
        B = planck_spectrum_wn(T_iso, ν)
        x, w = vSmartMOM.Scattering.gauleg(400, 0.0, 1.0)
        flux = 2 * sum(w .* x .* (1 .- exp.(-τ_tot ./ x)))
        μ = cosd.([0.0, 30.0, 60.0])
        expected = [B[n] * ((1 - exp(-τ_tot / μ[i])) + exp(-τ_tot / μ[i]) * a * flux)
                    for i in eachindex(μ), n in eachindex(B)]
        @test r0.toa[:, 1, :] ≈ expected rtol=1e-6
        @test r0.boa[:, 1, :] ≈ [B[n] * (1 - exp(-τ_tot / μ[i]))
                                 for i in eachindex(μ), n in eachindex(B)] rtol=1e-10
        # A non-scattering layer between scattering layers still emits.
        τ_abs = [0.1, 0.5, 0.1]
        r0 = _ns_run([0.05, 0.0, 0.05]; τ_abs, pol=Stokes_I(), sources=thermal)
        rε = _ns_run([0.05, 1e-12, 0.05]; τ_abs, pol=Stokes_I(), sources=thermal)
        @test r0.toa ≈ rε.toa rtol=1e-9
        @test r0.boa ≈ rε.boa rtol=1e-9
    end

    @testset "thermal: Float32 ≈ Float64 in the thermal IR" begin
        # At 1000 cm⁻¹ the per-layer Rayleigh optical depth falls below
        # 2eps(Float32) in most layers, so Float32 routes them through the
        # non-scattering shortcut while Float64 still scatters.
        function ir(FT)
            _ns_quiet() do
                p = parameters_from_yaml(_NS_YAML)
                p.architecture = vSmartMOM.Architectures.CPU()
                p.float_type = FT
                p.spec_bands = [[1000.0, 1001.0]]
                p.vza = [0.0, 45.0]
                p.vaz = [0.0, 0.0]
                p.polarization_type = Stokes_I()
                model = model_from_parameters(p)
                fill!(model.τ_abs[1], FT(0.05))
                ν = collect(CoreRT.get_spec_bands(model)[1])
                nsc = count(iz -> maximum(model.τ_rayl[1][:, iz]) <= 2eps(FT),
                            axes(model.τ_rayl[1], 2))
                r = rt_run(model; sources=ThermalEmission(FT.(model.profile.T), ν))
                (Float64.(r.toa), Float64.(r.boa), nsc)
            end
        end
        R64, T64, nsc64 = ir(Float64)
        R32, T32, nsc32 = ir(Float32)
        @test nsc64 == 0
        @test nsc32 > 0
        @test R32 ≈ R64 rtol=1e-4
        @test T32 ≈ T64 rtol=1e-4
    end

    # ---------------------------------------------------------------------
    @testset "linearized: forward parity and surface Jacobian with a zero layer" begin
        albedo = 0.3
        # Purely absorbing (non-scattering) interior and bottom layers.
        τ_rayl = [0.1, 0.0, 0.1, 0.0]
        τ_abs = [0.0, 0.05, 0.0, 0.05]
        lin(τ) = _ns_quiet() do
            model, lin_model = model_from_parameters(
                LinMode(), _ns_params(length(τ); albedo))
            _ns_set_optics!(model, τ, τ_abs)
            rt_run(model, lin_model, 0, size(lin_model.τ̇_abs[1], 1), 1)
        end
        res = lin(τ_rayl)
        fwd = _ns_run(τ_rayl; τ_abs, albedo)
        @test res.toa ≈ fwd.toa rtol=1e-10
        @test res.boa ≈ fwd.boa rtol=1e-10
        fwdε = _ns_run([0.1, 1e-12, 0.1, 1e-12]; τ_abs, albedo)
        @test res.toa ≈ fwdε.toa rtol=1e-9
        @test res.boa ≈ fwdε.boa rtol=1e-9
        col = first(CoreRT.surface_range(res.layout))
        δ = 1e-5
        plus = _ns_run(τ_rayl; τ_abs, albedo=albedo + δ)
        minus = _ns_run(τ_rayl; τ_abs, albedo=albedo - δ)
        @test res.toa_jacobian[:, :, :, col] ≈ (plus.toa .- minus.toa) ./ 2δ rtol=1e-6
        @test res.boa_jacobian[:, :, :, col] ≈ (plus.boa .- minus.boa) ./ 2δ rtol=1e-6
    end
end
