using Test
using LinearAlgebra
using vSmartMOM
using vSmartMOM.CoreRT
import Interpolations
import CUDA
import AtmosphericAbsorption

# Synthetic spectroscopy keeps these pressure regressions independent of
# downloaded HITRAN/ABSCO data. The line centre and both wings respond to p.
function _pressure_test_line_model(::Type{FT}; water=false,
        profile=AtmosphericAbsorption.Voigt(), mixing=0,
        cpf=AtmosphericAbsorption.HumlicekWeideman32(), advanced=false,
        architecture=AtmosphericAbsorption.CPU()) where {FT}
    lines = AtmosphericAbsorption.LineDatabase(
        mol=Int32[water ? 1 : 2], iso=Int32[1], ν0=FT[6200],
        S=FT[water ? 1e-25 : 1e-23], E_lower=FT[0], g_upper=FT[1],
        γ_air=FT[0.07], γ_self=FT[0.30], n_air=FT[0.7],
        # At 6200 cm⁻¹, a Float32 centre cannot resolve the tiny shift from
        # a 0.5-hPa FD step. Isolate broadening in Float32; Float64 and the
        # resolved Doppler tests below exercise pressure shifts separately.
        δ_air=FT[FT === Float32 ? 0 : -0.01], n_self=FT[0.5],
        δ_self=FT[FT === Float32 ? 0 : 0.04],
        Y_LM=FT[mixing], n_Y_LM=FT[0.5],
        γ2_air=FT[advanced ? 0.005 : 0], δ2_air=FT[advanced ? 0.001 : 0],
        νVC=FT[advanced ? 0.01 : 0], η=FT[advanced ? 0.2 : 0],
        molar_mass=FT[water ? 18 : 44],
        meta=AtmosphericAbsorption.SourceMetadata("synthetic pressure test", 296.0, 1013.25))
    partition = AtmosphericAbsorption.TabulatedPF(FT[100, 300, 500], FT[1, 1, 1])
    AtmosphericAbsorption.LineByLineModel(lines, partition;
        profile, cpf, wing_cutoff=FT(10), vmr=zero(FT), architecture)
end

function _pressure_test_absco(::Type{FT}; architecture=AtmosphericAbsorption.CPU()) where {FT}
    ν = FT[6199.8, 6200, 6200.2]
    p = FT[100, 600, 1100]
    # Temperature coordinates slide with pressure. A slope of raw table
    # entries would be wrong: first interpolate at the fixed physical T.
    T = FT[210 230 250; 290 310 330]
    broadener = FT[0, 0.05]
    σ = zeros(FT, length(ν), length(broadener), size(T, 1), length(p))
    for ip in eachindex(p), it in axes(T, 1), ib in eachindex(broadener), iν in eachindex(ν)
        σ[iν, ib, it, ip] = FT(1e-23) * FT(iν) *
            (1 + FT(0.001) * p[ip] + FT(0.000001) * p[ip]^2 +
             FT(0.004) * T[it, ip] + 2 * broadener[ib])
    end
    AtmosphericAbsorption.AbscoLUT(2, -1, ν, p, T, broadener, σ;
        architecture)
end

function _pressure_test_regular_lut(::Type{FT}) where {FT}
    ν, p, T = FT[6199.8, 6200, 6200.2], FT[100, 600, 1100], FT[210, 290, 330]
    σ = FT[FT(1e-23) * FT(i) *
        (1 + FT(0.001) * pp + FT(0.000001) * pp^2 + FT(0.004) * tt)
        for i in eachindex(ν), pp in p, tt in T]
    AtmosphericAbsorption.InterpolationModel(ν, p, T, σ, zero(FT), AtmosphericAbsorption.CPU())
end

@testset "Line profiles, CPF and line-mixing pressure response" begin
    rational = AtmosphericAbsorption.HumlicekWeideman32()
    cases = ((AtmosphericAbsorption.Doppler(), rational, false),
             (AtmosphericAbsorption.Lorentz(), rational, false),
             (AtmosphericAbsorption.Voigt(), rational, false),
             (AtmosphericAbsorption.Voigt(), AtmosphericAbsorption.ErfcxCPF(), false),
             (AtmosphericAbsorption.HartmannTran(), rational, true))
    for (profile, cpf, advanced) in cases
        model = _pressure_test_line_model(Float64; profile, cpf, advanced, mixing=0.15)
        # Resolve the narrow Doppler core; pressure shifts it even though its
        # Doppler width itself is independent of pressure at fixed T.
        grid = collect(6199.98:0.001:6200.02)
        p, T, h = 750.0, 280.0, 0.05
        σ(p) = collect(vSmartMOM.Absorption.absorption_cross_section(model, grid, p, T))
        fd = (σ(p + h) - σ(p - h)) / (2h)
        actual = CoreRT._absorption_pressure_derivative(model, grid, p, T, nothing)
        @test norm(fd) > 0
        @test actual ≈ fd rtol=2e-6
        @test all(iszero, CoreRT._absorption_pressure_derivative(
            model, [6100.0, 6300.0], p, T, nothing))
    end
end

function _pressure_test_parameters(::Type{FT}, absorber; species=:variable) where {FT}
    params = read_parameters(Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[6199.8, 6199.95, 6200.0, 6200.05, 6200.2]"],
            "surface" => ["LambertianSurfaceScalar(0.2)"],
            "nstreams" => 3, "polarization_type" => "Stokes_I()",
            "truncation" => "NoTruncation()", "depol" => -1,
            "float_type" => string(FT), "architecture" => "CPU()"),
        "geometry" => Dict("sza" => 35.0, "vza" => [25.0], "vaz" => [40.0], "obs_alt" => [0]),
        "atmospheric_profile" => Dict("T" => [270.0, 280.0],
            "p" => [100.0, 500.0, 1000.0], "profile_reduction" => -1)))
    fixed = species in (:fixed, :mixed) ? ["CO2"] : String[]
    variable = species in (:variable, :mixed) ? ["CH4"] : String[]
    vmr = Dict{String,Any}("CO2" => FT(4e-4), "CH4" => FT[3e-4, 5e-4])
    luts = [Any[absorber for _ in vcat(fixed, variable)]]
    h2o = species === :water ? Any[absorber] : Any[:disabled]
    params.absorption_params = CoreRT.AbsorptionParameters(
        [fixed], [variable], vmr, AtmosphericAbsorption.Voigt(),
        AtmosphericAbsorption.HumlicekWeideman32(), FT(10), luts, h2o, String[], "")
    params.q .= species === :water ? FT(0.01) : zero(FT)
    return params
end

# A downstream absorber can implement forward cross sections without yet
# supporting pressure differentiation. Dispatch must preserve that extension.
struct _ForwardOnlyPressureTestAbsorber{M}
    model::M
end
vSmartMOM.Absorption.absorption_cross_section(a::_ForwardOnlyPressureTestAbsorber,
    grid, p, T) = CoreRT._layer_absorption_cross_section(a.model, grid, p, T, nothing)

struct _FixedPressureTestFlavor <: AbstractJacobianFlavor end
CoreRT.requires_pressure_jacobians(::_FixedPressureTestFlavor) = false
function CoreRT.jacobian_plan(flavor::_FixedPressureTestFlavor, params, model, lin)
    keys = [ParameterKey(:gas, :vmr; component="CH4", layer=1),
            ParameterKey(:surface, :P0; component=1, band=1)]
    names = ["ch4_layer1", "albedo"]
    # Native pressure=1, q-H2O=2:3, variable CH4=4:5.
    layout = ActiveParameterLayout(keys, names, keys, names, [4], 1;
                                   surface_columns=2:2)
    JacobianPlan(flavor, keys, names, [layout])
end

@testset "Fixed pressure allows forward-only absorber extensions" begin
    absorber = _pressure_test_regular_lut(Float64)
    params = _pressure_test_parameters(Float64, _ForwardOnlyPressureTestAbsorber(absorber))
    @test_throws r"pressure derivatives are not implemented" model_from_parameters(LinMode(), params)
    model, planned = model_from_parameters(_FixedPressureTestFlavor(), params)
    @test planned.τ̇_abs_psurf === nothing
    @test planned.τ̇_rayl_psurf === nothing
    @test planned.τ̇_aer_psurf === nothing
    # Forward and linearized column integration differ in multiplication
    # order, so allow a few ulps here; selected RT/tangents below remain exact.
    @test only(model.τ_abs) ≈ only(model_from_parameters(params).τ_abs) rtol=4eps(Float64)
    # The raw full layout requests pressure and must not expose a false zero.
    @test_throws r"Surface-pressure derivatives are unavailable" rt_run(
        model, planned.base, 0, size(planned.τ̇_abs[1], 1), 1)
    supported = _pressure_test_parameters(Float64, absorber)
    reference, reference_lin = model_from_parameters(_FixedPressureTestFlavor(), supported;
        compute_pressure_jacobians=true)
    for basis in (:local, :physical)
        result = rt_run(model, planned; jacobian_basis=basis)
        expected = rt_run(reference, reference_lin; jacobian_basis=basis)
        @test result.toa == expected.toa
        @test result.toa_jacobian == expected.toa_jacobian
    end
end

@testset "Surface-pressure absorption full rebuild" begin
    for FT in (Float32, Float64), kind in (:line, :absco, :regular), species in (:fixed, :variable, :mixed, :water)
        absorber = kind === :line ? _pressure_test_line_model(FT; water=species === :water) :
                   kind === :absco ? _pressure_test_absco(FT) : _pressure_test_regular_lut(FT)
        params = _pressure_test_parameters(FT, absorber; species)
        model, lin = model_from_parameters(LinMode(), params)
        h = FT === Float32 ? FT(0.5) : FT(0.05)
        plus, minus = deepcopy(params), deepcopy(params)
        plus.p[end] += h
        minus.p[end] -= h
        mplus, mminus = model_from_parameters(plus), model_from_parameters(minus)
        fd = (mplus.τ_abs[1] - mminus.τ_abs[1]) / (2h)
        relative_error = norm(lin.τ̇_abs_psurf[1] - fd) / norm(fd)
        @test relative_error < (FT === Float32 ? 2e-3 : 2e-6)
        @test eltype(lin.τ̇_abs_psurf[1]) === FT
        @test all(iszero, lin.τ̇_abs_psurf[1][:, 1])
        @test mplus.profile.p_full[end] - mminus.profile.p_full[end] ≈ h
        @test model.τ_abs[1] ≈ model_from_parameters(params).τ_abs[1]

        # Disabling H2O VMR columns must retain the surface-pressure response.
        if species === :water
            _, lin_fixed_water = model_from_parameters(LinMode(), params;
                compute_h2o_jacobians=false)
            @test lin_fixed_water.τ̇_abs_psurf[1] ≈ lin.τ̇_abs_psurf[1]
        end
    end
end

@testset "Native ABSCO pressure intervals and clamping" begin
    for FT in (Float32, Float64)
        lut = _pressure_test_absco(FT)
        grid = FT[6199, 6199.8, 6200, 6200.2, 6201]
        rtol = FT === Float32 ? 2e-6 : 2e-13
        for (p, interval_sum) in ((300, 700), (600, 1700), (800, 1700))
            # At the interior knot 600 hPa, the forward interpolator selects
            # the interval on the right. Match that directional derivative.
            expected = FT[0, 1, 2, 3, 0] .* FT(1e-23) .*
                (FT(0.001) + FT(0.000001) * FT(interval_sum))
            actual = CoreRT._absorption_pressure_derivative(lut, grid, FT(p), FT(280), FT(0.02))
            @test actual ≈ expected rtol=rtol
            @test eltype(actual) === FT
        end
        for p in (50, 100, 1100, 1200)
            @test all(iszero, CoreRT._absorption_pressure_derivative(
                lut, grid, FT(p), FT(280), FT(0.02)))
        end
    end
end

@testset "LBL pressure derivative and helper accumulation" begin
    for FT in (Float32, Float64), water in (false, true)
        absorber = _pressure_test_line_model(FT; water)
        params = _pressure_test_parameters(FT, absorber; species=water ? :water : :variable)
        model = model_from_parameters(params)
        profile = model.profile
        grid = params.spec_bands[1]
        p, T = profile.p_full[end], profile.T[end]
        broadener = water ? profile.vmr_h2o[end] / (1 + profile.vmr_h2o[end]) : nothing
        σ(p) = collect(CoreRT._layer_absorption_cross_section(absorber, grid, p, T, broadener))
        h = FT === Float32 ? FT(0.5) : FT(0.05)
        fd = (σ(p + h) - σ(p - h)) / (2h)
        derivative = CoreRT._absorption_pressure_derivative(absorber, grid, p, T, broadener)
        @test derivative ≈ fd rtol=(FT === Float32 ? 2e-3 : 2e-6)
        @test eltype(derivative) === FT
        vmr = water ? profile.vmr_h2o : profile.vmr["CH4"]
        expected = zeros(FT, size(model.τ_abs[1]))
        expected[:, end] .= fd .* profile.vcd_dry[end] .* vmr[end] ./ 2
        for linearized in (false, true)
            τ = zeros(FT, size(expected))
            pressure_tangent = zero(τ)
            gas_tangent = zeros(FT, length(profile.T), size(τ)...)
            if water
                args = linearized ? (τ, gas_tangent, 1, absorber, grid, profile) :
                                    (τ, absorber, grid, profile)
                CoreRT.compute_h2o_absorption_profile!(args...; pressure_tangent)
                # Accumulation is essential when several gases share a band.
                CoreRT.compute_h2o_absorption_profile!(args...; pressure_tangent)
            else
                args = linearized ? (τ, gas_tangent, 1, absorber, grid, vmr, profile) :
                                    (τ, absorber, grid, vmr, profile)
                CoreRT.compute_absorption_profile!(args...; pressure_tangent)
                CoreRT.compute_absorption_profile!(args...; pressure_tangent)
            end
            @test pressure_tangent ≈ 2expected rtol=(FT === Float32 ? 2e-3 : 2e-6)
            @test τ ≈ 2model.τ_abs[1]
        end
    end
end

@testset "Legacy interpolator pressure derivative" begin
    for FT in (Float32, Float64)
        ν = range(FT(6199), FT(6201); length=3)
        p = range(FT(100), FT(1100); length=3)
        T = range(FT(200), FT(320); length=3)
        table = FT[FT(1e-23) * (1 + FT(0.001) * pp + FT(0.002) * tt + FT(i))
            for i in eachindex(ν), pp in p, tt in T]
        itp = Interpolations.interpolate(table, Interpolations.BSpline(Interpolations.Linear()))
        absorber = vSmartMOM.Absorption.InterpolationModel(itp, 2, 1, ν, p, T)
        grid = FT[6198, 6199, 6200, 6201, 6202]
        derivative = CoreRT._absorption_pressure_derivative(absorber, grid, FT(750), FT(280), nothing)
        @test derivative ≈ FT[0, 1, 1, 1, 0] .* FT(1e-26) rtol=(FT === Float32 ? 2e-6 : 2e-13)
    end
end

@testset "Surface-pressure radiance and fixed convolution" begin
    for kind in (:line, :absco), external_solar in (false, true)
        absorber = kind === :line ? _pressure_test_line_model(Float64) : _pressure_test_absco(Float64)
        params = _pressure_test_parameters(Float64, absorber; species=:mixed)
        model, lin = model_from_parameters(LinMode(), params; external_solar)
        analytic = rt_run(model, lin, 0, size(lin.τ̇_abs[1], 1), 1)
        h = 0.05
        plus, minus = deepcopy(params), deepcopy(params)
        plus.p[end] += h
        minus.p[end] -= h
        mplus = model_from_parameters(plus; external_solar)
        mminus = model_from_parameters(minus; external_solar)
        forward(m) = external_solar ? rt_run_toa(m) : rt_run(m).toa
        fd = (forward(mplus) - forward(mminus)) / (2h)
        pressure_column = analytic.toa_jacobian[:, :, :, 1]
        @test pressure_column ≈ fd rtol=5e-5

        # Explicit, normalized detector weights are independent of atmosphere
        # construction. First validate native radiance, then apply the same
        # fixed linear operator to the Jacobian and the rebuilt spectra.
        weights = [1.0 2 1 0 0; 0 0 1 2 1] ./ 4
        @test weights * vec(pressure_column) ≈ weights * vec(fd) rtol=5e-5
    end
end

@testset "CUDA pressure cross sections" begin
    if CUDA.functional()
        for FT in (Float32, Float64), kind in (:line, :absco)
            fixture = kind === :line ? _pressure_test_line_model : _pressure_test_absco
            cpu = fixture(FT)
            gpu = fixture(FT; architecture=AtmosphericAbsorption.GPU())
            grid = FT[6199.8, 6199.95, 6200, 6200.05, 6200.2]
            p, T = FT(750), FT(280)
            broadener = kind === :absco ? FT(0.02) : nothing
            cpu_derivative = CoreRT._absorption_pressure_derivative(cpu, grid, p, T, broadener)
            gpu_derivative = CoreRT._absorption_pressure_derivative(gpu, grid, p, T, broadener)
            @test gpu_derivative ≈ cpu_derivative rtol=(FT === Float32 ? 1e-5 : 2e-12)
            h = FT === Float32 ? FT(0.5) : FT(0.05)
            σ(p) = Array(CoreRT._layer_absorption_cross_section(gpu, grid, p, T, broadener))
            fd = (σ(p + h) - σ(p - h)) / (2h)
            @test gpu_derivative ≈ fd rtol=(FT === Float32 ? 2e-3 : 2e-6)
        end
    else
        @test_skip CUDA.functional()
    end
end
