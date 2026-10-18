using vSmartMOM
using Test
import YAML

@testset "Float32 Consistency" begin

    yaml64 = "test_parameters/PureRayleighParameters.yaml"
    yaml32 = "test_parameters/PureRayleighParameters_f32.yaml"

    @testset "Float64 baseline" begin
        params64 = vSmartMOM.parameters_from_yaml(yaml64)
        params64.max_m = 2
        params64.l_trunc = 10
        # Keep `truncation` in sync with the mutated `l_trunc` — model_from_parameters
        # no longer auto-rebuilds the default δBGE so that explicit user truncation
        # (e.g. NoTruncation()) isn't silently overwritten.
        params64.truncation = vSmartMOM.Scattering.δBGE{Float64}(10, params64.Δ_angle)
        params64.architecture = vSmartMOM.CPU()
        model64 = model_from_parameters(params64)
        R64, _, _, _, _, _, _ = rt_run(model64)
        @test eltype(R64) == Float64
        @test all(isfinite, R64)
        @test size(R64, 1) == length(params64.vza)
    end

    @testset "Float32 forward" begin
        params32 = vSmartMOM.parameters_from_yaml(yaml32)
        params32.max_m = 2
        params32.l_trunc = 10
        params32.truncation = vSmartMOM.Scattering.δBGE{Float32}(10, params32.Δ_angle)
        params32.architecture = vSmartMOM.CPU()
        model32 = model_from_parameters(params32)
        R32, _, _, _, _, _, _ = rt_run(model32)
        @test eltype(R32) == Float32
        @test all(isfinite, R32)
        @test size(R32, 1) == length(params32.vza)
    end

    @testset "Float32 vs Float64 agreement" begin
        params64 = vSmartMOM.parameters_from_yaml(yaml64)
        params64.max_m = 2
        params64.l_trunc = 10
        # Keep `truncation` in sync with the mutated `l_trunc` — model_from_parameters
        # no longer auto-rebuilds the default δBGE so that explicit user truncation
        # (e.g. NoTruncation()) isn't silently overwritten.
        params64.truncation = vSmartMOM.Scattering.δBGE{Float64}(10, params64.Δ_angle)
        params64.architecture = vSmartMOM.CPU()
        model64 = model_from_parameters(params64)
        R64, _, _, _, _, _, _ = rt_run(model64)
        
        params32 = vSmartMOM.parameters_from_yaml(yaml32)
        params32.max_m = 2
        params32.l_trunc = 10
        params32.truncation = vSmartMOM.Scattering.δBGE{Float32}(10, params32.Δ_angle)
        params32.architecture = vSmartMOM.CPU()
        model32 = model_from_parameters(params32)
        R32, _, _, _, _, _, _ = rt_run(model32)
        
        R64_f32 = Float32.(R64)
        rel_diff = abs.(R64_f32 .- R32) ./ max.(abs.(R64_f32), Float32(1e-10))
        max_rel = maximum(rel_diff)
        mean_rel = sum(rel_diff) / length(rel_diff)
        
        println("  Float32 vs Float64: max_rel=$max_rel, mean_rel=$mean_rel")
        @test max_rel < 0.01  # 1% relative tolerance
    end
end

@testset "Float32 LinMode tangents use the model float type" begin
    # A Float64 CO2 entry in absorption.vmr must not leave Float64 tangent
    # arrays in a Float32 model (they are coerced to the model float type).
    d = YAML.load_file("test_parameters/JacobianTestFast.yaml")
    d["radiative_transfer"]["float_type"] = "Float32"
    d["absorption"]["vmr"]["CO2"] = 400e-6
    model, lin = model_from_parameters(LinMode(), read_parameters(d))
    @test eltype(model.optics.τ_abs[1]) == Float32
    for field in (:τ̇_abs, :τ̇_aer, :τ̇_rayl_psurf, :τ̇_abs_psurf, :τ̇_aer_psurf)
        @test eltype(first(getfield(lin, field))) == Float32
    end
end
