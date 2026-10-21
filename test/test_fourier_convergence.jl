using Test
using vSmartMOM
using vSmartMOM.CoreRT
import YAML

@testset "Fourier convergence strategy guards" begin
    @test RTNumericalParameters{Float64}().fourier_convergence isa
          AllFourierMoments
    @test_throws ArgumentError IntensityConvergence(0.0)
    @test_throws ArgumentError IntensityConvergence(-1e-5)
    @test_throws ArgumentError IntensityConvergence(1e-5; min_m=2)
    @test_throws ArgumentError IntensityConvergence(1e-5; n_consecutive=0)
    @test IntensityConvergence(1e-5; n_consecutive=1).n_consecutive == 1
    @test IntensityConvergence(1e-5).min_m == 3
    @test IntensityConvergence(1e-5).n_consecutive == 2
    @test_throws ArgumentError StokesConvergence(0.0)
    @test_throws ArgumentError StokesConvergence(1e-5; min_m=2)
    @test_throws ArgumentError StokesConvergence(1e-5; n_consecutive=0)
    @test StokesConvergence(1e-5; n_consecutive=1).n_consecutive == 1
    @test StokesConvergence(1e-5).min_m == 3
    @test StokesConvergence(1e-5).n_consecutive == 2
    convergence = StokesConvergence(1e-5)
    @test CoreRT._fourier_guard_through_m(convergence, 0) == 0
    @test CoreRT._fourier_guard_through_m(convergence, 1) == 1
    @test CoreRT._fourier_guard_through_m(convergence, 2) == 2
    @test CoreRT._fourier_guard_through_m(convergence, 15) == 2
end

@testset "one passing moment exits at m=3 after protected m=2" begin
    convergence = StokesConvergence(1e-5; n_consecutive=1)
    output = zeros(Float64, 1, 3, 1)
    outputs = CoreRT._fourier_outputs(output)
    snapshots = CoreRT._fourier_snapshots(outputs)

    for m in 0:2
        m == 0 && (output[1, 1, 1] = 1)
        m == 2 && (output[1, 2, 1] = 0.125)
        stop, passes = CoreRT._fourier_convergence_step!(
            convergence, outputs, snapshots, 0, m, 15)
        @test !stop && passes == 0
    end
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, 0, 3, 15)
    @test stop && passes == 1
end

@testset "m>=3 guard retains scalar Rayleigh beta_2" begin
    convergence = IntensityConvergence(1e-5; n_consecutive=2)
    output = zeros(Float64, 1, 1, 1)
    outputs = CoreRT._fourier_outputs(output)
    snapshots = CoreRT._fourier_snapshots(outputs)

    output[1, 1, 1] = 1
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, 0, 0, 2)
    @test !stop && passes == 0
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 1, 2)
    @test !stop && passes == 0
    output[1, 1, 1] += 0.125
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 2, 2)
    @test !stop && passes == 0
    @test output[1, 1, 1] == 1.125
end

@testset "full Stokes convergence waits beyond m=2" begin
    convergence = StokesConvergence(1e-5; n_consecutive=2)
    output = zeros(Float64, 1, 3, 1)
    outputs = CoreRT._fourier_outputs(output)
    snapshots = CoreRT._fourier_snapshots(outputs)

    output[1, 1, 1] = 1
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, 0, 0, 15)
    @test !stop && passes == 0
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 1, 15)
    @test !stop && passes == 0
    output[1, 2, 1] = 0.125
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 2, 15)
    @test !stop && passes == 0
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 3, 15)
    @test !stop && passes == 1
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 4, 15)
    @test stop && passes == 2
end

@testset "Fourier convergence YAML parsing" begin
    input = YAML.load_file("test_parameters/JacobianTestFast.yaml")
    input["radiative_transfer"]["numerics"] = Dict(
        "fourier_convergence" => "intensity",
        "fourier_tolerance" => 2e-6,
        "fourier_min_m" => 4,
        "fourier_n_consecutive" => 2,
    )
    mktempdir() do directory
        path = joinpath(directory, "fourier.yaml")
        YAML.write_file(path, input)
        strategy = parameters_from_yaml(path).numerics.fourier_convergence
        @test strategy isa IntensityConvergence
        @test strategy.tolerance == 2e-6
        @test strategy.min_m == 4
        @test strategy.n_consecutive == 2

        input["radiative_transfer"]["numerics"] = Dict(
            "fourier_convergence" => "stokes",
            "fourier_tolerance" => 3e-6,
            "fourier_n_consecutive" => 2,
        )
        YAML.write_file(path, input)
        strategy = parameters_from_yaml(path).numerics.fourier_convergence
        @test strategy isa StokesConvergence
        @test strategy.tolerance == 3e-6
        @test strategy.min_m == 3
        @test strategy.n_consecutive == 2

        input["radiative_transfer"]["numerics"] = Dict(
            "fourier_convergence" => "all")
        YAML.write_file(path, input)
        @test parameters_from_yaml(path).numerics.fourier_convergence isa
              AllFourierMoments

        input["radiative_transfer"]["numerics"] = Dict(
            "fourier_convergence" => "bogus")
        YAML.write_file(path, input)
        @test_throws ArgumentError parameters_from_yaml(path)
    end
end

@testset "full-series and runtime-converged forward solves" begin
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    FT = params.float_type
    function set_convergence!(parameters, convergence)
        numerics = parameters.numerics
        parameters.numerics = RTNumericalParameters{FT}(
            dτ_max_threshold=numerics.dτ_max_threshold,
            dτ_min_floor=numerics.dτ_min_floor,
            blas_threads=numerics.blas_threads,
            verbose=numerics.verbose,
            fourier_convergence=convergence,
            ss_correction=numerics.ss_correction)
        return parameters
    end

    full_model = model_from_parameters(set_convergence!(params, AllFourierMoments()))
    R_full, T_full = rt_run(full_model)
    @test CoreRT._LAST_FOURIER_M_USED[] == full_model.solver.m_max_bands[1]

    tiny_model = model_from_parameters(set_convergence!(
        params, IntensityConvergence(1e-30)))
    R_tiny, T_tiny = rt_run(tiny_model)
    @test Array(R_tiny) == Array(R_full)
    @test Array(T_tiny) == Array(T_full)
    @test CoreRT._LAST_FOURIER_M_USED[] == tiny_model.solver.m_max_bands[1]

    tolerance = 1e-4
    converged_model = model_from_parameters(set_convergence!(
        params, IntensityConvergence(tolerance)))
    R_converged, T_converged = rt_run(converged_model)
    m_used = CoreRT._LAST_FOURIER_M_USED[]
    @test m_used <= converged_model.solver.m_max_bands[1]
    I_full = Array(R_full)[:, 1, :]
    I_converged = Array(R_converged)[:, 1, :]
    relative_error = maximum(abs.(I_converged .- I_full) ./
                             max.(abs.(I_full), 1e-12))
    @test relative_error <= 20 * tolerance
end
