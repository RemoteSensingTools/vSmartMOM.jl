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

@testset "one passing moment exits at m=3 after the protected m=2 term" begin
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

    # Scalar Rayleigh schematic: beta_0 is nonzero, beta_1 is structurally
    # zero, and beta_2 is nonzero. No convergence pass is allowed through m=2.
    output[1, 1, 1] = 1
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, 0, 0, 2)
    @test !stop
    @test passes == 0

    # m=1 would pass a naive test because its contribution is exactly zero.
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 1, 2)
    @test !stop
    @test passes == 0

    output[1, 1, 1] += 0.125
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 2, 2)
    @test !stop
    @test passes == 0
    @test output[1, 1, 1] == 1.125

    # m_max=2: the exact Rayleigh-supported series is now complete. No
    # convergence test was needed and no m=3/m=4 work was manufactured.
end

@testset "full Stokes convergence waits beyond the m=2 Q term" begin
    convergence = StokesConvergence(1e-5; n_consecutive=2)
    output = zeros(Float64, 1, 3, 1)
    outputs = CoreRT._fourier_outputs(output)
    snapshots = CoreRT._fourier_snapshots(outputs)

    output[1, 1, 1] = 1
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, 0, 0, 15)
    @test !stop && passes == 0

    # Neither structural zeros nor a nonzero m=2 term are tested.
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 1, 15)
    @test !stop && passes == 0

    output[1, 2, 1] = 0.125
    stop, passes = CoreRT._fourier_convergence_step!(
        convergence, outputs, snapshots, passes, 2, 15)
    @test !stop && passes == 0

    # Both m=3 and m=4 must be negligible in I, Q, and U.
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
        parameters = parameters_from_yaml(path)
        strategy = parameters.numerics.fourier_convergence
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
    end
end
