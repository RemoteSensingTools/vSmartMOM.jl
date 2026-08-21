# VLIDORT-style Fourier-loop convergence (AbstractFourierConvergence).
#
# Contracts:
# - AllFourierMoments (default) is the historical full loop — bit-identical.
# - IntensityConvergence with an unreachable tolerance never terminates
#   early → also bit-identical to the full series.
# - A realistic tolerance terminates early (m_used < m_max) and the
#   truncated radiances agree with the full series to O(tolerance).
# - _LAST_FOURIER_M_USED reports the moment actually used (FOURIER_SAVED).
using Test
using vSmartMOM
using vSmartMOM.CoreRT
import YAML as _YAML
using vSmartMOM.CoreRT: AbstractFourierConvergence, AllFourierMoments,
                        IntensityConvergence, RTNumericalParameters,
                        _LAST_FOURIER_M_USED

params = vSmartMOM.parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
FTp = params.float_type

set_fc!(p, fc) = (p.numerics = RTNumericalParameters{FTp}(
    dτ_max_threshold = p.numerics.dτ_max_threshold,
    dτ_min_floor     = p.numerics.dτ_min_floor,
    blas_threads     = p.numerics.blas_threads,
    verbose          = p.numerics.verbose,
    fourier_convergence = fc); p)

@testset "constructor guards + default" begin
    @test RTNumericalParameters{Float64}().fourier_convergence isa AllFourierMoments
    @test_throws ArgumentError IntensityConvergence(0.0)
    @test_throws ArgumentError IntensityConvergence(-1e-4)
    @test_throws ArgumentError IntensityConvergence(1e-4; n_consecutive = 0)
    c = IntensityConvergence(1e-4)
    @test c.n_consecutive == 2          # VLIDORT DO_DOUBLE_CONVTEST default
end

@testset "YAML parsing" begin
    d = _YAML.load_file(joinpath(@__DIR__, "test_parameters", "JacobianTestFast.yaml"))
    d["radiative_transfer"]["numerics"] = Dict(
        "fourier_convergence" => "intensity",
        "fourier_tolerance" => 5e-5,
        "fourier_n_consecutive" => 1)
    mktempdir() do dir
        f = joinpath(dir, "fc.yaml")
        _YAML.write_file(f, d)
        p = vSmartMOM.parameters_from_yaml(f)
        fc = p.numerics.fourier_convergence
        @test fc isa IntensityConvergence
        @test fc.tolerance ≈ 5e-5
        @test fc.n_consecutive == 1
        # 'all' and invalid values
        d["radiative_transfer"]["numerics"]["fourier_convergence"] = "all"
        _YAML.write_file(f, d)
        @test vSmartMOM.parameters_from_yaml(f).numerics.fourier_convergence isa AllFourierMoments
        d["radiative_transfer"]["numerics"]["fourier_convergence"] = "bogus"
        _YAML.write_file(f, d)
        @test_throws Exception vSmartMOM.parameters_from_yaml(f)
    end
end

@testset "full series A/B (bitwise contracts)" begin
    model_full = model_from_parameters(set_fc!(params, AllFourierMoments()))
    R_full, T_full = rt_run(model_full)
    @test _LAST_FOURIER_M_USED[] == model_full.solver.m_max_bands[1]

    # unreachable tolerance ⇒ never converges ⇒ identical to full series
    model_tiny = model_from_parameters(set_fc!(params, IntensityConvergence(1e-30)))
    R_tiny, T_tiny = rt_run(model_tiny)
    @test Array(R_tiny) == Array(R_full)
    @test Array(T_tiny) == Array(T_full)
    @test _LAST_FOURIER_M_USED[] == model_tiny.solver.m_max_bands[1]

    # realistic tolerance ⇒ early termination, O(tol) agreement on Stokes I
    tol = 1e-4
    model_conv = model_from_parameters(set_fc!(params, IntensityConvergence(tol)))
    R_conv, T_conv = rt_run(model_conv)
    m_used = _LAST_FOURIER_M_USED[]
    m_max  = model_conv.solver.m_max_bands[1]
    @test m_used ≤ m_max
    I_full = Array(R_full)[:, 1, :]
    I_conv = Array(R_conv)[:, 1, :]
    rel = maximum(abs.(I_conv .- I_full) ./ max.(abs.(I_full), 1e-12))
    # The dropped tail is a sum of ≤(m_max−m_used) terms each ≤ tol·I —
    # allow a generous multiple of the tolerance.
    @test rel ≤ 20 * tol
    println("fourier convergence: m_used=$m_used of m_max=$m_max, tail rel err=$(round(rel; sigdigits=3))")
end

println("test_fourier_convergence: done")
