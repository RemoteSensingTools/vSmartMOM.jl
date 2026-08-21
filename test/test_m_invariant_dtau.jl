# Bitwise A/B for the m-invariant max(τ·ϖ) cache in rt_run's Fourier loop.
#
# τ·ϖ does not depend on the Fourier moment, so caching one per-layer device
# reduction from m = 0 and reusing the scalar for every m feeds the EXACT
# same value into the scatter test and dτ/ndoubl derivation — the contract
# is bitwise equality of the radiances, not a tolerance.
using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.CoreRT: _M_INVARIANT_DTAU_CACHE_ENABLED

@testset "m-invariant dτ cache: bitwise rt_run A/B" begin
    @test _M_INVARIANT_DTAU_CACHE_ENABLED[] == true      # default on

    params = vSmartMOM.parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    model  = model_from_parameters(params)

    _M_INVARIANT_DTAU_CACHE_ENABLED[] = false
    R_off, T_off = rt_run(model)
    _M_INVARIANT_DTAU_CACHE_ENABLED[] = true
    R_on, T_on = rt_run(model)

    @test Array(R_on) == Array(R_off)
    @test Array(T_on) == Array(T_off)
    @test all(isfinite, Array(R_on))
end

println("test_m_invariant_dtau: done")
