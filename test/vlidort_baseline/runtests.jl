# VLIDORT baseline validation suite — standalone runner.
#
# Compares vSmartMOM against VLIDORT 2.8.3 shipped saved_results tables
# (no PyVLIDORT/Fortran runtime needed; reference data is committed to
# reference_data/).
#
# Usage:
#   cd test && julia --project=. vlidort_baseline/runtests.jl
#
# It is also included from test/runtests.jl as part of the default suite.

using Test

# Switch to this directory so relative paths in cases/ resolve.
const _ORIGINAL_PWD = pwd()
cd(@__DIR__)
try
    include(joinpath(@__DIR__, "harness.jl"))

    @testset "VLIDORT baseline" begin
        include(joinpath(@__DIR__, "cases", "case_0_definition_checks.jl"))
        include(joinpath(@__DIR__, "cases", "case_A_siewert2000.jl"))
        include(joinpath(@__DIR__, "cases", "case_B_solar_tester.jl"))
        include(joinpath(@__DIR__, "cases", "case_C_solar_tester_vector.jl"))
        include(joinpath(@__DIR__, "cases", "case_D_rayleigh_brdf_coupled.jl"))
    end
finally
    cd(_ORIGINAL_PWD)
end
