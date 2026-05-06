# VLIDORT baseline validation suite — standalone runner.
#
# Compares vSmartMOM against VLIDORT 2.8.3 shipped saved_results tables
# (no PyVLIDORT/Fortran runtime needed; reference data is committed to
# reference_data/).
#
# Usage:
#   cd test && julia --project=. vlidort_baseline/runtests.jl
#
# This file is intentionally NOT included from test/runtests.jl. It can be
# wired in later when we want CI integration.

using Test

# Switch to this directory so relative paths in cases/ resolve.
const _ORIGINAL_PWD = pwd()
cd(@__DIR__)
try
    include(joinpath(@__DIR__, "harness.jl"))

    @testset "VLIDORT baseline" begin
        include(joinpath(@__DIR__, "cases", "case_A_siewert2000.jl"))
        # Future cases (B1/B2 scalar, C1/C2 vector) will be added here.
    end
finally
    cd(_ORIGINAL_PWD)
end
