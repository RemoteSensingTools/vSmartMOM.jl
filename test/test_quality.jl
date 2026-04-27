using Aqua
using Test
using vSmartMOM

@testset "Aqua" begin
    Aqua.test_all(vSmartMOM;
        # Keep method-ambiguity rollout separate; current branch still has
        # legacy broad signatures that make this too noisy for the first gate.
        ambiguities = false,
        # CUDA and lazy artifact precompilation can leave background tasks in
        # the test environment. Keep this check separate from the package
        # hygiene gate until those initialization paths are isolated.
        persistent_tasks = false,
    )
end
