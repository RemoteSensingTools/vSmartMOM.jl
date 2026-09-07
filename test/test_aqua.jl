using Test
using Aqua
using vSmartMOM

@testset "Aqua: package health" begin
    Aqua.test_all(vSmartMOM;
        # Recursive dispatch ambiguity coverage includes loaded CUDA methods.
        ambiguities = true,
        # CUDA initialization (ext/vSmartMOMCUDAExt.jl) and lazy artifact
        # precompilation can leave background tasks in the test environment on
        # some Julia versions. Disabled to prevent non-deterministic CI failures;
        # tracked as a separate clean-up item.
        persistent_tasks = false,
    )
end
