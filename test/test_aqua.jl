using Test
using Aqua
using vSmartMOM

@testset "Aqua: package health" begin
    Aqua.test_all(vSmartMOM;
        # Method-ambiguity check excluded: the package has pre-existing broad
        # signatures (e.g. AbstractRTModel dispatches) that produce false positives
        # and are tracked in a separate clean-up pass. Keeping the check in the
        # test_quality.jl Quality Gates testset with its own tracking.
        ambiguities = false,
        # The CUDA weak dependency (ext/vSmartMOMCUDAExt.jl) and lazy artifact
        # precompilation can leave background tasks in the test environment on
        # some Julia versions. Disabled to prevent non-deterministic CI failures;
        # tracked as a separate clean-up item.
        persistent_tasks = false,
    )
end
