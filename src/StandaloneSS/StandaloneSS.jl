"""
    StandaloneSS

Standalone exact single-scattering tools. Phase 1 implements scalar,
plane-parallel, Lambertian paths 1 and 2 with KernelAbstractions CPU kernels.
"""
module StandaloneSS

using KernelAbstractions

export AbsorptionSSContributor,
       ExactSSConfig,
       HGAerosolSSContributor,
       LambertianSSSurface,
       RayleighSSContributor,
       SSGeometry,
       TruncatedAndExactScatteringOpticalProperties,
       exact_phase_function,
       run_exact_ss

include("types.jl")
include("kernels.jl")
include("solver.jl")

end
