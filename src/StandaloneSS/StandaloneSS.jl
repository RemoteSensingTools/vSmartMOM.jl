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
       determine_required_l_aerosol,
       determine_required_l_from_moments,
       determine_required_nbrdf,
       determine_required_nbrdf_coxmunk,
       determine_required_nquad,
       determine_required_nquad_inner,
       determine_required_nstreams,
       exact_phase_function,
       run_exact_ss

include("types.jl")
include("kernels.jl")
include("quadrature_required.jl")
include("solver.jl")

end
