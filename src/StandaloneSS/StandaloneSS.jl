"""
    StandaloneSS

Standalone exact single-scattering tools. Phase 1 implements scalar,
plane-parallel exact paths with KernelAbstractions CPU kernels.
"""
module StandaloneSS

using FastGaussQuadrature: gausslegendre
using KernelAbstractions
using SpecialFunctions: erfc
using ..Architectures: AbstractArchitecture, CPU, GPU, MetalGPU, array_type, devi

export AbsorptionSSContributor,
       CoxMunkSSSurface,
       ExactSSConfig,
       HGAerosolSSContributor,
       LambertianSSSurface,
       RayleighSSContributor,
       SSMeasurementSelector,
       SSGeometry,
       TruncatedAndExactScatteringOpticalProperties,
       chain_rule_combine_dP,
       chain_rule_combine_dτ,
       chain_rule_combine_dϖ,
       chain_rule_combine_surface_brdf,
       determine_required_l_aerosol,
       determine_required_l_from_moments,
       determine_required_nbrdf,
       determine_required_nbrdf_coxmunk,
       determine_required_nquad,
       determine_required_nquad_inner,
       determine_required_nstreams,
       exact_phase_function,
       run_exact_ss,
       run_exact_ss_with_jacobians,
       selected_measurement_jacobian,
       selected_measurements

include("types.jl")
include("surfaces.jl")
include("kernels.jl")
include("quadrature_required.jl")
include("solver.jl")
include("measurement_selector.jl")
include("linearized_f2.jl")
include("chain_rule.jl")

end
