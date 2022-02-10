module InelasticScattering

    using DocStringExtensions
    using Parameters
    using OffsetArrays
    using Interpolations
    using ..Scattering
    using UnPack
    #import PhysicalConstants.CODATA2018:c_0, h, k_B

    include("src/raman_constants.jl")
    include("src/molecular_constructors.jl")
    include("src/inelastic_cross_section.jl")
    include("src/apply_lineshape.jl")
    include("types.jl")
    include("inelastic_helper.jl")
    include("raman_atmo_prop.jl")

    "Speed of light in vacuum in `[cm/s]`"
    global const c = 2.99792458e10
    "Planck constant in `[erg ̇s]`"
    global const h = 6.62607015e-27
    "Boltzmann constant in `[erg/K]`"
    global const k_B = 1.380649e-16

    export compute_effective_coefficents!
    export compute_σ_Rayl_coeff!
    export compute_σ_Rayl_VibRaman_coeff_hires!
    export compute_energy_levels!
    export apply_lineshape!

    #export compute_σ_Raman_coeff!
    export compute_σ_VibRaman_coeff!, compute_σ_RoVibRaman_coeff!
    export getRamanSSProp!
    export AbstractRamanType
    export noRS, RRS, VS_0to1, VS_1to0
    export noRS_plus, RRS_plus, VS_0to1_plus, VS_1to0_plus

end # module
