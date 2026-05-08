"""
    InelasticScattering

Molecular Raman and Cabannes/Rayleigh helper types and precomputations used by
CoreRT inelastic-scattering modes.
"""
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
    # stellar inelastic scattering
    include("stellar_inelastic_helper.jl")
    include("raman_stellar_prop.jl")

    "Speed of light in vacuum in `[cm/s]`"
    global const c = 2.99792458e10
    "Planck constant in `[erg ̇s]`"
    global const h = 6.62607015e-27
    "Boltzmann constant in `[erg/K]`"
    global const k_B = 1.380649e-16
    "Wavenumber-to-wavelength conversion: `λ_nm = nm_per_m / ν̄_cm⁻¹`.
    Despite the name (preserved for cross-module consistency with
    `Absorption.constants`), this is `1e7` (= nm per cm), the Raman
    band-grid code uses ν̄ in cm⁻¹ throughout."
    global const nm_per_m = 1.0e7

    export compute_γ_air_Cabannes! 
    export compute_γ_air_Rayleigh!
    export compute_effective_coefficents!
    export compute_σ_Rayl_coeff!, compute_stellar_Rayl
    export compute_σ_Rayl_VibRaman_coeff_hires!
    export compute_energy_levels!
    export apply_lineshape!, get_n₀_n₁
    export get_greek_raman, get_greek_raman_VS

    #export compute_σ_Raman_coeff!
    export compute_σ_VibRaman_coeff!, compute_σ_RoVibRaman_coeff!
    export getRamanSSProp!
    export AbstractRamanType
    export noRS, RRS, VS_0to1, VS_1to0
    export noRS_plus, VS_0to1_plus, VS_1to0_plus
    export has_inelastic, uses_cabannes_phase, needs_interaction_workspace,
           needs_rayleigh_expansion, normalize_raman_weights!


    # stellar inelastic scattering
    export sol_RRS, sol_VS_0to1, sol_VS_1to0
    export sol_VS_0to1_plus, sol_VS_1to0_plus

    @doc """
        apply_lineshape!(Δνᵢ, σᵢ, λ₀, Δν_out, σ_out, pressure, temperature, molMass; wavelength_flag=false)

    Broaden discrete Raman transitions onto `Δν_out`, writing cross sections
    into `σ_out` for incident wavelength `λ₀`.
    """ apply_lineshape!

    @doc """
        compute_energy_levels!(mol; vmax=2, Jmax=30)

    Populate molecular vibrational/rotational energy levels in cm⁻¹ for
    `v = 0:vmax` and `J = 0:Jmax`.
    """ compute_energy_levels!

    @doc """
        compute_stellar_Rayl(λ₀, h2)

    Return the H₂ Rayleigh cross section at stellar/solar wavelength `λ₀`.
    """ compute_stellar_Rayl

    @doc """
        compute_γ_air_Cabannes!(λ₀, rs_or_n2, [o2])

    Compute the effective Cabannes Greek coefficient `γ` and elastic Cabannes
    fraction for an air mixture at wavelength `λ₀` in nm.
    """ compute_γ_air_Cabannes!

    @doc """
        compute_γ_air_Rayleigh!(λ₀, rs_or_n2, [o2])

    Compute the effective Rayleigh Greek coefficient `γ` and Rayleigh cross
    section for an air mixture at wavelength `λ₀` in nm.
    """ compute_γ_air_Rayleigh!

    @doc """
        compute_σ_Rayl_VibRaman_coeff_hires!(T, mol; Jmax=30)

    Update high-resolution Rayleigh and vibrational Raman transition
    prefactors over rotational levels up to `Jmax`.
    """ compute_σ_Rayl_VibRaman_coeff_hires!

    @doc """
        compute_σ_Rayl_coeff!(mol)

    Update the molecular Rayleigh cross-section prefactor stored on `mol`.
    """ compute_σ_Rayl_coeff!

    @doc """
        compute_σ_VibRaman_coeff!(T, mol; vmax=2, Jmax=30)

    Update Stokes and anti-Stokes vibrational Raman cross-section prefactors
    for `mol` at temperature `T`.
    """ compute_σ_VibRaman_coeff!

    @doc """
        get_greek_raman(rs, n2, o2)

    Return Raman phase-function Greek coefficients for rotational or
    rovibrational Raman scattering in an N₂/O₂ atmosphere.
    """ get_greek_raman



end # module
