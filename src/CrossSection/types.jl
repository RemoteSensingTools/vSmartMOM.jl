"""
    type AbstractCrossSection
Abstract Cross Section type for generic cross section calculations
"""
abstract type AbstractCrossSection end

"""
    struct HitranCrossSection{FT}

An [`AbstractCrossSection`](@ref) type struct, which provides all HITRAN line
parameters needed to compute absorption cross sections

See https://hitran.org/docs/definitions-and-units/ for details

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct HitranTable{FT<:AbstractFloat} <: AbstractCrossSection
    "The molecular species identification (ID) number"
    mol::Array{Int,1}
    "The isotopologue ID number"
    iso::Array{Int,1}
    "The wavenumber of the spectral line transition (cm-1) in vacuum"
    νᵢ::Array{FT,1}
    "The spectral line intensity (cm−1/(molecule·cm−2)) at Tref=296K"
    Sᵢ::Array{FT,1}
    "The Einstein-A coefficient (s-1) of a transition"
    Aᵢ::Array{FT,1}
    "The air-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm"
    γ_air::Array{FT,1}
    "The self-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm"
    γ_self::Array{FT,1}
    "The lower-state energy of the transition (cm-1)"
    E″::Array{FT,1}
    "The coefficient of the temperature dependence of the air-broadened half width"
    n_air::Array{FT,1}
    "The pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to the vacuum transition wavenumber νij"
    δ_air::Array{FT,1}
    "The electronic and vibrational quantum numbers and labels for the upper state of a transition"
    global_upper_quanta::Array{String,1}
    "The electronic and vibrational quantum numbers and labels for the lower state of a transition"
    global_lower_quanta::Array{String,1}
    "Rotational, hyperfine and other quantum numbers and labels for the upper state of a transition"
    local_upper_quanta::Array{String,1}
    "Rotational, hyperfine and other quantum numbers and labels for the lower state of a transition"
    local_lower_quanta::Array{String,1}
    "Ordered list of indices corresponding to uncertainty estimates of transition parameters"
    ierr::Array{String,1}
    "Ordered list of reference identifiers for transition parameters"
    iref::Array{String,1}
    "A flag indicating the presence of additional data and code relating to line-mixing"
    line_mixing_flag::Array{String,1}
    "The upper state degeneracy"
    g′::Array{FT,1}
    "The lower state degeneracy"
    g″::Array{FT,1}
end

@enum BroadeningFunction doppler=1 lorentz=2 voigt=3




abstract type AbstractComplexErrorFunction end

"Humlicek only formulation for Complex Error Function"
struct HumlicekErrorFunction <: AbstractComplexErrorFunction end
"Mix of Humlicek and Weidemann (N=32) Error Function, suggested for Voigt function"
struct HumlicekWeidemann32VoigtErrorFunction <: AbstractComplexErrorFunction end
"Mix of Humlicek and Weidemann (N=32) Error Function, suggested for Speed Dependent Voigt function"
struct HumlicekWeidemann32SDErrorFunction <: AbstractComplexErrorFunction end

struct CPF12ErrorFunction <: AbstractComplexErrorFunction end

"Mix of Humplicek and erfc Special Function for Voigt"
struct ErfcHumliErrorFunctionVoigt  <: AbstractComplexErrorFunction end

"Mix of Humplicek and erfc Special Function for Voigt"
struct ErfcHumliErrorFunctionSD  <: AbstractComplexErrorFunction end

"erfc Special Function for Voigt"
struct ErfcErrorFunction  <: AbstractComplexErrorFunction end


struct HitranEmptyError <: Exception end
Base.showerror(io::IO, e::HitranEmptyError) = print(io, e, "No HITRAN records match the parameters")
