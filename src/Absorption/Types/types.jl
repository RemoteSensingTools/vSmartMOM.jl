#####
##### Types for holding input transition states, intensities and related
##### variables. (Currently intended for storing HITRAN database)
#####

using CUDA
using ..Architectures

"""
    struct HitranTable{FT}

A struct, which provides all HITRAN line parameters needed to compute 
absorption cross sections

See https://hitran.org/docs/definitions-and-units/ for details

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct HitranTable{FT<:AbstractFloat} 
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

#####
##### Types of Broadening Functions
##### Currently: Doppler, Lorentz, and Voigt
#####

"""
    type AbstractBroadeningFunction
Abstract Broadening Function type for generic line broadening function
"""
abstract type AbstractBroadeningFunction end

"Doppler line broadening"
struct Doppler <: AbstractBroadeningFunction end

"Lorentz line broadening"
struct Lorentz <: AbstractBroadeningFunction end

"Voigt line broadening"
struct Voigt <: AbstractBroadeningFunction end


#####
##### Types of Complex Error Functions
#####

"""
    type AbstractComplexErrorFunction
Abstract Complex Error Function type for generic complex error functions
"""
abstract type AbstractComplexErrorFunction end

"Humlicek only formulation for Complex Error Function"
struct HumlicekErrorFunction <: AbstractComplexErrorFunction end

"Mix of Humlicek and Weidemann (N=32) Error Function, suggested for Voigt function"
struct HumlicekWeidemann32VoigtErrorFunction <: AbstractComplexErrorFunction end

"Mix of Humlicek and Weidemann (N=32) Error Function, suggested for Speed Dependent Voigt function"
struct HumlicekWeidemann32SDErrorFunction <: AbstractComplexErrorFunction end

"Humlicek with a single rational approximation."
struct CPF12ErrorFunction <: AbstractComplexErrorFunction end

"Mix of Humplicek and erfc Special Function for Voigt"
struct ErfcHumliErrorFunctionVoigt  <: AbstractComplexErrorFunction end

"Mix of Humplicek and erfc Special Function for Voigt"
struct ErfcHumliErrorFunctionSD  <: AbstractComplexErrorFunction end

"erfc Special Function for Voigt"
struct ErfcErrorFunction  <: AbstractComplexErrorFunction end

#####
##### Types of models that can be used to calculate an absorption cross-
##### section. 
#####

"""
    type AbstractCrossSectionModel
Abstract Cross Section Model type for generic cross section calculations
Currently either a HitranModel or an InterpolationModel
"""
abstract type AbstractCrossSectionModel end

"""
    struct HitranModel{FT}

A struct which provides all model parameters needed for cross-section 
calculations using HITRAN data

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct HitranModel <: AbstractCrossSectionModel

    "Struct with hitran data"
    hitran::HitranTable
    "Broadening function (Doppler/Lorentz/Voigt)"
    broadening::AbstractBroadeningFunction
    "Wing cutoff [cm-1]"
    wing_cutoff::Real
    "VMR of gas itself [0-1]"
    vmr::Union{Real, Vector}
    "Complex Error Function to Use"
    CEF::AbstractComplexErrorFunction
    "Computer `Architecture` on which `Model` is run"
    architecture::AbstractArchitecture 
end

# @with_kw struct HitranModel{H, A<:AbstractArchitecture, AB, FT, CF } <: AbstractCrossSectionModel
#         hitran :: H     # Struct with hitran data
#     broadening :: AB    # Broadening function (Doppler/Lorentz/Voigt)
#    wing_cutoff :: FT    # Wing cutoff [cm-1]
#            vmr :: FT    # VMR of gas itself [0-1]
#            CEF :: CF    # "Complex Error Function to Use"
#            mol :: Int   # Molecule number
#            iso :: Int   # Isotope number 
#          ν_min :: FT    # min wavenumber (cm-1)
#          ν_max :: FT    # max wavenumber (cm-1)
#   architecture :: A     # Computer `Architecture` on which `Model` is run
# end

"""
    struct InterpolationModel{FT}

A struct which provides all model parameters needed for cross-section 
calculations using an Interpolator

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw struct InterpolationModel <: AbstractCrossSectionModel 

    "The interpolator"
    itp

    # Everything below is metadata associated with the interpolator

    "The molecular species identification (ID) number"
    mol::Int
    "The isotopologue ID number"
    iso::Int
    "Broadening function (Doppler/Lorentz/Voigt), as a string"
    broadening::AbstractBroadeningFunction
    "Wavelength grids"
    ν_grid::AbstractRange{<:Real}
    "Wavelength grids"
    p_grid::AbstractRange{<:Real}
    "Wavelength grids"
    t_grid::AbstractRange{<:Real}
    "Wing cutoff [cm-1]"
    wing_cutoff::Real
    "VMR of gas itself [0-1]"
    vmr::Real
    "Complex Error Function to Use"
    CEF::AbstractComplexErrorFunction
    "Computer `Architecture` on which `Model` is run"
    architecture::AbstractArchitecture     

end


#####
##### Types of Errors that may be thrown
#####

struct HitranEmptyError <: Exception end
Base.showerror(io::IO, e::HitranEmptyError) = print(io, e, "No HITRAN records match the parameters")
