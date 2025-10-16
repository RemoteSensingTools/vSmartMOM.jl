#=
 
This file contains all types that are used in the vSmartMOM module:

- `AtmosphericProfile` stores all relevant atmospheric profile information 
- `AbstractObsGeometry` specifies the RT geometry
- `RT_Aerosol` holds an Aerosol with additional RT parameters
- `AbstractQuadratureType` specifies the quadrature type to use
- `AbstractSourceType` specifies the source type
- `CompositeLayer` and `AddedLayer` specify the layer properties
- `AbstractScatteringInterface` specifies the scattering interface type
- `AbstractSurfaceType` specify the type of surface in the RT simulation
- `AbsorptionParameters`, `ScatteringParameters`, and `vSmartMOM_Model` hold model parameters
- `QuadPoints` holds quadrature points, weights, etc. 
- `ComputedAtmosphereProperties` and `ComputedLayerProperties` hold intermediate computed properties

=#
#=
"Struct for an atmospheric profile"
struct AtmosphericProfile{FT, VMR <: Union{Real, Vector}}
    "Temperature Profile"
    T::Array{FT,1}
    "Pressure Profile (Full)"
    p_full::Array{FT,1}
    "Specific humidity profile"
    q::Array{FT,1}
    "Pressure Levels"
    p_half::Array{FT,1}
    "H2O Volume Mixing Ratio Profile"
    vmr_h2o::Array{FT,1}
    "Vertical Column Density (Dry)"
    vcd_dry::Array{FT,1}
    "Vertical Column Density (H2O)"
    vcd_h2o::Array{FT,1}
    "Volume Mixing Ratio of Constituent Gases"
    vmr::Dict{String, VMR}
end=#

#"Types for describing atmospheric parameters"
#abstract type AbstractObsGeometry end
#=
mutable struct RT_Aerosol_Lin{}#FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "Aerosol"
    aerosol_lin::AerosolLin#{FT}
    #"Reference τ"
    #τ_ref#::FT
    #"Mode z (km)"
    #z₀#::FT
    #"Peak width"
    #σ₀#::FT
    #"Pressure peak (Pa)"
    #p₀#::FT
    #"Pressure peak width (Pa)"
    #σp#::FT
end
=#
"Abstract Type for Layer Ṙ,Ṫ and J̇ matrices"
abstract type AbstractLayerLin end

"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct CompositeLayerLin{FT} <: AbstractLayerLin 
    "Composite layer Reflectance matrix R (from + -> -)"
    Ṙ⁻⁺::AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    Ṙ⁺⁻::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from + -> +)"
    Ṫ⁺⁺::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from - -> -)"
    Ṫ⁻⁻::AbstractArray{FT,4}
    "Composite layer source matrix J (in + direction)"
    J̇₀⁺::AbstractArray{FT,4}
    "Composite layer source matrix J (in - direction)"
    J̇₀⁻::AbstractArray{FT,4}
end

"Added (Single) Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct AddedLayerLin{FT} <: AbstractLayerLin 
    # Derivatives with respect to (layer) τ, ϖ and Z only
    "Added layer Reflectance matrix R (from + -> -)"
    ṙ⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    ṫ⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    ṙ⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    ṫ⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    J̇₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    J̇₀⁻::AbstractArray{FT,4}
    # Derivatives with respect to all state parameters:
    "Added layer Reflectance matrix R (from + -> -)"
    ap_ṙ⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    ap_ṫ⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    ap_ṙ⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    ap_ṫ⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    ap_J̇₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    ap_J̇₀⁻::AbstractArray{FT,4}
end

# Multisensor Composite layers 
# Elastic
#=
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct CompositeLayerMS{M} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    topR⁻⁺::M#AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    topR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    topT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    topT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    topJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    topJ₀⁻::M
    "Composite layer Reflectance matrix R (from + -> -)"
    botR⁻⁺::M
    "Composite layer Reflectance matrix R (from - -> +)"
    botR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    botT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    botT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    botJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    botJ₀⁻::M
end
=#

#=
"Abstract Type for Surface Types" 
abstract type AbstractSurfaceType end

"Lambertian Surface (scalar per band)"
mutable struct LambertianSurfaceScalar{FT} <: AbstractSurfaceType
    "Albedo (scalar)"
    albedo::FT
end

"Defined as Array (has to have the same length as the band!)"
mutable struct LambertianSurfaceSpectrum{FT} <: AbstractSurfaceType
    "Albedo (vector)"
    albedo::AbstractArray{FT,1}
end

"Defined by Legendre polynomial terms as function of spectral grid, which is scaled to [-1,1] (degree derived from length of `a_coeff`)"
struct LambertianSurfaceLegendre{FT} <: AbstractSurfaceType
    "albedo = legendre_coeff[1] * P₀ + legendre_coeff[2]*P₁ + legendre_coeff[3]*P₂ + ... "
    legendre_coeff::AbstractArray{FT,1}
end

"""
    struct AbsorptionParameters

A struct which holds all absorption-related parameters (before any computations)
"""
mutable struct AbsorptionParameters
    "Molecules to use for absorption calculations (`nBand, nMolecules`)"
    molecules::AbstractArray
    "Volume-Mixing Ratios"
    vmr::Dict
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::AbstractBroadeningFunction
    "Complex Error Function to use in Voigt calculations"
    CEF::AbstractComplexErrorFunction
    "Wing cutoff to use in cross-section calculation (cm⁻¹)"
    wing_cutoff::Integer
    "Lookup table type"
    luts::AbstractArray 
end

"""
    struct ScatteringParameters

A struct which holds all scattering-related parameters (before any computations)
"""
mutable struct ScatteringParametersLin{FT<:AbstractFloat}#{FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "List of scattering aerosols and their properties"
    rt_aerosols_lin::Vector{RT_Aerosol_Lin}
    #="Maximum aerosol particle radius for quadrature points/weights (µm)"
    r_max::FT
    "Number of quadrature points for integration of size distribution"
    nquad_radius::Integer
    "Reference wavelength (µm)"
    λ_ref::FT
    "Algorithm to use for fourier decomposition (NAI2/PCW)"
    decomp_type::AbstractFourierDecompositionType=#
end

"""
    struct vSmartMOM_Parameters

A struct which holds all initial model parameters (before any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Lin_Parameters{FT<:AbstractFloat}#{FT<:Union{AbstractFloat, ForwardDiff.Dual}} 

    # absorption group
    "Optional struct that holds all absorption-related parameters"
    absorption_params_lin::Union{AbsorptionParametersLin, Nothing}

    # scattering group
    "Optional struct that holds all aerosol scattering-related parameters"
    scattering_params_lin::Union{ScatteringParametersLin, Nothing}
    
end
=#
"""
    struct vSmartMOM_Lin

A struct which holds all derived model parameters 

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Lin
    #"Struct with all individual parameters"
    #params_lin::vSmartMOM_Lin_Parameters # for example to include information like 
                                         # whether and which gases/aerosols are to 
                                         # be linearized, Nparams, etc.
    "Array to hold cross-sections over entire atmospheric profile"
    τ̇_abs::AbstractArray{AbstractArray} # w.r.t. psurf, mol. conc.
    #"Rayleigh optical thickness"
    #τ̇_rayl::AbstractArray{AbstractArray} # w.r.t. psurf
    #"Aerosol optical thickness"
    τ̇_aer::AbstractArray{AbstractArray} # w.r.t. τ_ref, nᵣ, nᵢ, r₀, σᵣ, z₀, σz        
    "Truncated aerosol optics"
    lin_aerosol_optics::AbstractArray{AbstractArray{linAerosolOptics}}
    #Nparams::Int16 # total number of state parameters (also consider surface parameters)
end
#=
"""
    struct ComputedAtmosphereProperties

A struct which holds (for the entire atmosphere) all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ComputedAtmospherePropertiesLin

    "Absorption optical depth vectors (wavelength dependent)"
    τ̇_λ_all
    "Albedo vectors (wavelength dependent)"
    ϖ̇_λ_all
    "Absorption optical depth scalars (not wavelength dependent)"
    τ̇_all
    "Albedo scalars (not wavelength dependent)"
    ϖ̇_all
    "Combined Z moments (forward)"
    Ż⁺⁺_all
    "Combined Z moments (backward)"
    Ż⁻⁺_all
    #"Maximum dτs"
    #dτ_max_all
    "dτs"
    dτ̇_all
    #"Number of doublings (for all layers)"
    #ndoubl_all
    "dτs (wavelength dependent)"
    dτ̇_λ_all
    "All expk"
    expk_lin_all
    #"Scattering flags"
    #scatter_all
    "Sum of optical thicknesses of all layers above the current layer"
    τ̇_sum_all
    #"elastic (Cabannes) scattering fraction of Rayleigh (Cabannes+Raman) scattering per layer"
    #ϖ_Cabannes_all
    #"Rayleigh fraction of scattering cross section per layer"
    #fscattRayl_all
    #"Scattering interface type for each layer"
    #scattering_interfaces_all
end

# TODO SUNITI: write a function to compute these properties and create this structure
#Base.@kwdef struct RamanAtmosphereProperties
    #"band spectral grid"
    #grid_in
    #"inelastic scattering SSA"
    #ϖ_λ₀λ₁
    #"inelastic scattering index"
    #i_λ₀λ₁
    #"inelastic (vibrational) scattering SSA: split later for each molecule"
    #ϖ_vib_λ₀λ₁
    #"inelastic (vibrational) scattering index: split later for each molecule"
    #i_vib_λ₀λ₁
    #"Greek coefs in Rayleigh calculations" 
    #greek_raman::GreekCoefs
    #"Combined o2 and n2 Z moments for rotational/rovibrational RS  (forward)"
    #Z⁺⁺_RRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for rotational/rovibrational RS  (backward)"
    #Z⁻⁺_RRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for vibrational RS (forward): split later for each molecule"
    #Z⁺⁺_VRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for vibrational RS (backward): split later for each molecule"
    #Z⁻⁺_VRS #same for rotational and rovibrational scattering
#end


"""
    struct ComputedLayerProperties

A struct which holds all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ComputedLayerPropertiesLin

    "Absorption optical depth vector (wavelength dependent)"
    τ̇_λ 
    "Albedo vector (wavelength dependent)"
    ϖ̇_λ 
    "Absorption optical depth scalar (not wavelength dependent)"
    τ̇ 
    "Albedo scalar (not wavelength dependent)"
    ϖ̇  
    "Combined Z moment (forward)"
    Ż⁺⁺ 
    "Combined Z moment (backward)"
    Ż⁻⁺ 
    #"Maximum dτ"
    #dτ_max 
    "dτ"
    dτ̇     
    #"Number of doublings"
    #ndoubl
    "dτ (wavelength dependent)"
    dτ̇_λ 
    "expk"
    expk_lin 
    #"Scattering flag"
    #scatter 
    "Sum of optical thicknesses of all layers above the current layer"
    τ̇_sum
    #"Fraction of scattering caused by Rayleigh"
    #fscattRayl
    #"Elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering"
    #ϖ_Cabannes 
    #"Scattering interface type for current layer"
    #scattering_interface
end
=#
abstract type AbstractOpticalPropertiesLin end

# Core optical Properties COP
Base.@kwdef struct CoreScatteringOpticalPropertiesLin{FT} <:  AbstractOpticalPropertiesLin
    "Absorption optical depth (scalar or wavelength dependent)"
    τ̇::Union{AbstractArray{FT,1}, AbstractArray{FT,2}}#FT3 
    "Single scattering albedo"
    ϖ̇::Union{AbstractArray{FT,1}, AbstractArray{FT,2}}#FT4   
    "Z scattering matrix (forward)"
    Ż⁺⁺::Union{AbstractArray{FT,3}, AbstractArray{FT,4}}#FT5 
    "Z scattering matrix (backward)"
    Ż⁻⁺::Union{AbstractArray{FT,3}, AbstractArray{FT,4}}#FT5
end

Base.@kwdef struct CoreAbsorptionOpticalPropertiesLin{FT} <:  AbstractOpticalPropertiesLin
    "Absorption optical depth (scalar or wavelength dependent)"
    τ̇::Union{AbstractArray{FT,1}, AbstractArray{FT,2}}
end

Base.@kwdef struct UmbrellaCoreScatteringOpticalProperties{FT} <:  AbstractOpticalPropertiesLin
    fwd::CoreScatteringOpticalProperties{FT}
    lin::Union{Nothing, CoreScatteringOpticalPropertiesLin{FT}}
end

Base.@kwdef struct UmbrellaCoreAbsorptionOpticalProperties{FT} <:  AbstractOpticalPropertiesLin
    fwd::CoreAbsorptionOpticalProperties{FT}
    lin::Union{Nothing, CoreAbsorptionOpticalPropertiesLin{FT}}
end
#=
function include_rayl!(combo::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
                    combo_lin::CoreScatteringOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3}, 
                    rayl::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
                    rayl_lin::CoreScatteringOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3})
    combo_lin.τ̇[1,:] = rayl_lin.τ̇[1,:]
    combo_lin.ϖ̇[1,:] .= 0.0
    combo_lin.Ż⁺⁺[1,:] .= 0.0
    combo_lin.Ż⁻⁺[1,:] .= 0.0
    combo_lin;
end

function include_aer!(iaer::Int16,
    combo::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
    combo_lin::CoreScatteringOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3}, 
    aer::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
    aer_lin::CoreScatteringOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3})
    
    iparam0 = 1+7*(iaer-1)
    for i=1:7
        iparam = iparam0+i
        combo_lin.τ̇[iparam,:] = aer_lin.τ̇[i,:]
        combo_lin.ϖ̇[iparam,:] = (aer_lin.τ̇[i,:].*(aer.ϖ-combo.ϖ) .+ 
                                aer.τ.*aer_lin.ϖ̇[i,:])./combo.τ
        combo_lin.Ż⁺⁺[iparam,:] .= ((aer_lin.τ̇[i,:].*aer.ϖ .+ aer.τ.*aer_lin.ϖ̇[i,:]).*
            (aer.Z⁺⁺-combo.Z⁺⁺) +
            aer.τ.*aer.ϖ.*aer_lin.Ż⁺⁺[i,:,:])./(combo.τ*combo.ϖ)
        combo_lin.Ż⁻⁺[iparam,:] .= ((aer_lin.τ̇[i,:].*aer.ϖ .+ aer.τ.*aer_lin.ϖ̇[i,:]).*
            (aer.Z⁻⁺-combo.Z⁻⁺) +
            aer.τ.*aer.ϖ.*aer_lin.Ż⁻⁺[i,:,:])./(combo.τ*combo.ϖ)
    end
    
    combo_lin;
end

function include_gas!(NAer::Int16, igas::Int16,
    combo::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
    combo_lin::CoreScatteringOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3}, 
    gas_lin::CoreAbsorptionOpticalPropertiesLin{ẋFT, ẋFT2, ẋFT3})
    
    iparam0 = 1+7*NAer+(igas-1)*2
    for i=1:2
        iparam = iparam0+i
        combo_lin.τ̇[iparam,:] = gas_lin.τ̇[i,:]
        combo_lin.ϖ̇[iparam,:] = gas_lin.τ̇[i,:].*combo.ϖ./combo.τ
        combo_lin.Ż⁺⁺[iparam,:] .= 0.0

        combo_lin.Ż⁻⁺[iparam,:] .= 0.0
    end
    
    combo_lin;
end
=#


# Adding Core Optical Properties, can have mixed dimensions!
function Base.:+(a::UmbrellaCoreScatteringOpticalProperties,
                 b::UmbrellaCoreScatteringOpticalProperties)

    x, ẋ = a.fwd, a.lin
    y, ẏ = b.fwd, b.lin

    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    yZ⁺⁺ = y.Z⁺⁺
    yZ⁻⁺ = y.Z⁻⁺

    if ẋ==nothing # Rayleigh    
        τ  = x.τ .+ y.τ
        τ̇  = ẏ.τ̇ #vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = (ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇ .- ϖ'.*ẏ.τ̇)./τ'#vcat((ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇ .- ϖ'.*ẋ.τ̇)./τ', 
            #        (ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇ .- ϖ'.*ẏ.τ̇)./τ')
        #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)), nothing : nothing, nothing
        #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)), nothing : nothing, nothing

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)
    
        nμ = size(xZ⁺⁺,1)
        n1 = 0
        n2 = size(ẏ.τ̇,1)

#=
        @show reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*
            reshape(yZ⁺⁺,1, nμ, nμ, 1)  
        @show 1
        @show reshape(y.τ.*y.ϖ, 1, 1, 1, n).*
            reshape(ẏ.Ż⁺⁺, n2, nμ, nμ, 1) 
        @show 2
        
        @show size(τ'.*ϖ̇ .+ τ̇.*ϖ')
        @show 2.1
        @show reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n)
        @show 2.2
        @show size(Z⁺⁺)
        @show 2.3
        @show reshape(Z⁺⁺,1, nμ, nμ, n)
        @show 2.4

        @show reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n).*
            reshape(Z⁺⁺,1, nμ, nμ, n)
        @show 3
        @show (reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*
            reshape(yZ⁺⁺,1, nμ, nμ, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, 1, n).*
            reshape(ẏ.Ż⁺⁺, n2, nμ, nμ, 1) .- 
            reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n).*
            reshape(Z⁺⁺,1, nμ, nμ, n))
        @show 4    
        @show (reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*
            reshape(yZ⁺⁺,1, nμ, nμ, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, 1, n).*
            reshape(ẏ.Ż⁺⁺, n2, nμ, nμ, 1) .- 
            reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n).*
            reshape(Z⁺⁺,1, nμ, nμ, n))./
            reshape(τ.*ϖ, 1, 1, 1, n)
        @show 5
=#
        Ż⁺⁺ = (reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*
            reshape(yZ⁺⁺,1, nμ, nμ, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, 1, n).*
            reshape(ẏ.Ż⁺⁺, n2, nμ, nμ, 1) .- 
            reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n).*
            reshape(Z⁺⁺,1, nμ, nμ, n))./
            reshape(τ.*ϖ, 1, 1, 1, n)

        Ż⁻⁺ = (reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*
            reshape(yZ⁻⁺, 1, nμ, nμ, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, 1, n).*
            reshape(ẏ.Ż⁻⁺, n2, nμ, nμ, 1) .- 
            reshape(τ'.*ϖ̇ .+ τ̇.*ϖ', n2, 1, 1, n).*
            reshape(Z⁻⁺, 1, nμ, nμ, n))./
            reshape(τ.*ϖ, 1, 1, 1, n)

        #=Ż⁺⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁺⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁺⁺,
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁺⁺,1,nμ,nμ,n) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*ẏ.Ż⁺⁺
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁺⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)


        Ż⁻⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁻⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁻⁺,
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁻⁺,1,nμ,nμ,n) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*ẏ.Ż⁻⁺
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁻⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)
        =#
    else
        τ  = x.τ .+ y.τ
        τ̇  = vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = vcat((ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇ .- ϖ'.*ẋ.τ̇)./τ', 
                    (ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇ .- ϖ'.*ẏ.τ̇)./τ')
        #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)), nothing : nothing, nothing
        #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)), nothing : nothing, nothing

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)
    
        nμ = size(xZ⁺⁺,1)
        n1 = size(ẋ.τ̇,1)
        n2 = size(ẏ.τ̇,1)
        Ż⁺⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁺⁺,1,nμ,nμ,1) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*reshape(ẋ.Ż⁺⁺,n1,nμ,nμ,1),
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁺⁺,1,nμ,nμ,1) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*reshape(ẏ.Ż⁺⁺,n2,nμ,nμ,1)
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁺⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)


        Ż⁻⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁻⁺,1,nμ,nμ,1) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*reshape(ẋ.Ż⁻⁺,n1,nμ,nμ,1),
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁻⁺,1,nμ,nμ,1) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*reshape(ẏ.Ż⁻⁺,n2,nμ,nμ,1)
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁻⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)
    end
    return UmbrellaCoreScatteringOpticalProperties(CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), CoreScatteringOpticalPropertiesLin(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺))    
end

function Base.:+(a::UmbrellaCoreScatteringOpticalProperties,
                 b::UmbrellaCoreAbsorptionOpticalProperties)

    x, ẋ = a.fwd, a.lin
    y, ẏ = b.fwd, b.lin

    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    #yZ⁺⁺ = y.Z⁺⁺
    #yZ⁻⁺ = y.Z⁻⁺

    if ẋ==nothing # Rayleigh    
        τ  = x.τ .+ y.τ
        τ̇  = ẏ.τ̇ #vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ #.* x.ϖ 
        wy = zero(wx) #y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = (- ϖ'.*ẏ.τ̇)./τ'#vcat((ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇ .- ϖ'.*ẋ.τ̇)./τ', 
            #        (ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇ .- ϖ'.*ẏ.τ̇)./τ')
        #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)), nothing : nothing, nothing
        #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)), nothing : nothing, nothing

        n = length(w);
        
        Z⁺⁺ = xZ⁺⁺  
        Z⁻⁺ = xZ⁻⁺ 
    
        nμ = size(xZ⁺⁺,1)
        n1 = 0
        n2 = size(ẏ.τ̇,1)
        Ż⁺⁺ = zeros(n2, nμ, nμ, n)
        Ż⁻⁺ = zeros(n2, nμ, nμ, n)

        #=Ż⁺⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁺⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁺⁺,
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁺⁺,1,nμ,nμ,n) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*ẏ.Ż⁺⁺
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁺⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)


        Ż⁻⁺ = (vcat(
            reshape(ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁻⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁻⁺,
            reshape(ẏ.τ̇.*y.ϖ' .+ y.τ'.*ẏ.ϖ̇, n2, 1, 1, n).*reshape(yZ⁻⁺,1,nμ,nμ,n) .+ reshape(y.τ.*y.ϖ,1,1,1,n).*ẏ.Ż⁻⁺
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁻⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)
        =#
    else
        τ  = x.τ .+ y.τ
        τ̇  = vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = zero(wx) #y.τ .* y.ϖ
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = vcat((ẋ.τ̇.*x.ϖ' .+ x.τ'.*ẋ.ϖ̇ .- ϖ'.*ẋ.τ̇)./τ', 
                (- ϖ'.*ẏ.τ̇)./τ')
        #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)), nothing : nothing, nothing
        #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)), nothing : nothing, nothing

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = xZ⁺⁺ 
        Z⁻⁺ = xZ⁻⁺ 
    
        nμ = size(xZ⁺⁺,1)
        n1 = size(ẋ.τ̇,1)
        n2 = size(ẏ.τ̇,1)
        
        Ż⁺⁺ = (vcat(
            reshape(ẋ.τ̇ .* x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁺⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁺⁺,
            zeros(n2, nμ, nμ, n)
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁺⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)

        Ż⁻⁺ = (vcat(
            reshape(ẋ.τ̇ .* x.ϖ' .+ x.τ'.*ẋ.ϖ̇, n1, 1, 1, n).*reshape(xZ⁻⁺,1,nμ,nμ,n) .+ reshape(x.τ.*x.ϖ,1,1,1,n).*ẋ.Ż⁻⁺,
            zeros(n2, nμ, nμ, n)
            ) .- reshape(τ'.*ϖ̇ .+ τ̇.*ϖ' ,n1+n2,1,1,n).*reshape(Z⁻⁺,1,nμ,nμ,n))./reshape(τ.*ϖ,1,1,1,n)
    end
    return UmbrellaCoreScatteringOpticalProperties(CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), CoreScatteringOpticalPropertiesLin(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺))
end

#=
function Base.:+(a::UmbrellaCoreAbsorptionOpticalProperties,
                 b::UmbrellaCoreAbsorptionOpticalProperties)

    x, ẋ = a.fwd, a.lin
    y, ẏ = b.fwd, b.lin

    # gaseous absorption only
    
    τ  = x.τ .+ y.τ
    τ̇ = vcat(ẋ.τ̇, ẏ.τ̇')
            ϖ  =  zeros(size(τ))
            ϖ̇ = vcat(zeros(size(ẋ.τ̇)), zeros(size(ẏ.τ̇)))
            Z⁺⁺ = zeros(1,1,length(τ))
            Z⁻⁺ = zeros(1,1,length(τ))
            Ż⁺⁺ = zeros(size(τ̇,1),1,1,length(τ))
            Ż⁻⁺ = zeros(size(τ̇,1),1,1,length(τ))

    return UmbrellaCoreAbsorptionOpticalProperties[CoreAbsorptionOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), 
            CoreAbsorptionOpticalPropertiesLin(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺)]
end=#


# Concatenate Core Optical Properties, can have mixed dimensions!
function Base.:*(ẋ::CoreScatteringOpticalPropertiesLin, ẏ::CoreScatteringOpticalPropertiesLin) 
    arr_type  = array_type(architecture(ẋ.τ̇))

    ẋ = expandOpticalProperties(ẋ, arr_type);
    ẏ = expandOpticalProperties(ẏ, arr_type);
    CoreScatteringOpticalPropertiesLin([ẋ.τ̇; ẏ.τ̇],
        [ẋ.ϖ̇; ẏ.ϖ̇],
        cat(ẋ.Ż⁺⁺,ẏ.Ż⁺⁺, dims=3), 
        cat(ẋ.Ż⁻⁺,ẏ.Ż⁻⁺, dims=3))
end

#=
function Base.:+( x::CoreScatteringOpticalProperties, y::CoreAbsorptionOpticalProperties ) 
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    #@show size(wx), size(τ)
    ϖ  = (wx) ./ τ
    CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)
end
=#
function Base.:+(a::UmbrellaCoreAbsorptionOpticalProperties,
                 b::UmbrellaCoreScatteringOpticalProperties)
    return b+a
end

#=
function Base.:*( x::FT, y::CoreScatteringOpticalProperties{FT} ) where FT
    CoreScatteringOpticalProperties(y.τ * x, y.ϖ, y.Z⁺⁺, y.Z⁻⁺)
end
=#
#=
# From https://gist.github.com/mcabbott/80ac43cca3bee8f57809155a5240519f
function _repeat(x::AbstractArray, counts::Integer...)
    N = max(ndims(x), length(counts))
    size_y = ntuple(d -> size(x,d) * get(counts, d, 1), N)
    size_x2 = ntuple(d -> isodd(d) ? size(x, 1+d÷2) : 1, 2*N)

    ## version without mutation
    # ignores = ntuple(d -> reshape(Base.OneTo(counts[d]), ntuple(_->1, 2d-1)..., :), length(counts))
    # y = reshape(broadcast(first∘tuple, reshape(x, size_x2), ignores...), size_y)

    # ## version with mutation
    size_y2 = ntuple(d -> isodd(d) ? size(x, 1+d÷2) : get(counts, d÷2, 1), 2*N)
    y = similar(x, size_y)
    reshape(y, size_y2) .= reshape(x, size_x2)
    y
end
=#