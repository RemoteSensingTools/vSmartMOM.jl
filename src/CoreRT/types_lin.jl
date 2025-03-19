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

"""
    struct vSmartMOM_Model

A struct which holds all derived model parameters (including any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Model_Lin
    "Struct with all individual parameters"
    params_lin::vSmartMOM_Lin_Parameters # for example to include information like 
                                         # whether and which gases/aerosols are to 
                                         # be linearized, Nparams, etc.
    
    "Truncated aerosol optics"
    aerosol_optics_lin::AbstractArray{AbstractArray{AerosolOpticsLin}}
    
    "Array to hold cross-sections over entire atmospheric profile"
    τ̇_abs::AbstractArray{AbstractArray} # w.r.t. psurf, mol. conc.
    "Rayleigh optical thickness"
    τ̇_rayl::AbstractArray{AbstractArray} # w.r.t. psurf
    "Aerosol optical thickness"
    τ̇_aer::AbstractArray{AbstractArray} # w.r.t. τ_ref, nᵣ, nᵢ, r₀, σᵣ, z₀, σz        

    Nparams::Int16 # total number of state parameters (also consider surface parameters)
end

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

abstract type AbstractOpticalProperties end

# Core optical Properties COP
Base.@kwdef struct CoreScatteringOpticalPropertiesLin{FT,FT2,FT3} <:  AbstractOpticalPropertiesLin
    "Absorption optical depth (scalar or wavelength dependent)"
    τ̇::FT 
    "Single scattering albedo"
    ϖ̇::FT2   
    "Z scattering matrix (forward)"
    Ż⁺⁺::FT3 
    "Z scattering matrix (backward)"
    Ż⁻⁺::FT3
end

Base.@kwdef struct CoreAbsorptionOpticalPropertiesLin{FT} <:  AbstractOpticalPropertiesLin
    "Absorption optical depth (scalar or wavelength dependent)"
    τ̇::FT 
end

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



# Adding Core Optical Properties, can have mixed dimensions!
function Base.:+( x::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
                  y::CoreScatteringOpticalProperties{yFT, yFT2, yFT3} 
                ) where {xFT, xFT2, xFT3, yFT, yFT2, yFT3} 
    # Predefine some arrays:            
    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    yZ⁺⁺ = y.Z⁺⁺
    yZ⁻⁺ = y.Z⁻⁺

    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    wy = y.τ .* y.ϖ  
    w  = wx .+ wy
    ϖ  =  w ./ τ
    #@show xFT, xFT2, xFT3
    all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)) : nothing
    all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)) : nothing

    # A bit more tedious for Z matrices:
    #@show size(wx), size(wy), size(w), wy, length(unique(wx))
    #length(unique(w))  == 1 ? w  = unique(w) : nothing
    #length(unique(wx)) == 1 ? wx = unique(wx) : nothing
    #length(unique(wy)) == 1 ? wy = unique(wy) : nothing
    n = length(w);
    #@show size(wx), size(wy), size(w)
    #@show n
    #if xFT <: AbstractFloat && yFT <: AbstractFloat
    #    Z⁺⁺ = ((wx .* xZ⁺⁺) .+ (wy .* yZ⁺⁺)) ./ w 
    #    Z⁻⁺ = ((wx .* xZ⁻⁺) .+ (wy .* yZ⁻⁺)) ./ w
    #else
        #@show xFT, yFT, length(wx), length(wy), length(w)
        #!(xFT <: AbstractFloat) ? wx = reshape(wx,1,1,n) : nothing
        #!(yFT <: AbstractFloat) ? wy = reshape(wy,1,1,n) : nothing
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        #w = reshape(w,1,1,n)
        #println("Block!")
        #@show size(wx), size(wy), size(xZ⁺⁺), size(yZ⁺⁺)
        #(xFT3 <: Matrix)  ? xZ⁺⁺ = _repeat(xZ⁺⁺,1,1,n) : nothing
        #(xFT3 <: Matrix)  ? xZ⁻⁺ = _repeat(xZ⁻⁺,1,1,n) : nothing
        #(yFT3 <: Matrix)  ? yZ⁺⁺ = _repeat(yZ⁺⁺,1,1,n) : nothing
        #(yFT3 <: Matrix)  ? yZ⁻⁺ = _repeat(yZ⁻⁺,1,1,n) : nothing
        #@show size(wx), size(xZ⁺⁺), size(wy), size(yZ⁺⁺)
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)
        #@show size(Z⁺⁺), size(Z⁻⁺)
        #println("###")
        
    #end
    CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)  
end

# Concatenate Core Optical Properties, can have mixed dimensions!
function Base.:*( x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties ) 
    arr_type  = array_type(architecture(x.τ))

    x = expandOpticalProperties(x,arr_type);
    y = expandOpticalProperties(y,arr_type);
    CoreScatteringOpticalProperties([x.τ; y.τ],[x.ϖ; y.ϖ],cat(x.Z⁺⁺,y.Z⁺⁺, dims=3), cat(x.Z⁻⁺,y.Z⁻⁺, dims=3))
end

function Base.:+( x::CoreScatteringOpticalProperties, y::CoreAbsorptionOpticalProperties ) 
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    #@show size(wx), size(τ)
    ϖ  = (wx) ./ τ
    CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)
end

function Base.:+(  y::CoreAbsorptionOpticalProperties, x::CoreScatteringOpticalProperties ) 
    x + y
end


function Base.:*( x::FT, y::CoreScatteringOpticalProperties{FT} ) where FT
    CoreScatteringOpticalProperties(y.τ * x, y.ϖ, y.Z⁺⁺, y.Z⁻⁺)
end

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