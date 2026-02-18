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

Base.@kwdef struct UmbrellaCoreScatteringOpticalProperties{FWD <: CoreScatteringOpticalProperties, LIN} <:  AbstractOpticalPropertiesLin
    fwd::FWD
    lin::LIN  # Union{Nothing, CoreScatteringOpticalPropertiesLin}
end

Base.@kwdef struct UmbrellaCoreAbsorptionOpticalProperties{FWD <: CoreAbsorptionOpticalProperties, LIN} <:  AbstractOpticalPropertiesLin
    fwd::FWD
    lin::LIN  # Union{Nothing, CoreAbsorptionOpticalPropertiesLin}
end


# Adding Core Optical Properties, can have mixed dimensions!
"""
    Base.:+(a::UmbrellaCoreScatteringOpticalProperties, b::UmbrellaCoreScatteringOpticalProperties)

Combine two scattering layers' optical properties and propagate linearized derivatives.

When combining Rayleigh (with `lin=nothing`) and aerosol optical properties, or two
scattering components, this operator computes the mixed effective properties:

```math
\\tau = \\tau_a + \\tau_b, \\quad
\\varpi = \\frac{\\tau_a \\varpi_a + \\tau_b \\varpi_b}{\\tau}, \\quad
\\mathbf{Z} = \\frac{\\tau_a \\varpi_a \\mathbf{Z}_a + \\tau_b \\varpi_b \\mathbf{Z}_b}{\\tau \\varpi}
```

The derivatives are propagated via the quotient/product rule. When `a.lin === nothing`
(e.g., Rayleigh without derivatives), only the derivatives from `b` are retained.
When both have derivatives, they are `vcat`-ed along the parameter dimension.
"""
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

"""
    Base.:+(a::UmbrellaCoreScatteringOpticalProperties, b::UmbrellaCoreAbsorptionOpticalProperties)

Add gas absorption to combined scattering optical properties and propagate derivatives.

Gas absorption only contributes to the total optical depth (no scattering):
```math
\\tau = \\tau_\\text{scat} + \\tau_\\text{gas}, \\quad
\\varpi = \\frac{\\tau_\\text{scat} \\varpi_\\text{scat}}{\\tau}, \\quad
\\mathbf{Z} = \\mathbf{Z}_\\text{scat}
```

The gas VMR derivative enters ``\\dot{\\varpi}`` as ``-\\varpi \\dot{\\tau}_\\text{gas}/\\tau``
(diluting the scattering fraction), while ``\\dot{\\mathbf{Z}}_\\text{gas} = 0``.
"""
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

function Base.:+(a::UmbrellaCoreAbsorptionOpticalProperties,
                 b::UmbrellaCoreScatteringOpticalProperties)
    return b+a
end
