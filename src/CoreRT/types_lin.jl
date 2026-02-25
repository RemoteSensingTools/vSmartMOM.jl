"Abstract Type for Layer Ṙ,Ṫ and J̇ matrices"
abstract type AbstractLayerLin end

"""
    CompositeLayerLin{FT} <: AbstractLayerLin

Linearized (Jacobian) counterpart of [`CompositeLayer`](@ref).  Each field
is a 4-D array whose extra (last) dimension spans the number of retrieval
parameters, storing ∂R/∂x, ∂T/∂x, and ∂J/∂x for the accumulated composite
layer.

# Fields
- `Ṙ⁻⁺::AbstractArray{FT,4}`: ∂R⁻⁺/∂x
- `Ṙ⁺⁻::AbstractArray{FT,4}`: ∂R⁺⁻/∂x
- `Ṫ⁺⁺::AbstractArray{FT,4}`: ∂T⁺⁺/∂x
- `Ṫ⁻⁻::AbstractArray{FT,4}`: ∂T⁻⁻/∂x
- `J̇₀⁺::AbstractArray{FT,4}`: ∂J₀⁺/∂x
- `J̇₀⁻::AbstractArray{FT,4}`: ∂J₀⁻/∂x
"""
Base.@kwdef struct CompositeLayerLin{FT} <: AbstractLayerLin 
    "Composite layer Reflectance matrix R (from + -> -)"
    Ṙ⁻⁺::AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    Ṙ⁺⁻::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from + -> +)"
    Ṫ⁺⁺::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from - -> -)"
    Ṫ⁻⁻::AbstractArray{FT,4}
    "Composite layer source matrix J (in + direction)"
    J̇₀⁺::AbstractArray{FT,4}
    "Composite layer source matrix J (in - direction)"
    J̇₀⁻::AbstractArray{FT,4}
end

"""
    AddedLayerLin{FT} <: AbstractLayerLin

Linearized (Jacobian) counterpart of [`AddedLayer`](@ref).  Stores
derivatives of the single-layer reflectance, transmission, and source
matrices with respect to two groups of parameters:

1. **Layer-intrinsic** derivatives (`ṙ`, `ṫ`, `J̇`) w.r.t. the layer's own
   τ, ϖ, and Z.
2. **All-parameter** derivatives (`ap_ṙ`, `ap_ṫ`, `ap_J̇`) w.r.t. the full
   state vector (surface albedo, VMR profiles, aerosol parameters, etc.).

# Fields
- `ṙ⁻⁺`, `ṫ⁺⁺`, `ṙ⁺⁻`, `ṫ⁻⁻`, `J̇₀⁺`, `J̇₀⁻`: layer-intrinsic Jacobians (4-D)
- `ap_ṙ⁻⁺`, `ap_ṫ⁺⁺`, `ap_ṙ⁺⁻`, `ap_ṫ⁻⁻`, `ap_J̇₀⁺`, `ap_J̇₀⁻`: full state-vector Jacobians (4-D)
"""
Base.@kwdef struct AddedLayerLin{FT} <: AbstractLayerLin 
    # Derivatives with respect to (layer) τ, ϖ and Z only
    "Added layer Reflectance matrix R (from + -> -)"
    ṙ⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    ṫ⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    ṙ⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    ṫ⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    J̇₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    J̇₀⁻::AbstractArray{FT,4}
    # Derivatives with respect to all state parameters:
    "Added layer Reflectance matrix R (from + -> -)"
    ap_ṙ⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    ap_ṫ⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    ap_ṙ⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    ap_ṫ⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    ap_J̇₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    ap_J̇₀⁻::AbstractArray{FT,4}
    # Doubling workspace (pre-allocated, reused across layers to avoid per-call allocations)
    "Doubling workspace: linearized geometric-progression reflectance [nμ × nμ × nSpec × Nparams]"
    dbl_gp_refl_lin::Union{AbstractArray{FT,4}, Nothing} = nothing
    "Doubling workspace: linearized T⁺⁺·gp_refl product [nμ × nμ × nSpec × Nparams]"
    dbl_tt_gp_refl_lin::Union{AbstractArray{FT,4}, Nothing} = nothing
    "Doubling workspace: per-parameter beam attenuation derivative [nSpec × Nparams]"
    dbl_ap_expk_lin::Union{AbstractArray{FT,2}, Nothing} = nothing
    "Doubling workspace: forward source temp J₁⁺ [nμ × 1 × nSpec]"
    dbl_J₁⁺::Union{AbstractArray{FT,3}, Nothing} = nothing
    "Doubling workspace: forward source temp J₁⁻ [nμ × 1 × nSpec]"
    dbl_J₁⁻::Union{AbstractArray{FT,3}, Nothing} = nothing
    "Doubling workspace: linearized source temp J̇₁⁺ [nμ × 1 × nSpec × Nparams]"
    dbl_ap_J̇₁⁺::Union{AbstractArray{FT,4}, Nothing} = nothing
    "Doubling workspace: linearized source temp J̇₁⁻ [nμ × 1 × nSpec × Nparams]"
    dbl_ap_J̇₁⁻::Union{AbstractArray{FT,4}, Nothing} = nothing
    "Doubling workspace: geometric progression (forward) [nμ × nμ × nSpec]"
    dbl_gp_refl::Union{AbstractArray{FT,3}, Nothing} = nothing
    "Doubling workspace: T⁺⁺·gp_refl (forward) [nμ × nμ × nSpec]"
    dbl_tt_gp_refl::Union{AbstractArray{FT,3}, Nothing} = nothing
end

"""
    struct vSmartMOM_Lin{A,B,C}

Holds linearized (Jacobian) model parameters: derivatives of optical depths
and aerosol properties w.r.t. physical state-vector elements.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Lin{A,B,C}
    "∂τ_abs/∂x per band: Vector of arrays [NGas × nSpec × nLayers]"
    τ̇_abs::A
    "∂τ_aer/∂x per band: Vector of arrays [NAer × 7 × nSpec × nLayers]"
    τ̇_aer::B
    "Linearized aerosol optics per band per aerosol: Vector{Vector{linAerosolOptics}}"
    lin_aerosol_optics::C
end
abstract type AbstractOpticalPropertiesLin end

"""
    CoreScatteringOpticalPropertiesLin{T1,T2,T3}

Per-layer Jacobian of the four core optical properties `(τ, ϖ, Z⁺⁺, Z⁻⁺)`
with respect to the retrieval state vector `x`.

This struct is the **AD boundary**: everything upstream of it (Mie code,
absorption cross-sections, atmospheric profiles) may use `ForwardDiff.Dual`
numbers, but by the time values reach the RT kernels they are extracted
into plain `Float64`/`Float32` arrays stored here.

The chain rule `lin_added_layer_all_params!` maps these derivatives into
the full `∂R/∂x` via:

```math
\\frac{\\partial R}{\\partial x_j} =
  \\frac{\\partial R}{\\partial \\tau} \\frac{\\partial \\tau}{\\partial x_j} +
  \\frac{\\partial R}{\\partial \\varpi} \\frac{\\partial \\varpi}{\\partial x_j} +
  \\frac{\\partial R}{\\partial \\mathbf{Z}} \\frac{\\partial \\mathbf{Z}}{\\partial x_j}
```

See also: [`OpticalPropertyJacobian`](@ref) (alias).
"""
Base.@kwdef struct CoreScatteringOpticalPropertiesLin{T1,T2,T3} <: AbstractOpticalPropertiesLin
    "∂τ/∂x — [Nparams] or [nSpec × Nparams]"
    τ̇::T1
    "∂ϖ/∂x — [Nparams] or [nSpec × Nparams]"
    ϖ̇::T2
    "∂Z⁺⁺/∂x — [nμ × nμ × nSpec] or [nμ × nμ × nSpec × Nparams]"
    Ż⁺⁺::T3
    "∂Z⁻⁺/∂x — [nμ × nμ × nSpec] or [nμ × nμ × nSpec × Nparams]"
    Ż⁻⁺::T3
end

"""
    OpticalPropertyJacobian

Type alias for [`CoreScatteringOpticalPropertiesLin`](@ref).  Use this name
in new code to emphasise the physical meaning: the Jacobian of the four
core optical properties `(τ, ϖ, Z⁺⁺, Z⁻⁺)` w.r.t. the state vector.
"""
const OpticalPropertyJacobian = CoreScatteringOpticalPropertiesLin

Base.@kwdef struct CoreAbsorptionOpticalPropertiesLin{T1} <: AbstractOpticalPropertiesLin
    "∂τ/∂x — [Nparams] or [Nparams × nSpec]"
    τ̇::T1
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
When both have derivatives, they are `hcat`/`cat`-ed along the parameter dimension.
"""
function Base.:+(a::UmbrellaCoreScatteringOpticalProperties,
                 b::UmbrellaCoreScatteringOpticalProperties)

    x, ẋ = a.fwd, a.lin
    y, ẏ = b.fwd, b.lin

    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    yZ⁺⁺ = y.Z⁺⁺
    yZ⁻⁺ = y.Z⁻⁺

    if ẋ==nothing # Rayleigh    
        τ  = x.τ .+ y.τ
        τ̇  = ẏ.τ̇ #vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = (ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇ .- ϖ.*ẏ.τ̇)./τ#vcat((ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇ .- ϖ.*ẋ.τ̇)./τ, 

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)
    
        nμ = size(xZ⁺⁺,1)
        n1 = 0
        n2 = size(ẏ.τ̇,2)

        Ż⁺⁺ = (reshape(ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇, 1, 1, n, n2).*
            reshape(yZ⁺⁺, nμ, nμ, 1, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, n, 1).*
            reshape(ẏ.Ż⁺⁺, nμ, nμ, 1, n2) .- 
            reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n2).*
            reshape(Z⁺⁺, nμ, nμ, n, 1))./
            reshape(τ.*ϖ, 1, 1, n, 1)

        Ż⁻⁺ = (reshape(ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇, 1, 1, n, n2).*
            reshape(yZ⁻⁺, nμ, nμ, 1, 1) .+ 
            reshape(y.τ.*y.ϖ, 1, 1, n, 1).*
            reshape(ẏ.Ż⁻⁺, nμ, nμ, 1, n2) .- 
            reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n2).*
            reshape(Z⁻⁺, nμ, nμ, n, 1))./
            reshape(τ.*ϖ, 1, 1, n, 1)

    else
        τ  = x.τ .+ y.τ
        τ̇  = hcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = hcat((ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇ .- ϖ.*ẋ.τ̇)./τ, 
                    (ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇ .- ϖ.*ẏ.τ̇)./τ)

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)
    
        nμ = size(xZ⁺⁺,1)
        n1 = size(ẋ.τ̇,2)
        n2 = size(ẏ.τ̇,2)
        Ż⁺⁺ = (cat(
            reshape(ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇, 1, 1, n, n1).*reshape(xZ⁺⁺,nμ,nμ,1,1) .+ reshape(x.τ.*x.ϖ,1,1,n,1).*reshape(ẋ.Ż⁺⁺,nμ,nμ,1,n1),
            reshape(ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇, 1, 1, n, n2).*reshape(yZ⁺⁺,nμ,nμ,1,1) .+ reshape(y.τ.*y.ϖ,1,1,n,1).*reshape(ẏ.Ż⁺⁺,nμ,nμ,1,n2),
                dims=4) .- reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n1+n2).*reshape(Z⁺⁺,nμ,nμ,n,1))./reshape(τ.*ϖ,1,1,n,1)


        Ż⁻⁺ = (cat(
            reshape(ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇, 1, 1, n, n1).*reshape(xZ⁻⁺,nμ,nμ,1,1) .+ reshape(x.τ.*x.ϖ,1,1,n,1).*reshape(ẋ.Ż⁻⁺,nμ,nμ,1,n1),
            reshape(ẏ.τ̇.*y.ϖ .+ y.τ.*ẏ.ϖ̇, 1, 1, n, n2).*reshape(yZ⁻⁺,nμ,nμ,1,1) .+ reshape(y.τ.*y.ϖ,1,1,n,1).*reshape(ẏ.Ż⁻⁺,nμ,nμ,1,n2),
                dims=4) .- reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n1+n2).*reshape(Z⁻⁺,nμ,nμ,n,1))./reshape(τ.*ϖ,1,1,n,1)
    end
    return UmbrellaCoreScatteringOpticalProperties(CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), CoreScatteringOpticalPropertiesLin(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺))    
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

    x, ẋ = a.fwd, a.lin
    y, ẏ = b.fwd, b.lin

    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺

    if ẋ==nothing # Rayleigh    
        τ  = x.τ .+ y.τ
        τ̇  = ẏ.τ̇ #vcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ #.* x.ϖ 
        wy = zero(wx) #y.τ .* y.ϖ  
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = (- ϖ.*ẏ.τ̇)./τ#vcat((ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇ .- ϖ.*ẋ.τ̇)./τ, 

        n = length(w);
        
        Z⁺⁺ = xZ⁺⁺  
        Z⁻⁺ = xZ⁻⁺ 
    
        nμ = size(xZ⁺⁺,1)
        n1 = 0
        n2 = size(ẏ.τ̇,2)
        Ż⁺⁺ = zeros(nμ, nμ, n, n2)
        Ż⁻⁺ = zeros(nμ, nμ, n, n2)

    else
        τ  = x.τ .+ y.τ
        τ̇  = hcat(ẋ.τ̇, ẏ.τ̇)
        wx = x.τ .* x.ϖ 
        wy = zero(wx) #y.τ .* y.ϖ
        w  = wx .+ wy
        ϖ  =  w ./ τ

        ϖ̇ = hcat((ẋ.τ̇.*x.ϖ .+ x.τ.*ẋ.ϖ̇ .- ϖ.*ẋ.τ̇)./τ, 
                (- ϖ.*ẏ.τ̇)./τ)

        n = length(w);
        
        wy = wy ./ w
        wx = wx ./ w
        wx = reshape(wx,1,1,n)
        wy = reshape(wy,1,1,n)
        
        Z⁺⁺ = xZ⁺⁺ 
        Z⁻⁺ = xZ⁻⁺ 
    
        nμ = size(xZ⁺⁺,1)
        n1 = size(ẋ.τ̇,2)
        n2 = size(ẏ.τ̇,2)
        
        Ż⁺⁺ = (cat(
            reshape(ẋ.τ̇ .* x.ϖ .+ x.τ.*ẋ.ϖ̇, 1, 1, n, n1).*reshape(xZ⁺⁺,nμ,nμ,n,1) .+ reshape(x.τ.*x.ϖ,1,1,n,1).*ẋ.Ż⁺⁺,
            zeros(nμ, nμ, n, n2),
                dims=4) .- reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n1+n2).*reshape(Z⁺⁺,nμ,nμ,n,1))./reshape(τ.*ϖ,1,1,n,1)

        Ż⁻⁺ = (cat(
            reshape(ẋ.τ̇ .* x.ϖ .+ x.τ.*ẋ.ϖ̇, 1, 1, n, n1).*reshape(xZ⁻⁺,nμ,nμ,n,1) .+ reshape(x.τ.*x.ϖ,1,1,n,1).*ẋ.Ż⁻⁺,
            zeros(nμ, nμ, n, n2),
                dims=4) .- reshape(τ.*ϖ̇ .+ τ̇.*ϖ, 1, 1, n, n1+n2).*reshape(Z⁻⁺,nμ,nμ,n,1))./reshape(τ.*ϖ,1,1,n,1)
    end
    return UmbrellaCoreScatteringOpticalProperties(CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), CoreScatteringOpticalPropertiesLin(τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺))
end



# Concatenate Core Optical Properties, can have mixed dimensions!
function Base.:*(ẋ::CoreScatteringOpticalPropertiesLin, ẏ::CoreScatteringOpticalPropertiesLin) 
    arr_type  = array_type(architecture(ẋ.τ̇))

    ẋ = expandOpticalProperties(ẋ, arr_type);
    ẏ = expandOpticalProperties(ẏ, arr_type);
    CoreScatteringOpticalPropertiesLin([ẋ.τ̇; ẏ.τ̇],
        [ẋ.ϖ̇; ẏ.ϖ̇],
        cat(ẋ.Ż⁺⁺,ẏ.Ż⁺⁺, dims=3), 
        cat(ẋ.Ż⁻⁺,ẏ.Ż⁻⁺, dims=3))
end

function Base.:+(a::UmbrellaCoreAbsorptionOpticalProperties,
                 b::UmbrellaCoreScatteringOpticalProperties)
    return b+a
end
