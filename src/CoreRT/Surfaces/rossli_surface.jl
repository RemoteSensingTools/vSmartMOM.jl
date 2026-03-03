"""
    reflectance(RossLi::RossLiSurfaceScalar, n, μᵢ, μᵣ, dϕ)

Ross-Li BRDF model (Ross & Li, 2000): linear combination of isotropic, volumetric, and geometric kernels.

``\\rho = f_{\\text{iso}} K_{\\text{iso}} + f_{\\text{vol}} K_{\\text{vol}} + f_{\\text{geo}} K_{\\text{geo}}``

- **fiso**: Isotropic coefficient.
- **fvol**: Volumetric (Ross Thick) coefficient — models dense canopy scattering.
- **fgeo**: Geometric (Li Sparse) coefficient — models shadowing effects.

For scalar RT (n=1) only; returns zero for polarized. See Ross & Li (2000), RSE.
"""
function reflectance(RossLi::RossLiSurfaceScalar{FT},n, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    (; fiso, fvol, fgeo) = RossLi
    # Function was defined for RAMI definition, have to reverse here:
    dϕ = π - dϕ
    # Convert cosines to angles for Ross-Li kernels
    θᵢ   = acos(μᵢ) #assert 0<=θᵢ<=π/2
    θᵣ   = acos(μᵣ) #assert 0<=θᵣ<=π/2

    if n==1
        return fiso * RossLi_K_iso() + 
            fvol * RossLi_K_vol(θᵢ, θᵣ, dϕ) + 
            fgeo * RossLi_K_geo(θᵢ, θᵣ, dϕ)
    else
        return FT(0)
    end
end
"""Isotropic kernel: ``K_{\\text{iso}} = 1``."""
function RossLi_K_iso()
    return 1.0;
end

"""
    RossLi_K_vol(θᵢ, θᵣ, dϕ)

Ross Thick volumetric kernel: models dense vegetation scattering.
``K_{\\text{vol}} = ((\\pi/2 - \\xi)\\cos\\xi + \\sin\\xi)/(\\cos\\theta_i + \\cos\\theta_r) - \\pi/4``
"""
function RossLi_K_vol(θᵢ::FT, θᵣ::FT, dϕ::FT) where FT
    ξ = acos(cos(θᵢ)*cos(θᵣ) + sin(θᵢ)*sin(θᵣ)*cos(dϕ))
    #@show ξ
    K_vol = ((π/FT(2) - ξ)*cos(ξ) + sin(ξ))/(cos(θᵢ)+cos(θᵣ)) - (π/FT(4))
    #@show K_vol
    return K_vol
end

"""
    RossLi_K_geo(θᵢ, θᵣ, dϕ)

Li Sparse geometric kernel: models shadowing by sparsely distributed objects.
Uses RAMI default: h/b = 2, b/r = 1. Returns overlap area ``O``.
"""
function RossLi_K_geo(θᵢ::FT, θᵣ::FT, dϕ::FT) where FT
    h_by_b = FT(2) #for RAMI
    b_by_r = FT(1) #for RAMI

    θᵢᵖ = atan(tan(θᵢ)*b_by_r)
    θᵣᵖ = atan(tan(θᵣ)*b_by_r)
    #@show θᵢᵖ, θᵣᵖ
    ξᵖ = acos(cos(θᵢᵖ)*cos(θᵣᵖ) + sin(θᵢᵖ)*sin(θᵣᵖ)*cos(dϕ))
    #@show ξᵖ
    D = sqrt(tan(θᵢᵖ)^2 + tan(θᵣᵖ)^2 - 2tan(θᵢᵖ)*tan(θᵣᵖ)*cos(dϕ))
    #@show  D
    ct = h_by_b * sqrt(D^2 + (tan(θᵢᵖ)*tan(θᵣᵖ)*sin(dϕ))^2)/(sec(θᵢᵖ)+sec(θᵣᵖ))
    if ct>FT(1) 
        ct=FT(1)
    elseif ct<FT(-1) 
        ct=FT(-1)
    end
    t = acos(ct)
    
    #@show  t
    _O = (1/π)*(t-sin(t)*cos(t))*(sec(θᵢᵖ)+sec(θᵣᵖ))
    #@show  _O 
    _O = _O - (sec(θᵢᵖ)+sec(θᵣᵖ)) + 0.5*(1+cos(ξᵖ))*sec(θᵢᵖ)*sec(θᵣᵖ)
    #@show  _O 
    return _O
end