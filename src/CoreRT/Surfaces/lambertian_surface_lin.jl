#=

This file specifies how to create surface layers and their linearized properties for
Lambertian surfaces in the linearized RT framework.

The Lambertian surface is characterized by a single albedo parameter ``\\alpha``. For the
0th Fourier moment (``m=0``), the surface reflection matrix is:
```math
\\mathbf{R}_\\text{surf} = 2\\alpha \\, \\mathbf{E} \\, \\mathbf{w}^T
```
where ``\\mathbf{E}`` is the identity-like Stokes vector selector and ``\\mathbf{w}`` are
the quadrature weights. For ``m > 0``, the surface contribution vanishes.

The derivative with respect to surface albedo is:
```math
\\frac{\\partial \\mathbf{R}_\\text{surf}}{\\partial \\alpha} = 2 \\, \\mathbf{E} \\, \\mathbf{w}^T
```

This is placed at parameter index `iparam = Nparams` (the last parameter in the state vector).

=#

"""
    create_surface_layer!(RS_type, lambertian, added_layer, added_layer_lin, iparam, 
                          SFI, m, pol_type, quad_points, τ_sum, τ̇_sum, architecture)

Compute surface reflection/transmission matrices and their derivatives for a scalar
Lambertian surface.

Sets the `added_layer` forward matrices and `added_layer_lin` derivative matrices for
the surface "layer", including the source function contribution from the solar beam
attenuated through the full atmosphere.

# Arguments
- `RS_type::noRS`: No Raman scattering.
- `lambertian::LambertianSurfaceScalar`: Surface with scalar albedo ``\\alpha``.
- `added_layer::AddedLayer`: Output forward matrices (modified in-place).
- `added_layer_lin::AddedLayerLin`: Output linearized matrices (modified in-place).
- `iparam::Int`: Parameter index for the surface albedo derivative.
- `SFI::Bool`: Source Function Integration flag.
- `m::Int`: Fourier moment (only ``m=0`` has nonzero surface contribution).
- `pol_type`: Polarization type.
- `quad_points`: Quadrature points and weights.
- `τ_sum`: Total optical depth from TOA to surface `[nSpec]`.
- `τ̇_sum`: Derivative of total τ `[nSpec × Nparams]`.
- `architecture`: CPU or GPU.
""" 
function create_surface_layer!(RS_type::noRS, 
                            lambertian::LambertianSurfaceScalar{FT}, 
                            #lambertian_lin::LambertianSurfaceScalarLin{FT}, 
                            added_layer::AddedLayer,
                            added_layer_lin::AddedLayerLin,
                            iparam::Int,
                            SFI,
                            m::Int,
                            pol_type,
                            quad_points,
                            τ_sum, τ̇_sum, F₀,
                            architecture) where {FT}
    
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
    (; n) = pol_type
    arr_type = array_type(architecture)
    
    nparams = size(τ̇_sum,2) # nparams ≠ Nparams (in general)
    nspec = length(τ_sum)
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    tmp = arr_type(ones(pol_type.n*Nquad))
    T_surf = Diagonal(tmp)
    i₀ = iμ₀Nstart:iμ₀Nstart+n-1
    # Surface source Jacobians are Fourier-local work arrays.  In particular,
    # Lambertian m>0 contributes exactly zero; clear the complete state-vector
    # slabs here so values from m=0 (or uninitialized storage) cannot leak into
    # later moments.
    added_layer_lin.ap_J̇₀⁺ .= zero(FT)
    added_layer_lin.ap_J̇₀⁻ .= zero(FT)
    if m == 0
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        ρ = FT(2) * lambertian.albedo#/FT(π)
        
        R_surf = Matrix(Diagonal(vcat(ρ, zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)
        Ṙ_surf = Matrix(Diagonal(vcat(FT(2), zeros(FT,pol_type.n-1))))
        Ṙ_surf = repeat(Ṙ_surf',Nquad)
        Ṙ_surf = repeat(Ṙ_surf',Nquad)
        #Ṙ_surf = Array{FT}(undef, Nparams, size(R_surf,1), size(R_surf,2))
        #Ṙ_surf[end,:,:] .= Ṙ_surf_tmp
        
        
        # Move to architecture:
        R_surf = arr_type(R_surf)
        Ṙ_surf = arr_type(Ṙ_surf)
        
        # Source function of surface:
        if SFI
            unweight = FT(2π) #this is multiplied to all non-solar, isotropic source functions to exclude them from the azimuthal weighting applied in the postprocessing step
            if quad_points.external_solar
                direct = _direct_solar_at_surface(
                    F₀, FT, pol_type, quad_points, τ_sum, architecture)
                incident_I = reshape(@view(direct[1, :]), 1, :)
                rows_I = 1:n:size(added_layer.r⁻⁺, 1)
                added_layer.j₀⁺ .= zero(FT)
                added_layer.j₀⁻ .= zero(FT)
                @views added_layer.j₀⁻[rows_I,1,:] .=
                    μ₀ .* FT(2) .* lambertian.albedo .* incident_I
                for ctr in 1:nparams
                    @views added_layer_lin.ap_J̇₀⁻[rows_I,1,:,ctr] .=
                        -added_layer.j₀⁻[rows_I,1,:] .*
                        reshape(τ̇_sum[:,ctr], 1, :) ./ μ₀
                end
                @views added_layer_lin.ap_J̇₀⁻[rows_I,1,:,iparam] .=
                    μ₀ .* FT(2) .* incident_I
            else
                F₀_NquadN = arr_type(zeros(length(qp_μN),length(τ_sum)))
                Ḟ₀_NquadN = arr_type(zeros(length(qp_μN),length(τ_sum),nparams+1))
                tmpF = F₀ .* arr_type(exp.(-τ_sum/μ₀))'
                F₀_NquadN[i₀,:] .= tmpF
                Ḟ₀_NquadN[i₀,:,1:nparams] .= -reshape(tmpF,n,nspec,1) .* reshape(τ̇_sum, 1, nspec, nparams) / μ₀
                added_layer.j₀⁺ .= zero(FT)
                added_layer.j₀⁻[:,1,:] .= μ₀ * (R_surf * F₀_NquadN)
                for ii=1:nspec
                    for ctr=1:nparams
                        added_layer_lin.ap_J̇₀⁻[:,1,ii,ctr] .= μ₀ * R_surf * Ḟ₀_NquadN[:,ii,ctr]
                    end
                    added_layer_lin.ap_J̇₀⁻[:,1,ii,iparam] .= μ₀ * Ṙ_surf * F₀_NquadN[:,ii]
                end
            end
        # for SIF
            #reinstate the following line after linearization works
            #added_layer.j₀⁻[:,1,:] .+= (1/π)*repeat(arr_type(RS_type.SIF₀),Nquad) * unweight
        end

        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        Ṙ_surf[:,:] = Ṙ_surf[:,:] * Diagonal(qp_μN.*wt_μN)
        
        #R_surf = 2R_surf * Diagonal(qp_μN.*wt_μN)
        #R_surf = R_surf * Diagonal(qp_μN.*wt_μN)/π

        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;#1. #0.0; #T_surf;
        added_layer.t⁻⁻ .= zero(FT); #T_surf;

        added_layer_lin.ap_ṙ⁻⁺[:,:,:,iparam] .= Ṙ_surf;
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT);
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT);#1. #0.0; #T_surf;
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT); #T_surf;

    else
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= zero(FT); #T_surf;
        added_layer.j₀⁺ .= zero(FT);
        added_layer.j₀⁻ .= zero(FT);

        added_layer_lin.ap_ṙ⁻⁺ .= zero(FT);
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT);
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT);#1. #0.0; #T_surf;
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT);
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        added_layer_lin.ap_J̇₀⁻ .= zero(FT)
    end
    synchronize_if_gpu()
end

"Legendre basis used by `LambertianSurfaceLegendre`, with columns P₀…Pₙ."
function _lambertian_legendre_basis(surface::LambertianSurfaceLegendre{FT},
                                    nSpec::Int) where {FT}
    ncoeff = length(surface.legendre_coeff)
    ncoeff > 0 || throw(ArgumentError(
        "LambertianSurfaceLegendre requires at least one coefficient"))
    ncoeff == 1 && return ones(FT, nSpec, 1)
    x = collect(range(FT(-1), FT(1), length=nSpec))
    return Scattering.compute_legendre_poly(x, ncoeff)[1]
end

"""
Linearized Lambertian-Legendre surface. Each coefficient is one consecutive
surface state-vector element, so `iparams[k]` stores ∂y/∂ρₖ₋₁.
"""
function create_surface_layer!(RS_type::noRS,
                               surface::LambertianSurfaceLegendre{FT},
                               added_layer::AddedLayer,
                               added_layer_lin::AddedLayerLin,
                               iparams::AbstractUnitRange{Int},
                               SFI, m::Int, pol_type, quad_points,
                               τ_sum, τ̇_sum, F₀, architecture) where {FT}
    nSpec = length(τ_sum)
    ncoeff = length(surface.legendre_coeff)
    length(iparams) == ncoeff || throw(DimensionMismatch(
        "LambertianSurfaceLegendre has $ncoeff coefficients but received " *
        "$(length(iparams)) Jacobian columns"))

    # Use the common forward scaffold, then populate its analytic tangent.
    create_surface_layer!(surface, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture; F₀)
    added_layer_lin.ap_ṙ⁻⁺ .= zero(FT)
    added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
    added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
    added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)
    added_layer_lin.ap_J̇₀⁺ .= zero(FT)
    added_layer_lin.ap_J̇₀⁻ .= zero(FT)
    m == 0 || return synchronize_if_gpu()

    arr_type = array_type(architecture)
    (; qp_μN, wt_μN, μ₀, iμ₀Nstart) = quad_points
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
    basis = arr_type(_lambertian_legendre_basis(surface, nSpec))
    albedo = arr_type(surface_albedo(surface, τ_sum))

    R_unit = Matrix(Diagonal(vcat(one(FT), zeros(FT, pol_type.n - 1))))
    R_unit = repeat(R_unit', Nquad)
    R_unit = arr_type(repeat(R_unit', Nquad))
    R_unit_weighted = R_unit .* reshape(qp_μN .* wt_μN, 1, :)

    # Atmospheric columns enter the reflected direct beam through attenuation.
    if SFI
        incident_I = if quad_points.external_solar
            direct = _direct_solar_at_surface(
                F₀, FT, pol_type, quad_points, τ_sum, architecture)
            reshape(@view(direct[1, :]), 1, :)
        else
            beam = _surface_beam_at_surface(F₀, FT, pol_type, quad_points,
                                            τ_sum, architecture)
            reshape(@view(beam[iμ₀Nstart, :]), 1, :)
        end
        rows_I = 1:pol_type.n:size(added_layer.r⁻⁺, 1)
        for p in axes(τ̇_sum, 2)
            @views added_layer_lin.ap_J̇₀⁻[rows_I, 1, :, p] .=
                -μ₀ .* reshape(FT(2) .* albedo, 1, :) .* incident_I .*
                reshape(τ̇_sum[:, p], 1, :) ./ μ₀
        end
        for (k, p) in enumerate(iparams)
            @views added_layer_lin.ap_J̇₀⁻[rows_I, 1, :, p] .=
                μ₀ .* reshape(FT(2) .* basis[:, k], 1, :) .* incident_I
        end
    end

    for (k, p) in enumerate(iparams)
        @views added_layer_lin.ap_ṙ⁻⁺[:, :, :, p] .=
            reshape(R_unit_weighted, size(R_unit_weighted, 1),
                    size(R_unit_weighted, 2), 1) .*
            reshape(FT(2) .* basis[:, k], 1, 1, nSpec)
    end
    return synchronize_if_gpu()
end
