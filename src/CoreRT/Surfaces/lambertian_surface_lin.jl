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
- `τ̇_sum`: Derivative of total τ `[Nparams × nSpec]`.
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
    
    nparams = size(τ̇_sum,1) # nparams ≠ Nparams (in general)
    nspec = length(τ_sum)
#@show nparams, nspec
#@show iμ₀Nstart, iμ₀
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    tmp = arr_type(ones(pol_type.n*Nquad))
    T_surf = Diagonal(tmp)
    i₀ = iμ₀Nstart:iμ₀Nstart+n-1
    #@show i₀
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
            F₀_NquadN = arr_type(zeros(length(qp_μN),length(τ_sum)));
            Ḟ₀_NquadN = arr_type(zeros(nparams+1,length(qp_μN),length(τ_sum)));
            #F₀_NquadN[:] .=0;
            #@show size(F₀), size(μ₀), size(R_surf), size(F₀_NquadN), size(added_layer.j₀⁻[:,1,:]), size(F₀ .* (exp.(-τ_sum/μ₀))'), size(F₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀,:])
            #@show iμ₀Nstart, pol_type.n*iμ₀
            #@show size(F₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀,:]), size(F₀), size(exp.(-τ_sum/μ₀)')
            tmpF = (F₀ .* arr_type(exp.(-τ_sum/μ₀))');
            F₀_NquadN[i₀,:] .= tmpF 
            #@show size(Ḟ₀_NquadN[:,i₀,:]), size(-reshape(tmpF,1,n,nspec).*reshape(τ̇_sum,nparams, 1, nspec)/μ₀)
            Ḟ₀_NquadN[1:nparams,i₀,:] .= -reshape(tmpF,1,n,nspec).*reshape(τ̇_sum, nparams, 1, nspec)/μ₀ # , arr_type(zeros(1, n, nspec)); dims=1)

            added_layer.j₀⁺[:,:,:] .= zero(FT);#
            added_layer.j₀⁻[:,1,:] .= μ₀*(R_surf*F₀_NquadN)#/FT(π);
            #added_layer_lin.J̇₀⁺[:,:,:,:] .= 0.;#
            

            added_layer_lin.ap_J̇₀⁺[:,:,1,:] .= zero(FT)
            for ii=1:nspec
                for ctr=1:nparams
                    added_layer_lin.ap_J̇₀⁻[ctr,:,1,ii] .= μ₀*R_surf*Ḟ₀_NquadN[ctr,:,ii]#/FT(π);
                end
                added_layer_lin.ap_J̇₀⁻[iparam,:,1,ii]  .= μ₀*Ṙ_surf[:,:]*F₀_NquadN[:,ii]#/FT(π);
            end    
                #@show size(added_layer.j₀⁻[:,1,1])
        #@show added_layer.j₀⁻[:,1,1]
        # for SIF
            #reinstate the following line after linearization works
            #added_layer.j₀⁻[:,1,:] .+= (1/π)*repeat(arr_type(RS_type.SIF₀),Nquad) * unweight
        end

        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        Ṙ_surf[:,:] = Ṙ_surf[:,:] * Diagonal(qp_μN.*wt_μN)
        
        #R_surf = 2R_surf * Diagonal(qp_μN.*wt_μN)
        #R_surf = R_surf * Diagonal(qp_μN.*wt_μN)/π

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.j₀⁻)
        #@show size(added_layer.r⁻⁺), size(R_surf)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;#1. #0.0; #T_surf;
        added_layer.t⁻⁻ .= zero(FT); #T_surf;

        added_layer_lin.ap_ṙ⁻⁺[iparam,:,:,:] .= Ṙ_surf;
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
    #@show size(T_surf), size(R_surf)
    #@show T_surf
    #@show R_surf
    synchronize_if_gpu()
end