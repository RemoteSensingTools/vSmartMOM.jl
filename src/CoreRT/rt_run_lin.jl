#=

This file contains the entry point for running the linearized RT simulation, `rt_run`.

The linearized RT computes both the forward radiance and its analytic Jacobians with
respect to physical parameters (aerosol properties, gas absorption, surface albedo)
using the Matrix Operator Method (MOM) of Plass, Hansen & Kattawar (1973) with the
linearization approach of Sanghavi & Stephens (2013), Sanghavi, Davis & Eldering (2014).

**References:**
- Sanghavi, S. & Stephens, G. (2013). "Adaptation of the delta-M and Successive
  Order of Scattering methods to the matrix operator method." *JQSRT*, 117, 1–12.
- Sanghavi, S., Davis, A. & Eldering, A. (2014). "vSmartMOM: A linearized discrete 
  ordinate radiative transfer implementation." *JQSRT*, 146, 182–207.
- de Haan, J.F., Bosma, P.B. & Hovenier, J.W. (1987). "The adding method for
  multiple scattering calculations of polarized light." *A&A*, 183, 371–391.

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users.

=#

"""
    rt_run(model::RTModel, lin_model, NAer, NGas, NSurf; i_band=1)

Perform linearized Radiative Transfer and return both radiances and their Jacobians.

Computes the reflected and transmitted Stokes vectors at the top and bottom of the
atmosphere, along with their analytic derivatives with respect to `Nparams = NAer×7 + NGas + NSurf`
physical parameters via the linearized Matrix Operator Method.

# Arguments
- `model::RTModel`: Forward model containing optical properties, geometry, etc.
- `lin_model`: Linearized model containing derivatives of optical properties.
- `NAer::Int`: Number of aerosol types.
- `NGas::Int`: Number of trace gas species with variable VMR.
- `NSurf::Int`: Number of surface parameters (typically 1 for Lambertian albedo).
- `i_band::Integer=1`: Spectral band index.

# Returns
- `R`: Reflected Stokes vector `[nVZA × nStokes × nSpec]`
- `T`: Transmitted Stokes vector `[nVZA × nStokes × nSpec]`
- `dR`: Jacobian of R, `[nVZA × nStokes × nSpec × Nparams]`
- `dT`: Jacobian of T, `[nVZA × nStokes × nSpec × Nparams]`

# Parameter Layout in `dR` / `dT`
The `Nparams` derivative dimension is ordered as:
1. **Aerosol sub-parameters** (7 per aerosol type):
   `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]` for each aerosol, so indices `1:7*NAer`
2. **Gas VMR parameters**: indices `7*NAer+1 : 7*NAer+NGas`
3. **Surface parameters**: indices `7*NAer+NGas+1 : Nparams`

# Theory
The forward model solves the vector radiative transfer equation via the discrete ordinate
method with the Matrix Operator Method (MOM). For each atmospheric layer ``k``, the
elemental reflection ``\\mathbf{r}`` and transmission ``\\mathbf{t}`` matrices are computed
from single-scattering, then doubled ``n_d`` times to obtain the full-layer matrices.
Layers are then combined via the adding (interaction) method from TOA to surface.

The linearized version simultaneously propagates derivatives with respect to three core
optical properties per layer:
```math
\\frac{\\partial \\mathbf{R}}{\\partial p_j} = 
  \\frac{\\partial \\mathbf{R}}{\\partial \\tau_k} \\frac{\\partial \\tau_k}{\\partial p_j} +
  \\frac{\\partial \\mathbf{R}}{\\partial \\varpi_k} \\frac{\\partial \\varpi_k}{\\partial p_j} +
  \\frac{\\partial \\mathbf{R}}{\\partial \\mathbf{Z}_k} \\frac{\\partial \\mathbf{Z}_k}{\\partial p_j}
```
where ``p_j`` is any physical parameter in the state vector.
"""

# Mockup if no Raman type is chosen:
function rt_run(model,
        lin_model,
        NAer::Int, NGas::Int, NSurf::Int;
        i_band::Integer = 1)
    rt_run(noRS(), model, lin_model, NAer, NGas, NSurf, i_band)
end

"""
    rt_run_lin(model, lin_model, NAer, NGas, NSurf; i_band=1)

Convenience alias for the linearized `rt_run` overload.  Equivalent to
`rt_run(model, lin_model, NAer, NGas, NSurf; i_band)`.
"""
rt_run_lin(model, lin_model,
           NAer::Int, NGas::Int, NSurf::Int; i_band::Integer = 1) =
    rt_run(model, lin_model, NAer, NGas, NSurf; i_band)

# Just to make sure we still have it:
function rt_run_test(RS_type::AbstractRamanType,
        model,
        lin_model,
        NAer, NGas, NSurf,
        iBand)
    rt_run(RS_type, model, lin_model,
        NAer, NGas, NSurf,
        iBand)
end

# Full multiple scattering
function rt_run(RS_type::AbstractRamanType,
                    model,
                    lin_model,
                    NAer::Int, NGas::Int, NSurf::Int,
                    iBand)
    (; obs_alt, sza, vza, vaz) = model.obs_geom   # Observational geometry properties
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad) = model.quad_points # All quadrature points
    pol_type = CoreRT.polarization_type(model)
    #@unpack max_m = model.max_m #params
    quad_points = model.quad_points
    max_m = model.max_m
    (; τ̇_abs, τ̇_aer, lin_aerosol_optics) = lin_model

    lin = LinMode()
    #@show iBand, sza, vza, vaz, model.params.brdf[iBand].albedo
    # Also to be changed!!
    #brdf = model.params.brdf[iBand[1]]
    #@show size(iBand)
    #@show iBand
    #@show iBand[1]
    #@show size(iBand[1])
    #bla
    brdf = get_surface(model, iBand)
    #brdf_lin = model_lin.brdf_lin[iBand]
    (; F₀) = RS_type
    # no Raman
    #if (typeof(RS_type)<:Union{RRS,RRS_plus})
    #    RS_type.ϖ_λ₁λ₀ .*=  (1. - model.ϖ_Cabannes[iBand])/sum(RS_type.ϖ_λ₁λ₀) # RS_type.ϖ_λ₁λ₀ .*=  (1. - model.ϖ_Cabannes[iBand[1]])/sum(RS_type.ϖ_λ₁λ₀) 
    #end   
    
    FT = eltype(sza)                    # Get the float-type to use

    Nz = length(model.profile.p_full)   # Number of vertical slices

    # For noRS: F₀ is the solar irradiance Stokes vector per spectral point.
    # Initialize to ones (unit solar flux, Stokes I only) if still at default 1×1 size.
    RS_type.bandSpecLim = UnitRange{Int}[]
    #put this code in model_from_parameters
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(model.τ_abs[iB], 1); # Number of spectral points
        push!(RS_type.bandSpecLim,nSpec0:nSpec);             
    end

    arr_type = CoreRT.array_type(model) # Type of array to use
    SFI = true                          # SFI flag
    NquadN = Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims = (NquadN,NquadN)              # nxn dims
    layout = ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                              n_gases=NGas, n_surface=NSurf)
    Nparams = n_total(layout)

    # For noRS: F₀ should be the solar irradiance Stokes vector per spectral point.
    # If still at its default 1×1 size, initialize to unit solar flux (Stokes I only).
    if size(F₀) != (pol_type.n, nSpec)
        F₀ = zeros(FT, pol_type.n, nSpec)
        F₀[1,:] .= one(FT)  # Only Stokes I = 1
        RS_type.F₀ = F₀
    end

    #FT = length(model.τ_aer[1][1]) > 0 ? typeof(model.τ_aer[1][1]) : FT
    # RT kernels always use pure FT (Float32/Float64), never Dual types
    # Output variables: Reflected and transmitted solar irradiation at TOA and BOA respectively # Might need Dual later!!
    #Suniti: consider adding a new dimension (iBand) to these arrays. The assignment of simulated spectra to their specific bands will take place after batch operations, thereby leaving the computational time unaffected 
    R       = zeros(FT, length(vza), pol_type.n, nSpec)
    T       = zeros(FT, length(vza), pol_type.n, nSpec)
    #R_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    #T_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    Ṙ       = zeros(FT, length(vza), pol_type.n, nSpec, Nparams)
    Ṫ       = zeros(FT, length(vza), pol_type.n, nSpec, Nparams)
    # Notify user of processing parameters
    msg = 
    """
    Processing on: $(architecture)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    # Create arrays
    @timeit "Creating layers" added_layer, added_layer_lin          = 
        make_added_layer(lin, RS_type, FT, arr_type, Nparams, dims, nSpec)
    # Just for now, only use noRS here

    @timeit "Creating layers" added_surface_layer, added_surface_layer_lin = 
        make_added_layer(lin, RS_type, FT, arr_type, Nparams, dims, nSpec)
    @timeit "Creating layers" composite_layer, composite_layer_lin  = 
        make_composite_layer(lin, RS_type, FT, arr_type, Nparams, dims, nSpec)
    @timeit "Creating arrays" I_static = 
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    #TODO: if RS_type!=noRS, create ϖ_λ₁λ₀, i_λ₁λ₀, fscattRayl, Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ (for input), and ieJ₀⁺, ieJ₀⁻, ieR⁺⁻, ieR⁻⁺, ieT⁻⁻, ieT⁺⁺, ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺ (for output)
    #getRamanSSProp(RS_type, λ, grid_in)

    # Loop over fourier moments
    for m = 0:max_m[iBand] - 1

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        #InelasticScattering.computeRamanZλ!(RS_type, pol_type,Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, layer_opt_props_lin, fScattRayleigh   = 
            constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model);
        #@show size(fScattRayleigh)
        #@show size(fScattRayleigh[1])
            # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all, τ̇_sum_all = 
            extractEffectiveProps(layer_opt_props, layer_opt_props_lin);
        
        
        # Loop over vertical layers: 
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA
            
            # Construct the atmospheric layer
            # From Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ
            # Suniti: modified to return fscattRayl as the last element of  computed_atmosphere_properties
            #if !(typeof(RS_type) <: noRS)
            #    RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz]) 
            #end
            
            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt, layer_opt_lin = 
                expandOpticalProperties(layer_opt_props[iz], layer_opt_props_lin[iz], arr_type)
            #@show size(layer_opt.Z⁺⁺[:,:,1]), size(RS_type.Z⁺⁺_λ₁λ₀)
            #@show typeof(layer_opt.Z⁺⁺[:,:,1]), typeof(RS_type.Z⁺⁺_λ₁λ₀)
            #aa = Array(layer_opt.Z⁺⁺[:,:,1]) #Array(RS_type.ϖ_Cabannes[1]*layer_opt.Z⁺⁺[:,:,1]) .+ (sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁺⁺_λ₁λ₀)
            #bb = Array(layer_opt.Z⁻⁺[:,:,1]) #Array(RS_type.ϖ_Cabannes[1]*layer_opt.Z⁻⁺[:,:,1]) .+ (sum(RS_type.ϖ_λ₁λ₀)*RS_type.Z⁻⁺_λ₁λ₀)
            #@show iz,(layer_opt), (layer_opt_lin)
            #@show iz, size(τ_sum_all), size(τ̇_sum_all)
            # Perform Core RT (doubling/elemental/interaction)
            rt_kernel!(RS_type::noRS, pol_type, SFI, 
                        #bandSpecLim, 
                        added_layer, added_layer_lin, 
                        composite_layer, composite_layer_lin,
                        layer_opt, layer_opt_lin,
                        scattering_interfaces_all[iz], 
                        τ_sum_all[:,iz], τ̇_sum_all[:,:,iz], 
                        m, quad_points, 
                        I_static, 
                        CoreRT.architecture(model), 
                        qp_μN, iz) 
        end 

        # Create surface matrices:
        iparam = surface_index(layout, iBand)
        create_surface_layer!(RS_type, brdf, #brdf_lin,
                            added_surface_layer, 
                            added_surface_layer_lin,
                            iparam,
                            SFI, m, 
                            pol_type, 
                            quad_points, 
                            arr_type(τ_sum_all[:,end]),
                            arr_type(τ̇_sum_all[:,:,end]), 
                            arr_type(F₀),
                            CoreRT.architecture(model));
        
        
        #@show F₀[:,1]
        #@show scattering_interfaces_all[end]
                            #@show scattering_interfaces_all[end]
        #blapl
        # One last interaction with surface:
        #@show composite_layer.J₀⁻[:,1,1] 
        
        @timeit "interaction" interaction!(#RS_type,
                                    #bandSpecLim,
                                    scattering_interfaces_all[end], 
                                    SFI, 
                                    #computed_layer_properties, computed_layer_properties_lin, 
                                    composite_layer, composite_layer_lin,
                                    added_surface_layer, added_surface_layer_lin,
                                    I_static)

        # Postprocess and weight according to vza
        postprocessing_vza!(RS_type, 
                            iμ₀, pol_type, 
                            composite_layer, 
                            composite_layer_lin,
                            vza, qp_μ, m, vaz, μ₀, 
                            weight, nSpec, 
                            SFI, 
                            R, 
                            T,
                            Ṙ, Ṫ)
        #@show R_SFI[:,1,1]
        #bla
    end
    
    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return R, T, Ṙ, Ṫ
    #return Array(added_layer.ieJ₀⁻), Array(composite_layer.ieJ₀⁻)#
end