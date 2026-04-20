#=

This file contains the entry point for running the RT simulation, rt_run. 

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users. 

=#




"""
    rt_run(model::RTModel; i_band=1) -> (R, T, ...)

Run the forward radiative transfer solver for one or more spectral bands.

Performs polarized adding-doubling RT through the atmosphere defined in `model`,
computing top-of-atmosphere (TOA) reflectance and bottom-of-atmosphere (BOA)
transmittance.  The solver iterates over azimuthal Fourier moments
`m = 0, …, max_m-1`, building layer R/T/J matrices via elemental → doubling →
interaction steps, then applies surface coupling and postprocessing.

Equivalent to `rt_run(noRS(), model, i_band)` (no Raman scattering).

# Arguments
- `model::RTModel`: Pre-built model from [`model_from_parameters`](@ref).
- `i_band::Integer=1`: Spectral band index (or vector of indices) to compute.

# Returns
A tuple `(R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)` when using
source-function integration (default), where:
- `R_SFI::Array{FT,3}`: TOA reflectance `[nVZA × nStokes × nSpec]`
- `T_SFI::Array{FT,3}`: BOA transmittance `[nVZA × nStokes × nSpec]`

For most use cases, only the first two elements are needed:
```julia
R, T = rt_run(model)
```

# Example
```julia
params = parameters_from_yaml("config/my_scene.yaml")
model = model_from_parameters(params)
R, T = rt_run(model)
R[1, 1, :]  # Stokes-I reflectance at first VZA across the spectrum
```

# See also
- [`rt_run(model, lin_model, NAer, NGas, NSurf)`](@ref) for linearized RT with Jacobians.
- [`model_from_parameters`](@ref) to build the model from parameters.
"""
function rt_run(model; i_band::Integer = 1)
    rt_run(InelasticScattering.noRS{float_type(model)}(), model, i_band)
end

"""
    rt_run_test(RS_type, model, iBand)

Test entry point for RT calculations with explicit Raman type.
"""
function rt_run_test(RS_type::AbstractRamanType, model, iBand)
    rt_run(RS_type, model, iBand)
end

"""
    rt_run(RS_type, model::RTModel, iBand)

Perform Radiative Transfer calculations with explicit Raman type.
"""
function rt_run(RS_type::AbstractRamanType, model, iBand)
    (; obs_alt, sza, vza, vaz) = model.obs_geom   # Observational geometry properties
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad) = model.quad_points # All quadrature points
    pol_type = CoreRT.polarization_type(model)
    max_m    = get_max_m(model)
    (; quad_points) = model
    FT       = CoreRT.float_type(model)

    n_aer = CoreRT.n_aerosols(model)

    # Also to be changed if more than 1 band is used!!
    brdf = get_surface(model, iBand[1])
    if length(iBand) > 1
        @info "More than one band has been chosen, be aware that multiple BRDFs are not yet implemented and only the first one will be used!"
    end

    (; ϖ_Cabannes) = RS_type

    Nz = length(model.profile.p_full)   # Number of vertical slices

    RS_type.bandSpecLim = UnitRange{Int}[]
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(model.τ_abs[iB], 1); # Number of spectral points
        push!(RS_type.bandSpecLim,nSpec0:nSpec);
    end

    arr_type = CoreRT.array_type(model) # Type of array to use
    arch     = CoreRT.architecture(model)
    SFI = true                          # SFI flag
    NquadN =  Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims   = (NquadN,NquadN)              # nxn dims

    # Output arrays for reflected and transmitted solar irradiation at TOA and BOA
    @timeit "Arrays"  R       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  R_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieR_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieT_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  hdr     = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  bhr_dw     = zeros(FT, pol_type.n, nSpec)
    @timeit "Arrays"  bhr_uw     = zeros(FT, pol_type.n, nSpec)
    @timeit "Arrays"  hdr_J₀⁻    = zeros(FT, length(vza), pol_type.n, nSpec)
    # Notify user of processing parameters
    msg =
    """
    Processing on: $(arch)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    # Create arrays
    @timeit "Creating layers" added_layer         =
        make_added_layer(RS_type, FT, arr_type, dims, nSpec)
    # Just for now, only use noRS here
    @timeit "Creating layers" added_layer_surface =
        make_added_layer(RS_type, FT, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     =
        make_composite_layer(RS_type, FT, arr_type, dims, nSpec)
    @timeit "Creating arrays" I_static =
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));

        # Note: Raman SS properties (ϖ_λ₁λ₀, Z matrices, etc.) are set up
    # in model_from_parameters when RS_type != noRS.

    # Build concatenated wavenumber grid for canopy spectral features
    _canopy_spec_wn = nothing
    if brdf isa CanopySurface
        _canopy_spec_wn = vcat([get_spec_bands(model)[iB] for iB in iBand]...)
    end

    # Pre-initialize canopy cache before the Fourier loop (Zazi precomputation is expensive)
    if brdf isa CanopySurface && brdf._cache === nothing
        @timeit "Canopy cache init" _init_canopy_cache!(
            brdf, added_layer_surface, pol_type, quad_points, arch;
            spec_bands_wn=_canopy_spec_wn, max_m=max_m)
    end

    # Pre-compute within-canopy atmospheric optical depth if requested
    if brdf isa CanopySurface && brdf.include_atm && brdf.canopy_dp !== nothing
        @timeit "Canopy atm tau" _compute_canopy_atm_tau!(brdf, model, _canopy_spec_wn)
    end

    # Initialize F₀ (solar irradiance Stokes vector) if still at default size.
    # Default: unit Stokes I across all spectral points; users can set actual
    # solar spectral irradiance for broadband instruments.
    if size(RS_type.F₀) != (pol_type.n, nSpec)
        F₀ = zeros(FT, pol_type.n, nSpec)
        F₀[1,:] .= one(FT)  # Only Stokes I = 1
        RS_type.F₀ = F₀
    end

    # Cumulative optical depth (m-independent, saved for TMS correction)
    τ_sum_all = nothing

    # Loop over fourier moments
    for m = 0:max_m - 1

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5) : FT(1.0)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        @timeit "IE"  InelasticScattering.computeRamanZλ!(RS_type, pol_type,collect(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, fScattRayleigh   =
            constructCoreOpticalProperties(RS_type,iBand,m,model);
        # Determine the scattering interface definitions:
        @timeit "Extract Optical Properties" scattering_interfaces_all, τ_sum_all =
            extractEffectiveProps(layer_opt_props,quad_points);

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            if !(typeof(RS_type) <: noRS)
                @timeit "Expand Bands" RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz])
            end

            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt =
                expandOpticalProperties(layer_opt_props[iz], arr_type)

            # Perform Core RT (doubling/elemental/interaction)
            @timeit "RT Kernel" rt_kernel!(RS_type, pol_type, SFI,
                        added_layer, composite_layer,
                        layer_opt,
                        scattering_interfaces_all[iz],
                        τ_sum_all[:,iz],
                        m, quad_points,
                        I_static,
                        arch,
                        qp_μN, iz)
        end

        # Create surface matrices:
        if brdf isa CanopySurface
            @timeit "Create Surface" create_surface_layer!(brdf,
                                added_layer_surface,
                                SFI, m,
                                pol_type,
                                quad_points,
                                arr_type(τ_sum_all[:,end]),
                                arch;
                                spec_bands_wn=_canopy_spec_wn,
                                max_m=max_m)
        else
            @timeit "Create Surface" create_surface_layer!(brdf,
                                added_layer_surface,
                                SFI, m,
                                pol_type,
                                quad_points,
                                arr_type(τ_sum_all[:,end]),
                                arch)
        end
        
        #@show composite_layer.J₀⁺[iμ₀,1,1:3]
        # One last interaction with surface:
        @timeit "interaction" interaction!(RS_type,
                                    #bandSpecLim,
                                    scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    I_static)
       #@show composite_layer.J₀⁺[iμ₀,1,1:3]
        hdr_J₀⁻ = similar(composite_layer.J₀⁻)
        # One last interaction with surface:
        @timeit "interaction_HDRF" interaction_hdrf!(#RS_type,
                                    #bandSpecLim,
                                    #scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    m, pol_type, quad_points,
                                    hdr_J₀⁻, bhr_uw, bhr_dw)
        
        # Postprocess and weight according to vza
        @timeit "Postprocessing VZA" postprocessing_vza!(RS_type, 
                            iμ₀, pol_type, 
                            composite_layer, 
                            vza, qp_μ, m, vaz, μ₀, 
                            weight, nSpec, 
                            SFI, 
                            R, R_SFI, 
                            T, T_SFI,
                            ieR_SFI, ieT_SFI)

        @timeit "Postprocessing HDRF" postprocessing_vza_hdrf!(RS_type, 
            iμ₀, pol_type, 
            hdr_J₀⁻, 
            vza, qp_μ, m, vaz, μ₀, 
            weight, nSpec, 
            hdr)
            
    end

    # Single-scattering correction for Cox-Munk specular hotspot (TMS)
    if brdf isa CoxMunkSurface && SFI
        @timeit "SS Correction" apply_ss_correction!(
            R_SFI, brdf, pol_type, vza, vaz, μ₀,
            Array(τ_sum_all[:,end]), max_m, nSpec)
    end

    # Show timing statistics
    print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    #if RAMI
    #@show size(hdr), size(bhr_dw)
    #hdr = hdr[:,1,:] ./ bhr_dw[1,:]
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw[1,:], bhr_dw[1,:]) : (R, T)
    #else
    #return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI) : (R, T)
    #end
end