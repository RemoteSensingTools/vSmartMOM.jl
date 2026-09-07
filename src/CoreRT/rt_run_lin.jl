#=

This file contains the entry point for running the linearized RT simulation, `rt_run`.

The linearized RT computes both the forward radiance and its analytic Jacobians with
respect to physical parameters (aerosol properties, gas absorption, surface albedo)
using the Matrix Operator Method (MOM) of Plass, Hansen & Kattawar (1973) with the
linearization approach of Sanghavi, Davis & Eldering (2014), following the
scalar formulation of Sanghavi, Martonchik, Davis & Diner (2013).

**References:**
- Sanghavi, S., Davis, A.B. & Eldering, A. (2014). "vSmartMOM: A vector matrix
  operator method-based radiative transfer model linearized with respect to
  aerosol properties." *JQSRT*, 133, 412–433; Eqs. (23)–(28), App. C.
- Sanghavi, S. & Frankenberg, C. (2023). "Raman scattering in the Earth's
  atmosphere, Part II: Radiative transfer modeling for remote sensing
  applications." *JQSRT*, 311, 108791; elastic Eqs. (10)–(12).
- de Haan, J.F., Bosma, P.B. & Hovenier, J.W. (1987). "The adding method for
  multiple scattering calculations of polarized light." *A&A*, 183, 371–391.

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users.

=#

"""
    rt_run(model::RTModel, lin_model, NAer, NGas, NSurf; i_band=1)

Perform linearized Radiative Transfer and return both radiances and their Jacobians.

Computes the reflected and transmitted Stokes vectors at the top and bottom of the
atmosphere, along with analytic derivatives with respect to
`Nparams = 1 + NAer×7 + NGasLayer + NSurf + NSIF` parameters, where
`NGasLayer = NGasSpecies×Nz`.
physical parameters via the linearized Matrix Operator Method.

# Arguments
- `model::RTModel`: Forward model containing optical properties, geometry, etc.
- `lin_model`: Linearized model containing derivatives of optical properties.
- `NAer::Int`: Number of aerosol types.
- `NGas::Int`: Number of layer-resolved gas VMR parameters (`NGasSpecies×Nz`).
- `NSurf::Int`: Number of surface parameters (typically 1 for Lambertian albedo).
- `i_band::Integer=1`: Spectral band index.
- `jacobian_adding=:matrix`: Use established matrix-tangent adding. Opt in to
  `:source` with a local basis for solar/Lambertian endpoint solves, including
  prescribed or retrievable surface SIF; atmospheric
  adding then propagates equivalent-source vectors instead of matrix tangents.
- `jacobian_basis=:auto`: Double a local optical basis when smaller than the
  atmospheric retrieval layout. `:physical` propagates retrieval columns
  directly; `:local` forces the local basis for a single-band solve.

# Returns
An [`ObserverRTResultLin`](@ref). It remains iterable as the historical
`(R, T, dR, dT)` endpoint tuple, while `result.levels` exposes forward
radiances and Jacobians at requested strict-interior observer heights.

# Parameter Layout in `dR` / `dT`
The `Nparams` derivative dimension is ordered as:
1. **Surface pressure** `p_surf` (hPa): column 1.
2. **Aerosol sub-parameters** (7 per aerosol type):
   `[τ_ref, nᵣ, nᵢ, μ_logr, σ_logr, profile_location, profile_width]` for each
   aerosol. The profile pair is `(p₀, σp)` for
   `Normal` and `(z₀, σ₀)` for `LogNormal`.
3. **Gas VMR parameters**, species-major with all `Nz` layers for each gas.
4. **Surface parameters**.
5. **SIF parameters**, when `SurfaceSIF(SIF755=...)` is present:
   `[SIF755, slope]`, referenced to 755 nm.

Use `result.layout` with `psurf_index`, `aerosol_range`, `gas_profile_range`,
`gas_layer_index`, `surface_range`, and `sif_range`; do not hard-code offsets.

# Theory
The forward model solves the vector radiative transfer equation via the discrete ordinate
method with the Matrix Operator Method (MOM). For each atmospheric layer ``k``, the
elemental reflection ``\\mathbf{r}`` and transmission ``\\mathbf{t}`` matrices are computed
from single-scattering, then doubled ``n_d`` times to obtain the full-layer matrices.
Layers are then combined via the adding (interaction) method from TOA to surface.

At the elemental boundary the code contracts core-optics partials with the
supplied tangent directions, then propagates those directional derivatives
through doubling. A local optical basis is contracted to retrieval columns
before adding, as described in `build_doubled_layer_lin!`. Schematically, for all layers:
```math
\\frac{\\partial \\mathbf{R}}{\\partial p_j} = 
  \\sum_k \\left[\\frac{\\partial \\mathbf{R}}{\\partial \\tau_k} \\frac{\\partial \\tau_k}{\\partial p_j} +
  \\frac{\\partial \\mathbf{R}}{\\partial \\varpi_k} \\frac{\\partial \\varpi_k}{\\partial p_j} +
  \\frac{\\partial \\mathbf{R}}{\\partial \\mathbf{Z}_k} \\frac{\\partial \\mathbf{Z}_k}{\\partial p_j}\\right]
```
where ``p_j`` is any physical parameter in the state vector. The Z term
denotes contraction over both forward/backward phase-matrix indices; after
doubling it cannot be represented by an elementwise scalar partial.
"""

# Mockup if no Raman type is chosen:
function rt_run(model,
        lin_model,
        NAer::Int, NGas::Int, NSurf::Int;
        i_band::Integer = 1,
        sources::Union{Nothing, AbstractSource} = nothing,
        jacobian_basis::Symbol=:auto,
        jacobian_adding::Symbol=:matrix)
    rt_run(InelasticScattering.noRS{float_type(model)}(), model, lin_model, NAer, NGas, NSurf, i_band; sources, jacobian_basis, jacobian_adding)
end

"""
    rt_run_lin(model, lin_model, NAer, NGas, NSurf; i_band=1)

Convenience alias for the linearized `rt_run` overload.  Equivalent to
`rt_run(model, lin_model, NAer, NGas, NSurf; i_band)`.
"""
rt_run_lin(model, lin_model,
           NAer::Int, NGas::Int, NSurf::Int;
           i_band::Integer = 1,
           sources::Union{Nothing,AbstractSource} = nothing,
           jacobian_basis::Symbol=:auto,
           jacobian_adding::Symbol=:matrix) =
    rt_run(model, lin_model, NAer, NGas, NSurf; i_band, sources, jacobian_basis, jacobian_adding)

"""
    rt_run(model, lin_model::PlannedRTModelLin; i_band=1, sources=nothing)

Run a retrieval-selected linearized calculation. The selected band layout is
compiled by the retrieval flavour and its compact tangent dimension is used
throughout the MOM kernels. Use [`globalize_jacobian`](@ref) to scatter the
returned band-local Jacobian into the shared retrieval state.
"""
function rt_run(model, lin_model::PlannedRTModelLin;
                i_band::Integer=1,
                sources::Union{Nothing,AbstractSource}=nothing,
                jacobian_basis::Symbol=:auto,
                jacobian_adding::Symbol=:matrix)
    layout = band_layout(lin_model.plan, i_band)
    NAer = CoreRT.n_aerosols(model)
    NGas = size(lin_model.base.τ̇_abs[i_band], 1)
    NSurf = surface_parameter_count(get_surface(model, i_band))
    return rt_run(InelasticScattering.noRS{float_type(model)}(),
                  model, lin_model.base, NAer, NGas, NSurf, i_band;
                  sources, active_layout=layout, jacobian_basis, jacobian_adding)
end

rt_run_lin(model, lin_model::PlannedRTModelLin;
           i_band::Integer=1,
           sources::Union{Nothing,AbstractSource}=nothing,
           jacobian_basis::Symbol=:auto,
           jacobian_adding::Symbol=:matrix) =
    rt_run(model, lin_model; i_band, sources, jacobian_basis, jacobian_adding)

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
                    iBand;
                    sources::Union{Nothing, AbstractSource} = nothing,
                    active_layout::Union{Nothing,ActiveParameterLayout} = nothing,
                    jacobian_basis::Symbol = :auto,
                    jacobian_adding::Symbol = :matrix)
    if InelasticScattering.has_inelastic(RS_type)
        throw(ArgumentError(
            "Linearized Raman-active RT is intentionally unsupported. " *
            "Use forward Raman RT only; do not call rt_run with lin_model for RRS/VS modes."))
    end

    # Per-model BLAS thread cap (see `rt_run` body for rationale).
    if model.numerics.blas_threads !== nothing
        LinearAlgebra.BLAS.set_num_threads(model.numerics.blas_threads)
    end

    (; obs_alt, sza, vza, vaz) = model.obs_geom   # Observational geometry properties
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad) = model.quad_points # All quadrature points
    pol_type = CoreRT.polarization_type(model)
    quad_points = model.quad_points
    # Numerics knobs threaded into the lin rt_kernel! (matches rt_run forward path).
    dτ_max_threshold = model.numerics.dτ_max_threshold
    dτ_min_floor     = model.numerics.dτ_min_floor
    # Per-band Fourier loop bound (order). Phase B unifies forward and
    # lin paths through `m_max_bands(model)`, fixing the lin-only
    # precedence bug at the previous lin formula.
    m_max = m_max_bands(model)[iBand]
    (; τ̇_abs, τ̇_aer, lin_aerosol_optics) = lin_model

    lin = LinMode()
    brdf = get_surface(model, iBand)
    expected_nsurf = surface_parameter_count(brdf)
    NSurf == expected_nsurf || throw(ArgumentError(
        "NSurf=$NSurf does not match $(typeof(brdf)), which has " *
        "$expected_nsurf linearized surface parameter(s)."))
    (; F₀) = RS_type

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
    # Resolve sources (v0.6 source-term refactor). Resolution: kwarg >
    # model.sources > pre-set `RS_type.F₀` for back-compat. See `rt_run`.
    # `prepared_sources` is now hoisted outside the conditional so it's
    # in scope for the surface step (`surface_source_contribute!`) that
    # routes SurfaceSIF / future per-source surface contributions.
    effective_sources = sources === nothing ? model.sources : sources
    validate_sif_solar_spectrum(effective_sources)
    NSIF = surface_sif_parameter_count(effective_sources)
    layout = active_layout === nothing ?
        ParameterLayout(aerosol_params=7, n_aerosols=NAer,
                        n_gases=NGas, n_surface=NSurf, n_sif=NSIF) :
        active_layout
    if active_layout !== nothing
        length(surface_range(layout)) == expected_nsurf || throw(ArgumentError(
            "active layout has $(length(surface_range(layout))) surface columns; " *
            "$expected_nsurf are required by $(typeof(brdf))"))
        length(sif_range(layout)) == NSIF || throw(ArgumentError(
            "active layout has $(length(sif_range(layout))) SIF columns but " *
            "the selected sources expose $NSIF"))
    end
    Nparams = n_total(layout)
    if !isempty(model.obs_geom.sensor_levels)
        _multisensor_source_supported(effective_sources) || throw(ArgumentError(
            "linearized interior-height radiances currently support " *
            "SolarBeam/NoSource only; thermal and surface-emission sources " *
            "require linearized multisensor source-slot propagation"))
        _require_multisensor_surfaces(model, (iBand,))
    end
    # External-solar SFI is TOA-only (see the forward driver's guard): the
    # interior/BOA direct carrier selects views via the embedded-solar iμ₀
    # node, which is a deliberate sentinel (0) in external mode. Interior
    # heights therefore fail loudly here, and the BOA endpoint is reported
    # as unavailable (`nothing`) at result assembly below instead of
    # returning a radiance that silently lacks the direct field.
    model.quad_points.external_solar && !isempty(model.obs_geom.sensor_levels) &&
        throw(ArgumentError(
            "external-solar SFI does not support interior sensor heights; " *
            "build the model with external_solar=false (the default)"))
    # Forward/lin parity policy: the TMS correction is not implemented in
    # the linearized driver yet, so it is rejected — never silently ignored
    # (the forward and linearized radiances of one configured model must
    # not diverge). See docs/dev_notes/forward_lin_parity_policy.md.
    model.numerics.ss_correction isa TMSCorrection && throw(ArgumentError(
        "TMSCorrection is not implemented in the linearized driver; " *
        "use NoSSCorrection for LinMode runs (parity policy: features are " *
        "matched or rejected, never silently divergent)"))
    prepared_sources = prepare_sources(effective_sources, FT, pol_type.n, nSpec, arr_type)
    if sources === nothing && size(F₀) == (pol_type.n, nSpec) && !iszero(F₀)
        # User pre-set F₀ honored — F₀ already in scope from RS_type unpack.
    else
        F₀_dev = extract_solar_F₀(prepared_sources, FT, pol_type.n, nSpec, arr_type)
        F₀ = Array{FT, 2}(F₀_dev)
        RS_type.F₀ = F₀
    end

    R       = zeros(FT, length(vza), pol_type.n, nSpec)
    T       = zeros(FT, length(vza), pol_type.n, nSpec)
    Ṙ       = zeros(FT, length(vza), pol_type.n, nSpec, Nparams)
    Ṫ       = zeros(FT, length(vza), pol_type.n, nSpec, Nparams)
    sensor_levels = model.obs_geom.sensor_levels
    level_uw = [zeros(FT, length(vza), pol_type.n, nSpec) for _ in sensor_levels]
    level_dw = [zeros(FT, length(vza), pol_type.n, nSpec) for _ in sensor_levels]
    level_direct_dw = [zeros(FT, length(vza), pol_type.n, nSpec) for _ in sensor_levels]
    level_uw_lin = [zeros(FT, length(vza), pol_type.n, nSpec, Nparams) for _ in sensor_levels]
    level_dw_lin = [zeros(FT, length(vza), pol_type.n, nSpec, Nparams) for _ in sensor_levels]
    level_direct_dw_lin = [zeros(FT, length(vza), pol_type.n, nSpec, Nparams) for _ in sensor_levels]
    # Notify user of processing parameters
    msg = 
    """
    Processing on: $(CoreRT.architecture(model))
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    local_basis = use_local_jacobian(jacobian_basis, layout, iBand, NAer)
    source_adding = use_source_adding(jacobian_adding,model,effective_sources,SFI,local_basis,brdf)
    local_workspace = local_basis && !source_adding ? make_local_jacobian_workspace(
        RS_type, FT, arr_type, local_jacobian_size(NAer), dims, nSpec,
        quad_points, pol_type) : nothing

    # Create arrays
    @timeit "Creating layers" added_layer, added_layer_lin          = 
        make_added_layer(lin, RS_type, FT, arr_type,
                         source_adding ? local_jacobian_size(NAer) : Nparams, dims, nSpec;
                         external_solar=quad_points.external_solar,
                         nStokes=pol_type.n, doubling_scratch=source_adding || !local_basis,
                         core_partials=source_adding || !local_basis)
    # Just for now, only use noRS here

    @timeit "Creating layers" added_surface_layer, added_surface_layer_lin = 
        make_added_layer(lin, RS_type, FT, arr_type,
                         source_adding ? NSurf : Nparams, dims, nSpec;
                         external_solar=quad_points.external_solar,
                         nStokes=pol_type.n, doubling_scratch=false,
                         core_partials=!source_adding)
    source_workspace = source_adding ? make_source_adding_workspace(
        added_layer,added_layer_lin,Nparams,NSurf,NSIF,Nz,arr_type) : nothing
    @timeit "Creating layers" composite_layer, composite_layer_lin = if source_adding
        (make_composite_layer(RS_type,FT,arr_type,dims,nSpec),source_workspace.result)
    else
        make_composite_layer(lin,RS_type,FT,arr_type,Nparams,dims,nSpec)
    end
    # Each interior observer owns two ordinary forward/tangent composites:
    # the atmosphere above it and the atmosphere below it. Reusing the
    # production CompositeLayer types means every per-layer Jacobian continues
    # through the same analytic interaction kernels as the endpoint column.
    top_pairs = [make_composite_layer(lin, RS_type, FT, arr_type,
                                      Nparams, dims, nSpec)
                 for _ in sensor_levels]
    bottom_pairs = [make_composite_layer(lin, RS_type, FT, arr_type,
                                         Nparams, dims, nSpec)
                    for _ in sensor_levels]
    top_layers = first.(top_pairs)
    top_layers_lin = last.(top_pairs)
    bottom_layers = first.(bottom_pairs)
    bottom_layers_lin = last.(bottom_pairs)
    bottom_interfaces = AbstractScatteringInterface[
        ScatteringInterface_00() for _ in sensor_levels]
    @timeit "Creating arrays" I_static = 
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));
    # Linearized RT intentionally supports pure-elastic noRS only.
    τ_sum_endpoint = nothing
    τ̇_sum_endpoint = nothing

    # τ, ϖ, their physical-parameter tangents, and gas objects do not depend
    # on Fourier order. Build them once; each moment attaches only Z(m), Ż(m),
    # Z₀(m), and Ż₀(m) before the analytic RT propagation.
    @timeit "OpticalProps invariant" m_invariant_cache =
        build_m_invariant_cache_lin(iBand, model, lin_model;
                                    active_layout)
    if local_basis
        @timeit "OpticalProps invariant" m_invariant_cache =
            build_local_jacobian_cache(iBand, model, lin_model, m_invariant_cache)
    end

    # The combined forward/analytic-linearization path uses the same forward
    # convergence decision as `rt_run` (I only or all Stokes components,
    # according to the selected strategy). Jacobian fields are accumulated
    # through the accepted moment and stop at that identical m; they do not
    # select a separate, potentially inconsistent Fourier range. Interior
    # sensor solves retain the complete series, matching the forward
    # multisensor contract.
    fourier_convergence = model.numerics.fourier_convergence
    convergence_active = _fourier_convergence_active(fourier_convergence) &&
                         SFI && isempty(sensor_levels)
    convergence_outputs = convergence_active ? _fourier_outputs(R, T) : ()
    convergence_snapshots = convergence_active ?
        _fourier_snapshots(convergence_outputs) : ()
    convergence_passes = 0
    m_used = m_max

    # Loop over fourier moments
    for m = 0:m_max

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        fill!(bottom_interfaces, ScatteringInterface_00())
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        #InelasticScattering.computeRamanZλ!(RS_type, pol_type,Array(qp_μ), m, arr_type)
        # Compute the core layer optical properties:
        @timeit "OpticalProps" layer_opt_props, layer_opt_props_lin, fScattRayleigh   = 
            (local_basis ? construct_local_optical_jacobians : constructCoreOpticalProperties)(
                RS_type, iBand, m, model, lin_model, m_invariant_cache);
        if active_layout !== nothing
            actual = size(layer_opt_props_lin[1].τ̇, 2)
            expected = n_layer_params(active_layout)
            actual == expected || throw(DimensionMismatch(
                "selective optical-property assembly returned $actual layer " *
                "columns; the active layout requires $expected"))
        end
            # Determine the scattering interface definitions:
        scattering_interfaces_all, τ_sum_all, τ̇_sum_all = 
            extractEffectiveProps(layer_opt_props, layer_opt_props_lin);

        # The collimated beam is outside the diffuse MOM source field. Report
        # it once at m=0, together with its atmospheric optical-depth tangent.
        if m == 0
            # The optical-depth sums are Fourier-independent. Retain the m=0
            # references explicitly because the loop introduces a local scope.
            τ_sum_endpoint = τ_sum_all
            τ̇_sum_endpoint = τ̇_sum_all
            if !isempty(sensor_levels)
                _set_unscattered_downwelling_lin!(
                    level_direct_dw, level_direct_dw_lin,
                    sensor_levels, τ_sum_all, τ̇_sum_all,
                    F₀, vza, qp_μ, iμ₀, μ₀, pol_type)
            end
        end
        
        
        # Loop over vertical layers: 
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA
            
            # Construct the atmospheric layer:
            # from Rayleigh and aerosol τ, ϖ, compute overall layer τ, ϖ.
            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt, layer_opt_lin = 
                expandOpticalProperties(layer_opt_props[iz], layer_opt_props_lin[iz], arr_type)
            # Perform Core RT (doubling/elemental/interaction)
            if source_adding
                source_adding_layer!(source_workspace,iz,composite_layer,
                    added_layer,added_layer_lin,RS_type,pol_type,layer_opt,layer_opt_lin,
                    τ_sum_all[:,iz],m,quad_points,I_static,CoreRT.architecture(model),
                    scattering_interfaces_all[iz];dτ_max_threshold,dτ_min_floor)
            else
                rt_kernel!(RS_type::noRS, pol_type, SFI,
                    added_layer, added_layer_lin, composite_layer, composite_layer_lin,
                    layer_opt, layer_opt_lin, scattering_interfaces_all[iz],
                    τ_sum_all[:,iz], τ̇_sum_all[:,:,iz], m, quad_points,
                    I_static, CoreRT.architecture(model), qp_μN, iz;
                    dτ_max_threshold, dτ_min_floor, local_workspace)
            end

            # Reuse the completed current added layer to grow the two
            # subcolumns belonging to every interior observer. The top
            # subcolumn follows the ordinary TOA-down interface dispatch.
            # The bottom subcolumn tracks its own scattering state because it
            # starts below the observer rather than at TOA.
            if !isempty(sensor_levels)
                layer_scatter = maximum(layer_opt.τ .* layer_opt.ϖ) > 2eps(FT)
                for ims in eachindex(sensor_levels)
                    boundary = sensor_levels[ims]
                    if iz == 1
                        seed_composite_from_added!(
                            top_layers[ims], top_layers_lin[ims],
                            added_layer, added_layer_lin)
                    elseif iz <= boundary
                        interaction!(scattering_interfaces_all[iz], SFI,
                                     top_layers[ims], top_layers_lin[ims],
                                     added_layer, added_layer_lin, I_static)
                    end

                    local_iz = iz - boundary
                    if local_iz == 1
                        seed_composite_from_added!(
                            bottom_layers[ims], bottom_layers_lin[ims],
                            added_layer, added_layer_lin)
                        bottom_interfaces[ims] = get_scattering_interface(
                            ScatteringInterface_00(), layer_scatter, 1)
                    elseif local_iz > 1
                        bottom_interfaces[ims] = get_scattering_interface(
                            bottom_interfaces[ims], layer_scatter, local_iz)
                        interaction!(bottom_interfaces[ims], SFI,
                                     bottom_layers[ims], bottom_layers_lin[ims],
                                     added_layer, added_layer_lin, I_static)
                    end
                end
            end
        end 

        # Create surface matrices. `surface_index(layout, i)` indexes WITHIN
        # the surface block (1..n_surface). For per-band rt_run calls the
        # layout has n_surface=1, so the surface-albedo Jacobian always lands
        # at the first surface slot, regardless of which atmospheric band this
        # rt_run call is processing. (Earlier code passed `iBand` here, which
        # blew through the surface block when iBand>1.)
        surface_columns = source_adding ? (1:NSurf) : surface_range(layout)
        iparam = brdf isa LambertianSurfaceLegendre ? surface_columns : first(surface_columns)
        create_surface_layer!(RS_type, brdf, #brdf_lin,
                            added_surface_layer,
                            added_surface_layer_lin,
                            iparam,
                            SFI, m,
                            pol_type,
                            quad_points,
                            arr_type(τ_sum_all[:,end]),
                            source_adding ? source_workspace.surface_zero_above :
                                arr_type(τ̇_sum_all[:,:,end]),
                            arr_type(F₀),
                            CoreRT.architecture(model));

        # Surface source contributions (v0.6 source-term refactor): mirror
        # the rt_run / rt_run_ss surface step. NoSource → no-op;
        # PreparedSurfaceSIF on a Lambertian BRDF → factor-2 SIF₀ injection
        # at m=0. SolarBeam contributions to surface j₀⁻ stay inside
        # `create_surface_layer!` for back-compat; Phase 5c will move them
        # out into the same dispatch table.
        # Source adding retains the solar-only attenuation source and forms
        # SIF derivative vectors before emission is added to the forward field.
        source_adding && prepare_source_adding_surface!(source_workspace,
            added_surface_layer,prepared_sources,brdf,m,pol_type,CoreRT.architecture(model))
        surface_source_contribute!(prepared_sources, brdf, added_surface_layer,
                                   m, pol_type, CoreRT.architecture(model))
        source_adding || surface_source_contribute_lin!(prepared_sources, brdf,
                                       added_surface_layer, added_surface_layer_lin,
                                       m, pol_type, CoreRT.architecture(model), sif_range(layout))

        # Close every lower subcolumn with the same linearized surface layer.
        # Treat the surface as a potentially reflecting added layer even when
        # its current reflectance is zero: a surface-parameter tangent can be
        # nonzero at that point.
        for ims in eachindex(sensor_levels)
            surface_interface = get_scattering_interface(
                bottom_interfaces[ims], true, Nz - sensor_levels[ims] + 1)
            interaction!(surface_interface, SFI,
                         bottom_layers[ims], bottom_layers_lin[ims],
                         added_surface_layer, added_surface_layer_lin, I_static)
        end
        
        if source_adding
            # Surface closure uses the same forward operator and physical
            # surface tangents as matrix adding; atmosphere adds only vectors.
            interaction!(scattering_interfaces_all[end],SFI,composite_layer,
                added_surface_layer,I_static)
            source_incident_fields!(source_workspace,added_surface_layer,I_static)
            source_adding_tangents!(source_workspace,added_surface_layer,
                added_surface_layer_lin,layer_opt_props_lin,τ̇_sum_all,quad_points,arr_type,
                surface_range(layout),sif_range(layout))
        else
            @timeit "interaction" interaction!(scattering_interfaces_all[end], SFI,
                composite_layer, composite_layer_lin,
                added_surface_layer, added_surface_layer_lin, I_static)
        end

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

        postprocessing_vza_ms_lin!(
            RS_type, top_layers, top_layers_lin,
            bottom_layers, bottom_layers_lin,
            pol_type, vza, qp_μ, m, vaz, weight, nSpec,
            level_uw, level_dw, level_uw_lin, level_dw_lin,
            I_static, arr_type)

        if convergence_active
            stop, convergence_passes = _fourier_convergence_step!(
                fourier_convergence, convergence_outputs,
                convergence_snapshots, convergence_passes, m, m_max)
            if stop
                m_used = m
                @info "Linearized Fourier series converged" m_used=m m_max=m_max tolerance=fourier_convergence.tolerance guard_through_m=_fourier_guard_through_m(fourier_convergence, m_max) n_consecutive=fourier_convergence.n_consecutive
                break
            end
        end
    end
    _LAST_FOURIER_M_USED[] = m_used

    # `J₀⁺` is the diffuse SFI field; the historical BOA forward output also
    # carries the attenuated collimated beam when a requested VZA resolves to
    # the solar ordinate. Add that carrier and its optical-depth tangent once
    # (outside the Fourier sum) so linearized endpoint radiances remain equal
    # to the production forward result. Interior records keep the same carrier
    # in their explicit `unscattered_downwelling` fields instead.
    if model.obs_geom.include_boa
        boa_direct = [zeros(FT, length(vza), pol_type.n, nSpec)]
        boa_direct_lin = [zeros(FT, length(vza), pol_type.n, nSpec, Nparams)]
        _set_unscattered_downwelling_lin!(
            boa_direct, boa_direct_lin, [Nz],
            τ_sum_endpoint, τ̇_sum_endpoint,
            F₀, vza, qp_μ, iμ₀, μ₀, pol_type)
        T .+= only(boa_direct)
        Ṫ .+= only(boa_direct_lin)
    end
    
    # Show timing statistics (gated on numerics.verbose; default off).
    model.numerics.verbose && print_timer()
    reset_timer!()

    # Preserve the historical four iterable slots while exposing named
    # per-height radiance/Jacobian records.
    R_out, T_out = _select_observer_endpoints(model, R, T)
    Ṙ_out, Ṫ_out = _select_observer_endpoints(model, Ṙ, Ṫ)
    if model.quad_points.external_solar
        # TOA-only representation: the BOA field would be missing its direct
        # component (external μ₀ has no operator node) — report unavailable.
        T_out = nothing
        Ṫ_out = nothing
    end
    level_type = isempty(sensor_levels) ?
        LevelRadianceLin{FT,Nothing,Nothing,Nothing,Nothing} :
        typeof(LevelRadianceLin(
            model.obs_geom.sensor_altitudes[1], sensor_levels[1],
            level_uw[1], level_dw[1], level_direct_dw[1],
            level_uw_lin[1], level_dw_lin[1], level_direct_dw_lin[1]))
    levels = Vector{level_type}()
    for (i, (height, boundary)) in enumerate(zip(
            model.obs_geom.sensor_altitudes, sensor_levels))
        push!(levels, LevelRadianceLin(
            height, boundary,
            level_uw[i], level_dw[i], level_direct_dw[i],
            level_uw_lin[i], level_dw_lin[i], level_direct_dw_lin[i]))
    end
    return ObserverRTResultLin(
        R_out, T_out, Ṙ_out, Ṫ_out, levels,
        FT(model.obs_geom.toa_altitude), layout)
end
