#=

This file contains the `model_from_parameters` function, which computes all derived information
like optical thicknesses, from the input parameters. Produces an RTModel object.

=#

"Generate default set of parameters for Radiative Transfer calculations (from ModelParameters/)"
default_parameters() = vSmartMOM.IO.parameters_from_yaml(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))

"""
    _resolved_truncation(params, FT)

Return `params.truncation` typed for `FT`. If a downstream caller has
mutated the legacy `l_trunc` / `Δ_angle` fields after construction (an
older idiom — e.g. `test_float32.jl`), they must also mutate
`params.truncation` to keep the two in sync; this function does NOT
silently rebuild from the legacy fields, because that would discard a
user-set explicit truncation (e.g. `NoTruncation()` or a custom δBGE).
"""
@inline _resolved_truncation(params, ::Type{FT}) where {FT} = params.truncation

_has_analytic_phase_function(c_aero::RT_Aerosol) =
    c_aero.phase_function !== nothing

function _finite_truncation_lmax(params, truncation_type)
    fallback = max(1, Int(params.l_trunc), 2 * Int(params.max_m) - 1)
    hasproperty(truncation_type, :l_max) || return fallback
    trunc_lmax = Int(getproperty(truncation_type, :l_max))
    return trunc_lmax == typemax(Int) ? fallback : max(1, trunc_lmax)
end

function _analytic_phase_lmax(params, truncation_type)
    trunc_lmax = _finite_truncation_lmax(params, truncation_type)
    base_lmax = max(1, trunc_lmax, 2 * Int(params.max_m) - 1)
    if truncation_type isa Scattering.δBGE
        return max(base_lmax + 1, 2 * trunc_lmax)
    end
    return base_lmax
end

function _analytic_aerosol_optics(c_aero::RT_Aerosol, params, truncation_type,
                                  ::Type{FT}) where {FT}
    lmax = _analytic_phase_lmax(params, truncation_type)
    return Scattering.analytic_aerosol_optics(
        c_aero.phase_function;
        single_scattering_albedo = convert(FT, c_aero.ϖ),
        extinction_cross_section = one(FT),
        l_max = lmax,
        nquad = max(2lmax + 1, 64))
end

"""
    model_from_parameters(params::vSmartMOM_Parameters) -> RTModel

Construct an [`RTModel`](@ref) from user-supplied [`vSmartMOM_Parameters`](@ref).

Computes all derived quantities needed by the RT solver:
- Observation geometry and quadrature points
- Rayleigh scattering coefficients (Greek/Cabannes)
- Aerosol optics via Mie theory with δ-M truncation
- Gas absorption cross-sections and optical depth profiles
- Per-band Fourier truncation limits

# Arguments
- `params::vSmartMOM_Parameters`: Configuration struct (typically from `parameters_from_yaml`).

# Returns
- `model::RTModel{ARCH, FT}`: Hierarchical model ready for `rt_run(model)`.

# See also
- `model_from_parameters(LinMode(), params)` for the linearized (Jacobian) variant.
- `parameters_from_yaml(path)` to load parameters from a YAML file.
"""
function model_from_parameters(params::vSmartMOM_Parameters;
                               sources::AbstractSource = SolarBeam())
    FT = params.float_type
    #@show FT
    # Number of total bands and aerosols (for convenience)
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Create observation geometry
    obs_geom = ObsGeometry{FT}(params.sza, params.vza, params.vaz, params.obs_alt)

    # Truncation method (typed; NoTruncation, δBGE, ...). The legacy
    # `params.Δ_angle` is only consulted via the default δBGE built in
    # `parameters_from_dict` when the user has not set `truncation`
    # explicitly. If a downstream caller mutates `l_trunc` / `Δ_angle`
    # after construction, they must also mutate `params.truncation` —
    # see `_resolved_truncation`.
    truncation_type = _resolved_truncation(params, FT)
    # Set quadrature points for streams
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Get AtmosphericProfile from parameters.
    # Convert T, p, q to the requested FT so that profile arrays (vcd_dry,
    # p_half, etc.) are consistent with the model float type.  Without this,
    # mutating params.float_type = Float32 on a Float64-constructed
    # vSmartMOM_Parameters causes getRayleighLayerOptProp to receive mixed
    # Float32/Float64 arguments and raises a MethodError.
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    T_ft, p_ft, q_ft = convert(Vector{FT}, params.T), convert(Vector{FT}, params.p), convert(Vector{FT}, params.q)
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(T_ft, p_ft, q_ft, vmr)

    profile = AtmosphericProfile(T_ft, p_full, q_ft, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz)

    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    # Rayleigh optical depth per spectral point per layer (uses reduced profile size).
    #
    # Depolarization sourcing rule (2026-04-24):
    #   params.depol < 0  → "auto": derive both Rayleigh (greek_rayleigh, τ_rayl)
    #                       and Cabannes (greek_cabannes) depolarizations per
    #                       band from the N₂/O₂ molecular constants (sanghavi
    #                       convention; appropriate for Earth atmospheres).
    #   params.depol ≥ 0  → "explicit": use params.depol uniformly for
    #                       greek_rayleigh, greek_cabannes, AND τ_rayl. Required
    #                       for idealizations such as Natraj (1980) which fix
    #                       depol = 0 — no molecular path can reproduce that.
    # ϖ_Cabannes is always taken from the molecular path (it relates Cabannes
    # single-scattering albedo to N₂/O₂ line ratios, not depol).
    τ_rayl = [zeros(FT,length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands];

    # Per-band Cabannes / Rayleigh greek coefs (depolarizations from molecular constants)
    greek_rayleigh_arr = Vector{Scattering.GreekCoefs{FT}}()
    greek_cabannes = Vector{Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(FT, n_bands)
    
    τ_abs     = [zeros(FT, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    
    # Track per-band l_max from aerosol greek coef lengths
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)
    
    # Loop over all bands:
    for i_band=1:n_bands

        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])
        
        # Per-band molecular-constant depolarizations (always needed for ϖ_Cab;
        # used as the depol values when params.depol < 0).
        νₘ = FT(0.5) * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = FT(1.0e7) / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(300))
        ϖ_Cab = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        γ_air_Ray, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
        ϖ_Cabannes[i_band] = FT(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        depol_air_Ray = 2γ_air_Ray / (1 + γ_air_Ray)

        # Apply the auto vs. explicit depol rule (see comment above the loop).
        depol_use_Cab = params.depol < 0 ? FT(depol_air_Cab) : FT(params.depol)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Ray) : FT(params.depol)

        push!(greek_cabannes,     Scattering.get_greek_rayleigh(depol_use_Cab))
        push!(greek_rayleigh_arr, Scattering.get_greek_rayleigh(depol_use_Ray))

        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end],
                                curr_band_λ,
                                depol_use_Ray, profile.vcd_dry);
        #@show τ_rayl[i_band]
        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        ap = params.absorption_params

        # Loop over fixed_molecules ∪ variable_molecules in this band; H2O is
        # handled separately below (driven by q, not by these lists).
        all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])
        for (molec_i, mol_name) in enumerate(all_species)
            if isempty(ap.luts)
                @timeit "Read HITRAN" hitran_data = read_hitran(artifact(mol_name), iso=-1)
                println("Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)")
                absorption_model = make_hitran_model(hitran_data,
                    ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    CEF = ap.CEF,
                    architecture = params.architecture,
                    vmr = 0)
                @timeit "Absorption Coeff" compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band], profile.vmr[mol_name], profile)
            else
                compute_absorption_profile!(τ_abs[i_band], ap.luts[i_band][molec_i], params.spec_bands[i_band], profile.vmr[mol_name], profile)
            end
        end

        # H₂O line absorption (driven by q). Use the band's H2O LUT if the
        # parser found one inside LUTfiles; otherwise fall back to artifact.
        if any(!iszero, params.q)
            if ap.h2o_lut[i_band] !== nothing
                @timeit "Absorption Coeff H2O" compute_absorption_profile!(
                    τ_abs[i_band], ap.h2o_lut[i_band],
                    params.spec_bands[i_band], profile.vmr_h2o, profile)
            else
                @timeit "Read HITRAN H2O" hitran_data = read_hitran(artifact("H2O"), iso=-1)
                println("Computing profile for H2O (q-driven) for band #$(i_band)")
                h2o_model = make_hitran_model(hitran_data,
                    ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    CEF = ap.CEF,
                    architecture = params.architecture,
                    vmr = 0)
                @timeit "Absorption Coeff H2O" compute_absorption_profile!(τ_abs[i_band], h2o_model, params.spec_bands[i_band], profile.vmr_h2o, profile)
            end
        end

        # Collision-induced absorption (HITRAN .cia files), if any.
        for cia_path in ap.cia_files
            @timeit "CIA $(basename(cia_path))" begin
                cia_table = Absorption.load_cia_table(cia_path,
                                                     params.spec_bands[i_band];
                                                     FT = FT)
                Absorption.compute_τ_cia!(τ_abs[i_band], cia_table, profile,
                                           ap.vmr)
            end
        end

        # MT_CKD H₂O continuum (self + foreign), if a reference table is configured.
        if !isempty(ap.mtckd_file)
            @timeit "MT_CKD H2O continuum" begin
                mtckd_table = Absorption.load_mtckd(ap.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(τ_abs[i_band], mtckd_table,
                                                     params.spec_bands[i_band], profile,
                                                     profile.vmr_h2o)
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];

    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]

        if _has_analytic_phase_function(c_aero)
            τ_profile = getAerosolLayerOptProp(one(FT), c_aero.profile, profile)
            for i_band=1:n_bands
                aerosol_optics_raw =
                    _analytic_aerosol_optics(c_aero, params, truncation_type, FT)
                β_len = length(aerosol_optics_raw.greek_coefs.β)
                aerosol_optics[i_band][i_aer] =
                    if truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                        Scattering.truncate_phase_lowconf(truncation_type,
                                                          aerosol_optics_raw; reportFit=false)
                    else
                        Scattering.truncate_phase(Scattering.NoTruncation(),
                                                  aerosol_optics_raw)
                    end
                l_max_aer[i_aer, i_band] =
                    length(aerosol_optics[i_band][i_aer].greek_coefs.β)
                τ_aer[i_band][i_aer, :] =
                    c_aero.τ_ref * τ_profile
                @info "AOD at band $i_band : $(sum(τ_aer[i_band][i_aer,:])), analytic phase function = $(typeof(c_aero.phase_function)), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"
            end
            continue
        end

        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #@show typeof(curr_aerosol.nᵣ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)   
        mie_model.aerosol.nᵣ = real(params.scattering_params.n_ref)
        mie_model.aerosol.nᵢ = -imag(params.scattering_params.n_ref)
        # k for reference wavelength
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)
        
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = FT.(1e4 ./ params.spec_bands[i_band])

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)
            n_ref = params.scattering_params.n_ref
            k = compute_ref_aerosol_extinction(mie_model,  params.float_type)
            
            #@show k
            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 
            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT);
            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            #@show i_aer, i_band
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)
            #aerosol_optics[i_band][i_aer] =  aerosol_optics_raw
            
            # Track greek coef length for l_max computation
            l_max_aer[i_aer, i_band] = min(length(aerosol_optics[i_band][i_aer].greek_coefs.β), 
                                            truncation_type.l_max)

            #@show aerosol_optics[i_band][i_aer].fᵗ
            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                getAerosolLayerOptProp(1, c_aero.profile, profile)
            @info "AOD at band $i_band : $(sum(τ_aer[i_band][i_aer,:])), truncation factor = $(aerosol_optics[i_band][i_aer].fᵗ)"
        end 
    end

    # Compute per-band l_max and max_m from aerosol greek coefficient lengths
    l_max = zeros(Int, n_bands)
    max_m_bands = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
        max_m_bands[i_band] = Int(ceil((l_max[i_band] + 1) / 2))
        # Clamp to user-specified max_m if it's lower (hard cutoff)
        max_m_bands[i_band] = min(max_m_bands[i_band], params.max_m)
    end

    # Build the hierarchical RTModel
    FT_ = FT
    solver = SolverConfig{FT_, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        params.max_m,
        max_m_bands,
        l_max,
        params.l_trunc,
        FT_(params.Δ_angle),
        FT_(params.depol),
    )
    spec_bands_ft = [convert(Vector{FT_}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft)
    rayleigh = RayleighScattering(greek_rayleigh_arr, greek_cabannes, FT_.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh, aerosols_s, τ_abs, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT_)
    return RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
end

"Re-type the user-supplied `RTNumericalParameters` to the resolved
RTModel float type. No-op when types already match."
@inline _convert_numerics(n::RTNumericalParameters{FT}, ::Type{FT}) where {FT} = n
@inline function _convert_numerics(n::RTNumericalParameters, ::Type{FT}) where {FT}
    RTNumericalParameters{FT}(
        dτ_max_threshold = FT(n.dτ_max_threshold),
        dτ_min_floor     = FT(n.dτ_min_floor),
        blas_threads     = n.blas_threads,
    )
end


#=

Modified version for vibrational Raman scattering

=#


"Take the parameters specified in the vSmartMOM_Parameters struct, and calculate derived attributes into an RTModel"
function model_from_parameters(RS_type::Union{VS_0to1_plus, VS_1to0_plus},
                    λ₀,
                    params::vSmartMOM_Parameters;
                    sources::AbstractSource = SolarBeam())
    # Number of total bands and aerosols (for convenience)
    n_bands = 3 #length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)

    # Create observation geometry
    obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt)

    # Create truncation type
    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    # Set quadrature points for streams
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    # Get AtmosphericProfile from parameters.
    # Convert T, p, q to params.float_type so profile arrays are type-consistent.
    FT_vrs_early = params.float_type
    vmr = isnothing(params.absorption_params) ? Dict() : params.absorption_params.vmr
    T_ft, p_ft, q_ft = convert(Vector{FT_vrs_early}, params.T), convert(Vector{FT_vrs_early}, params.p), convert(Vector{FT_vrs_early}, params.q)
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(T_ft, p_ft, q_ft, vmr)
    profile = AtmosphericProfile(T_ft, p_full, q_ft, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz)

    # Reduce the profile to the number of target layers (if specified)
    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile);
    end

    effT = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry);
    # Define RS type
    # Compute N2 and O2
    RS_type.n2, RS_type.o2 = 
        InelasticScattering.getRamanAtmoConstants(1.e7/λ₀,effT);
    InelasticScattering.getRamanSSProp!(RS_type, λ₀);
    n_bands = length(RS_type.iBand)
    params.spec_bands = RS_type.grid_in

    # Rayleigh optical properties calculation
    greek_rayleigh = Scattering.get_greek_rayleigh(params.depol)
    τ_rayl = [zeros(params.float_type,1, length(profile.p_full)) for i=1:n_bands];

    # Per-band Cabannes / Rayleigh depolarization (for inelastic scattering support)
    FT_vrs = params.float_type
    greek_cabannes = typeof(greek_rayleigh)[]
    ϖ_Cabannes = zeros(FT_vrs, n_bands)
    l_max_aer = zeros(Int, max(n_aer, 1), n_bands)

    # Pre-allocated absorption arrays
    τ_abs     = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    # Loop over all bands:
    for i_band=1:n_bands
        # i'th spectral band (convert from cm⁻¹ to μm)
        curr_band_λ = 1e4 ./ params.spec_bands[i_band]

        # Compute per-band Cabannes properties.
        # Explicit (λ₀, n2, o2) form (effT = 300 K, Earth atmospheres).
        νₘ = 0.5 * (params.spec_bands[i_band][1] + params.spec_bands[i_band][end])
        λₘ = 1.0e7 / νₘ
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(1.0e7 / λₘ, 300.0)
        ϖ_Cab = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        ϖ_Cabannes[i_band] = FT_vrs(ϖ_Cab)
        depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
        push!(greek_cabannes, Scattering.get_greek_rayleigh(FT_vrs(depol_air_Cab)))

        # Compute Rayleigh properties per layer for `i_band` band center
        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end], 
                            (maximum(curr_band_λ) + minimum(curr_band_λ))/2, 
                            params.depol, profile.vcd_dry);

        # If no absorption, continue to next band
        isnothing(params.absorption_params) && continue
        ap = params.absorption_params

        # Loop over fixed_molecules ∪ variable_molecules in this band; H2O is
        # handled separately below (driven by q, not by these lists).
        all_species = vcat(ap.fixed_molecules[i_band], ap.variable_molecules[i_band])
        for (molec_i, mol_name) in enumerate(all_species)
            if isempty(ap.luts)
                @timeit "Read HITRAN" hitran_data = read_hitran(artifact(mol_name), iso=-1)
                println("Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)")
                absorption_model = make_hitran_model(hitran_data,
                    ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    CEF = ap.CEF,
                    architecture = params.architecture,
                    vmr = 0)
                @timeit "Absorption Coeff" compute_absorption_profile!(τ_abs[i_band], absorption_model, params.spec_bands[i_band], profile.vmr[mol_name], profile)
            else
                compute_absorption_profile!(τ_abs[i_band], ap.luts[i_band][molec_i], params.spec_bands[i_band], profile.vmr[mol_name], profile)
            end
        end

        # H₂O line absorption (driven by q). Use the band's H2O LUT if the
        # parser found one inside LUTfiles; otherwise fall back to artifact.
        if any(!iszero, params.q)
            if ap.h2o_lut[i_band] !== nothing
                @timeit "Absorption Coeff H2O" compute_absorption_profile!(
                    τ_abs[i_band], ap.h2o_lut[i_band],
                    params.spec_bands[i_band], profile.vmr_h2o, profile)
            else
                @timeit "Read HITRAN H2O" hitran_data = read_hitran(artifact("H2O"), iso=-1)
                println("Computing profile for H2O (q-driven) for band #$(i_band)")
                h2o_model = make_hitran_model(hitran_data,
                    ap.broadening_function,
                    wing_cutoff = ap.wing_cutoff,
                    CEF = ap.CEF,
                    architecture = params.architecture,
                    vmr = 0)
                @timeit "Absorption Coeff H2O" compute_absorption_profile!(τ_abs[i_band], h2o_model, params.spec_bands[i_band], profile.vmr_h2o, profile)
            end
        end

        # Collision-induced absorption (HITRAN .cia files), if any.
        for cia_path in ap.cia_files
            @timeit "CIA $(basename(cia_path))" begin
                cia_table = Absorption.load_cia_table(cia_path,
                                                     params.spec_bands[i_band];
                                                     FT = FT)
                Absorption.compute_τ_cia!(τ_abs[i_band], cia_table, profile,
                                           ap.vmr)
            end
        end

        # MT_CKD H₂O continuum (self + foreign), if a reference table is configured.
        if !isempty(ap.mtckd_file)
            @timeit "MT_CKD H2O continuum" begin
                mtckd_table = Absorption.load_mtckd(ap.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(τ_abs[i_band], mtckd_table,
                                                     params.spec_bands[i_band], profile,
                                                     profile.vmr_h2o)
            end
        end
    end

    # aerosol_optics[iBand][iAer]
    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands];
        
    # τ_aer[iBand][iAer,iZ]
    τ_aer = [zeros(FT, n_aer, length(profile.p_full)) for i=1:n_bands];

    # Loop over aerosol type
    for i_aer=1:n_aer

        # Get curr_aerosol
        c_aero = params.scattering_params.rt_aerosols[i_aer]
        curr_aerosol = c_aero.aerosol
        
        # Create Aerosol size distribution for each aerosol species
        size_distribution = curr_aerosol.size_distribution

        # Create a univariate aerosol distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)
        #mie_aerosol = make_mie_aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ, params.scattering_params.r_max, params.scattering_params.nquad_radius) #Suniti: why is the refractive index needed here?

        # Create the aerosol extinction cross-section at the reference wavelength:
        mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                        mie_aerosol, 
                                        params.scattering_params.λ_ref, 
                                        params.polarization_type, 
                                        truncation_type, 
                                        params.scattering_params.r_max, 
                                        params.scattering_params.nquad_radius)
        mie_model.aerosol.nᵣ = real(params.scattering_params.n_ref)
        mie_model.aerosol.nᵢ = -imag(params.scattering_params.n_ref)
        k_ref          = compute_ref_aerosol_extinction(mie_model, params.float_type)

        #params.scattering_params.rt_aerosols[i_aer].p₀, params.scattering_params.rt_aerosols[i_aer].σp
        # Loop over bands
        for i_band=1:n_bands
            
            # i'th spectral band (convert from cm⁻¹ to μm)
            curr_band_λ = 1e4 ./ params.spec_bands[i_band]

            # Create the aerosols:
            mie_model      = make_mie_model(params.scattering_params.decomp_type, 
                                            mie_aerosol, 
                                            (maximum(curr_band_λ)+minimum(curr_band_λ))/2, 
                                            params.polarization_type, 
                                            truncation_type, 
                                            params.scattering_params.r_max, 
                                            params.scattering_params.nquad_radius)

            # Compute raw (not truncated) aerosol optical properties (not needed in RT eventually) 

            @timeit "Mie calc"  aerosol_optics_raw = compute_aerosol_optical_properties(mie_model, FT2);

            # Compute truncated aerosol optical properties (phase function and fᵗ), consistent with Ltrunc:
            aerosol_optics[i_band][i_aer] = Scattering.truncate_phase(truncation_type, 
                                                    aerosol_optics_raw; reportFit=false)

            # Compute nAer aerosol optical thickness profiles
            τ_aer[i_band][i_aer,:] = 
                params.scattering_params.rt_aerosols[i_aer].τ_ref * 
                (aerosol_optics[i_band][i_aer].k/k_ref) * 
                CoreRT.getAerosolLayerOptProp(1.0, c_aero.p₀, c_aero.σp, profile)
        end 
    end

    # Compute per-band l_max and max_m
    l_max = zeros(Int, n_bands)
    max_m_bands = zeros(Int, n_bands)
    for i_band = 1:n_bands
        if n_aer > 0
            l_max[i_band] = maximum(l_max_aer[:, i_band])
        else
            l_max[i_band] = params.l_trunc
        end
        max_m_bands[i_band] = Int(ceil((l_max[i_band] + 1) / 2))
        max_m_bands[i_band] = min(max_m_bands[i_band], params.max_m)
    end

    # Build the hierarchical RTModel
    FT_vrs2 = params.float_type
    solver = SolverConfig{FT_vrs2, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        params.max_m,
        max_m_bands,
        l_max,
        params.l_trunc,
        FT_vrs2(params.Δ_angle),
        FT_vrs2(params.depol),
    )
    spec_bands_ft2 = [convert(Vector{FT_vrs2}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft2)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT_vrs2.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh_s, aerosols_s, τ_abs, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT_vrs2)
    return RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
end

function loadAbsco(file; scale=(1.0))
    absco = Dataset(file)
    mol = absco["Gas_Index"][1]
    
    cs_name = "Gas_"* mol * "_Absorption"
    # Loading cross sections:
    σ = Float32(scale)*absco[cs_name][:]
    # Temperature
    T = absco["Temperature"][:]
    p = absco["Pressure"][:]/100
    ν = absco["Wavenumber"][:]
    return Absorption.AbscoTable(parse(Int,mol), -1, ν, σ, p, T )
end
