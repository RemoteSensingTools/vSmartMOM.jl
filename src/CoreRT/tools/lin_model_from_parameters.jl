"Set a common l_max in a given band for each aerosol type"
function set_uniform_lmax!(lmax::Vector{Int}, aerosol_optics, lin_aerosol_optics=nothing)
    n_bands = length(lmax)
    n_aer = length(aerosol_optics[1])
    for i_band = 1:n_bands
        pref_lmax = lmax[i_band]
        for i_aer = 1:n_aer
            aer_lmax = size(aerosol_optics[i_band][i_aer].greek_coefs.β, 1)
            if aer_lmax < pref_lmax
                for fname in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    old = getfield(aerosol_optics[i_band][i_aer].greek_coefs, fname)
                    new_arr = zeros(eltype(old), pref_lmax, size(old)[2:end]...)
                    view(new_arr, 1:aer_lmax, ntuple(_ -> Colon(), ndims(old)-1)...) .= old
                    setfield!(aerosol_optics[i_band][i_aer].greek_coefs, fname, new_arr)
                end
                if lin_aerosol_optics !== nothing
                    for fname in (:α̇, :β̇, :γ̇, :δ̇, :ϵ̇, :ζ̇)
                        old = getfield(lin_aerosol_optics[i_band][i_aer].lin_greek_coefs, fname)
                        new_arr = zeros(eltype(old), size(old,1), pref_lmax, size(old)[3:end]...)
                        view(new_arr, :, 1:aer_lmax,
                             ntuple(_ -> Colon(), ndims(old)-2)...) .= old
                        setfield!(lin_aerosol_optics[i_band][i_aer].lin_greek_coefs,
                                  fname, new_arr)
                    end
                end
            end
        end
    end
end

@inline function _zero_lin_greek(greek, ::Type{FT}) where {FT<:AbstractFloat}
    make_zero(field) = zeros(FT, 4, size(getfield(greek, field))...)
    return linGreekCoefs(make_zero(:α), make_zero(:β), make_zero(:γ),
                         make_zero(:δ), make_zero(:ϵ), make_zero(:ζ))
end

@inline function _zero_microphysics_field(value, ::Type{FT}) where {FT<:AbstractFloat}
    return value isa Number ? zeros(FT, 4) : zeros(FT, 4, size(value)...)
end

"Create structurally valid zero Mie tangents when microphysics is fixed."
function _zero_aerosol_microphysics_jacobian(aerosol, ::Type{FT}) where {FT<:AbstractFloat}
    phase_lin = aerosol.phase_greek === nothing ? nothing :
        [_zero_lin_greek(greek, FT) for greek in aerosol.phase_greek]
    return linAerosolOptics(
        lin_greek_coefs=_zero_lin_greek(aerosol.greek_coefs, FT),
        ω̃̇=_zero_microphysics_field(aerosol.ω̃, FT),
        k̇=_zero_microphysics_field(aerosol.k, FT),
        ḟᵗ=_zero_microphysics_field(aerosol.fᵗ, FT),
        phase_ν=aerosol.phase_ν === nothing ? nothing : FT.(aerosol.phase_ν),
        phase_lin_greek=phase_lin)
end

"""
    model_from_parameters(::LinMode, params::vSmartMOM_Parameters)

Construct both the forward `RTModel` and the linearized `RTModelLin` objects
from the input parameters, for use in linearized (Jacobian) RT computations.

This is the **linearized** counterpart of `model_from_parameters(params)`. It computes:

1. **Forward optical properties** (same as the forward-only model):
   - Rayleigh optical depth per layer and spectral point.
   - Aerosol optical properties (via Mie theory) with δ-M truncation.
   - Trace gas absorption cross-sections and optical depths.

2. **Linearized optical properties** (derivatives of the above):
   - `τ̇_aer[iB][iaer, 7, nSpec, nLayers]`: Derivatives of aerosol τ w.r.t. 7 sub-parameters
     `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]` per aerosol type.
   - `τ̇_abs[iB][NGas*Nz, nSpec, Nz]`: Layer-resolved derivatives of gas
     absorption τ w.r.t. `VMR(gas, z)`, flattened species-major.
   - `lin_aerosol_optics[iB][iaer]`: Derivatives of Mie properties (ω̃, fᵗ, greek coefficients)
     w.r.t. Mie parameters `[nᵣ, nᵢ, rₘ, σᵣ]`.

# Returns
- `model::RTModel`: Forward model (optical properties, geometry, quadrature).
- `lin_model::RTModelLin`: Linearized model (all derivative arrays).

# Notes
- Strict-interior observer heights are inserted as exact atmospheric
  interfaces during model construction and are supported by the elastic
  analytic-linearized solver.
- Aerosol Mie calculations use `ForwardDiff.Dual` numbers to simultaneously obtain
  derivatives of the extinction cross-section, single-scattering albedo, truncation
  factor, and greek coefficients with respect to `[nᵣ, nᵢ, rₘ, σᵣ]`.
- Set `compute_aerosol_microphysics_jacobians=false` when all four Mie
  parameters are fixed. The forward aerosol optics are retained and
  structurally valid zero tangents are supplied to the downstream mixer.
- Set `compute_h2o_jacobians=false` when the q-driven H2O profile is fixed.
  Its forward absorption is retained without generating H2O tangent entries.
- `aerosol_anchor_bands` has the same meaning as in the forward constructor:
  it fixes aerosol spectral endpoint calculations to canonical full-band
  grids while evaluating the interpolated optics on chunked or
  shoulder-expanded solve grids. Forward and tangent optics always share the
  same anchors.
"""
function model_from_parameters(lin::LinMode,
    params::vSmartMOM_Parameters;
    sources::AbstractSource = SolarBeam(),
    external_solar::Bool = true,
    aerosol_anchor_bands=nothing,
    compute_aerosol_microphysics_jacobians::Bool = true,
    compute_h2o_jacobians::Bool = true)
    FT = params.float_type
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
    aerosol_bands = aerosol_anchor_bands === nothing ? params.spec_bands :
                     aerosol_anchor_bands
    length(aerosol_bands) == n_bands || throw(DimensionMismatch(
        "aerosol_anchor_bands must have one grid per spectral band"))
    all(!isempty, aerosol_bands) || throw(ArgumentError(
        "aerosol interpolation anchor grids must be nonempty"))
    scat = params.scattering_params
    abs_params = params.absorption_params
    if scat !== nothing && any(_has_analytic_phase_function, scat.rt_aerosols)
        throw(ArgumentError(
            "model_from_parameters(LinMode(), ...) currently supports Mie aerosols only; " *
            "analytic phase-function aerosols do not yet define Mie-parameter Jacobians."))
    end

    truncation_type = _resolved_truncation(params, params.float_type)

    vmr = isnothing(abs_params) ? Dict() : abs_params.vmr
    T_ft, p_ft, q_ft = convert(Vector{FT}, params.T), convert(Vector{FT}, params.p), convert(Vector{FT}, params.q)
    obs_alt = params.obs_alt isa Real ? FT(params.obs_alt) : convert(Vector{FT}, params.obs_alt)
    profile, observation = prepare_observer_profile(
        T_ft, p_ft, q_ft, vmr, obs_alt, params.profile_reduction_n)
    input_p_full = (p_ft[1:end-1] .+ p_ft[2:end]) ./ FT(2)
    sources = reframe_vertical_sources(sources, input_p_full, profile.p_full)
    obs_geom = ObsGeometry{FT}(
        FT(params.sza), convert(Vector{FT}, params.vza), convert(Vector{FT}, params.vaz), obs_alt,
        observation.sensor_levels, observation.interior_altitudes,
        observation.include_toa, observation.include_boa, observation.toa_altitude)
    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom,
                                 params.polarization_type, array_type(params.architecture);
                                 external_solar)
    rayleigh_molecular_T = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry)

    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(FT, n_bands)
    τ_rayl = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]

    FT2 = isnothing(abs_params) || !haskey(abs_params.vmr, "CO2") ? params.float_type : eltype(abs_params.vmr["CO2"])
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    N_fix_gas = isnothing(abs_params) ? 0 : length(unique(Iterators.flatten(abs_params.fixed_molecules)))
    N_var_gas = isnothing(abs_params) ? 0 : length(unique(Iterators.flatten(abs_params.variable_molecules)))
    # Gas Jacobian columns are species-major and layer-resolved:
    # (gas 1, z=1:Nz), (gas 2, z=1:Nz), ... .  At layer iz only the
    # corresponding gas/layer column is nonzero at the optics boundary.
    N_gas_species = 1 + N_var_gas
    Nz = length(profile.p_full)
    τ̇_abs     = [zeros(FT2, N_gas_species * Nz,
                        length(params.spec_bands[i]), Nz) for i in 1:n_bands]
    τ̇_rayl_psurf = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    τ̇_abs_psurf = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    τ̇_aer_psurf = [zeros(FT2, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    psurf_tangents = psurf_profile_tangents(profile)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)

    for i_band=1:n_bands

        # `params` may have been parsed as Float64 and subsequently switched
        # to Float32 by a caller. Keep the wavelength grid consistent with the
        # reframed profile and Rayleigh arrays, as the forward constructor does.
        curr_band_λ = FT.(FT(1e4) ./ params.spec_bands[i_band])
        νₘ = FT(0.5)*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = FT(1.e7)/νₘ
        # Per-band molecular-constant depolarizations. ϖ_Cabannes is always
        # taken from the molecular path; the depol values feed greek coefs and
        # τ_rayl only when params.depol < 0 (auto). See model_from_parameters.jl
        # for the full rule.
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(rayleigh_molecular_T))
        ϖ_Cabannes[i_band] = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
        γ_air_Cabannes, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
        γ_air_Rayleigh, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
        depol_air_Cabannes = 2γ_air_Cabannes/(1+γ_air_Cabannes)
        depol_air_Rayleigh = 2γ_air_Rayleigh/(1+γ_air_Rayleigh)

        depol_use_Cab = params.depol < 0 ? FT(depol_air_Cabannes) : FT(params.depol)
        depol_use_Ray = params.depol < 0 ? FT(depol_air_Rayleigh) : FT(params.depol)

        push!(greek_cabannes, Scattering.get_greek_rayleigh(depol_use_Cab))
        push!(greek_rayleigh, Scattering.get_greek_rayleigh(depol_use_Ray))

        τ_rayl[i_band]   .= getRayleighLayerOptProp(profile.p_half[end],
                                curr_band_λ,
                                depol_use_Ray, profile.vcd_dry)
        # τ_rayl = τ_total(p_surf) * vcd_dry/sum(vcd_dry).
        dry = profile.vcd_dry
        drydot = psurf_tangents.vcd_dry_dot
        total = vec(sum(τ_rayl[i_band], dims=2))
        totaldot = total ./ profile.p_half[end]
        # Use the factored quotient rule: the direct form squares the total
        # dry column (~1e25), which overflows in Float32 and makes the entire
        # surface-pressure Jacobian NaN.
        frac, fracdot = _normalized_column_fraction_tangent(dry, drydot)
        τ̇_rayl_psurf[i_band] .= totaldot * frac' .+ total * fracdot'

        (isnothing(abs_params) && isnothing(params.q)) && continue

        if !isnothing(params.q) && any(!iszero, params.q)
            jac_idx = 1
            h2o_setting = isnothing(abs_params) || isempty(abs_params.h2o_lut) ?
                          nothing : abs_params.h2o_lut[i_band]
            if h2o_setting === :disabled
                nothing
            elseif h2o_setting !== nothing
                if compute_h2o_jacobians
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], τ̇_abs[i_band], jac_idx,
                        h2o_setting, params.spec_bands[i_band], profile)
                else
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], h2o_setting,
                        params.spec_bands[i_band], profile)
                end
            else
                @timeit "Read HITRAN" lines_h2o = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact("H2O")); FT)
                @debug "Computing profile for water vapor (q-driven) in band #$(i_band)"
                bf  = isnothing(abs_params) ? AtmosphericAbsorption.Voigt() : abs_params.broadening_function
                cef = isnothing(abs_params) ? AtmosphericAbsorption.HumlicekWeideman32() : abs_params.CEF
                wc  = isnothing(abs_params) ? 150 : abs_params.wing_cutoff
                absorption_model = AtmosphericAbsorption.LineByLineModel(lines_h2o;
                    profile = bf,
                    wing_cutoff = wc,
                    cpf = cef,
                    architecture = _to_aa_arch(params.architecture),
                    vmr = 0)
                if compute_h2o_jacobians
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], τ̇_abs[i_band], jac_idx,
                        absorption_model, params.spec_bands[i_band], profile)
                else
                    @timeit "Absorption Coeff H2O" compute_h2o_absorption_profile!(
                        τ_abs[i_band], absorption_model,
                        params.spec_bands[i_band], profile)
                end
            end
        end

        if !isnothing(abs_params)
            if !isempty(abs_params.fixed_molecules[i_band])
                for molec_i in 1:length(abs_params.fixed_molecules[i_band])
                    mol_name = abs_params.fixed_molecules[i_band][molec_i]
                    if isempty(abs_params.luts)
                        @timeit "Read HITRAN" lines = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact(mol_name)); FT)

                        @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"
                        absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                            profile = abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            cpf = abs_params.CEF,
                            architecture = _to_aa_arch(params.architecture),
                            vmr = 0)
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band],
                                absorption_model,
                                params.spec_bands[i_band],
                                profile.vmr[mol_name],
                                profile)
                    else
                        compute_absorption_profile!(τ_abs[i_band],
                            abs_params.luts[i_band][molec_i],
                            params.spec_bands[i_band],
                            profile.vmr[mol_name],
                            profile)
                    end
                end
            end
            if !isempty(abs_params.variable_molecules[i_band])
                # luts[i_band] is parallel to vcat(fixed_molecules[i_band],
                # variable_molecules[i_band]); offset variable indices past fixed.
                lut_offset = length(abs_params.fixed_molecules[i_band])
                for molec_i in 1:length(abs_params.variable_molecules[i_band])
                    mol_name = abs_params.variable_molecules[i_band][molec_i]
                    jac_idx = molec_i + 1
                    if isempty(abs_params.luts)
                        @timeit "Read HITRAN" lines = AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(artifact(mol_name)); FT)
                        @debug "Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)"

                        absorption_model = AtmosphericAbsorption.LineByLineModel(lines;
                            profile = abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            cpf = abs_params.CEF,
                            architecture = _to_aa_arch(params.architecture),
                            vmr = 0)
                        @timeit "Absorption Coeff"  compute_absorption_profile!(
                                τ_abs[i_band],
                                τ̇_abs[i_band],
                                jac_idx,
                                absorption_model,
                                params.spec_bands[i_band],
                                profile.vmr[mol_name],
                                profile)
                    else
                        compute_absorption_profile!(
                            τ_abs[i_band],
                            τ̇_abs[i_band],
                            jac_idx,
                            abs_params.luts[i_band][lut_offset + molec_i],
                            params.spec_bands[i_band],
                            profile.vmr[mol_name],
                            profile)
                    end
                end
            end
        end

        # The pressure tangent holds cross sections, T, q, and VMR fixed for
        # ordinary line absorption, so its bottom-layer dependence is through
        # molecular column only. CIA and MT_CKD also scale with midpoint
        # pressure and receive that additional analytic factor below.
        dry_ratio_dot =
            psurf_tangents.vcd_dry_dot[end] / profile.vcd_dry[end]
        τ̇_abs_psurf[i_band][:, end] .=
            τ_abs[i_band][:, end] .* dry_ratio_dot
        midpoint_pressure_ratio_dot = FT(0.5) / profile.p_full[end]
        binary_ratio_dot = dry_ratio_dot + midpoint_pressure_ratio_dot

        # Collision-induced absorption (HITRAN .cia files), if any. CIA is
        # fixed with respect to the gas-profile Jacobian τ̇_abs, but its
        # surface-pressure tangent is exact at fixed T and composition:
        # τ_CIA ∝ n_midpoint * N_column.
        if !isnothing(abs_params)
            for (cia_i, cia_path) in enumerate(abs_params.cia_files)
                @timeit "CIA $(basename(cia_path))" begin
                    cia_table = _load_configured_cia_table(
                        abs_params, cia_i, params.spec_bands[i_band], FT)
                    τ_cia = zeros(eltype(τ_abs[i_band]), size(τ_abs[i_band]))
                    Absorption.compute_τ_cia!(τ_cia, cia_table, profile,
                                               profile.vmr)
                    τ_abs[i_band] .+= τ_cia
                    @views τ̇_abs_psurf[i_band][:, end] .+=
                        τ_cia[:, end] .* binary_ratio_dot
                end
            end

            # MT_CKD H₂O continuum, if a reference table is configured.
            # It is fixed with respect to the H₂O-profile Jacobian for now,
            # but τ_cont ∝ p_midpoint * N_H2O gives the same exact pressure
            # factor as CIA when q and T are held fixed.
            if !isempty(abs_params.mtckd_file)
                @timeit "MT_CKD H2O continuum" begin
                    mtckd_table = Absorption.load_mtckd(abs_params.mtckd_file)
                    τ_continuum = zeros(
                        eltype(τ_abs[i_band]), size(τ_abs[i_band]))
                    Absorption.compute_τ_h2o_continuum!(
                        τ_continuum, mtckd_table, params.spec_bands[i_band],
                        profile, profile.vmr_h2o)
                    τ_abs[i_band] .+= τ_continuum
                    @views τ̇_abs_psurf[i_band][:, end] .+=
                        τ_continuum[:, end] .* binary_ratio_dot
                end
            end
        end
    end

    aerosol_optics = [Array{AerosolOptics}(undef, (n_aer)) for i=1:n_bands]
    lin_aerosol_optics = [Array{linAerosolOptics}(undef, (n_aer)) for i=1:n_bands]

    FT2 = params.float_type

    τ_aer = [zeros(FT2, n_aer, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]
    τ̇_aer = [zeros(FT2, n_aer, 7, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]

    for i_aer=1:n_aer

        c_aero = scat.rt_aerosols[i_aer]
        τ_ref = c_aero.τ_ref
        curr_aerosol = c_aero.aerosol

        size_distribution = curr_aerosol.size_distribution
        mie_aerosol = Aerosol(size_distribution, curr_aerosol.nᵣ, curr_aerosol.nᵢ)

        # Linearized Mie has NO GPU variant: ForwardDiff/analytic-derivative Greek
        # coefficients are only implemented on the CPU analytic kernel. Force CPU
        # architecture here regardless of params.architecture so a GPU run still
        # builds Jacobians correctly via compute_aerosol_optical_properties(lin, …).
        _mie(λ) = make_mie_model(scat.decomp_type, mie_aerosol, λ,
            params.polarization_type, truncation_type, scat.r_max, scat.nquad_radius;
            architecture = Architectures.CPU())

        # k_ref normalization at λ_ref uses the aerosol's own nᵣ/nᵢ so the Jacobian
        # chain w.r.t. nᵣ/nᵢ stays intact — the FD test perturbs the aerosol's nᵣ, and a
        # fixed-n_ref k_ref would zero ∂k_ref/∂nᵣ. For the common case (n_ref unset →
        # defaults to the aerosol's own index) this equals the forward path's k_ref value.
        # KNOWN LIMITATION: for an explicit n_ref that differs from the aerosol's index,
        # the linearized forward/Jacobian normalization uses the aerosol index here while
        # the non-linear forward uses n_ref — a small inconsistency for that niche config.
        if compute_aerosol_microphysics_jacobians
            mie_model_ref_lin = _mie(scat.λ_ref)
            k_ref, k̇_ref = compute_ref_aerosol_extinction(
                lin, mie_model_ref_lin, params.float_type)
        else
            # The selective OCO plan fixes aerosol microphysics. Match the
            # nonlinear forward constructor exactly: its AOD reference
            # normalization uses the configured common `n_ref`, not each
            # aerosol's in-band refractive index. This distinction matters for
            # species 2/3 and is required for truth/linearized-forward closure.
            ref_aerosol = Aerosol(size_distribution,
                                  real(scat.n_ref), -imag(scat.n_ref))
            mie_model_ref = make_mie_model(
                scat.decomp_type, ref_aerosol, scat.λ_ref,
                params.polarization_type, truncation_type,
                scat.r_max, scat.nquad_radius;
                architecture=Architectures.CPU())
            k_ref = compute_ref_aerosol_extinction(
                mie_model_ref, params.float_type)
            k̇_ref = zeros(FT2, 4)
        end

        function aerosol_pair(λ)
            mie_model = _mie(λ)
            if compute_aerosol_microphysics_jacobians
                return compute_aerosol_optical_properties(lin, mie_model, FT2)
            end
            aerosol = compute_aerosol_optical_properties(mie_model)
            return aerosol, _zero_aerosol_microphysics_jacobian(aerosol, FT2)
        end

        for i_band=1:n_bands

            curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]
            anchor_ν = FT2.(aerosol_bands[i_band])
            anchor_band_λ = FT2(1e4) ./ anchor_ν
            already_truncated = false

            if length(curr_band_λ) == 1 && length(anchor_band_λ) == 1
                @timeit "Mie calc"  aerosol_optics_raw, lin_aerosol_optics_raw =
                                aerosol_pair((maximum(curr_band_λ)+minimum(curr_band_λ))/2)

            else
                # Multi-spectral band: compute Mie properties at band edges (λ[1], λ[end]) and
                # linearly interpolate k_ext, k_sca and their derivatives in wavenumber across
                # the band. Greek coefficients are averaged. This avoids a full Mie calculation
                # at every spectral point (Mie varies smoothly with λ).
                n_spec = length(curr_band_λ)
                length(anchor_band_λ) > 1 || throw(ArgumentError(
                    "a multi-wavelength solve grid requires at least two aerosol anchors"))
                @timeit "Mie calc"  aerosol_optics_raw_0, lin_aerosol_optics_raw_0 =
                                aerosol_pair(anchor_band_λ[1])

                @timeit "Mie calc"  aerosol_optics_raw_1, lin_aerosol_optics_raw_1 =
                                aerosol_pair(anchor_band_λ[end])

                function truncate_endpoint(a, la)
                    if truncation_type isa Scattering.δBGE &&
                       length(a.greek_coefs.β) > truncation_type.l_max
                        return Scattering.truncate_phase(
                            truncation_type, a, la; reportFit=false)
                    end
                    return Scattering.truncate_phase(
                        Scattering.NoTruncation(), a, la)
                end
                endpoint₀, lin_endpoint₀ = truncate_endpoint(
                    aerosol_optics_raw_0, lin_aerosol_optics_raw_0)
                endpoint₁, lin_endpoint₁ = truncate_endpoint(
                    aerosol_optics_raw_1, lin_aerosol_optics_raw_1)
                ν_spec = FT2.(params.spec_bands[i_band])
                ν_ref_phase = FT2(1e4) / scat.λ_ref
                endpoint_ref = lin_endpoint_ref = nothing
                if first(extrema(anchor_ν)) < ν_ref_phase < last(extrema(anchor_ν))
                    @timeit "Mie calc" aerosol_optics_raw_ref, lin_aerosol_optics_raw_ref =
                        aerosol_pair(scat.λ_ref)
                    endpoint_ref, lin_endpoint_ref = truncate_endpoint(
                        aerosol_optics_raw_ref, lin_aerosol_optics_raw_ref)
                end
                aerosol_optics_raw, lin_aerosol_optics_raw =
                    _spectralize_truncated_endpoints(
                        endpoint₀, lin_endpoint₀, endpoint₁, lin_endpoint₁, ν_spec;
                        reference=endpoint_ref, lin_reference=lin_endpoint_ref,
                        ν_ref=ν_ref_phase,
                        ν_endpoints=(first(anchor_ν), last(anchor_ν)))
                already_truncated = true
            end

            # Always go through `truncate_phase`. For NoTruncation this
            # is the cheap passthrough that resets the raw Mie sentinel
            # `fᵗ = 1` to 0 (and `ḟᵗ` to zero) so downstream
            # `delta_m_truncation_lin` doesn't silently kill the
            # aerosol SSA. For δBGE the existing least-squares fit runs
            # only when the Greek expansion is long enough to truncate.
            β_len = size(aerosol_optics_raw.greek_coefs.β, 1)
            if already_truncated
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                lin_aerosol_optics[i_band][i_aer] = lin_aerosol_optics_raw
                l_max_aer[i_aer, i_band] = β_len
            elseif truncation_type isa Scattering.δBGE && β_len > truncation_type.l_max
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] =
                    Scattering.truncate_phase(truncation_type,
                                aerosol_optics_raw, lin_aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
            else
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] =
                    Scattering.truncate_phase(Scattering.NoTruncation(),
                                aerosol_optics_raw, lin_aerosol_optics_raw)
                l_max_aer[i_aer, i_band] = β_len
            end

            k_band = aerosol_optics[i_band][i_aer].k
            k̇_band = lin_aerosol_optics[i_band][i_aer].k̇
            ν_spec = FT2.(params.spec_bands[i_band])
            ν_ref = FT2(1e4) / scat.λ_ref
            aod_scale = _aod_spectral_scale(ν_spec, k_band, k_ref, ν_ref)
            aod_scale_dot = zeros(FT2, 4, length(ν_spec))
            if first(extrema(ν_spec)) < ν_ref < last(extrema(ν_spec))
                for ctr in 1:4
                    q̇₀ = k̇_band[ctr, 1] / k_ref -
                          k_band[1] * k̇_ref[ctr] / k_ref^2
                    q̇₁ = k̇_band[ctr, end] / k_ref -
                          k_band[end] * k̇_ref[ctr] / k_ref^2
                    aod_scale_dot[ctr, :] .= _natural_cubic_three(
                        ν_spec, [first(ν_spec), ν_ref, last(ν_spec)],
                        [q̇₀, zero(FT2), q̇₁])
                end
            else
                for ctr in 1:4
                    aod_scale_dot[ctr, :] .= k̇_band[ctr, :] ./ k_ref .-
                        k_band .* k̇_ref[ctr] ./ k_ref^2
                end
            end

            # Match the production forward profile discretization exactly.
            # Normal profiles use pressure; altitude-form LogNormal profiles
            # use exact CDF differences at geometric layer interfaces.
            τₚ, dτₚdp₀, dτₚdσp =
                getAerosolLayerOptProp(lin, 1, c_aero.profile, profile)
            dτₚdpsurf = aerosol_profile_psurf_tangent(
                c_aero.profile, profile, psurf_tangents.Δz_dot)

            # ────────────────────────────────────────────────────────────────
            # Aerosol optical depth per layer:
            #   τ_aer(λ,z) = (τ_ref / k_ref) · k(λ) · τₚ(z)
            #
            # Jacobian parameters [1..7]:
            #   1: τ_ref    → ∂τ_aer/∂τ_ref = (k/k_ref) · τₚ
            #   2-5: nᵣ,nᵢ,rₘ,σᵣ → quotient rule on k(λ)/k_ref:
            #        ∂τ_aer/∂p = τ_ref·(dk/dp/k_ref − k·dk_ref/dp/k_ref²)·τₚ
            #   6-7: p₀,σₚ  → only τₚ(z) depends on these:
            #        ∂τ_aer/∂p = (τ_ref/k_ref)·k(λ)·∂τₚ/∂p
            # ────────────────────────────────────────────────────────────────
            τ_aer[i_band][i_aer,:,:] = τ_ref .* aod_scale .* τₚ'
            τ̇_aer_psurf[i_band][i_aer, :, :] .=
                τ_ref .* aod_scale .* dτₚdpsurf'

            τ̇_aer[i_band][i_aer,1,:,:] .=
                aod_scale .* τₚ'

            for ctr=1:4
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    τ_ref .* aod_scale_dot[ctr, :] .* τₚ'
            end
            for ctr=5:6
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    τ_ref .* aod_scale .*
                    (ctr==5 ? dτₚdp₀' : dτₚdσp')
            end

        end
    end
    # Mirror forward path: with no aerosols, fall back to the legacy l_trunc.
    for i_band = 1:n_bands
        l_max[i_band] = n_aer > 0 ? maximum(l_max_aer[:,i_band]) : params.l_trunc
    end
    if n_aer > 0
        set_uniform_lmax!(l_max, aerosol_optics, lin_aerosol_optics)
    end

    # Per-band Fourier loop bound (order). Phase C: trait-based aggregator
    # via `component_m_max(c, ctx)`. Same helper as the forward path so
    # forward and lin can never silently disagree.
    components_per_band = [_band_components(params, aerosol_optics, sources, i_band;
                                rayleigh_active=_rayleigh_active(τ_rayl[i_band]))
                            for i_band in 1:n_bands]
    m_max_bands = _derive_m_max_bands_via_traits(l_max, params.max_m,
                                                  components_per_band,
                                                  quad_points.Nstreams;
                                                  greek_beta_cutoff=params.greek_beta_cutoff)
    n_fourier_moments_bands = m_max_bands .+ 1

    # Build the hierarchical RTModel
    solver = SolverConfig{FT, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        m_max_bands,
        n_fourier_moments_bands,
        l_max,
        params.l_trunc,
        FT(params.Δ_angle),
        FT(params.depol),
        true,   # use_component_traits — flipped on in Phase C
    )
    spec_bands_ft = [convert(Vector{FT}, b) for b in params.spec_bands]
    atm = Atmosphere(profile, spec_bands_ft)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    # τ_abs may have been allocated at FT2 (= eltype of the VMR for AD)
    # which can differ from FT (= params.float_type, e.g. Float32 with
    # Float64 VMR). Optics requires both τ-vectors to share the same
    # element type, so coerce τ_abs back to FT here. Matches the forward
    # path's behavior of always allocating τ_abs at params.float_type.
    τ_abs_ft = eltype(eltype(τ_abs)) === FT ? τ_abs : [convert(Matrix{FT}, m) for m in τ_abs]
    optics = Optics(rayleigh_s, aerosols_s, τ_abs_ft, τ_rayl)
    numerics = _convert_numerics(params.numerics, FT)
    model = RTModel(params.architecture, solver, numerics, obs_geom, quad_points, atm, optics, params.brdf, sources)
    return model, RTModelLin(τ̇_abs, τ̇_aer, lin_aerosol_optics,
                             τ̇_rayl_psurf, τ̇_aer_psurf, τ̇_abs_psurf)
end

"""
    model_from_parameters_lin(params)

Convenience alias for `model_from_parameters(LinMode(), params)`.
Returns `(model, lin_model)`.
"""
model_from_parameters_lin(params::vSmartMOM_Parameters) =
    model_from_parameters(LinMode(), params)
