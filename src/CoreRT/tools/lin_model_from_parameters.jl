"Set a common l_max in a given band for each aerosol type"
function set_uniform_lmax!(lmax::Vector{Int}, aerosol_optics)
    n_bands = length(lmax)
    n_aer = length(aerosol_optics[1])
    for i_band = 1:n_bands
        pref_lmax = lmax[i_band]
        for i_aer = 1:n_aer
            aer_lmax = length(aerosol_optics[i_band][i_aer].greek_coefs.β)
            if aer_lmax < pref_lmax
                for fname in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    old = getfield(aerosol_optics[i_band][i_aer].greek_coefs, fname)
                    new_arr = zeros(eltype(old), pref_lmax)
                    new_arr[1:aer_lmax] .= old
                    setfield!(aerosol_optics[i_band][i_aer].greek_coefs, fname, new_arr)
                end
            end
        end
    end
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
   - `τ̇_abs[iB][NGas, nSpec, nLayers]`: Derivatives of gas absorption τ w.r.t. VMR.
   - `lin_aerosol_optics[iB][iaer]`: Derivatives of Mie properties (ω̃, fᵗ, greek coefficients)
     w.r.t. Mie parameters `[nᵣ, nᵢ, rₘ, σᵣ]`.

# Returns
- `model::RTModel`: Forward model (optical properties, geometry, quadrature).
- `lin_model::RTModelLin`: Linearized model (all derivative arrays).

# Notes
- The atmospheric profile may be truncated to the observer altitude for tower/airborne sensors.
- Aerosol Mie calculations use `ForwardDiff.Dual` numbers to simultaneously obtain
  derivatives of the extinction cross-section, single-scattering albedo, truncation
  factor, and greek coefficients with respect to `[nᵣ, nᵢ, rₘ, σᵣ]`.
"""
function model_from_parameters(lin::LinMode,
    params::vSmartMOM_Parameters)
    FT = params.float_type
    n_bands = length(params.spec_bands)
    n_aer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
    scat = params.scattering_params
    abs_params = params.absorption_params

    obs_geom = ObsGeometry(params.sza, params.vza, params.vaz, params.obs_alt)

    truncation_type = Scattering.δBGE{params.float_type}(params.l_trunc, params.Δ_angle)

    quad_points = rt_set_streams(params.quadrature_type, params.l_trunc, obs_geom, params.polarization_type, array_type(params.architecture))

    vmr = isnothing(abs_params) ? Dict() : abs_params.vmr
    p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz = compute_atmos_profile_fields(params.T, params.p, params.q, vmr)

    profile = AtmosphericProfile(params.T, p_full, params.q, p_half, vmr_h2o, vcd_dry, vcd_h2o, new_vmr, Δz)

    if params.profile_reduction_n != -1
        profile = reduce_profile(params.profile_reduction_n, profile)
    end

    greek_cabannes = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    greek_rayleigh = Vector{vSmartMOM.Scattering.GreekCoefs{FT}}()
    ϖ_Cabannes = zeros(n_bands)
    τ_rayl = [zeros(params.float_type, length(params.spec_bands[i]), length(profile.p_full)) for i=1:n_bands]

    FT2 = isnothing(abs_params) || !haskey(abs_params.vmr, "CO2") ? params.float_type : eltype(abs_params.vmr["CO2"])
    τ_abs     = [zeros(FT2, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    N_fix_gas = length(unique(Iterators.flatten(abs_params.fixed_molecules)))
    N_var_gas = length(unique(Iterators.flatten(abs_params.variable_molecules)))
    τ̇_abs     = [zeros(FT2, 1+N_var_gas, length(params.spec_bands[i]), length(profile.p_full)) for i in 1:n_bands]
    max_m = zeros(Int, n_bands)
    l_max = zeros(Int, n_bands)
    l_max_aer = zeros(Int, n_aer, n_bands)

    for i_band=1:n_bands

        curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]
        νₘ = FT(0.5)*(params.spec_bands[i_band][1]+params.spec_bands[i_band][end])
        λₘ = FT(1.e7)/νₘ
        # Per-band molecular-constant depolarizations. ϖ_Cabannes is always
        # taken from the molecular path; the depol values feed greek coefs and
        # τ_rayl only when params.depol < 0 (auto). See model_from_parameters.jl
        # for the full rule.
        _n2, _o2 = InelasticScattering.getRamanAtmoConstants(FT(1.0e7) / λₘ, FT(300))
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

        (isnothing(abs_params) && isnothing(params.q)) && continue

        if !isnothing(params.q) && any(!iszero, params.q)
            jac_idx = 1
            if !isnothing(abs_params) && !isempty(abs_params.h2o_lut) && abs_params.h2o_lut[i_band] !== nothing
                @timeit "Absorption Coeff H2O"  compute_absorption_profile!(
                    τ_abs[i_band],
                    τ̇_abs[i_band],
                    jac_idx,
                    abs_params.h2o_lut[i_band],
                    params.spec_bands[i_band],
                    profile.vmr_h2o,
                    profile)
            else
                @timeit "Read HITRAN" hitran_data = read_hitran(artifact("H2O"), iso=-1)
                println("Computing profile for water vapor (q-driven) in band #$(i_band)")
                bf  = isnothing(abs_params) ? vSmartMOM.Absorption.Voigt() : abs_params.broadening_function
                cef = isnothing(abs_params) ? vSmartMOM.Absorption.HumlicekWeidemann32SDErrorFunction() : abs_params.CEF
                wc  = isnothing(abs_params) ? 150 : abs_params.wing_cutoff
                absorption_model = make_hitran_model(hitran_data, bf,
                    wing_cutoff = wc,
                    CEF = cef,
                    architecture = params.architecture,
                    vmr = 0)
                @timeit "Absorption Coeff H2O"  compute_absorption_profile!(
                    τ_abs[i_band],
                    τ̇_abs[i_band],
                    jac_idx,
                    absorption_model,
                    params.spec_bands[i_band],
                    profile.vmr_h2o,
                    profile)
            end
        end
        if !isnothing(abs_params)
            if !isempty(abs_params.fixed_molecules[i_band])
                for molec_i in 1:length(abs_params.fixed_molecules[i_band])
                    mol_name = abs_params.fixed_molecules[i_band][molec_i]
                    if isempty(abs_params.luts)
                        @timeit "Read HITRAN" hitran_data =
                            read_hitran(artifact(mol_name), iso=-1)

                        println("Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)")
                        absorption_model = make_hitran_model(hitran_data,
                            abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            CEF = abs_params.CEF,
                            architecture = params.architecture,
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
                        @timeit "Read HITRAN" hitran_data =
                            read_hitran(artifact(mol_name), iso=-1)
                        println("Computing profile for $(mol_name) with vmr $(profile.vmr[mol_name]) for band #$(i_band)")

                        absorption_model = make_hitran_model(hitran_data,
                            abs_params.broadening_function,
                            wing_cutoff = abs_params.wing_cutoff,
                            CEF = abs_params.CEF,
                            architecture = params.architecture,
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

        # Collision-induced absorption (HITRAN .cia files), if any.
        # CIA is treated as a fixed contribution — no Jacobian (τ̇_abs unchanged).
        for cia_path in abs_params.cia_files
            @timeit "CIA $(basename(cia_path))" begin
                cia_table = Absorption.load_cia_table(cia_path,
                                                     params.spec_bands[i_band];
                                                     FT = FT)
                Absorption.compute_τ_cia!(τ_abs[i_band], cia_table, profile,
                                           abs_params.vmr)
            end
        end

        # MT_CKD H₂O continuum, if a reference table is configured.
        # Treated as fixed (no Jacobian wrt H₂O VMR for now).
        if !isempty(abs_params.mtckd_file)
            @timeit "MT_CKD H2O continuum" begin
                mtckd_table = Absorption.load_mtckd(abs_params.mtckd_file)
                Absorption.compute_τ_h2o_continuum!(τ_abs[i_band], mtckd_table,
                                                     params.spec_bands[i_band], profile,
                                                     profile.vmr_h2o)
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

        _mie(λ) = make_mie_model(scat.decomp_type, mie_aerosol, λ,
            params.polarization_type, truncation_type, scat.r_max, scat.nquad_radius)

        mie_model = _mie(scat.λ_ref)
        k_ref, k̇_ref   = compute_ref_aerosol_extinction(lin, mie_model, params.float_type)

        for i_band=1:n_bands

            curr_band_λ = params.float_type(1e4) ./ params.spec_bands[i_band]

            if length(curr_band_λ)==1
                mie_model = _mie((maximum(curr_band_λ)+minimum(curr_band_λ))/2)
                @timeit "Mie calc"  aerosol_optics_raw, lin_aerosol_optics_raw =
                                compute_aerosol_optical_properties(lin, mie_model, FT2)

            else
                # Multi-spectral band: compute Mie properties at band edges (λ[1], λ[end]) and
                # linearly interpolate k_ext, k_sca and their derivatives in wavenumber across
                # the band. Greek coefficients are averaged. This avoids a full Mie calculation
                # at every spectral point (Mie varies smoothly with λ).
                n_spec = length(curr_band_λ)
                mie_model_0 = _mie(curr_band_λ[1])
                @timeit "Mie calc"  aerosol_optics_raw_0, lin_aerosol_optics_raw_0 =
                                compute_aerosol_optical_properties(lin, mie_model_0, FT2)

                mie_model_1 = _mie(curr_band_λ[end])
                @timeit "Mie calc"  aerosol_optics_raw_1, lin_aerosol_optics_raw_1 =
                                compute_aerosol_optical_properties(lin, mie_model_1, FT2)

                Nl  = length(aerosol_optics_raw_1.greek_coefs.α)
                Nl_ = length(aerosol_optics_raw_0.greek_coefs.α)
                gc0 = aerosol_optics_raw_0.greek_coefs
                gc1 = aerosol_optics_raw_1.greek_coefs
                lgc0 = lin_aerosol_optics_raw_0.lin_greek_coefs
                lgc1 = lin_aerosol_optics_raw_1.lin_greek_coefs

                greek_arrs = Dict{Symbol,Vector{Float64}}()
                lin_greek_arrs = Dict{Symbol,Matrix{Float64}}()
                for (fn, lfn) in zip((:α,:β,:γ,:δ,:ϵ,:ζ), (:α̇,:β̇,:γ̇,:δ̇,:ϵ̇,:ζ̇))
                    g0 = getfield(gc0, fn); g1 = getfield(gc1, fn)
                    arr = zeros(Nl)
                    arr[1:Nl_]     .= 0.5 .* (g0[1:Nl_] .+ g1[1:Nl_])
                    arr[1+Nl_:Nl]  .= 0.5 .* g1[1+Nl_:Nl]
                    greek_arrs[fn] = arr

                    lg0 = getfield(lgc0, lfn); lg1 = getfield(lgc1, lfn)
                    darr = zeros(4, Nl)
                    darr[:, 1:Nl_]    .= 0.5 .* (lg0[:, 1:Nl_] .+ lg1[:, 1:Nl_])
                    darr[:, 1+Nl_:Nl] .= 0.5 .* lg1[:, 1+Nl_:Nl]
                    lin_greek_arrs[lfn] = darr
                end

                greek_coeffs = GreekCoefs(greek_arrs[:α], greek_arrs[:β], greek_arrs[:γ],
                                          greek_arrs[:δ], greek_arrs[:ϵ], greek_arrs[:ζ])
                lin_greek_coeffs = linGreekCoefs(lin_greek_arrs[:α̇], lin_greek_arrs[:β̇], lin_greek_arrs[:γ̇],
                                  lin_greek_arrs[:δ̇], lin_greek_arrs[:ϵ̇], lin_greek_arrs[:ζ̇])
                ν_grid = [1e4/curr_band_λ[1], 1e4/curr_band_λ[end]]
                kext_grid = [aerosol_optics_raw_0.k, aerosol_optics_raw_1.k]
                ksca_grid = [aerosol_optics_raw_0.k*aerosol_optics_raw_0.ω̃, aerosol_optics_raw_1.k*aerosol_optics_raw_1.ω̃]
                interp_linear_kext = LinearInterpolation(ν_grid, kext_grid)
                interp_linear_ksca = LinearInterpolation(ν_grid, ksca_grid)
                k = zeros(n_spec)
                ω̃ = zeros(n_spec)
                fᵗ = zeros(n_spec)
                for i = 1:n_spec
                    k[i] = interp_linear_kext(1e4/curr_band_λ[i])
                    ω̃[i] = interp_linear_ksca(1e4/curr_band_λ[i])/k[i]
                end
                k̇ = zeros(4, n_spec)
                ω̃̇ = zeros(4, n_spec)
                ḟᵗ= zeros(4, n_spec)
                for ctr=1:4
                    k̇ext_grid = [lin_aerosol_optics_raw_0.k̇[ctr], lin_aerosol_optics_raw_1.k̇[ctr]]
                    # Interpolate dk_sca/dp using product rule: d(k·ω̃)/dp = dk/dp·ω̃ + k·dω̃/dp
                    k̇sca_grid = [lin_aerosol_optics_raw_0.k̇[ctr] * aerosol_optics_raw_0.ω̃ +
                                  aerosol_optics_raw_0.k * lin_aerosol_optics_raw_0.ω̃̇[ctr],
                                  lin_aerosol_optics_raw_1.k̇[ctr] * aerosol_optics_raw_1.ω̃ +
                                  aerosol_optics_raw_1.k * lin_aerosol_optics_raw_1.ω̃̇[ctr]]
                    interp_linear_k̇ext = LinearInterpolation(ν_grid, k̇ext_grid)
                    interp_linear_k̇sca = LinearInterpolation(ν_grid, k̇sca_grid)

                    for i = 1:n_spec
                        k̇[ctr,i] = interp_linear_k̇ext(1e4/curr_band_λ[i])
                        # Recover dω̃/dp via quotient rule: d(k_sca/k_ext)/dp = (dk_sca/dp − ω̃·dk_ext/dp)/k_ext
                        ω̃̇[ctr,i] = (interp_linear_k̇sca(1e4/curr_band_λ[i]) - ω̃[i] * k̇[ctr,i]) / k[i]
                    end
                end

                aerosol_optics_raw = AerosolOptics(greek_coefs=greek_coeffs, ω̃=ω̃, k=k, fᵗ=fᵗ)
                lin_aerosol_optics_raw = linAerosolOptics(lin_greek_coefs=lin_greek_coeffs, ω̃̇=ω̃̇, k̇=k̇, ḟᵗ=ḟᵗ)
            end

            if length(aerosol_optics_raw.greek_coefs.β) > truncation_type.l_max
                aerosol_optics[i_band][i_aer], lin_aerosol_optics[i_band][i_aer] =
                    Scattering.truncate_phase(truncation_type,
                                aerosol_optics_raw, lin_aerosol_optics_raw; reportFit=false)
                l_max_aer[i_aer, i_band] = truncation_type.l_max
            else
                aerosol_optics[i_band][i_aer] = aerosol_optics_raw
                lin_aerosol_optics[i_band][i_aer] = lin_aerosol_optics_raw
                l_max_aer[i_aer, i_band] = length(aerosol_optics_raw.greek_coefs.β)
            end

            k_band = aerosol_optics[i_band][i_aer].k
            k̇_band = lin_aerosol_optics[i_band][i_aer].k̇

            aer_p₀ = mean(c_aero.profile)
            aer_σp = std(c_aero.profile)
            τₚ, dτₚdp₀, dτₚdσp = getAerosolLayerOptProp(lin, 1, aer_p₀, aer_σp, profile.p_half)

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
            τ_aer[i_band][i_aer,:,:] =
                (τ_ref/k_ref) * k_band * τₚ'

            τ̇_aer[i_band][i_aer,1,:,:] .=
                (k_band/k_ref) * τₚ'

            for ctr=1:4
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    ((τ_ref/k_ref) * k̇_band[ctr,:] .-
                    (τ_ref/k_ref^2) * k̇_ref[ctr] .* k_band)*τₚ'
            end
            for ctr=5:6
                τ̇_aer[i_band][i_aer,ctr+1,:,:] =
                    (τ_ref/k_ref) * k_band *
                    (ctr==5 ? dτₚdp₀' : dτₚdσp')
            end

        end
    end
    for i_band = 1:n_bands
        l_max[i_band] = maximum(l_max_aer[:,i_band])
        max_m[i_band] = Int(ceil(l_max[i_band] + 1)/2)
    end
    set_uniform_lmax!(l_max, aerosol_optics)

    # Build the hierarchical RTModel
    solver = SolverConfig{FT, typeof(params.polarization_type), typeof(params.quadrature_type)}(
        params.polarization_type,
        params.quadrature_type,
        params.max_m,
        max_m,
        l_max,
        params.l_trunc,
        FT(params.Δ_angle),
        FT(params.depol),
    )
    atm = Atmosphere(profile, params.spec_bands)
    rayleigh_s = RayleighScattering(greek_rayleigh, greek_cabannes, FT.(ϖ_Cabannes))
    aerosols_s = AerosolState(aerosol_optics, τ_aer)
    optics = Optics(rayleigh_s, aerosols_s, τ_abs, τ_rayl)
    model = RTModel(params.architecture, solver, obs_geom, quad_points, atm, optics, params.brdf)
    return model, RTModelLin(τ̇_abs, τ̇_aer, lin_aerosol_optics)
end

"""
    model_from_parameters_lin(params)

Convenience alias for `model_from_parameters(LinMode(), params)`.
Returns `(model, lin_model)`.
"""
model_from_parameters_lin(params::vSmartMOM_Parameters) =
    model_from_parameters(LinMode(), params)
