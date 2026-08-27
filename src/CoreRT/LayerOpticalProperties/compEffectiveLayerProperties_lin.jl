#=
AD Architecture Boundary
========================
This file is the interface between the "physical parameter → optical property" chain rule
and the "optical property → TOA radiance" adding-doubling propagation.

Physical parameters (VMR, aerosol τ_ref/n_r/n_i/size/profile, surface albedo) produce
derivatives ∂(τ,ϖ,Z)/∂x_phys via Mie theory (ForwardDiff), Beer-Lambert law, and
mixing rules (quotient rule). These derivatives are computed HERE and stored in
`CoreScatteringOpticalPropertiesLin`.

The RT kernel (elemental!, doubling!, interaction!) then only needs to propagate
w.r.t. the 3 core variables (τ, ϖ, Z) using the chain rule:
    ∂R/∂x = ∂R/∂τ · ∂τ/∂x + ∂R/∂ϖ · ∂ϖ/∂x + ∂R/∂Z · ∂Z/∂x
This clean separation makes it straightforward to add new physical parameters
without touching the RT propagation code.
=#

function _compute_phase_blocks_lin(model, greek, lin_greek, m, arr_type)
    q = model.quad_points
    pol = CoreRT.polarization_type(model)
    if q.external_solar
        μ = collect(q.qp_μ)
        Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺ = Scattering.compute_Z_moments(
            pol, μ, greek, lin_greek, m; arr_type=arr_type)
        Z₀⁺, Z₀⁻, Ż₀⁺, Ż₀⁻ = Scattering.compute_Z_source_moments(
            pol, μ, q.μ₀, greek, lin_greek, m; arr_type=arr_type)
        return Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺, Z₀⁺, Z₀⁻, Ż₀⁺, Ż₀⁻
    end
    Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺ = Scattering.compute_Z_moments(
        pol, collect(q.qp_μ), greek, lin_greek, m; arr_type=arr_type)
    return Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺, nothing, nothing, nothing, nothing
end

"Evaluate forward/tangent aerosol phase blocks at 2–3 retained nodes and interpolate."
function _compute_aerosol_phase_blocks_lin(model, optics, lin_optics, ν_spec,
                                            m, arr_type)
    optics.phase_ν === nothing && return _compute_phase_blocks_lin(
        model, optics.greek_coefs, lin_optics.lin_greek_coefs, m, arr_type)
    length(optics.phase_ν) == length(optics.phase_greek) ==
        length(lin_optics.phase_lin_greek) || throw(DimensionMismatch(
        "aerosol phase-node values and tangents must have matching lengths"))
    blocks = [_compute_phase_blocks_lin(model, g, lg, m, Array) for
              (g, lg) in zip(optics.phase_greek, lin_optics.phase_lin_greek)]
    interp(k) = arr_type(_interpolate_phase_nodes(ν_spec, optics.phase_ν,
                                                   [b[k] for b in blocks]))
    first_four = ntuple(interp, 4)
    if blocks[1][5] === nothing
        return first_four..., nothing, nothing, nothing, nothing
    end
    return first_four..., interp(5), interp(6), interp(7), interp(8)
end

"Moment-independent, δ-M-scaled aerosol quantities for one layer."
struct LinAerosolInvariant{T1,T2}
    τ::T1
    ϖ::T1
    τ̇::T2
    ϖ̇::T2
end

"""
    NativeLayerSelection

Internal component-local decomposition of an active native atmospheric
Jacobian basis. `include_pressure` controls the shared leading pressure
column; `aerosol_columns[iaer]` contains indices in that aerosol's historical
seven-column block; `gas_columns` contains indices in the flattened native gas
block. Splitting the selection before mixing prevents inactive columns from
entering the combined phase-Jacobian tensors.
"""
struct NativeLayerSelection
    include_pressure::Bool
    aerosol_columns::Vector{Vector{Int}}
    gas_columns::Vector{Int}
end

function _native_layer_selection(layout::ActiveParameterLayout,
                                 n_aerosol::Int, n_gas::Int)
    columns = native_layer_columns(layout)
    selected_pressure = 1 in columns
    aerosol_columns = Vector{Vector{Int}}(undef, n_aerosol)
    reconstructed = Int[]
    selected_pressure && push!(reconstructed, 1)
    for iaer in 1:n_aerosol
        native = (2 + 7 * (iaer - 1)):(1 + 7 * iaer)
        local_columns = [column - first(native) + 1
                         for column in columns if column in native]
        aerosol_columns[iaer] = local_columns
        append!(reconstructed, first(native) .+ local_columns .- 1)
    end
    gas_native = (2 + 7 * n_aerosol):(1 + 7 * n_aerosol + n_gas)
    gas_columns = [column - first(gas_native) + 1
                   for column in columns if column in gas_native]
    append!(reconstructed, first(gas_native) .+ gas_columns .- 1)
    reconstructed == columns || throw(ArgumentError(
        "active native layer columns must follow pressure, aerosol-component, " *
        "then gas order; got $columns"))
    return NativeLayerSelection(selected_pressure, aerosol_columns, gas_columns)
end

"""
Moment-invariant cache used by the linearized Fourier loop. When `selection`
is non-`nothing`, aerosol and gas tangent arrays already carry only the active
component-local columns; all forward arrays remain complete.
"""
struct LinMInvariantCache
    rayl_τ_dev::Vector{Vector}
    aerosol::Vector{Vector{Vector}}
    gas::Vector{Vector}
    lin_gas::Vector{Vector}
    selection::Union{Nothing,NativeLayerSelection}
end

_same_greek(a, b) = all(name -> getfield(a, name) == getfield(b, name),
                         fieldnames(typeof(a)))

function _createAero_invariant(τAer, aerosol_optics, τ̇Aer,
                               lin_aerosol_optics, arr_type,
                               columns::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    (; fᵗ, ω̃) = aerosol_optics
    (; ḟᵗ, ω̃̇) = lin_aerosol_optics
    n = size(τAer, 1)
    τ̇Aer = _to_device(arr_type, collect(τ̇Aer'))
    ω̃ = ω̃ isa Number ? arr_type(fill(ω̃, n)) : _to_device(arr_type, ω̃)
    fᵗ = fᵗ isa Number ? arr_type(fill(fᵗ, n)) : _to_device(arr_type, fᵗ)
    ω̃̇_block = _lift_mie_param_to_n_x_4(ω̃̇, n, arr_type)
    ḟᵗ_block = _lift_mie_param_to_n_x_4(ḟᵗ, n, arr_type)
    fω = fᵗ .* ω̃
    τ_mod = (1 .- fω) .* τAer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fω)
    τ̇_mod = arr_type(zeros(eltype(τAer), n, 7))
    ϖ̇_mod = arr_type(zeros(eltype(τAer), n, 7))
    τ̇_mod[:,1] .= (1 .- fω) .* τ̇Aer[:,1]
    tmp = fᵗ .* ω̃̇_block .+ ω̃ .* ḟᵗ_block
    τ̇_mod[:,2:5] .= (1 .- fω) .* τ̇Aer[:,2:5] .- tmp .* τAer
    ϖ̇_mod[:,2:5] .= (ω̃̇_block .* (1 .- fᵗ) .-
        ḟᵗ_block .* (ω̃ .* (1 .- ω̃))) ./ (1 .- fω).^2
    τ̇_mod[:,6:7] .= (1 .- fω) .* τ̇Aer[:,6:7]
    if columns !== nothing
        τ̇_mod = τ̇_mod[:, columns]
        ϖ̇_mod = ϖ̇_mod[:, columns]
    end
    return LinAerosolInvariant(τ_mod, ϖ_mod, τ̇_mod, ϖ̇_mod)
end

"""
    build_m_invariant_cache_lin(iBand, model, lin_model;
                                active_layout=nothing)

Construct Fourier-independent Rayleigh, aerosol, and gas optical properties.
With an active layout, split its native atmospheric columns by component and
cache only those aerosol/gas tangents. This is the earliest safe selection
point: the forward mixing weights are already known, but no Fourier-dependent
phase tensor or MOM operator workspace has been allocated.
"""
function build_m_invariant_cache_lin(iBand, model, lin_model;
                                     active_layout::Union{Nothing,ActiveParameterLayout}=nothing)
    bands = iBand isa Integer ? (iBand,) : iBand
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics) = model
    (; τ̇_aer, τ̇_abs, lin_aerosol_optics) = lin_model
    arr_type = CoreRT.array_type(model)
    nZ = size(τ_rayl[1], 2)
    nAero = size(τ_aer[first(bands)], 1)
    nGas = size(τ̇_abs[first(bands)], 1)
    selection = active_layout === nothing ? nothing :
        _native_layer_selection(active_layout, nAero, nGas)
    rayl = Vector{Vector}(undef, length(bands))
    aeros = Vector{Vector{Vector}}(undef, length(bands))
    gas = Vector{Vector}(undef, length(bands))
    lin_gas = Vector{Vector}(undef, length(bands))
    for (iBi, iB) in enumerate(bands)
        rayl[iBi] = [_to_device(arr_type, τ_rayl[iB][:,iz]) for iz in 1:nZ]
        aeros[iBi] = [[_createAero_invariant(
            _to_device(arr_type, τ_aer[iB][iaer,:,iz]), aerosol_optics[iB][iaer],
            τ̇_aer[iB][iaer,:,:,iz], lin_aerosol_optics[iB][iaer], arr_type,
            selection === nothing ? nothing : selection.aerosol_columns[iaer])
            for iz in 1:nZ] for iaer in 1:nAero]
        gas[iBi] = [CoreAbsorptionOpticalProperties(
            _to_device(arr_type, τ_abs[iB][:,iz])) for iz in 1:nZ]
        gas_columns = selection === nothing ? axes(τ̇_abs[iB], 1) :
                      selection.gas_columns
        lin_gas[iBi] = [CoreAbsorptionOpticalPropertiesLin(
            _to_device(arr_type, collect(τ̇_abs[iB][gas_columns,:,iz]'))) for iz in 1:nZ]
    end
    return LinMInvariantCache(rayl, aeros, gas, lin_gas, selection)
end

function _attach_aerosol_phase(inv::LinAerosolInvariant,
        Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺, Z₀⁺, Z₀⁻, Ż₀⁺, Ż₀⁻, arr_type,
        columns::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    n = length(inv.τ)
    active_columns = columns === nothing ? collect(1:7) : collect(columns)
    function lift_phase_dot(Z, nrow, ncol)
        FT = Z === nothing ? eltype(inv.τ) : eltype(Z)
        if Z === nothing || ndims(Z) == 3
            out = arr_type(zeros(FT, nrow, ncol, length(active_columns)))
            if Z !== nothing
                for (iout, column) in enumerate(active_columns)
                    2 <= column <= 5 || continue
                    out[:,:,iout] .= @view Z[column - 1, :, :]
                end
            end
        else
            out = arr_type(zeros(FT, nrow, ncol, n, length(active_columns)))
            for (iout, column) in enumerate(active_columns)
                2 <= column <= 5 || continue
                out[:,:,:,iout] .= @view Z[column - 1, :, :, :]
            end
        end
        return out
    end
    Ż7⁺⁺ = lift_phase_dot(Ż⁺⁺, size(Z⁺⁺,1), size(Z⁺⁺,2))
    Ż7⁻⁺ = lift_phase_dot(Ż⁻⁺, size(Z⁻⁺,1), size(Z⁻⁺,2))
    Ż7₀⁺ = Z₀⁺ === nothing ? nothing : lift_phase_dot(Ż₀⁺, size(Z₀⁺,1), size(Z₀⁺,2))
    Ż7₀⁻ = Z₀⁻ === nothing ? nothing : lift_phase_dot(Ż₀⁻, size(Z₀⁻,1), size(Z₀⁻,2))
    return CoreScatteringOpticalProperties(inv.τ, inv.ϖ, Z⁺⁺, Z⁻⁺, Z₀⁺, Z₀⁻),
        CoreScatteringOpticalPropertiesLin(inv.τ̇, inv.ϖ̇,
                                            Ż7⁺⁺, Ż7⁻⁺, Ż7₀⁺, Ż7₀⁻)
end

"""
    constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model)

Construct the combined effective layer optical properties and their linearized derivatives
for all atmospheric layers, merging Rayleigh scattering, aerosol scattering (with δ-M
truncation), and trace gas absorption.

For each atmospheric layer, this function:
1. Creates Rayleigh scattering properties (τ_rayl, ϖ=1, Z_rayl).
2. For each aerosol type, computes δ-M scaled optical properties via [`createAero`](@ref),
   combining Mie scattering with the aerosol vertical profile.
3. Adds Rayleigh and aerosol properties using the `+` operator on
   `UmbrellaCoreScatteringOpticalProperties`, which correctly propagates derivatives through
   the mixing formulas for combined ϖ and Z.
4. Adds gas absorption via `UmbrellaCoreAbsorptionOpticalProperties`, updating ϖ and τ.

The output derivative dimension has `Nparams = 7×NAer + NGasLayer` entries per layer,
where `NGasLayer = NGasSpecies×Nz`.
(surface parameters are handled separately in the interaction step).

# Returns
- `layer_opt`: Array of `CoreScatteringOpticalProperties` (one per layer).
- `layer_opt_lin`: Array of `CoreScatteringOpticalPropertiesLin` (derivatives, one per layer).
- `fscat_opt`: Rayleigh scattering fraction per layer (for inelastic scattering weight).
"""
function constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model) #where {FT<:Real}
    if InelasticScattering.has_inelastic(RS_type)
        throw(ArgumentError(
            "Linearized Raman-active optical properties are intentionally unsupported. " *
            "Use forward Raman RT only."))
    end

    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, 
            greek_rayleigh) = model
    (; τ̇_aer, τ̇_abs, lin_aerosol_optics,
       τ̇_rayl_psurf, τ̇_aer_psurf, τ̇_abs_psurf) = lin_model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"
    FT = eltype(τ_rayl[1])
    
    arr_type = CoreRT.array_type(model)

    pol_type = CoreRT.polarization_type(model)
    # Do this in CPU space only first:
    
    # Quadrature points:
    μ = collect(model.quad_points.qp_μ )
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    # Rayleigh Z matrix. Linearized RT is pure-elastic only; Raman-active
    # Cabannes/RRS modes are rejected above.
    
    band_layer_props     = [];
    band_layer_props_lin = [];
    band_fScattRayleigh  = [];
    for iB in iBand
        # Linearized RT is pure-elastic only — Raman/Cabannes types are rejected
        # at the lin entry points — so both noRS and noRS_plus use the full
        # Rayleigh phase matrix. greek_cabannes is intentionally NOT destructured
        # for this path; a Cabannes branch here would be dead code (and an
        # UndefVarError for noRS_plus). See project_raman_not_linearized.
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻ =
            _compute_phase_blocks(model, greek_rayleigh[iB], m, arr_type)

        rayl = [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), 1.0,
                Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻) for i=1:nZ]
        # Initiate combined properties with rayleigh
        combrella = [UmbrellaCoreScatteringOpticalProperties(rayl[i],nothing) for i=1:nZ]
        aer_ps_components = [Tuple[] for _ in 1:nZ]
        # Loop over all aerosol types:
        for iaer=1:nAero
            # Precomute Z matrices per type (constant per layer)
            AerZ⁺⁺, AerZ⁻⁺, AerŻ⁺⁺, AerŻ⁻⁺,
            AerZ₀⁺, AerZ₀⁻, AerŻ₀⁺, AerŻ₀⁻ = _compute_aerosol_phase_blocks_lin(
                model, aerosol_optics[iB][iaer], lin_aerosol_optics[iB][iaer],
                get_spec_bands(model)[iB], m, arr_type)
            # Generate Core optical properties for Aerosols iaer
            aer = []
            lin_aer = []
            for iz=1:nZ
                t_aer, t_lin_aer =  createAero(
                            arr_type(τ_aer[iB][iaer,:,iz]), 
                            aerosol_optics[iB][iaer], 
                            arr_type(AerZ⁺⁺), arr_type(AerZ⁻⁺), 
                            arr_type(τ̇_aer[iB][iaer,:,:,iz]), 
                            lin_aerosol_optics[iB][iaer], 
                            arr_type(AerŻ⁺⁺), arr_type(AerŻ⁻⁺),
                            AerZ₀⁺, AerZ₀⁻, AerŻ₀⁺, AerŻ₀⁻,
                            arr_type) 
                push!(aer, t_aer)
                push!(lin_aer, t_lin_aer)
                rawdot = arr_type(τ̇_aer_psurf[iB][iaer, :, iz])
                push!(aer_ps_components[iz], (t_aer, rawdot))
            end 
            aer_combrella = [UmbrellaCoreScatteringOpticalProperties(aer[i],lin_aer[i]) for i=1:nZ]
            # Mix with previous Core Optical Properties
            combrella = [combrella[i]+aer_combrella[i] for i=1:nZ]
        end


        # Create Core Optical Properties merged with trace gas absorptions:
        gas = [CoreAbsorptionOpticalProperties(arr_type(τ_abs[iB][:,iz])) for iz=1:nZ]
        lin_gas = [CoreAbsorptionOpticalPropertiesLin(arr_type(collect(τ̇_abs[iB][:,:,iz]'))) for iz=1:nZ] 
        gas_combrella = [UmbrellaCoreAbsorptionOpticalProperties(gas[iz],lin_gas[iz]) for iz=1:nZ]   
        tmp = combrella .+ gas_combrella            
        combrella = tmp
        fScattRayleigh =
            [_rayleigh_fraction_of_total_extinction(rayl[iz], combrella[iz].fwd)
             for iz=1:nZ]
        #for i=1:nZ
        #end
        #combo_lin = [include_rayl!(combo[iz], combo_lin[iz], rayl[iz], rayl_lin[iz]) for iz=1:nZ]
        #for iaer=1:nAero
        #    aer_lin =  [createAeroLin(arr_type(τ̇_aer[iB][iaer,1:7,:,iz]), 
        #                aerosol_optics_lin[iB][iaer], 
        #                AerŻ⁺⁺, AerŻ⁻⁺, arr_type) for iz=1:nZ]
        #    combo_lin = [include_aer!(iaer, combo[iz], combo_lin[iz], aer[iz], aer_lin[iz]) for i=1:nZ]
        #end
        # Use the following two lines for every gas to be included in the state vector
        #igas=1
        #combo2_lin = [include_gas(nAero, igas[iz], combo[iz], combo_lin[iz], gas_lin[iz]) for i=1:nZ]
        combo =     [combrella[iz].fwd for iz=1:nZ]
        # `Any` is intentional here: legacy gas-only mixing may allocate its Z
        # tangent on the host, while prepending p_surf below moves the completed
        # tangent to the selected device.
        combo_lin = Any[combrella[iz].lin for iz=1:nZ]

        # One shared p_surf column precedes all aerosol and gas columns. Its
        # component tangents are combined at the final mixed-optics boundary.
        for iz in 1:nZ
            fwd = combo[iz]
            nSpec = length(fwd.τ)
            raydot = arr_type(τ̇_rayl_psurf[iB][:, iz])
            absdot = arr_type(τ̇_abs_psurf[iB][:, iz])
            τdot = raydot .+ absdot
            Wdot = copy(raydot)  # Rayleigh has ϖ=1.
            Ndot⁺⁺ = reshape(raydot, 1, 1, nSpec) .* Rayl𝐙⁺⁺
            Ndot⁻⁺ = reshape(raydot, 1, 1, nSpec) .* Rayl𝐙⁻⁺
            Ndot₀⁺ = fwd.Z₀⁺ === nothing ? nothing :
                reshape(raydot, 1, 1, nSpec) .* _ensure_3d(RaylZ₀⁺)
            Ndot₀⁻ = fwd.Z₀⁻ === nothing ? nothing :
                reshape(raydot, 1, 1, nSpec) .* _ensure_3d(RaylZ₀⁻)
            for (iaer, (aer_fwd, rawdot)) in enumerate(aer_ps_components[iz])
                optics = aerosol_optics[iB][iaer]
                trunc_factor = arr_type(one(FT) .- optics.fᵗ .* optics.ω̃)
                moddot = trunc_factor .* rawdot
                τdot .+= moddot
                scatterdot = moddot .* aer_fwd.ϖ
                Wdot .+= scatterdot
                Ndot⁺⁺ .+= reshape(scatterdot, 1, 1, nSpec) .* aer_fwd.Z⁺⁺
                Ndot⁻⁺ .+= reshape(scatterdot, 1, 1, nSpec) .* aer_fwd.Z⁻⁺
                if Ndot₀⁺ !== nothing
                    Ndot₀⁺ .+= reshape(scatterdot, 1, 1, nSpec) .* aer_fwd.Z₀⁺
                    Ndot₀⁻ .+= reshape(scatterdot, 1, 1, nSpec) .* aer_fwd.Z₀⁻
                end
            end
            W = fwd.τ .* fwd.ϖ
            ϖdot = (Wdot .- fwd.ϖ .* τdot) ./ fwd.τ
            Zdot⁺⁺ = (Ndot⁺⁺ .- reshape(Wdot, 1, 1, nSpec) .* fwd.Z⁺⁺) ./
                      reshape(W, 1, 1, nSpec)
            Zdot⁻⁺ = (Ndot⁻⁺ .- reshape(Wdot, 1, 1, nSpec) .* fwd.Z⁻⁺) ./
                      reshape(W, 1, 1, nSpec)
            Zdot₀⁺ = Ndot₀⁺ === nothing ? nothing :
                (Ndot₀⁺ .- reshape(Wdot, 1, 1, nSpec) .* fwd.Z₀⁺) ./
                reshape(W, 1, 1, nSpec)
            Zdot₀⁻ = Ndot₀⁻ === nothing ? nothing :
                (Ndot₀⁻ .- reshape(Wdot, 1, 1, nSpec) .* fwd.Z₀⁻) ./
                reshape(W, 1, 1, nSpec)
            old = combo_lin[iz]
            combo_lin[iz] = CoreScatteringOpticalPropertiesLin(
                hcat(τdot, old.τ̇), hcat(ϖdot, old.ϖ̇),
                cat(reshape(Zdot⁺⁺, size(Zdot⁺⁺)..., 1), arr_type(old.Ż⁺⁺); dims=4),
                cat(reshape(Zdot⁻⁺, size(Zdot⁻⁺)..., 1), arr_type(old.Ż⁻⁺); dims=4),
                Zdot₀⁺ === nothing ? nothing : cat(
                    reshape(Zdot₀⁺, size(Zdot₀⁺)..., 1), arr_type(old.Ż₀⁺); dims=4),
                Zdot₀⁻ === nothing ? nothing : cat(
                    reshape(Zdot₀⁻, size(Zdot₀⁻)..., 1), arr_type(old.Ż₀⁻); dims=4))
        end

        push!(band_layer_props, combo)
        push!(band_layer_props_lin, combo_lin)
        push!(band_fScattRayleigh,fScattRayleigh)
    end
    layer_opt     = [prod([band_layer_props[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    layer_opt_lin = [prod([band_layer_props_lin[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    fscat_opt     = [[band_fScattRayleigh[i][iz] for i=1:length(iBand)] for iz=1:nZ]
    return layer_opt, layer_opt_lin, fscat_opt
end

"Cache-aware linearized construction: only phase blocks and phase mixing vary with m."
function constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model,
                                        cache::LinMInvariantCache)
    bands = iBand isa Integer ? (iBand,) : iBand
    InelasticScattering.has_inelastic(RS_type) && throw(ArgumentError(
        "Linearized Raman-active optical properties are intentionally unsupported."))
    (; τ_rayl, τ_aer, aerosol_optics, greek_rayleigh) = model
    (; lin_aerosol_optics, τ̇_rayl_psurf, τ̇_aer_psurf,
       τ̇_abs_psurf) = lin_model
    arr_type = CoreRT.array_type(model)
    nAero = size(τ_aer[first(bands)], 1)
    nZ = size(τ_rayl[1], 2)
    FT = eltype(τ_rayl[1])
    selected_rayleigh = [greek_rayleigh[iB] for iB in bands]
    shared_rayleigh = length(selected_rayleigh) > 1 &&
                      all(g -> _same_greek(g, selected_rayleigh[1]), selected_rayleigh)
    shared_blocks = shared_rayleigh ?
        _compute_phase_blocks(model, selected_rayleigh[1], m, arr_type) : nothing

    band_layer_props = Any[]
    band_layer_props_lin = Any[]
    band_fScattRayleigh = Any[]
    for (iBi, iB) in enumerate(bands)
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻ = shared_blocks === nothing ?
            _compute_phase_blocks(model, greek_rayleigh[iB], m, arr_type) : shared_blocks
        rayl = [CoreScatteringOpticalProperties(cache.rayl_τ_dev[iBi][iz], 1.0,
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻) for iz in 1:nZ]
        combrella = [UmbrellaCoreScatteringOpticalProperties(rayl[iz], nothing)
                     for iz in 1:nZ]
        aer_ps_components = [Tuple[] for _ in 1:nZ]
        for iaer in 1:nAero
            aerosol_columns = cache.selection === nothing ? nothing :
                              cache.selection.aerosol_columns[iaer]
            needs_mie_phase = aerosol_columns === nothing ||
                              any(column -> 2 <= column <= 5, aerosol_columns)
            if needs_mie_phase
                AerZ⁺⁺, AerZ⁻⁺, AerŻ⁺⁺, AerŻ⁻⁺,
                AerZ₀⁺, AerZ₀⁻, AerŻ₀⁺, AerŻ₀⁻ =
                    _compute_aerosol_phase_blocks_lin(
                        model, aerosol_optics[iB][iaer], lin_aerosol_optics[iB][iaer],
                        get_spec_bands(model)[iB], m, arr_type)
            else
                AerZ⁺⁺, AerZ⁻⁺, AerZ₀⁺, AerZ₀⁻ =
                    _compute_aerosol_phase_blocks(
                        model, aerosol_optics[iB][iaer],
                        get_spec_bands(model)[iB], m, arr_type)
                AerŻ⁺⁺ = AerŻ⁻⁺ = AerŻ₀⁺ = AerŻ₀⁻ = nothing
            end
            aer = Any[]
            lin_aer = Any[]
            for iz in 1:nZ
                afwd, alin = _attach_aerosol_phase(cache.aerosol[iBi][iaer][iz],
                    AerZ⁺⁺, AerZ⁻⁺, AerŻ⁺⁺, AerŻ⁻⁺,
                    AerZ₀⁺, AerZ₀⁻, AerŻ₀⁺, AerŻ₀⁻, arr_type,
                    aerosol_columns)
                push!(aer, afwd); push!(lin_aer, alin)
                push!(aer_ps_components[iz],
                      (afwd, _to_device(arr_type, τ̇_aer_psurf[iB][iaer,:,iz])))
            end
            combrella = [combrella[iz] +
                UmbrellaCoreScatteringOpticalProperties(aer[iz], lin_aer[iz])
                for iz in 1:nZ]
        end
        combrella = [combrella[iz] + UmbrellaCoreAbsorptionOpticalProperties(
            cache.gas[iBi][iz], cache.lin_gas[iBi][iz]) for iz in 1:nZ]
        fScattRayleigh = [_rayleigh_fraction_of_total_extinction(
            rayl[iz], combrella[iz].fwd) for iz in 1:nZ]
        combo = [combrella[iz].fwd for iz in 1:nZ]
        combo_lin = Any[combrella[iz].lin for iz in 1:nZ]

        # Prepend the shared surface-pressure column when selected. Its phase
        # numerators are m-dependent; τ/ϖ component tangents are cached above.
        include_pressure = cache.selection === nothing ||
                           cache.selection.include_pressure
        for iz in (include_pressure ? (1:nZ) : (1:0))
            fwd = combo[iz]; nSpec = length(fwd.τ)
            raydot = _to_device(arr_type, τ̇_rayl_psurf[iB][:,iz])
            absdot = _to_device(arr_type, τ̇_abs_psurf[iB][:,iz])
            τdot = raydot .+ absdot
            Wdot = copy(raydot)
            Ndot⁺⁺ = reshape(raydot,1,1,nSpec) .* Rayl𝐙⁺⁺
            Ndot⁻⁺ = reshape(raydot,1,1,nSpec) .* Rayl𝐙⁻⁺
            Ndot₀⁺ = fwd.Z₀⁺ === nothing ? nothing :
                reshape(raydot,1,1,nSpec) .* _ensure_3d(RaylZ₀⁺)
            Ndot₀⁻ = fwd.Z₀⁻ === nothing ? nothing :
                reshape(raydot,1,1,nSpec) .* _ensure_3d(RaylZ₀⁻)
            for (iaer, (aer_fwd, rawdot)) in enumerate(aer_ps_components[iz])
                optics = aerosol_optics[iB][iaer]
                raw_factor = one(FT) .- optics.fᵗ .* optics.ω̃
                trunc_factor = raw_factor isa AbstractArray ?
                    _to_device(arr_type, raw_factor) : raw_factor
                moddot = trunc_factor .* rawdot
                τdot .+= moddot
                scatterdot = moddot .* aer_fwd.ϖ
                Wdot .+= scatterdot
                Ndot⁺⁺ .+= reshape(scatterdot,1,1,nSpec) .* aer_fwd.Z⁺⁺
                Ndot⁻⁺ .+= reshape(scatterdot,1,1,nSpec) .* aer_fwd.Z⁻⁺
                if Ndot₀⁺ !== nothing
                    Ndot₀⁺ .+= reshape(scatterdot,1,1,nSpec) .* aer_fwd.Z₀⁺
                    Ndot₀⁻ .+= reshape(scatterdot,1,1,nSpec) .* aer_fwd.Z₀⁻
                end
            end
            W = fwd.τ .* fwd.ϖ
            ϖdot = (Wdot .- fwd.ϖ .* τdot) ./ fwd.τ
            zdot(Ndot, Z) = (Ndot .- reshape(Wdot,1,1,nSpec) .* Z) ./
                             reshape(W,1,1,nSpec)
            Zdot⁺⁺, Zdot⁻⁺ = zdot(Ndot⁺⁺, fwd.Z⁺⁺), zdot(Ndot⁻⁺, fwd.Z⁻⁺)
            Zdot₀⁺ = Ndot₀⁺ === nothing ? nothing : zdot(Ndot₀⁺, fwd.Z₀⁺)
            Zdot₀⁻ = Ndot₀⁻ === nothing ? nothing : zdot(Ndot₀⁻, fwd.Z₀⁻)
            old = combo_lin[iz]
            combo_lin[iz] = CoreScatteringOpticalPropertiesLin(
                hcat(τdot,old.τ̇), hcat(ϖdot,old.ϖ̇),
                cat(reshape(Zdot⁺⁺,size(Zdot⁺⁺)...,1),old.Ż⁺⁺;dims=4),
                cat(reshape(Zdot⁻⁺,size(Zdot⁻⁺)...,1),old.Ż⁻⁺;dims=4),
                Zdot₀⁺ === nothing ? nothing : cat(reshape(Zdot₀⁺,size(Zdot₀⁺)...,1),old.Ż₀⁺;dims=4),
                Zdot₀⁻ === nothing ? nothing : cat(reshape(Zdot₀⁻,size(Zdot₀⁻)...,1),old.Ż₀⁻;dims=4))
        end
        push!(band_layer_props,combo); push!(band_layer_props_lin,combo_lin)
        push!(band_fScattRayleigh,fScattRayleigh)
    end
    layer_opt = [prod([band_layer_props[i][iz] for i in eachindex(bands)]) for iz in 1:nZ]
    layer_opt_lin = [prod([band_layer_props_lin[i][iz] for i in eachindex(bands)]) for iz in 1:nZ]
    fscat_opt = [[band_fScattRayleigh[i][iz] for i in eachindex(bands)] for iz in 1:nZ]
    return layer_opt, layer_opt_lin, fscat_opt
end
 

"""
    createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺, arr_type)

!!! note "Legacy function"
    Consider using [`delta_m_truncation_lin`](@ref) directly for new code.
    It provides the same δ-M chain rule in a cleaner, composable form
    with explicit `RawAerosolJacobian` as the AD boundary.

Compute **δ-M scaled** aerosol optical properties and their derivatives for one aerosol
type in one atmospheric layer.

Applies the δ-M truncation correction (Nakajima & Tanaka, 1988; Sanghavi & Stephens, 2013)
to the aerosol optical depth and single-scattering albedo, and computes derivatives with
respect to 7 aerosol sub-parameters: `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`.

# δ-M scaling formulas

```math
\\tau_\\text{mod} = (1 - f^t \\tilde{\\omega}) \\cdot \\tau_\\text{aer}
```
```math
\\varpi_\\text{mod} = \\frac{(1 - f^t) \\tilde{\\omega}}{1 - f^t \\tilde{\\omega}}
```

# Derivative chain rule

For each physical parameter ``p_j``:
```math
\\frac{\\partial \\tau_\\text{mod}}{\\partial p_j} = 
  (1 - f^t \\tilde{\\omega}) \\frac{\\partial \\tau_\\text{aer}}{\\partial p_j}
  - \\tau_\\text{aer} \\left(f^t \\frac{\\partial \\tilde{\\omega}}{\\partial p_j} + 
    \\tilde{\\omega} \\frac{\\partial f^t}{\\partial p_j}\\right)
```

For `τ_ref, p₀, σp`: only the `τ_aer` chain contributes (Mie properties are independent).

!!! note "Bug 18 fix"
    Corrected index mapping: `τ̇Aer[k,:]` now correctly maps to parameter `k`
    (was previously off-by-one, mixing τ_ref with Mie parameter derivatives).

# Returns
- `CoreScatteringOpticalProperties`: Forward δ-M scaled properties.
- `CoreScatteringOpticalPropertiesLin`: Linearized properties (7 sub-params).
"""
# ----------------------------------------------------------------------
# Helper: normalise the Mie-parameter linearisation arrays into a
# uniform `(n_spec, 4)` block for the v0.6+ Jacobian assembly path.
#
# `lin_aerosol_optics.ω̃̇` and `.ḟᵗ` carry derivatives w.r.t. the four Mie
# parameters (n_r, n_i, r_m, sigma_r). Their concrete shape depends on
# the Fourier truncation chosen at parameters_from_yaml time:
#   ()                 -- scalar (NoTruncation: f^t = 0 everywhere)
#   (nparams,)         -- spectrally flat (legacy delta-BGE pre-v0.7)
#   (nparams, 1)       -- degenerate matrix form
#   (nparams, n_spec)  -- fully resolved (v0.7 auto + Mie Dual chain)
# `createAero` then writes columns 2:5 of the 7-parameter layout
# [tau_ref, n_r, n_i, r_m, sigma_r, p_0, sigma_p] uniformly.
function _lift_mie_param_to_n_x_4(x, n, arr_type)
    eltyp = eltype(x)
    if x isa Number
        return arr_type(fill(eltyp(x), n, 4))
    elseif ndims(x) == 1
        # (nparams,) -> broadcast to (n, nparams), pad/truncate to 4 cols
        x4 = length(x) == 4 ? collect(x) : _pad4(collect(x), eltyp)
        return arr_type(repeat(reshape(x4, 1, 4), n, 1))
    elseif ndims(x) == 2
        sz = size(x)
        # Accept (nparams, 1) (degenerate) -- broadcast across spectral.
        if sz[2] == 1
            x4 = sz[1] == 4 ? collect(x[:, 1]) : _pad4(collect(x[:, 1]), eltyp)
            return arr_type(repeat(reshape(x4, 1, 4), n, 1))
        end
        # (nparams, n_spec) -- transpose to (n_spec, nparams), pad to 4 cols.
        if sz[1] in (1, 2, 3, 4) && sz[2] == n
            xt = collect(transpose(x))   # (n, nparams)
            return sz[1] == 4 ? arr_type(xt) : arr_type(_pad4_cols(xt, eltyp))
        end
        # (n_spec, nparams) already in the assembly orientation.
        if sz[1] == n && sz[2] in (1, 2, 3, 4)
            return sz[2] == 4 ? arr_type(collect(x)) : arr_type(_pad4_cols(collect(x), eltyp))
        end
    end
    error("createAero: unexpected shape ", size(x), " for Mie parameter linearisation array; expected scalar, (nparams,), or 2-D with one dim equal to n_spec=", n)
end

_pad4(v::AbstractVector, T) = (out = zeros(T, 4); out[1:length(v)] .= v; out)
_pad4_cols(M::AbstractMatrix, T) = (out = zeros(T, size(M, 1), 4); out[:, 1:size(M, 2)] .= M; out)

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, 
                    τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺,
                    AerZ₀⁺, AerZ₀⁻, AerŻ₀⁺, AerŻ₀⁻,
                    arr_type)

    (; fᵗ, ω̃) = aerosol_optics
    (; ḟᵗ, ω̃̇) = lin_aerosol_optics

    n  = size(τAer,1)
    τ̇Aer = arr_type(collect(τ̇Aer'))  # (7, n) → (n, 7) — transpose upstream input
    #fᵗ = arr_type(fᵗ)
    #ω̃  = arr_type(ω̃)
    #ḟᵗ = arr_type(ḟᵗ)
    #ω̃̇  = arr_type(ω̃̇)

    tmpŻ⁺⁺ = AerŻ⁺⁺
    tmpŻ⁻⁺ = AerŻ⁻⁺
    if ndims(tmpŻ⁺⁺) == 3
        AerŻ⁺⁺ = arr_type(zeros(size(AerZ⁺⁺,1), size(AerZ⁺⁺,2), 7))
        AerŻ⁻⁺ = arr_type(zeros(size(AerZ⁻⁺,1), size(AerZ⁻⁺,2), 7))
        AerŻ⁺⁺[:,:,2:5] .= permutedims(tmpŻ⁺⁺, (2, 3, 1))
        AerŻ⁻⁺[:,:,2:5] .= permutedims(tmpŻ⁻⁺, (2, 3, 1))
    elseif ndims(tmpŻ⁺⁺) == 4
        AerŻ⁺⁺ = arr_type(zeros(size(AerZ⁺⁺,1), size(AerZ⁺⁺,2), n, 7))
        AerŻ⁻⁺ = arr_type(zeros(size(AerZ⁻⁺,1), size(AerZ⁻⁺,2), n, 7))
        AerŻ⁺⁺[:,:,:,2:5] .= permutedims(tmpŻ⁺⁺, (2, 3, 4, 1))
        AerŻ⁻⁺[:,:,:,2:5] .= permutedims(tmpŻ⁻⁺, (2, 3, 4, 1))
    else
        error("createAero: expected phase derivatives (nParam,nμ,nμ) or (nParam,nμ,nμ,nSpec), got $(size(tmpŻ⁺⁺))")
    end
    if AerZ₀⁺ === nothing
        AerŻ₀⁺ = AerŻ₀⁻ = nothing
    elseif ndims(AerŻ₀⁺) == 3
        tmp⁺, tmp⁻ = AerŻ₀⁺, AerŻ₀⁻
        AerŻ₀⁺ = arr_type(zeros(size(AerZ₀⁺,1), size(AerZ₀⁺,2), 7))
        AerŻ₀⁻ = arr_type(zeros(size(AerZ₀⁻,1), size(AerZ₀⁻,2), 7))
        AerŻ₀⁺[:,:,2:5] .= permutedims(tmp⁺, (2, 3, 1))
        AerŻ₀⁻[:,:,2:5] .= permutedims(tmp⁻, (2, 3, 1))
    else
        tmp⁺, tmp⁻ = AerŻ₀⁺, AerŻ₀⁻
        AerŻ₀⁺ = arr_type(zeros(size(AerZ₀⁺,1), size(AerZ₀⁺,2), n, 7))
        AerŻ₀⁻ = arr_type(zeros(size(AerZ₀⁻,1), size(AerZ₀⁻,2), n, 7))
        AerŻ₀⁺[:,:,:,2:5] .= permutedims(tmp⁺, (2, 3, 4, 1))
        AerŻ₀⁻[:,:,:,2:5] .= permutedims(tmp⁻, (2, 3, 4, 1))
    end

    # Ensure arrays are in the right memory space (CPU or GPU)
    ω̃  = (ω̃ isa Number) ? arr_type(fill(ω̃,n)) : arr_type(ω̃)
    fᵗ  = (fᵗ isa Number) ? arr_type(fill(fᵗ,n)) : arr_type(fᵗ)
    
    ω̃̇ = arr_type(ω̃̇)
    ḟᵗ = arr_type(ḟᵗ)
    
    # ω̃̇ and ḟᵗ from Mie are derivatives w.r.t. the 4 Mie
    # parameters (n_r, n_i, r_m, sigma_r). They arrive in any of:
    #   ()                 -- scalar (e.g. NoTruncation case)
    #   (nparams,)         -- spectrally flat (legacy delta-BGE pre-v0.7)
    #   (nparams, 1)       -- degenerate matrix form
    #   (nparams, n_spec)  -- fully resolved (v0.7 Mie Dual chain)
    # `_lift_mie_param_to_n_x_4` normalises all four to (n_spec, 4) so
    # we can write columns 2:5 of the 7-parameter
    # [tau_ref, n_r, n_i, r_m, sigma_r, p0, sigma_p] layout uniformly.
    # Replaces an `eachindex` loop on ḟᵗ that overflowed when the
    # matrix form contained n_spec*4 > 6 entries (BoundsError on the
    # docs build's `_jacobian_aod_plot` after the v0.7 YAML migration).
    ω̃̇_block = _lift_mie_param_to_n_x_4(ω̃̇, n, arr_type)
    ḟᵗ_block = _lift_mie_param_to_n_x_4(ḟᵗ, n, arr_type)
    ω̃̇ = arr_type(zeros(n, 7));  ω̃̇[:, 2:5] .= ω̃̇_block
    ḟᵗ = arr_type(zeros(n, 7));  ḟᵗ[:, 2:5] .= ḟᵗ_block

    # Forward modified properties
    fω = fᵗ .* ω̃
    τ_mod = (1 .- fω) .* τAer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fω)

    # Allocate linearized outputs
    τ̇_mod = arr_type(zeros(n,7)) #similar(τ̇Aer)
    ϖ̇_mod = arr_type(zeros(n,7)) #similar(ω̃̇)

    # Bug 18 fix: derivatives w.r.t. 7 aerosol sub-params [τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]
    τ̇_mod[:,1] .= (1 .- fω) .* τ̇Aer[:,1]
    ϖ̇_mod[:,1] .= 0.0
    # Vectorized form over iparam dimension
    # Dimensions: iparam × spectral
    tmp = fᵗ .* ω̃̇[:,2:5] .+ ω̃ .* ḟᵗ[:,2:5]  # (nSpec, iparam)
    τ̇_mod[:,2:5] .= (1 .- fω) .* τ̇Aer[:,2:5] .- tmp .* τAer  # (nSpec, iparam)
    ϖ̇_mod[:,2:5] .= (ω̃̇[:,2:5] .* (1 .- fᵗ) .- ḟᵗ[:,2:5] .* (ω̃ .* (1 .- ω̃))) ./ (1 .- fω).^2

    τ̇_mod[:,6:7] .= (1 .- fω) .* τ̇Aer[:,6:7]
    ϖ̇_mod[:,6:7] .= 0.0
    
    return CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺,
                                           AerZ₀⁺, AerZ₀⁻),
        CoreScatteringOpticalPropertiesLin(τ̇_mod, ϖ̇_mod, AerŻ⁺⁺, AerŻ⁻⁺,
                                           AerŻ₀⁺, AerŻ₀⁻)
end


# Extract scattering definitions and integrated absorptions for the source function!
"""
    extractEffectiveProps(lods, lods_lin)

Extract effective scattering properties from combined layer optical property arrays.

Computes the cumulative optical depth sums and determines scattering interfaces for
each layer. The scattering interface type (`ScatteringInterface_00`, `01`, `10`, `11`)
controls which adding method variant is dispatched during the interaction step.

# Returns
- `scattering_interfaces_all`: Array of `ScatteringInterface` types per layer.
- `τ_sum_all`: Cumulative optical depth `[nSpec × (Nz+1)]`.
- `τ̇_sum_all`: Cumulative τ derivative `[nSpec × Nparams × (Nz+1)]`.
"""
function extractEffectiveProps(
                lods::Array,#{CoreScatteringOpticalProperties{FT},1}
                lods_lin::Array) #where FT

    FT    = eltype(lods[1].τ)
    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    nParams = size(lods_lin[1].τ̇, 2)
    # First the Scattering Interfaces:
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []
    τ_sum_all = similar(lods[1].τ,(nSpec,nZ+1))
    τ_sum_all[:,1] .= 0
    τ̇_sum_all = similar(lods_lin[1].τ̇,(nSpec,nParams,nZ+1))
    τ̇_sum_all[:,:,1] .= 0
    for iz =1:nZ
        # Need to check max entries in Z matrices here as well later!
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + lods[iz].τ 
        for ip=1:nParams
            τ̇_sum_all[:,ip,iz+1] = τ̇_sum_all[:,ip,iz] + lods_lin[iz].τ̇[:,ip] 
        end
    end
    return scattering_interfaces_all, τ_sum_all, τ̇_sum_all
end

"""
    expandOpticalProperties(in, in_lin, arr_type; expand_Z=false)

Convert optical properties and their linearized derivatives to the target
device array type.  When `expand_Z=false` (default for RT-kernel callers)
singleton-Z arrays (`size(Z,3)==1`) are **not** replicated to `nSpec` copies —
the elemental and fused linearized kernels already branch on `size(Z,3)>1` and
use index `n2=1` when the phase matrix is spectrally flat, so no replication
is needed.  Set `expand_Z=true` only when the result will be concatenated
along the spectral dimension (e.g. the `Base.:*` multi-layer combiner).
"""
function expandOpticalProperties(in::CoreScatteringOpticalProperties,
                                  in_lin::CoreScatteringOpticalPropertiesLin,
                                  arr_type;
                                  expand_Z::Bool = false)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, Z₀⁺, Z₀⁻) = in
    (; τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, Ż₀⁺, Ż₀⁻) = in_lin
    solar_plus = Z₀⁺ === nothing ? nothing : _to_device(arr_type, _ensure_3d(Z₀⁺))
    solar_minus = Z₀⁻ === nothing ? nothing : _to_device(arr_type, _ensure_3d(Z₀⁻))
    solar_dot_plus = Ż₀⁺ === nothing ? nothing : _to_device(arr_type,
        ndims(Ż₀⁺) == 3 ? reshape(Ż₀⁺, size(Ż₀⁺,1), size(Ż₀⁺,2), 1, size(Ż₀⁺,3)) : Ż₀⁺)
    solar_dot_minus = Ż₀⁻ === nothing ? nothing : _to_device(arr_type,
        ndims(Ż₀⁻) == 3 ? reshape(Ż₀⁻, size(Ż₀⁻,1), size(Ż₀⁻,2), 1, size(Ż₀⁻,3)) : Ż₀⁻)
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    @assert length(τ̇) == length(ϖ̇) "τ̇ and ϖ̇ sizes need to match"

    if size(Z⁺⁺, 3) == 1
        if expand_Z
            Z⁺⁺ = _repeat(Z⁺⁺, 1, 1, length(τ))
            Z⁻⁺ = _repeat(Z⁻⁺, 1, 1, length(τ))
            Ż⁺⁺ = _repeat(reshape(Ż⁺⁺, size(Ż⁺⁺,1), size(Ż⁺⁺,2), 1, size(Ż⁺⁺,3)), 1, 1, length(τ), 1)
            Ż⁻⁺ = _repeat(reshape(Ż⁻⁺, size(Ż⁻⁺,1), size(Ż⁻⁺,2), 1, size(Ż⁻⁺,3)), 1, 1, length(τ), 1)
            return CoreScatteringOpticalProperties(
                        _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                        _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
                        solar_plus, solar_minus),
                   CoreScatteringOpticalPropertiesLin(
                        _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇),
                        _to_device(arr_type, Ż⁺⁺), _to_device(arr_type, Ż⁻⁺),
                        solar_dot_plus, solar_dot_minus)
        else
            # Fast path: fused linearized kernels branch on size(Z,3) and
            # size(Ż,3), using n2=1 / n2_lin=1 for singleton spectral dim.
            # Ensure 3-D / 4-D shapes so Z[i,j,1] and Ż[i,j,1,p] are valid.
            Z3⁺⁺ = _to_device(arr_type, _ensure_3d(Z⁺⁺))
            Z3⁻⁺ = _to_device(arr_type, _ensure_3d(Z⁻⁺))
            Ż4⁺⁺ = _to_device(arr_type, ndims(Ż⁺⁺) == 3 ? reshape(Ż⁺⁺, size(Ż⁺⁺,1), size(Ż⁺⁺,2), 1, size(Ż⁺⁺,3)) : Ż⁺⁺)
            Ż4⁻⁺ = _to_device(arr_type, ndims(Ż⁻⁺) == 3 ? reshape(Ż⁻⁺, size(Ż⁻⁺,1), size(Ż⁻⁺,2), 1, size(Ż⁻⁺,3)) : Ż⁻⁺)
            return CoreScatteringOpticalProperties(
                        _to_device(arr_type, τ), _to_device(arr_type, ϖ), Z3⁺⁺, Z3⁻⁺,
                        solar_plus, solar_minus),
                   CoreScatteringOpticalPropertiesLin(
                        _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇), Ż4⁺⁺, Ż4⁻⁺,
                        solar_dot_plus, solar_dot_minus)
        end
    else
        @assert size(Z⁺⁺, 3) == length(τ) "Z and τ dimensions need to match"
        return CoreScatteringOpticalProperties(
                    _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                    _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
                    solar_plus, solar_minus),
               CoreScatteringOpticalPropertiesLin(
                    _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇),
                    _to_device(arr_type, Ż⁺⁺), _to_device(arr_type, Ż⁻⁺),
                    solar_dot_plus, solar_dot_minus)
    end
end
