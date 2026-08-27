"""
    jacobian_plan(flavor, params, model, lin_model)

Compile a retrieval-specific Jacobian flavour into a global named state and
compact band-local tangent layouts. New retrievals extend this function by
multiple dispatch; the radiative-transfer kernels consume only the resulting
generic [`JacobianPlan`](@ref).
"""
function jacobian_plan end

"""
    requires_aerosol_microphysics_jacobians(flavor)

Upstream-work trait. Return `false` only when refractive index and size
distribution are fixed: forward Mie optics are still computed, but their
ForwardDiff tangents and Fourier phase derivatives are omitted.
"""
requires_aerosol_microphysics_jacobians(::AbstractJacobianFlavor) = true
requires_aerosol_microphysics_jacobians(::OCO_RRS_synth) = false

"""
    requires_h2o_jacobians(flavor)

Upstream-work trait controlling the q-driven H2O tangent block. Returning
`false` retains the complete forward H2O absorption while leaving its native
gas derivative columns zero and therefore eligible for omission by the plan.
"""
requires_h2o_jacobians(::AbstractJacobianFlavor) = true
requires_h2o_jacobians(::OCO_RRS_synth) = false

const _OCO_RRS_AEROSOL_NAMES =
    ("sulfate", "organic_carbon", "utls_sulfate")

@inline _oco_psurf_key() =
    ParameterKey(:atmosphere, :surface_pressure)
@inline _oco_gas_key(layer) =
    ParameterKey(:gas, :vmr; component="CO2", layer)
@inline _oco_aerosol_key(iaer, field) =
    ParameterKey(:aerosol, field; component=iaer)
@inline _oco_surface_key(i_band, icoeff) =
    ParameterKey(:surface, Symbol("P$(icoeff - 1)");
                 component=i_band, band=i_band)
@inline _oco_sif_key(field) =
    ParameterKey(:source, field; component="SIF", band=1)

function _oco_parameter_name(key::ParameterKey)
    key.kind === :atmosphere && return "surface_pressure"
    if key.kind === :gas
        return "co2_vmr_layer$(lpad(key.layer, 2, '0'))"
    elseif key.kind === :aerosol
        name = _OCO_RRS_AEROSOL_NAMES[key.component]
        return "$(name)_$(key.field)"
    elseif key.kind === :surface
        bands = ("o2a", "weak_co2", "strong_co2")
        return "$(bands[key.band])_surface_$(key.field)"
    elseif key.kind === :source
        return String(key.field)
    end
    return string(key.kind, '_', key.field)
end

function _oco_active_co2_layers(model)
    z_half = half_level_altitudes(model.profile)
    z_center = (z_half[1:end-1] .+ z_half[2:end]) ./ 2
    # The OCO_RRS_synth prior fixes CO2 above 10 km (zero variance). Those
    # layers remain in the forward absorption but need no tangent column.
    return findall(z -> z < 10, z_center)
end

function _oco_global_keys(n_aerosol::Int, n_band::Int,
                          active_co2_layers::AbstractVector{<:Integer})
    keys = ParameterKey[_oco_psurf_key()]
    append!(keys, (_oco_gas_key(iz) for iz in active_co2_layers))
    append!(keys, (_oco_aerosol_key(iaer, :tau_ref)
                   for iaer in 1:n_aerosol))
    append!(keys, (_oco_aerosol_key(iaer, :z0)
                   for iaer in 1:n_aerosol))
    for i_band in 1:n_band, icoeff in 1:3
        push!(keys, _oco_surface_key(i_band, icoeff))
    end
    push!(keys, _oco_sif_key(:SIF760), _oco_sif_key(:mSIF))
    return keys
end

"""
Compile the selective three-band synthetic OCO/RRS layout. The global basis
contains 30 active physical precursors for the current 16-layer profile; each
band receives a smaller kernel-ordered layout. CO2 layers are selected from
actual geometric center altitudes (`z < 10 km`), rather than from hard-coded
layer numbers.
"""
function jacobian_plan(flavor::OCO_RRS_synth, params, model, lin_model)
    n_band = length(get_spec_bands(model))
    n_band == 3 || throw(ArgumentError(
        "OCO_RRS_synth requires exactly three bands; found $n_band"))
    n_aerosol = CoreRT.n_aerosols(model)
    n_aerosol == 3 || throw(ArgumentError(
        "OCO_RRS_synth requires exactly three aerosols; found $n_aerosol"))
    all(surface_parameter_count(get_surface(model, ib)) == 3 for ib in 1:3) ||
        throw(ArgumentError(
            "OCO_RRS_synth requires three surface coefficients in every band"))

    Nz = length(model.profile.p_full)
    active_co2_layers = _oco_active_co2_layers(model)
    isempty(active_co2_layers) && throw(ArgumentError(
        "OCO_RRS_synth found no active CO2 layers below 10 km"))

    abs_params = params.absorption_params
    abs_params === nothing && throw(ArgumentError(
        "OCO_RRS_synth requires configured CO2 absorption"))
    variable_species = unique(Iterators.flatten(abs_params.variable_molecules))
    ico2 = findfirst(==("CO2"), variable_species)
    ico2 === nothing && throw(ArgumentError(
        "OCO_RRS_synth requires CO2 in variable_molecules"))
    # The native gas block reserves species 1 for q-driven H2O.
    native_co2_species = 1 + ico2

    global_keys = _oco_global_keys(n_aerosol, n_band, active_co2_layers)
    global_names = _oco_parameter_name.(global_keys)
    layouts = ActiveParameterLayout[]

    for i_band in 1:n_band
        n_native_gas = size(lin_model.τ̇_abs[i_band], 1)
        n_native_gas % Nz == 0 || throw(DimensionMismatch(
            "native gas derivative count $n_native_gas is not divisible by Nz=$Nz"))
        native_gas_species = n_native_gas ÷ Nz
        native_co2_species <= native_gas_species || throw(ArgumentError(
            "native CO2 gas block is absent in band $i_band"))
        native = ParameterLayout(aerosol_params=7, n_aerosols=n_aerosol,
                                 n_gases=n_native_gas, n_surface=3)

        keys = ParameterKey[_oco_psurf_key()]
        native_columns = Int[psurf_index(native)]
        for iaer in 1:n_aerosol
            # Native aerosol order is [tau_ref,nr,ni,r_m,sigma_g,z0,sigma0].
            append!(keys, (_oco_aerosol_key(iaer, :tau_ref),
                           _oco_aerosol_key(iaer, :z0)))
            arange = aerosol_range(native, iaer)
            append!(native_columns, (arange[1], arange[6]))
        end

        co2_active_in_band = "CO2" in abs_params.variable_molecules[i_band]
        if co2_active_in_band
            for iz in active_co2_layers
                push!(keys, _oco_gas_key(iz))
                push!(native_columns,
                      gas_layer_index(native, native_co2_species, iz, Nz))
            end
        end
        n_layer = length(keys)

        surface_start = n_layer + 1
        append!(keys, (_oco_surface_key(i_band, icoeff) for icoeff in 1:3))
        surface_columns = surface_start:(surface_start + 2)
        if i_band == 1
            sif_start = length(keys) + 1
            push!(keys, _oco_sif_key(:SIF760), _oco_sif_key(:mSIF))
            sif_columns = sif_start:(sif_start + 1)
        else
            sif_columns = (length(keys) + 1):length(keys)
        end
        names = _oco_parameter_name.(keys)
        push!(layouts, ActiveParameterLayout(
            keys, names, global_keys, global_names, native_columns, n_layer;
            surface_columns, sif_columns))
    end
    return JacobianPlan(flavor, global_keys, global_names, layouts)
end

"""
    model_from_parameters(flavor::AbstractJacobianFlavor, params; kwargs...)

Build the ordinary forward/linearized optical models, compile `flavor`, and
return a `PlannedRTModelLin` that activates compact propagation in `rt_run`.
"""
function model_from_parameters(flavor::AbstractJacobianFlavor, params; kwargs...)
    defaults = (
        compute_aerosol_microphysics_jacobians =
            requires_aerosol_microphysics_jacobians(flavor),
        compute_h2o_jacobians = requires_h2o_jacobians(flavor),
    )
    options = merge(defaults, (; kwargs...))
    model, base_lin = model_from_parameters(LinMode(), params; options...)
    plan = jacobian_plan(flavor, params, model, base_lin)
    return model, PlannedRTModelLin(base_lin, plan)
end

"""
    globalize_jacobian(jacobian, layout)

Scatter a band-local Jacobian into the plan's global named physical basis.
The trailing dimension is replaced by the global dimension, inactive
parameters are initialized to exact zero, and active columns are assigned
through `layout.local_to_global`. This operation reorders/zero-fills only; it
does not apply unit, AOD-reference, or logarithmic-coordinate chain rules.
"""
function globalize_jacobian(jacobian::AbstractArray,
                            layout::ActiveParameterLayout)
    size(jacobian, ndims(jacobian)) == n_total(layout) ||
        throw(DimensionMismatch(
            "Jacobian trailing dimension does not match its active layout"))
    dims = (size(jacobian)[1:end-1]..., length(layout.global_keys))
    result = similar(jacobian, eltype(jacobian), dims)
    fill!(result, zero(eltype(result)))
    dest = (ntuple(_ -> Colon(), ndims(result) - 1)...,
            layout.local_to_global)
    result[dest...] .= jacobian
    return result
end
